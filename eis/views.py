#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext
from time import time
from hashlib import md5
import os
from django.conf import settings
from django.core.urlresolvers import reverse
from django.contrib import messages
import socket
from stat import ST_MTIME, S_IRWXU, S_IRWXG
from sdftools.moleculeformats import smiles_to_sdf, InputError
from eis import InvalidInputError, first_mol
from tempfile import mkstemp
from subprocess import Popen
import sys
from copy import copy
from guest.decorators import guest_allowed, login_required
from pugsearch.views import search as pug_search
from simplejson import dumps


# try:
# ....from renderer import do_post as sdf2png
# except:
# ....sys.path.append('/srv/renderer/bin')
# ....from renderer import do_post as sdf2png

@guest_allowed
def read(request, s):
    fp = os.path.join(settings.WORK_DIR, s) + '.j'
    if not os.path.exists(fp):
        msg = {'status': 'wait'}
    else:
        f = file(fp)
        c = f.read()
        f.close()
        if not c:
            msg = {'status': 'wait'}
        else:
            msg = {'status': 'ok', 'content': c}
    return HttpResponse(dumps(msg), mimetype='text/plain')


@guest_allowed
def search(request):
    if request.method == 'GET':
        smi = ''
        if 'smi' in request.GET:
            smi = str(request.GET['smi'])
        return render_to_response('eis/search.html', dict(mode='form',
                                  smi=smi),
                                  context_instance=RequestContext(request))
    else:
        if 'smiles' in request.POST:
            try:
                sdf = smiles_to_sdf(str(request.POST['smiles']))
            except InputError:
                messages.error(request, 'Error: Invalid SMILES string!')
                return HttpResponseRedirect(reverse('compound_search'))
        elif 'sdf' in request.FILES:
            try:
                sdf = first_mol(request.FILES['sdf'])
            except:
                messages.error(request, 'Error: Invalid input data!')
                return HttpResponseRedirect(reverse('compound_search'))
        elif 'sdf' in request.POST:
            try:
                sdf = request.POST['sdf']
                sdf = first_mol(request.POST['sdf'])
            except:
                messages.error(request, 'Error: Invalid input data!')
                return HttpResponseRedirect(reverse('compound_search'))
        else:
            messages.error(request, 'Error: Invalid input data!')
            return HttpResponseRedirect(reverse('compound_search'))

        # dos2unix, and remove 2nd line (fpcdb complains JME editor SDF)

        sdf = sdf.split('\r\n')
        sdf = '\n'.join(sdf)
        sdf = sdf.split('\n')
        sdf[1] = ''
        sdf = '\n'.join(sdf)

        if not sdf.strip():
            return HttpResponseRedirect(request.get_full_path())
        else:
            s = md5(sdf).hexdigest()
            dir = os.path.join(settings.WORK_DIR, s)
            if not os.path.exists(dir):
                os.mkdir(dir)
                os.chmod(dir, S_IRWXU | S_IRWXG)
            query_sdf = os.path.join(dir, 'query.sdf')
            if os.path.exists(query_sdf):
                os.unlink(query_sdf)
            f = file(query_sdf, 'w')
            f.write(sdf)
            f.close()

            # s_ = sdf2png('file://' + os.path.join(dir, 'query.sdf'),
            # ............local=True)
            # assert s == s_

            # remove status file if none exists (note this could yield race condition)

            filename = os.path.join(dir, 'status')
            if os.path.exists(filename):
                os.unlink(filename)
            filename = os.path.join(dir, 'out')
            if os.path.exists(filename):
                os.unlink(filename)

            # if pubchem PUG search

            if 'pug' in request.POST:
                try:
                    return pug_search(request, s)
                except:
                    messages.error(request, 'No hits!')
                    return HttpResponseRedirect(reverse('compound_search'
                            ))

            # submit to server

            (HOST, PORT) = settings.ASSEI_SERVER
            try:
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                sock.connect((HOST, PORT))
                hello = sock.recv(3)
                assert hello == 'OK\n'
                sock.send(dir + '\n')
            except:
                messages.error(request, 'Query server not reachable!')
                return HttpResponseRedirect(reverse('compound_search'))

            return HttpResponseRedirect(reverse('search_wait',
                    kwargs=dict(sid=s)))


@guest_allowed
def wait(request, sid):
    return render_to_response('eis/wait.html', dict(sid=sid),
                              context_instance=RequestContext(request))


def session_check(f):

    def func(
        request,
        sid,
        *args,
        **kargs
        ):

        dir = os.path.join(settings.WORK_DIR, sid)
        query = os.path.join(dir, 'query.sdf')
        if not os.path.exists(dir):
            messages.error(request,
                           'No such job found. Please try again.')
            return HttpResponseRedirect(reverse('compound_search'))

        status_fp = os.path.join(dir, 'status')
        if os.path.exists(status_fp):
            status_f = file(status_fp)
            status = status_f.read().strip()
            status_f.close()
        if not os.path.exists(status_fp) or status != 'done' and status \
            != 'failed':
            if time() - os.stat(query)[ST_MTIME] \
                > settings.QUERY_TIMEOUT:
                if 'ajax' in request.GET:
                    return HttpResponse('{"status":"ok"}',
                            mimetype='text/plain')
                return render_to_response('eis/result.html',
                        dict(sid=sid, mode='perm_failure'),
                        context_instance=RequestContext(request))
            else:
                if 'ajax' in request.GET:
                    return HttpResponse('{"status":"wait"}',
                            mimetype='text/plain')
                return HttpResponseRedirect(reverse('search_wait',
                        kwargs=dict(sid=sid)))

        if 'ajax' in request.GET:
            return HttpResponse('{"status":"ok"}', mimetype='text/plain'
                                )

        if status == 'failed':
            error_fp = os.path.join(dir, 'error')
            error_f = file(error_fp)
            error = error_f.read()
            error_f.close()
            return render_to_response('eis/result.html', dict(sid=sid,
                    mode='error', error=error),
                    context_instance=RequestContext(request))
        elif status == 'done':

            return f(request, sid, *args, **kargs)

    return func


@guest_allowed
@session_check
def result(request, sid):
    dir = os.path.join(settings.WORK_DIR, sid)
    out_fp = os.path.join(dir, 'out')
    out_f = file(out_fp)
    line = out_f.next()
    assert line.startswith('#')
    timing = [float(i) for i in line[1:].strip().split()]
    total_t = reduce(lambda x, y: x + y, timing, 0)
    r_timing = [int(i * 1000) for i in timing]
    t_table_width = int(total_t * 1000) + 300
    compounds = []
    for line in out_f:
        (c, d) = line.split()
        compounds.append((c, 1 - float(d)))
    out_f.close()

    return render_to_response('eis/result.html', dict(
        sid=sid,
        mode='ready',
        timing=timing,
        compounds=compounds,
        r_timing=r_timing,
        t_table_width=t_table_width,
        total_t=total_t,
        slow=total_t > 1.5,
        ), context_instance=RequestContext(request))


@guest_allowed
@session_check
def download(request, sid, format):
    dir = os.path.join(settings.WORK_DIR, sid)
    out_fp = os.path.join(dir, 'out')
    out_f = file(out_fp)
    out_f.next()
    cids = [int(line.split()[0]) for line in out_f]
    out_f.close()
    cids_str = '\n'.join(['%s' % i for i in cids])
    s = md5(cids_str + format).hexdigest()
    fp = os.path.join(settings.WORK_DIR, s) + '.j'
    (tmp_fd, tmp_fp) = mkstemp(dir=settings.WORK_DIR)
    fo = os.fdopen(tmp_fd, 'w')
    fo.write(cids_str)
    fo.close()
    cmd = '%s %s %s %s' % (settings.PUBCHEM_DOWNLOADER, tmp_fp, fp,
                           format)
    Popen(cmd, close_fds=True, shell=True)

    return render_to_response('compound_download.html', dict(sid=s),
                              context_instance=RequestContext(request))

