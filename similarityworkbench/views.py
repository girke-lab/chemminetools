#!/usr/bin/python
# -*- coding: utf-8 -*-

from builtins import str
from django.shortcuts import render
from django.template import RequestContext
from django.http import Http404, HttpResponseBadRequest, HttpResponse
from simplejson import dumps
from sdftools.moleculeformats import batch_smiles_to_smiles, \
    batch_sdf_to_smiles
from similarityworkbench import funcs
from myCompounds.views import getMyCompounds
from guest.decorators import guest_allowed, login_required
from django.views.decorators.cache import cache_page
from django.utils.http import urlunquote, urlquote
import pybel
import re

import os
import tempfile
import sys
import traceback
from django.conf import settings


@guest_allowed
@cache_page(60 * 120)
def renderer(request, smiles):
    try:
        smiles = urlunquote(smiles)
        smiles = re.match(r"^(\S{1,2000})", smiles).group(1)
        mymol = pybel.readstring('smi', smiles)
        #this simple method below does not work as it returns the image as a unicode string
        # which cannot be converted back to raw bytes
        ##png = mymol.write(format='_png2')

        # so we write the image to a temp file and then read it back as bytes
        fp = tempfile.NamedTemporaryFile()
        fp.close() # deletes the file, but we still use the unique file name
        mymol.write(format='_png2',filename=fp.name)
        imageFile = open(fp.name,"rb")
        image = imageFile.read()
        imageFile.close()
        os.remove(fp.name)
        return HttpResponse(image, content_type='image/png')
    except:
        print("error rendering "+smiles)
        print("Unexpected error:", sys.exc_info())
        traceback.print_tb(sys.exc_info()[2])
        raise Http404


def get_workbench_compounds(request):
    print("getting workbench compounds")
    ret = []
    from hashlib import md5
    (page, compounds) = getMyCompounds(request)
    for compound in compounds:
        try:
            md5id = md5(compound.cid.encode()).hexdigest()
            img = '/compounds/' + str(compound.id) + '/png'
            ret.append(dict(img=img, md5=md5id, title=compound.cid,
                       smiles=compound.smiles))
        except:
            print("error while processing compound: ")
            print("Unexpected error:", sys.exc_info())
            traceback.print_tb(sys.exc_info()[2])
            raise
    return ret


def add_compounds(smiles):
    """given SMILES, generate a list of dictionaries, one for each compound,
....in format that is accepted at the client side"""

    ret = []
    from hashlib import md5  # for md5 ID of compounds
    for (i, c) in enumerate(smiles.splitlines()):
        md5id = md5(c.encode()).hexdigest()
        try:
            (s, name) = c.split(None, 1)
        except:
            s = c
            name = 'Unnamed:' + md5id[:5]
        img = '/similarity/renderer/' + s
        ret.append(dict(img=img, md5=md5id, title=name, smiles=c))
    return ret


@guest_allowed
def ui(request):
    if request.method == 'GET':

        # on GET, show the UI page
        print("loading similarity page")

        try:
            preload = dumps(get_workbench_compounds(request))
        except:
            preload = ''
        return render(request,'similarityworkbench.html',
                                  dict(preload=preload))
    else:

        # on POST, do actual work; the form tells what function to invoke

        func = request.POST.get('func', None)

        # there can be overrring function in the form. mostly javascript uses
        # this to override form-specified function

        if 'func-override' in request.POST:
            func = request.POST.get('func-override')

        # first type of function is compound upload

        if func in ['sdf', 'smiles']:
            if func == 'smiles':

                # process SMILES upload
                # first, format the SMILES

                (smiles, err) = \
                    batch_smiles_to_smiles(request.POST.get('smiles', ''
                        ))
            elif func == 'sdf':

                # process SDF upload

                sdff = request.FILES.get('sdffile', None)
                if not sdff or not sdff.size:
                    return HttpResponseBadRequest()
                sdf = sdff.read()
                (smiles, err) = batch_sdf_to_smiles(sdf)

            # for each compound passed the format test, let the client know

            status = add_compounds(smiles)

            # also notify about the failed ones

            failed = err
            return HttpResponse(dumps(dict(failed=failed,
                                status=status)), content_type='text/json')
        elif hasattr(funcs, func):

        # the other type is actually processing request. Check whether the
        # the requested function is supported
            # there is correponding processor function in the funcs module

            f = getattr(funcs, func)
            ret = f(*request.POST.getlist('compound'))
            return HttpResponse(dumps(ret), content_type='text/json')
        else:

        # we don't know how to process this request

            return HttpResponseBadRequest()
