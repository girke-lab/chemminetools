import sys
import os
from pugsearch.PUG_services import *
from time import sleep, time
from django.conf import settings
from django.http import Http404, HttpResponseBadRequest, HttpResponseRedirect
from compounddb.models import Compound
from django.shortcuts import render_to_response
import logging
from django.template import RequestContext
from django.core.urlresolvers import reverse

loc = PUGLocator()
port = loc.getPUGSoap()

def search(request, sid):
    try:
        cutoff = int(request.POST['cutoff'])
    except:
        if not 'session_msg' in request.session:
        	request.session['session_msg'] = []
        request.session['session_msg'].append("Invalid cutoff!")
        request.session.save()
        return HttpResponseRedirect(reverse('compound_search'))
    if (cutoff < 90) or (cutoff > 99):
	if not 'session_msg' in request.session:
                request.session['session_msg'] = []
        request.session['session_msg'].append("Cutoff must be between .90 and .99")
        request.session.save()
        return HttpResponseRedirect(reverse('compound_search'))
    dir = os.path.join(settings.WORK_DIR, sid)
    query = os.path.join(dir, 'query.sdf')
    if not os.path.exists(query):
        raise Http404, 'no such job: %s' % sid
    try:
        f = file(os.path.join(dir, 'pugsearch.results.%d' % cutoff))
        elapse = float(f.next())
        url = f.next().strip()
        id_list = [int(i) for i in f]
        f.close()
        logging.info("successfully loading cached results for sid:%s"%sid)
    except:
        req = InputStructureBase64SoapIn()
        req.set_element_structure(file(query).read())
        req.set_element_format('eFormat_SDF')

        start = time()
        strKey = port.InputStructureBase64(req).get_element_StrKey()
        
        req = SimilaritySearch2DSoapIn();
        req.set_element_StrKey(strKey)
        req.set_element_simOptions(req.new_simOptions())
        req.get_element_simOptions().set_element_threshold(cutoff)
        listKey = port.SimilaritySearch2D(req).get_element_ListKey()
        logging.info("PubChem job key: %s" % listKey)
        
        # call GetOperationStatus until the operation is finished
        req = GetOperationStatusSoapIn()
        req.set_element_AnyKey(listKey)
        status = port.GetOperationStatus(req).get_element_status()
        while (status == 'eStatus_Queued' or status == 'eStatus_Running'):
            sleep(1)
            status = port.GetOperationStatus(req).get_element_status()
            
        logging.info("PUG search done for %s" % listKey)
        elapse = time() - start
        warning = ''
        # check status
        if (status == 'eStatus_Success'
                or status == 'eStatus_TimeLimit'
                or status == 'eStatus_HitLimit'): 

            if (status == 'eStatus_TimeLimit'):
                warning = 'Warning: time limit reached before entire db searched'
            elif (status == 'eStatus_HitLimit'):
                warning = 'Warning: hit limit reached before entire db searched'

            # get Entrez key
            req = GetEntrezKeySoapIn()
            req.set_element_ListKey(listKey);
            entrezKey = port.GetEntrezKey(req).get_element_EntrezKey()
            
            # get URL
            req = GetEntrezUrlSoapIn()
            req.set_element_EntrezKey(entrezKey)
            url = port.GetEntrezUrl(req).get_element_url()

            # get ID list
            req = GetIDListSoapIn()
            req.set_element_ListKey(listKey)
            id_list = port.GetIDList(req).get_element_IDList().get_element_int()[:200]
            f = file(os.path.join(dir, 'pugsearch.results.%d' % cutoff), 'w')
            f.write('%f\n' % elapse)
            f.write('%s\n' % url)
            for i in id_list: f.write('%d\n' % i)
            f.close()

        else:   # status indicates error
        	if not 'session_msg' in request.session:
                	request.session['session_msg'] = []
        	request.session['session_msg'].append("No hits! Please try the EI search, or select a lower cutoff.")
        	request.session.save()
        	return HttpResponseRedirect(reverse('compound_search'))
	
		# commented out, don't let users see these TERRIBLE error messages
        	   # see if there's some explanatory message
            	# req = GetStatusMessageSoapIn()
            	# req.set_element_AnyKey(listKey)
            	# error = port.GetStatusMessage(req).get_element_message()
            	# return render_to_response(
                # 'pugsearch/result.html',
                # dict(mode='error', error=error),
                # context_instance=RequestContext(request)
                # )
    # filter IDs
    # _id_list = ['%d'%i for i in id_list]
    # from django.db import connection
    # cursor = connection.cursor()
    # cursor.execute("SELECT cid FROM eic_compound WHERE cid IN %s" %
    #    str(tuple(_id_list)))

    # covered_ids = [int(row[0]) for row in cursor.fetchall()]
    #covered_ids = Compound.objects.filter(cid__IN=_id_list).values_list('cid')

    # write IDs to file for downloader to access
    filename = os.path.join(dir, 'out')
    if os.path.exists(filename):
        os.unlink(filename)
    f = file(filename, 'w')
    f.write('#\n')
    for cid in id_list:
        f.write('%s %s\n' % (cid, 0))
    f.close()

    # write status file for PUG download
    filename = os.path.join(dir, 'status')
    if os.path.exists(filename):
        os.unlink(filename)
    f = file(filename, 'w')
    f.write('done\n')
    f.close()

    # load EI-Search result
    # out_fp = os.path.join(dir, 'out')
    # out_f = file(out_fp)
    # line = out_f.next()
    # assert line.startswith('#')
    # eis_time = sum([float(i) for i in line[1:].split()])
    # hits = dict([(int(i.split()[0]), 1-float(i.split()[1])) for i in out_f])

    # compounds = [(i, hits.get(i, None), i not in covered_ids, i not in hits and i in covered_ids) for i in id_list]
    compounds = [(i,None,None,None,None) for i in id_list]
    return render_to_response(
        'pugsearch/result.html',
        dict(mode='ready', sid=sid, compounds=compounds, elapse=elapse,
            eis_time=0, e_url=url),
        context_instance=RequestContext(request)
        )
        

# vim:tabstop=4:shiftwidth=4:expandtab
