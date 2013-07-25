#!/usr/bin/python

import re
import time
from string import join
from django.shortcuts import redirect, render_to_response
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.contrib import messages
from guest.decorators import guest_allowed, login_required
from django.template import RequestContext
from tools.models import Application, Job
from tools.runapp import updateJob, createJob, getAppForm, parseToolForm
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_sdf, \
    InputError, sdf_to_smiles

@guest_allowed
def search(request):
    application = Application.objects.get(name='EI Search')
    AppFormSet = getAppForm(application.id, request.user)
    if request.method != 'POST':
        smi = ''
        if 'smi' in request.GET:
            smi = str(request.GET['smi'])
        form = AppFormSet()
        form = str(form)
        return render_to_response('search.html', dict(mode='form',
            smi=smi, form=form),
            context_instance=RequestContext(request)) 
    else:
        sdf = None
        smiles = None
        compid = u'query'
        if 'smiles' in request.POST:
            input_mode = 'smiles-input'
            sdf = u''
            try:
                smiles = request.POST['smiles']
                sdf = smiles_to_sdf(str(smiles))
            except:
                messages.error(request, 'Error: Invalid SMILES string!')
                sdf = None
        elif 'sdf' in request.FILES:
            input_mode = 'sdf-upload'
            try:
                sdf = request.FILES['sdf']
                sdf = sdf.read()
                smiles = sdf_to_smiles(sdf)
            except (InputError, InvalidInputError):
                messages.error(request, 'Invalid SDF!')
                sdf = None
        elif 'sdf' in request.POST:
            if 'draw' in request.POST:
                input_mode = 'draw'
                sdf = request.POST['sdf'] + '$$$$'
                try:
                    smiles = sdf_to_smiles(sdf)
                    smiles = re.match(r"^(\S+)", smiles).group(1)
                    smiles = smiles + ' ' + compid
                    sdf = smiles_to_sdf(smiles)
                except:
                    messages.error(request, 'Invalid drawing!')
                    sdf = None
            else:
                try:
                    input_mode = 'sdf-input'
                    sdf = request.POST['sdf']
                    smiles = sdf_to_smiles(sdf)
                except:
                    messages.error(request, 'Invalid input SDF!')
        form = AppFormSet(request.POST)
        if form.is_valid():
            commandOptions, optionsList = parseToolForm(form) 
        else:
            sdf = None
            messages.error(request, "Invalid form options!")
        if not sdf:
            return redirect('eisearch.views.search')
        smiles = re.search(r'(\S+)', smiles).group(1)
        newJob = createJob(request.user, 'EI Search', optionsList, commandOptions, sdf, smiles)
        time.sleep(2)
        return redirect('tools.views.view_job', job_id=newJob.id,resource='')

@guest_allowed
def getStructures(request, job_id, format):

    # takes a search job ID and returns an SDF of the compounds from this result 

    try:
        job = updateJob(request.user, job_id)
        f = open(job.output, 'r')
        result = f.read()
        f.close() 
        result = join(re.findall(r'^\S+', result, re.MULTILINE), sep='\n')
    except Job.DoesNotExist:
        raise Http404
    # note: add in here a regex to keep only cids
    newJob = createJob(request.user, 'pubchemID2SDF', '', '', result, format) 
    time.sleep(2)
    return redirect('tools.views.view_job', job_id=newJob.id,resource='')
