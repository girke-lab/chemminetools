#!/usr/bin/python

from builtins import str
import re
import time
#from string import join
from django.shortcuts import redirect, render
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.urls import reverse
from django.contrib import messages
from guest.decorators import guest_allowed, login_required
from django.template import RequestContext
from django.utils.http import urlunquote, urlquote
from tools.models import Application, Job, ApplicationCategories
from tools.runapp import updateJob, createJob, getAppForm, parseToolForm
import tools
import structure_search 
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_sdf, \
    InputError, sdf_to_smiles
from structure_search import first_mol
import traceback
from compounddb.models import Compound, SDFFile,Tag
import csv
import functools

from django.forms import Form,  ModelChoiceField

@guest_allowed
def search(request):
    if request.is_ajax():
        if 'application_id' in request.GET:
            form = getAppForm(request.GET['application_id'], request.user)()
            return HttpResponse(str(form))
        else:
            return HttpResponse("",status_code=404)
    elif request.method != 'POST':
        smi = ''
        if 'smi' in request.GET:
            smi = str(request.GET['smi'])
            smi = urlunquote(smi)
        allTags = Tag.allUserTagNames(request.user)


        fields = {}
        category = ApplicationCategories.objects.get(name="Search")
        apps = Application.objects.filter(category=category).exclude(name__exact="pubchemID2SDF")
        fields['application'] = ModelChoiceField(queryset=apps, empty_label='')
        form = type('%sForm' % 'choose application', (Form, ), fields)


        return render(request,'search.html', dict(mode='form',
            smi=smi, appsForm=form, tags=allTags))
    else:
        sdf = None
        smiles = None
        compid = u'query'
        form = None
        application = None
        application_id = None

        if 'application' in request.POST:
            application_id = request.POST['application']
            application = Application.objects.get(id= application_id)


        if 'tags' in request.POST:
            givenTags = request.POST.getlist("tags")
            compoundList = Compound.byTagNames(givenTags,request.user)
            if len(compoundList) == 0:
                messages.error(request,"Error: No compounds found with selected tags")
            else:
                sdf = u''
                for compound in compoundList:
                    sdf = sdf + compound.sdffile_set.all()[0].sdffile.rstrip() + '\n'
                smiles = sdf_to_smiles(sdf)
        elif 'smiles' in request.POST:
            input_mode = 'smiles-input'
            sdf = u''
            try:
                smiles = request.POST['smiles']
                sdf = smiles_to_sdf(smiles)
            except:
                messages.error(request, 'Error: Invalid SMILES string!')
                sdf = None
        elif 'sdf' in request.FILES:
            input_mode = 'sdf-upload'
            try:
                sdf = request.FILES['sdf']
                sdf = first_mol(sdf.read())
                smiles = sdf_to_smiles(sdf)
            except:
                print(traceback.format_exc())
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
                    print(traceback.format_exc())
                    messages.error(request, 'Invalid drawing!')
                    sdf = None
            else:
                try:
                    input_mode = 'sdf-input'
                    sdf = first_mol(request.POST['sdf'])
                    smiles = sdf_to_smiles(sdf)
                except:
                    print(traceback.format_exc())
                    messages.error(request, 'Invalid input SDF!')
                    sdf = None

        if application_id != None:
            AppFormSet = getAppForm(request.POST['application'], request.user)
            form = AppFormSet(request.POST)

        if form != None and form.is_valid():
            commandOptions, optionsList = parseToolForm(form) 
        else:
            sdf = None
            messages.error(request, "Invalid form options!")
        if not sdf:
            print("no sdf found")
            return redirect(structure_search.views.search)
        smiles = re.search(r'(\S+)', smiles).group(1)
        smiles = urlquote(smiles)


        #print("command options: "+str(commandOptions))

        #if algorithm == u'fp':
        if application.name == "PubChem Fingerprint Search":
            newJob = createJob(request.user, application.name, optionsList, commandOptions, sdf, smiles)
        #elif algorithm == "fp-chembl" :
        elif application.name == "ChEMBL Fingerprint Search":
            newJob = createJob(request.user, application.name, optionsList, commandOptions, sdf, smiles)
        #elif algorithm == "ei":
        elif application.name == "EI Search":
            newJob = createJob(request.user, application.name, optionsList, commandOptions, sdf, smiles)

        time.sleep(1)
        return redirect(tools.views.view_job, job_id=newJob.id,resource='')

@guest_allowed
def getStructures(request, job_id, format):

    # takes a search job ID and returns an SDF of the compounds from this result 

    try:
        job = updateJob(request.user, job_id)
        f = open(job.output, 'r')
        csvinput = csv.reader(f,delimiter=' ')
        #read each line, extrace the second column, and combine all values 
        #into a new-line separated string
        targetIds = "\n".join([line[0] for line in csvinput])
        print("targetIds: "+str(targetIds))
        f.close()
        #result = '\n'.join(re.findall(r'^\S+', result, re.MULTILINE))
    except Job.DoesNotExist:
        print(traceback.format_exc())
        raise Http404
    newJob = createJob(request.user, 'pubchemID2SDF', '', [], targetIds,
                       format, async=False) 
    if format == 'smiles':
        filename = 'search_result.smi'
    else:
        filename = 'search_result.sdf' 
    return redirect(tools.views.view_job, job_id=newJob.id,resource='other',filename=filename)
