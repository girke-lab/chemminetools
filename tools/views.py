#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from builtins import next
from builtins import str
import re
import os
import csv
import time
import traceback
import myCompounds
import tools
from collections import defaultdict, OrderedDict
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import redirect, render
from django.template import RequestContext
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.forms import ModelForm
from django.contrib import messages
from django.utils.http import urlquote
from guest.decorators import guest_allowed, login_required
from compounddb.models import Compound, Tag
from tools.runapp import *
from .models import *
from sdftools.moleculeformats import batch_sdf_to_smiles
from simplejson import dumps

@guest_allowed
def launch_job(request, category=None):
    if request.is_ajax():

        # for ajax requests, return HTML form for each app

        currentApp = request.GET['currentApp']
        try:
            application = Application.objects.get(id=currentApp)
            AppFormSet = getAppForm(application.id, request.user)
            form = AppFormSet()
            form = str(form)
            response = dict(form=form, desc=application.description)
        except Exception as e:
            print(traceback.format_exc())
            response = dict(form='ERROR')
        return HttpResponse(dumps(response), 'text/json')
    if request.method == 'POST':
        appForm = getAppForm(request.POST['application'], request.user)
        form = appForm(request.POST, request.FILES)
        if form.is_valid():
            try:
                appid = int(form.cleaned_data['application'])
                application = Application.objects.get(id=str(appid))
            except Application.DoesNotExist:
                messages.error(request, 'Application does not exist')
                return redirect(launch_job, category=category)
        else:
            messages.error(request, str(form.errors))
            return redirect(launch_job, category=category)
        commandOptions, optionsList = parseToolForm(form)

        # setup input

        if application.input_type == 'chemical/x-mdl-sdfile':
            input = 'chemical/x-mdl-sdfile'
        elif application.input_type == 'upload':
            input = request.FILES['File Upload'].read().decode("utf-8")
        else:
            input = ''
        tagNames = []
        if "tags" in request.POST:
            tagNames = request.POST.getlist("tags")
        #print("got tagNames: "+str(tagNames))
        newJob = createJob(request.user, application.name, optionsList,
                           commandOptions, input,tagNames=tagNames)
        messages.success(request, 'Success: job launched.')
        return redirect(view_job, job_id=newJob.id, resource='')
    else:
        if category:
            fromWorkbench = False
            try:
                category = ApplicationCategories.objects.get(name=category)
                compoundCount = Compound.objects.filter(user=request.user).count()
                if category.name != 'Upload':
                    fromWorkbench = True
                if category.name == 'Clustering' and compoundCount < 3:
                    messages.info(request,
                                  'Notice: you must have at least 3 compounds to perform clustering. Please use this form to add more compounds and then click "Cluster" again.'
                                  )
                    return redirect(myCompounds.views.uploadCompound,job_id=None, resource=None)
                if category.name == 'Properties' and compoundCount < 1:
                    messages.info(request,
                                  'Notice: you must have at least one compound to compute properties. Please use this form to add compounds and then click "Properties" again.'
                                  )
                    return redirect(myCompounds.views.uploadCompound,job_id=None, resource=None)
                title = 'Launch ' + category.name + ' Job'
                apps = Application.objects.filter(category=category)
            except:
                print(traceback.format_exc())
                raise Http404
        else:
            title = 'Launch Job'
            apps = Application.objects.filter()
        fields = {}
        fields['application'] = ModelChoiceField(queryset=apps, empty_label='')
        form = type('%sForm' % 'choose application', (Form, ), fields)
        allTags = Tag.allUserTagNames(request.user)
        return render(request,'submitForm.html', dict(title=title,
                                  form=form,
                                  tags = allTags,
                                  fromWorkbench=fromWorkbench,
                                  totalCompounds=Compound.objects.filter(user=request.user).count()))


@guest_allowed
def view_job(
    request,
    job_id,
    resource=None,
    filename=None,
    ):

    job = updateJob(request.user, job_id)
    if job == False:
        raise Http404
    if resource:
        if resource == 'delete':
            deleteJob(request.user, job.id)
            return HttpResponse('deleted', content_type='text/plain')
        if resource == 'download':
            f = open(job.output, 'r')
            result = f.read()
            f.close()
            if job.application.output_type == 'application/json.canvasxpress':
                result = re.search(r'smps":"(.*)"}},', result).group(1)
            return HttpResponse(result,
                             content_type=job.application.output_type)
    if request.is_ajax():
        if job.status == Job.RUNNING:
            response = dict(status='RUNNING',job_id=job_id)
        else:
            response = dict(status='DONE',job_id=job_id)
        return HttpResponse(dumps(response), 'text/json')
    if job.status == Job.FINISHED:
        finalResult = job.output
        finalResult = re.sub('.*/', '', finalResult, count=0)
        finalResult = '/working/' + finalResult

        job.start_time = str(job.start_time)

        # select correct viewer here based on output type

        if job.application.output_type == 'chemical/sdfUpload':
            f = open(job.output, 'r')
            sdf = f.read()
            f.close()
            result, error = batch_sdf_to_smiles(sdf) 
            quotedSmiles = []
            for smiles in result.splitlines():
                match = re.search(r'(\S+)\s+(\S+)', smiles)
                smiOnly = match.group(1)
                cid = match.group(2)
                quotedSmiles.append({'smiles': urlquote(smiOnly), 'cid':cid})
            return render(request,'sdfUpload.html',
                    dict(title=str(job.application) + ' Results',
                    job=job, compounds=quotedSmiles))
        if job.application.output_type == 'text/ei.search.result':
            f = open(job.output, 'r')
            csvinput = csv.reader(f, delimiter=' ')
            csvOutput = []
            queryColumnEmpty=True
            for line in csvinput:
                if line[0] != "":
                    queryColumnEmpty =False
                #print("csv line: "+str(line))
                csvOutput.append(line)
            #print("query column empty? "+str(queryColumnEmpty))
            f.close()
            allTags = Tag.allUserTagNames(request.user)
            #print("job input: "+str(job.input))
            if queryColumnEmpty:
                singleQuery = job.input
            else:
                singleQuery = None
            return render(request,'structure_search_result.html',
                    dict(title=str(job.application) + ' Results',
                    job=job,tags=allTags, compounds=csvOutput, 
                    singleQuery=singleQuery,resultType="ei"))
        if job.application.output_type == 'text/fp.search.result':
            f = open(job.output, 'r')
            csvinput = csv.reader(f, delimiter=' ')
            csvOutput = []
            queryColumnEmpty=True
            for line in csvinput:
                if line[0] != "":
                    queryColumnEmpty =False
                #print("csv line: "+str(line))
                csvOutput.append(line)
            #print("query column empty? "+str(queryColumnEmpty))
            f.close()

            allTags = Tag.allUserTagNames(request.user)

            if queryColumnEmpty:
                singleQuery = job.input
            else:
                singleQuery = None

            return render(request,'structure_search_result.html',
                    dict(title=str(job.application) + ' Results',
                    job=job, compounds=csvOutput, 
                    singleQuery=singleQuery,tags=allTags,
                    resultType="fp"))
        elif job.application.output_type == 'text/sdf.upload':
            f = open(job.output, 'r')
            message = f.read()
            f.close()
            deleteJob(request.user, job.id)
            if re.search(r"^ERROR:", message):
                messages.error(request, message)
                return redirect(myCompounds.views.uploadCompound,job_id=None, resource=None)
            else:
                messages.success(request, message)
                return redirect(myCompounds.views.showCompounds,
                                resource='')
        elif job.application.output_type \
            == 'application/json.cytoscape':
            f = open(job.output, 'r')
            plotJSON = f.read()
            f.close()
            return render(request,'view_network.html',
                    dict(title=str(job.application) + ' Results',
                    result=finalResult, job=job, plotJSON=plotJSON))
        elif job.application.output_type \
            == 'application/json.canvasxpress':
            f = open(job.output, 'r')
            plotJSON = f.read()
            f.close()
            if (job.application.name == 'Hierarchical Clustering')\
                or (job.application.name == 'Numeric Data Clustering'):
                tree = True 
                return render(request,'view_job.html',
                        dict(title=str(job.application) + ' Results',
                        result=finalResult, tree=tree, job=job, plotJSON=plotJSON))
            else:
                tree = False
                return render(request,'view_job.html',
                        dict(title=str(job.application) + ' Results',
                        result=finalResult, tree=tree, job=job, plotJSON=plotJSON))
        elif job.application.output_type == 'text/properties.table':
            f = open(job.output, 'r')
            csvinput = csv.reader(f)
            csvOutput = []
            for line in csvinput:
                csvOutput.append(line)
            f.close()
            return render(request,'view_csv.html',
                    dict(title=str(job.application) + ' Results',
                    result=finalResult, job=job, csv=csvOutput))
        elif job.application.output_type == 'text/bins.table':
            f = open(job.output, 'r')
            csvinput = csv.reader(f)
            bins = defaultdict(list)
            iterCsvInput = iter(csvinput)
            next(iterCsvInput)  # skip first line
            for line in iterCsvInput:
                cid = line[0]
                binSize = line[1]
                bin = line[2]
                bins[bin].append(cid)
            f.close()
            bins = OrderedDict(sorted(list(bins.items()), key=lambda t: \
                               int(t[0])))
            return render(request,'bins.html',
                    dict(title=str(job.application) + ' Results',
                    result=finalResult, job=job, bins=bins))
        elif job.application.output_type == 'chemical/x-mdl-sdfile':
            f = open(job.output, 'r')
            sdf = f.read()
            f.close()
            nextStep = job.input
            deleteJob(request.user, job.id)
            if nextStep == 'workbench':
                newJob = createJob(request.user, 'Upload Compounds', '',
                    ['--user=' + str(request.user.id)], sdf)
                time.sleep(2)
                return redirect(tools.views.view_job, job_id=newJob.id,
                    resource='')
            elif nextStep == 'smiles':
                result, error = batch_sdf_to_smiles(sdf) 
                return HttpResponse(result,
                                 content_type=job.application.output_type)
            else:
                return HttpResponse(sdf,
                                 content_type=job.application.output_type)
        else:

            # if mimetype is unknown, just send the file to the user

            return redirect(view_job, job_id=job.id, resource='download'
                            , filename='output')
    elif job.status == Job.RUNNING:
        return render(request,'wait.html',
                                  dict(title=job.application.name
                                  + ' Job Running <img src="/static/images/spinner.gif"/>'
                                  ))
    elif job.status == Job.FAILED:
        messages.error(request,
                       'Job Failed due to invalid input data and/or invalid selected options. Please double check that your uploaded data (compounds and/or numeric data), and input options are valid and try running the tool again.'
                       )
        return render(request,'view_job.html', dict(title='Error: '
                                  + str(job.application) + ' Job Failed'
                                  , job=job))


@guest_allowed
def list_jobs(request):
    matches = getJobList(request.user)
    for match in matches:
        match.start_time = str(match.start_time)
    return render(request,'list_jobs.html', dict(matches=matches) )
