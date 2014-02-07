#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.shortcuts import get_object_or_404, get_list_or_404, \
    render_to_response, redirect
from django.template import RequestContext
from django.http import HttpResponseRedirect, HttpResponse, Http404
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.core.exceptions import ObjectDoesNotExist
from pubchem_soap_interface.pubchemdl import download
from pubchem_soap_interface.SimilaritySearch import SimilaritySearch
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_sdf, \
    InputError, sdf_to_smiles
from django.views.decorators.csrf import csrf_exempt
from tools.models import Application, Job, ApplicationCategories
from tools.runapp import *
from django.contrib.auth.models import User
from guest.decorators import guest_allowed, login_required
import random, string, time, re
from converters import inputConverters, outputConverters
from bs4 import BeautifulSoup

@csrf_exempt
def listCMTools(request, url):
    if request.method == 'POST':
        toolList = u'Category\tName\tInput\tOutput\n'
        allTools = Application.objects.exclude(category = ApplicationCategories.objects.get(name="Internal"))
        for tool in allTools:
            inputObject = RObjectType(tool, 'input')
            outputObject = RObjectType(tool, 'output') 
            outputObject = re.sub('\n.*', '', outputObject, count=0) 
            toolList = toolList + tool.category.name + "\t" + tool.name + "\t"\
                    + inputObject + "\t" + outputObject + "\n"
        return HttpResponse(toolList,
                            mimetype='text/plain')
    else:
        return HttpResponse('ERROR: query must be an HTTP POST\n',
                            mimetype='text/plain')

def RObjectType(tool, type):
    if type == 'input':
        if tool.input_type in inputConverters.keys():
            object = inputConverters[tool.input_type]
        else:
            object = inputConverters['default']
    else:
        if tool.output_type in outputConverters.keys():
            object = outputConverters[tool.output_type]
        else:
            object = outputConverters['default']
    return re.sub('\n.*', '', object, count=0) 

@csrf_exempt
def toolDetails(request, url):
    if not request.method == 'POST':
        return HttpResponse('ERROR: query must be an HTTP POST\n',
                            mimetype='text/plain')
    tool_name = request.POST['tool_name']
    try:
        app = Application.objects.get(name__iexact=tool_name)
    except ObjectDoesNotExist:
        return HttpResponse('ERROR: tool name not in database.\n Check that the name matches exactly.\n',
                            mimetype='text/plain')
    details = 'Category:\t\t' + app.category.name + '\n' +\
              'Name:\t\t\t' + app.name + '\n' +\
              'Input R Object:\t\t' + RObjectType(app, 'input') + '\n' +\
              'Input mime type:\t' + app.input_type + '\n' +\
              'Output R Object:\t' + RObjectType(app, 'output') + '\n' +\
              'Output mime type:\t' + app.output_type + '\n' +\
              '###### BEGIN DESCRIPTION ######\n' +\
              BeautifulSoup(app.description).get_text() + '\n' +\
              '####### END DESCRIPTION #######\n' 
    optionCounter = 0
    for option in ApplicationOptions.objects.filter(application=app):
        optionCounter += 1
        details = details + "Option " + str(optionCounter) + ": '" +\
                option.name + "'\nAllowed Values: " 
        for value in ApplicationOptionsList.objects.filter(category=option):
            details = details + " '" + value.name + "'"
        details = details + "\n"
    details = details + "Example function call:\n\t" +\
            "job <- launchCMTool(\n\t\t'" +\
                app.name + "',\n\t\t" +\
                RObjectType(app, 'input')
    for option in ApplicationOptions.objects.filter(application=app):
        details = details + ",\n\t\t'" + option.name + "'='" +\
            ApplicationOptionsList.objects.filter(category=option)[0].name + "'"
    details = details + "\n\t)\n"
    return HttpResponse(details,
                        mimetype='text/plain')

@csrf_exempt
def getConverter(request, url):
    if not request.method == 'POST':
        return HttpResponse('ERROR: query must be an HTTP POST\n',
                            mimetype='text/plain')
    converterType = request.POST['converterType']
    toolName = request.POST['toolName']
    if converterType == 'input':
        mimeType = Application.objects.get(name=toolName).input_type
        if mimeType in inputConverters.keys():
            result = inputConverters[mimeType]
        else:
            result = inputConverters['default']
    else:
        mimeType = Application.objects.get(name=toolName).output_type
        if mimeType in outputConverters.keys():
            result = outputConverters[mimeType]
        else:
            result = outputConverters['default']
    return HttpResponse(result, mimetype='text/plain')

@csrf_exempt
def launchCMTool(request, url):
    if not request.method == 'POST':
        return HttpResponse('ERROR: query must be an HTTP POST\n',
                            mimetype='text/plain')
    # get ChemmineR user
    try:
        user = User.objects.get(username='ChemmineR')
    except ObjectDoesNotExist:
        # create user with a random password
        sysrand = random.SystemRandom()
        length = 64
        chars = string.ascii_letters + string.digits + '!@#$%^&*()'
        password = ''.join(sysrand.choice(chars) for _ in xrange(length))
        user = User.objects.create_user('ChemmineR', 'none', password)

    # create and validate job form
    tool_name = request.POST['tool_name']
    try:
        app = Application.objects.get(name__iexact=tool_name)
    except ObjectDoesNotExist:
        return HttpResponse('ERROR: tool name not in database.\n Check that the name matches exactly.\n',
                            mimetype='text/plain')

    # parse form options
    fields = {'application': app.id}
    for optionName in request.POST:
        if optionName in ('input', 'tool_name'):
            continue
        try:
            AppOption = ApplicationOptions.objects.get(name__iexact=optionName, application=app)
            listid = ApplicationOptionsList.objects.get(category=AppOption, name__iexact=request.POST[optionName]).id
            fields[AppOption.name] = listid
        except ObjectDoesNotExist:
            pass

    appForm = getAppForm(app.id, user)
    form = appForm(fields, auto_id=False)
    if not form.is_valid():
        return HttpResponse('ERROR: invalid or missing input options',
                           mimetype='text/plain')

    # launch job
    commandOptions, optionsList = parseToolForm(form)
    newJob = createJob(user, app.name, optionsList, commandOptions, request.POST['input'])

    # return task id token to user
    return HttpResponse(newJob.task_id, mimetype='text/plain')

@csrf_exempt
def jobStatus(request, url):
    if not request.method == 'POST':
        return HttpResponse('ERROR: query must be an HTTP POST\n',
                            mimetype='text/plain')
    time.sleep(2)
    task_id = request.POST['task_id']
    user = User.objects.get(username='ChemmineR')
    try:
        job = Job.objects.get(task_id=task_id, user=user)
    except ObjectDoesNotExist:
        return HttpResponse('ERROR: job not in database.\n', mimetype='text/plain')
    job = updateJob(user, job.id)
    if job.status == Job.RUNNING:
        return HttpResponse('RUNNING', mimetype='text/plain')
    if job.status == Job.FAILED:
        return HttpResponse('FAILED', mimetype='text/plain')
    return HttpResponse('FINISHED', mimetype='text/plain')

@guest_allowed
def showJob(request, task_id):
    cmuser = User.objects.get(username='ChemmineR')
    webuser = request.user
    try:
        job = Job.objects.get(task_id=task_id, user=cmuser)
    except ObjectDoesNotExist:
        raise Http404
    job.user = webuser
    job.save()
    return redirect('tools.views.view_job', job_id=job.id, resource='')

@csrf_exempt
def jobResult(request, url):
    if not request.method == 'POST':
        return HttpResponse('ERROR: query must be an HTTP POST\n',
                            mimetype='text/plain')
    time.sleep(2)
    task_id = request.POST['task_id']
    user = User.objects.get(username='ChemmineR')
    try:
        job = Job.objects.get(task_id=task_id, user=user)
    except ObjectDoesNotExist:
        return HttpResponse('ERROR: job not in database.\n', mimetype='text/plain')
    job = updateJob(user, job.id)
    if job.status == Job.RUNNING:
        return HttpResponse('RUNNING', mimetype='text/plain')
    if job.status == Job.FAILED:
        deleteJob(user, job.id)
        return HttpResponse('FAILED', mimetype='text/plain')
    f = open(job.output, 'r')
    result = f.read()
    f.close()
    # deleteJob(user, job.id)
    return HttpResponse(result, mimetype='text/plain')

# below this line are legacy tools
# to be eliminated in a future version

@csrf_exempt
def runapp(request, url):
    if request.method == 'POST':
        app = str(request.GET.get('app', '1'))
        if app == 'getIds':
            try:

                # parse and test input

                cids = str(request.POST['cids'])

                # launch app

                response = getIds(cids)
            except:
                response = \
                    HttpResponse('ERROR: no results or invalid query',
                                 mimetype='text/plain')
        elif app == 'searchString':
            try:

                # parse and test input

                smiles = str(request.POST['smiles'])

                # launch app

                response = searchString(smiles)
            except:
                response = \
                    HttpResponse('ERROR: no results or invalid query',
                                 mimetype='text/plain')
        elif app == 'sdf2smiles':
            try:

                # parse and test input

                sdf = str(request.POST['sdf'])

                # launch app

                response = sdf2smiles(sdf)
            except:
                response = \
                    HttpResponse('ERROR: no results or invalid query',
                                 mimetype='text/plain')
        elif app == 'smiles2sdf':
            try:

                # parse and test input

                smiles = str(request.POST['smiles'])

                # launch app

                response = smiles2sdf(smiles)
            except:
                response = \
                    HttpResponse('ERROR: no results or invalid query',
                                 mimetype='text/plain')
        else:
            response = HttpResponse('ERROR: invalid app specified',
                                    mimetype='text/plain')
        return response
    else:
        return HttpResponse('ERROR: query must be an HTTP POST\n',
                            mimetype='text/plain')

def getIds(cids):
    cids = cids.split(',')
    cids = [int(cid) for cid in cids]
    response = download(cids)
    return HttpResponse(response, mimetype='text/plain')


def searchString(smiles):
    id_list = SimilaritySearch(smiles)
    id_list = [str(cid) for cid in id_list]
    id_list = ','.join(id_list)
    return HttpResponse(id_list, mimetype='text/plain')


def smiles2sdf(smiles):
    sdf = smiles_to_sdf(smiles)
    return HttpResponse(sdf, mimetype='text/plain')


def sdf2smiles(sdf):
    sdf = sdf_to_smiles(sdf)
    return HttpResponse(sdf, mimetype='text/plain')
