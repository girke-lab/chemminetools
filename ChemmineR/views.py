#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.shortcuts import get_object_or_404, get_list_or_404, \
    render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect, HttpResponse
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.core.exceptions import ObjectDoesNotExist
from pubchem_soap_interface.pubchemdl import download
from pubchem_soap_interface.SimilaritySearch import SimilaritySearch
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_sdf, \
    InputError, sdf_to_smiles
from django.views.decorators.csrf import csrf_exempt
from tools.models import Application
from tools.runapp import *
from django.contrib.auth.models import User
import random, string

@csrf_exempt
def listCMTools(request, url):
    if request.method == 'POST':
        toolList = u'Category\tName\tInput\tOutput\n'
        allTools = Application.objects.all()
        for tool in allTools:
            toolList = toolList + tool.category.name + "\t" + tool.name + "\t"\
                    + tool.input_type + "\t" + tool.output_type + "\n"
        return HttpResponse(toolList,
                            mimetype='text/plain')
    else:
        return HttpResponse('ERROR: query must be an HTTP POST\n',
                            mimetype='text/plain')

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
        password = ''.join(foo.choice(chars) for _ in xrange(length))
        user = User.objects.create_user('ChemmineR', 'none', password)

    # create and validate job form
    tool_name = request.POST['tool_name']
    try:
        app = Application.objects.get(name=tool_name)
    except ObjectDoesNotExist:
        return HttpResponse(tool_name, mimetype='text/plain')
        return HttpResponse('ERROR: tool name not in database.\n Check that the name matches exactly (case sensitive).\n',
                            mimetype='text/plain')
    # appForm = getAppForm(app.id, user)
    # form = appForm(request.POST)
    # if not form.is_valid():
    #    return HttpResponse('ERROR: Invalid input options. Check that the names match exactly (case sensitive).\n',
    #                        mimetype='text/plain')

    # launch job
    # commandOptions, optionsList = parseToolForm(form)
    newJob = createJob(user, app.name, '', '', request.POST['input'])

    # return task id token to user
    return HttpResponse(newJob.task_id, mimetype='text/plain')

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


