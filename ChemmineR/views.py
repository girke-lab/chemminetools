#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.shortcuts import get_object_or_404, get_list_or_404, \
    render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect, HttpResponse
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from pubchem_soap_interface.pubchemdl import download
from pubchem_soap_interface.SimilaritySearch import SimilaritySearch
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_sdf, \
    InputError, sdf_to_smiles
from django.views.decorators.csrf import csrf_exempt
from tools.models import Application
from tools.runapp import createJob

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

    # create and validate job form

    # launch job

    # return task id token to user


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


