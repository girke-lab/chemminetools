from django.shortcuts import render, redirect
from django.http import HttpResponse, JsonResponse, Http404
from django.urls import reverse
from lockdown.decorators import lockdown
from compounddb.models import Compound, Tag
from django.contrib import messages
from collections import OrderedDict
from django.template import Context, Engine
import sys
import traceback


from .helpers import *

from django.conf import settings
from django.views.decorators.cache import cache_page

import os
import csv
import json

def readSources(type):
    sources = {}
    filename = os.path.join(settings.PROJECT_DIR,type+"_sources.txt")
    with open(filename) as f:
        for line in f:
            (key, val) = line.split("\t")
            sources[key.strip()] = val.strip()
    return sources

def newTS(request):
    # Default local variables
    query_submit = False
    message = None
    annotation_info = None
    annotation_matches = None
    drugind_json = None
    activity_info = None
    activity_matches = None
    allTags = Tag.allUserTagNames(request.user)
    compoundDbs = readSources("unichem")
    proteinDbs= readSources("uniprot")
    defaultCompoundDb = "1"
    defaultProteinDb = "ACC+ID"
    groupingCol = 0

    # Default GET request variables
    id_type = 'compound'
    ids = list()
    include_activity = False
    source_id = 1

    # Retrieve GET request variables
    if 'id_type' in request.GET:
        id_type = request.GET['id_type']
    if 'ids' in request.GET:
        if id_type == 'homolog-target':
            ids = request.GET.getlist('ids')
        else:
            ids = list(request.GET['ids'].split())
    if 'include_activity' in request.GET:
        include_activity = True
    if 'tags' in request.GET:
        for c in [compound.cid for compound in Compound.byTagNames(request.GET.getlist("tags"),request.user)]:
            ids.append(c)
    if 'source_id' in request.GET:
        source_id = request.GET['source_id']

    # Generate content
    try:
        idMapping={}
        if id_type == 'compound' and source_id != '1':
            idMapping= mapToChembl(ids, source_id)
            ids = list(idMapping.keys())
        elif id_type == 'target' and source_id != 'ACC':
            idMapping= mapToUniprot(ids, source_id)
            ids = list(idMapping.keys())
        elif id_type == 'homolog-target':
            #full_ids = ids # Full context for homolog handling
            #ids = [ i.split(',')[2] for i in full_ids ] # Bare Accession IDs for old target-search

            # Context dictionary for homolog handling
            # ex:
            # homolog_context = {
            #     'P29274': {'paralog': 'P30542'},
            #     'P29275': {'paralog': 'P30542'},
            #     'P0DMS8': {'paralog': 'P30542'},
            #     }
            homolog_context = dict()
            for i in ids:
                [relation, src_id, homolog_id] = i.split(',')
                if homolog_id not in homolog_context.keys():
                    homolog_context[homolog_id] = dict()
                if src_id == homolog_id:
                    continue # Don't bother with "X is a homolog of X"
                if relation not in homolog_context[homolog_id].keys():
                    homolog_context[homolog_id][relation] = set()
                homolog_context[homolog_id][relation].add(src_id)

            ids = list(homolog_context.keys()) # Bare Accession IDs for old target-search

            # Prepare homolog relation descriptions
            homolog_desc = dict()
            for homolog_id, relations in homolog_context.items():
                desc_parts = list()
                for relation, src_ids in sorted(relations.items()):
                    desc_parts.append("{} of {}".format(relation, ', '.join(sorted(list(src_ids)))))
                homolog_desc[homolog_id] = '; '.join(desc_parts)

        if len(ids) != 0:
            query_submit = True
            queryIdCol = {
                    "id": "query_id",
                    "sql": None,
                    "table": "Query ID",
                    "name": "Query ID",
                    "desc": "Original compound ID prior to ChEMBL conversion",
                    "visible": True,
                    }
            headerTextCol = {
                "id": "header_text",
                "sql": None,
                "table": "Header Text",
                "name": "Header Text",
                "desc": "Description text to show in row-group headers (i.e. source query IDs, translations, etc.)",
                "visible": False,
                }

            myAnnotationSearch = AnnotationWithDrugIndSearch(id_type, ids)
            annotation_info = myAnnotationSearch.table_info
            annotation_matches = myAnnotationSearch.get_grouped_results()

            # Generate Drug Indication JSON data
            drugind_json = dict()
            for chembl_id, drugind_obj in myAnnotationSearch.drugind_objs.items():
                colnames = [ {'title': c['name']} for c in drugind_obj.table_info ]

                data = list()
                for row in drugind_obj.get_results():
                    data.append(list(row.values()))

                drugind_json[chembl_id] = dict()
                drugind_json[chembl_id]['colnames'] = colnames
                drugind_json[chembl_id]['data'] = data

            # Exclude ActivitySearch from search-by-target by default
            if id_type in ['target', 'homolog-target'] and not include_activity:
                activity_info = None
                activity_matches = None
            else:
                myActivitySearch = ActivitySearch(id_type, ids)
                activity_info = myActivitySearch.table_info
                activity_matches = myActivitySearch.get_grouped_results()

            # Add column with original query info if necessary
            def addQueryCol(matches) :
                for key, rowGroup in matches.items():
                    queryKey = idMapping[key]
                    for i in range(len(rowGroup)):
                        matches[key][i] = OrderedDict([("query_id",queryKey)] + [ item for item in rowGroup[i].items()] )

            if len(idMapping) != 0:
                #groupingCol += 1
                annotation_info.insert(0,queryIdCol)

                addQueryCol(annotation_matches)
                if activity_matches != None:
                    activity_info.insert(0,queryIdCol)
                    addQueryCol(activity_matches)

            # Add row-group header column
            def addHeaderCol(matches):
                for key, rowGroup in matches.items():
                    if id_type in ['compound', 'target'] and len(idMapping) != 0:
                        header_text = "{} ({})".format(key, idMapping[key])
                    elif id_type == 'homolog-target':
                        if homolog_desc[key]:
                            header_text = "{} ({})".format(key, homolog_desc[key])
                        else:
                            header_text = key
                    else:
                        header_text = key
                    for i in range(len(rowGroup)):
                        matches[key][i] = OrderedDict([('header_text', header_text)] + [item for item in rowGroup[i].items()])

            #annotation_info.append(headerTextCol)
            annotation_info.insert(0, headerTextCol)
            addHeaderCol(annotation_matches)
            if activity_matches != None:
                #activity_info.append(headerTextCol)
                activity_info.insert(0, headerTextCol)
                addHeaderCol(activity_matches)

    except Exception as e:
        print("exception in newTS:", sys.exc_info())
        traceback.print_tb(sys.exc_info()[2])
        message = str(e)

    context = {
        'query_submit' : query_submit,
        'message' : message,
        'id_type' : id_type,
        'annotation_info' : annotation_info,
        'annotation_matches' : annotation_matches,
        'drugind_json' : json.dumps(drugind_json),
        'activity_info' : activity_info,
        'activity_matches' : activity_matches,
        'tags' : allTags,
        'compoundDbs' : compoundDbs,
        'defaultCompoundDb': defaultCompoundDb,
        'proteinDbs': proteinDbs,
        'defaultProteinDb': defaultProteinDb,
        'groupingCol' : groupingCol
        }

    return render(request, 'targetsearch/new_ts.html', context)

def tableHtml(search_obj, header=True, footer=False, table_class=""):
    """Takes a SearchBase object and returns an HTML table"""

    error = None
    table_info = None
    table_data = None

    engine = Engine.get_default()
    template = engine.get_template('targetsearch/table.html')

    try:
        table_info = search_obj.table_info
        table_data = search_obj.get_results()
    except Exception as e:
        error = str(e)

    context = {
        'error' : error,
        'header' : header,
        'footer' : footer,
        'table_info' : table_info,
        'table_data' : table_data,
        'table_class' : table_class,
    }

    return template.render(Context(context))

def drugIndTable(request):
    if 'molregno' in request.GET:
        molregno = request.GET['molregno']
    else:
        error = "drugIndTable: Missing molregno in GET variables"
        return render(request, 'targetsearch/table.html', {'error':error})

    myDrugIndSearch = DrugIndicationSearch(molregno)
    html = tableHtml(myDrugIndSearch)

    return HttpResponse(html)

def compoundNames(request,query):
    names = compoundNameAutocomplete(query)
    data = [ {'chembl_id': n.chembl_id, 'name': n.synonyms} for n in names]
    return JsonResponse(data,safe=False)

def targetNames(request, query):
    names = targetNameAutocomplete(query)
    data = [ {'accession_id': n.accession,
              'name': (n.description+' ('+n.organism+')')} for n in names ]
    return JsonResponse(data, safe=False)

@cache_page(60 * 120)
def chemblPNG(request, chembl_id):
    try:
        img = getChemblPNG(chembl_id, mwt_limit=2000, shrink=False)
        return HttpResponse(img, content_type='image/png')
    except Exception as e:
        raise Http404(str(e))

@cache_page(60 * 120)
def chemblSVG(request, chembl_id):
    try:
        img = getChemblSVG(chembl_id)
        return HttpResponse(img, content_type='image/svg+xml')
    except Exception as e:
        raise Http404(str(e))

def ajax(request, action):
    def die(msg):
        ajaxResponse = { 'success' : False, 'message' : msg }
        return JsonResponse(ajaxResponse)

    if action == 'homologs':
        try:
            ids = request.POST['ids']
            homolog_type = request.POST['homolog_type']
            result = homologsAjax(ids.split(), homolog_type)
            return JsonResponse(result)
        except Exception as e:
            return die(str(e))
    else:
        return die('Unknown action: {}'.format(action))
