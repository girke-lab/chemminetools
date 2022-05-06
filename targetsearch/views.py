# views.py - Put anything that returns an HTTP response here
# Non-trivial view functions should be split into a "debug" part that can be
# easily evaluated/tested in a shell, and a rendering stub here.

from django.shortcuts import render, redirect
from django.http import HttpResponse, JsonResponse, Http404
from django.urls import reverse
from lockdown.decorators import lockdown
from compounddb.models import Compound, Tag
from tools.models import *
from django.contrib import messages
#from collections import OrderedDict
#from django.template import Context, Engine
from guest.decorators import guest_allowed
import sys
import traceback

#from .helpers import *
from .views_debug import *

#from django.conf import settings
from django.views.decorators.cache import cache_page

#import os
#import csv
#import json

def detailPage(request,id):
    targetInfo=None
    compoundInfo=None
    id_type = None
    compoundDbs = readSources("unichem")
    proteinDbs= readSources("uniprot")
    activity_info= None
    activity_matches= None


    #print("rendering detail page for "+str(id))

    if 'id_type' in request.GET:
        id_type = request.GET['id_type']
    else:
        id_type = getEntityType(id)

    if id_type == None or id_type == "unknown":
        return render(request,'detail-select-database.html',{
            'given_id': id,
            'compoundDbs' : compoundDbs,
            'proteinDbs':proteinDbs,
        })


    try:
        if id_type == "compound" and 'compound_source_id' in request.GET:
            source_id = request.GET['compound_source_id']
            #get mapping dict, fetch keys, then get either first element, or None
            id=next(iter(mapToChembl([id], source_id).keys()), None)
            #print("new compound id: "+str(id)+", source_id: "+str(source_id))
        elif id_type == "target" and 'target_source_id' in request.GET:
            source_id = request.GET['target_source_id']
            #get mapping dict, fetch keys, then get either first element, or None
            id=next(iter(mapToUniprot([id], source_id).keys()), None)
            #print("new target id: "+str(id)+", source_id: "+str(source_id))

        if id == None:
            raise Exception("ID not found in selected database")
    except Exception as e:
        message = "Sorry, I could not convert that ID to a ChEMBL compound or target. Error was: "+str(e)
        return render(request,'detail-select-database.html',{
            'given_id': id,
            'compoundDbs' : compoundDbs,
            'proteinDbs':proteinDbs,
            'message': message,
        })


    myAnnotationSearch = AnnotationWithDrugIndSearch(id_type, [id])
    annotation_info = myAnnotationSearch.table_info
    annotation_matches = myAnnotationSearch.get_grouped_results()


    drugind_json = drugIndicationData(myAnnotationSearch.drugind_objs);

    if id_type == "target":
        targetInfo = targetSummaryInfo(id)

    if id_type == "compound":
        compoundInfo = compoundSummaryInfo(id);
        # activity search is often very slow, so disable for now
        #myActivitySearch = ActivitySearch(id_type, [id])
        #activity_info = myActivitySearch.table_info
        #activity_matches = myActivitySearch.get_grouped_results()


    return render(request,'detail.html', {
        'given_id':id,
        'type': id_type,
        'annotation_info' : annotation_info,
        'annotation_matches' : annotation_matches,
        'activity_info' : activity_info,
        'activity_matches' : activity_matches,
        'drugind_json' : json.dumps(drugind_json),
        'target_info': targetInfo,
        'compound_info': compoundInfo,
        })

@guest_allowed
def newTS(request,initial_ids=""):
    # Default local variables
    query_submit = False
    message = None
    annotation_info = None
    annotation_matches = None
    drugind_json = None
    activity_info = None
    activity_matches = None
    homolog_type_value= None
    allTags = Tag.allUserTagNames(request.user)
    compoundDbs = readSources("unichem")
    proteinDbs= readSources("uniprot")
    defaultCompoundDb = "1"
    defaultProteinDb = "ACC+ID"
    groupingCol = 0
    similarityJobs = [(job.id,str(job)) for job in Job.objects.filter(user=request.user,application__category_id=5)]
    ids_resubmit_set = set()

    # Default GET request variables
    id_type = 'compound'
    ids = list()
    include_activity = False
    source_id = 1
    similarity_job_id=-1

    # Retrieve GET request variables
    if 'id_type' in request.GET:
        id_type = request.GET['id_type']
    #see if we can identify the entity type of any initial_ids
    elif initial_ids != "": 
        ii_array = initial_ids.split(",") 
        if len(ii_array) > 0:
            type_guess = getEntityType(ii_array[0])
            if type_guess != "unknown":
                id_type = type_guess


    if 'ids' in request.GET:
        if id_type == 'homolog-target':
            ids = request.GET.getlist('ids')
        else:
            ids = list(request.GET['ids'].split())
        initial_ids = " ".join(ids)
    if 'include_activity' in request.GET:
        include_activity = True
    if 'tags' in request.GET:
        for c in [compound.cid for compound in Compound.byTagNames(request.GET.getlist("tags"),request.user)]:
            ids.append(c)
    if 'source_id' in request.GET:
        source_id = request.GET['source_id']
    if 'similarity_job_id' in request.GET:
        similarity_job_id = int(request.GET['similarity_job_id'])

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
            originalQueryCol = {
                    "id": "original_query_id",
                    "sql": None,
                    "table": "Original Query ID",
                    "name": "Original Query ID",
                    "desc": "The compound that the current query compound was originally derived from, based on similarity",
                    "visible": True,
                    }

            myAnnotationSearch = AnnotationWithDrugIndSearch(id_type, ids)
            annotation_info = myAnnotationSearch.table_info
            annotation_matches = myAnnotationSearch.get_grouped_results()

            drugind_json = drugIndicationData(myAnnotationSearch.drugind_objs);

            for v in annotation_matches.values():
                for row in v:
                    ids_resubmit_set.add(row['annotation__chembl_id_lookup__chembl_id'])

            # Exclude ActivitySearch from search-by-target by default
            if id_type in ['target', 'homolog-target'] and not include_activity:
                activity_info = None
                activity_matches = None
            else:
                myActivitySearch = ActivitySearch(id_type, ids)
                activity_info = myActivitySearch.table_info
                activity_matches = myActivitySearch.get_grouped_results()

            if len(idMapping) != 0:
                groupingCol += 1
                addMappedQueryColumn(idMapping,queryIdCol,
                        annotation_info,annotation_matches,
                        activity_info,activity_matches)

            if similarity_job_id != -1:
                similarityMapping = readSimilarityMappingData(request.user,similarity_job_id)
                if len(idMapping) != 0: # need to compose our mapping with previous mapping
                    similarityMapping = composeMaps(idMapping,similarityMapping)
                #print("similarity mapping: \n"+str(similarityMapping))
                if len(similarityMapping) != 0:
                    groupingCol += 1
                    addMappedQueryColumn(similarityMapping, originalQueryCol,
                            annotation_info,annotation_matches,
                            activity_info,activity_matches)

            #if similarity_job_id != -1, then read job_<similarity_job_id>
            # map chembl id (2nd column) back to the original compound id (first column)
            # insert new column to show original compound id.

        if ts_paralog_cache():
            homolog_type_value = 'paralog-cache'
        else:
            homolog_type_value = 'paralog'

        initial_ids = initial_ids.replace(","," ")

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
        'groupingCol' : groupingCol,
        'similarityJobs': similarityJobs,
        'homolog_type_value': homolog_type_value,
        'initial_ids': initial_ids,
        'ids_resubmit' : '+'.join(ids_resubmit_set),
        'chembl_version': getChemblVersion(),
        }

    return render(request, 'targetsearch/new_ts.html', context)

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

def extAnnoByChembl(request, chembl_id, db=None):
    check_avail = checkExtAnno(chembl_id)
    if db is None:
        data = None
    else:
        data = getExtAnno(chembl_id, db)

    context = {
        'check_avail': check_avail,
        'data': data,
        'chembl_id': chembl_id,
        'db': db,
    }

    return render(request, 'targetsearch/extannobychembl.html', context)

def tsFilter(request):
    if 'ids' in request.GET:
        ids = list(request.GET['ids'].split())
    #else:
    #    ids = None

    context = tsFilter_debug(ids)

    return render(request, 'targetsearch/filter.html', context)

def tsFilter2(request):
    if 'ids' in request.GET:
        ids = list(request.GET['ids'].split())
    #else:
    #    ids = None

    context = tsFilter2_debug(ids)

    return render(request, 'targetsearch/filter2.html', context)

def tsAnnoFilter1(request):
    context = {
        'chembl_version': getChemblVersion(),
        'search_type': 'chembl',
    }

    if 'ids' not in request.GET:
        return render(request, 'targetsearch/annofilter1_form.html', context)

    ids = list(request.GET['ids'].split())

    context.update(tsAnnoFilter1_debug(ids))

    if len(context['anno_matches']) == 0:
        return render(request, 'targetsearch/annofilter1_form.html', context)

    return render(request, 'targetsearch/annofilter1.html', context)

def tsAnnoFilter2(request):
    context = {
        'chembl_version': getChemblVersion(),
        'search_type': 'go',
    }

    if 'ids' not in request.GET:
        return render(request, 'targetsearch/annofilter1_form.html', context)

    ids = list(request.GET['ids'].split())

    context.update(tsAnnoFilter2_debug(ids))

    if len(context['anno_matches']) == 0:
        return render(request, 'targetsearch/annofilter1_form.html', context)

    return render(request, 'targetsearch/annofilter1.html', context)
