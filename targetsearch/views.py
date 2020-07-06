from django.shortcuts import render
from django.http import HttpResponse, JsonResponse, Http404
from lockdown.decorators import lockdown
from compounddb.models import Compound, Tag
from tools.models import *
from django.contrib import messages
from collections import OrderedDict
from django.template import Context, Engine
import sys
import traceback


from .helpers import (
    AnnotationSearch,
    ActivitySearch,
    DrugIndicationSearch,
    AnnotationWithDrugIndSearch,
    mapToChembl,
    mapToUniprot,
    compoundNameAutocomplete,
    targetNameAutocomplete,
    getChemblPNG,
    getChemblSVG,
    )

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

def detailPage(request,id):
    
    print("rendering detail page for "+str(id))

    myAnnotationSearch = AnnotationWithDrugIndSearch("compound", [id])
    annotation_info = myAnnotationSearch.table_info
    annotation_matches = myAnnotationSearch.get_grouped_results()

    myActivitySearch = ActivitySearch("compound", [id])
    activity_info = myActivitySearch.table_info
    activity_matches = myActivitySearch.get_grouped_results()

    drugind_json = drugIndicationData(myAnnotationSearch.drugind_objs);

    return render(request,'detail.html', {
        
        'annotation_info' : annotation_info,
        'annotation_matches' : annotation_matches,
        'activity_info' : activity_info,
        'activity_matches' : activity_matches,
        'drugind_json' : json.dumps(drugind_json),
        })

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
    similarityJobs = [(job.id,str(job)) for job in Job.objects.filter(user=request.user,application__category_id=5)]

    # Default GET request variables
    id_type = 'compound'
    ids = list()
    include_activity = False
    source_id = 1
    similarity_job_id=-1

    # Retrieve GET request variables
    if 'id_type' in request.GET:
        id_type = request.GET['id_type']
    if 'ids' in request.GET:
        ids = list(request.GET['ids'].split())
    if 'include_activity' in request.GET:
        include_activity = True
    if 'tags' in request.GET:
        for c in [compound.cid for compound in Compound.byTagNames(request.GET.getlist("tags"),request.user)]:
            ids.append(c)
    if 'source_id' in request.GET:
        source_id = request.GET['source_id']
    if 'similarity_job_id' in request.GET:
        similarity_job_id = request.GET['similarity_job_id']

    # Generate content
    try:
        # Only attempt conversion for compound search. Skip if already ChEMBL.
        idMapping={}
        if id_type == 'compound' and source_id != '1':
            idMapping= mapToChembl(ids, source_id)
            ids = list(idMapping.keys())
        elif id_type == 'target' and source_id != 'ACC':
            idMapping= mapToUniprot(ids, source_id)
            ids = list(idMapping.keys())

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

            # Exclude ActivitySearch from search-by-target by default
            if id_type == 'target' and not include_activity:
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
        'similarityJobs': similarityJobs
        }

    return render(request, 'targetsearch/new_ts.html', context)

def addMappedQueryColumn(mapping,columnDefinition,annotation_info, annotation_matches,activity_info,activity_matches):

    def addQueryCol(matches) :
        for key, rowGroup in matches.items():
            if key in mapping:
                queryKey = mapping[key]
                for i in range(len(rowGroup)):
                    matches[key][i] = OrderedDict([(columnDefinition["id"],queryKey)] + [ item for item in rowGroup[i].items()] )


    annotation_info.insert(0,columnDefinition)

    addQueryCol(annotation_matches)
    if activity_matches != None:
        activity_info.insert(0,columnDefinition)
        addQueryCol(activity_matches)


def readSimilarityMappingData(user,job_id):
    job = Job.objects.get(id=job_id, user=user)
    mapping = {}

    try:
        f = open(job.output, 'r')
        csvinput = csv.reader(f, delimiter=' ')
        for line in csvinput:
            #map from target back to query
            # multiple queries may map to the same target
            queries = mapping.setdefault(line[1],[]);
            queries.append(line[0]);
    finally:
        f.close()

    for key, value in mapping.items():
        mapping[key] = ",".join(value)

    return mapping


## map1: x->y and map2: y -> z, return x->z    
def composeMaps(map1,map2):
    composedMap = {}
    for key1,value1 in map1.items():
        if value1 in map2:
            composedMap[key1] = map2[value1]
    return composedMap


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

def drugIndicationData(drugind_objs):
    # Generate Drug Indication JSON data
    drugind_json = dict()
    for chembl_id, drugind_obj in drugind_objs.items():
        colnames = [ {'title': c['name']} for c in drugind_obj.table_info ]

        data = list()
        for row in drugind_obj.get_results():
            data.append(list(row.values()))

        drugind_json[chembl_id] = dict()
        drugind_json[chembl_id]['colnames'] = colnames
        drugind_json[chembl_id]['data'] = data
    return drugind_json


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
