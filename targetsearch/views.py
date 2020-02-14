from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from lockdown.decorators import lockdown
from compounddb.models import Compound, Tag
from django.contrib import messages
from collections import OrderedDict
from django.template import Context, Engine
import sys
import traceback


from .chembl_helpers import (
    AnnotationSearch,
    ActivitySearch,
    DrugIndicationSearch,
    AnnotationWithDrugIndSearch,
    mapToChembl,
    mapToUniprot,
    compoundNameAutocomplete,
    targetNameAutocomplete,
    )

from django.conf import settings

import os
import csv

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
    molregno_to_chembl = None
    drugind_tables = None
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

            myAnnotationSearch = AnnotationWithDrugIndSearch(id_type, ids)
            annotation_info = myAnnotationSearch.table_info
            annotation_matches = myAnnotationSearch.get_grouped_results()
            molregno_to_chembl = myAnnotationSearch.molregno_to_chembl

            # Generate Drug Indication tables
            drugind_tables = dict()
            for molregno, drugind_obj in myAnnotationSearch.drugind_objs.items():
                drugind_tables[molregno] = tableHtml(drugind_obj, table_class="table table-striped")

            # Exclude ActivitySearch from search-by-target by default
            if id_type == 'target' and not include_activity:
                activity_info = None
                activity_matches = None
            else:
                myActivitySearch = ActivitySearch(id_type, ids)
                activity_info = myActivitySearch.table_info
                activity_matches = myActivitySearch.get_grouped_results()

            def addQueryCol(matches) :
                for key, rowGroup in matches.items():
                    queryKey = idMapping[key]
                    for i in range(len(rowGroup)):
                        matches[key][i] = OrderedDict([("query_id",queryKey)] + [ item for item in rowGroup[i].items()] )


            if len(idMapping) != 0:
                groupingCol += 1
                annotation_info.insert(0,queryIdCol)

                addQueryCol(annotation_matches)
                if activity_matches != None:
                    activity_info.insert(0,queryIdCol)
                    addQueryCol(activity_matches)

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
        'molregno_to_chembl' : molregno_to_chembl,
        'drugind_tables' : drugind_tables,
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
