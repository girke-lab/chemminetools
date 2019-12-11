from django.shortcuts import render
from django.http import HttpResponse,JsonResponse
from lockdown.decorators import lockdown
from compounddb.models import Compound, Tag
from django.contrib import messages
from collections import OrderedDict
import sys
import traceback


from .chembl_helpers import (
    AnnotationWithMeshSearch,
    ActivitySearch,
    mapToChembl,
    mapToUniprot,
    compoundNameAutocomplete
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

@lockdown()
def newTS(request):
    # Default local variables
    query_submit = False
    message = None
    annotation_list = None
    annotation_matches = None
    activity_list = None
    activity_matches = None
    allTags = Tag.allUserTagNames(request.user)
    compoundDbs = readSources("unichem")
    proteinDbs= readSources("uniprot")
    defaultCompoundDb = "1"
    defaultProteinDb = "ACC+ID"
    groupingCol = 1
    
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
                    "name":"Query ID",
                    "sql" :"",
                    "id":"query_id"
                    }
            
            myAnnotationSearch = AnnotationWithMeshSearch(id_type)
            annotation_list = myAnnotationSearch.annotation_list
            annotation_matches = myAnnotationSearch.search_grouped(ids)
                
            
            myActivitySearch = ActivitySearch(id_type)
            activity_list = myActivitySearch.activity_list
            
            # Exclude ActivitySearch from search-by-target by default
            if id_type == 'target' and not include_activity:
                activity_matches = None
            else:
                activity_matches = myActivitySearch.search_grouped(ids)

            def addQueryCol(matches) :
                for key, rowGroup in matches.items():
                    queryKey = idMapping[key]
                    for i in range(len(rowGroup)):
                        matches[key][i] = OrderedDict([("query_id",queryKey)] + [ item for item in rowGroup[i].items()] )


            if len(idMapping) != 0:
                groupingCol = 2
                annotation_list.insert(0,queryIdCol)
                activity_list.insert(0,queryIdCol)
                
                addQueryCol(annotation_matches)
                if activity_matches != None:
                    addQueryCol(activity_matches)
                                           


            #extract table name for grouping columns
            for col in annotation_list:
                try:
                    col["table"] = col["sql"][0:col["sql"].index(".")].replace("_"," ")
                except:
                    col["table"] = ""
            for col in activity_list:
                try:
                    col["table"] = col["sql"][0:col["sql"].index(".")].replace("_"," ")
                except:
                    col["table"] = ""



    except Exception as e:
        print("exception in newTS:", sys.exc_info())
        traceback.print_tb(sys.exc_info()[2])
        message = str(e)
    
    context = {
        'query_submit' : query_submit,
        'message' : message,
        'id_type' : id_type,
        'annotation_list' : annotation_list,
        'annotation_matches' : annotation_matches,
        'annotation_child_rows': True,
        'activity_list' : activity_list,
        'activity_matches' : activity_matches,
        'activity_child_rows' : False,
        'tags' : allTags,
        'compoundDbs' : compoundDbs,
        'defaultCompoundDb': defaultCompoundDb,
        'proteinDbs': proteinDbs,
        'defaultProteinDb': defaultProteinDb,
        'groupingCol' : groupingCol
        }
    
    return render(request, 'targetsearch/new_ts.html', context)
def compoundNames(request,query):
    names = compoundNameAutocomplete(query)
    data = [ {'chembl_id': n.chembl_id, 'name': n.synonyms} for n in names]
    return JsonResponse(data,safe=False)
