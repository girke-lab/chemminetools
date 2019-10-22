from django.shortcuts import render
from django.http import HttpResponse
from lockdown.decorators import lockdown
from compounddb.models import Compound, Tag
from django.contrib import messages

from .chembl_helpers import (
    byActivity,
    byAnnotations,
    AnnotationWithMeshSearch,
    ActivitySearch,
    mapToChembl
    )

from django.conf import settings

import os
import csv

@lockdown()
def search(request):
    sourceDbs  = readSources()
    query_submit = False
    annotationMatches = None
    activityMatches = None
    message = None
    ids = ()
    search_type="compound"
    allTags = Tag.allUserTagNames(request.user)

    if 'id' in request.GET:
        search_type = request.GET['search_type']


        ids = tuple(request.GET['id'].split())

        print("ids: "+str(ids))
        if search_type == 'compound':
            print("doing compound search")
            source_id = request.GET['source_id']
            
            if 'tags' in request.GET:
                ids =ids + tuple( [compound.cid for compound in  Compound.byTagNames(request.GET.getlist("tags"),request.user)])

            try:
                if source_id != "1" : # not CHEMBL
                    print("mapping to source_id")
                    ids = mapToChembl(ids,source_id)


                if len(ids) > 0 :
                    query_submit = True
                    annotationMatches = byAnnotations(chemblIds=ids)
                    activityMatches = byActivity(chemblIds=ids)
                else:
                    message = "None of the given IDs could be converted to ChEMBL IDs"
            except Exception as e:
                print("failed to do compound search: "+str(e))
                message = e

        elif search_type == 'target':
            print("doing target search")
            query_submit = True
            annotationMatches = byAnnotations(accessionIds=ids)
            if 'include_activity' in request.GET:
                activityMatches = byActivity(accessionIds=ids)

    context = {
        'query_submit': query_submit,
        'ids': ids,
        'activityMatches': activityMatches,
        'annotationMatches':annotationMatches,
        'sources':sourceDbs,
        'message': message,
        'search_type':search_type,
        'tags':allTags,
        }

    #print("results: "+str(context))
    
    return render(request, 'targetsearch.html', context)


def readSources():
    sources = {}
    filename = os.path.join(settings.PROJECT_DIR,"unichem_sources.txt")
    with open(filename) as f:
        for line in f:
            (key, val) = line.split("\t")
            sources[int(key)] = val
    return sources

def bs4test(request):
    if 'id_type' in request.GET:
        id_type = request.GET['id_type']
    if 'ids' in request.GET:
        ids = list(request.GET['ids'].split())
    
    #annotation_column_defaults = list()
    #l = AnnotationWithMeshSearch.annotation_list
    #for i in range(0, len(l)):
    #    if l[i]['visible']:
    #        annotation_column_defaults.append(i)
    
    #activity_column_defaults = list()
    #l = ActivitySearch.activity_list
    #for i in range(0, len(l)):
    #    if l[i]['visible']:
    #        activity_column_defaults.append(i)
    
    context = {
        'annotation_list' : AnnotationWithMeshSearch.annotation_list,
        #'annotation_column_defaults' : str(annotation_column_defaults),
        'annotation_matches': AnnotationWithMeshSearch.search(id_type, ids),
        'activity_list' : ActivitySearch.activity_list,
        #'activity_column_defaults': str(activity_column_defaults),
        'activity_matches': ActivitySearch.search(id_type, ids),
        }
    
    return render(request, 'targetsearch/bs4test.html', context)


def newTS(request):
    # Default local variables
    query_submit = False
    message = None
    annotation_list = None
    annotation_matches = None
    activity_list = None
    activity_matches = None
    
    # Default GET request variables
    id_type = 'compound'
    ids = list()
    include_activity = False
    
    # Retrieve GET request variables
    if 'id_type' in request.GET:
        id_type = request.GET['id_type']
    if 'ids' in request.GET:
        ids = list(request.GET['ids'].split())
    if 'include_activity' in request.GET:
        include_activity = True
    
    # Generate content
    if id_type != None and len(ids) != 0:
        query_submit = True
        
        try:
            annotation_list = AnnotationWithMeshSearch.annotation_list
            annotation_matches = AnnotationWithMeshSearch.search(id_type, ids)
            activity_list = ActivitySearch.activity_list
            
            # Exclude ActivitySearch from search-by-target by default
            if id_type == 'target' and not include_activity:
                activity_matches = None
            else:
                activity_matches = ActivitySearch.search(id_type, ids)
        except Exception as e:
            message = str(e)
    
    context = {
        'query_submit' : query_submit,
        'message' : message,
        'id_type' : id_type,
        'annotation_list' : annotation_list,
        'annotation_matches': annotation_matches,
        'activity_list' : activity_list,
        'activity_matches': activity_matches,
        }
    
    return render(request, 'targetsearch/new_ts.html', context)
    
