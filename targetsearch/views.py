from django.shortcuts import render
from django.http import HttpResponse
from lockdown.decorators import lockdown

from .chembl_helpers import (
    byActivity,
    byAnnotations,
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
    ids = None
    search_type=None

    if 'id' in request.GET:
        search_type = request.GET['search_type']
        ids = tuple(request.GET['id'].split())
        print("ids: "+str(ids))
        if search_type == 'compound':
            source_id = request.GET['source_id']
            
            if source_id != "1" : # not CHEMBL
                print("mapping to source_id")
                ids = mapToChembl(ids,source_id)

            if len(ids) > 0 :
                query_submit = True
                annotationMatches = byAnnotations(chemblIds=ids)
                activityMatches = byActivity(chemblIds=ids)
            else:
                message = "None of the given IDs could be converted to ChEMBL IDs"

        elif search_type == 'target':
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




