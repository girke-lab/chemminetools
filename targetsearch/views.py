from django.shortcuts import render

from .chembl_helpers import (
    chemblTargetAccessions,
    accessionToChembl,
    chemblTargetAccessionsByAnnotations,
    mapToChembl
    )

from django.conf import settings
import os

def home(request):
    return render(request, 'targetsearch/home.html')

def getTargets(request):
    sourceDbs  = readSources()
    query_submit = False
    annotationMatches = None
    matches = None
    chembl_id_tuple = None
    message = None

    if 'chembl_id' in request.GET:
        source_id = request.GET['source_id']
        chembl_id_str = request.GET['chembl_id']
        chembl_id_tuple = tuple(chembl_id_str.split())
        print("source_id: "+str(source_id))
        if(source_id != "1"): # not CHEMBL
            print("mapping to source_id")
            chembl_id_tuple = mapToChembl(chembl_id_tuple,source_id)

        if len(chembl_id_tuple) > 0 :
            query_submit = True
            annotationMatches = chemblTargetAccessionsByAnnotations(chembl_id_tuple)
            matches = chemblTargetAccessions(chembl_id_tuple)
        else:
            message = "None of the given IDs could be converted to ChEMBL IDs"
    
    context = {
        'query_submit': query_submit,
        'chembl_id_tuple': chembl_id_tuple,
        'matches': matches,
        'annotationMatches':annotationMatches,
        'sources':sourceDbs,
        'message': message,
        }

    
    return render(request, 'targetsearch/getTargets.html', context)

def getChembl(request):
    if 'accession_id' in request.GET:
        query_submit = True
        accession_id_str = request.GET['accession_id']
        accession_id_tuple = tuple(accession_id_str.split())
        matches = accessionToChembl(accession_id_tuple)
    else:
        query_submit = False
        accession_id_tuple = None
        matches = None
    
    context = {
        'query_submit': query_submit,
        'accession_id_tuple': accession_id_tuple,
        'matches': matches,
        }
    
    return render(request, 'targetsearch/getChembl.html', context)

def readSources():
    sources = {}
    filename = os.path.join(settings.PROJECT_DIR,"unichem_sources.txt")
    with open(filename) as f:
        for line in f:
            (key, val) = line.split("\t")
            sources[int(key)] = val
    return sources




