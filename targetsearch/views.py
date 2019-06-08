from django.shortcuts import render
from django.http import HttpResponse

from .chembl_helpers import (
    chemblTargetAccessionsByActivity,
    accessionToChembl,
    chemblTargetAccessionsByAnnotations,
    mapToChembl
    )

from django.conf import settings

import os
import csv

def home(request):
    return render(request, 'targetsearch/home.html')

def getTargets(request):
    sourceDbs  = readSources()
    query_submit = False
    annotationMatches = None
    activityMatches = None
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
            activityMatches = chemblTargetAccessionsByActivity(chembl_id_tuple)
        else:
            message = "None of the given IDs could be converted to ChEMBL IDs"
    
    context = {
        'query_submit': query_submit,
        'chembl_id_tuple': chembl_id_tuple,
        'activityMatches': activityMatches,
        'annotationMatches':annotationMatches,
        'sources':sourceDbs,
        'message': message,
        }

    
    return render(request, 'targetsearch/getTargets.html', context)

def getTargetsTSV(request):
    # Set download to False to open directly in browser
    download = True
    
    if 'chembl_id' in request.GET:
        chembl_id = request.GET['chembl_id']
    else:
        message = "Error: Expected 'chembl_id' value."
        return HttpResponse(message, content_type='text/plain')
    
    if 'by' in request.GET:
        by = request.GET['by']
    else:
        message = "Error: Expected 'by' value."
        return HttpResponse(message, content_type='text/plain')
    
    if 'download' in request.GET:
        if request.GET['download'] in ['true', 'True']:
            download = True
        else:
            download = False    
    
    if by == 'annotation':
        matches = chemblTargetAccessionsByAnnotations(tuple([chembl_id]))
        
        fieldnames = ['accession', 'mechanism_of_action', 'tid', 'component_id', 'description', 'organism', 'mesh_id', 'mesh_heading']
        header = {'accession': 'Accession ID',
                  'mechanism_of_action': 'Mechanism of Action',
                  'tid': 'Target ID',
                  'component_id': 'Component ID',
                  'description': 'Description',
                  'organism': 'Organism',
                  'mesh_id': 'Mesh ID',
                  'mesh_heading': 'Mesh Heading'}
        
        response = HttpResponse(content_type='text/plain')
        if download:
            response['Content-Disposition'] = "attachment; filename=\"{chembl_id}_byAnnotation.tsv\"".format(chembl_id=chembl_id)
        writer = csv.DictWriter(response, fieldnames=fieldnames, dialect='excel-tab', extrasaction='ignore')
        
        writer.writerow(header)
        if len(matches) != 0:
            for m in matches[chembl_id]:
                writer.writerow(m._asdict())
        
        return response
    elif by == 'activity':
        matches = chemblTargetAccessionsByActivity(tuple([chembl_id]))
        
        fieldnames = ['accession', 'pref_name', 'activity_id', 'chembl_assay_id', 'description', 'organism', 'standard_value', 'standard_units', 'standard_flag', 'standard_type']
        header = {'accession': 'Accession ID',
                  'pref_name': 'Name',
                  'activity_id': 'Activity ID',
                  'chembl_assay_id': 'Assay ID',
                  'description': 'Description',
                  'organism': 'Organism',
                  'standard_value': 'Value',
                  'standard_units': 'Units',
                  'standard_flag': 'Flag',
                  'standard_type': 'Type'}
        
        response = HttpResponse(content_type='text/plain')
        if download:
            response['Content-Disposition'] = "attachment; filename=\"{chembl_id}_byActivity.tsv\"".format(chembl_id=chembl_id)
        writer = csv.DictWriter(response, fieldnames=fieldnames, dialect='excel-tab', extrasaction='ignore')
        
        writer.writerow(header)
        if len(matches) != 0:
            for m in matches[chembl_id]:
                writer.writerow(m._asdict())
        
        return response
    else:
        message = "Error: Invalid 'by' value. Should be 'annotation' or 'activity'."
        return HttpResponse(message, content_type='text/plain')

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

def getChemblTSV(request):
    download = True
    
    if 'accession_id' in request.GET:
        accession_id = request.GET['accession_id']
    else:
        message = "Error: Expected 'accession_id' value."
        return HttpResponse(message, content_type='text/plain')
    
    if 'download' in request.GET:
        if request.GET['download'] in ['true', 'True']:
            download = True
        else:
            download = False  
    
    matches = accessionToChembl(tuple([accession_id]))
    
    fieldnames = ['chembl_id']
    header = {'chembl_id': 'ChEMBL ID'}
    
    response = HttpResponse(content_type='text/plain')
    if download:
        response['Content-Disposition'] = "attachment; filename=\"{accession_id}_toChembl.tsv\"".format(accession_id=accession_id)
    writer = csv.DictWriter(response, fieldnames=fieldnames, dialect='excel-tab', extrasaction='ignore')
    
    writer.writerow(header)
    if len(matches) != 0:
        for m in matches[accession_id]:
            writer.writerow(m._asdict())
    
    return response

def readSources():
    sources = {}
    filename = os.path.join(settings.PROJECT_DIR,"unichem_sources.txt")
    with open(filename) as f:
        for line in f:
            (key, val) = line.split("\t")
            sources[int(key)] = val
    return sources




