from django.shortcuts import render

from .chembl_helpers import (
    chemblTargetAccessions,
    accessionToChembl
    )

def home(request):
    return render(request, 'targetsearch/home.html')

def getTargets(request):
    if 'chembl_id' in request.GET:
        query_submit = True
        chembl_id_str = request.GET['chembl_id']
        chembl_id_tuple = tuple(chembl_id_str.split())
        matches = chemblTargetAccessions(chembl_id_tuple)
    else:
        query_submit = False
        chembl_id_tuple = None
        matches = None
    
    context = {
        'query_submit': query_submit,
        'chembl_id_tuple': chembl_id_tuple,
        'matches': matches,
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
