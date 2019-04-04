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
        chembl_id = request.GET['chembl_id']
        result = chemblTargetAccessions(chembl_id)
        if len(result) == 0:
            matches = None
        else:
            matches = result[chembl_id]
    else:
        query_submit = False
        chembl_id = None
        matches = None
    
    context = {
        'query_submit': query_submit,
        'chembl_id': chembl_id,
        'matches': matches,
        }
    
    return render(request, 'targetsearch/getTargets.html', context)

def getChembl(request):
    if 'accession_id' in request.GET:
        query_submit = True
        accession_id = request.GET['accession_id']
        result = accessionToChembl(accession_id)
        if len(result) == 0:
            matches = None
        else:
            matches = result[accession_id]
    else:
        query_submit = False
        accession_id = None
        matches = None
    
    context = {
        'query_submit': query_submit,
        'accession_id': accession_id,
        'matches': matches,
        }
    
    return render(request, 'targetsearch/getChembl.html', context)
