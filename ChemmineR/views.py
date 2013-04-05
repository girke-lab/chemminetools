#Create your views here.
from django.shortcuts import get_object_or_404, get_list_or_404, render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect, HttpResponse
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from eis.pubchemdl import download
from SimilaritySearch import SimilaritySearch
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_sdf, InputError, sdf_to_smiles
from django.views.decorators.csrf import csrf_exempt

@csrf_exempt
def runapp(request, url):	
	if request.method == 'POST':
		app = str(request.GET.get('app', '1'))
		if app == 'getIds':
			try:
				# parse and test input
				cids = str(request.POST['cids'])
				# launch app
				response = getIds(cids)
			except:
				response = HttpResponse("ERROR: no results or invalid query", mimetype="text/plain")
		elif app == 'searchString':
			try:
				# parse and test input
				smiles = str(request.POST['smiles'])
				# launch app
				response = searchString(smiles)
			except:
				response = HttpResponse("ERROR: no results or invalid query", mimetype="text/plain")
		elif app == 'sdf2smiles':
			try:
				# parse and test input
				sdf = str(request.POST['sdf'])
				# launch app
				response = sdf2smiles(sdf)
			except:
				response = HttpResponse("ERROR: no results or invalid query", mimetype="text/plain")
		elif app == 'smiles2sdf':
			try:
				# parse and test input
				smiles = str(request.POST['smiles'])
				# launch app
				response = smiles2sdf(smiles)
			except:
				response = HttpResponse("ERROR: no results or invalid query", mimetype="text/plain")
		else:
			response = HttpResponse("ERROR: invalid app specified", mimetype="text/plain")
		return response
	else:
		return HttpResponse("ERROR: query must be an HTTP POST\n", mimetype="text/plain")

def getIds(cids):
	cids = cids.split(',')
	cids = [int(cid) for cid in cids]
	response = download(cids)
	return HttpResponse(response, mimetype="text/plain")

def searchString(smiles):
	id_list = SimilaritySearch(smiles)
	id_list = [str(cid) for cid in id_list]
	id_list = ','.join(id_list)
	return HttpResponse(id_list, mimetype="text/plain")

def smiles2sdf(smiles):
	sdf = smiles_to_sdf(smiles)
	return HttpResponse(sdf, mimetype="text/plain")

def sdf2smiles(sdf):
	sdf = sdf_to_smiles(sdf)
	return HttpResponse(sdf, mimetype="text/plain")
