from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.shortcuts import render_to_response
from django.template import RequestContext
from compounddb.models import Compound
from guest.decorators import guest_allowed, login_required
from django.contrib.auth import authenticate, login
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
import pybel

@guest_allowed
@cache_page(60 * 120)
@vary_on_cookie         # hide private compounds from other users 
def compound_detail(request, id, resource):
	username = request.user.username
	try:
		compound = Compound.objects.get(id__iexact=id, username=username)
	except Compound.DoesNotExist:
		raise Http404

	inchi = compound.inchi
	smiles = compound.smiles
	sdf = compound.sdffile_set.all()[0].sdffile
	
	if resource:
		if resource == 'smiles':
			return HttpResponse(smiles, mimetype='text/plain')
		elif resource == 'inchi':
			return HttpResponse(inchi, mimetype='text/plain')
		elif resource == 'sdf':
			return HttpResponse(sdf, mimetype='text/plain')
		elif resource == 'png':
			mymol = pybel.readstring("smi", str(smiles))
			png = mymol.write(format='png')
			return HttpResponse(png, mimetype='image/png') 

	return render_to_response('compound.html', dict(
		compound=compound,
		sdf=sdf,
		smiles=smiles,
		inchi=inchi,
		),
		context_instance=RequestContext(request))
