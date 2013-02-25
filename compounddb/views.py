from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.shortcuts import render_to_response
from django.template import RequestContext
from compounddb.models import Compound
from guest.decorators import guest_allowed, login_required
from django.contrib.auth import authenticate, login
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
from django.forms import ModelForm
from django.contrib import messages
import pybel

class compoundForm(ModelForm):
	class Meta:
		model = Compound
		fields = ('cid','name',)

@guest_allowed
@cache_page(60 * 120)
@vary_on_cookie         # hide private compounds from other users 
def render_image(request, id):
	username = request.user.username
	try:
		compound = Compound.objects.get(id__iexact=id, username=username)
	except Compound.DoesNotExist:
		raise Http404
	smiles = compound.smiles
	mymol = pybel.readstring("smi", str(smiles))
	png = mymol.write(format='png')
	return HttpResponse(png, mimetype='image/png') 

@guest_allowed
def compound_detail(request, id, resource):
	username = request.user.username
	try:
		compound = Compound.objects.get(id__iexact=id, username=username)
	except Compound.DoesNotExist:
		raise Http404

	if request.method == 'POST':
		form = compoundForm(request.POST)
		if form.is_valid():
			compound.cid = form['cid'].value()	
			compound.name = form['name'].value()
			compound.save()
			messages.success(request, 'Success: details updated.')
		else:
			messages.error(request, 'Error: invalid form data.')

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
		elif resource == 'delete':
			compound.delete()
			return HttpResponse("deleted", mimetype='text/plain')
		elif resource == 'editform':
			form = compoundForm(instance=compound)
			return render_to_response('genericForm.html', dict(
				title='Edit Compound \'' + compound.cid + '\'',
				form=form,
			),
			context_instance=RequestContext(request))

	return render_to_response('compound.html', dict(
		compound=compound,
		sdf=sdf,
		smiles=smiles,
		inchi=inchi,
		),
		context_instance=RequestContext(request))
