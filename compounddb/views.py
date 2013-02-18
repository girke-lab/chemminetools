from django.http import Http404, HttpResponse, HttpResponseRedirect
from compounddb.models import Compound

def compound_detail(request, cid, resource):
	try:
		compound = Compound.objects.get(cid__iexact=cid, username=username)
	except Compound.DoesNotExist:
		raise Http404

	inchi = compound.inchi
	smiles = compound.smiles
	sdf = compound.sdffile_set.all()[0].sdffile

	return render_to_response('compound.html', dict(
		compound=compound,
		sdf=sdf,
		smiles=smiles,
		inchi=inchi,
		),
		context_instance=RequestContext(request))
