from django.http import HttpResponse, Http404, HttpResponseRedirect
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.core.urlresolvers import reverse

from tile import index, png, full, data, SupportError
def viewer(request, reference):
	try:
		context = index(request.session, reference, mode=request.GET.get('mode',
			''), dref=request.GET.get('dref'))
		return render_to_response("treeviewer.html", context,
			context_instance=RequestContext(request))
	except SupportError:
		return render_to_response("treeviewer.html", dict(error=
			"Sorry. The system encountered a problem while rendering your tree."),
			context_instance=RequestContext(request))

def serv_tile(request, img_ref):
    zoom = request.GET['zoom']
    tx = request.GET['tx']
    ty = request.GET['ty']
    pixels = png(request.session, img_ref, zoom, tx, ty)
    return HttpResponse(pixels, mimetype="image/png")
    
def serv_fullsize(request, img_ref):
    zoom = request.GET['zoom']
    pixels = full(request.session, img_ref, zoom)
    return HttpResponse(pixels, mimetype="image/png")

def userdata(request, reference):
    try:
        dref = data(request.session, reference, request.FILES['data'],
            name=request.POST['name'])
    except:
        raise 

    return HttpResponseRedirect(reverse('tree-viewer', args=[reference]) + '?dref=' + 
        dref + '&mode=' + request.GET.get('mode', ''))


