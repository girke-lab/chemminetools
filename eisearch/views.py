#!/usr/bin/python

from django.shortcuts import redirect, render_to_response
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.contrib import messages
from guest.decorators import guest_allowed, login_required
from django.template import RequestContext

@guest_allowed
def search(request):
    if request.method == 'GET':
        smi = ''
        if 'smi' in request.GET:
            smi = str(request.GET['smi'])
        return render_to_response('search.html', dict(mode='form',
            smi=smi),
            context_instance=RequestContext(request)) 
