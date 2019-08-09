#!/usr/bin/python
# -*- coding: utf-8 -*-

from builtins import str
from builtins import object
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from django.template import RequestContext
from compounddb.models import Compound, Tag
from guest.decorators import guest_allowed, login_required
from django.contrib.auth import authenticate, login
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
from django.forms import ModelForm
from django.contrib import messages
import pybel
import re
from simplejson import dumps
import tempfile
import os
import json


class compoundForm(ModelForm):

    class Meta(object):

        model = Compound
        fields = ('cid', 'name')


                        # hide private compounds from other users

@guest_allowed
@cache_page(60 * 120)
@vary_on_cookie
def render_image(request, filename,id=None,cid=None ):
    compound = None
    try:
        if id is not None:
            compound = Compound.objects.get(id__iexact=id, user=request.user)
        elif cid is not None:
            compound = Compound.objects.get(cid__iexact=cid, user=request.user)

    except Compound.DoesNotExist:
        raise Http404
    if compound.weight > 2000:
        raise Http404
    smiles = re.match(r"^(\S+)", compound.smiles).group(1)
    mymol = pybel.readstring('smi', str(smiles))
    #this simple method below does not work as it returns the image as a unicode string
    # which cannot be converted back to raw bytes
    ##png = mymol.write(format='_png2')

    # so we write the image to a temp file and then read it back as bytes
    fp = tempfile.NamedTemporaryFile()
    fp.close() # deletes the file, but we still use the unique file name
    mymol.write(format='_png2',filename=fp.name)
    imageFile = open(fp.name,"rb")
    image = imageFile.read()
    imageFile.close()
    os.remove(fp.name)
    return HttpResponse(image, content_type='image/png')


@guest_allowed
def cid_lookup(request):
    if request.is_ajax():
        try:
            cid = request.GET['cid']
            compound = Compound.objects.get(cid=cid, user=request.user)
            response = dict(id=str(compound.id))
        except:
            response = dict(id='ERROR')
        return HttpResponse(dumps(response), 'text/json')
    else:
        raise Http404


@guest_allowed
def compound_detail(
    request,
    id,
    resource,
    filename,
    ):

    try:
        compound = Compound.objects.get(id__iexact=id,
                user=request.user)
    except Compound.DoesNotExist:
        raise Http404

    if request.method == 'POST':
        form = compoundForm(request.POST)
        if form.is_valid():
            compound.cid = form['cid'].value()
            compound.name = form['name'].value()
            compound.smiles = re.match(r"^(\S*)",
                    compound.smiles).group(1) + ' ' + compound.cid
            compound.save()
            sdf = compound.sdffile_set.all()[0]
            sdf.sdffile = compound.cid + '\n' + re.match(r"^.*?\n(.*)$"
                    , sdf.sdffile, flags=re.M | re.S).group(1)
            sdf.save()
            messages.success(request, 'Success: details updated.')
        else:
            messages.error(request, 'Error: invalid form data.')

    inchi = compound.inchi
    smiles = compound.smiles
    sdf = compound.sdffile_set.all()[0].sdffile

    if resource:
        if resource == 'smiles':
            return HttpResponse(smiles, content_type='text/plain')
        elif resource == 'inchi':
            return HttpResponse(inchi, content_type='text/plain')
        elif resource == 'sdf':
            return HttpResponse(sdf, content_type='text/plain')
        elif resource == 'delete':
            compound.delete()
            return HttpResponse('deleted', content_type='text/plain')
        elif resource == 'editform':
            form = compoundForm(instance=compound)
            return render(request,'genericForm.html',
                    dict(title='Edit Compound \'' + compound.cid + '\''
                    , form=form))

    return render(request,'compound.html', dict(compound=compound,
                              sdf=sdf, smiles=smiles, inchi=inchi))
@guest_allowed
def tagCompounds(request,action):
    if request.is_ajax() and request.method == 'POST':
        rawJson = request.body
        data = json.loads(rawJson.decode("utf-8"))
        compoundIds = []
        tags = []
        try:
            tags = data["tags"]
            compoundIds = [int(id) for id in data["ids"]]
        except:
            print("failed to parse json for tags and compound ids: "+str(rawJson))
            return HttpResponse("failed to parse JSON",status=404)

        #print("adding tags "+str(tags)+" to compounds "+str(compoundIds))

        Tag.ensureAllExist(tags,request.user)

        tagObjects = Tag.objects.filter(user=request.user,name__in=tags)
        compounds = Compound.objects.filter(user=request.user,id__in=compoundIds)
        for compound in compounds:
            if action == 'add':
                compound.tags.add(*tagObjects)
            elif action == 'remove':
                compound.tags.remove(*tagObjects)
            else:
                return HttpResponse("Unknown action given",status=404)
            compound.save()
    

    return HttpResponse('')

@guest_allowed
def batchOperation(request,action):
    #print("batch operation, action: "+action)
    if request.method == "POST" and action == "delete" :
        compoundsNotFound=[]
        compoundIds = request.POST.getlist("id[]")
        #print("deleteing compounds "+str(compoundIds))
        for compoundId in compoundIds:
            try:
                compound = Compound.objects.get(id__iexact=compoundId, user=request.user)
                compound.delete()
                #print("compound deleted")
            except Compound.DoesNotExist:
                print("compound "+str(compoundId)+" not found")
                compoundsNotFound.append(compoundId)
        if len(compoundsNotFound) != 0:
            return HttpResponse("Not all compounds could be deleted: "+str(compoundsNotFound),status=404)

    return HttpResponse('')

@guest_allowed
def countCompoundsWithTags(request,tags):
    return HttpResponse(Compound.byTagNames(tags.split(","),request.user).count())


