from django.shortcuts import render
from builtins import str
from django.http import HttpResponse, HttpResponseRedirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from guest.decorators import guest_allowed, login_required
from django.template import RequestContext
from django.shortcuts import redirect, render
from compounddb import first_mol, InvalidInputError


from compounddb.tools import parse_annotation, insert_single_compound
from compounddb.models import Compound, SDFFile,Tag
from pubchem_rest_interface.Pubchem_pug import DownloadCIDs
from django.contrib import messages


from django.http import Http404, HttpResponseBadRequest, HttpResponse
from django.utils.http import urlunquote, urlquote
import pybel
import re
import sys
import traceback
import os
import tempfile
from django.views.decorators.cache import cache_page
from simplejson import dumps,loads
from similarityworkbench import funcs

@guest_allowed
def showCompoundGrid(request ):

    # perform query for existing myCompounds

    matches = None
    givenTags = []
    allTags = [tag.name for tag in Tag.objects.filter(user = request.user)]

    if 'tags' in request.POST:
        givenTags = request.POST.getlist("tags")

    if "all" in givenTags or len(givenTags) == 0: 
        matches = Compound.objects.filter(user=request.user)
    else:
        tags = Tag.objects.filter(name__in = givenTags)
        matches = Compound.objects.filter(user=request.user,tags__in=tags)

    return render(request,'compoundGrid.html', dict(matches=matches,tags=allTags,currentTags=givenTags))


def ap(request):
    result = funcs.ap(request.POST.get("compound1"),
                      request.POST.get("compound2"));
    return HttpResponse(dumps(result), content_type='text/json')

def mcs(request):
    result = funcs.mcs(request.POST.get("compound1"),
                       request.POST.get("compound2"));
    return HttpResponse(dumps(result), content_type='text/json')

@guest_allowed
@cache_page(60 * 120)
def renderer(request, smiles):
    try:
        smiles = urlunquote(smiles)
        smiles = re.match(r"^(\S{1,2000})", smiles).group(1)
        mymol = pybel.readstring('smi', smiles)
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
    except:
        print("error rendering "+smiles)
        print("Unexpected error:", sys.exc_info())
        traceback.print_tb(sys.exc_info()[2])
        raise Http404


