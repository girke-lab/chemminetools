#!/usr/bin/python
# -*- coding: utf-8 -*-

from django.http import HttpResponse, HttpResponseRedirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from guest.decorators import guest_allowed, login_required
from django.template import RequestContext
from django.shortcuts import redirect, render_to_response
from compounddb import first_mol, InvalidInputError

# from compounddb.search import search
# from compounddb.views import get_library_by_name

from compounddb.tools import parse_annotation, insert_single_compound
from compounddb.models import Compound, SDFFile
from pubchem_soap_interface.DownloadCIDs import DownloadCIDs
from django.contrib import messages
import random
import openbabel
import re
import string
import time
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_sdf, \
    InputError, sdf_to_smiles
from django.conf import settings
from tools.runapp import createJob, updateJob
from django.utils.http import urlquote


@guest_allowed
def showCompounds(request, resource):

    # perform query for existing myCompounds

    (page, matches) = getMyCompounds(request)
    if resource:
        if resource == 'deleteAll':
            deleteMyCompounds(request)
            matches = None
            messages.error(request, 'All Compounds Deleted!')
	if resource == 'downloadSMILES.smi':
            smiles = makeSMILES(request.user)
            return HttpResponse(smiles, mimetype='text/plain')
        if resource == 'downloadSDF.sdf':
            sdf = makeSDF(request.user)
            return HttpResponse(sdf, mimetype='text/plain')
    if matches:
        for match in matches:
            match.smiles = re.match(r"^(\S+)", match.smiles).group(1)
            match.smiles = urlquote(match.smiles)
    return render_to_response('showCompounds.html', dict(p=page,
                              matches=matches),
                              context_instance=RequestContext(request))


@guest_allowed
def uploadCompound(request, resource = None, job_id = None):
    if (request.method == 'GET') and (resource != u'job'):
        return render_to_response('addCompounds.html',
                                  dict(input_mode='smiles-input'),
                                  context_instance=RequestContext(request))
    else:
        sdf = None
        name = None
        compid = None
        smiles = None
        if 'smiles' in request.POST:
            input_mode = 'smiles-input'
            sdf = u''
            try:
                smiles = request.POST['smiles'].split('\n')
                for line in smiles:
                    if re.match(r"^\S+", line):
                        sdf = sdf + smiles_to_sdf(str(line))
            except:
                messages.error(request, 'Error: Invalid SMILES string!')
                sdf = None
        elif resource == 'job':
            input_mode = 'sdf-upload'
            job = updateJob(request.user, job_id)
            f = open(job.output, 'r')
            sdf = f.read()
            f.close()
        elif 'sdf' in request.FILES:
            input_mode = 'sdf-upload'
            try:
                sdf = request.FILES['sdf']
                sdf = sdf.read()
            except (InputError, InvalidInputError):
                messages.error(request, 'Invalid SDF!')
                sdf = None
        elif 'sdf' in request.POST:
            if 'draw' in request.POST:
                input_mode = 'draw'
                sdf = request.POST['sdf'] + '$$$$'
                compid = str(request.POST['id'])
                compid = re.match(r"^(\S{0,20})", compid).group(1)
                try:
                    smiles = sdf_to_smiles(sdf)
                    smiles = re.match(r"^(\S+)", smiles).group(1)
                    smiles = smiles + ' ' + compid
                    sdf = smiles_to_sdf(smiles)
                except:
                    messages.error(request, 'Invalid drawing!')
                    sdf = None
            else:
                input_mode = 'sdf-input'
                sdf = request.POST['sdf']
                if not sdf:
                    messages.error(request, 'No input found!')
        elif 'pubchem' in request.POST:
            cids = request.POST['pubchem']
            cids = cids.split()
            input_mode = 'pubchem'
            filteredCIDs = []
            for cid in cids[:]:
                match = re.search("(\d{1,200})", cid)
                if match:
                    filteredCIDs.append(int(match.group(1)))

            if len(filteredCIDs) > 0:
                try:
                    sdf = DownloadCIDs(filteredCIDs)
                except:
                    messages.error(request,
                                   'Invalid CIDs or no response from PubChem!'
                                   )
                    sdf = None
            else:
                messages.error(request, 'Error: No valid CIDs entered!')
                sdf = None

        if not sdf:
            return render_to_response('addCompounds.html',
                    dict(input_mode=input_mode,
                    post_data=request.POST),
                    context_instance=RequestContext(request))
        newJob = createJob(request.user, 'Upload Compounds', '',
                           ['--user=' + str(request.user.id)], sdf)
        time.sleep(2)
        return redirect('tools.views.view_job', job_id=newJob.id,
                        resource='')


def makeSDF(user):
    compoundList = Compound.objects.filter(user=user)
    sdf = u''
    for compound in compoundList:
        sdf = sdf + compound.sdffile_set.all()[0].sdffile.rstrip() \
            + '\n'
    return sdf


def makeSMILES(user):
    compoundList = Compound.objects.filter(user=user)
    smiles = u''
    for compound in compoundList:
        smiles = smiles + compound.smiles.rstrip() + '\n'
    return smiles


def getMyCompounds(request):
    page = request.GET.get('page')
    result = []
    matches = Compound.objects.filter(user=request.user)
    paginator = Paginator(matches, 100)

    try:
        result = paginator.page(page)
    except PageNotAnInteger:

        # If page is not an integer, deliver first page.

        result = paginator.page(1)
    except EmptyPage:

        # If page is out of range (e.g. 9999), deliver last page of results.

        result = paginator.page(paginator.num_pages)

    if len(result) == 0:
        result = None
    return (paginator, result)


def getMW(sdf):
    if isinstance(sdf, unicode):
        sdf = sdf.encode('ascii', 'ignore')

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('sdf', 'smi')
    mol = openbabel.OBMol()
    if not obConversion.ReadString(mol, sdf):
        raise InputError
    return str(mol.GetMolWt())


def getFormula(sdf):
    if isinstance(sdf, unicode):
        sdf = sdf.encode('ascii', 'ignore')

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('sdf', 'smi')
    mol = openbabel.OBMol()
    if not obConversion.ReadString(mol, sdf):
        raise InputError
    return str(mol.GetFormula())


def getInChI(sdf):
    if isinstance(sdf, unicode):
        sdf = sdf.encode('ascii', 'ignore')

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('sdf', 'InChI')
    mol = openbabel.OBMol()
    if not obConversion.ReadString(mol, sdf):
        raise InputError
    return obConversion.WriteString(mol)


def deleteMyCompounds(request, cids=None):
    compoundsToDelete = Compound.objects.filter(user=request.user)
    if cids:
        for cid in cids:
            Compound.objects.filter(cid=cid).delete()
    else:
        compoundsToDelete.delete()


def addToWorkbench(user, cids):
    library = get_library_by_name('myCompounds')
    for cid in cids:
        compound = \
            WorkbenchCompounds(compound=Compound.objects.get(cid__iexact=cid,
                               library=library), user=user)
        compound.save()


