#!/usr/bin/python
# -*- coding: utf-8 -*-

from builtins import str
from django.http import HttpResponse, HttpResponseRedirect, JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from guest.decorators import guest_allowed, login_required
from django.template import RequestContext
from django.shortcuts import redirect, render
from compounddb import first_mol, InvalidInputError

# from compounddb.search import search
# from compounddb.views import get_library_by_name

from compounddb.tools import parse_annotation, insert_single_compound
from compounddb.models import Compound, SDFFile,Tag
from pubchem_rest_interface.Pubchem_pug import DownloadCIDs
from django.contrib import messages
import random
import openbabel
import re
import string
import sys
import traceback
import time
import tools
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_sdf, \
    InputError, sdf_to_smiles
from django.conf import settings
from tools.models import Job
from tools.runapp import createJob, updateJob
from django.utils.http import urlquote
from targetsearch.chembl_helpers import get_chembl_sdfs,get_chembl_smiles




@guest_allowed
def showCompounds(request, resource):

    # perform query for existing myCompounds


    #(page, matches) = getMyCompounds(request)
    matches = Compound.objects.filter(user=request.user)
    #print("page count: "+str(page.count))
    #print("page num pages: "+str(page.num_pages))
    if resource:
        if resource == 'deleteAll':
            deleteMyCompounds(request)
            matches = None
            messages.error(request, 'All Compounds Deleted!')
        if resource == 'downloadSMILES.smi':
            smiles = makeSMILES(request.user)
            return HttpResponse(smiles, content_type='text/plain')
        if resource == 'downloadSDF.sdf':
            sdf = makeSDF(request.user)
            return HttpResponse(sdf, content_type='text/plain')
    if matches:
        for match in matches:
            match.smiles = re.match(r"^(\S+)", match.smiles).group(1)
            match.smiles = urlquote(match.smiles)
    allTags = Tag.allUserTagNames(request.user)
    return render(request,'showCompounds.html', dict(matches=matches,tags=allTags))

@guest_allowed
def downloadCompounds(request, outputFormat):
    params = None
    compounds = None
    if request.method == "GET":
        params = request.GET
    elif request.method == "POST":
        params = request.POST
    else:
        return HttpResponse("dowload only supported via GET or POST",status=501)
        
    if 'chembl_id' in params: #tested
        cids = params.getlist('chembl_id')
        try:
            if outputFormat == "sdf":
                sdfs = get_chembl_sdfs(tuple(cids))
                sdf = "\n$$$$\n".join(sdfs)+"\n$$$$\n"
                return fileDownloadResponse(sdf,"chembl.sdf")
            elif outputFormat == "smi":
                smiles = get_chembl_smiles(tuple(cids))
                data = "\n".join(smiles)+"\n"
                return fileDownloadResponse(data,"chembl.smi")
        except:
            print("Unexpected error:", sys.exc_info())
            traceback.print_tb(sys.exc_info()[2])
            return HttpResponse("Failed to fetch ChEMBL compounds",status=500)
    elif 'cid' in params or 'tag' in params: # NOT YET TESTED
        if 'cid' in params:
            cids = params.getlist('cid')
            compounds=Compound.objects.filter(user=request.user, cid__in=cids)
        else:
            tagNames= params.getlist('tag')
            tags = Tag.objects.filter(user=request.user, name__in=tagNames)
            compounds=Compound.objects.filter(user=request.user, tags__in=tags)

        if outputFormat == "sdf":
            sdf ="\n".join( 
                    [ compound.sdffile_set.all()[0].sdffile.rstrip() for compound in compounds])+"\n"
            return fileDownloadResponse(sdf,"download.sdf")
        elif outputFormat == "smi":
            smiles ="\n".join( 
                    [ compound.smiles.rstrip() for compound in compounds])+"\n"
            return fileDownloadResponse(smiles,"download.smi")
    elif 'pubchem_id' in params: # NOT YET TESTED
        cids = params.getlist('pubchem_id')
        job = createJob(request.user,'pubchemID2SDF','',[],cids,outputFormat,wait=True)
        return redirect(tools.views.view_job,job_id=job.id,resource="other",filename="download."+outputFormat)
    else:
        return HttpResponse("unknown search paramaters",status=400)
    

def fileDownloadResponse(content,filename):
    response = HttpResponse(content,content_type="application/force-download")
    response["Content-Disposition"]="attachement; filename=\""+filename+"\""
    return response


   


@guest_allowed
def uploadCompound(request, resource = None, job_id = None):
    allTags = Tag.allUserTagNames(request.user)


    if (request.method == 'GET') and (resource != u'job'):
        return render(request,'addCompounds.html',
                                  dict(input_mode='smiles-input',
                                      tags=allTags))
    else:
        sdf = None
        name = None
        compid = None
        smiles = None
        compoundTags = []
        input_mode='smiles-input'

        if 'tags' in request.POST:

            compoundTags = set(request.POST.getlist('tags'))

            Tag.ensureAllExist(compoundTags,request.user)

            #existingTags = set(allTags)
            #print("compound tags: "+str(compoundTags))

            #for newTag in compoundTags.difference(existingTags):
            #    print("creating new tag: "+newTag+" for user "+request.user.username)
            #    Tag.objects.create(name = newTag, user=request.user)


        if 'smiles' in request.POST:
            input_mode = 'smiles-input'
            sdf = u''
            try:
                smiles = request.POST['smiles'].split('\n')
                for line in smiles:
                    if re.match(r"^\S+", line):
                        sdf = sdf + smiles_to_sdf(str(line))
            except:
                print("Unexpected error:", sys.exc_info())
                traceback.print_tb(sys.exc_info()[2])
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
                sdf = sdf.read().decode("utf-8")
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
                    print("Unexpected error:", sys.exc_info())
                    traceback.print_tb(sys.exc_info()[2])
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
                    sdf = DownloadCIDs(cids)
                except:
                    print("Unexpected error:", sys.exc_info())
                    traceback.print_tb(sys.exc_info()[2])
                    messages.error(request,
                                   'Invalid CIDs or no response from PubChem!'
                                   )
                    sdf = None
            else:
                messages.error(request, 'Error: No valid CIDs entered!')
                sdf = None
        elif 'chembl' in request.POST:
            cids = tuple(request.POST['chembl'].split())

            if len(cids) > 0:
                try:
                    sdfs = get_chembl_sdfs(cids)
                    sdf = "\n$$$$\n".join(sdfs)+"\n$$$$\n"
                except:
                    print("Unexpected error:", sys.exc_info())
                    traceback.print_tb(sys.exc_info()[2])
                    messages.error(request,
                                   'Invalid CIDs or no response from ChEMBL!'
                                   )
                    sdf = None
            else:
                print("no chembl cids given")
                messages.error(request, 'Error: No valid ChEMBL CIDs entered!')
                sdf = None


        if not sdf:
            return render('addCompounds.html',
                    dict(input_mode=input_mode,
                    post_data=request.POST,
                    tags=compoundTags))
        newJob = createJob(request.user, 'Upload Compounds', '',
                           ['--user=' + str(request.user.id),"--tags="+(",".join(compoundTags))], sdf)
        time.sleep(2)
        return redirect(tools.views.view_job, job_id=newJob.id,
                        resource='')

@guest_allowed
def ajax(request):
    compoundTags = set()
    
    def die(msg):
        ajaxResponse = { 'success' : False, 'message' : msg }
        return JsonResponse(ajaxResponse)
    
    if 'tags' in request.POST:
        compoundTags = set(request.POST.getlist('tags'))
        Tag.ensureAllExist(compoundTags, request.user)
    
    if 'action' not in request.POST:
        return die('Expected "action" in request.')
    action = request.POST['action']
    
    if action == 'add':
        if 'source_id' not in request.POST:
            return die('Expected "source_id" in request.')
        source_id = request.POST['source_id']
        
        if 'ids' not in request.POST:
            return die('No "ids" given in request.')
        ids = tuple(request.POST['ids'].split())
        
        if len(ids) == 0:
            return die('No "ids" given in request.')
        
        if source_id == 'chembl':
            try:
                sdfs = get_chembl_sdfs(ids)
                sdf = '\n$$$$\n'.join(sdfs) + '\n$$$$\n'
                
                newJob = createJob(request.user, 'Upload Compounds', '',
                                   ['--user=' + str(request.user.id), '--tags=' + (','.join(compoundTags))], sdf)
                time.sleep(2)
                newJob = updateJob(request.user, newJob.id)
                if newJob.status == Job.RUNNING:
                    ajaxResponse = { 'success' : True, 'message' : 'Compound upload in progress. Check "Past Jobs" for status.' }
                elif newJob.status == Job.FINISHED:
                    ajaxResponse = { 'success' : True, 'message' : 'Compounds uploaded successfully.'}
                else:
                    ajaxResponse = { 'success' : False, 'message' : 'An error occurred while uploading your compounds.' }
                return JsonResponse(ajaxResponse)
            except Exception as e:
                return die(str(e))
        else:
            return die('Unknown source_id: {}'.format(source_id))
    else:
        return die('Unknown action: {}'.format(action))


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
    #if isinstance(sdf, str):
        #sdf = sdf.encode('ascii', 'ignore')

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('sdf', 'smi')
    mol = openbabel.OBMol()
    if not obConversion.ReadString(mol, sdf):
        raise InputError
    return str(mol.GetMolWt())


def getFormula(sdf):
    #if isinstance(sdf, str):
        #sdf = sdf.encode('ascii', 'ignore')

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('sdf', 'smi')
    mol = openbabel.OBMol()
    if not obConversion.ReadString(mol, sdf):
        raise InputError
    return str(mol.GetFormula())


def getInChI(sdf):
    #if isinstance(sdf, str):
        #sdf = sdf.encode('ascii', 'ignore')

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


