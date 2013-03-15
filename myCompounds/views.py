from django.http import HttpResponse, HttpResponseRedirect
from guest.decorators import guest_allowed, login_required
from django.template import RequestContext
from django.shortcuts import redirect, render_to_response
from compounddb import first_mol, InvalidInputError
# from compounddb.search import search
# from compounddb.views import get_library_by_name
from compounddb.tools import parse_annotation, insert_single_compound
from compounddb.models import Compound, SDFFile
from myCompounds.DownloadCIDs import DownloadCIDs
from django.contrib import messages
import random
import openbabel
import re
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_sdf, InputError, sdf_to_smiles
from django.conf import settings

MAX_COMPOUND_LIMIT = settings.MAX_COMPOUND_LIMIT
MAX_SDF_LENGTH = settings.MAX_SDF_LENGTH

@guest_allowed
def showCompounds(request, resource):
    # perform query for existing myCompounds
    page, matches = getMyCompounds(request)
    if resource:
	if resource == 'deleteAll':
		deleteMyCompounds(request)
		matches = None
		messages.error(request, 'All Compounds Deleted!')
	if resource == 'downloadSMILES':
		smiles = makeSMILES(request.user)
		return HttpResponse(smiles, mimetype='text/plain')
	if resource == 'downloadSDF':
		sdf = makeSDF(request.user)
		return HttpResponse(sdf, mimetype='text/plain')
    if matches:
        for match in matches:
	    match.smiles = re.match(r"^(\S+)", match.smiles).group(1) 
    return render_to_response('showCompounds.html', dict(p=page, matches=matches,), context_instance=RequestContext(request))    

@guest_allowed
def uploadCompound(request, *args, **kargs):
    if request.method == 'GET':
		return render_to_response('addCompounds.html', dict(input_mode='smiles-input'), context_instance=RequestContext(request))
    else:
		sdf = None
		name = None
		compid = None
		smiles = None
		if 'smiles' in request.POST:
			input_mode = 'smiles-input'
			sdf = u''
			try:
				smiles = request.POST['smiles'].split("\n")
				for line in smiles:
					if re.match(r"^\S+", line):
						sdf = sdf + smiles_to_sdf(str(line))
			except:
				messages.error(request, 'Error: Invalid SMILES string!')
				sdf = None
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
					smiles = smiles + " " + compid 
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

		    if (len(filteredCIDs) > 0):
			try:
				sdf = DownloadCIDs(filteredCIDs)
			except:
				messages.error(request, 'Invalid CIDs or no response from PubChem!')
				sdf = None
		    else:
				messages.error(request, 'Error: No valid CIDs entered!')
				sdf = None
		    
		if not sdf:
			return render_to_response('addCompounds.html', dict(
				input_mode=input_mode,
				post_data=request.POST,
				),
				context_instance=RequestContext(request))
		message = addMyCompounds(sdf, request.user)
		if re.search(r"^ERROR:", message):
			messages.error(request, message)
			return render_to_response('addCompounds.html', dict(
				input_mode=input_mode,
				post_data=request.POST,
				),
				context_instance=RequestContext(request))
		messages.success(request, message)
		return redirect(showCompounds, resource='')

def makeSDF(user):
	compoundList = Compound.objects.filter(user=user)
	sdf = u''
	for compound in compoundList:
		sdf = sdf + compound.sdffile_set.all()[0].sdffile.rstrip() + '\n'	
	return sdf 

def makeSMILES(user):
	compoundList = Compound.objects.filter(user=user)
	smiles = u''
	for compound in compoundList:
		smiles = smiles + compound.smiles.rstrip() + '\n'
	return smiles

def getMyCompounds(request):
	page = int(request.GET.get('p', '1'))
	matches = []
	# pure_query, matches = search('library: myCompounds', page, request)
	matches = Compound.objects.filter(user=request.user) 
	if len(matches) == 0:
		matches = None
	return page, matches

def addMyCompounds(sdf, user):
	sdffile = u''
	counter = 0
	linecounter = 0
	namekey = 'PUBCHEM_IUPAC_NAME'
	message = 'ERROR: bad input data.'
	added_ids = []
	# try:
	if 1 == 1:
		sdf = sdf.split("\n")
		for line in sdf:
			linecounter += 1
			if linecounter > MAX_SDF_LENGTH:
				message = "ERROR: an input sdf exceeds " + str(MAX_SDF_LENGTH) + " lines."
				raise Exception
			try:
				line = unicode(line, 'utf-8')
			except:
				pass
			sdffile += line
			sdffile += '\n'
			if line.startswith("$$$$"):
				try: 
					moldata = parse_annotation(sdffile, namekey)			
				except:
					message = "ERROR: invalid input format."	
					raise Exception
				if re.search("^unspecified_", moldata['id']):
					sdffile = moldata['id'] + sdffile
				counter += 1
				linecounter = 0
				if counter > MAX_COMPOUND_LIMIT:
					message = "ERROR: upload exceeds " + str(MAX_COMPOUND_LIMIT) + " compounds."
					raise Exception
				try:
					newid = insert_single_compound(moldata, sdffile, namekey, 'id', user)
				except:
					message = "ERROR: Database error, possibly excessively large compound?"
					raise Exception
				added_ids.append(newid)
				sdffile = u''
		if counter > 0:
			return "Success: Added " + str(counter) + " compounds."
		else:
			return "ERROR: No valid input found."
	try:
		pass
	except:
		for id in added_ids:
			try:
				Compound.objects.get(id=id).delete()	
			except:
				pass
		return message	
			
def getMW(sdf):
	if isinstance(sdf, unicode):
		sdf = sdf.encode('ascii', 'ignore')
	
	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats("sdf", "smi")
	mol = openbabel.OBMol()
	if not obConversion.ReadString(mol, sdf):
		raise InputError
	return str(mol.GetMolWt())
    
def getFormula(sdf):
	if isinstance(sdf, unicode):
		sdf = sdf.encode('ascii', 'ignore')
	
	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats("sdf", "smi")
	mol = openbabel.OBMol()
	if not obConversion.ReadString(mol, sdf):
		raise InputError
	return str(mol.GetFormula())
	
def getInChI(sdf):
	if isinstance(sdf, unicode):
		sdf = sdf.encode('ascii', 'ignore')
	
	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats("sdf", "InChI")
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
        compound = WorkbenchCompounds(compound = Compound.objects.get(cid__iexact=cid, library=library), user = user)
        compound.save()
