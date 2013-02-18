from django.db import models
from django.contrib.auth.models import User
from logging import root, basicConfig
import openbabel
import sys
import re
import tempfile
import os
import codecs
import md5
import compounddb.sdfiterator

cur_dir = os.path.dirname(__file__)

from compounddb.models import *

basicConfig()
inchiconv = openbabel.OBConversion()

############################
# library-related operations
############################
def get_library(name, version=0):
	""" only return library w/ given name and version
	    if not exist, return nothing w/o creating new library"""

	exist = False
	try:
		header = LibraryHeader.objects.get(name=name)
	except:
		#raise
		# header not exist
		return exist

	version = int(version)
	try:
		library = Library.objects.get(header=header, version=version)
	except:
		#raise
		# library not exist
		return exist

	return library

def create_library_w_header(header):
	""" given header, create library with latest version + 1"""

	# library_set w/ given header
	libs = header.library_set.all()

	# no library yet
	if len(libs) == 0:
		version = 0
	else:
		version = libs.latest().version + 1

	library = Library(header=header, version=version)
	library.save()

	root.warning('  --> : New library created: %s' % library)
	return library

def get_or_create_library(name, lib_ver=None):
	
	try:
		header = LibraryHeader.objects.get(name=name)
	except:
		header = LibraryHeader(name=name)
		header.save()

	try:
		library = Library.objects.get(header=header, version=lib_ver)
	except:
		try:
			cur_ver = header.library_set.latest().version
		except:
			cur_ver = -1
		library = Library(header=header, version=cur_ver+1)
		library.save()

	return library

#########################
# sdf-related processings
#########################
def get_sdf_tags(sdf):
	"""parse the sdf tags"""

	tag_pattern = re.compile(">\s+<([^>]+)>[^\n]*\n([^>$]+)")
	tags = tag_pattern.findall(sdf)
	tagdict = dict()

	# process each tag
	for name, value in tags:
		tagdict[name.strip()] = value.strip()
	return tagdict

def parse_annotation(sdf, namekey, idkey):
	""" parse annotation from SDF file """

	# parse the sdf tags
	moldata = get_sdf_tags(sdf)
	moldata[idkey] = moldata[idkey].replace(' ', '-')

	# --- inchi
	inchiconv.SetInAndOutFormats("sdf", "Inchi")
	mol = openbabel.OBMol()
	if '(noMol)' not in moldata[idkey]:
		res = inchiconv.ReadString(mol, codecs.encode(sdf, 'utf-8'))
		if mol.Empty():
			root.warning("  -->  ERROR on sdf")
			raise Exception

	# standard data generated
	# --- inchi/formula/weight
	moldata['inchi'] = inchiconv.WriteString(mol).strip()
	moldata['formula'] = mol.GetFormula()
	mol.AddHydrogens()
	moldata['weight'] = str(mol.GetMolWt())

	# if the name is not in sdf:
	if namekey == '__NA__':
		moldata[namekey] = ''
	
	# smiles
	inchiconv.SetInAndOutFormats("sdf", "smi")
	mol = openbabel.OBMol()
	if '(noMol)' not in moldata[idkey]:
		res = inchiconv.ReadString(mol, codecs.encode(sdf, 'utf-8'))
		moldata['smiles'] = inchiconv.WriteString(mol).strip()
	else:
		moldata['smiles'] = ''

	# check if necessary keys are present: the name and key
	if not moldata.has_key(idkey):
		print "compound data error: the idkey is not found"
		raise Exception
	if not moldata.has_key(namekey):
		# print "warning: no name key"
		moldata[namekey] = ""

	# do any necessary fixes for tags
	return moldata

############################
# single compound operations
############################
def _update_single_compound(moldata, sdf, library, nameky, idkey):

	# sdf file
	s = SDFFile(sdffile=sdf)
	s.save()
	sdfid = s.id

def insert_single_compound(moldata, sdf, library, namekey, idkey, username=''):
	""" insert single compound into database """

	cid = moldata[idkey]
	name = moldata[namekey]
	if '\n' in name:
		name = name.split('\n')[0]

	# compound
	c = Compound( cid=cid,
					name=name, 
					formula=moldata['formula'],
					weight=moldata['weight'],
					inchi=moldata['inchi'],
					smiles=moldata['smiles'],
					username=username )
					#sdf_file=s)
	c.save()
	c.library.add(library)
	c_id = c.id
	root.warning("  -->new compound inserted: c_id=%s, cid=%s" % (c_id, cid))

	# sdf file
	s = SDFFile(sdffile=sdf, compound=c)
	s.save()
	sdfid = s.id

	# annotation
	a_ids = []
	for name, value in moldata.items():
		if name in (namekey, idkey):
			continue
		a = Annotation(name=name, value=value, compound=c)
		a.save()
		aid = a.id
		a_ids.append(aid)

	return False

def link_compound_to_library(compound, library):
	compound.library.add(library)
	print compound
	print library

############################
# batch compounds operations
############################

def check_or_update_compounds(filename, library, namekey, idkey, check_only=True):
	"""Before updating a library, manually check which compounds have update """

	need_update = False

	# if update, create a new library first
	if check_only == False:
		new_lib = create_library_w_header(header=library.header)

	for sdf in sdfiterator.sdf_iter(codecs.open(filename, 'r', 'latin1')):
	#for sdf in sdfiterator.sdf_iter(filename):
		moldata = parse_annotation(sdf, namekey, idkey)
		cid = moldata[idkey]

		try:
			c = library.compound_set.get(cid=cid)
			sdf_old = SDFFile.objects.get(compound=c).sdffile
			
			if sdf != sdf_old:
				need_update = True
				
				if check_only:
					return need_update
				else:
					# update cmp: insert cmp actually
					insert_single_compound(moldata, sdf, new_lib, namekey, idkey)
			else:
				# same compound, link it to new library
				if not check_only:
					link_compound_to_library(c, new_lib)
		except Compound.DoesNotExist:
			insert_single_compound(moldata, sdf, new_lib, namekey, idkey)
			
	return False
	
def update_compounds(filename, library_name, version, namekey, idkey):
	"""update an existing library:
		 library version += 1 """

	library = get_library(library_name, version)

	if check_or_update_compounds(filename, library, namekey, idkey, check_only=True):
		"""if the 1st different compound found, update whole lib """
		root.warning('  --> : compound update found. Will update now.')
		check_or_update_compounds(filename, library, namekey, idkey, check_only=False)
	else:
		root.warning('  --> : No compound update found.')

def insert_compounds(filename, lib_name, namekey, idkey, lib_ver=0):
	""" Read SDF file, fetch each compound and process it.  
		All compounds will be inserted into a new library,
		If you are updating an existing library, 
		please use update_compounds()"""

	library = get_or_create_library(lib_name, lib_ver)

	f = open(filename)
	sdffile = u''
	count = 0
	for line in f.xreadlines():
		line = unicode(line, 'latin1')
		sdffile += line

		if line.startswith("$$$$"):
			moldata = parse_annotation(sdffile, namekey, idkey)
			insert_single_compound(moldata, sdffile, library, namekey, idkey)

			sdffile = u''
			count += 1
	root.warning('  -->  count:%s'%count)
	f.close()

###################################
# Physical Chemical Property Fields
###################################
def insert_property_fields(jfile):
	"""file format: NAME, DESCRIPTION, IS_INTEGER?, SOURCE_TAG"""
	"""e.g.: """
	"""HA1, Number of Hydrogen Bond Acceptors (HBA) 1, 1, count.HBA1"""
	fp = file(jfile)
	counter = 0
	for i in fp.readlines():
		counter += 1
		try:
			name, description, is_integer, source_tag = i.strip().split(',')[0:4]
		except:
			root.warning("  -->  line %s"%counter)
			raise
		insert_single_property_field(name, description, is_integer, source_tag)

	fp.close()

	return False

def insert_single_property_field(name, description, is_integer, source_tag):
	
	if is_integer:
		is_integer = True
	else:
		is_integer = False

	field = PropertyField(name=name, description=description, is_integer=is_integer, source_tag=source_tag)
	field.save()
	root.warning('  -->  id=%s'%field.id)

	return False

#####################################
# Physical Chemical Property - QSAR
#####################################
def prepare_qsar_property_file(library_name, library_ver, qsar_raw_file, output_file):
	library = get_library(library_name, library_ver)
	index = dict()
	fields = dict()

	fp = file(qsar_raw_file)
	tag_line = fp.readline()
	tags = tag_line.rstrip('\r').split('\t')
	id = tags.pop(0)
	c = 0
	for t in tags:
		fname = t.strip('\n').strip('"')
		try:
			fid = PropertyField.objects.get(name=fname).id
		except:
			print 'fname:',fname
			raise

		index[c] = fname
		fields[fname] = fid
		c += 1

	out = file(output_file, 'w')
	counter = 0
	for line in fp:
		values = line.strip().split('\t')
		cid = values.pop(0).strip('"')
		try:
			c_id = Compound.objects.get(cid=cid,library=library).id
		except:
			print cid
			print "counter: %s" % counter
			raise

		query = "INSERT INTO compounddb_property (compound_id, value, field_id)"
		c1 = 0
		for value in values:
			try:
				name = index[c1]
			except:
				print cid
				print name
				print value
				print index
				raise
			fid = fields[name]

			if c1 == 0:
				query += " VALUES (%s, %s, %s)" % (c_id, value, fid)
			else:
				query += ", (%s, %s, %s)" % (c_id, value, fid)

			c1 += 1
			counter += 1
		query += ";"
		#print query
		out.write('%s\n' % query)

	out.close()
	fp.close()

def format_qsar_property_file(lib_name, lib_ver, input, output, sdffile, id_tag):
	"""cerius2 program not giving correct cid for some
	   libraries, need to find out them manually
	   note: input file and sdffile have same number of compounds"""

	library = get_library(lib_name, lib_ver)
	fp = file(input)
	result = file(output, 'w')
	line1 = fp.readline()
	result.write(line1)
	for sdf in sdfiterator.sdf_iter(sdffile):
		tagdict = get_sdf_tags(sdf)
		cid = tagdict[id_tag]
		line = fp.readline()
		#line = line.split('\t',1)[1]
		result.write("%s\t%s" % (cid,line))

	fp.close()
	result.close()
	return

#def format_qsar_property_file_from_db(lib_id, input, output, file_tag):
#	"""cerius2 program not giving correct cid for some
#	   libraries, need to find out them manually
#	   note: input file and sdffile have diff. number of compounds """
#
#	from django.db import connection
#	c = connection.cursor()
#
#	fp = file(input)
#	out_fp = file(output, 'w')
#	line1 = fp.readline()
#	out_fp.write(line1)
#
#	for line in fp:
#		id_old = line.strip().split('\t')[0].strip('"')
#		if file_tag == 'name':
#			query = """select cid from compounddb_compound, compounddb_compound_library where compounddb_compound.name='%s' and compounddb_compound.id=compounddb_compound_library.compound_id and compounddb_compound_library.library_id=%s; """ % (id_old, lib_id)
#			c.execute(query)
#			raw = c.fetchall()
#			if len(raw) == 1:
#				cid = raw[0][0]
#			else:
#				print 'multiple compounds, id_old=%s' % id_old
#				print query
#				continue
#			out_fp.write('%s\t%s' % (cid, line.split('\t', 1)[1]))
#		else:
#			print 'to do'
#			return
#
#	fp.close()
#	out_fp.close()
#	return

def remove_property_in_file(lib_name, lib_ver, filename):
	"""in case wrong properties were inserted """

	library = get_library(lib_name, lib_ver)
	fp = file(filename, 'w')
	
	for c in library.compound_set.all():
			query = """delete from compounddb_property where compound_id=%s and field_id>42;""" % c.id
			fp.write('%s\n' % query)

	fp.close()
	return

def remove_duplicat_mw(lib_name, lib_ver):
	"""specifically MW, sometime inserted twice """
	
	library = get_library(lib_name, lib_ver)

	field = PropertyField.objects.get(name = 'MW')

	for c in library.compound_set.all():
		if c.property_set.filter(field=field).count() == 2:
			c.property_set.filter(field=field)[1].delete()
		elif c.property_set.filter(field=field).count() != 1:
			print '========exception: c_id=%s=====' % c.id

	return

#####################################
# Physical Chemical Property - JOELib
#####################################
def gen_joelib_property(sdf):
	"""run and parse the property output """

	# save the input in FS
	t = tempfile.NamedTemporaryFile(suffix=".sdf")
	t.write(codecs.encode(sdf, 'utf-8'))
	t.flush()

	# prepare the output file
	(f, out) = tempfile.mkstemp(suffix='.sdf')
	os.close(f)

	# convert
	cmd = """JAVA_HOME=%s/support/jre1.6.0_06/ JOELIB2=%s/support/JOELib2-alpha-20070303/ %s/support/JOELib2-alpha-20070303/moleculeConversion.sh +d +h -iSDF -osdf "%s" "%s" > /dev/null""" % (cur_dir, cur_dir, cur_dir, t.name, out)
	if os.uname()[1] == 'chemmineweb': 
		cmd = """JAVA_HOME=/usr/lib/jvm/java-6-sun JOELIB2=%s/support/JOELib2-alpha-20070303/ %s/support/JOELib2-alpha-20070303/moleculeConversion.sh +d +h -iSDF -osdf "%s" "%s" > /dev/null""" % (cur_dir, cur_dir, t.name, out)
	elif os.uname()[1] == 'biocluster': 
		cmd = """JAVA_HOME=/usr/lib/jvm/java-1.5.0-sun-1.5.0.14/ JOELIB2=%s/support/JOELib2-alpha-20070303/ %s/support/JOELib2-alpha-20070303/moleculeConversion.sh +d +h -iSDF -osdf "%s" "%s" > /dev/null""" % (cur_dir, cur_dir, t.name, out)
	root.warning("  -->   running:%s" % cmd)
	if os.system(cmd) != 0:
		os.unlink(out)
		raise "cannot run JOELib"
	
	# read and parse
	f = file(out)
	tags = get_sdf_tags(codecs.decode(f.read(), 'utf-8'))
	f.close()

	# clean
	os.unlink(out)
	return tags

def _insert_property(compound, propertydata, tags):
	"""insert JOELib property into DB"""

	count = 0
	for key in propertydata.keys():
		if key in tags:
			try:
				j = Property(field=PropertyField.objects.get(source_tag=key), 
					value=propertydata[key],
					compound=compound)
				j.save()
			except:
				#root.warning("!!!!!!!!!!!key=%s" % key)
				raise
			count += 1
	return False
			
def insert_property(library_name, library_ver, offset, limit, filename):
	"""for JOELib property
	   good for small library as it is slow
	   for library with large amount of compounds:
	    1. use extra scripts to generate Property files
		2. use insert_property_from_file() """
	
	def do_process(compounds, filep, tags):
		for cmp in compounds:
			sdf = SDFFile.objects.get(compound=cmp).sdffile
			print (' running compound:%s' % cmp.id)
			try:
				propertydata = gen_joelib_property(sdf)
				_insert_property(cmp, propertydata, tags)
			except:
				root.warning('  -->  cannot run property on %s' % cmp.id)
				raise

	library = get_library(library_name, library_ver)

	try:
		offset = int(offset)
		limit = int(limit)
	except:
		raise

	# number of compounds
	filep = file(filename, 'w')
	compounds = library.compound_set.all()
	count = len(compounds)
	root.warning('  -->  %s compounds to be processed' % count)

	tags = []
	for i in PropertyField.objects.values_list('source_tag'):
		tags.append(i[0])

	limit = limit + offset
	count = 0
	while(offset < limit):
		cur_limit = 1000
		#cur_limit = 7 # for testing only
		if offset + cur_limit < limit:
			cur_limit = limit - offset
		root.warning('  --> cur_limit:%s'%cur_limit)

		filep_progress = file(filename+'.progress', 'w')
		filep_progress.write('working on %s - %s\n' % (offset, offset + cur_limit))
		filep_progress.flush()

		# process
		cmps = compounds[offset:offset+cur_limit]
		do_process(cmps, filep, tags)

		filep_progress.close()
		count += 1
		offset += cur_limit

	filep.close()
	return False

def insert_property_from_file(library_name, library_ver, ID_tag, filename):
	"""for JOELib Property
	   input files: parsed name-value pairs in SDF format
	   name-value pairs: ID and property only """

	import datetime
	begin = datetime.datetime.now()
	print 'starts at: %s' % begin

	library = get_library(library_name, library_ver)
	for sdf in sdfiterator.sdf_iter(filename):
		cid = ''
		tagdict = get_sdf_tags(sdf)

		cid = tagdict.pop(ID_tag)
		compound = Compound.objects.get(library=library,cid=cid)

		for tag, value in tagdict.items():
			try:
				field = PropertyField.objects.get(source_tag=tag)
			except:
				raise
			p = Property(field=field,
						value=value,
						compound=compound)
			p.save()

	end = datetime.datetime.now()
	print 'finished at: %s' % end
	print 'time lapsed: %s' % (end-begin)

	return

def prepare_property_file(library_name, library_ver, ID, input, output):
	"""for JOELib property
	   input: parsed name-value pairs in SDF format
	   output: result file to write SQL sentences
	   name-value pair: ID and value only"""
	import datetime
	begin = datetime.datetime.now()
	print 'now begin at: %s' % begin

	library = get_library(library_name, library_ver)
	fields = {}
	for f in PropertyField.objects.all():
		fields[f.source_tag] = f.id

	fp = file(output, 'w')

	for sdf in sdfiterator.sdf_iter(input):
		cid = ''
		tagdict = get_sdf_tags(sdf)
		cid = tagdict.pop(ID)
		try:
			c_id = Compound.objects.get(library=library,
								cid=cid).id
		except:
			print '-------exception!!!-----------'
			print Compound.objects.filter(library=library, cid=cid)
			continue

		for tag,value in tagdict.items():
			try:
				field_id = fields[tag]
			except:
				print '-------field exception!!!-----------'
				print 'cid:%s, c_id:%s' % (cid, c_id)
				print 'tag:%s\n\n' % tag
				raise
			query = """INSERT INTO compounddb_property (compound_id, value, field_id) VALUES (%s, %s, %s);""" % (c_id, value, field_id)
			fp.write("%s\n" % query)

	fp.close()
	end = datetime.datetime.now()
	print 'finished at: %s' % end
	print 'time lapsed: %s' % (end-begin)

	return
################
# plate and well
################
def fast_insert_plate_test(library_name, library_ver, platefile, format):

	format = int(format)
	library = get_library(library_name, library_ver)

	info = dict()
	fp = file(platefile)

	for i in fp:
		cid = ''
		plate = ''
		well = ''

		cid, plate, well = i.strip().split(',')[0:3]
		cid = str(cid)
		if cid == 'NA':
			continue
		plate = int(plate)
		well = str(well)

		info[cid] = '%s,%s' % (plate, well)

	counter = 0
	for cmp in library.compound_set.all():
		try:
			plate, well = info[cmp.cid].split(',')[0:2]
		except:
			raise
		p = Plate(format=format,
					plate=plate,
					well=well,
					compound=cmp)
		p.save()
		counter += 1
		print counter

	fp.close()

def fast_insert_plate_from_file(library_name, library_ver, platefile, format):
	"""if we have a file containing plate mapping info,
		use plate info from the file.
		'library': library name/library id
		'platefile' format:
			1st column: compound id
			2nd column: plate
			3rd column: well
		'format': 24/48/96/384 """

	format = int(format)
	library = get_library(library_name, library_ver)

	fp = file(platefile)
	count = 0
	for i in fp.xreadlines():
		cid = ''
		plate = ''
		well =  ''

		cid, plate, well = i.strip().split(',')[0:3]
		cid = str(cid)
		if cid == 'NA':
			continue
		plate = int(plate)
		well = str(well)
		#sys.stdout.write('\r  -->  cid:%s' % (cid, ))
		compound = Compound.objects.get(library=library, cid=cid)
		try:
			p = Plate(format=format, 
					plate=int(plate),
					well=well,
					compound=compound)
			p.save()
			count += 1
		except:
			print("  -->  count:%s, cid:%s" % (count, cid) )
			raise
	root.warning("  -->  count:%s" % (count, ) )
	fp.close()
	return

def insert_plate_from_db(library_name, library_ver, format, plate_key, plate_base, well_key, row_key, col_key):
	"""if no file about plate mapping info,
		use plate info from sdf in db
		'format': 24/48/96/384
		'library': library name/id
		'plate_key': required
		'plate_base': '__NA__'/required, e.g. 'LAT', to be stripped from plate
		'well_key': optional, if no, must be '__NA__'
		'row_key' 'col_key': if no well_key, these 2 keys are required,
		                     if has well_key, will be ignored
		e.g.:
		insert_plate_from_db('Latca', 96, 'PLATE', 'LAT', '__NA__','ROW', 'COL')"""

	library = get_library(library_name, library_ver)
	compounds = Compound.objects.filter(library=library)
	format = int(format)

	count = 0
	for cmp in compounds:
		plate = Annotation.objects.get(compound=cmp, name=plate_key).value
		if plate_base is not '__NA__':
			try:
				plate = int(plate.lstrip(plate_base))
			except:
				root.warning('  -->  plate:%s, base:%s' % (plate, plate_base))
		if well_key is not "__NA__":
			well = Annotation.objects.get(compound=cmp, name=well_key).value
		else:
			row = Annotation.objects.get(compound=cmp, name=row_key).value
			col = Annotation.objects.get(compound=cmp, name=col_key).value
			if len(col) == 1:
				col = "0%s" % col
			well = "%s%s" % (row, col)

		try:
			p = Plate(format=format,
					plate=int(plate),
					well=well,
					compound=cmp)
			p.save()
			count += 1
			root.warning('  -->  count:%s, compound:%s' % (count, cmp) )
		except:
			raise
			root.warning('  -->  count:%s, compound:%s' % (count, cmp) )

	root.warning('  -->  total count:%s' % count )
	return

def insert_plate_from_file(filename, library_name, library_ver, format, id_col, plate_col, well_col):
	"""for libraries using UCR own plate info, stored in file 
	   columbe number starts from 0"""

	import csv

	library = get_library(library_name, library_ver)

	id_col = int(id_col)
	format = int(format)
	plate_col = int(plate_col)
	well_col = int(well_col)

	rd = csv.reader(open(filename), delimiter=',')
	rd.next()

	counter = 0
	for data in rd:
		cid = data[id_col]
		if cid == 'NA':
			continue
		plate = data[plate_col]
		well = data[well_col]

		try:
			cmp = Compound.objects.get(library=library, cid=cid)
		except:
			print cid
			raise

		try:
			p = Plate(format=format,
					plate=int(plate),
					well=well,
					compound=cmp)
			p.save()
		except:
			print "cid:%s, plate:%s, well:%s" % cid, plate, well
			raise
		counter += 1
		print "\r		counter: %s" % counter

	print "in total %s plates inserted" % counter
	return


def insert_fingerprint(library_name, library_ver, namekey=None, idkey = None):
	"""calculate fingerprint and store in DB """
	import pybel
	#fpt = '00000000 00000000 00000100 00010200 00000040 0000000000000000 00000000 00000100 00080840 08000000 4000900011008000 00000000 40040000 00080008 08000000 0200000001000000 00000000 00000000 08000a00 01000000 0040001000000000 80000000 00000000 00000000 00000000 000200020000a000 00000400'
	hash = md5.new() # '57e27aea58c3ec6abcc3c77963ce88f9' for this fingerprint
	id = None

	library = get_library(library_name, library_ver)
	compounds = library.compound_set.all()
	count = 0
	for cmp in compounds:
		count += 1
		if len(cmp.fingerprint_set.all()) > 0:
			print 'skipping #%s'%count
			continue
		try:
			sdf = SDFFile.objects.get(compound=cmp).sdffile
			if idkey:
				moldata = parse_annotation(sdf, namekey, idkey)
				id = moldata[idkey]
			if id and '(noMol)' in id:
				continue
			mol = pybel.readstring('sdf', codecs.encode(sdf, 'utf-8'))
			fpt = mol.calcfp()
			fpt_string = ' '.join([str(i) for i in fpt.bits])
			hash.update(fpt_string)
			fpt_md5 = hash.hexdigest()
			#print fpt_md5

			fp = Fingerprint(compound=cmp,
					md5 = fpt_md5)
			fp.save()
			print '\r    %s'%count
		except:
			print cmp
			raise

	return

######
# MISC
######

def update_mw(lib_name, lib_ver, input, rev=False):
	"""goal: to update MW value with hydrogen added 
	   when calculating JOELib
	   'rev': in ChemMineV2, some libraries got compound ID and compound name switched, like 'Aurora'"""

	import datetime
	begin = datetime.datetime.now()
	print 'starts at: %s' % begin

	library = get_library(lib_name, lib_ver)
	mw = PropertyField.objects.get(name='MW')

	fp = file(input)
	line1 = fp.readline()

	count = 1
	for line in fp:
		cid, weight = line.strip().split('\t')
		try:
			if rev:
				c = Compound.objects.get(library=library,
										name=cid)
			else:
				c = Compound.objects.get(library=library,
									cid=cid)
		except Compound.DoesNotExist:
			print 'not found: line %s, cid=%s' % (count, cid)
			pass

		try:
			p = Property.objects.get(compound=c, field=mw)
			p.value = weight
			p.save()
		except Property.DoesNotExist:
			p = Property(field=mw, 
						compound=c, 
						value=weight)
			p.save()
			print 'new p for %s, line %s' % (cid, count)
		except:
			print '----->line %s, cid=%s' % (count, cid)
			pass

		count += 1
		#print "%s: %s -> %s", (cid, old, weight)

	fp.close()

	end = datetime.datetime.now()
	print 'ends at: %s' % end

	return

def del_duplicate_mw(lib_name, lib_ver):
	"""some libraries has 2 mw """

	library = get_library(lib_name, lib_ver)
	mw = PropertyField.objects.get(name='MW')

	for c in library.compound_set.all():
		if c.property_set.filter(field=mw).count() == 2:
			c.property_set.filter(field=mw)[1].delete()

	return

def fix_kegg_cid():
	"""some cid in KEGG still has '(noMol)', fix them"""

	library = get_library('KEGG', 0)
	
	count = 0
	for c in library.compound_set.all():
		if '(noMol)' in c.cid:
			old = c.cid
			print old
			c.cid = old.strip('(noMol)')
			c.save()
			count += 1

	print '%s compounds updated with new cid' % count
	return


def format_sdf_for_qsar(sdffile, output, ID_tag):
	"""Cerius2 uses 1st line in SDF as ID tag
	   some sdf has blank 1st line, so we need to format SDF
	   by filling cid to 1st line in SDF"""
	
	fp = file(output, 'w')
	for sdf in sdfiterator.sdf_iter(sdffile):
		tagdict = get_sdf_tags(sdf)
		cid = tagdict[ID_tag]

		fp.write('%s\n' % cid)
		fp.write(sdf.split('\n',1)[1].split('M  END')[0])
		fp.write('M  END\n')
		fp.write('> <%s>\n%s\n\n' % (ID_tag, cid))
		fp.write('$$$$\n')

	fp.close()
	return

def list_all_cid_from_sdf(sdffile, ID_tag, outfile):
	fp = file(outfile, 'w')
	for sdf in sdfiterator.sdf_iter(sdffile):
		tagdict = get_sdf_tags(sdf)
		cid = tagdict[ID_tag]
		fp.write('%s\n' % cid)

	fp.close()
	return
