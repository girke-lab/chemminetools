#!/usr/bin/python
# -*- coding: utf-8 -*-

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
import string
import random

cur_dir = os.path.dirname(__file__)

from compounddb.models import *

basicConfig()
inchiconv = openbabel.OBConversion()


#########################
# sdf-related processings
#########################

def get_sdf_tags(sdf):
    """parse the sdf tags"""

    tag_pattern = re.compile(""">\s+<([^>]+)>[^
]*
([^>$]+)""")
    tags = tag_pattern.findall(sdf)
    tagdict = dict()

    # process each tag

    for (name, value) in tags:
        tagdict[name.strip()] = value.strip()
    return tagdict


def parse_annotation(sdf, namekey):
    """ parse annotation from SDF file """

    # parse the sdf tags

    moldata = get_sdf_tags(sdf)

    # --- inchi

    inchiconv.SetInAndOutFormats('sdf', 'Inchi')
    mol = openbabel.OBMol()
    res = inchiconv.ReadString(mol, codecs.encode(sdf, 'utf-8'))
    if mol.Empty():
        root.warning('  -->  ERROR on sdf')
        raise Exception

    # standard data generated
    # --- inchi/formula/weight

    moldata['inchi'] = inchiconv.WriteString(mol).strip()
    moldata['formula'] = mol.GetFormula()
    moldata['id'] = mol.GetTitle()
    if moldata['id'] == '':
        moldata['id'] = 'unspecified_' \
            + ''.join(random.sample(string.digits, 6))
    mol.AddHydrogens()
    moldata['weight'] = str(mol.GetMolWt())

    # if the name is not in sdf:

    if not moldata.has_key(namekey):
        moldata[namekey] = ''

    # smiles

    inchiconv.SetInAndOutFormats('sdf', 'smi')
    mol = openbabel.OBMol()
    res = inchiconv.ReadString(mol, codecs.encode(sdf, 'utf-8'))
    if mol.Empty():
        root.warning('  -->  ERROR on sdf')
        raise Exception
    moldata['smiles'] = inchiconv.WriteString(mol).strip()

    return moldata


############################
# single compound operations
############################

def _update_single_compound(
    moldata,
    sdf,
    library,
    nameky,
    idkey,
    ):

    # sdf file

    s = SDFFile(sdffile=sdf)
    s.save()
    sdfid = s.id


def insert_single_compound(
    moldata,
    sdf,
    namekey,
    idkey,
    user,
    ):
    """ insert single compound into database """

    cid = moldata[idkey]
    name = moldata[namekey]
    if '\n' in name:
        name = name.split('\n')[0]

    # compound

    c = Compound(
        cid=cid,
        name=name,
        formula=moldata['formula'],
        weight=moldata['weight'],
        inchi=moldata['inchi'],
        smiles=moldata['smiles'],
        user=user,
        )

                    # sdf_file=s)

    c.save()
    c_id = c.id
    root.warning('  -->new compound inserted: c_id=%s, cid=%s' % (c_id,
                 cid))

    # sdf file

    s = SDFFile(sdffile=sdf, compound=c)
    s.save()
    sdfid = s.id

    return c.id


#####################################
# Physical Chemical Property - JOELib
#####################################

def gen_joelib_property(sdf):
    """run and parse the property output """

    # save the input in FS

    t = tempfile.NamedTemporaryFile(suffix='.sdf')
    t.write(codecs.encode(sdf, 'utf-8'))
    t.flush()

    # prepare the output file

    (f, out) = tempfile.mkstemp(suffix='.sdf')
    os.close(f)

    # convert

    cmd = \
        """JAVA_HOME=/opt/jre/ JOELIB2=/opt/JOELib2-alpha-20070303/ /opt/JOELib2-alpha-20070303/moleculeConversion.sh +d +h -iSDF -osdf "%s" "%s" > /dev/null""" \
        % (t.name, out)
    root.warning('  -->   running:%s' % cmd)
    if os.system(cmd) != 0:
        os.unlink(out)
        raise 'cannot run JOELib'

    # read and parse

    f = file(out)
    tags = get_sdf_tags(codecs.decode(f.read(), 'utf-8'))
    f.close()

    # clean

    os.unlink(out)
    return tags


######
# MISC
######

def update_mw(
    lib_name,
    lib_ver,
    input,
    rev=False,
    ):
    """goal: to update MW value with hydrogen added 
....   when calculating JOELib
....   'rev': in ChemMineV2, some libraries got compound ID and compound name switched, like 'Aurora'"""

    import datetime
    begin = datetime.datetime.now()
    print 'starts at: %s' % begin

    library = get_library(lib_name, lib_ver)
    mw = PropertyField.objects.get(name='MW')

    fp = file(input)
    line1 = fp.readline()

    count = 1
    for line in fp:
        (cid, weight) = line.strip().split('\t')
        try:
            if rev:
                c = Compound.objects.get(library=library, name=cid)
            else:
                c = Compound.objects.get(library=library, cid=cid)
        except Compound.DoesNotExist:
            print 'not found: line %s, cid=%s' % (count, cid)
            pass

        try:
            p = Property.objects.get(compound=c, field=mw)
            p.value = weight
            p.save()
        except Property.DoesNotExist:
            p = Property(field=mw, compound=c, value=weight)
            p.save()
            print 'new p for %s, line %s' % (cid, count)
        except:
            print '----->line %s, cid=%s' % (count, cid)
            pass

        count += 1

        # print "%s: %s -> %s", (cid, old, weight)

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
....   some sdf has blank 1st line, so we need to format SDF
....   by filling cid to 1st line in SDF"""

    fp = file(output, 'w')
    for sdf in sdfiterator.sdf_iter(sdffile):
        tagdict = get_sdf_tags(sdf)
        cid = tagdict[ID_tag]

        fp.write('%s\n' % cid)
        fp.write(sdf.split('\n', 1)[1].split('M  END')[0])
        fp.write('M  END\n')
        fp.write('''> <%s>
%s

''' % (ID_tag, cid))
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


