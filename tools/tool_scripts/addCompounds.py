#!/usr/bin/env python
# -*- coding: utf-8 -*-

# make it possible to run as standalone program

import sys
import re
import random
import string
import traceback
sys.path.append('/srv/chemminetools')
sys.path.append('/srv/chemminetools/sdftools')  # allow tools.py to be imported

import os
import django
#from chemmminetools import models


os.environ.setdefault("DJANGO_SETTINGS_MODULE","chemminetools.settings")
django.setup()

#import chemminetools.settings
#from django.core.management import setup_environ
#setup_environ(chemminetools.settings)

import argparse
from django.contrib.auth.models import User
from django.conf import settings
from compounddb.tools import parse_annotation, insert_single_compound
from compounddb.models import Compound, SDFFile, Tag
MAX_COMPOUND_LIMIT = settings.MAX_COMPOUND_LIMIT
MAX_SDF_LENGTH = settings.MAX_SDF_LENGTH

parser = \
    argparse.ArgumentParser(description='Add new compounds to a users workbench'
                            )
parser.add_argument('-u', '--user', help='numeric user id',
                    required=True)
parser.add_argument('-t', '--tags', help='comma seperated string of tags',
                    required=False)
parser.add_argument('-o', '--outfile', help='output file',
                    required=True)
parser.add_argument('-d', '--deduplicate', help='skip inserting any compounds that already exist for given user',
                    action='store_true',required=False)
args = vars(parser.parse_args())

def addMyCompounds(sdf, user,tags,dedup):
    sdffile = ''
    counter = 0
    linecounter = 0
    namekeys = ['NAME', 'PUBCHEM_IUPAC_NAME']
    added_ids = []
    initial_compound_count = Compound.objects.filter(user=user).count()

    try:
        if not isinstance(sdf, str):
            sdf = str(sdf, 'utf-8')
        sdf = sdf.split('\n')
        for line in sdf:
            linecounter += 1
            if linecounter > MAX_SDF_LENGTH:
                raise Exception( 'ERROR: an input sdf exceeds ' \
                    + str(MAX_SDF_LENGTH) + ' lines.')
            if linecounter == 1:

                # clean up cid with regexes

                line = re.match(r"^\W*(.*?)\W*$", line).group(1)
                if line == '':
                    line = 'unspecified_' \
                        + ''.join(random.sample(string.digits, 6))
                line = re.sub(r"[^a-zA-Z_0-9-]", '_', line, count=0)

                # loop adding postfix numbers to the cid until we find a unique cid in the database
                # disabling unique CIDs, since adding the postfix will cause
                # problems when searching with CIDs.
                #appendNumber = 1
                #oldCid = line
                #while len(Compound.objects.filter(cid=line, user=user)) \
                #    > 0:
                #    appendNumber += 1
                #    line = str(oldCid) + '_' + str(appendNumber)
            sdffile += line
            sdffile += '\n'
            if line.startswith('$$$$'):
                try:
                    moldata = parse_annotation(sdffile, namekeys)
                except Exception as e:
                    print("error while processing sdf: ")
                    print("Unexpected error:", sys.exc_info())
                    traceback.print_tb(sys.exc_info()[2])
                    raise Exception( 'ERROR: invalid input format.') from e
                counter += 1
                linecounter = 0
                if initial_compound_count + counter > MAX_COMPOUND_LIMIT:
                    raise Exception('ERROR: upload exceeds ' \
                        + str(MAX_COMPOUND_LIMIT) + ' compounds.')
                try:
                    newid = insert_single_compound(moldata, sdffile,
                            'id', user, tags,dedup)
                except Exception as e:
                    print("error while inserting compound: ")
                    print("Unexpected error:", sys.exc_info())
                    traceback.print_tb(sys.exc_info()[2])
                    raise Exception( \
                        'ERROR: failed to insert compound with id '+str(moldata['id'])) from e
                added_ids.append(newid)
                sdffile = ''
        if counter > 0:
            return 'Success: Added ' + str(counter) + ' compounds.'
        else:
            return 'ERROR: No valid input found.'
    except Exception as e:
        print("error while reading sdf: ")
        print("Unexpected error:", sys.exc_info())
        traceback.print_tb(sys.exc_info()[2])
        for id in added_ids:
            try:
                Compound.objects.get(id=id).delete()
            except:
                pass
        message = str(sys.exc_info()[1])
        if e.__cause__ != None:
            message = message + " Caused by: "+str(e.__cause__)
        return message


def main():
    sdf = sys.stdin.read()
    user = User.objects.get(id=args['user'])
    tags = []
    if 'tags' in args and args['tags'] is not None:
        tags = Tag.objects.filter(user=user,name__in=args['tags'].split(','))

    output = addMyCompounds(sdf, user,tags, args['deduplicate'])

    f = open(args['outfile'], 'w')
    f.write(output)
    f.close()


if __name__ == '__main__':
    main()

