#!/usr/bin/python
# -*- coding: utf-8 -*-

# make it possible to run as standalone program

from __future__ import unicode_literals
from builtins import str
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
args = vars(parser.parse_args())

def addMyCompounds(sdf, user,tags):
    sdffile = ''
    counter = 0
    linecounter = 0
    namekey = 'PUBCHEM_IUPAC_NAME'
    message = 'ERROR: bad input data.'
    added_ids = []
    initial_compound_count = Compound.objects.filter(user=user).count()

    try:
        if not isinstance(sdf, str):
            sdf = str(sdf, 'utf-8')
        sdf = sdf.split('\n')
        for line in sdf:
            linecounter += 1
            if linecounter > MAX_SDF_LENGTH:
                message = 'ERROR: an input sdf exceeds ' \
                    + str(MAX_SDF_LENGTH) + ' lines.'
                raise Exception
            if linecounter == 1:

                # clean up cid with regexes

                line = re.match(r"^\W*(.*?)\W*$", line).group(1)
                if line == '':
                    line = 'unspecified_' \
                        + ''.join(random.sample(string.digits, 6))
                line = re.sub(r"[^a-zA-Z_0-9-]", '_', line, count=0)

                # loop adding postfix numbers to the cid until we find a unique cid in the database

                appendNumber = 1
                oldCid = line
                while len(Compound.objects.filter(cid=line, user=user)) \
                    > 0:
                    appendNumber += 1
                    line = str(oldCid) + '_' + str(appendNumber)
            sdffile += line
            sdffile += '\n'
            if line.startswith('$$$$'):
                try:
                    moldata = parse_annotation(sdffile, namekey)
                except:
                    print("error while processing sdf: ")
                    print("Unexpected error:", sys.exc_info())
                    traceback.print_tb(sys.exc_info()[2])
                    message = 'ERROR: invalid input format.'
                    raise Exception
                counter += 1
                linecounter = 0
                if initial_compound_count + counter > MAX_COMPOUND_LIMIT:
                    message = 'ERROR: upload exceeds ' \
                        + str(MAX_COMPOUND_LIMIT) + ' compounds.'
                    raise Exception
                try:
                    newid = insert_single_compound(moldata, sdffile,
                            namekey, 'id', user,tags)
                except:
                    print("error while inserting compound: ")
                    print("Unexpected error:", sys.exc_info())
                    traceback.print_tb(sys.exc_info()[2])
                    message = \
                        'ERROR: Database error, possibly excessively large compound?'
                    raise Exception
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
        print("exception: "+message)
        return message


def main():
    sdf = sys.stdin.read()
    user = User.objects.get(id=args['user'])
    tags = Tag.objects.filter(name__in=args['tags'].split(','))

    output = addMyCompounds(sdf, user,tags)

    f = open(args['outfile'], 'w')
    f.write(output)
    f.close()


if __name__ == '__main__':
    main()

