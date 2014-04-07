#!/usr/bin/python
# -*- coding: utf-8 -*-

# make it possible to run as standalone program

import sys
import re
import random
import string
sys.path.append('/srv/chemminetools')
sys.path.append('/srv/chemminetools/sdftools')  # allow tools.py to be imported
from django.core.management import setup_environ
import chemminetools.settings
import argparse
setup_environ(chemminetools.settings)
from django.contrib.auth.models import User
from django.conf import settings
from compounddb.tools import parse_annotation, insert_single_compound
from compounddb.models import Compound, SDFFile
MAX_COMPOUND_LIMIT = settings.MAX_COMPOUND_LIMIT
MAX_SDF_LENGTH = settings.MAX_SDF_LENGTH

parser = \
    argparse.ArgumentParser(description='Add new compounds to a users workbench'
                            )
parser.add_argument('-u', '--user', help='numeric user id',
                    required=True)
parser.add_argument('-o', '--outfile', help='output file',
                    required=True)
args = vars(parser.parse_args())


def addMyCompounds(sdf, user):
    sdffile = u''
    counter = 0
    linecounter = 0
    namekey = 'PUBCHEM_IUPAC_NAME'
    message = 'ERROR: bad input data.'
    added_ids = []
    try:
        if not isinstance(sdf, unicode):
            sdf = unicode(sdf, 'utf-8')
        sdf = sdf.encode('ascii', 'ignore')
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
                    line = oldCid + '_' + str(appendNumber)
            sdffile += line
            sdffile += '\n'
            if line.startswith('$$$$'):
                try:
                    moldata = parse_annotation(sdffile, namekey)
                except:
                    message = 'ERROR: invalid input format.'
                    raise Exception
                counter += 1
                linecounter = 0
                if counter > MAX_COMPOUND_LIMIT:
                    message = 'ERROR: upload exceeds ' \
                        + str(MAX_COMPOUND_LIMIT) + ' compounds.'
                    raise Exception
                try:
                    newid = insert_single_compound(moldata, sdffile,
                            namekey, 'id', user)
                except:
                    message = \
                        'ERROR: Database error, possibly excessively large compound?'
                    raise Exception
                added_ids.append(newid)
                sdffile = u''
        if counter > 0:
            return 'Success: Added ' + str(counter) + ' compounds.'
        else:
            return 'ERROR: No valid input found.'
    except:
        for id in added_ids:
            try:
                Compound.objects.get(id=id).delete()
            except:
                pass
        return message


def main():
    sdf = sys.stdin.read()
    user = User.objects.get(id=args['user'])

    output = addMyCompounds(sdf, user)

    f = open(args['outfile'], 'w')
    f.write(output)
    f.close()


if __name__ == '__main__':
    main()

