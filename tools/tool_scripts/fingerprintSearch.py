#!/usr/bin/python
# -*- coding: utf-8 -*-

# make it possible to run as standalone program
from builtins import str
import sys
sys.path.append('/srv/chemminetools')
sys.path.append('/srv/chemminetools/sdftools')  # allow tools.py to be imported

import os 
import django
os.environ.setdefault("DJANGO_SETTINGS_MODULE","chemminetools.settings")
django.setup()
#from django.core.management import setup_environ

#import chemminetools.settings
#setup_environ(chemminetools.settings)
import argparse
from django.contrib.auth.models import User
from django.conf import settings

# load functions used
from pubchem_rest_interface.Pubchem_pug import SimilaritySearch
from sdftools.moleculeformats import sdf_to_smiles

# parse command line arguments
parser = argparse.ArgumentParser(description='perform pubchem fingerprint search')
parser.add_argument('-o', '--outfile', help='output file',required=True)
parser.add_argument('-s', '--similarity', help='similarity cutoff',required=True)
parser.add_argument('-c', '--compounds', help='max compounds returned',required=True)
args = vars(parser.parse_args())

def main():
    sdf = sys.stdin.read()
    similarity = args['similarity']
    similarity = int(float(similarity) * 100)
    compounds = int(args['compounds'])
    query = sdf_to_smiles(sdf)

    try:
        output = SimilaritySearch(query, similarity, compounds)
        #print("id_list: {}".format(id_list))
       # if id_list == 'error':
       #     raise Exception
       # id_list = [str(cid) for cid in id_list]
       # output = '\n'.join(id_list)
    except:
        output = ''

    f = open(args['outfile'], 'w')
    f.write(output)
    f.close()

if __name__ == '__main__':
    main()
