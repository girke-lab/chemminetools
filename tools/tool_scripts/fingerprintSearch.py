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

from compounddb.sdfiterator import sdf_iter
from io import StringIO
import re

# parse command line arguments
parser = argparse.ArgumentParser(description='perform pubchem fingerprint search')
parser.add_argument('-o', '--outfile', help='output file',required=True)
parser.add_argument('-s', '--similarity', help='similarity cutoff',required=True)
parser.add_argument('-c', '--compounds', help='max compounds returned',required=True)
args = vars(parser.parse_args())

def main():
    sdfs = sys.stdin.read()
    similarity = args['similarity']
    similarity = int(float(similarity) * 100)
    compounds = int(args['compounds'])

    with open(args['outfile'], 'w') as f:
        for sdf in sdf_iter(StringIO(str(sdfs))):
            query_id =""
            parts = re.split('\s+',sdf_to_smiles(sdf).strip())
            smile_query = parts[0]
            if len(parts) > 1: 
                query_id =  parts[1]
            output = SimilaritySearch(smile_query, similarity, compounds)
            for match in output.strip().split('\n'):
                f.write(query_id +" "+match+"\n")



    #try:
    #    output = SimilaritySearch(query, similarity, compounds)
    #    #print("id_list: {}".format(id_list))
    #   # if id_list == 'error':
    #   #     raise Exception
    #   # id_list = [str(cid) for cid in id_list]
    #   # output = '\n'.join(id_list)
    #except:
    #    output = ''


if __name__ == '__main__':
    main()
