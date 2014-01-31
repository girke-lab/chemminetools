#!/usr/bin/python
# -*- coding: utf-8 -*-

# make it possible to run as standalone program
import sys
sys.path.append('/srv/chemminetools')
sys.path.append('/srv/chemminetools/sdftools')  # allow tools.py to be imported
from django.core.management import setup_environ
import chemminetools.settings
import argparse
setup_environ(chemminetools.settings)
from django.contrib.auth.models import User
from django.conf import settings

# load functions used
from pubchem_soap_interface.SimilaritySearch import SimilaritySearch
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
        id_list = SimilaritySearch(query, similarity, compounds)
        if id_list == 'error':
            raise Exception
        id_list = [str(cid) for cid in id_list]
        output = '\n'.join(id_list)
    except:
        output = ''

    f = open(args['outfile'], 'w')
    f.write(output)
    f.close()

if __name__ == '__main__':
    main()
