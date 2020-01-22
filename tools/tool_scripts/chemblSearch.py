#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings
Settings.Instance().CACHING = False
#Settings.Instance().CACHE_NAME = "/tmp/chembl_webresource_cache"
from io import StringIO
import re

sys.path.append('/srv/chemminetools')
sys.path.append('/srv/chemminetools/sdftools')  # allow tools.py to be imported
from sdftools.moleculeformats import sdf_to_smiles
from compounddb.sdfiterator import sdf_iter

parser = argparse.ArgumentParser(
            description='Compute similarity to CHEMBL compounds using REST API')
parser.add_argument('-o', '--outfile', help='output file', required=True)
parser.add_argument('-s', '--similarity', help='similarity cutoff', required=True)
parser.add_argument('-c', '--compounds', help='max number of matches returned',required=True)
args = vars(parser.parse_args())



def chembl_similarity(smile,threshold,numResults):
    similarity_query = new_client.similarity
    print("smile query: "+smile)
    print("similarity: "+str(threshold))
    results = similarity_query.filter(smiles=smile, similarity=threshold). \
                only(['molecule_chembl_id','similarity'])
    #print(str(results))
    if results: 
        return [ (results[i]['molecule_chembl_id'],results[i]['similarity']) for i in range(numResults) if i < len(results)]
    else:
        return None

def main():

    sdfs= sys.stdin.read()

    similarity = args['similarity']
    similarity = int(float(similarity) * 100)
    numResults = int(args['compounds'])

    with open(args['outfile'], 'w') as f:
        for sdf in sdf_iter(StringIO(str(sdfs))):
            query_id =""
            parts = re.split('\s+',sdf_to_smiles(sdf).strip())
            smile_query = parts[0]
            if len(parts) > 1: 
                query_id =  parts[1]
            matches = chembl_similarity(smile_query, similarity, numResults)
            if matches != None:
                for match in matches:
                    print("match: "+str(match))
                    f.write(query_id+" "+match[0]+" "+match[1]+"\n")



if __name__ == '__main__':
    main()
