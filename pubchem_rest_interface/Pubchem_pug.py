#!/usr/bin/python
# -*- coding: utf-8 -*-


from builtins import str
import requests
from requests.compat import quote



def DownloadCIDs(cids):

    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+(",".join([str(x) for x in cids]))+"/SDF"
    req = requests.get(url)
    if req.status_code == requests.codes.ok:
        return req.text
    else:
        print("failed to download sdf from pubchem. status: {}, result: {}".format(req.status_code,req.text))
        
def SimilaritySearch(smilesQuery,threshold,maxResults):
    #print("SimilaritySearch. query: {}, threshold: {}, maxResults: {}".format(smilesQuery,threshold,maxResults))
    #url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/"+
     #       quote(smilesQuery.rstrip(),safe='')+"/cids/TXT?Threshold="+str(threshold)+"&MaxRecords="+str(maxResults))
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/"+
            "cids/TXT?Threshold="+str(threshold)+"&MaxRecords="+str(maxResults))

    print("SimilaritySearch. url: "+url)
    req = requests.post(url,data={'smiles':smilesQuery})
    if req.status_code == requests.codes.ok:
        return req.text
    else:
        print("failed to run similarity search on pubchem. status: {}, result: {}".format(req.status_code,req.text))
        raise Exception("failed to run similarity search on pubchem. status: {}, result: {}".format(req.status_code,req.text))

     
if __name__ == '__main__':
    DownloadCIDs()
