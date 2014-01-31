#!/usr/bin/python
# -*- coding: utf-8 -*-

# make it possible to run as standalone program
import sys
import string
import re
sys.path.append('/srv/chemminetools')
from django.core.management import setup_environ
import chemminetools.settings
setup_environ(chemminetools.settings)
from django.contrib.auth.models import User
from django.conf import settings
from django.contrib import messages

# load functions used
from pubchem_soap_interface.DownloadCIDs import DownloadCIDs 

def main():
    cids = sys.stdin.read()
    cids = cids.split()
    filteredCIDs = []
    for cid in cids[:]:
        match = re.search("(\d{1,200})", cid)
        if match:
            filteredCIDs.append(int(match.group(1)))

    if len(filteredCIDs) > 0:
        try:
            sdf = DownloadCIDs(filteredCIDs)
        except:
            messages.error(request, 'Invalid CIDs or no response from PubChem!')
            sdf = None
    else:
        messages.error(request, 'Error: No valid CIDs entered!')
        sdf = None
    print sdf

if __name__ == '__main__':
    main()
