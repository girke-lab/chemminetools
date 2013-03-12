#!/usr/bin/python
# make it possible to run as standalone program
import sys
sys.path.append('/srv/chemminetools')
sys.path.append('/srv/chemminetools/sdftools')	# allow tools.py to be imported
import os

from django.core.management import setup_environ
import chemminetools.settings
import argparse
setup_environ(chemminetools.settings)

from django.contrib.auth.models import User

from compounddb import tools
from compounddb.sdfiterator import sdf_iter
from csv import writer
from cStringIO import StringIO

parser = argparse.ArgumentParser(description='Add new compounds to a users workbench')
parser.add_argument('-f','--format', help='input format smiles, sdf, or pubchem', required=True)
parser.add_argument('-u','--user', help='numeric user id', required=True)
parser.add_argument('-o','--outfile', help='output file', required=True)
args = vars(parser.parse_args())

def main():
	# check if user account exists
	try:
		user = User.objects.get(id=user)	
	except:
		return False

		
	sdf = sys.stdin.read()
	output = app_processor(sdf)
	f = open(args['outfile'], 'w')
	f.write(output)
	f.close()

if __name__ == "__main__":
	main()
