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

from compounddb import tools
from compounddb.sdfiterator import sdf_iter
from csv import writer
from cStringIO import StringIO

# giant messy fields table for making sense of joelib results:
fields = [(u'joelib2.feature.types.MolecularWeight', u'MW'), (u'joelib2.feature.types.LogP', u'LGP'), (u'joelib2.feature.types.count.HBA1', u'HA1'), (u'joelib2.feature.types.count.HBA2', u'HA2'), (u'joelib2.feature.types.count.HBD1', u'HD1'), (u'joelib2.feature.types.count.HBD2', u'HD2'), (u'joelib2.feature.types.count.AcidicGroups', u'ACG'), (u'joelib2.feature.types.count.AliphaticOHGroups', u'AOH'), (u'joelib2.feature.types.count.AromaticBonds', u'ARB'), (u'joelib2.feature.types.count.BasicGroups', u'BAG'), (u'joelib2.feature.types.FractionRotatableBonds', u'FRB'), (u'joelib2.feature.types.RotatableBonds', u'ROT'), (u'joelib2.feature.types.HeavyBonds', u'HB'), (u'joelib2.feature.types.HeteroCycles', u'HCY'), (u'joelib2.feature.types.HydrophobicGroups', u'HPG'), (u'joelib2.feature.types.MolarRefractivity', u'MR'), (u'joelib2.feature.types.count.NumberOfAtoms', u'AT'), (u'joelib2.feature.types.count.NumberOfHal', u'HAL'), (u'joelib2.feature.types.count.NumberOfB', u'B'), (u'joelib2.feature.types.count.NumberOfBr', u'Br'), (u'joelib2.feature.types.count.NumberOfBonds', u'BON'), (u'joelib2.feature.types.count.NumberOfCl', u'Cl'), (u'joelib2.feature.types.count.NumberOfI', u'I'), (u'joelib2.feature.types.count.NumberOfF', u'F'), (u'joelib2.feature.types.count.NumberOfN', u'N'), (u'joelib2.feature.types.count.NumberOfO', u'O'), (u'joelib2.feature.types.count.NumberOfP', u'P'), (u'joelib2.feature.types.count.NumberOfS', u'S'), (u'joelib2.feature.types.count.NO2Groups', u'NO2'), (u'joelib2.feature.types.count.OSOGroups', u'OSO'), (u'joelib2.feature.types.PolarSurfaceArea', u'POL'), (u'joelib2.feature.types.count.SOGroups', u'SO'), (u'joelib2.feature.types.count.SO2Groups', u'SO2'), (u'joelib2.feature.types.GeometricalDiameter', u'GDI'), (u'joelib2.feature.types.GeometricalRadius', u'GRA'), (u'joelib2.feature.types.Geometrical diameter', u'GSH'), (u'joelib2.feature.types.TopologicalDiameter', u'TDI'), (u'joelib2.feature.types.TopologicalRadius', u'TRA'), (u'joelib2.feature.types.KierShape1', u'KS1'), (u'joelib2.feature.types.KierShape2', u'KS2'), (u'joelib2.feature.types.ZagrebIndex1', u'ZI1'), (u'NA', u'IC'), (u'NA', u'BIC'), (u'NA', u'CIC'), (u'NA', u'SIC'), (u'NA', u'IAC-Total'), (u'NA', u'V-ADJ-mag'), (u'NA', u'V-DIST-mag'), (u'NA', u'E-ADJ-mag'), (u'NA', u'E-DIST-mag'), (u'NA', u'Rotlbonds'), (u'NA', u'Hbond acceptor'), (u'NA', u'Hbond donor'), (u'NA', u'AlogP98'), (u'NA', u'JX'), (u'NA', u'Kappa-1'), (u'NA', u'Kappa-2'), (u'NA', u'Kappa-3'), (u'NA', u'Kappa-1-AM'), (u'NA', u'Kappa-2-AM'), (u'NA', u'Kappa-3-AM'), (u'NA', u'PHI'), (u'NA', u'SC-0'), (u'NA', u'SC-1'), (u'NA', u'SC-2'), (u'NA', u'SC-3_P'), (u'NA', u'SC-3_C'), (u'NA', u'SC-3_CH'), (u'NA', u'CHI-0'), (u'NA', u'CHI-1'), (u'NA', u'CHI-2'), (u'NA', u'CHI-3_P'), (u'NA', u'CHI-3_C'), (u'NA', u'CHI-3_CH'), (u'NA', u'CHI-V-0'), (u'NA', u'CHI-V-1'), (u'NA', u'CHI-V-2'), (u'NA', u'CHI-V-3_P'), (u'NA', u'CHI-V-3_C'), (u'NA', u'CHI-V-3_CH'), (u'NA', u'Wiener'), (u'NA', u'Zagreb'), (u'joelib2.feature.types.Geometrical diameter', u'GSC')]

parser = argparse.ArgumentParser(description='Compute joelib properties and return table')
parser.add_argument('-o','--outfile', help='output file', required=True)
args = vars(parser.parse_args())

def float_format(x):
	try:
		x = float(x)
		x = '%.6g' % x
	except: pass
	return x

def table_format(fields_list, joelibs_list):
	"""given a list of fieldsets, find the common fieldset; then format the list
	of joelibs into a table, whose columns correpond to the common fieldset"""
	common = reduce(lambda x,y:x.intersection(set(y)), fields_list,
			set(fields_list[0]))
	common = list(common)
	joelibs_list = [
		[float_format(j[i]) for i in common]
		for j in joelibs_list]
	columns = [i[1] for i in fields if i[0] in common]
	return columns, joelibs_list

# from sample import smiles_to_sdf
# @smiles_to_sdf
def app_processor(x):
	fields_list = []
	joelibs_list = []
	cids = []
	for sdf in sdf_iter(iter([i + '\n'
				for i in x.splitlines()])):
		cids.append(sdf.splitlines()[0].strip())
		j = tools.gen_joelib_property(sdf)
		available_fields = [i[0] for i in fields if i[0] in j]
		fields_list.append(available_fields)
		joelibs_list.append(j)
	columns, joelibs_list = table_format(fields_list, joelibs_list)
	storage = StringIO()
	w = writer(storage)
	columns.insert(0, 'ID')
	w.writerow(columns)
	for i, joelib in enumerate(joelibs_list):
		joelib.insert(0, cids[i])
		w.writerow(joelib)
	return storage.getvalue()

def main():
	sdf = sys.stdin.read()
	output = app_processor(sdf)
	f = open(args['outfile'], 'w')
	f.write(output)
	f.close()

if __name__ == "__main__":
	main()
