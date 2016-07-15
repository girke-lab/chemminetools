#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import tempfile
import re
import pybel
import openbabel as ob

parser = \
    argparse.ArgumentParser(description='Compute InChI string and return table'
                            )
parser.add_argument('-o', '--outfile', help='output file',
                    required=True)
args = vars(parser.parse_args())

def main():
    fa = sys.stdin.read()
    fa = fa.rstrip()
    of = open(args['outfile'], 'w')
    of.write('cid\t')
    labels = 'InChI\tInChIkey'
    of.write(labels + '\n')
    inputTemp = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False)
    inputTemp.write(fa)
    inputTempName = inputTemp.name
    inputTemp.close()
    for mol in pybel.readfile('sdf', inputTempName):
        info = ''
        myid = mol.title
        mol.addh()
        mol.make3D()
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("sdf", "inchi")
        inchi = conv.WriteString(mol.OBMol)
        conv.SetInAndOutFormats("sdf", "inchikey")
        inchikey = conv.WriteString(mol.OBMol)
        info = re.match(r"^InChI=(.*)\n", inchi).group(1)
        of.write(myid.strip() + '\t' + info + '\t' + inchikey + '\n')
    of.close()
    os.unlink(inputTempName)


if __name__ == '__main__':
    main()
