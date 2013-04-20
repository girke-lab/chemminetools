#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import tempfile
import re
import pybel

parser = \
    argparse.ArgumentParser(description='Compute openbabel properties and return table'
                            )
parser.add_argument('-o', '--outfile', help='output file',
                    required=True)
args = vars(parser.parse_args())

properties = [
    'abonds',
    'atoms',
    'bonds',
    'dbonds',
    'HBA1',
    'HBA2',
    'HBD',
    'logP',
    'MR',
    'MW',
    'nF',
    'sbonds',
    'tbonds',
    'TPSA',
    ]


def main():
    fa = sys.stdin.read()
    fa = fa.rstrip()
    of = open(args['outfile'], 'w')
    of.write('cid,')
    labels = ''
    for i in properties:
        labels = labels + i + ','
    labels = re.match(r"^(.*),", labels).group(1)
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
        desc = mol.calcdesc()
        for thisdesc in properties:
            info = info + str(desc[thisdesc]) + ','
        info = re.match(r"^(.*),", info).group(1)
        of.write(myid.strip() + ',' + info + '\n')
    of.close()
    os.unlink(inputTempName)


if __name__ == '__main__':
    main()
