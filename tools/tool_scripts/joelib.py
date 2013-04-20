#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import pickle
import urllib2
import urllib
from StringIO import StringIO
import tempfile
import re
import pybel

# dir = '/home/ycao/xmlrpc/'

columns = [
    'Molecular_weight',
    'LogP',
    'Number_of_HBA_1',
    'Number_of_HBA_2',
    'Number_of_HBD_1',
    'Number_of_HBD_2',
    'Number_of_acidic_groups',
    'Number_of_aliphatic_OH_groups',
    'Number_of_basic_groups',
    'Fraction_of_rotatable_bonds',
    'Number_of_heavy_bonds',
    'Number_of_heterocycles',
    'Number_of_hydrophobic_groups',
    'MolarRefractivity',
    'Number_of_atoms',
    'Number_of_halogen_atoms',
    'Number_of_B_atoms',
    'Number_of_Br_atoms',
    'Number_of_Cl_atoms',
    'Number_of_I_atoms',
    'Number_of_F_atoms',
    'Number_of_N_atoms',
    'Number_of_O_atoms',
    'Number_of_P_atoms',
    'Number_of_S_atoms',
    'Number_of_bonds',
    'Number_of_NO2_groups',
    'Number_of_SO_groups',
    'Number_of_OSO_groups',
    'Number_of_SO2_groups',
    'PolarSurfaceArea',
    'Geometrical_diameter',
    'Geometrical_radius',
    'Geometrical_shape_coefficient',
    'Kier_shape_1',
    'Kier_shape_2',
    'Zagreb_group_index_1',
    'Zagreb_group_index_2',
    ]

parser = \
    argparse.ArgumentParser(description='Compute joelib properties and return table'
                            )
parser.add_argument('-o', '--outfile', help='output file',
                    required=True)
args = vars(parser.parse_args())


def getinfo(f):
    flag = 1

    # skip to annotation

    while flag == 1:
        line = f.readline()
        if line == '':
            return ''
        if line[0:6] == 'M  END':
            flag = 0
    flag = 1
    data = dict()
    tag = ''
    tagdata = ''
    while flag == 1:
        line = f.readline()
        if line == '':
            break
        if line == None:
            break
        if re.match('>\s+<[A-Za-z0-9_]+>', line) is not None:

            # save the previous tag if it is ther

            if tag != '':
                if tag in columns:
                    data[tag] = tagdata.strip()
                    tagdata = ''
            tag = line.replace('>', '')
            tag = tag.replace('<', '')
            tag.strip()
            tag = re.sub('\s+', '', tag)
            tagdata = ''
        elif line[0:5] == '$$$$':
            if tag != '':
                if tag in columns:
                    data[tag] = tagdata.strip()
            flag = 0
            break
        else:
            tagdata = tagdata + line
    r = ''
    for i in columns:
        r = r + str(data[i].strip()) + ','
    r = re.match(r"^(.*),", r).group(1)
    return r


def main():
    fa = sys.stdin.read()
    fa = fa.rstrip()

    # fa = open('assofile')

    of = open(args['outfile'], 'w')
    of.write('cid,')
    labels = ''
    for i in columns:
        labels = labels + i + ','
    labels = re.match(r"^(.*),", labels).group(1)
    of.write(labels + '\n')
    inputTemp = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False)
    inputTemp.write(fa)
    inputTempName = inputTemp.name
    inputTemp.close()
    for mol in pybel.readfile('sdf', inputTempName):
        myid = mol.title
        smi = mol.write(format='smi')
        sdf = pybel.readstring('smi', smi)
        sdf.addh()
        sdf.make3D()
        sdf = sdf.write(format='sdf')
        sdfTemp = tempfile.NamedTemporaryFile(suffix='.sdf',
                delete=False)
        sdfTemp.write(sdf)
        sdfName = sdfTemp.name
        sdfTemp.close()
        fh = tempfile.NamedTemporaryFile(suffix='.sdf')

        # for every compound run the conversion

        cmd = \
            'JAVA_HOME=/opt/jre/ JOELIBDIR=/opt/JOELib/ CLASSPATH=/opt/JOELib/build /opt/JOELib/convert.sh +d -iSDF -oSDF ' \
            + sdfName + ' ' + fh.name + ' > /dev/null'
        os.system(cmd)
        os.unlink(sdfName)
        info = getinfo(fh)
        of.write(myid.strip() + ',' + info + '\n')
        fh.close()
    of.close()
    os.unlink(inputTempName)


if __name__ == '__main__':
    main()
