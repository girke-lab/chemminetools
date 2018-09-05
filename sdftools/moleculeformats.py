#!/usr/bin/python
# -*- coding: utf-8 -*-

"""format conversion"""
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
from builtins import str
import openbabel
import pybel
import re
import sys
import traceback


class InputError(Exception):

    pass



def smiles_to_sdf(smiles):
    try:
        mymol = pybel.readstring('smi', smiles)
    except:
        print("Unexpected error:", sys.exc_info())
        traceback.print_tb(sys.exc_info()[2])
        raise InputError
    try:
        mymol.draw(show=False, update=True)
    except:
        pass
    return mymol.write(format='sdf')

def sdf_to_smiles(sdf):

    # dos2unix

    sdf = sdf.split('\r\n')
    sdf = '\n'.join(sdf)
    sdf = sdf.split('\n')
    sdf[1] = ''
    sdf = '\n'.join(sdf)

    # remove unicode

    #if isinstance(sdf, str):
        #sdf = sdf.encode('ascii', 'ignore')

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('sdf', 'smi')

    mol = openbabel.OBMol()
    if not obConversion.ReadString(mol, sdf):
        raise InputError

    return obConversion.WriteString(mol)


def sdf_to_sdf(sdf):
    return smiles_to_sdf(sdf_to_smiles(sdf))


def smiles_to_smiles(smiles):
    return sdf_to_smiles(smiles_to_sdf(smiles))


def batch_sdf_to_smiles(sdfs):
    from compounddb.sdfiterator import sdf_iter
    from io import StringIO
    buf = ''
    err = 0
    for sdf in sdf_iter(StringIO(str(sdfs))):
        try:
            buf += sdf_to_smiles(sdf)
        except InputError:
            err += 1
    return (buf, err)


def batch_smiles_to_sdf(smiles):
    buf = ''
    err = 0
    for s in smiles.splitlines():
        try:
            buf += smiles_to_sdf(s)
        except InputError:
            err += 1
    return (buf, err)


def batch_sdf_to_sdf(sdfs):
    (smiles, _) = batch_sdf_to_smiles(sdfs)
    (sdf, __) = batch_smiles_to_sdf(smiles)
    return (sdf, _ + __)


def batch_smiles_to_smiles(smiles):
    (sdf, __) = batch_smiles_to_sdf(smiles)
    (smiles, _) = batch_sdf_to_smiles(sdf)
    return (smiles, _ + __)
