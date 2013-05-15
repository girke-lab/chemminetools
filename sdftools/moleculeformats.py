#!/usr/bin/python
# -*- coding: utf-8 -*-

"""format conversion"""

import openbabel
import pybel
import re


class InputError(Exception):

    pass


def unicode_wrapper(f):

    def wrapped(inp, *args, **kwargs):
        if isinstance(inp, unicode):
            inp = unicode.encode(inp, 'iso8859-1')
        return f(inp, *args, **kwargs)

    return wrapped


@unicode_wrapper
def smiles_to_sdf(smiles):
    try:
        mymol = pybel.readstring('smi', str(smiles))
    except:
        raise InputError
    try:
        mymol.draw(show=False, update=True)
    except:
        pass
    return mymol.write(format='sdf')

@unicode_wrapper
def sdf_to_smiles(sdf):

    # dos2unix

    sdf = sdf.split('\r\n')
    sdf = '\n'.join(sdf)
    sdf = sdf.split('\n')
    sdf[1] = ''
    sdf = '\n'.join(sdf)

    # remove unicode

    if isinstance(sdf, unicode):
        sdf = sdf.encode('ascii', 'ignore')

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('sdf', 'smi')

    mol = openbabel.OBMol()
    if not obConversion.ReadString(mol, sdf):
        raise InputError

    return obConversion.WriteString(mol)


@unicode_wrapper
def sdf_to_sdf(sdf):
    return smiles_to_sdf(sdf_to_smiles(sdf))


@unicode_wrapper
def smiles_to_smiles(smiles):
    return sdf_to_smiles(smiles_to_sdf(smiles))


@unicode_wrapper
def batch_sdf_to_smiles(sdfs):
    from sdfiterator import sdf_iter
    from cStringIO import StringIO
    buf = ''
    err = 0
    for sdf in sdf_iter(StringIO(sdfs)):
        try:
            buf += sdf_to_smiles(sdf)
        except InputError:
            err += 1
    return (buf, err)


@unicode_wrapper
def batch_smiles_to_sdf(smiles):
    buf = ''
    err = 0
    for s in smiles.splitlines():
        try:
            buf += smiles_to_sdf(s)
        except InputError:
            err += 1
    return (buf, err)


@unicode_wrapper
def batch_sdf_to_sdf(sdfs):
    (smiles, _) = batch_sdf_to_smiles(sdfs)
    (sdf, __) = batch_smiles_to_sdf(smiles)
    return (sdf, _ + __)


@unicode_wrapper
def batch_smiles_to_smiles(smiles):
    (sdf, __) = batch_smiles_to_sdf(smiles)
    (smiles, _) = batch_sdf_to_smiles(sdf)
    return (smiles, _ + __)
