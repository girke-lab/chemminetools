#!/usr/bin/python
# -*- coding: utf-8 -*-

"""define backend functions for similarity workbench here. At the frontend,
you pass "func=XXX" via ajax POST request, and XXX will be used to look up
a function in this module. This function must take SMILES as input and
return a dictionary. The returned value is encoded using JSON and passed to the
client side via ajax response. For client code example, see similarity.js and
search for "// option for ajaxified MCS calculation" or "// option for ajaxified
ap similarity calculation" to see how to define AJAX response handling. And
search for "ajaxSubmit" to see how to invoke the POST request.
"""

from tempfile import mkstemp
from subprocess import Popen, PIPE
import os
cur_dir = os.path.dirname(os.path.abspath(__file__))
ap_cmd = os.path.join(cur_dir, 'bin', 'compare cmp-cmp')
mcs_cmd = os.path.join(cur_dir, 'bin', 'calc_mcs.py')
from sdftools.moleculeformats import smiles_to_sdf, sdf_to_smiles


class ExternalProgramError(Exception):

    pass


def popen(cmd):
    """run the command <cmd>, and return the STDOUT. On non-zero return code, 
....raise an exception of type ExternalProgramError"""

    p = Popen(cmd, shell=True, stdout=PIPE, close_fds=True)
    sts = os.waitpid(p.pid, 0)[1]
    if sts != 0:
        raise ExternalProgramError()
    return p.stdout.read()


def smiles_to_sdf_input_adaptor(func):
    """a decorator to wrap a function that takes SDF file path input to make it
....   look like it takes SMILES input instead. We assume each input to be one
....   compound. We only convert positioned arguments and not keyword arguments
...."""

    def wrapped(*args, **kwargs):
        _args = []
        try:
            for i in args:

                # i is a SMILES string

                sdf = mkstemp()
                _args.append(sdf[1])
                sdf = os.fdopen(sdf[0], 'w')
                sdf.write(smiles_to_sdf(i))
                sdf.close()
        except:
            for i in _args:
                try:
                    os.unlink(i)
                except:
                    pass
            raise

        try:
            ret = func(*_args, **kwargs)
        finally:
            for i in _args:
                os.unlink(i)
        return ret

    return wrapped


@smiles_to_sdf_input_adaptor
def mcs(sdf1, sdf2):
    calc_cmd = "%s '%s' '%s'" % (mcs_cmd, sdf1, sdf2)
    ret = popen(calc_cmd)
    (status, mcs) = ret.split('\n', 1)
    if status == 'error':
        raise ExternalProgramError()
    (status_code, m_size, size1, size2) = status.split()
    m_size = int(m_size)
    min_size = min(int(size1), int(size2))
    max_size = max(int(size1), int(size2))
    sim_min = m_size * 1. / min_size
    sim_max = m_size * 1. / max_size
    sim_tanimoto = m_size * 1. / (min_size + max_size - m_size)

    mcs = sdf_to_smiles(mcs)
    (s, title) = mcs.split(None, 1)
    from hashlib import md5
    return dict(
        m_size=m_size,
        sim_min='%.4f' % sim_min,
        sim_max='%.4f' % sim_max,
        sim_tanimoto='%.4f' % sim_tanimoto,
        img='/similarity/renderer/' + s,
        md5=md5(mcs).hexdigest(),
        title=title,
        smiles=mcs,
        )


@smiles_to_sdf_input_adaptor
def ap(sdf1, sdf2):
    calc_cmd = "%s '%s' '%s'" % (ap_cmd, sdf1, sdf2)
    sim = popen(calc_cmd)
    sim = float(sim.strip())
    return dict(sim=sim)
