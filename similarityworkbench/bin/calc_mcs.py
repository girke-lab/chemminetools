#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import mcs

openbabel_import_warned = False


def convert_aromatic(graph, flexible=False):
    """detect aromatic bonds and set it to a unique type (type 4, or type 103).
....When flexible is set, it will be set to type 103 to allow it to match
....to both single- and double-bonds
...."""

    assert isinstance(graph, mcs.SDFGraph)
    assert graph.content

    if flexible:
        repl = '103'
    else:
        repl = '  4'

    try:
        from pybel import readstring
    except:
        global openbabel_import_warned
        if not openbabel_import_warned:
            openbabel_import_warned = True
            sys.stderr.write('WARNING: Cannot import OpenBabel. Will skip aromatic detection. This means aromatic bonds will be treated as single- or double-bonds as is hinted by the SDF.\n'
                             )
        return graph

    mol = readstring('sdf', graph.content)
    mol2_lines = mol.write('mol2').splitlines()

    indices = None
    for line in mol2_lines:
        if line.strip() == '@<TRIPOS>BOND':
            indices = []
        if indices is not None and len(line.split()) >= 3 \
            and line.split()[3] == 'ar':
            indices.append(int(line.split()[0]))

    sdf = graph.content.splitlines()
    n_atoms = int((sdf[3])[:3])
    for i in indices:
        row = 4 + n_atoms + i - 1
        sdf[row] = (sdf[row])[:6] + repl + (sdf[row])[9:]

    g = mcs.SDFGraph('\n'.join(sdf))
    return g


def sdf_mcs(sdffile1, sdffile2, aromatic=0):
    """calculate MCS from two input SDF files. 
....@param sdffile1 path to the first SDF file
....@param sdffile2 path to the second SDF file
....@param aromatic bonds processing:
........0: as is represented in SDF files
........1: treat as aromatic bonds
........2: allow it to match to both single and double bonds
...."""

    g1 = mcs.SDFGraph.from_file(sdffile1)
    g1_size = int((g1.content.splitlines()[3])[:3])
    g2 = mcs.SDFGraph.from_file(sdffile2)
    g2_size = int((g2.content.splitlines()[3])[:3])
    if aromatic == 0:
        m = mcs.mcs(g1, g2)
    elif aromatic == 1:

        # detect aromatic bonds first

        _g1 = convert_aromatic(g1)
        _g2 = convert_aromatic(g2)
        m = mcs.mcs(_g1, _g2)
    elif aromatic == 2:
        _g1 = convert_aromatic(g1, flexible=True)
        _g2 = convert_aromatic(g2, flexible=True)
        m = mcs.mcs(_g1, _g2)
    else:
        raise 'wrong aromatic processing hint'
    m_size = len(m)
    sub = g1.sub(m.keys())
    return (g1_size, g2_size, m_size, sub)


if __name__ == '__main__':
    if len(sys.argv) != 3 and len(sys.argv) != 4:
        sys.stderr.write('Usage: calc_mcs 1.sdf 2.sdf [1]\n')
        sys.stderr.write('''The optional argument at the end will enable flexiable matching of aromatic bonds

''')
        sys.exit(1)
    try:
        (sdffile1, sdffile2) = sys.argv[1:3]
        flex = 1
        if len(sys.argv) == 4 and sys.argv[3] == '1':
            flex = 2
        (g1_size, g2_size, m_size, sub) = sdf_mcs(sdffile1, sdffile2,
                aromatic=flex)
        print 'ok %d %d %d' % (m_size, g1_size, g2_size)
        print sub
    except:
        print 'error'
        raise
