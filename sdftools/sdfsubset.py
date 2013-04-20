#!/usr/bin/python
# -*- coding: utf-8 -*-

"""subset a big sdf"""

import sys
if len(sys.argv) != 4:
    sys.stderr.write('Usage: %s input.sdf indices.file output[.sdf]\n'
                     % sys.argv[0])
    sys.stderr.write('Indices are 1-based\n')
    sys.stderr.write('''If second arguments ends in .sdf, a single file will be created. Otherwise, a
folder will be created and each SDF will be placed in a separate file
''')

    sys.exit(1)
inp = sys.argv[1]
indices = set([int(i) for i in file(sys.argv[2]).readlines()])
out = sys.argv[3]
if out.endswith('.sdf'):
    outf = file(out, 'w')
    mode = 'single-file'
else:
    import shutil
    import os
    shutil.os.mkdir(out)
    mode = 'directory'

from sdfiterator import sdf_iter
iter = sdf_iter(inp)

cntr = 0
to_capture = len(indices)
for i in iter:
    cntr += 1
    if not cntr in indices:
        continue
    if mode == 'single-file':
        outf.write(i)
    elif mode == 'directory':

        # name = i.split('\n')[0].strip()

        name = str(cntr)
        file(os.path.join(out, name), 'w').write(i)
    to_capture -= 1
    if to_capture == 0:
        break
