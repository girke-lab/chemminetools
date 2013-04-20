#!/usr/bin/python
# -*- coding: utf-8 -*-

"""make sdf smaller by keeping only the MOL. Also, compounds with bonds
fewer than the limit will be skipped"""

limit = 0

import sys
inp = sys.argv[1]
if sys.argv[2] == '-':
    outp = sys.stdout
else:
    outp = file(sys.argv[2], 'w')
if len(sys.argv) == 4:
    limit = int(sys.argv[3])

from sdfiterator import sdf_iter
iter = sdf_iter(inp, True)

cnt = 0
for i in iter:
    bonds_cnt = int((i.split('\n')[3])[3:6])
    if bonds_cnt >= limit:
        outp.write(i)
    else:
        sys.stderr.write('skipping %s\n' % cnt)
    cnt += 1
