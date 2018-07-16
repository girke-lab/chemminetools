#!/usr/bin/python
# -*- coding: utf-8 -*-

"""remove sdf with 0 bonds"""

import sys
inp = sys.argv[1]
outp = file(sys.argv[2], 'w')

from sdfiterator import sdf_iter
iter = sdf_iter(inp)

cnt = 0
for i in iter:
    bonds_cnt = int((i.split('\n')[3])[3:6])
    if bonds_cnt != 0:
        outp.write(i)
    else:
        print 'skipping ' + cnt
    cnt += 1
