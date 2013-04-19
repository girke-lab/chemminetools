#!/usr/bin/python
# -*- coding: utf-8 -*-

"""remove sdf with 0 bonds"""

import sys
inp = sys.argv[1]

from sdfiterator import sdf_iter
iter = sdf_iter(inp)

asum = 0
bsum = 0
cntr = 0
for i in iter:
    asum += int((i.split('\n', 5)[3])[0:3])
    bsum += int((i.split('\n', 5)[3])[3:6])
    cntr += 1
print asum / cntr
print bsum / cntr
