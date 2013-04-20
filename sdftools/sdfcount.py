#!/usr/bin/python
# -*- coding: utf-8 -*-

"""remove sdf with 0 bonds"""

import sys
inp = sys.argv[1]

from sdfiterator import sdf_iter
iter = sdf_iter(inp)

cntr = 0
for i in iter:
    cntr += 1

print cntr
