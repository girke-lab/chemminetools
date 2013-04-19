#!/usr/bin/python
# -*- coding: utf-8 -*-

"""clean lines over 80 characters"""

import sys
input = file(sys.argv[1])
output = file(sys.argv[2], 'w')

for i in input.xreadlines():
    if len(i) > 80:
        i = i[:80] + '\n'
    output.write(i)

input.close()
output.close()
