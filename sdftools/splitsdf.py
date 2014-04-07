#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
cid = ''
for i in sys.stdin.xreadlines():
    if not cid:
        if not i.strip():
            continue
        assert i.startswith('NSC')
        cid = i.strip()[len('NSC '):]

        # f = file('%s.sdf' % cid, 'w')

        print '%s' % cid

    # f.write(i)

    if i.startswith('$$$$'):
        cid = ''

        # f.close()
