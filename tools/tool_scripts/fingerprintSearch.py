#!/usr/bin/python
# -*- coding: utf-8 -*-

# make it possible to run as standalone program

import sys
sys.path.append('/srv/chemminetools')
sys.path.append('/srv/chemminetools/sdftools')  # allow tools.py to be imported
from django.core.management import setup_environ
import chemminetools.settings
import argparse
setup_environ(chemminetools.settings)
from django.contrib.auth.models import User
from django.conf import settings

def main():
    sdf = sys.stdin.read()
    user = User.objects.get(id=args['user'])

    output = addMyCompounds(sdf, user)

    f = open(args['outfile'], 'w')
    f.write(output)
    f.close()


if __name__ == '__main__':
    main()
