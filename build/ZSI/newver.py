#! /usr/bin/env python

# This script does not use 2.x features so that you can build an RPM
# on a system that doesn't have 2.x installed (e.g., basic RedHat).

import ConfigParser
import os
import sys
import time

cf = ConfigParser.ConfigParser()
cf.read('setup.cfg')

major = cf.getint('version', 'major')
minor = cf.getint('version', 'minor')
patchlevel = cf.getint('version', 'patchlevel')

if '--incr' in sys.argv[1:]:
    patchlevel = patchlevel + 1
    cf.set('version', 'patchlevel', patchlevel)
    f = open('setup.cfg', 'w')
    cf.write(f)
    f.close()

def write_changes(filename, text):
    """Write text to the file only if the text changed."""
    oldtext = None
    if os.path.exists(filename):
        f = open(filename)
        oldtext = f.read()
        f.close()
    if oldtext != text:
        f = open(filename, 'w')
        f.write(text)
        f.close()

write_changes('ZSI/version.py',
              '# Auto-generated file; do not edit\n'
              'Version = (%d, %d, %d)\n' % (major, minor, patchlevel))

write_changes('doc/version.tex',
              '% Auto-generated file; do not edit\n'
              + time.strftime('\\date{%B %d, %Y}\n', time.localtime())
              + ('\\release{%d.%d.%d}\n' % (major, minor, patchlevel))
              + ('\\setshortversion{%d.%d}\n' % (major, minor)))
