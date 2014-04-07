#!/usr/bin/python
# -*- coding: utf-8 -*-

"""create a db to handle random access to large SDF file"""

from anydbm import open
from sdfiterator import sdf_iter
import logging
from logging import info, debug, NOTSET, DEBUG, INFO

logging.basicConfig(level=NOTSET)


class SDFDB(object):

    def __init__(self, dbpath=None):
        if dbpath is None:
            self.db = None
        else:
            info('opening %s' % dbpath)
            self.db = open(dbpath)
            info('%s opened' % dbpath)

    def newdb(self, filepath, dbpath=None):
        """take a filepath, and store all SDFs inside to a database"""

        if dbpath is None:
            dbpath = filepath + '.db'
        info('opening %s for writing' % dbpath)
        db = open(dbpath, 'n')
        cntr = 1
        for sdf in sdf_iter(filepath):
            info(str(cntr))
            db[str(cntr)] = sdf
            cntr += 1
        info('writing %s finished' % dbpath)
        if self.db:
            info('rebinding to %s' % dbpath)
            self.db.close()
        else:
            info('binding to %s' % dbpath)
        self.db = db

    def close(self):
        if self.db:
            info('closing bindings')
            self.db.close()
            self.db = None
        else:
            debug('nothing is opened')

    def __getitem__(self, index):
        if self.db is None:
            critical("access SDFDB instance when it's not bound")
            raise Exception('SDFDB is not bound')
        info('accessing with index %s' % index)
        if not isinstance(index, str):
            index = str(index)
        return self.db[index]


if __name__ == '__main__':
    logging.basicConfig(level=DEBUG)
    import sys
    db = SDFDB()
    db.newdb(sys.argv[1])
    db.close()
