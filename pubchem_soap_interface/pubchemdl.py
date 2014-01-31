#!/usr/bin/python
# -*- coding: utf-8 -*-

# bring in the ZSI-generated interface

from django.utils.daemonize import become_daemon
import os
import sys
sys.path.append('/srv/chemminetools')
sys.path.append('/srv/chemminetools/pubchem_soap_interface')
from django.core.management import setup_environ
import chemminetools.settings
setup_environ(chemminetools.settings)
from pubchem_soap_interface.PUG_services import *

# other modules/functions

from time import sleep, time
from urllib import urlopen


def asynced(_f):

    def func(*args, **kargs):
        try:
            outfile = kargs.pop('asynced_outfile')
        except:
            return _f(*args, **kargs)

        # the following is only executed when in async mode

        try:
            if os.path.exists(outfile):
                t = os.stat(outfile).st_mtime
                if time() - t < 3600:
                    f = file(outfile)
                    c = f.read()
                    f.close()
                    if c and not c.startswith('Error'):
                        return
        except:
            pass

        try:
            if os.path.exists(outfile):
                os.unlink(outfile)
            f = file(outfile, 'w')
            f.close()
        except Exception, e:
            sys.stderr.write('Cannot write to output %s\n' % outfile)
            sys.stderr.write(e.message)

        # try:
        #    if os.fork() > 0:
        #        return
        # except Exception, e:
        #        sys.stderr.write('Cannot fork!')
        #        sys.stderr.write(e.message)
        # become_daemon()

        try:
            out = _f(*args, **kargs)
            err = None
        except Exception, e:
            raise
            err = e.message
        try:
            f = file(outfile, 'w')
            if err:
                f.write('Error:' + err)
            else:
                f.write(out)
            f.close()
        except Exception, e:
            sys.stderr.write('Cannot write to output %s\n' % outfile)
            sys.stderr.write(e.message)

    return func


@asynced
def download(cids, format='sdf'):
    """download SDF/SMILES from PubChem website given the CIDS. If async is
set, it gives a path to a file to write the output URL"""

    if format == 'sdf':
        _format = 'eFormat_SDF'
        _compress = 'eCompress_GZip'
    elif format == 'smiles':
        _format = 'eFormat_SMILES'
        _compress = 'eCompress_None'
    elif format == 'webtools':
        _format = 'eFormat_SMILES'
        _compress = 'eCompress_None'
    else:
        raise Exception('Unknown format: %s' % format)

    # get a PUG SOAP port instance

    loc = PUGLocator()
    port = loc.getPUGSoap()

    # start with a list of CIDs

    req = InputListSoapIn()
    req.set_element_ids(req.new_ids())
    req.get_element_ids().set_element_int(cids)
    req.set_element_idType('eID_CID')
    listKey = port.InputList(req).get_element_ListKey()

    # request download

    req = DownloadSoapIn()
    req.set_element_ListKey(listKey)
    req.set_element_eFormat(_format)
    req.set_element_eCompress(_compress)
    downloadKey = port.Download(req).get_element_DownloadKey()

    # call GetOperationStatus until the operation is finished

    req = GetOperationStatusSoapIn()
    req.set_element_AnyKey(downloadKey)
    status = port.GetOperationStatus(req).get_element_status()
    while status == 'eStatus_Queued' or status == 'eStatus_Running':
        sleep(5)
        status = port.GetOperationStatus(req).get_element_status()

    # check status

    if status == 'eStatus_Success':

        # get the url of the prepared file

        req = GetDownloadUrlSoapIn()
        req.set_element_DownloadKey(downloadKey)
        url = port.GetDownloadUrl(req).get_element_url()
        if format == 'webtools':
            try:
                f = urlopen(url)
                content = f.read()
                f.close()
                return content
            except:
                raise Exception('Error: download %s failed' % url)
        return url
    else:

            # status indicates error

        # see if there's some explanatory message

        req = GetStatusMessageSoapIn()
        req.set_element_AnyKey(downloadKey)
        err = 'Error:' \
            + port.GetStatusMessage(req).get_element_message()
        raise Exception(err)

if __name__ == '__main__':
    f = file(sys.argv[1])
    cids = [int(i) for i in f.read().splitlines()]
    f.close()
    if sys.argv[1].startswith('tmp'):
        os.unlink(sys.argv[1])
    download(cids, asynced_outfile=sys.argv[2], format=sys.argv[3])

# vim: noai:ts=4:sw=4:expandtab
