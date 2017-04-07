#!/usr/bin/python
# -*- coding: utf-8 -*-

# bring in the ZSI-generated interface

from PUG_services import *
import re

# other modules/functions

from time import sleep
from urllib2 import urlopen
import zlib


def DownloadCIDs(cids):

    # get a PUG SOAP port instance

    loc = PUGLocator()
    port = loc.getPUGSoap()

    # start with a list of CIDs

    req = InputListSoapIn()
    req.set_element_ids(req.new_ids())
    req.get_element_ids().set_element_int(cids)
    req.set_element_idType('eID_CID')
    listKey = port.InputList(req).get_element_ListKey()
    # print 'ListKey =', listKey

    # request download in SDF format, gzip-compressed

    req = DownloadSoapIn()
    req.set_element_ListKey(listKey)
    req.set_element_eFormat('eFormat_SDF')
    req.set_element_eCompress('eCompress_GZip')
    downloadKey = port.Download(req).get_element_DownloadKey()
    # print 'DownloadKey =', downloadKey

    # call GetOperationStatus until the operation is finished

    req = GetOperationStatusSoapIn()
    req.set_element_AnyKey(downloadKey)
    status = port.GetOperationStatus(req).get_element_status()
    while status == 'eStatus_Queued' or status == 'eStatus_Running':
        # print 'Waiting for operation to finish...'
        sleep(10)
        status = port.GetOperationStatus(req).get_element_status()

    # check status

    if status == 'eStatus_Success':

        # get the url of the prepared file

        req = GetDownloadUrlSoapIn()
        req.set_element_DownloadKey(downloadKey)
        url = port.GetDownloadUrl(req).get_element_url()
        url = re.sub(r'^ftp://', 'http://', url)
        # print 'Success! URL =', url

        # download to a local file

        sdfCompressed = urlopen(url, timeout=10)
        sdf = zlib.decompress(sdfCompressed.read(), 31)
        sdfCompressed.close()
        return sdf
    else:

            # status indicates error

        # see if there's some explanatory message

        req = GetStatusMessageSoapIn()
        req.set_element_AnyKey(downloadKey)
        # print 'Error:', port.GetStatusMessage(req).get_element_message()


if __name__ == '__main__':
    DownloadCIDs()

# $Id: DownloadCIDs.py 127842 2008-05-15 21:52:43Z thiessen $
