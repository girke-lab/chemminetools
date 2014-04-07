#!/usr/bin/python
# -*- coding: utf-8 -*-

# bring in the ZSI-generated interface

from PUG_services import *

# other modules/functions

from time import sleep


def SimilaritySearch(smiles, similaritythreshold=95, maxReturned=200):

    # get a PUG SOAP port instance

    loc = PUGLocator()
    port = loc.getPUGSoap()

    # input o-bromopyridine in SMILES format

    req = InputStructureSoapIn()
    req.set_element_structure(smiles)
    req.set_element_format('eFormat_SMILES')
    strKey = port.InputStructure(req).get_element_StrKey()

    # print 'StrKey =', strKey

    # Initialize similarity search 

    req = SimilaritySearch2DSoapIn()
    req.set_element_StrKey(strKey)
    req.set_element_simOptions(req.new_simOptions())
    req.get_element_simOptions().set_element_threshold(similaritythreshold)
    req.set_element_limits(req.new_limits())
    req.get_element_limits().set_element_maxRecords(maxReturned)
    listKey = port.SimilaritySearch2D(req).get_element_ListKey()


    # print 'ListKey =', listKey

    # call GetOperationStatus until the operation is finished

    req = GetOperationStatusSoapIn()
    req.set_element_AnyKey(listKey)
    status = port.GetOperationStatus(req).get_element_status()
    while status == 'eStatus_Queued' or status == 'eStatus_Running':

        # print 'Waiting for operation to finish...'

        sleep(10)
        status = port.GetOperationStatus(req).get_element_status()

    # check status

    if status == 'eStatus_Success' or status == 'eStatus_TimeLimit' \
        or status == 'eStatus_HitLimit':

        # if (status == 'eStatus_TimeLimit'):
            # print 'Warning: time limit reached before entire db searched'
        # elif (status == 'eStatus_HitLimit'):
            # print 'Warning: hit limit reached before entire db searched'

        # get Entrez key

        req = GetEntrezKeySoapIn()
        req.set_element_ListKey(listKey)
        entrezKey = port.GetEntrezKey(req).get_element_EntrezKey()

        # get URL

        req = GetEntrezUrlSoapIn()
        req.set_element_EntrezKey(entrezKey)
        url = port.GetEntrezUrl(req).get_element_url()

    # get ID list

        req = GetIDListSoapIn()
        req.set_element_ListKey(listKey)
        id_list = \
            port.GetIDList(req).get_element_IDList().get_element_int()
        return id_list
    else:

            # status indicates error

        # see if there's some explanatory message

        req = GetStatusMessageSoapIn()
        req.set_element_AnyKey(listKey)

        # print 'Error:', port.GetStatusMessage(req).get_element_message()

        return 'error'


if __name__ == '__main__':
    SimilaritySearch()

# $Id: SimilaritySearch.py 127842 2008-05-15 21:52:43Z thiessen $
