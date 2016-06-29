#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys
from ZSI.ServiceContainer import AsServer
from EchoServer_services_server import EchoServer

"""
EchoServer example service

WSDL:  ../../samples/Echo/Echo.wsdl

"""


class Service(EchoServer):
    def soap_Echo(self, ps):
        response = EchoServer.soap_Echo(self, ps)
        response.EchoResult = self.request.EchoIn
        return response


if __name__ == "__main__" :
    port = int(sys.argv[1])
    AsServer(port, (Service('test'),))
