#! /usr/bin/env python
from ZSI import *
from ZSI import _copyright, resolvers, _child_elements, _textprotect
import sys, time, cStringIO as StringIO
from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer

from sclasses import Operation, WSDL_DEFINITION, TC_SOAPStruct

class InteropRequestHandler(BaseHTTPRequestHandler):
    server_version = 'ZSI/1.2 OS390/VM5.4 ' + BaseHTTPRequestHandler.server_version

    def send_xml(self, text, code=200):
        '''Send some XML.'''
        self.send_response(code)
        self.send_header('Content-type', 'text/xml; charset="utf-8"')
        self.send_header('Content-Length', str(len(text)))
        self.end_headers()
        self.wfile.write(text)
        self.trace(text, 'SENT')
        self.wfile.flush()

    def send_fault(self, f):
        '''Send a fault.'''
        self.send_xml(f.AsSOAP(), 500)

    def trace(self, text, what):
        '''Log a debug/trace message.'''
        F = self.server.tracefile
        if not F: return
        print >>F, '=' * 60, '\n%s %s %s %s:' % \
            (what, self.client_address, self.path, time.ctime(time.time()))
        print >>F, text
        print >>F, '=' * 60, '\n'
        F.flush()

    def do_QUIT(self):
        '''Quit.'''
        self.server.quitting = 1
        self.log_message('Got QUIT command')
        sys.stderr.flush()
        raise SystemExit

    def do_GET(self):
        '''The GET command.  Always returns the WSDL.'''
        self.send_xml(WSDL_DEFINITION.replace('>>>URL<<<', self.server.url))

    def do_POST(self):
        '''The POST command.'''

        try:
            # SOAPAction header.
            action = self.headers.get('soapaction', None)
            if not action:
                self.send_fault(Fault(Fault.Client,
                                    'SOAPAction HTTP header missing.'))
                return
            if action != Operation.SOAPAction:
                self.send_fault(Fault(Fault.Client,
                    'SOAPAction is "%s" not "%s"' % \
                    (action, Operation.SOAPAction)))
                return

            # Parse the message.
            ct = self.headers['content-type']
            if ct.startswith('multipart/'):
                cid = resolvers.MIMEResolver(ct, self.rfile)
                xml = cid.GetSOAPPart()
                ps = ParsedSoap(xml, resolver=cid.Resolve)
            else:
                cl = int(self.headers['content-length'])
                IN = self.rfile.read(cl)
                self.trace(IN, 'RECEIVED')
                ps = ParsedSoap(IN)
        except ParseException, e:
            self.send_fault(FaultFromZSIException(e))
            return
        except Exception, e:
            # Faulted while processing; assume it's in the header.
            self.send_fault(FaultFromException(e, 1, sys.exc_info()[2]))
            return

        try:
            # Actors?
            a = ps.WhatActorsArePresent()
            if len(a):
                self.send_fault(FaultFromActor(a[0]))
                return

            # Is the operation defined?
            root = ps.body_root
            if root.namespaceURI != Operation.ns:
                self.send_fault(Fault(Fault.Client,
                    'Incorrect namespace "%s"' % root.namespaceURI))
                return
            n = root.localName
            op = Operation.dispatch.get(n, None)
            if not op:
                self.send_fault(Fault(Fault.Client,
                    'Undefined operation "%s"' % n))
                return

            # Scan headers.  First, see if we understand all headers with
            # mustUnderstand set. Then, get the ones intended for us (ignoring
            # others since step 1 insured they're not mustUnderstand).
            for mu in ps.WhatMustIUnderstand():
                if mu not in op.headers:
                    uri, localname = mu
                    self.send_fault(FaultFromNotUnderstood(uri, localname))
                    return
            headers = [ e for e in ps.GetMyHeaderElements()
                        if (e.namespaceURI, e.localName) in op.headers ]
            nsdict={ 'Z': Operation.ns }
            if headers:
                nsdict['E'] = Operation.hdr_ns
                self.process_headers(headers, ps)
            else:
                self.headertext = None

            try:
                results = op.TCin.parse(ps.body_root, ps)
            except ParseException, e:
                self.send_fault(FaultFromZSIException(e))
            self.trace(str(results), 'PARSED')
            if op.convert:
                results = op.convert(results)
            if op.nsdict: nsdict.update(op.nsdict)
            reply = StringIO.StringIO()
            sw = SoapWriter(reply, nsdict=nsdict, header=self.headertext)
            sw.serialize(results, op.TCout,
                    name = 'Z:' + n + 'Response', inline=1)
            sw.close()
            self.send_xml(reply.getvalue())
        except Exception, e:
            # Fault while processing; now it's in the body.
            self.send_fault(FaultFromException(e, 0, sys.exc_info()[2]))
            return

    def process_headers(self, headers, ps):
        '''Process headers, set self.headertext to be what to output.
        '''
        self.headertext = ''
        for h in headers:
            if h.localName == 'echoMeStringRequest':
                s = TC.String().parse(h, ps)
                self.headertext += \
    '<E:echoMeStringResponse>%s</E:echoMeStringResponse>\n' % _textprotect(s)
            elif h.localName == 'echoMeStructRequest':
                tc = TC_SOAPStruct('echoMeStructRequest', inline=1)
                data = tc.parse(h, ps)
                s = StringIO.StringIO()
                sw = SoapWriter(s, envelope=0)
                tc.serialize(sw, data, name='E:echoMeStructResponse')
                sw.close()
                self.headertext += s.getvalue()
            else:
                raise TypeError('Unhandled header ' + h.nodeName)
        pass


class InteropHTTPServer(HTTPServer):
    def __init__(self, me, url, **kw):
        HTTPServer.__init__(self, me, InteropRequestHandler)
        self.quitting = 0
        self.tracefile = kw.get('tracefile', None)
        self.url = url
    def handle_error(self, req, client_address):
        if self.quitting: sys.exit(0)
        HTTPServer.handle_error(self, req, client_address)


import getopt
try:
    (opts, args) = getopt.getopt(sys.argv[1:], 'l:p:t:u:',
                        ('log=', 'port=', 'tracefile=', 'url=') )
except getopt.GetoptError, e:
    print >>sys.stderr, sys.argv[0] + ': ' + str(e)
    sys.exit(1)
if args:
    print >>sys.stderr, sys.argv[0] + ': Usage error.'
    sys.exit(1)

portnum = 1122
tracefile = None
url = None
for opt, val in opts:
    if opt in [ '-l', '--logfile' ]:
        sys.stderr = open(val, 'a')
    elif opt in [ '-p', '--port' ]:
        portnum = int(val)
    elif opt in [ '-t', '--tracefile' ]:
        if val == '-':
            tracefile = sys.stdout
        else:
            tracefile = open(val, 'a')
    elif opt in [ '-u', '--url' ]:
        url = val
ME = ( '', portnum )

if not url:
    import socket
    url = 'http://' + socket.getfqdn()
    if portnum != 80: url += ':%d' % portnum
    url += '/interop.wsdl'

try:
    InteropHTTPServer(ME, url, tracefile=tracefile).serve_forever()
except SystemExit:
    pass
sys.exit(0)
