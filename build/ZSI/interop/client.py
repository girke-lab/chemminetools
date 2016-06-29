#! /usr/bin/env python
import getopt, socket, sys
from cpackets import testlist

try:
    (opts, args) = getopt.getopt(sys.argv[1:],
                    'h:lp:qst:w',
                    ( 'host=', 'list', 'port=',
                        'quit', 'statusonly', 'test=', 'wsdl', 'help'))
except getopt.GetoptError, e:
    print >>sys.stderr, sys.argv[0] + ': ' + str(e)
    sys.exit(1)
if args:
    print sys.argv[0] + ': Usage error; try --help.'
    sys.exit(1)

hostname, portnum, tests, quitting, getwsdl, verbose = \
        'localhost', 1122, [0,1], 0, 0, 1
for opt, val in opts:
    if opt in [ '--help' ]:
        print '''Options include:
    --host HOST (-h HOST)       Name of server host
    --port PORT (-p PORT)       Port server is listening on
    --quit (-q)                 Send server a QUIT command
    --testnum 1,2,3 (-t ...)    Run comma-separated tests; use * or all for all
    --list (-l)                 List tests (brief description)
    --statusonly (-s)           Do not output reply packets; just status code
    --wsdl (-w)                 Get the WSDL file
Default is -h%s -p%d -t%s''' % \
    (hostname, portnum, ','.join([str(x) for x in tests]))
        sys.exit(0)
    if opt in [ '-h', '--host' ]:
        hostname = val
    elif opt in [ '-p', '--port' ]:
        portnum = int(val)
    elif opt in [ '-s', '--statusonly' ]:
        verbose = 0
    elif opt in [ '-q', '--quit' ]:
        quitting = 1
    elif opt in [ '-t', '--testnum' ]:
        if val in [ '*', 'all' ]:
            tests = range(len(testlist))
        else:
            tests = [ int(t) for t in val.split(',') ]
    elif opt in [ '-l', '--list' ]:
        for i in range(len(testlist)):
            print i, testlist[i][0]
        sys.exit(0)
    elif opt in [ '-w', '--wsdl' ]:
        getwsdl = 1

if quitting:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.connect((hostname, portnum))
    except socket.error, e:
        if e.args[1] == 'Connection refused': sys.exit(0)
        raise
    f = s.makefile('r+')
    f.write('QUIT / HTTP/1.0\r\n')
    f.flush()
    sys.exit(0)

if getwsdl:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((hostname, portnum))
    f = s.makefile('r+')
    f.write('GET /wsdl HTTP/1.0\r\n\r\n')
    f.flush()
    status = f.readline()
    print status,
    while 1:
        l = f.readline()
        if l == '': break
        print l,
    sys.exit(0)

for T in tests:
    descr, IN, header = testlist[T]
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((hostname, portnum))
    f = s.makefile('r+')

    print '-' * 60, '\n\n\n', T, descr
    f.write('POST / HTTP/1.0\r\n')
    f.write('SOAPAction: "http://soapinterop.org/"\r\n')
    if header == None:
        f.write('Content-type: text/xml; charset="utf-8"\r\n')
        f.write('Content-Length: %d\r\n\r\n' % len(IN))
    else:
        f.write(header)
    f.write(IN)
    f.flush()

    status = f.readline()
    print status,
    while 1:
        l = f.readline()
        if l == '': break
        if verbose: print l,

    f.close()
