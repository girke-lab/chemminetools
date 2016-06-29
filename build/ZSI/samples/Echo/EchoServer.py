#
# How to build an echo server using the extended code generation
#

# Import the ZSI stuff you'd need no matter what
from ZSI.ServiceContainer import ServiceContainer

# This is a new method imported to show it's value
from ZSI.ServiceContainer import GetSOAPContext

# Import the generated Server Object
import EchoServer_interface

# Create a Server implementation

# Here we use a Tie approach
class EchoServer:
    def authorize(self, auth_info, post, action):
        print "NOT Authorizing TIE Echo"
        ctx = GetSOAPContext()
        print dir(ctx)
        print "Container: ", ctx.connection
        print "Parsed SOAP: ", ctx.parsedsoap
        print "Container: ", ctx.container
        print "HTTP Headers:\n", ctx.httpheaders
        print "----"
        print "XML Data:\n", ctx.xmldata
        return 0

    def Echo(self, input):
        return input * 3
    
# This using a derived server instead
class EchoServIn(EchoServer_interface.EchoServer):
    def __init__(self, post='', **kw):
        EchoServer_interface.EchoServer.__init__(self, post, kw=kw)
        
    def authorize(self, auth_info, post, action):
        print "Authorizing INHERIT Echo"
        ctx = GetSOAPContext()
        print dir(ctx)
        print "Container: ", ctx.connection
        print "Parsed SOAP: ", ctx.parsedsoap
        print "Container: ", ctx.container
        print "HTTP Headers:\n", ctx.httpheaders
        print "----"
        print "XML Data:\n", ctx.xmldata
        return 1

    def Echo(self, input):
        return input * 3

# Here we set up the server
serviceContainer = ServiceContainer(('localhost', 9999))

# Create the TIE version of the server
hws = EchoServer()
hwsi = EchoServer_interface.EchoServer(impl=hws,
                                                auth_method_name="authorize")
serviceContainer.setNode(hwsi, url="/EchoServer")

# Create the Inherited version of the server
hws2 = EchoServIn()
serviceContainer.setNode(hws2, url="/EchoServIn")

# Run the service container
serviceContainer.serve_forever()
