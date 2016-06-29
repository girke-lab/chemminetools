#
# Exampe echo client, to show extended code generation in ZSI
#
# Import the client proxy object
from EchoServer_client import EchoServerSOAP

# Instantiate a client proxy object, then call it
hw = EchoServerSOAP("http://localhost:9999/EchoServer")
print hw.Echo("Test TIE String")

# Instantiate a client proxy object, then call it
hw2 = EchoServerSOAP("http://localhost:9999/EchoServIn")
try:
    print hw2.Echo("Test INHERIT String")
except Exception, e:
    print "Failed to echo (Inherited): ", e

