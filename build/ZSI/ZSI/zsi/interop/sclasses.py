#! /usr/bin/env python

WSDL_DEFINITION = '''<?xml version="1.0"?>

<definitions name="InteropTest"
    targetNamespace="http://soapinterop.org/" 
    xmlns="http://schemas.xmlsoap.org/wsdl/" 
    xmlns:soap="http://schemas.xmlsoap.org/wsdl/soap/" 
    xmlns:tns="http://soapinterop.org/">

  <import
      location="http://www.whitemesa.com/interop/InteropTest.wsdl"
      namespace="http://soapinterop.org/xsd"/>
  <import
      location="http://www.whitemesa.com/interop/InteropTest.wsdl"
      namespace="http://soapinterop.org/"/>
  <import
      location="http://www.whitemesa.com/interop/InteropTestB.wsdl"
      namespace="http://soapinterop.org/"/>
  <import
      location="http://www.whitemesa.com/interop/echoHeaderBindings.wsdl"
      namespace="http://soapinterop.org/"/>
  <import
      location="http://www.whitemesa.com/interop/InteropTestMap.wsdl"
      namespace="http://soapinterop.org/"/>
<!--  DOCSTYLE; soon.
  <import
      location="http://www.whitemesa.com/interop/interopdoc.wsdl"
      namespace="http://soapinterop.org/"/>
-->

  <service name="interop">
    <port name="TestSoap" binding="tns:InteropTestSoapBinding">
      <soap:address location=">>>URL<<<"/>
    </port>
    <port name="TestSoapB" binding="tns:InteropTestSoapBindingB">
      <soap:address location=">>>URL<<<"/>
    </port>
    <port name="EchoHeaderString" binding="tns:InteropEchoHeaderStringBinding">
      <soap:address location=">>>URL<<<"/>
    </port>
    <port name="EchoHeaderStruct" binding="tns:InteropEchoHeaderStructBinding">
      <soap:address location=">>>URL<<<"/>
    </port>
    <port name="TestSoapMap" binding="tns:InteropTestSoapBindingMap">
      <soap:address location=">>>URL<<<"/>
    </port>
<!--  DOCSTYLE; soon.
    <port name="TestDoc" binding="tns:doc_test_binding">
      <soap:address location=">>>URL<<<"/>
    </port>
-->
  </service>

</definitions>
'''

from ZSI import *
from ZSI import _copyright, _seqtypes
import types

class SOAPStruct:
    def __init__(self, name):
        pass
    def __str__(self):
        return str(self.__dict__)

class TC_SOAPStruct(TC.Struct):
    def __init__(self, pname=None, **kw):
        TC.Struct.__init__(self, SOAPStruct, [
            TC.String('varString', strip=0, inline=1),
            TC.Iint('varInt'),
            TC.FPfloat('varFloat', format='%.18g'),
        ], pname, **kw)

class TC_SOAPStructStruct(TC.Struct):
    def __init__(self, pname=None, **kw):
        TC.Struct.__init__(self, SOAPStruct, [
            TC.String('varString', strip=0),
            TC.Iint('varInt'),
            TC.FPfloat('varFloat', format='%.18g'),
            TC_SOAPStruct('varStruct'),
        ], pname, **kw)

class TC_SOAPArrayStruct(TC.Struct):
    def __init__(self, pname=None, **kw):
        TC.Struct.__init__(self, SOAPStruct, [
            TC.String('varString', strip=0),
            TC.Iint('varInt'),
            TC.FPfloat('varFloat', format='%.18g'),
            TC.Array('xsd:string', TC.String(string=0), 'varArray'),
        ], pname, **kw)

class TC_ArrayOfstring(TC.Array):
    def __init__(self, pname=None, **kw):
        TC.Array.__init__(self, 'xsd:string', TC.String(string=0), pname, **kw)

class TC_ArrayOfint(TC.Array):
    def __init__(self, pname=None, **kw):
        TC.Array.__init__(self, 'xsd:int', TC.Iint(), pname, **kw)

class TC_ArrayOffloat(TC.Array):
    def __init__(self, pname=None, **kw):
        TC.Array.__init__(self, 'xsd:float', TC.FPfloat(format='%.18g'),
                    pname, **kw)

class TC_ArrayOfSOAPStruct(TC.Array):
    def __init__(self, pname=None, **kw):
        TC.Array.__init__(self, 'Za:SOAPStruct', TC_SOAPStruct(), pname, **kw)

#class TC_ArrayOfstring2D(TC.Array):
#    def __init__(self, pname=None, **kw):
#       TC.Array.__init__(self, 'xsd:string', TC.String(string=0), pname, **kw)

class RPCParameters:
    def __init__(self, name):
        pass
    def __str__(self):
        t = str(self.__dict__)
        if hasattr(self, 'inputStruct'):
            t += '\ninputStruct\n'
            t += str(self.inputStruct)
        if hasattr(self, 'inputStructArray'):
            t += '\ninputStructArray\n'
            t += str(self.inputStructArray)
        return t
    def frominput(self, arg):
        self.v = s = SOAPStruct(None)
        self.v.varString = arg.inputString
        self.v.varInt = arg.inputInteger
        self.v.varFloat = arg.inputFloat
        return self

class Operation:
    dispatch = {}
    SOAPAction = '''"http://soapinterop.org/"'''
    ns = "http://soapinterop.org/"
    hdr_ns = "http://soapinterop.org/echoheader/"

    def __init__(self, name, tcin, tcout, **kw):
        self.name = name
        if type(tcin) not in _seqtypes: tcin = tcin,
        self.TCin = TC.Struct(RPCParameters, tuple(tcin), name)
        if type(tcout) not in _seqtypes: tcout = tcout,
        self.TCout = TC.Struct(RPCParameters, tuple(tcout), name + 'Response')
        self.convert = kw.get('convert', None)
        self.headers = kw.get('headers', [])
        self.nsdict = kw.get('nsdict', {})
        Operation.dispatch[name] = self

Operation("echoString",
    TC.String('inputString', strip=0),
    TC.String('inputString', oname='return', strip=0)
)
Operation("echoStringArray",
    TC_ArrayOfstring('inputStringArray'),
    TC_ArrayOfstring('inputStringArray', oname='return')
)
Operation("echoInteger",
    TC.Iint('inputInteger'),
    TC.Iint('inputInteger', oname='return'),
)
Operation("echoIntegerArray",
    TC_ArrayOfint('inputIntegerArray'),
    TC_ArrayOfint('inputIntegerArray', oname='return'),
)
Operation("echoFloat",
    TC.FPfloat('inputFloat', format='%.18g'),
    TC.FPfloat('inputFloat', format='%.18g', oname='return'),
)
Operation("echoFloatArray",
    TC_ArrayOffloat('inputFloatArray'),
    TC_ArrayOffloat('inputFloatArray', oname='return'),
)
Operation("echoStruct",
    TC_SOAPStruct('inputStruct'),
    TC_SOAPStruct('inputStruct', oname='return'),
)
Operation("echoStructArray",
    TC_ArrayOfSOAPStruct('inputStructArray'),
    TC_ArrayOfSOAPStruct('inputStructArray', oname='return'),
    nsdict={'Za': 'http://soapinterop.org/xsd'}
)
Operation("echoVoid",
    [],
    [],
    headers=( ( Operation.hdr_ns, 'echoMeStringRequest' ),
                ( Operation.hdr_ns, 'echoMeStructRequest' ) )
)
Operation("echoBase64",
    TC.Base64String('inputBase64'),
    TC.Base64String('inputBase64', oname='return'),
)
Operation("echoDate",
    TC.gDateTime('inputDate'),
    TC.gDateTime('inputDate', oname='return'),
)
Operation("echoHexBinary",
    TC.HexBinaryString('inputHexBinary'),
    TC.HexBinaryString('inputHexBinary', oname='return'),
)
Operation("echoDecimal",
    TC.Decimal('inputDecimal'),
    TC.Decimal('inputDecimal', oname='return'),
)
Operation("echoBoolean",
    TC.Boolean('inputBoolean'),
    TC.Boolean('inputBoolean', oname='return'),
)
Operation("echoStructAsSimpleTypes",
    TC_SOAPStruct('inputStruct'),
    ( TC.String('outputString', strip=0), TC.Iint('outputInteger'),
        TC.FPfloat('outputFloat', format='%.18g') ),
    convert=lambda s: (s.v.varString, s.v.varInt, s.v.varFloat),
)
Operation("echoSimpleTypesAsStruct",
    ( TC.String('inputString', strip=0), TC.Iint('inputInteger'),
        TC.FPfloat('inputFloat') ),
    TC_SOAPStruct('v', opname='return'),
    convert=lambda arg: RPCParameters(None).frominput(arg),
)
#Operation("echo2DStringArray",
#    TC_ArrayOfstring2D('input2DStringArray'),
#    TC_ArrayOfstring2D('return')
#),
Operation("echoNestedStruct",
    TC_SOAPStructStruct('inputStruct'),
    TC_SOAPStructStruct('inputStruct', oname='return'),
)
Operation("echoNestedArray",
    TC_SOAPArrayStruct('inputStruct'),
    TC_SOAPArrayStruct('inputStruct', oname='return'),
)
