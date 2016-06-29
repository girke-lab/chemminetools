#! /usr/bin/env python

from ZSI import *

"""Polymorphic containers using TC.Choice

Based on code and text from Dan Gunter <dkgunter@lbl.gov>.

The TC.Choice typecode can be used to create a polymorphic type for
a specified set of object types.  Here's how:
    1.  Define all your data classes (D1, D2, ...). If they derive from
        a common base class (Base), some feel you get "cleaner" code.
    2.  Create a typecode class (TC_D1, TC_D2, ...) for each data class.
    3.  Create a "base" typecode that uses TC.Choice to do the actual
        parsing and serializing.  Two versions are shown, below.
Then you can instantiate, e.g., an Array that can handle multiple
datatypes with: 
    TC.Array("base-class-name", TC_Base(), "MyArray") 
"""

class Base: pass
class D1(Base): pass
class D2(Base): pass

class TC_D1(TC.TypeCode): pass
class TC_D2(TC.TypeCode): pass
D1.typecode = TC_D1()
D2.typecode = TC_D2()

# A simple version of TC_Base that is "hardwired" with the types of
# objects it can handle.  We defer setting the choice attribute because
# with nested containers you could get too much recursion.
class TC_Base(TC.TypeCode):
    def parse(self, elt, ps):
        return self.choice.parse(elt, ps)[1]

    def serialize(self, sw, pyobj, **kw):
        if not isinstance(pyobj, Base):
            raise TypeError(str(pyobj.__class__) + " not in type hierarchy")
        if isinstance(pyobj, D1):
            self.choice.serialize(sw, ('D1', pyobj), **kw)
        elif isinstance(pyobj, D2):
            self.choice.serialize(sw, ('D2', pyobj), **kw)
        else:
            raise TypeError(str(pyobj.__class__) + " unknown type")
        return

    def __getattr__(self, attr):
        if attr == 'choice':
            choice = TC.Choice((D1.typecode, D2.typecode), 'Item')
            self.__dict__['choice'] = choice
            return choice
        raise AttributeError(attr)

## Another version that takes a dictionary that maps element names to
## the python class.

class TC_Polymorphic(TC.TypeCode):
    def __init__(self, name2class, pname=None, **kw):
        TC.TypeCode.__init__(self, pname, **kw)
        self.name2class = name2class

    def parse(self, elt, ps):
        return self.choice.parse(elt, ps)[1]

    def serialize(self, sw, pyobj, **kw):
        self.choice.serialize(sw,
            (self.class2name[pyobj.__class__], pyobj), **kw)

    def __getattr__(self, attr):
        if attr == 'choice':
            choice = TC.Choice(
                [getattr(v, 'typecode') for k,v in self.name2class.items()],
                'Item')
            self.__dict__['choice'] = choice
            return choice
        if attr == 'class2name':
            class2name = {}
            for k,v in self.name2class.items(): class2name[v] = k
            self.__dict__['class2name'] = class2name
            return class2name
        raise AttributeError(attr)

class P1: pass
class P2: pass
P1.typecode = TC.String('s')
P2.typecode = TC.Integer('i')
myTC = TC.Array("Base", TC_Polymorphic({'i': P2, 's': P1}))

test = '''<SOAP-ENV:Envelope
  xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
  xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
  xmlns:xsd="http://www.w3.org/2001/XMLSchema"
  SOAP-ENV:encodingStyle="http://schemas.xmlsoap.org/soap/encoding/">
    <SOAP-ENV:Body xmlns='test-uri'>
        <array SOAP-ENC:arrayType="Base">
            <i>34</i>
            <s>hello</s>
            <s>34</s>
            <i>12</i>
            <s>goodbye</s>
        </array>
    </SOAP-ENV:Body>
</SOAP-ENV:Envelope>'''

ps = ParsedSoap(test)
a = myTC.parse(ps.body_root, ps)
print a

if 0:
    # XXX.  Does not work. :(
    b = [ P1(), P1(), P2(), P1() ]
    b[0].s = 'string'
    b[1].s = '34'
    b[2].i = 34
    b[3].s = 'adios'

    import sys
    sw = SoapWriter(sys.stdout)
    myTC.serialize(sw, b)
