#! /usr/bin/env python

pocketsoap_struct_test = '''<SOAP-ENV:Envelope
 xmlns="http://www.example.com/schemas/TEST"
 xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
 xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"
 xmlns:ZSI="http://www.zolera.com/schemas/ZSI/"
xmlns:ps42='http://soapinterop.org/xsd'
xmlns:xsd='http://www.w3.org/1999/XMLSchema'
xmlns:xsi='http://www.w3.org/1999/XMLSchema-instance'>
<SOAP-ENV:Header>
<ww:echoMeStringRequest xmlns:ww="http://soapinterop.org/echoheader/">
Please release &lt;me</ww:echoMeStringRequest>
</SOAP-ENV:Header>
<SOAP-ENV:Body>
    <m:echoStruct xmlns:m='http://soapinterop.org/'>
      <inputStruct xsi:type='ps42:SOAPStruct'>
        <varInt xsi:type='xsd:int'>1073741824</varInt>

        <varFloat xsi:type='xsd:float'>-42.24</varFloat>

        <varString xsi:type='xsd:string'>pocketSOAP
        rocks!&lt;g&gt;</varString>
      </inputStruct>
    </m:echoStruct>
</SOAP-ENV:Body>
</SOAP-ENV:Envelope>'''

phalanx_b64_test = '''<SOAP-ENV:Envelope xmlns:SOAP-ENV='http://schemas.xmlsoap.org/soap/envelope/' xmlns:m='http://soapinterop.org/' xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xmlns:SOAP-ENC='http://schemas.xmlsoap.org/soap/encoding/' xmlns:xsd='http://www.w3.org/2001/XMLSchema' SOAP-ENV:encodingStyle='http://schemas.xmlsoap.org/soap/encoding/'>
  <SOAP-ENV:Body>
    <m:echoBase64>
      <inputBase64 xsi:type='xsd:base64Binary'>Ty4rY6==</inputBase64>
    </m:echoBase64>
  </SOAP-ENV:Body>
</SOAP-ENV:Envelope>'''

hexbin_test = '''<SOAP-ENV:Envelope xmlns:SOAP-ENV='http://schemas.xmlsoap.org/soap/envelope/' xmlns:m='http://soapinterop.org/' xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xmlns:SOAP-ENC='http://schemas.xmlsoap.org/soap/encoding/' xmlns:xsd='http://www.w3.org/2001/XMLSchema' SOAP-ENV:encodingStyle='http://schemas.xmlsoap.org/soap/encoding/'>
  <SOAP-ENV:Body>
    <m:echoHexBinary>
      <inputHexBinary xsi:type='xsd:hexBinary'>656174206d792073686f72747321</inputHexBinary>
    </m:echoHexBinary>
  </SOAP-ENV:Body>
</SOAP-ENV:Envelope>'''

phalanx_badhref_test = '''<SOAP-ENV:Envelope xmlns:SOAP-ENV='http://schemas.xmlsoap.org/soap/envelope/' xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xmlns:ns1='http://soapinterop.org/xsd' xmlns:xsd='http://www.w3.org/2001/XMLSchema' xmlns:ns='http://soapinterop.org/' xmlns:SOAP-ENC='http://schemas.xmlsoap.org/soap/encoding/' SOAP-ENV:encodingStyle='http://schemas.xmlsoap.org/soap/encoding/'>
  <SOAP-ENV:Body>
    <ns:echoStructArray>
      <inputStructArray xsi:type='SOAP-ENC:Array' SOAP-ENC:arrayType='ns1:SOAPStruct[3]'>
        <item xsi:type='ns1:SOAPStruct' href='#2'>invalid value</item>
        <item xsi:type='ns1:SOAPStruct'>
          <varFloat xsi:type='xsd:float'>21.02</varFloat>
          <varString xsi:type='xsd:string'>c</varString>
          <varInt xsi:type='xsd:int'>3</varInt>
        </item>
        <item xsi:type='ns1:SOAPStruct' href='#1'/>
      </inputStructArray>
    </ns:echoStructArray>
    <mrA xsi:type='ns1:SOAPStruct' id='1'>
      <varInt xsi:type='xsd:int'>-33</varInt>
      <varFloat xsi:type='xsd:float'>33.33</varFloat>
      <varString xsi:type='xsd:string'>test 1</varString>
    </mrA>
    <mrB xsi:type='ns1:SOAPStruct' id='2'>
      <varFloat xsi:type='xsd:float'>11.11</varFloat>
      <varString xsi:type='xsd:string'>test 2</varString>
      <varInt xsi:type='xsd:int'>-11</varInt>
    </mrB>
  </SOAP-ENV:Body>
</SOAP-ENV:Envelope>'''

someones_b64_test = '''<S:Envelope
      S:encodingStyle='http://schemas.xmlsoap.org/soap/encoding/'
            xmlns:S='http://schemas.xmlsoap.org/soap/envelope/'
                  xmlns:E='http://schemas.xmlsoap.org/soap/encoding/'
                        xmlns:a='http://soapinterop.org/'
                              xmlns:b='http://www.w3.org/2001/XMLSchema-instance'>
                              <S:Body>
                              <a:echoBase64><inputBase64
                              b:type='E:base64'>AAECAwQFBgcICQoLDA0ODxAREhMUFRYXGBkaGxwdHh8gISIjJCUmJygpKissLS4vMDEyMzQ1Njc4OTo7PD0+P0BBQkNERUZHSElKS0xNTk9QUVJTVFVWV1hZWltcXV5fYGFiY2RlZmdoaWprbG1ub3BxcnN0dXZ3eHl6e3x9fn+AgYKDhIWGh4iJiouMjY6PkJGSk5SVlpeYmZqbnJ2en6ChoqOkpaanqKmqq6ytrq+wsbKztLW2t7i5uru8vb6/wMHCw8TFxsfIycrLzM3Oz9DR0tPU1dbX2Nna29zd3t/g4eLj5OXm5+jp6uvs7e7v8PHy8/T19vf4+fr7/P3+AAECAwQFBgcICQoLDA0ODxAREhMUFRYXGBkaGxwdHh8gISIjJCUmJygpKissLS4vMDEyMzQ1Njc4OTo7PD0+P0BBQkNERUZHSElKS0xNTk9QUVJTVFVWV1hZWltcXV5fYGFiY2RlZmdoaWprbG1ub3BxcnN0dXZ3eHl6e3x9fn+AgYKDhIWGh4iJiouMjY6PkJGSk5SVlpeYmZqbnJ2en6ChoqOkpaanqKmqq6ytrq+wsbKztLW2t7i5uru8vb6/wMHCw8TFxsfIycrLzM3Oz9DR0tPU1dbX2Nna29zd3t/g4eLj5OXm5+jp6uvs7e7v8PHy8/T1</inputBase64>
                              </a:echoBase64>
                              </S:Body></S:Envelope>'''

phalanx_void_test = '''<SOAP-ENV:Envelope
 xmlns="http://www.example.com/schemas/TEST"
 xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
 xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"
 xmlns:ZSI="http://www.zolera.com/schemas/ZSI/"
xmlns:ps42='http://soapinterop.org/xsd'
xmlns:xsd='http://www.w3.org/1999/XMLSchema'
xmlns:xsi='http://www.w3.org/1999/XMLSchema-instance'>
<SOAP-ENV:Header>
<ww:echoMeStringRequest SOAP-ENV:mustUnderstand="1" xmlns:ww="http://soapinterop.org/echoheader/">
Please release me</ww:echoMeStringRequest>
<ww:echoMeStructRequest xmlns:ww="http://soapinterop.org/echoheader/">
        <varInt>111</varInt>
        <varFloat>-42.24</varFloat>
        <varString>Header text string.</varString>
      </ww:echoMeStructRequest>
</SOAP-ENV:Header>
<SOAP-ENV:Body>
    <m:echoVoid xmlns:m='http://soapinterop.org/'/>
</SOAP-ENV:Body>
</SOAP-ENV:Envelope>'''

multipart_test = '''
Ignore this
--sep
Content-Type: text/xml

<SOAP-ENV:Envelope
 xmlns="http://www.example.com/schemas/TEST"
 xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
 xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"
 xmlns:ZSI="http://www.zolera.com/schemas/ZSI/"
xmlns:ps42='http://soapinterop.org/xsd'
xmlns:xsd='http://www.w3.org/1999/XMLSchema'
xmlns:xsi='http://www.w3.org/1999/XMLSchema-instance'>
<SOAP-ENV:Body>
    <m:echoStruct xmlns:m='http://soapinterop.org/'>
      <inputStruct xsi:type='ps42:SOAPStruct'>
        <varInt xsi:type='xsd:int'>1073741824</varInt>

        <varFloat xsi:type='xsd:float'>-42.24</varFloat>

        <varString xsi:type='xsd:string' href="cid:123@456"/>
      </inputStruct>
    </m:echoStruct>
</SOAP-ENV:Body>
</SOAP-ENV:Envelope>
--sep
Content-ID: otherwise

yehll

--sep
Content-ID: 123@456

this is a very long string
it is separated over several lines.
hwlleasdfasdf 
asdfad
--sep--
'''
phalanx_badstructtype_test = '''<SOAP-ENV:Envelope xmlns:SOAP-ENV='http://schemas.xmlsoap.org/soap/envelope/' xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xmlns:ns1='http://soapinterop.org/xsd' xmlns:xsd='http://www.w3.org/2001/XMLSchema' xmlns:m='http://soapinterop.org/' xmlns:SOAP-ENC='http://schemas.xmlsoap.org/soap/encoding/' SOAP-ENV:encodingStyle='http://schemas.xmlsoap.org/soap/encoding/'>
  <SOAP-ENV:Body>
    <m:echoStruct>
      <inputStruct xsi:type='ns1:roastbeef'>
        <varString xsi:type='xsd:string'>easy test</varString>
        <varInt xsi:type='xsd:int'>11</varInt>
        <varFloat xsi:type='xsd:float'>22.33</varFloat>
      </inputStruct>
    </m:echoStruct>
  </SOAP-ENV:Body>
</SOAP-ENV:Envelope>'''

phalanx_int_href_test = '''<SOAP-ENV:Envelope xmlns:SOAP-ENV='http://schemas.xmlsoap.org/soap/envelope/' xmlns:m='http://soapinterop.org/' xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xmlns:SOAP-ENC='http://schemas.xmlsoap.org/soap/encoding/' xmlns:xsd='http://www.w3.org/2001/XMLSchema' SOAP-ENV:encodingStyle='http://schemas.xmlsoap.org/soap/encoding/'>
  <SOAP-ENV:Body>
    <m:echoInteger>
      <inputInteger href='#a1'/>
    </m:echoInteger>
    <multiref xsi:type='xsd:int' id='a1'>13</multiref>
  </SOAP-ENV:Body>
</SOAP-ENV:Envelope>
'''

wm_simple_as_struct_test = '''<?xml version="1.0" encoding="UTF-8"?>
<SOAP-ENV:Envelope
SOAP-ENV:encodingStyle="http://schemas.xmlsoap.org/soap/encoding/"
xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/"
xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
xmlns:xsd="http://www.w3.org/2001/XMLSchema"
xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<SOAP-ENV:Body>
<m:echoSimpleTypesAsStruct xmlns:m="http://soapinterop.org/">
<inputString>White Mesa Test</inputString>
<inputInteger>42</inputInteger>
<inputFloat>0.0999</inputFloat>
</m:echoSimpleTypesAsStruct>
</SOAP-ENV:Body>
</SOAP-ENV:Envelope>
'''

apache_float_test = '''<?xml version="1.0" encoding="UTF-8"?>
<SOAP-ENV:Envelope SOAP-ENV:encodingStyle="http://schemas.xmlsoap.org/soap/encoding/" xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
 <SOAP-ENV:Body>
   <ns1:echoFloat xmlns:ns1="http://soapinterop.org/">
      <inputFloat xsi:type="xsd:float">3.7</inputFloat>
        </ns1:echoFloat>
         </SOAP-ENV:Body>
</SOAP-ENV:Envelope>'''

testlist = (
    ( 'struct test', pocketsoap_struct_test, None),
    ( 'base64 test', phalanx_b64_test, None),
    ( 'hexBinary', hexbin_test, None),
    ( 'big base64 test', someones_b64_test, None),
    ( 'echovoid', phalanx_void_test, None),
    ( 'simple2struct', wm_simple_as_struct_test, None),
    ( 'multipart', multipart_test,
        'Content-type: multipart/related; boundary="sep"\r\n' ),
    ( 'int href test', phalanx_int_href_test, None),
    ( 'apache float', apache_float_test, None),
    ( 'bad href test', phalanx_badhref_test, None),
    ( 'bad type attr on struct', phalanx_badstructtype_test, None),
)
