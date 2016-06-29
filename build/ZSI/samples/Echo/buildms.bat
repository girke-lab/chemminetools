wsdl Echo.wsdl
csc /t:library EchoService.cs
csc Echo.cs /r:EchoService.dll