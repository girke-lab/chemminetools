//
// To use this:
//  wsdl Echo.wsdl
//  csc /t:library EchoService.cs
//  csc Echo.cs /r:EchoService.dll

using System;

class Echo {
  public static void Main(string[] args) {
   
   EchoService server = new EchoService();

   string request = "Test C#";

   server.Echo(ref request);

   Console.WriteLine("Return: " + request);
  }
}