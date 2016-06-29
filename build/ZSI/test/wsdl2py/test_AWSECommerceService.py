#!/usr/bin/env python

############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite, TestException
from ZSI.schema import ElementDeclaration, GED
from ZSI import ParsedSoap

"""
Unittest for contacting the Amazon ECommerce Service

WSDL: 

"""
# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(AmazonTestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(AmazonTestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(AmazonTestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(AmazonTestCase, 'test_'))
    return suite


class AmazonTestCase(ServiceTestCase):
    """Test case for Amazon ECommerce Web service
    """
    name = "test_AWSECommerceService"
    client_file_name = "AWSECommerceService_services.py"
    types_file_name  = "AWSECommerceService_services_types.py"
    server_file_name = "AWSECommerceService_services_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')
        self.wsdl2py_args.append('--lazy')

    def test_local_bug_1525567(self):
        element = GED("http://webservices.amazon.com/AWSECommerceService/2006-09-18", 'Items')
        # Make sure this is a GED
        self.failUnless(isinstance(element, ElementDeclaration), '"%s" not a GED' %element)
    
    def test_local_parse_ItemSearch(self):
        msg = self.client_module.ItemSearchResponseMsg()
        ps = ParsedSoap(ItemSearchResponseMsg)
        response = ps.Parse(msg.typecode)
        response.OperationRequest.Arguments
        for i in response.OperationRequest.Arguments.Argument: 
             i.get_attribute_Name()
             i.get_attribute_Value()

        for i in response.OperationRequest.HTTPHeaders.Header or []:
             i.get_attribute_Name()
             i.get_attribute_Value()
             
        response.OperationRequest.RequestId
        response.OperationRequest.RequestProcessingTime
        for its in response.Items:
            self.failUnless(its.TotalResults == 55, '')
            self.failUnless(its.TotalPages == 6, '')
            for it in its.Item:
                it.ASIN; 
                it.Accessories; 
                #it.AlternateVersions; 
                it.BrowseNodes
                #it.Collections; 
                it.CustomerReviews ;it.DetailPageURL
                it.EditorialReviews; it.Errors; it.ImageSets; it.ItemAttributes
                it.LargeImage; it.ListmaniaLists; it.MediumImage; it.MerchantItemAttributes
                it.OfferSummary; it.Offers; 
                #it.ParentASIN; 
                it.SalesRank; it.SearchInside
                it.SimilarProducts; it.SmallImage; it.Subjects; it.Tracks;


    def test_net_ItemSearch(self):
        loc = self.client_module.AWSECommerceServiceLocator()
        port = loc.getAWSECommerceServicePortType(**self.getPortKWArgs())

        msg = self.client_module.ItemSearchRequestMsg()
        msg.SubscriptionId = '0HP1WHME000749APYWR2'
        request = msg.new_Request()
        msg.Request = [request]

        # request
        request.ItemPage = 1
        request.SearchIndex = "Books"
        request.Keywords = 'Tamerlane'
        request.ResponseGroup = ['Medium',]

        response = port.ItemSearch(msg)

        response.OperationRequest
        self.failUnless(response.OperationRequest.Errors is None, 'ecommerce site reported errors')

        response.OperationRequest.Arguments
        for i in response.OperationRequest.Arguments.Argument: 
             i.get_attribute_Name()
             i.get_attribute_Value()

        for i in response.OperationRequest.HTTPHeaders.Header or []:
             i.get_attribute_Name()
             i.get_attribute_Value()
             
        response.OperationRequest.RequestId
        response.OperationRequest.RequestProcessingTime
        for its in response.Items:
            for it in its.Item:
                it.ASIN; 
                it.Accessories; 
                #it.AlternateVersions; 
                it.BrowseNodes
                #it.Collections; 
                it.CustomerReviews ;it.DetailPageURL
                it.EditorialReviews; it.Errors; it.ImageSets; it.ItemAttributes
                it.LargeImage; it.ListmaniaLists; it.MediumImage; it.MerchantItemAttributes
                it.OfferSummary; it.Offers; 
                #it.ParentASIN; 
                it.SalesRank; it.SearchInside
                it.SimilarProducts; it.SmallImage; it.Subjects; it.Tracks;
                it.VariationSummary; it.Variations


ItemSearchResponseMsg="""<?xml version="1.0" encoding="UTF-8"?>
<SOAP-ENV:Envelope xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/" 
xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/" 
xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" 
xmlns:xsd="http://www.w3.org/2001/XMLSchema"><SOAP-ENV:Body><ItemSearchResponse xmlns="http://webservices.amazon.com/AWSECommerceService/2006-09-18"><OperationRequest><HTTPHeaders><Header Name="UserAgent"></Header></HTTPHeaders><RequestId>1C167XDF9BX253MEYAF2</RequestId><Arguments><Argument Name="Service" Value="AWSECommerceService"></Argument></Arguments><RequestProcessingTime>0.987805128097534</RequestProcessingTime></OperationRequest><Items><Request><IsValid>True</IsValid><ItemSearchRequest><ItemPage>1</ItemPage><Keywords>Tamerlane</Keywords><ResponseGroup>Medium</ResponseGroup><SearchIndex>Books</SearchIndex></ItemSearchRequest></Request><TotalResults>55</TotalResults><TotalPages>6</TotalPages><Item><ASIN>030681465X</ASIN><DetailPageURL>http://www.amazon.com/gp/redirect.html%3FASIN=030681465X%26tag=ws%26lcode=sp1%26cID=2025%26ccmID=165953%26location=/o/ASIN/030681465X%253FSubscriptionId=0HP1WHME000749APYWR2</DetailPageURL><SalesRank>135340</SalesRank><SmallImage><URL>http://images.amazon.com/images/P/030681465X.01._SCTHUMBZZZ_V66860320_.jpg</URL><Height Units="pixels">75</Height><Width Units="pixels">50</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/030681465X.01._SCMZZZZZZZ_V66860320_.jpg</URL><Height Units="pixels">160</Height><Width Units="pixels">108</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/030681465X.01._SCLZZZZZZZ_V66860320_.jpg</URL><Height Units="pixels">500</Height><Width Units="pixels">337</Width></LargeImage><ImageSets><ImageSet Category="primary"><SmallImage><URL>http://images.amazon.com/images/P/030681465X.01._SCTHUMBZZZ_V66860320_.jpg</URL><Height Units="pixels">75</Height><Width Units="pixels">50</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/030681465X.01._SCMZZZZZZZ_V66860320_.jpg</URL><Height Units="pixels">160</Height><Width Units="pixels">108</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/030681465X.01._SCLZZZZZZZ_V66860320_.jpg</URL><Height Units="pixels">500</Height><Width Units="pixels">337</Width></LargeImage></ImageSet></ImageSets><ItemAttributes><Author>Justin Marozzi</Author><Binding>Hardcover</Binding><DeweyDecimalNumber>920</DeweyDecimalNumber><EAN>9780306814655</EAN><Edition>New Ed</Edition><ISBN>030681465X</ISBN><Label>Da Capo Press</Label><ListPrice><Amount>2695</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$26.95</FormattedPrice></ListPrice><Manufacturer>Da Capo Press</Manufacturer><NumberOfItems>1</NumberOfItems><NumberOfPages>368</NumberOfPages><PackageDimensions><Height Units="hundredths-inches">150</Height><Length Units="hundredths-inches">904</Length><Weight Units="hundredths-pounds">174</Weight><Width Units="hundredths-inches">640</Width></PackageDimensions><ProductGroup>Book</ProductGroup><PublicationDate>2006-02-22</PublicationDate><Publisher>Da Capo Press</Publisher><Studio>Da Capo Press</Studio><Title>Tamerlane: Sword of Islam, Conqueror of the World</Title></ItemAttributes><OfferSummary><LowestNewPrice><Amount>831</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$8.31</FormattedPrice></LowestNewPrice><LowestUsedPrice><Amount>832</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$8.32</FormattedPrice></LowestUsedPrice><TotalNew>43</TotalNew><TotalUsed>26</TotalUsed><TotalCollectible>0</TotalCollectible><TotalRefurbished>0</TotalRefurbished></OfferSummary><EditorialReviews><EditorialReview><Source>Book Description</Source><Content>A powerful account of the life of Tamerlane the Great (1336-1405), the last great Mongol conqueror of Central Asia, ruler of a vast empire, and one of history's most brutal tyrants &lt;P&gt; Tamerlane, aka Temur-the Mongol successor to Genghis Khan-ranks with Alexander the Great as one of the world's great conquerors, yet the details of his life are scarcely known in the West. Born in obscurity and poverty, he rose to become a fierce tribal leader, and with that his dominion and power grew with astonishing speed. He blazed through Asia, razing cities to the ground. He tortured conquered inhabitants without mercy, sometimes ordering them buried alive, at other times decapitating them. Over the ruins of conquered Baghdad, Tamerlane had his soldiers erect a pyramid of 90,000 enemy heads. As he and his armies swept through Central Asia, sacking, and then rebuilding cities, Tamerlane gradually imposed an iron rule and a refined culture over a vast territory-from the steppes of Asia to the Syrian coastline. &lt;P&gt; Justin Marozzi traveled in the footsteps of this fearsome emperor of Samarkand (modern-day Uzbekistan) to write this book, which is part history, part travelogue. He carefully follows the path of this infamous and enigmatic conqueror, recounting the history and the story of this cruel, cultivated, and indomitable warrior.</Content></EditorialReview></EditorialReviews></Item><Item><ASIN>1853141046</ASIN><DetailPageURL>http://www.amazon.com/gp/redirect.html%3FASIN=1853141046%26tag=ws%26lcode=sp1%26cID=2025%26ccmID=165953%26location=/o/ASIN/1853141046%253FSubscriptionId=0HP1WHME000749APYWR2</DetailPageURL><SalesRank>366445</SalesRank><ItemAttributes><Author>David Nicolle</Author><Author>Richard Hook</Author><Binding>Hardcover</Binding><DeweyDecimalNumber>950.20922</DeweyDecimalNumber><EAN>9781853141041</EAN><ISBN>1853141046</ISBN><Label>Firebird</Label><ListPrice><Amount>2495</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$24.95</FormattedPrice></ListPrice><Manufacturer>Firebird</Manufacturer><NumberOfItems>1</NumberOfItems><NumberOfPages>208</NumberOfPages><PackageDimensions><Height Units="hundredths-inches">1000</Height><Length Units="hundredths-inches">75</Length><Weight Units="hundredths-pounds">160</Weight><Width Units="hundredths-inches">775</Width></PackageDimensions><ProductGroup>Book</ProductGroup><PublicationDate>1990-09</PublicationDate><Publisher>Firebird</Publisher><Studio>Firebird</Studio><Title>The Mongol Warlords: Ghengis Khan, Kublai Khan, Hulegu, Tamerlane (Heroes &amp; Warriors)</Title></ItemAttributes><OfferSummary><LowestUsedPrice><Amount>1095</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$10.95</FormattedPrice></LowestUsedPrice><TotalNew>0</TotalNew><TotalUsed>3</TotalUsed><TotalCollectible>0</TotalCollectible><TotalRefurbished>0</TotalRefurbished></OfferSummary></Item><Item><ASIN>0521633842</ASIN><DetailPageURL>http://www.amazon.com/gp/redirect.html%3FASIN=0521633842%26tag=ws%26lcode=sp1%26cID=2025%26ccmID=165953%26location=/o/ASIN/0521633842%253FSubscriptionId=0HP1WHME000749APYWR2</DetailPageURL><SalesRank>429712</SalesRank><SmallImage><URL>http://images.amazon.com/images/P/0521633842.01._SCTHUMBZZZ_V1114821525_.jpg</URL><Height Units="pixels">60</Height><Width Units="pixels">39</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/0521633842.01._SCMZZZZZZZ_V1114821525_.jpg</URL><Height Units="pixels">140</Height><Width Units="pixels">90</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/0521633842.01._SCLZZZZZZZ_V1114821525_.jpg</URL><Height Units="pixels">475</Height><Width Units="pixels">306</Width></LargeImage><ImageSets><ImageSet Category="primary"><SmallImage><URL>http://images.amazon.com/images/P/0521633842.01._SCTHUMBZZZ_V1114821525_.jpg</URL><Height Units="pixels">60</Height><Width Units="pixels">39</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/0521633842.01._SCMZZZZZZZ_V1114821525_.jpg</URL><Height Units="pixels">140</Height><Width Units="pixels">90</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/0521633842.01._SCLZZZZZZZ_V1114821525_.jpg</URL><Height Units="pixels">475</Height><Width Units="pixels">306</Width></LargeImage></ImageSet></ImageSets><ItemAttributes><Author>Beatrice Forbes Manz</Author><Binding>Paperback</Binding><DeweyDecimalNumber>950.2</DeweyDecimalNumber><EAN>9780521633840</EAN><Edition>Reprint</Edition><ISBN>0521633842</ISBN><Label>Cambridge University Press</Label><ListPrice><Amount>2299</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$22.99</FormattedPrice></ListPrice><Manufacturer>Cambridge University Press</Manufacturer><NumberOfItems>1</NumberOfItems><NumberOfPages>246</NumberOfPages><PackageDimensions><Height Units="hundredths-inches">67</Height><Length Units="hundredths-inches">850</Length><Weight Units="hundredths-pounds">74</Weight><Width Units="hundredths-inches">562</Width></PackageDimensions><ProductGroup>Book</ProductGroup><PublicationDate>1999-03-28</PublicationDate><Publisher>Cambridge University Press</Publisher><Studio>Cambridge University Press</Studio><Title>The Rise and Rule of Tamerlane (Canto original series)</Title></ItemAttributes><OfferSummary><LowestNewPrice><Amount>1298</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$12.98</FormattedPrice></LowestNewPrice><LowestUsedPrice><Amount>1054</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$10.54</FormattedPrice></LowestUsedPrice><LowestCollectiblePrice><Amount>2299</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$22.99</FormattedPrice></LowestCollectiblePrice><TotalNew>25</TotalNew><TotalUsed>21</TotalUsed><TotalCollectible>1</TotalCollectible><TotalRefurbished>0</TotalRefurbished></OfferSummary><EditorialReviews><EditorialReview><Source>Book Description</Source><Content>This is the first serious study of Tamerlane, the great nomad conqueror who rose to power in 1370 on the ruins of the Mongol Empire and led his armies on campaigns of unprecedented scope, ranging from Moscow to Delhi.  As the last nomad ruler to unite the steppe regions of Eurasia, Tamerlane marks the transition from the era of nomad conquest and rule to the modern ascendency of the settled world.</Content></EditorialReview></EditorialReviews></Item><Item><ASIN>1885221770</ASIN><DetailPageURL>http://www.amazon.com/gp/redirect.html%3FASIN=1885221770%26tag=ws%26lcode=sp1%26cID=2025%26ccmID=165953%26location=/o/ASIN/1885221770%253FSubscriptionId=0HP1WHME000749APYWR2</DetailPageURL><SalesRank>474520</SalesRank><SmallImage><URL>http://images.amazon.com/images/P/1885221770.01._SCTHUMBZZZ_V1056534986_.jpg</URL><Height Units="pixels">60</Height><Width Units="pixels">40</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/1885221770.01._SCMZZZZZZZ_V1056534986_.jpg</URL><Height Units="pixels">140</Height><Width Units="pixels">93</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/1885221770.01._SCLZZZZZZZ_V1056534986_.jpg</URL><Height Units="pixels">475</Height><Width Units="pixels">317</Width></LargeImage><ImageSets><ImageSet Category="primary"><SmallImage><URL>http://images.amazon.com/images/P/1885221770.01._SCTHUMBZZZ_V1056534986_.jpg</URL><Height Units="pixels">60</Height><Width Units="pixels">40</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/1885221770.01._SCMZZZZZZZ_V1056534986_.jpg</URL><Height Units="pixels">140</Height><Width Units="pixels">93</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/1885221770.01._SCLZZZZZZZ_V1056534986_.jpg</URL><Height Units="pixels">475</Height><Width Units="pixels">317</Width></LargeImage></ImageSet></ImageSets><ItemAttributes><Author>Roy Stier</Author><Binding>Paperback</Binding><EAN>9781885221773</EAN><ISBN>1885221770</ISBN><Label>Bookpartners</Label><ListPrice><Amount>1695</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$16.95</FormattedPrice></ListPrice><Manufacturer>Bookpartners</Manufacturer><NumberOfItems>1</NumberOfItems><NumberOfPages>304</NumberOfPages><PackageDimensions><Height Units="hundredths-inches">79</Height><Length Units="hundredths-inches">900</Length><Weight Units="hundredths-pounds">92</Weight><Width Units="hundredths-inches">606</Width></PackageDimensions><ProductGroup>Book</ProductGroup><PublicationDate>1998-09</PublicationDate><Publisher>Bookpartners</Publisher><Studio>Bookpartners</Studio><Title>Tamerlane: The Ultimate Warrior</Title></ItemAttributes><OfferSummary><LowestUsedPrice><Amount>975</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$9.75</FormattedPrice></LowestUsedPrice><LowestCollectiblePrice><Amount>2295</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$22.95</FormattedPrice></LowestCollectiblePrice><TotalNew>0</TotalNew><TotalUsed>3</TotalUsed><TotalCollectible>1</TotalCollectible><TotalRefurbished>0</TotalRefurbished></OfferSummary><EditorialReviews><EditorialReview><Source>Book Description</Source><Content>From humble beginnings, Tamerlane, the ancient Turki-Mongol conqueror, rose to become the scourge of his time and changed the course of history.   &lt;P&gt;The name Tamerlane runs the gamut of human emotions, evoking in many a revulsion for the devil incarnate, in others, an appreciation for the benefactor of millions.   &lt;P&gt;By using accounts from Tamerlane's detractors and his admirers, Roy Stier has captured an amazing story that gives credence to the old adage, "truth is stranger than fiction."   &lt;P&gt;Tamerlane, the Ultimate Warrior is presented as a fascinating series of events and captures the reader in the first comprehensive view of this historical figure who dominated Asia and made Europe tremble.</Content></EditorialReview></EditorialReviews></Item><Item><ASIN>0850459494</ASIN><DetailPageURL>http://www.amazon.com/gp/redirect.html%3FASIN=0850459494%26tag=ws%26lcode=sp1%26cID=2025%26ccmID=165953%26location=/o/ASIN/0850459494%253FSubscriptionId=0HP1WHME000749APYWR2</DetailPageURL><SalesRank>480359</SalesRank><SmallImage><URL>http://images.amazon.com/images/P/0850459494.01._SCTHUMBZZZ_V1128022797_.jpg</URL><Height Units="pixels">75</Height><Width Units="pixels">56</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/0850459494.01._SCMZZZZZZZ_V1128022797_.jpg</URL><Height Units="pixels">160</Height><Width Units="pixels">119</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/0850459494.01._SCLZZZZZZZ_V1128022797_.jpg</URL><Height Units="pixels">500</Height><Width Units="pixels">372</Width></LargeImage><ImageSets><ImageSet Category="primary"><SmallImage><URL>http://images.amazon.com/images/P/0850459494.01._SCTHUMBZZZ_V1128022797_.jpg</URL><Height Units="pixels">75</Height><Width Units="pixels">56</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/0850459494.01._SCMZZZZZZZ_V1128022797_.jpg</URL><Height Units="pixels">160</Height><Width Units="pixels">119</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/0850459494.01._SCLZZZZZZZ_V1128022797_.jpg</URL><Height Units="pixels">500</Height><Width Units="pixels">372</Width></LargeImage></ImageSet></ImageSets><ItemAttributes><Author>David Nicolle</Author><Binding>Paperback</Binding><Creator Role="Illustrator">Angus Mcbride</Creator><EAN>9780850459494</EAN><ISBN>0850459494</ISBN><Label>Osprey Publishing</Label><ListPrice><Amount>1595</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$15.95</FormattedPrice></ListPrice><Manufacturer>Osprey Publishing</Manufacturer><NumberOfItems>1</NumberOfItems><NumberOfPages>48</NumberOfPages><PackageDimensions><Height Units="hundredths-inches">15</Height><Length Units="hundredths-inches">978</Length><Weight Units="hundredths-pounds">36</Weight><Width Units="hundredths-inches">722</Width></PackageDimensions><ProductGroup>Book</ProductGroup><PublicationDate>1990-07-26</PublicationDate><Publisher>Osprey Publishing</Publisher><ReleaseDate>1990-07-26</ReleaseDate><Studio>Osprey Publishing</Studio><Title>The Age of Tamerlane (Men-at-Arms)</Title></ItemAttributes><OfferSummary><LowestNewPrice><Amount>920</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$9.20</FormattedPrice></LowestNewPrice><LowestUsedPrice><Amount>815</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$8.15</FormattedPrice></LowestUsedPrice><TotalNew>12</TotalNew><TotalUsed>9</TotalUsed><TotalCollectible>0</TotalCollectible><TotalRefurbished>0</TotalRefurbished></OfferSummary><EditorialReviews><EditorialReview><Source>Book Description</Source><Content>Tamerlane or Timur-i-Lenk ('Timur the Lame') is one of the most extraordinary conquerors in history. In the late 14th century his armies seized huge territories from the borders of Mongolia to Palestine and Anatolia. His passage was marked by massacres that outdid even those of the Mongols for sheer savagery. Timur's career was unequalled since Alexander the Great in terms of constant battlefield success. Only in his youth, while recovering his family estates south of Samarqand, did he face occasional defeat. This title tells the remarkable story of Timur and details the organisation, tactics, arms and armour of his all-conquering army.</Content></EditorialReview></EditorialReviews></Item><Item><ASIN>1851684573</ASIN><DetailPageURL>http://www.amazon.com/gp/redirect.html%3FASIN=1851684573%26tag=ws%26lcode=sp1%26cID=2025%26ccmID=165953%26location=/o/ASIN/1851684573%253FSubscriptionId=0HP1WHME000749APYWR2</DetailPageURL><SalesRank>600493</SalesRank><SmallImage><URL>http://images.amazon.com/images/P/1851684573.01._SCTHUMBZZZ_V63318467_.jpg</URL><Height Units="pixels">75</Height><Width Units="pixels">49</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/1851684573.01._SCMZZZZZZZ_V63318467_.jpg</URL><Height Units="pixels">160</Height><Width Units="pixels">104</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/1851684573.01._SCLZZZZZZZ_V63318467_.jpg</URL><Height Units="pixels">500</Height><Width Units="pixels">324</Width></LargeImage><ImageSets><ImageSet Category="primary"><SmallImage><URL>http://images.amazon.com/images/P/1851684573.01._SCTHUMBZZZ_V63318467_.jpg</URL><Height Units="pixels">75</Height><Width Units="pixels">49</Width></SmallImage><MediumImage><URL>http://images.amazon.com/images/P/1851684573.01._SCMZZZZZZZ_V63318467_.jpg</URL><Height Units="pixels">160</Height><Width Units="pixels">104</Width></MediumImage><LargeImage><URL>http://images.amazon.com/images/P/1851684573.01._SCLZZZZZZZ_V63318467_.jpg</URL><Height Units="pixels">500</Height><Width Units="pixels">324</Width></LargeImage></ImageSet></ImageSets><ItemAttributes><Author>Robert Rand</Author><Binding>Paperback</Binding><DeweyDecimalNumber>320</DeweyDecimalNumber><EAN>9781851684571</EAN><ISBN>1851684573</ISBN><Label>Oneworld Publications</Label><ListPrice><Amount>1495</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$14.95</FormattedPrice></ListPrice><Manufacturer>Oneworld Publications</Manufacturer><NumberOfItems>1</NumberOfItems><NumberOfPages>224</NumberOfPages><PackageDimensions><Height Units="hundredths-inches">66</Height><Length Units="hundredths-inches">784</Length><Weight Units="hundredths-pounds">54</Weight><Width Units="hundredths-inches">512</Width></PackageDimensions><ProductGroup>Book</ProductGroup><PublicationDate>2006-09-11</PublicationDate><Publisher>Oneworld Publications</Publisher><Studio>Oneworld Publications</Studio><Title>Tamerlane's Children: Dispatches from Contemporary Uzbekistan</Title></ItemAttributes><OfferSummary><LowestNewPrice><Amount>915</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$9.15</FormattedPrice></LowestNewPrice><LowestUsedPrice><Amount>950</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$9.50</FormattedPrice></LowestUsedPrice><TotalNew>31</TotalNew><TotalUsed>6</TotalUsed><TotalCollectible>0</TotalCollectible><TotalRefurbished>0</TotalRefurbished></OfferSummary><EditorialReviews><EditorialReview><Source>Book Description</Source><Content>In the central park of Tashkent, in a place the Uzbeks call the square, a magnificent statue of a mounted warrior dominates the surroundings - Tamerlane - national hero of post-Soviet Uzbekistan. And yet how does this 14th century conqueror reflect one of the world's most diverse and politically intriguing countries?Having spent three years in the region, renowned journalist Robert Rand seeks to answer this question, covering an assortment of fascinating topics, ranging from the effect of 9/11 to the clash of culture in Uzbek pop music.  Overflowing with charming anecdotes and loveable personalities, Rand gives the reader a real sense of the country's confused identity and  the challenges which it and its people will face in generations to come.</Content></EditorialReview></EditorialReviews></Item><Item><ASIN>B000GKT9AQ</ASIN><DetailPageURL>http://www.amazon.com/gp/redirect.html%3FASIN=B000GKT9AQ%26tag=ws%26lcode=sp1%26cID=2025%26ccmID=165953%26location=/o/ASIN/B000GKT9AQ%253FSubscriptionId=0HP1WHME000749APYWR2</DetailPageURL><SalesRank>619234</SalesRank><ItemAttributes><Author>Harold Lamb</Author><Binding>Hardcover</Binding><Label>Garden City, N.Y., U.S.A.: Doubleday &amp; Company, Inc.</Label><Manufacturer>Garden City, N.Y., U.S.A.: Doubleday &amp; Company, Inc.</Manufacturer><ProductGroup>Book</ProductGroup><ProductTypeName>BOOKS_1973_AND_LATER</ProductTypeName><PublicationDate>1941</PublicationDate><Publisher>Garden City, N.Y., U.S.A.: Doubleday &amp; Company, Inc.</Publisher><Studio>Garden City, N.Y., U.S.A.: Doubleday &amp; Company, Inc.</Studio><Title>Earth Shakers: The: The March of the Barbarians and Tamerlane</Title></ItemAttributes><OfferSummary><LowestUsedPrice><Amount>850</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$8.50</FormattedPrice></LowestUsedPrice><TotalNew>0</TotalNew><TotalUsed>2</TotalUsed><TotalCollectible>0</TotalCollectible><TotalRefurbished>0</TotalRefurbished></OfferSummary></Item><Item><ASIN>B00087SKA2</ASIN><DetailPageURL>http://www.amazon.com/gp/redirect.html%3FASIN=B00087SKA2%26tag=ws%26lcode=sp1%26cID=2025%26ccmID=165953%26location=/o/ASIN/B00087SKA2%253FSubscriptionId=0HP1WHME000749APYWR2</DetailPageURL><SalesRank>855185</SalesRank><ItemAttributes><Author>Harold Lamb</Author><Binding>Unknown Binding</Binding><Label>Bantam Books</Label><Manufacturer>Bantam Books</Manufacturer><NumberOfPages>216</NumberOfPages><ProductGroup>Book</ProductGroup><PublicationDate>1955</PublicationDate><Publisher>Bantam Books</Publisher><Studio>Bantam Books</Studio><Title>Tamerlane: Conqueror of the earth</Title></ItemAttributes><OfferSummary><LowestUsedPrice><Amount>650</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$6.50</FormattedPrice></LowestUsedPrice><TotalNew>0</TotalNew><TotalUsed>5</TotalUsed><TotalCollectible>0</TotalCollectible><TotalRefurbished>0</TotalRefurbished></OfferSummary></Item><Item><ASIN>0689304463</ASIN><DetailPageURL>http://www.amazon.com/gp/redirect.html%3FASIN=0689304463%26tag=ws%26lcode=sp1%26cID=2025%26ccmID=165953%26location=/o/ASIN/0689304463%253FSubscriptionId=0HP1WHME000749APYWR2</DetailPageURL><SalesRank>1028876</SalesRank><ItemAttributes><Author>Barbara Corcoran</Author><Binding>Hardcover</Binding><Creator Role="Illustrator">Charles Robinson</Creator><EAN>9780689304460</EAN><Edition>[1st ed.]</Edition><ISBN>0689304463</ISBN><Label>Macmillan Pub Co</Label><ListPrice><Amount>695</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$6.95</FormattedPrice></ListPrice><Manufacturer>Macmillan Pub Co</Manufacturer><NumberOfItems>1</NumberOfItems><NumberOfPages>152</NumberOfPages><ProductGroup>Book</ProductGroup><PublicationDate>1975-02</PublicationDate><Publisher>Macmillan Pub Co</Publisher><ReadingLevel>Young Adult</ReadingLevel><Studio>Macmillan Pub Co</Studio><Title>Meet Me at Tamerlane's Tomb</Title></ItemAttributes><OfferSummary><LowestUsedPrice><Amount>30</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$0.30</FormattedPrice></LowestUsedPrice><LowestCollectiblePrice><Amount>1000</Amount><CurrencyCode>USD</CurrencyCode><FormattedPrice>$10.00</FormattedPrice></LowestCollectiblePrice><TotalNew>0</TotalNew><TotalUsed>20</TotalUsed><TotalCollectible>2</TotalCollectible><TotalRefurbished>0</TotalRefurbished></OfferSummary></Item></Items></ItemSearchResponse></SOAP-ENV:Body></SOAP-ENV:Envelope>"""

if __name__ == '__main__':
    main()
