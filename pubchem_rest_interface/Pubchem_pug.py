#!/usr/bin/python
# -*- coding: utf-8 -*-


#from builtins import str
import requests
#from requests.compat import quote
import gzip
import time
import urllib.request
import xml.dom.minidom


def DownloadCIDs(cids):

    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"+(",".join([str(x) for x in cids]))+"/SDF"
    req = requests.get(url)
    if req.status_code == requests.codes.ok:
        return req.text
    else:
        print("failed to download sdf from pubchem. status: {}, result: {}".format(req.status_code,req.text))

def SimilaritySearch(smilesQuery,threshold,maxResults):
    #print("SimilaritySearch. query: {}, threshold: {}, maxResults: {}".format(smilesQuery,threshold,maxResults))
    #url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/"+
     #       quote(smilesQuery.rstrip(),safe='')+"/cids/TXT?Threshold="+str(threshold)+"&MaxRecords="+str(maxResults))
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/"+
            "cids/TXT?Threshold="+str(threshold)+"&MaxRecords="+str(maxResults))

    print("SimilaritySearch. url: "+url)
    req = requests.post(url,data={'smiles':smilesQuery})
    if req.status_code == requests.codes.ok:
        return req.text
    else:
        print("failed to run similarity search on pubchem. status: {}, result: {}".format(req.status_code,req.text))
        raise Exception("failed to run similarity search on pubchem. status: {}, result: {}".format(req.status_code,req.text))

def makeQueryXML(cid_list, dl_format):
    xml_template = """
    <PCT-Data>
      <PCT-Data_input>
        <PCT-InputData>
          <PCT-InputData_download>
            <PCT-Download>
              <PCT-Download_uids>
                <PCT-QueryUids>
                  <PCT-QueryUids_ids>
                    <PCT-ID-List>
                      <PCT-ID-List_db>pccompound</PCT-ID-List_db>
                      <PCT-ID-List_uids>
                        {xml_query_data}
                      </PCT-ID-List_uids>
                    </PCT-ID-List>
                  </PCT-QueryUids_ids>
                </PCT-QueryUids>
              </PCT-Download_uids>
              <PCT-Download_format value="{dl_format}"/>
              <PCT-Download_compression value="gzip"/>
              <PCT-Download_use-3d value="false"/>
            </PCT-Download>
          </PCT-InputData_download>
        </PCT-InputData>
      </PCT-Data_input>
    </PCT-Data>"""

    xml_query_list = [ "<PCT-ID-List_uids_E>{}</PCT-ID-List_uids_E>".format(cid) for cid in cid_list ]

    return xml_template.format(xml_query_data="".join(xml_query_list), dl_format=dl_format)

def makePollXML(reqid, action):
    xml_template = """
    <PCT-Data>
      <PCT-Data_input>
        <PCT-InputData>
          <PCT-InputData_request>
            <PCT-Request>
              <PCT-Request_reqid>{reqid}</PCT-Request_reqid>
              <PCT-Request_type value="{action}"/>
            </PCT-Request>
          </PCT-InputData_request>
        </PCT-InputData>
      </PCT-Data_input>
    </PCT-Data>"""

    return xml_template.format(reqid=reqid, action=action)

def checkQueryResult(xml_text):
    result = dict()
    dom = xml.dom.minidom.parseString(xml_text)
    status_node = dom.getElementsByTagName("PCT-Status")
    if len(status_node) == 0:
        raise Exception("PCT-Status node not found")

    status = status_node[0].getAttribute("value")
    if status == "success":
        node = dom.getElementsByTagName("PCT-Download-URL_url")[0]
        result["action"] = "download"
        result["url"] = node.firstChild.data
        return result
    elif status in ["running", "queued"]:
        node = dom.getElementsByTagName("PCT-Waiting_reqid")[0]
        result["action"] = "wait"
        result["reqid"] = node.firstChild.data
        msg_node = dom.getElementsByTagName("PCT-Waiting_message")
        if len(msg_node) == 0:
            result["message"] = None
        else:
            result["message"] = msg_node[0].firstChild.data
        return result
    else:
        msg = "PCT-Status value: {}. Should be one of {}\n".format(status, ["success", "running", "queued"])
        msg += "Response:\n" + xml_text
        raise Exception(msg)

def sendPUGQuery(xml_text):
    url = "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"
    req = requests.post(url, data=xml_text)

    if req.status_code != requests.codes.ok:
        msg = "Status code: {}. Expected code: {}\n".format(req.status_code, requests.codes.ok)
        msg += "Response:\n" + req.text
        raise Exception(msg)

    return req.text

def pubchemDownload(cid_list, dl_format, poll_interval=20, timeout=None, verbose=False):
    query_xml = makeQueryXML(cid_list, dl_format)
    query_response = sendPUGQuery(query_xml)
    query_result = checkQueryResult(query_response)

    runtime = 0
    while query_result["action"] == "wait":
        if timeout is not None and runtime >= timeout:
            cancel_xml = makePollXML(query_result["reqid"], "cancel")
            cancel_response = sendPUGQuery(cancel_xml)
            raise Exception("runtime exceeded timeout ({} sec)".format(timeout))
        if verbose:
            print("reqid: " + query_result["reqid"])
            print("message: " + query_result["message"])
            print("runtime: " + str(runtime))
            print("waiting for {} seconds...".format(poll_interval))
        time.sleep(poll_interval)
        runtime += poll_interval
        query_xml = makePollXML(query_result["reqid"], "status")
        query_response = sendPUGQuery(query_xml)
        query_result = checkQueryResult(query_response)

    assert query_result["action"] == "download"
    if verbose:
        print("url: " + query_result["url"])
    req = urllib.request.urlopen(query_result["url"])
    gzip_data = req.file.file.read()
    data = gzip.decompress(gzip_data)

    return data.decode()

#if __name__ == '__main__':
#    DownloadCIDs()
