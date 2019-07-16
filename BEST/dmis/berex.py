"""
.. module:: dmis.berex
    :platform: Unix, linux, Windows
.. moduleauthor:: Minji Jeon <minji.jeon1@gmail.com>

=====================================
Biomedical entity network query API
=====================================

.. note:: Usage: Biomedical entity network query API

>>> from dmis import berex
>>> 
>>> berexQuery = berex.BEReXQuery({"keywordA":["chronic myeloid leukemia"], "keywordB":["ABL1", "imatinib"], "mode":"subnet"})
>>> relevantRelations = getRelevantBioRelations(berexQuery)
>>> 
>>> print(relevantRelations)


"""

from functools import reduce

import json
import urllib.request

#berex query object
class BEReXQuery():
    """
    dmis.BEReXQuery class is basic query object for BEReX API.
    
    """
        
    berexurl = "http://berex.korea.ac.kr/api/"
    
    #query has a mode
    #mode 1 is to get edges about query entities
    #mode 7 is to get enriched gene ontology
        
    def __init__(self, queryObj={"keywordA":[], "keywordB":[], "mode":1}):
        """BEReXQuery
        
        :param queryObj: keywordA (list), keywordB (list), mode (["subnet", "GOTerms"]) dict return.
        
        >>> query = BEReXQuery({"keywordA":["chronic myeloid leukemia", "BCR"], "keywordB":["ABL1", "imatinib"], "mode":"subnet"})
        """
        
        if queryObj["keywordA"] == None:
            queryObj["keywordA"] = [""]
            
        if len(queryObj["keywordA"]) == 0:
            queryObj["keywordA"] = [""]
        
        if ("keywordA" not in queryObj) or (type(queryObj["keywordA"]) is not list) or ("keywordB" not in queryObj) or(type(queryObj["keywordB"]) is not list) or ("mode" not in queryObj):
            print ("Initialize error: invalid query object, query object should contains 'keywordA (list of string)', 'keywordB (list of string)', 'mode (['subnet', 'GOTerms'])'")
            print (queryObj)
        
        for keya in queryObj["keywordA"] :
            if type(keya) is not str :
                print ("Initialize error: invalid keywordA. keywordA should be either None, empty list or list of string")
                print (queryObj["keywordA"])
            
        for keyb in queryObj["keywordB"] :
            if type(keyb) is not str :
                print ("Initialize error: invalid keywordB. keywordB should be list of string")
                print (queryObj["keywordB"])
                
        if queryObj["mode"] not in ["subnet", "GOTerms"] :
            print("Initialize error: invalid search mode. The search mode should be either 'subnet' or 'GOTerms'")
        
        self.keywordA = queryObj["keywordA"]
        self.keywordB = queryObj["keywordB"]
        self.mode = queryObj["mode"]
    
    def setKeywordA (self, keywords):
        """Setting the primary keywords (Keyword A)
        
        :param keyword: primary keywords, which must be a list of str
        
        >>> query.setKeywordA(["cancer"])
        """
       
        for keya in keywords :
            if type(keya) is not str :
                print ("Initialize error : invalid keywordA. keywordA should be list of string")
                print (keywords)
                return
  
        if len(keywords) == 0:
            keywords = [""]
            return
            
        self.keywordA = keywords
        
    def getKeywordA (self):
        """Getting the primary keyword (Keyword A)
        
        :return: keyword A list
        
        >>> keywordA = query.getKeywordA()
        >>> print (keywordA)
        ['chronic myeloid leukemia']
        """
        return self.keywordA
    
    def addKeywordtoA (self, keyword):
        """Adding a keyword to the primary keyword list (Keyword A)
        
        :param keyword: the keyword to be added to the primary keyword list
        
        >>> print (query.getKeywordA())
        ['chronic myeloid leukemia', 'BCR']
        >>> query.addKeywordtoA("EGFR")
        >>> print (query.getKeywordA())
        ['chronic myeloid leukemia', 'BCR', "EGFR']
        """
        self.keywordA.append(keyword)
    
    def removeKeywordfromA(self, keyword):
        """Removing a keyword from the primary keyword list (Keyword A)
        
        :param keyword: the keyword to be removed from the primary keyword list
        
        >>> print (query.getKeywordA())
        ['chronic myeloid leukemia', 'BCR', 'EGFR']
        >>> query.removeKeywordfromA("EGFR")
        >>> print (query.getKeywordA())
        ['chronic myeloid leukemia', 'BCR']
        """
        self.keywordA.remove(keyword)
        
    
    def setKeywordB (self, keywords):
        """Setting the secondary keywords (Keyword B)
        
        :param keywords: the secondary keywords, which must be a list of str
        
        >>> keywordB = ["ABL1", "imatinib"]
        >>> query.setKeywordB(keywordB)
        """
        for keyb in keywords :
            if type(keyb) is not str :
                print ("Initialize error : invalid keywordB. keywordB should be list of string")
                print (keywords)
                return
        
        if type(keywords) is list:
            self.keywordB = keywords
        else :
            print ("Warning! keywords should be list type : " + str(keywords))
    
    def isValid(self):
        if type(self.keywordA) is not list:
            return False
        
        for keya in self.keywordA :
            if type(keya) is not str :
                return False
        
        for keyb in self.keywordB :
            if type(keyb) is not str :
                return False
        
        if len(self.keywordB) == 0:
            return False
        
        #if self.mode != 1 and self.mode != 7:
#            return False

        if self.mode not in ["subnet", "GOTerms"]:
            return False
        
        return True
    
    def getKeywordB (self):
        """Getting the secondary keywords (Keyword B)
        
        :return: list of keyword B string
        
        >>> keywordB = query.getKeywordB()
        >>> print (keywordB)
        ['ABL1', 'imatinib']
        """
        return self.keywordB
    
    def addKeywordtoB (self, keyword):
        """Adding a keyword to the secondary keyword list (Keyword B)
        
        :param keyword: the keyword to be added to the secondary keyword list
        
        >>> print (query.getKeywordB())
        ['ABL1', 'imatinib']
        >>> query.addKeywordtoB("EGFR")
        >>> print (query.getKeywordB())
        ['ABL1', 'imatinib', 'EGFR']
        """
        self.keywordB.append(keyword)
    
    def removeKeywordfromB(self, keyword):
        """Removing a keyword from the secondary keyword list (Keyword B)
        
        :param keyword: the keyword to be removed from the secondary keyword list
        
        >>> print (query.getKeywordB())
        ['ABL1', 'imatinib', 'EGFR']
        >>> query.removeKeywordfromB("EGFR")
        >>> print (query.getKeywordB())
        ['ABL1', 'imatinib']
        """
        self.keywordB.remove(keyword)
    
    def setMode (self, mode):
        """ Setting a searching mode
        
        :param mode: searching mode. This should be either "subnet" or "GOTerms"
        
        >>> query.setMode("subnet")
        """
        self.mode = mode
        
    def getMode (self):
        """ Getting the searching mode
        
        :return: the searching mode
        
        >>> print (query.getMode())
        subnet
        """
        return self.mode
    
    def makeQueryString(self):        
 
        #paramKeywordA = self.keywordA
        paramKeywordA = reduce(lambda x, y: x +"@@" + y , self.keywordA)
        paramKeywordB = reduce(lambda x, y: x +"@@" + y , self.keywordB)
        queryKeywords = paramKeywordB
        
        if paramKeywordA != "":
            queryKeywords = queryKeywords + "@@" + paramKeywordA
        queryKeywords = queryKeywords.replace(" ","%20")
        
        mode = 1
        if self.mode == "subnet":
            mode = 1
        elif self.mode == "GOTerms" :
            mode = 7
        
        strQuery = self.berexurl +"server.php?query=" + queryKeywords +"&mode=" + str(mode)
        return strQuery
    
#get raw data from berex
def getRelevantBioRelations(berexQuery):
    """ Function for retrieval from BOSS
    
    :param berexQuery: BEReXQuery
    
    :return: subnetwork object (mode "subnet") or list of enriched GO terms (mode "GOTerms")
    
    * subnetwork object (list): [BIO_EDGE]
        - BIO_EDGE: {"source":str, "target":str, "dbsource":str, "interaction":str, "reference":str}
    * list of GO terms (object): {"biological process":[str], "molecular function":[str], "cellular component":[str]}
    
    >>> berexQuery = BEReXQuery({"keywordA":["chronic myeloid leukemia"], "keywordB":["ABL1", "imatinib"], "mode":"subnet"})
    >>> relevantRelations = getRelevantBioRelations(bossQuery)
    """
    if not (type(berexQuery) is BEReXQuery):
        print ("query is invalid! please check your query object.")
       
    if not berexQuery.isValid() :
        print ("Query object is invalid. Please check the query")
        print ("Query: ")
        print ("   keywordA: " + str(berexQuery.keywordA))
        print ("   keywordB: " + str(berexQuery.keywordB))
        print ("   mode: " + str(berexQuery.mode))
            
        return {}
    
    Query = berexQuery.makeQueryString()
    
   
    request = urllib.request.Request(Query)
    request.add_header('User-Agent', 'Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1)')
                
    Url = urllib.request.urlopen(request)
    
    print (Query)
    ResultStr = Url.read().decode('utf-8')
    
    
    if berexQuery.getMode() == "subnet":
        Result = makeDataFromBEReXQueryResult(ResultStr)
    elif berexQuery.getMode() == "GOTerms":
        Result = makeGODataFromBEReXQueryResult(ResultStr)
    return Result

#make data for mode 1 (entity search mode)
#resultDataArr is array containing each edge
#get source node name : resultDataArr[i]['source']
#edge contains source, target, dbsource, interaction,and reference

def makeDataFromBEReXQueryResult(resultStr):    
    if resultStr == "":
        return None
    resultStr = resultStr.replace("<img src=\\\"file:imgs/pubmed.png\\\" border=\\\"0\\\" width=\\\"70\\\"/>","")
    resultStr = resultStr.replace("<br/> ","")
    resultObject = json.loads(resultStr, strict=False)
    edges = resultObject['data']['edges']
    
    linesCnt = len(edges)
    resultDataArr = []
    curData = {}
    for i in range(0, linesCnt) :
        edge = edges[i]
 
        curData={"source":edge['source'], "target":edge['target'], "dbsource":edge['dbsource'], "interaction":edge['interaction'],"reference":edge['reference']}
        resultDataArr.append(curData)
    
    return resultDataArr

#make data for mode 7 (go term search mode)
#There are three types of ontology, biological process, molecular function, cellular component
#get biological process: resultData['biological process'] 
def makeGODataFromBEReXQueryResult(resultStr):    
    if resultStr == "":
        return None
    resultObject = json.loads(resultStr, strict=False)
    
    bplist = []
    mflist = []
    cclist = []
    
    for i in range(0, 3):
        if resultObject[i]['mode']=="bp":
            bplist = resultObject[i]['values']
        elif resultObject[i]['mode']=="mf":
            mflist = resultObject[i]['values']
        elif resultObject[i]['mode']=="cc":
            cclist = resultObject[i]['values']
    
    resultData = {"biological process":bplist, "molecular function":mflist, "cellular component":cclist}
    return resultData