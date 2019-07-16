import best

bestQuery = best.BESTQuery({"keywordA":["cancer"], "keywordB":["BRCA1", "EGFR"], "filterObjectName":"", "topN":30})
relevantEntities = best.getRelevantBioEntities(bestQuery)
print(relevantEntities)
