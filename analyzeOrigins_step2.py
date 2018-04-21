
# takes in (1) the metadata, and (2) the output from the first step script

import sys

# get taxon of each sample
taxonDict = {}
with open(sys.argv[1], 'r') as infile:
    for line in infile:
        sample, group = line.strip().split()
        if group not in taxonDict:
            taxonDict[group] = []
        taxonDict[group].append(sample)

# read in the gene list for reference
geneList = []
with open(sys.argv[2], 'r') as infile:
    infile.readline()
    for line in infile:
        newline = line.strip().split()
        geneList.append(newline)

# read in the results table
sampleDict = {}
with open(sys.argv[3], 'r') as infile:
    sampleIDs = infile.readline().strip().split()
    for s in sampleIDs:
        sampleDict[s] = []
    numTranscripts = 0
    for line in infile:
        numTranscripts += 1
        newline = line.strip().split()
        for d in range(len(newline)):
            sampleDict[sampleIDs[d]].append(newline[d])

# reading in a list of genes that look suspicious/cool by looking at their box plots
coolDict = {}
with open(sys.argv[4], 'r') as infile:
    for line in infile:
        coolDict[line.strip()] = 0
            

            
# now ask questions!

# first. when the domesticated splice pattern is not found in wild annuus, where IS it found?
otherSpecies = ["deblis", "petiolaris", "paradoxus", "deserticola", "praecox", "anomalus", "argophyllus"]
totalWilds = len(taxonDict["wild"])
totalDoms = len(taxonDict["domesticated"])
totalLandraces = len(taxonDict["landraces"])
for transcript in range(numTranscripts):
    geneId = geneList[transcript]

    # checking in landraces
    countDomTypes = 0
    numWithExpression = 0
    countWildTypes = 0
    sampsWithDomType = []
    for sample in taxonDict["landraces"]:
        version = sampleDict[sample][transcript]
        if version != "NA":
            numWithExpression += 1
        if version == "domesticated":
            countDomTypes += 1
            sampsWithDomType.append(sample)
        elif version == "wild":
            countWildTypes += 1
    if numWithExpression > 0:
        #print "landrace:", countDomTypes, numWithExpression, geneId, sampsWithDomType
        print "landrace:", countWildTypes, numWithExpression, geneId, sampsWithDomType
        
    # checking in wilds
    numWithExpression = 0
    countDomTypes = 0 
    countWildTypes = 0
    for sample in taxonDict["wild"]:
        version = sampleDict[sample][transcript]
        if version != "NA":
            numWithExpression += 1     
        if version == "domesticated":
            countDomTypes += 1
        elif version == "wild":
            countWildTypes += 1
    if numWithExpression > 0:
#        print "wilds", countDomTypes, numWithExpression, geneId
        print "wilds", countWildTypes, numWithExpression, geneId
        
    # new way, giving only a list of genes that look good in the wilds, seeing if in the other species 
    if geneId[0] in coolDict:
        for taxon in otherSpecies:
            for sample in taxonDict[taxon]:
                if sample in sampleDict:
                    theirVersion = sampleDict[sample][transcript]
                    if theirVersion == "domesticated":
                        print "gene", geneId[0], "not found in wilds, but taxon", taxon, "has the domesticated version"





                        
