
# input:
# 1. list of transcripts, and versions
# 2. RSEM table

import sys

minTPM = 0.25

outputType = sys.argv[3]

# find out which isoforms to look at
transcriptDict = {}
with open(sys.argv[1], "r") as infile:
    infile.readline() # header
    for line in infile:
        newline = line.strip().split()
        gene, midPt, isos, domLE, wildLE = newline[0], float(newline[1]), newline[2:4], float(newline[4]), float(newline[5])
        transcriptDict[gene] = [midPt, isos, domLE, wildLE]

# look at those isoforms
isoDict = {}
with open(sys.argv[2], "r") as infile:
    infile.readline()
    for line in infile:
        newline = line.strip().split()
        splitID = newline[0].split("_")
        gene, iso = "_".join(splitID[0:2]), splitID[2]
        if gene in transcriptDict:
            if iso in transcriptDict[gene][1]:
                isoDict[ "_".join([gene, iso]) ] = float(newline[5])
                
# classify as dom, wild, intermediate, or NA
for gene in transcriptDict:
    isos = transcriptDict[gene][1]
    domVersion = isos[0]
    wildVersion = isos[1]
    domTPM = isoDict["_".join([gene, domVersion])]
    wildTPM = isoDict["_".join([gene, wildVersion])]
    tot = domTPM + wildTPM
    version = "NA"
    if tot >= minTPM and tot > 0: # minimum TPM filter
        mid = transcriptDict[gene][0]
        propDom = domTPM / tot
        
        # check which side of the parental midpoint they fall on
        if outputType == "1":
            mid = transcriptDict[gene][0]
            if propDom > mid:
                version = "domesticated"
            elif propDom < mid:
                version = "wild"
            else:
                version = "tie"
                
            
        # check if they have proportions as extreme as the parents
        if outputType == "2":
            domLimit = transcriptDict[gene][2]
            wildLimit = transcriptDict[gene][3]
            if propDom >= domLimit:
                version = "domesticated"
            elif propDom <= wildLimit:
                version = "wild"
            else:
                version = "intermediate"

    # print
    print '\t'.join([gene, version]) 




    
