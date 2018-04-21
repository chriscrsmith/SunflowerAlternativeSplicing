
import sys
from Bio.Seq import Seq


hitDict = {}
def process(newline):
    gene = "_".join(newline[0].split("_")[0:2])
    iso = newline[0].split("_")[2]
    subjectID = newline[1]
    if gene not in hitDict:
        hitDict[gene] = {}
    if iso not in hitDict[gene]:
        hitDict[gene][iso] = {}
    if subjectID not in hitDict[gene][iso]:
        hitDict[gene][iso][subjectID] = []
    hitDict[gene][iso][subjectID].append(newline)


    

# read in data
with open(sys.argv[1], "r") as infile:
    infile.readline() # header
    for line in infile:
        newline = line.strip().split()
        process(newline)

with open(sys.argv[2], "r") as infile: # reading in both blast files
    infile.readline()
    for line in infile:
        newline = line.strip().split()
        process(newline)

isoDict = {}
with open(sys.argv[3], "r") as infile:
    for line in infile:
        newline = line.strip().split()
        if newline[0][0] == ">":
            isoName = newline[0].split(">")[1]
        else:
            isoDict[isoName] = newline[0]

sangerDict = {}
with open(sys.argv[4], "r") as infile:
    for line in infile:
        newline = line.strip().split()
        if newline[0][0] == ">":
            subName = newline[0].split(">")[1]
        else:
            sangerDict[subName] = newline[0]
with open(sys.argv[5], "r") as infile:
    for line in infile:
        newline = line.strip().split()
        if newline[0][0] == ">":
            subName = newline[0].split(">")[1]
        else:
            sangerDict[subName] = newline[0]






# now analyze
for gene in hitDict:
    for iso in hitDict[gene]:
        for sub in hitDict[gene][iso]:
            for hit in hitDict[gene][iso][sub]:
                queryLength = float(hit[3])
                subjectLength = float(hit[4])
                alignmentLength = float(hit[5])
            if (alignmentLength / queryLength) > 0.75:
                # this way, we're simply looking for clean, one-chunk blast hits that cover >90% of the transcript
                # can filter the output for genes, with "good" hits for each isoform, on different sanger reads
                sangerRead = sangerDict[sub]
                transcript = isoDict["_".join([gene, iso])]
                print
                print ">goodHit", gene, iso, sub
                print sangerRead
                print transcript









        
