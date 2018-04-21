
# this will take in a list of fpkm files from wild and dom lines and output proportions of isoform 1 and 2 for ILR transformation in R
# Input:
# 1. list of all isoforms that we're interested in- need the exact gene and isoform versions
# 2. comma separated list of wild fpkm files
# 3. comma separated list of dom fpkm files

import sys

addedToZero = 0.01
maxNA = 1
minTPM = 0.25


# go through and get isos you want
isoDict = {}
wildIso1Dict = {}
wildIso2Dict = {}
domIso1Dict = {}
domIso2Dict = {}
with open(sys.argv[1], 'r') as infile:
    for line in infile:
        newline = line.strip().split()
        gene, iso = newline[0:2]
        if gene not in isoDict:
            isoDict[gene] = []
            wildIso1Dict[gene] = []
            wildIso2Dict[gene] = []
            domIso1Dict[gene] = []
            domIso2Dict[gene] = []
        isoDict[gene].append(iso)

# go through and get tpm's for each iso you want, in the wild files
for wildFile in sys.argv[2].split(','):
    with open(wildFile, 'r') as infile:
        infile.readline() # header
        for line in infile:
            newline = line.strip().split()
            gene, iso, tpm = newline[1], newline[0].split('_')[2], float(newline[5])+addedToZero
            if gene in isoDict:
                if iso in isoDict[gene]:
                    isoNum = isoDict[gene].index(iso)
                    if isoNum == 0:
                        wildIso1Dict[gene].append(tpm)
                    elif isoNum == 1:
                        wildIso2Dict[gene].append(tpm)
                    else:
                        sys.stdout.write("\n\nMore than 2 isoforms for a gene\n\n")
                        1/0

# go through dom files next
for domFile in sys.argv[3].split(','):
    with open(domFile, 'r') as infile:
        infile.readline() # header   
        for line in infile:
            newline = line.strip().split()
            gene, iso, tpm = newline[1], newline[0].split('_')[2], float(newline[5])+addedToZero
            if gene in isoDict:
                if iso in isoDict[gene]:
                    isoNum = isoDict[gene].index(iso)
                    if isoNum == 0:
                        domIso1Dict[gene].append(tpm)
                    elif isoNum == 1:
                        domIso2Dict[gene].append(tpm)
                    else:
                        sys.stdout.write("\n\nMore than 2 isoforms for a gene\n\n")
                        1/0


# output header
header = ["gene", "isoVersion"]
for wildInd in range(len(sys.argv[2].split(','))):
    header.append("wild"+str(wildInd+1))
for domInd in range(len(sys.argv[3].split(','))):
    header.append("dom"+str(domInd+1))
print '\t'.join(header)


# go back through the isos you want and output proportions 
for gene in isoDict:
    props = [gene, isoDict[gene][0]] # outputting the isoform version of the first isoform listed in the table
    for wildInd in range(len(wildIso1Dict[gene])):
        total = wildIso1Dict[gene][wildInd] # for the first iso, just appending to totals, next adding 
        total += wildIso2Dict[gene][wildInd]
        if total <= (addedToZero*2):
            props.append("NA")
        elif total < minTPM:
            props.append("NA")
        else:
            propIso1 = wildIso1Dict[gene][wildInd] / total
            props.append(propIso1)
    for domInd in range(len(domIso1Dict[gene])):
        total = domIso1Dict[gene][domInd]
        total += domIso2Dict[gene][domInd]
        if total <= (addedToZero*2):
            props.append("NA")
        elif total < minTPM:
            props.append("NA")
        else:
            propIso1 = domIso1Dict[gene][domInd] / total
            props.append(propIso1)
    
    if props[2:7].count("NA") <= maxNA and props[7:12].count("NA") <= maxNA: # line checking for underrepresentation within populations
        print '\t'.join(map(str, props))


















