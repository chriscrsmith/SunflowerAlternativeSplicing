

# this script will be for taking in a tpm table with two* isoforms per gene, and outputting wild and dom proportions

import sys, numpy

geneDict = {}
with open(sys.argv[1], 'r') as tpmFile:
    for line in tpmFile:
        newline = line.strip().split()
        gene = newline[0]
        if gene not in geneDict:
            geneDict[gene] = [newline]
        else:
            geneDict[gene].append(newline)

# quick check that you gave it a good table
for gene in geneDict:
    if len(geneDict[gene]) != 2:
        print "\n\n  bad table: need exactly two isoforms per gene\n\n"
        print gene
        1 / 0

# now go through and get proportions
print '\t'.join([ 'domProp', 'wildProp', 'isoVersion', 'gene'   ])
for gene in geneDict:
    # get wild and domestic mean props of each isoform
    dversion = ['w','w'] # wild default
    tots = [0,0,0,0,0,0]
    means = []
    for isoNum in range(2):
        iso = geneDict[gene][isoNum]
        for p in range(6):
            tpm = float(iso[p+2])
            tots[p] += tpm
    for isoNum in range(2):
        iso = geneDict[gene][isoNum]
        props = [0,0,0,0,0,0]
        for p in range(6):
            prop = float(iso[p+2]) / tots[p]
            props[p] = prop
        wtot = 0
        for w in [0,2,3]:
            wtot += props[w]
        wmean = wtot / 3
        dtot = 0
        for d in [1,4,5]:
            dtot += props[d]
        dmean = dtot / 3
        if dmean > wmean:
            dversion[isoNum] = 'd'
        means.append([dmean, wmean])
    # what's complicated here, is if there is not* a clear domestic isoform
    # e.g., if both wild and domestic have more* isoform 1, but domestic has a disproportionate increase
    # so. first check if this happens or not
    if dversion == ['d', 'd'] or dversion == ['w', 'w']:
        print "\n\n     there's a gene without a clear domestic isoform\n\n"
        1/0
    else:
        for isoNum in range(2):
            print '\t'.join(map(str, [ means[isoNum][0], means[isoNum][1], dversion[isoNum], gene ]))



