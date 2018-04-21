


# this will make bands from a qtl file, and a gene locs file
# NOTE: for the qtl file, need to give only trans qtl, or else won't be able to output cis regulated genes


# e.g. python ../Scripts/makeCircosBands.py ../InterpolateMapPos/lg.ALL.bronze14.path.txt ../InterpolateMapPos/lgLengths.txt ../Clean_gene_lists/Significant/varExplained_1000perm_08162017_greaterThan1varExp_genesAndQtlInfo_transOnly_sameChromOnly_plus1.txt ../Clean_gene_lists/PlacingOntoGenome/sigQvals_pt05_1fpkm_ILR_04062016_cMpositions_plus1.txt ../Clean_gene_lists/Significant/varExplained_1000perm_08162017_greaterThan1varExp_propCis.txt 

import sys
sys.path.insert(0, '/data5/AltSplicing09182015/Scripts/')
from interpolateBpPos import fillGeneticMap
from interpolateBpPos import interpolateBP

# fill genetic map
geneticMap = fillGeneticMap(sys.argv[1])

qtlDict = {}
for i in range(1,18):
    qtlDict[str(i)] = []
with open(sys.argv[3], 'r') as qtlFile:
    for line in qtlFile:
        newline = line.strip().split()
        gene, chrom, start, end, varExplained = newline[0], newline[2], newline[3], newline[4], newline[5]
        start = interpolateBP(chrom, float(start), geneticMap )
        end = interpolateBP(chrom, float(end), geneticMap )
        if [start, end] not in qtlDict[chrom]:
            qtlDict[chrom].append([start, end, gene])
            
cisDict = {}
with open(sys.argv[5], 'r') as cisFile:
    for line in cisFile:
        newline = line.strip().split()
        gene, prop = newline[0], newline[1]
        if prop == '1.0':
            cisDict[gene] = 0

geneDict = {}
for i in range(1,18):
    geneDict[str(i)] = []
with open(sys.argv[4], 'r') as geneFile:
    for line in geneFile:
        newline = line.strip().split()
        gene, chrom, start, end = newline[7], newline[0], newline[1], newline[1]
        if chrom[0] == '0':
            chrom = chrom[1:]
        if chrom != "0_73Ns":
            if gene in cisDict:
                geneDict[chrom].append([start, end, gene])


chromLengths = {}
with open(sys.argv[2], 'r') as chromFile:
    for line in chromFile:
        newline = line.strip().split()
        chrom, length = newline[0], newline[1]
        if chrom[0] == '0':
            chrom = chrom[1:]
        chromLengths[chrom] = float(length)

# make qtl bands
radius = 5000000
color = "black"
for chrom in qtlDict:
    bandCounter = 0
    for qtl in qtlDict[chrom]:
        bandCounter += 1
        start, end, gene = qtl[:]
        newMid = float(start) + ((float(end) - float(start)) / 2)
        newStart = int(newMid - radius)
        newEnd = int(newMid + radius)
        if newStart < 0:
            newStart = 0
            newEnd = (radius*2)
        if newEnd > chromLengths[chrom]:
            newStart = chromLengths[chrom] - (radius*2)
            newEnd = chromLengths[chrom]            

        print '\t'.join(["band", 
                         "lg"+chrom, 
                         "q"+chrom+"."+str(bandCounter), 
                         "q"+chrom+"."+str(bandCounter), 
                         str(newStart), 
                         str(newEnd),
                         color])
        

# make cis-regulated gene bands
for chrom in geneDict:
    bandCounter = 0
    for gene in geneDict[chrom]:
        color = "white"
        bandCounter += 1
        start, end, g = gene[:]
        midPoint = float(start) + ((float(end) - float(start)) / 2)
        newMid = float(midPoint)
        newStart = int(newMid - radius)
        newEnd = int(newMid + radius)

        # color specific bands                                                                                                   
        if g == "c107290_g1":
            color = "yellow"
            newStart = int(midPoint - radius*2)
            newEnd = int(midPoint + radius*2)
        if g == "c106588_g1":
            color = "purple"
            newStart = int(midPoint - radius*2)
            newEnd = int(midPoint + radius*2)
        pList = ["c94209_g2", "c107987_g3","c90415_g2","c43713_g1","c95147_g1","c93396_g1","c95232_g2","c107251_g1","c97577_g2","c97044_g1","c99599_g1","c107740_g1","c88028_g1"]
        if g in pList and chrom == "5":
            color = "green"
            newStart = int(midPoint - radius*2)
            newEnd = int(midPoint + radius*2)
        #
        
        if newStart < 0:
            newStart = 0
            newEnd = (radius*2)
        if newEnd > chromLengths[chrom]:
            newStart = chromLengths[chrom] - (radius*2)
            newEnd = chromLengths[chrom]
        print '\t'.join(["band",
                         "lg"+chrom,
                         "g"+chrom+"."+str(bandCounter),
                         "g"+chrom+"."+str(bandCounter),
                         str(newStart),
                         str(newEnd),
                         color])
