


# e.g. python /data5/AltSplicing09182015/Scripts/makeCircosLinks.py /data5/AltSplicing09182015/InterpolateMapPos/lg.ALL.bronze14.path.txt /data5/AltSplicing09182015/InterpolateMapPos/lgLengths.txt /data5/AltSplicing09182015/Clean_gene_lists/PlacingOntoGenome/sigQvals_pt05_1fpkm_ILR_04062016_cMpositions.txt /data5/AltSplicing09182015/Clean_gene_lists/Significant/varExplained_1000perm_08162017_greaterThan1varExp_genesAndQtlInfo_transOnly.txt 

import sys, numpy
sys.path.insert(0, '/data5/AltSplicing09182015/Scripts/')
from interpolateBpPos import fillGeneticMap
from interpolateBpPos import interpolateBP

# fill genetic map
geneticMap = fillGeneticMap(sys.argv[1])

# get lg lengths
lengthDict = {}
with open(sys.argv[2], 'r') as lgLengthFile:
    for line in lgLengthFile:
        newline = line.strip().split()
        chrom, length = newline[0],newline[1]
        if chrom[0] == '0':
            chrom = chrom[1:]
        lengthDict[chrom] = float(length)

                                                            
# read in the gene locations
locDict = {}
with open(sys.argv[3], 'r') as geneLocsFile:
    for line in geneLocsFile:
        newline = line.strip().split()
        geneName = newline[7]
        chrom = newline[0]
        if chrom[0] == '0':
            chrom = chrom[1:]
        if chrom != '0_73Ns':
            pos = float(newline[1])
            locDict[geneName] = [chrom, pos]

#def blowUpLinkWidth():
    

            
# finally, read in the qtl data, and output links
with open(sys.argv[4], 'r') as qtlFile:
    c = 2500000
    coef = (float(1)/float(2))
    for line in qtlFile:
        newline = line.strip().split()
        gene = newline[0]
        if gene in locDict:
            chromLenQTL = lengthDict[newline[2]]
            chromLenGene = lengthDict[locDict[gene][0]]
            varExplained = float(newline[5])
            qtlMid = ((float(newline[4]) - float(newline[3]))/2) + float(newline[3])
            pos1 = interpolateBP(newline[2], qtlMid, geneticMap ) - (c * ((varExplained+1)**coef))
            pos2 = interpolateBP(newline[2], qtlMid, geneticMap) + (c * ((varExplained+1)**coef))
            geneMid = locDict[gene][1]
            pos3 = geneMid - (c * ((varExplained+1)**coef))
            pos4 = geneMid + (c * ((varExplained+1)**coef))
            if pos1 < 0 and pos2 > chromLenQTL: # check if the link is bigger than the chromosome
                print "          huge", pos1, pos2, chromLenQTL
                1/0
            if pos3 < 0 and pos4 > chromLenGene:
                print "          huge", pos3, pos4, chromLenGene
                1/0
            if pos1 < 0: # check if the links run off the end of the chromosome
                pos1 = 0
                pos2 = 2* (c * ((varExplained+1)**coef))
            if pos3 < 0:
                pos3 = 0
                pos4 = 2* (c * ((varExplained+1)**coef))
            if pos2 > chromLenQTL:
                pos2 = chromLenQTL
                pos1 = chromLenQTL - (2* (c * ((varExplained+1)**coef)))
            if pos4 > chromLenGene:
                pos4 = chromLenGene
                pos3 = chromLenGene -(2* (c * ((varExplained+1)**coef)))
            if (pos2 >= pos3 and pos2 <= pos4): # if the start and end are overlapping on same chrom  
                pass
            elif (pos1 >= pos3 and pos1 <= pos4): 
                pass
            link = ' '.join([ "lg"+newline[2], str(int(pos1)), str(int(pos2)), "lg"+locDict[gene][0], str(int(pos3)), str(int(pos4))    ])
            if gene == "c107196_g1" or gene == "c94264_g1":
                
                link = ' '.join([ "lg"+newline[2], str(int(pos1)), str(int(pos2)), "lg"+locDict[gene][0], str(int(pos3)), str(int(pos4)), "color=yellow_a1" ])
                print link
            elif gene == "c110619_g2" or (gene == "c110699_g2" and newline[2] == "4"):
                link = ' '.join([ "lg"+newline[2], str(int(pos1)), str(int(pos2)), "lg"+locDict[gene][0], str(int(pos3)), str(int(pos4)), "color=blue_a1" ])
                print link
            elif gene in ["c94209_g2","c107987_g3","c90415_g2","c43713_g1","c95147_g1","c93396_g1","c95232_g2","c107251_g1","c97577_g2","c97044_g1","c99599_g1","c107740_g1","c88028_g1"] and newline[2] == "5":
                link = ' '.join([ "lg"+newline[2], str(int(pos1)), str(int(pos2)), "lg"+locDict[gene][0], str(int(pos3)), str(int(pos4)), "color=green_a1" ])
                print link
            elif newline[2] != locDict[gene][0]: # skipping same-chromosome links
                print link





