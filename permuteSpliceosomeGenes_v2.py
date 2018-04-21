

# this script differs from the first version, because it will permute positions of QTL regions, instead of spliceosome genes
# find what proportion contain spliceosome genes, and output that proportion.

# Input:
# 1. qtl locations, e.g. varExplained_1000perm_08162017_greaterThan1varExp_indepQtl.txt
# 2. spliceosome gene locations (only for getting the number to permute), e.g. tomatoSpliceosomeGene_homologues.txt
# 3. lg lengths, because the lgs have different lengths, and the probability of a gene falling on each one is therefore different

import sys, random

# get lg lengths
cumulativeLgLength = 0
cumulativeLgLengths = []
with open(sys.argv[3], "r") as lengthsFile:
     for line in lengthsFile:
          theLength = float(line.strip().split()[1])
          cumulativeLgLength += theLength
          cumulativeLgLengths.append(cumulativeLgLength)



          
# get the sizes of QTL to permute
QTLsizes = []
with open(sys.argv[1], "r") as qtls:
     for line in qtls:
          start, end = line.strip().split()[1:3]
          size = float(end)-float(start)
          QTLsizes.append(size)
QTLsizes.append(size)


          



# permute QTL locations
QTLlocations = [None] * len(QTLsizes)
for q in range(len(QTLsizes)):
     # the lg is non random because they are different lengths. Using jerry-rigged (clever!) way to find the lg
     randLgNum = random.uniform(0,1) * cumulativeLgLength
     tempList = list(cumulativeLgLengths)
     tempList.append(randLgNum)
     tempList.sort()
     lg = tempList.index(randLgNum) + 1
     # and a uniform cM position
     cMstart = random.uniform(0,100)
     cMend = cMstart + QTLsizes[q]
     while cMend > 100:
          cMstart = random.uniform(0,100)
          cMend = cMstart + QTLsizes[q]
     QTLlocations[q] = [lg, cMstart, cMend]




# read in spliceosome gene locations
spliceosomeLocs = {}
for c in range(17):
     spliceosomeLocs[str(c+1)] = []
with open(sys.argv[2], "r") as spliceosomeGenes:
     for line in spliceosomeGenes:
          newline = line.strip().split()
          chrom = newline[0]
          if chrom[0] == "0":
               chrom = chrom[1:]
          pos = newline[2]
          spliceosomeLocs[chrom].append(pos)


# find the proportion that overlap
overlapping = 0
for q in QTLlocations:
     overlap = False
     lg, start, end = q[0:3]
     for s in spliceosomeLocs[str(lg)]:
          if (float(s) >= start) and (float(pos) <= end):
               overlap = True
     if overlap == True:
          overlapping += 1


# print proportion that overlapped
print float(overlapping) / float(len(QTLlocations))
