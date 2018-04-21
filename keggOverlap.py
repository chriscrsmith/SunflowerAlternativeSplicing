

import sys

spliceoDict = {}
with open(sys.argv[1], "r") as spliceosomeLocs:
    for line in spliceosomeLocs:
        chrom, pos, cM, contig, length, nuc1, nuc2, gene = line.strip().split()
        if chrom not in spliceoDict:
            spliceoDict[chrom] = [[chrom, pos, cM, contig, length, nuc1, nuc2, gene]]
        else:
            spliceoDict[chrom].append([chrom, pos, cM, contig, length, nuc1, nuc2, gene])

            
with open(sys.argv[2], "r") as qtls:
    for line in qtls:
        chrom, start, end, numGenes = line.strip().split()
        if len(chrom) ==1:
            chrom = "0" + chrom
        for spliceoGene in spliceoDict[chrom]:
            pos = spliceoGene[2]
            if (float(pos) >= float(start)) and (float(pos) <= float(end)):
                print chrom, start, end, numGenes, '\t'.join(spliceoGene)

