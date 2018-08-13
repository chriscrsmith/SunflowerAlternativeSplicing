


#    output types
# I'm redoig this to make it correct and more simple
# 1. First and foremost, I want a list of independent qtl
# 2. Second, I want a list of genes and which (numbered/named) qtl are associated with them for go analysis
# 3. Third, I want a list of genes and which qtl are associated with them and their locations and info
# 4. Fourth, for circos links, the same list as above but only trans qtl

# this script will take in
# 1. list of LOD hits inside QTL regions identidied by rQTL.R
# 2. splice type file, two cols, first call gene, second call splicing model
# 3. output type (see below for options)

import sys
from operator import itemgetter


outputType = sys.argv[3]

# read in all LOD data
lodDict = {}
chromDict = {} # for constructing independent regions (below)
with open(sys.argv[1], 'r') as lodFile:
    for line in lodFile:
        newline = line.strip().split()
        gene, chrom, pos, lod, start, end, varExplained, cisTran = newline
        if gene not in lodDict:
            lodDict[gene] = {}
            lodDict[gene][chrom] = [[float(pos), float(lod), float(start), float(end), float(varExplained), cisTran, chrom]]
        else:
            if chrom not in lodDict[gene]:
                lodDict[gene][chrom] = [[float(pos), float(lod), float(start), float(end), float(varExplained), cisTran, chrom]]
            else:
                lodDict[gene][chrom].append( [float(pos), float(lod), float(start), float(end), float(varExplained), cisTran, chrom] )
        if chrom not in chromDict:
            chromDict[chrom] = [[float(start), float(end)]]
        else:
            chromDict[chrom].append( [float(start), float(end)] )




# read in splice type
typeDict = {}
with open(sys.argv[2]) as spliceTypeFile:
    for line in spliceTypeFile:
        newline = line.strip().split()
        typeDict[newline[0]] = newline[1]







# trying to identify independent overlapping regions 
# make a list of clusters, one cluster per region. Next, will try to add clusters together recursively              
indepList = {}
for chrom in chromDict:
    chromList = list(chromDict[chrom])
    removeList = []
    combineList = []
    repeat = True
    while repeat == True:
        repeat = False
        for hit1 in range(0, len(chromList)-1):
            for hit2 in range(hit1+1, len(chromList)):
                if (round(chromList[hit1][0],5) <= round(chromList[hit2][1],5)) and (round(chromList[hit1][1],5) >= round(chromList[hit2][0],5)): # I checked and they always have the same orientation
                    combined = [min(chromList[hit1][0],chromList[hit2][0]), max(chromList[hit1][1],chromList[hit2][1])]
                    if combined not in combineList:
                        combineList.append( combined ) # combining into 1
                    if chromList[hit1] not in removeList:
                        removeList.append(chromList[hit1])
                    if chromList[hit2] not in removeList:
                        removeList.append(chromList[hit2])
                    repeat = True
        for r in removeList:
            count = chromList.count(r)
            for i in range(count):
                chromList.remove(r)
        removeList = []
        for c in combineList:
            chromList.append(c)
        combineList = []
    indepList[chrom] = chromList






# print the independent regions and the number of genes regulated by each region
qtlID = 0
whichQtls = {} # keeping track of which qtl are associated with which gene
for chrom in range(1,17+1):
    for ind in indepList[str(chrom)]:
        geneCount = 0
        qtlID += 1
        for gene in lodDict:
            if str(chrom) in lodDict[gene]:
                for hit in lodDict[gene][str(chrom)]:
                    if (round(float(hit[0]),5) >= round(float(ind[0]),5)) and (round(float(hit[0]),5) <= round(float(ind[1]),5)): 
                    # if the hit falls within the region 
                        geneCount += 1
                        if gene not in whichQtls:
                            whichQtls[gene] = [str(qtlID)]
                        else:
                            whichQtls[gene].append(str(qtlID))
                        if outputType == '3':
                            newOut = [ gene, qtlID, chrom, hit[2], hit[3], hit[4] ]
                            print '\t'.join(map(str,newOut))
                        if outputType == '4':
                            ct = hit[5]
                            if ct == 'trans':
                                newOut = [ gene, qtlID, chrom, hit[2], hit[3], hit[4] ]
                                print '\t'.join(map(str,newOut))
        if outputType == '1':
            print '\t'.join( map(str,[chrom, ind[0], ind[1], geneCount]) )
if outputType == '2':
    for gene in whichQtls:
        outline = gene + '\t' + ','.join( whichQtls[gene] )
        print outline


    
        

# # how often cis, trans, or both?    
# typeKey = {'1':'exonSkip', '3':'altSpliceSite', '5':'exonSkip', '9':'substitution', '10':'intronRetention'}
# if outputType == '4':
#     print '\t'.join([ 'geneName', 'spliceType', 'cisVtrans', 'varCis', 'varTrans', 'totVarExpl', 'proportionCis' ])
#     for ge in lodDict:
#         # grab splice type
#         st = typeKey[ typeDict[ge] ]
#         tr = False
#         ci = False
#         trVar = 0
#         ciVar = 0
#         for ch in lodDict[ge]:
#             for hi in lodDict[ge][ch]:
#                 if hi[5] == 'trans':
#                     tr = True
#                     trVar += float(hi[4])
#                 if hi[5] == 'cis':
#                     ci = True
#                     ciVar += float(hi[4])
#         tot = ciVar + trVar
#         if tr == True and ci == True:
#             propCis = ciVar/(ciVar+trVar)
#             print '\t'.join( map(str, [ge, st, 'both', ciVar, trVar, tot, propCis]) )
#         elif tr == True:
#             propCis = 0
#             print '\t'.join( map(str, [ge, st, 'trans', ciVar, trVar, tot, propCis]) )
#         elif ci == True:
#             propCis = 1
#             print '\t'.join( map(str, [ge, st, 'cis', ciVar, trVar, tot, propCis]) )
#         else:
#             propCis = 'NA'
#             print '\t'.join( map(str, [ge, st, 'NA', ciVar, trVar, tot, propCis]) )










