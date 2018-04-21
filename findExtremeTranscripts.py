
# this one will take a tpm table for the significant genes, and output ones that are extreme in the parents

import sys, numpy

print '\t'.join(['geneID','domIso','wildIso'])

with open(sys.argv[1], "r") as infile:
    firstIso = True
    for line in infile:
        newline = line.strip().split()
        gene = newline[0]
        if firstIso == True:
            iso1 = newline[1]
            dom1 = map(float, [ newline[3], newline[6], newline[7] ])
            wild1 = map(float, [ newline[2], newline[4], newline[5] ])
            firstIso = False
        else:
            iso2 = newline[1]
            dom2 = map(float, [ newline[3], newline[6], newline[7] ])
            wild2 = map(float, [ newline[2], newline[4], newline[5] ])
            firstIso = True

            # calculate proportions
            domTots = numpy.sum([dom1, dom2], axis = 0)
            domProp1s = numpy.divide(dom1, domTots)
            domProp1 = numpy.mean( domProp1s )
            wildTots = numpy.sum([wild1, wild2], axis = 0)
            wildProp1s = numpy.divide(wild1, wildTots)
            wildProp1 = numpy.mean( wildProp1s )

            # outputting midpoint, parent ranges, iso versions
            if domProp1 > wildProp1:
                leastExtremeDom = min(domProp1s)
                leastExtremeWild = max(wildProp1s)
                mid = ((leastExtremeDom-leastExtremeWild) /2) + leastExtremeWild
                print '\t'.join(map(str, ([ gene, mid, iso1, iso2, leastExtremeDom, leastExtremeWild])))
            else:
                domProp2s = numpy.subtract([1,1,1], domProp1s)
                wildProp2s = numpy.subtract([1,1,1], wildProp1s)
                leastExtremeDom = min(domProp2s)
                leastExtremeWild = max(wildProp2s)
                mid = ((leastExtremeDom-leastExtremeWild) /2) +leastExtremeWild
                print '\t'.join(map(str, ([ gene, mid, iso2, iso1, leastExtremeDom, leastExtremeWild])))
                


                
