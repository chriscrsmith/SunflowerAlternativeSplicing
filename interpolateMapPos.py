



# this script will take in a list of snps with chrom and bp and the lg map file, and output chrom and genetic map pos
# NOTE: the lenght of the output will differ from the length of the input if the input contains lines with LG "Ha073Ns"


def main():

    import sys
    import decimal 
    decimal.getcontext().prec = 22 # setting decimal precission to 22. R's limit is 22 digits
    snpPath = sys.argv[1] # snp bp positions on genome, chrom tab bp
    mapPath = sys.argv[2] # rieseberg genetic map
    
    # let's do this by linkage group, to make list-searching quick
    lgs = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17']
    for lg in range(17):
        currentChrom = lgs[lg]

        # first, make a list of snps (list, so it's ordered by position)
        with open(mapPath, 'r') as mapFile:
            mapList = []
            for line in mapFile:
                newline = line.strip().split()
                cM = newline[2]
                chrom = newline[0]
                if chrom == currentChrom:
                    if (cM != 'cM') and (cM != 'NA'):
                        bp = decimal.Decimal(newline[1])
                        info = [ bp, decimal.Decimal(cM), chrom ]
                        mapList.append(info)
        mapList = sorted(mapList)

        # then go through your snp list and find positions for each one
        with open(snpPath, 'r') as snpFile:
            for line in snpFile:
                Hachrom = line.strip().split()[0]
                chrom = Hachrom.split('Ha')[1]
                contig = line.strip().split()[2]
                contig_pos = line.strip().split()[3]
                pop1 = line.strip().split()[4]
                pop2 = line.strip().split()[5]
                geneID = line.strip().split()[6] # this one sometimes necessary 
                if len(chrom) < 2:
                    chrom = '0'+chrom
                if chrom == currentChrom:
                    bp = decimal.Decimal(line.strip().split()[1])
                    
                    # search for the position in the genetic map
                    listIndex = 0
                    beforePosition = False
                    afterPosition = False
                    exactPosition = False
                    mapSize = len(mapList)
                    while (afterPosition == False) and (listIndex < mapSize):
                        currentPos = mapList[listIndex]
                        if bp > currentPos[0]:  # snp is greater than the current map index
                            beforePosition = currentPos # the "beforePosition" gets updated as long as it's before the SNP
                        if bp < currentPos[0]:
                            afterPosition = currentPos
                        if bp == currentPos[0]:
                            exactPosition = currentPos
                        listIndex += 1 # increment map index

                    # output
                    if exactPosition:
                        print chrom, bp, exactPosition[1], contig, contig_pos, pop1, pop2, geneID # comment out geneID sometimes
                    elif (beforePosition) and (afterPosition):
                        # if no exact match, interpolate
                        beforeBP = beforePosition[0]
                        afterBP = afterPosition[0]
                        beforecM = beforePosition[1]
                        aftercM = afterPosition[1]
                        refdiff = afterBP - beforeBP
                        snpdiff = bp - beforeBP
                        ratio = snpdiff / refdiff
                        cMdiff = aftercM - beforecM
                        finalcM = beforecM + (cMdiff*ratio)
                        print chrom, bp, finalcM, contig, contig_pos, pop1, pop2, geneID # comment out geneID sometimes
                    else:
                        print chrom, bp, 'NA', contig, contig_pos, pop1, pop2, geneID # comment out geneID sometimes


    snpFile.close()

main()
