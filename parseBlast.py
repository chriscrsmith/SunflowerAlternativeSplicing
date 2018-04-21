

# NOTE: format requirement: the blast file MUST end with a line that reads "# BLAST processed ..."
# this will parse (1) a blast output file, and (2) the fixed snps file, and output the best hits with chrom and position, 
#     but only if there's a good hit, 
#     and only if there aren't 2 equally best hits,
#     AND only if it's not on one of the 17 chromosomes and not cp, mt, or that other strange one

def main():

    import sys
    infile = open(sys.argv[1], 'r')
    snpfile = open(sys.argv[2], 'r')
    lgs = ['Ha1', 'Ha10', 'Ha11', 'Ha12', 'Ha13', 'Ha14', 'Ha15','Ha16','Ha17','Ha2','Ha3','Ha4','Ha5','Ha6','Ha7','Ha8', 'Ha9']

    ##### dictionary of snps #####
    snpDict = {}
    for line in snpfile:
        newline = line.strip().split()
        snp = newline[0]
        if snp != 'contig':
            pos = int(newline[1])
            pop1 = newline[2]
            pop2 = newline[3]
            info = [pos,pop1,pop2]
            if snp not in snpDict:
                snpDict[snp] = [info]
            else:
                snpDict[snp].append(info)

    ##### parsing blast output #####
    firstQuery = True
    inside_header = False
    bestHit = False
    for line in infile:
        newline = line.strip().split()
        if newline[0] == '#': # header lines
            if (newline[1] == "BLAST") and (newline[2] == "processed"): # if you reach the end of the file do this bit
                contig = bestHit[0]
                snps = snpDict[contig]
                chromo = bestHit[1]
                qstart = int(bestHit[6])
                qend = int(bestHit[7])
                rstart = int(bestHit[8])
                rend = int(bestHit[9])
                for snp in snps: # for each snp, find position in new reference                                                              
                    bp = snp[0]
                    if rstart < rend:
                        ref_pos = rstart + (bp - qstart)
                    else:
                        ref_pos = rstart - (bp - qstart)
                    print '\t'.join([ chromo, str(ref_pos), contig, str(bp), snp[1], snp[2] ])

            if inside_header == False: # this is the first commented line of each new query
                inside_header = True
                if firstQuery == False: # output results from previous query, now, unless you're looking at the first query
                    if bestHit:
                        if skipSNP == False: # if there was a good hit for the transcriptome contig, then output snps
                            contig = bestHit[0]
                            snps = snpDict[contig]
                            chromo = bestHit[1]
                            qstart = int(bestHit[6])
                            qend = int(bestHit[7])
                            rstart = int(bestHit[8])
                            rend = int(bestHit[9])
                            for snp in snps: # for each snp, find position in new reference
                                bp = snp[0]
                                if rstart < rend:
                                    ref_pos = rstart + (bp - qstart)
                                else: # here, the ref is reverse complimented in the blast output
                                    ref_pos = rstart - (bp - qstart)
                                print '\t'.join([ chromo, str(ref_pos), contig, str(bp), snp[1], snp[2] ])
                firstQuery = False
                bestHit = False # restart best hit

        else: # data lines
            chrom = newline[1]
            score = float(newline[2])
            length = int(newline[3])
            ref_start = int(newline[8])
            if inside_header == True: # first line of data
                inside_header = False
                bestScore = False
                skipSNP = False
                if chrom in lgs: # one of the 17 chromosomes
                    if score >= 95: # at least 90% ID
                        if length >= 100: # at least 100bp match
                            bestHit = newline # then, current best score is this
                            bestScore = [ score,length ]
            else: # subsequent lines of data
                if bestHit:
                    if chrom != bestHit[1]: # if different chrom than the current best hit
                        if score >= 95:
                            if length >= 100:
                                if chrom in lgs:
                                    skipSNP = True # there are multiple, equal best hits on different chromosomes, so not keeping this snp
                    else: # if same chrom as current best hit
                        if abs(ref_start - int(bestHit[8])) > 1000000: # if the hits are greater than 1kb apart, skip this contig
                            skipSNP = True
                        
                else:
                    if chrom in lgs:
                        if score >= 95:
                            if length >= 100: 
                                bestHit = newline # then, current best score is this
                                bestScore = [ score,length ]

    infile.close()
    snpfile.close()

main()



