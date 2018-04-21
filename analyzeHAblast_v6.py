


# this script will take in (1) the blast output from the command:
# blastn -query retention_deletion_1TPM_01052016.fa -db ../../InterpolateMapPos/Ha412v1r1_genome.fasta -outfmt "6 qseqid sseqid pident qlen length qstart qend sstart send evalue gaps" -perc_identity 90 > blastOut3.txt
# (2) the filtered TPM file to compare means and get wild v dom versions
# and (3) the output type: 1 for gene list, 2 for visualization, 3 list with splice type per blast hit
# we are looking for "good" looking genes. 
# don't forget that for SOME of these genes, we will be able to use the HA genome to confidently identify splicing
# for the other genes (will be interesting to see if this is a disproportionate number or not), we will not be able to rule out allelic variation 

# We are interested in several different possible scenarios in this blast analysis:
# some will clearly indicate alternative splicing, some will rule out splicing, and others will remain ambiguous because we don't have enough information.
# I have come up with 8 models to categorize what we find, here.

# Different number of blast hits per isoform
# (1) (given HA reference) HA89 isoform has more hits and 1238 has fewer, longer hits with same start/end reference positions
# splicing: intron retention in 1238 (model 1) / exon skipping in HA89 (model 10)
# (2) HA89 isoform has fewer hits and 1238 has more hits
# ambiguous: could be genomic deletion in 1238, we don't have that genome

# Same number of blast hits, different reference positions
# (3) 1238 version longer, HA89 version shorter; both iso versions share start and end reference positions***
# splicing: alternative splice site in HA89      edit 03132017: it doesn't seem like this is right? could be insertion in 1238?
# (4) 1238 version shorter, HA89 version longer
# ambiguous: could be genomic deletion in 1238, we don't have that genome

# Same reference positions, different query lengths
# (5) HA89 version shorter and aligns to larger portion of genome, 1238 version longer
# splicing: intron retention in 1238
# (6) HA89 version longer, 1238 version short and aligns to larger portion of the genome
# ambiguous: could be genomic deletion in 1238, can't tell without that genome
# (7) HA89 version shorter, 1238 version longer and aligns to smaller portion of the genome
# NOT splicing: 1238 genomic insertion
# (8) HA89 version long and aligns to smaller portion of the genome, 1238 version shorter
# (would not expect to see this) NOT splicing: genomic insertion in HA89

# Same number of blast hits, at least one of the hits are completely non-overlapping
# (9) this should be checked BEFORE models 3 and 4, because it makes more sense    # NEED TO MAKE SURE ALL EXON/INTRON SEQUENCES ARE FOUND IN THE REFERENCE
# splicing: substitution-type 




def main():

    import sys, numpy

    outputType = sys.argv[3]
    lgs = ['Ha1', 'Ha10', 'Ha11', 'Ha12', 'Ha13', 'Ha14', 'Ha15','Ha16','Ha17','Ha2','Ha3','Ha4','Ha5','Ha6','Ha7','Ha8', 'Ha9']

    maxIntronSize = 500000
    minSpliceSiteDifference = 10
    maxQueryGap = 10 # maximum gap OR OVERLAP between two query start/ends of two blast hits
    maxRefStartEndOffset = 1 # min 1 # the two isoforms may map to the reference at difference start/end positions
                               # surprisingly, this hurts your numbers. It may be because it allows not-good exons to get added onto the growing transcript, and preventing other exons from being added on

    #########################################
    # read in and organize all the blast data
    #########################################
    genes = {} # sorted list of genes with corresponding isos
    isos = {} # sorted list of isos with corresponding blast hits
    with open(sys.argv[1], 'r') as infile:
        for line in infile:
            newline = line.strip().split()
            iso = newline[0]
            gene = '_'.join(iso.split('_')[:-1])
            if gene not in genes:
                genes[gene] = [iso]
                isos[iso] = [newline]
            else:
                if iso not in genes[gene]:
                    genes[gene].append(iso)
                    isos[iso] = [newline]
                else:
                    isos[iso].append(newline)





    ###########################################################################
    # read in the TPM file, compare wild vs domestic to get wild v dom versions
    ###########################################################################
    # first, dictionary of genes and isos in the filtered tpm table
    import sys, numpy
    tpmDict = {}
    with open(sys.argv[2], 'r') as infile:
        infile.readline() # getting rid of header
        for line in infile:
            newline = line.strip().split()
            geneNAME = newline[0]
            if geneNAME not in tpmDict:
                tpmDict[geneNAME] = [newline]
            else:
                tpmDict[geneNAME].append(newline)

    # quickly filter for only genes with two isoforms
    removeListy = []
    for potential in tpmDict:
        if len(tpmDict[potential]) != 2:
            removeListy.append(potential)
    for removed in removeListy:
        del tpmDict[removed]
        if removed in genes: # not all genes had blast hits
            del genes[removed]

    # calculate proportions to distinguish wild versus dom versions
    deleteGeneList = []
    for gene in genes:
        deleteGene = False
        totals = [0,0,0,0,0,0]
        props = [ [0,0,0,0,0,0], [0,0,0,0,0,0] ]
        if gene in tpmDict:
            for iso in range(2): # getting totals
                isoInfo = tpmDict[gene][iso]
                tpms = isoInfo[2:]
                for i in range(6):
                    totals[i] += float(tpms[i]) 
            for i in range(6):
                if totals[i] == 0:
#                    totals[i] = None
                    deleteGene = True
                    deleteGeneList.append(gene)
            if deleteGene == False:
                for iso in range(2): # getting proportions
                    isoInfo = tpmDict[gene][iso]
                    tpms = isoInfo[2:]
                    for i in range(6):
                        if totals[i] != None:
                            prop = (float(tpms[i]) + 0.00001) / float(totals[i])
                            props[iso][i] = prop
            # now, check which is the wild version
            if deleteGene == False:
                version = None
                x1238_iso1 = [ props[0][0], props[0][2], props[0][3] ]
                #x1238_iso2 = [ props[1][0], props[1][2], props[1][3] ]
                mean1238_iso1 = numpy.mean(x1238_iso1)
                #mean1238_iso2 = numpy.mean(x1238_iso2)
                ha89_iso1 = [ props[0][1], props[0][4], props[0][5] ] 
                #ha89_iso2 = [ props[1][1], props[1][4], props[1][5] ]
                meanHA89_iso1 = numpy.mean(ha89_iso1)
                #meanHA89_iso2 = numpy.mean(ha89_iso2)
                if mean1238_iso1 > meanHA89_iso1:
                    tpmDict[gene][0].append('1238')
                    tpmDict[gene][1].append('HA89')
                elif mean1238_iso1 < meanHA89_iso1:
                    tpmDict[gene][0].append('HA89')
                    tpmDict[gene][1].append('1238')
                else:
                    pass # tie

    # lastly, go through and delete genes where at least one parent had 0 expression
    changeToSet = set(deleteGeneList)
    for gene in changeToSet:
        del genes[gene]





    ############################################################################
    # now go through the genes individually and find alternatively spliced genes
    ############################################################################
    outputDict = {}
    otherDict= {}
    for gene in genes:

        ##### make sure both isoforms have hits #####
        if len(genes[gene]) == 2: # there are a handful of genes, 15, with only one isoform hitting the genome

            ##### make dictionary of hits on each chromosome #####
            hitDict = {}
            isoNum = 0
            for iso in genes[gene]:
                for hit in isos[iso]:
                    chrom = hit[1]
                    if chrom not in hitDict:
                        hitDict[chrom] = [ [],[] ] # list of two lists- one for each isoform
                        hitDict[chrom][isoNum].append(hit)
                    else:
                        hitDict[chrom][isoNum].append(hit)
                isoNum += 1 # only 0 or 1

            ##### go through the hits, and put together exons #####
            transcript = None
            for chromosome in hitDict:
                for isoform in range(2):
                    for currentHit in hitDict[chromosome][isoform]:
                        hits = list( hitDict[chromosome] ) # use list() or else python links the variables
                        if len(hits) == 2: # making sure we have hits for both isoforms, not just one. This doesn't mean that both isoforms have hits corresponding to EACH 'best hit'
                            # see which isoform the best hit came from, and start with that isoform
                            if isoform == 0:
                                isoList = list(hits[0])
                                isoList.pop(isoList.index(currentHit)) # remove hits from the list as they are added to the growing transcript
                                currentIso = 'iso1' # to be used further down 
                            else:
                                isoList = list(hits[1])
                                isoList.pop(isoList.index(currentHit))
                                currentIso = 'iso2'
                            transcript = list( [currentHit] ) # this is the transcript we will start building
                            if ( int(transcript[-1][8]) - int(transcript[0][7])) > 0:             
                                strand = '+'
                            else:
                                strand = '-'
                            ##### recursively add on exons to the growing transcript #####
                            removeList = ['space_holder'] # need to define this before the following while loop
                            while (len(removeList) > 0) and (len(isoList) > 0): # things get added to remove list as they are added to the growing transcript
                                removeList = []
                                for hit in isoList:
                                    exonRefStart = int(hit[7])
                                    exonRefEnd = int(hit[8])
                                    if ( int(exonRefEnd) - int(exonRefStart) ) > 0: # hit is positive (+) strand
                                        hitStrand = '+'
                                    else:
                                        hitStrand = '-'
                                    if hitStrand == strand:
                                        exonQueryStart = int(hit[5])
                                        exonQueryEnd = int(hit[6])
                                        transcriptRefStart = int(transcript[0][7])
                                        transcriptRefEnd = int(transcript[-1][8])
                                        transcriptQueryStart = int(transcript[0][5])
                                        transcriptQueryEnd = int(transcript[-1][6])
                                        if abs(exonRefStart - transcriptRefEnd) < maxIntronSize: # making sure hits do not have too-large of gaps between them 
                                            if abs(exonQueryStart - transcriptQueryEnd) <= maxQueryGap: # if exon start is near transcript end
                                                if strand == '+':
                                                    if exonRefStart > transcriptRefEnd: # ref positions do NOT overlap, because separated by intron
                                                        transcript.append(hit)
                                                        removeList.append(hit)
                                                else:
                                                    if exonRefStart < transcriptRefEnd:
                                                        transcript.append(hit)
                                                        removeList.append(hit)
                                            elif abs(exonQueryEnd - transcriptQueryStart) <= maxQueryGap: # if exon end is near transcript start
                                                if abs(exonQueryEnd - transcriptQueryStart) < maxIntronSize:
                                                    if strand == '+':
                                                        if exonRefEnd < transcriptRefStart:
                                                            transcript.insert(0,hit)
                                                            removeList.append(hit)
                                                    else:
                                                        if exonRefEnd > transcriptRefStart:
                                                            transcript.insert(0,hit)
                                                            removeList.append(hit)
                                for item in removeList:
                                    isoList.remove(item)
                                    

                            ##### build second isoform, beginning with same start as first transcript #####                                                      
                            transcriptStart = int(transcript[0][7])
                            transcript2 = None
                            if currentIso == 'iso1':
                                isoList = list(hits[1])
                            else:
                                isoList = list(hits[0])
                            for hit in isoList:                            
                                if (int(hit[8]) - int(hit[7])) > 0:
                                    hitStrand = '+'
                                else:
                                    hitStrand = '-'
                                if (int(hit[7]) < (transcriptStart + maxRefStartEndOffset)) and (int(hit[7]) > (transcriptStart - maxRefStartEndOffset)): # exon close to reference start pos
                                    if hitStrand == strand: # same strand as other transcript
                                        transcript2 = list([hit])
                                        isoList.pop(isoList.index(hit))
                            if transcript2: # found exon with same start pos
                                removeList = ['space_holder']
                                while (len(removeList) > 0) and (len(isoList) > 0): 
                                    removeList = []
                                    for hit in isoList:
                                        if (int(hit[8]) - int(hit[7])) > 0:
                                            hitStrand = '+'
                                        else:
                                            hitStrand = '-'
                                        if hitStrand == strand: # if hit is on same strand
                                            exonQueryStart = int(hit[5])
                                            exonRefStart = int(hit[7])
                                            transcriptRefEnd = int(transcript2[-1][8])
                                            transcriptQueryEnd = int(transcript2[-1][6])                            
                                            if abs(exonQueryStart - transcriptQueryEnd) <= maxQueryGap: # if exon start is near transcript end
                                                if abs(exonRefStart - transcriptRefEnd) < maxIntronSize:
                                                    if strand == '+':
                                                        if exonRefStart > transcriptRefEnd: # ref positions do NOT overlap, because separated by intron
                                                            transcript2.append(hit)
                                                            removeList.append(hit)
                                                    else:
                                                        if exonRefStart < transcriptRefEnd: 
                                                            transcript2.append(hit)
                                                            removeList.append(hit)
                                    for item in removeList:
                                        isoList.remove(item)






                            #####################################
                            # compare alternative splicing models
                            #####################################
                            keepGene = False
                            splicingInfo = None
                            spliceModel = 'NA'
                            if transcript and transcript2: # both transcripts assembled
                                if ( abs(int(transcript[-1][8]) - int(transcript2[-1][8])) < maxRefStartEndOffset ): # reference end the same (and start)
                                    isoVersion = tpmDict[gene][isoform][-1] # at this point get wild vs dom version
                                    queryMatches = 0 # next, checking if the hits are all the same query start and end positions 
                                    refMatches = 0
                                    overlappingHits = 0
                                    transcript1length = 0
                                    transcript2length = 0
                                    for exon in range(len(transcript)):
                                        transcript1length += abs(int(transcript[exon][6]) - int(transcript[exon][5]))
                                    for exon in range(len(transcript2)):
                                        transcript2length += abs(int(transcript2[exon][6]) - int(transcript2[exon][5]))
                                    if len(transcript) != len(transcript2): # different number of hits; model 1 or 2 
                                        keepGene = True
                                        # regardless of currentIso, isoVersion corresponds to transcript, not transcript2
                                        if transcript1length > transcript2length: # current transcript has longer alignment length                                            
                                            if isoVersion == '1238':
                                                if len(transcript) < len(transcript2): # intron retention
                                                    spliceModel = 1
                                                if len(transcript) > len(transcript2): # exon skipping
                                                    spliceModel = 10
                                            else:
                                                spliceModel = 2
                                        else:
                                            if isoVersion == '1238':
                                                spliceModel = 2
                                            else:
                                                if len(transcript) > len(transcript2):
                                                    spliceModel = 1
                                                if len(transcript) < len(transcript2): 
                                                    spliceModel = 10

                                    else: # same number of hits
                                        for exon in range(len(transcript)):
                                            if abs(int(transcript[exon][6]) - int(transcript[exon][5])) == abs(int(transcript2[exon][6]) - int(transcript2[exon][5])):
                                                queryMatches += 1
                                            if (transcript[exon][7] == transcript2[exon][7]) and (transcript[exon][8] == transcript2[exon][8]):
                                                refMatches += 1
                                            if (max(transcript[exon][7:9]) < min(transcript2[exon][7:9])) or (min(transcript[exon][7:9]) > max(transcript2[exon][7:9])): # checking, haphazardly, if the exons DON'T overlap
                                                overlappingHits += 1
                                            transcript1length += abs(int(transcript[exon][6]) - int(transcript[exon][5]))
                                            transcript2length += abs(int(transcript2[exon][6]) - int(transcript2[exon][5]))
                                        if refMatches != len(transcript): # different reference positions; model 3 or 4, or 9
                                            if overlappingHits > 0: # if there is a non-overlapping hit (substitution). Regarless of wild v dom version, at least one of the exon skipping is real splicing
                                                spliceModel = 9
                                                keepGene = True
                                            elif queryMatches != len(transcript): # if the query lengths are not all the same (and not substitution)
                                                keepGene = True
                                                if transcript1length > (transcript2length + minSpliceSiteDifference): # current transcript alignment longer
                                                    if isoVersion == '1238':
                                                        spliceModel = 3
                                                    else:
                                                        spliceModel = 4
                                                elif transcript1length < (transcript2length - minSpliceSiteDifference):
                                                    if isoVersion == '1238':
                                                        spliceModel = 4
                                                    else:
                                                        spliceModel = 3
                                                else:
                                                    spliceModel = 3.5 # same alingment lengths. confusing (probably not good splicing)
                                        else: # same reference positions; models 5-8 
                                            if queryMatches != len(transcript): 
                                                models = [] # may have different splicing model for each blast hit  
                                                for exon in range(len(transcript)):                                                 
                                                    refAlnLength = abs(int(transcript[exon][7]) - int(transcript[exon][8])) # same in both isoforms
                                                    queryLength1 = abs(int(transcript[exon][5]) - int(transcript[exon][6]))
                                                    queryLength2 =  abs(int(transcript2[exon][5]) - int(transcript2[exon][6]))
                                                    # NOTE: query end - query start is never equal to 'alignment length' in the blastn output
                                                    # this is because of gaps in the alignment, I'm assuming
                                                    # therefore we want to use the query start and end positions for this part
                                                    newMod = 'NA'
                                                    if isoVersion == 'HA89':
                                                        if (queryLength1 < (queryLength2 - minSpliceSiteDifference)) and ( abs(queryLength2-refAlnLength) <= (float(queryLength2)/50) ): # arbitrary max difference between longer hit and the reference alignment length
                                                            if queryLength1 < (refAlnLength - minSpliceSiteDifference): # a criterion for this model is that one transcript is shorter than reference alignment length
                                                                newMod = 5
                                                            if queryLength2 > (refAlnLength + minSpliceSiteDifference):
                                                                newMod = 7
                                                        elif queryLength1 > (queryLength2 + minSpliceSiteDifference) and ( abs(queryLength1-refAlnLength) <= (float(queryLength1)/50) ):
                                                            if queryLength2 < (refAlnLength - minSpliceSiteDifference):
                                                                newMod = 6
                                                            if queryLength1 > (refAlnLength + minSpliceSiteDifference):
                                                                newMod = 8
                                                    elif isoVersion == '1238':
                                                        if queryLength1 > (queryLength2 + minSpliceSiteDifference) and ( abs(queryLength1-refAlnLength) <= (float(queryLength1)/50) ):
                                                            if queryLength2 < (refAlnLength - minSpliceSiteDifference):
                                                                newMod = 5
                                                            if queryLength1 > (refAlnLength + minSpliceSiteDifference):
                                                                newMod = 7
                                                        elif queryLength1 < (queryLength2 - minSpliceSiteDifference) and ( abs(queryLength2-refAlnLength) <= (float(queryLength2)/50) ):
                                                            if queryLength1 < (refAlnLength - minSpliceSiteDifference):
                                                                newMod = 6
                                                            if queryLength2 > (refAlnLength + minSpliceSiteDifference):
                                                                newMod = 8
                                                    models.append(newMod)
                                                while 'NA' in models: models.remove('NA') # getting rid of NA's 
                                                if len(models) > 0:
                                                    combined = ''.join(map(str,models))
                                                    if combined == '5':
                                                        spliceModel = 5
                                                        keepGene = True
                                                    if combined[0] == '6':
                                                        spliceModel = 6
                                                        keepGene = True                                                
                                                    if len(combined) > 1: # this transcript's exons match multiple models
                                                        if models.count(5) > (len(combined)/2):                                                 
                                                            spliceModel = 5                                                    
                                                            keepGene = True
                                                        if models.count(6) > (len(combined)/2):
                                                            spliceModel = 6
                                                            keepGene = True

                            # if outputting this gene, add all exons to the output list
                            if keepGene == True: 
                                if gene not in outputDict:
                                    outputDict[gene] = [ [transcript, transcript2, spliceModel] ]
                                else:
                                    # lastly, want to make sure the current transcript doesn't share any blast hits with other transcripts in the dictionary 
                                    # in other words, if we successfully contructed multiple different transcripts sharing ref positions, we assume the blast is poor, not outputting
                                    newAS = True        
                                    for tran in [transcript, transcript2]:
                                        for newHit in tran:
                                            for AS in outputDict[gene]:
                                                for oldTranscript in AS[0:2]:
                                                    for oldHit in oldTranscript:
                                                        if newHit == oldHit:
                                                            newAS = False
                                    if newAS == True:
                                        outputDict[gene].append( [transcript, transcript2, spliceModel] )







    ##############################
    #            OUTPUT  
    ############################## 
    for gene in outputDict:
        if outputType == '1':
            outLine = [gene]
            for scenario in range(len(outputDict[gene])):
                mod = str(outputDict[gene][scenario][2])
                outLine.append(mod)
            print '\t'.join(outLine)

        if outputType == '2':
            for trans in outputDict[gene]:
                transcript = trans[0]
                transcript2 = trans[1]
                print
                print
                print
                print gene, trans[2]
                print '    iso1'
                print '\t'.join(['qseqid', 'sseqid', 'pident', 'qlen' ,'length', 'qstart', 'qend', 'sstart', 'send', 'evalue'])
                for exon in transcript:
                    print '\t'.join(exon)
                print '    iso2'
                for exon in transcript2:
                    print '\t'.join(exon)




    
main()
