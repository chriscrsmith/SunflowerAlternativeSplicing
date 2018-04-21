
# with this script I want to give input 1 the interpolated cM SNPs txt file, and input 2 a RIL vcf, and output a table of SNPs and wild vs dom for each
# output a single list (for now) wild vs dom for each SNP, or none, or hetero?

def main():

    import sys

    snps_path = sys.argv[1]
    RIL_path = sys.argv[2]
    SNP_file = open(snps_path, "r")
    RIL_file = open(RIL_path, "r")
    

    SNP_dict = {}
    SNP_file.readline() # getting rid of header
    headerLine = ['RIL']
    for line in SNP_file:
        newline = line.strip().split()
        chrom = newline[0]
        if chrom[0] == '0':
            chrom = chrom[1]
        cM = newline[2]
        contig = newline[3]
        contig_pos = newline[4]
        key = '_'.join( [contig, contig_pos] )
        HA89version = newline[5]
        version1238 = newline[6]
        info = [HA89version, version1238, chrom, cM]
        SNP_dict[key] = list(info) # use list() to avoid linking variables
        
#    for snp in SNP_dict: # doing it like this to maintain the order in the dictionary
 #       headerLine.append( '_'.join([SNP_dict[snp][2], SNP_dict[snp][3]]) )
  #  print '\t'.join(headerLine)

#    for snp in SNP_dict: # doing it like this to maintain the order in the dictionary                     
 #       headerLine.append(SNP_dict[snp][2])
  #  print '\t'.join(headerLine)

    for snp in SNP_dict: # doing it like this to maintain the order in the dictionary
        headerLine.append(SNP_dict[snp][3])
    print '\t'.join(headerLine)

    # now, genotype the RIL for each SNP
    count1238 = 0
    HA89count = 0
    heterocount = 0
    tot = 0
    nuc_list = ['A', 'T', 'C', 'G']
    RIL_line = []
    for line in RIL_file:
        parent_type = None
        if line[0:2] != "##": # skipping headers
            if line.split()[0] == '#CHROM':
                RIL_name = line.split()[9].split('.')[0]
                RIL_line.append(RIL_name)
            else:
                new_line = line.strip().split("\t")
                location = "_".join(new_line[0:2])
                if location in SNP_dict: # if this snp is fixed between parents
                    new_ref = new_line[3]
                    new_alt = new_line[4]
                    if (new_ref in nuc_list) and (new_alt in nuc_list): # filtering indels
                        if int(new_line[9].split(':')[2]) < 20: # filterig bad quality snps for the current RIL
                            parentType = 'NA'
                        else:
                            new_type = new_line[9].split(":")[0]
                            HA89version = SNP_dict[location][0]
                            version1238 = SNP_dict[location][1]
                            if new_type == "0/0":
                                genotype = new_ref
                            elif new_type == "1/1":
                                genotype = new_alt
                            elif new_type == "0/1":
                                genotype = new_ref + '/' + new_alt
                            if (new_type == "0/0") or (new_type == "1/1"): #homozygote
                                if genotype == HA89version:
                                    parentType = 'HA89'
                                elif genotype == version1238:
                                    parentType = '1238'
                                else:
                                    parentType = 'NA' # non parent genotype
                            else: # heterozygote
                                alleleList = genotype.split('/')
                                parentList = [HA89version, version1238]
                                if (alleleList[0] in parentList) and (alleleList[1] in parentList):
                                    parentType = 'hetero'
                                else:
                                    parentType = 'NA' # non parent genotype
                        SNP_dict[location].append(parentType)

    for snp in SNP_dict:
        if len(SNP_dict[snp]) == 5:
            RIL_line.append(SNP_dict[snp][4])
        else:
            RIL_line.append('NA')

    percentHetero = float(RIL_line.count('hetero')) / float(len(RIL_line))
    percentHA89 = float(RIL_line.count('HA89')) / float(len(RIL_line))
    percent1238 = float(RIL_line.count('1238')) / float(len(RIL_line))
    percentNA = float(RIL_line.count('NA')) / float(len(RIL_line))
    # other output
#    print '\t'.join( ['hetero', 'HA89', '1238', 'NA', 'total'] )
 #   print '\t'.join( [str(percentHetero), str(percentHA89), str(percent1238), str(percentNA), str(percentHetero +  percentHA89 +percent1238+  percentNA)] )

    # genotyping output
#    if (percentHetero < .4) and (percent1238 > .2):
 #       print '\t'.join(RIL_line)    

    SNP_file.close()
    RIL_file.close()

main()
