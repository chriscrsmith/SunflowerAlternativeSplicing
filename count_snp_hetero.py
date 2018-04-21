
# NOTE: format requirement: the input file must contain all 3 header lines
# it will take in a preRQTL-formatted file (rows of RILs, cols of SNPs)
# it will output a similarly formatted file, but with only snps that pass a heterozygosity threshold


def main():

    import sys

    data_path = sys.argv[1]
    data_file = open(data_path, "r")

    SNPmatrix = []
    header = data_file.readline().strip().split()
    for snp in header[1:]: # skipping one column, which are RIL id's
        SNPmatrix.append( [0,0,0,0,snp,0,0] ) # snp matrix has length equal to number of snps. 
                                      # And each entry has 4 values: hetero, HA89, 1238, NA
    secondLine = data_file.readline().strip().split()
    thirdLine = data_file.readline().strip().split()

    ril_count = 0
    for line in data_file:
        ril_count += 1
        newline = line.strip().split()
        for snp in range(len(newline)-1): # skipping first col
            value = newline[snp+1] 
            if value == 'hetero':
                SNPmatrix[snp][0] += 1
            if value == 'HA89':
                SNPmatrix[snp][1] += 1
            if value == '1238':
                SNPmatrix[snp][2] += 1
            if value == 'NA':
                SNPmatrix[snp][3] += 1

    # now filtering snps based on heterozygosity, and parent type makeup, and NA's
    filtered_snps = {}
    for snp in SNPmatrix:
        if float(snp[0])/ril_count <= .2: # at most 20% heterozygotes
            if float(snp[1])/ril_count <= .95: # at most 95% HA89
                if float(snp[2])/ril_count >= .05: # at least 5% 1238
                    if float(snp[3])/ril_count < .5: # less than 50% NA
                        filtered_snps[snp[4]] = 0
    data_file.close()

    # printing filtered genotype table
    data_file = open(data_path, "r")
    for line in data_file:
        currentLine = line.strip().split()
        newline = [currentLine[0]] #ril name
        for snp in range( len(SNPmatrix) ):
            if SNPmatrix[snp][4] in filtered_snps:
                newline.append( currentLine[snp+1] ) # this works for the header line, too
        print '\t'.join(newline)
    data_file.close()


main()
