

# filtering genes by ones found in all parents to some degree, that differ in population mean in each isoform in opposite directions
# and finally outputting proportions of each isoform per parent


def main():

    from scipy import stats
    import sys
    import numpy
    import math

    #open file
    enterfile = sys.argv[1] # fpkm table (filtered) 
    infile = open(enterfile, "r")
    header = infile.readline()
    assemblypath = sys.argv[2]
    assembly = open(assemblypath, 'r') # give it the single line trinity assmebly fasta
    clean_gene_path = sys.argv[3]
    clean_gene_list = open(clean_gene_path, 'r')

    #init some variables
    count_1 = 0
    count_2 = 0
    X1238_count = 0
    HA89_count = 0
    count_5 = 0
    count_6 = 0
    gene_list = {}
    gene_allParents = {}
    X1238A_dict = {}
    HA89b_dict = {}
    X1238b_dict = {}
    X1238e_dict = {}
    HA89e_dict = {}
    HA89A_dict = {}
    clean_isoforms = {}

    # here making a list (in the form of a dictionary) of clean genes
    clean_genes = {}
    clean_gene_list.readline() # to get header out of the way
    for line in clean_gene_list:
        gene = line.strip().split()[0]
        clean_genes[gene] = 1

    #first, we're looking at >0 vs =0 reads, filling up dictionaries of who has what gene, and which genes are found in all parents
    for line in infile:
        newline = line.strip().split("\t")
        gene = newline[0]

        #check if parents have it
        if newline[1+1] != "0.00": #plus one because of col organization, for first parent 1238A
            X1238A = 1 #X because it doesn't like you to init a variable starting with number
            if gene not in X1238A_dict:
                X1238A_dict[gene] = []
                X1238A_dict[gene].append(newline)
            else:
                X1238A_dict[gene].append(newline)
        else:
            X1238A = 0
        if newline[1+2] != "0.00":
            HA89b = 1
            if gene not in HA89b_dict:
                HA89b_dict[gene] = []
                HA89b_dict[gene].append(newline)
            else:
                HA89b_dict[gene].append(newline)
        else:
            HA89b = 0
        if newline[1+3] != "0.00":
            X1238b = 1
            if gene not in X1238b_dict:
                X1238b_dict[gene] = []
                X1238b_dict[gene].append(newline)
            else:
                X1238b_dict[gene].append(newline)
        else:
            X1238b = 0
        if newline[1+4] != "0.00":
            X1238e = 1
            if gene not in X1238e_dict:
                X1238e_dict[gene] = []
                X1238e_dict[gene].append(newline)
            else:
                X1238e_dict[gene].append(newline)
        else:
            X1238e = 0
        if newline[1+5] != "0.00":
            HA89e = 1
            if gene not in HA89e_dict:
                HA89e_dict[gene] = []
                HA89e_dict[gene].append(newline)
            else:
                HA89e_dict[gene].append(newline)
        else:
            HA89e = 0
        if newline[1+6] != "0.00":
            HA89A = 1
            if gene not in HA89A_dict:
                HA89A_dict[gene] = []
                HA89A_dict[gene].append(newline)
            else:
                HA89A_dict[gene].append(newline)
        else:
            HA89A = 0

        #total wild and dom scores                                                                                  
        dom = X1238A + X1238b + X1238e
        wild = HA89b + HA89e + HA89A

        # counting total genes                                                                                        
        if gene not in gene_list:
            gene_list[gene] = []
            gene_list[gene].append(newline)
        else:
            gene_list[gene].append(newline)

    for gene in gene_list:
        if (gene in X1238A_dict) and (gene in X1238b_dict) and (gene in X1238e_dict) and (gene in HA89b_dict) and (gene in HA89e_dict) and (gene in HA89A_dict): #this line SHOULD reduce analysis to genes that all parents have
            if (len(gene_list[gene]) > 1) and (len(gene_list[gene]) < 3):  #this line reduces the analysis to only genes with more than 1 isoform, but at most 2
                if gene in clean_genes:
                    gene_allParents[gene] = list(gene_list[gene]) # use list() to avoid linking the variables
                    
    # get a dictionary of isoforms lengths from the assembly
    assembly_dict = {}
    for line in assembly:
        if line[0] == '>':
            fasta_iso = line.split(' ')[0][1:]
            length = line.split(' ')[1].split('=')[1]
            assembly_dict[fasta_iso] = length

    # here, we are recalculating fpkm values, because we have to calculate fpkm in order to calcualte tpm,
    # and we need the value of 1 reads tpm in order to get rid of zeros, to avoid infinite values in the ILR transformation
    # we have already filtered for genes where every parent has expression; this has to do with 100% isoform 1 situaitons
    # output will be in the form of:
    # gene    number_of_isoforms
    # ind1    iso1 iso2 iso3 iso4 etc.
    # ind2    iso1 iso2 NA NA (if only two isoforms)
    # ind3
    # ind4
    # ind5
    # ind6
    tot1238A = 1046660 # these are total fpkm's for the parents. Using them to calculate TPM
    tot1238b = 1063381
    tot1238e = 986468
    totHA89A = 1078733
    totHA89b = 1101558
    totHA89e = 1124123
    total_fpkms = [tot1238A, totHA89b, tot1238b, tot1238e, totHA89e, totHA89A] # reordered corresponding to TPM table column order
    
    mappedReads = [ 181273670, 95039498, 229715226, 154401412, 117856890, 193329542 ] # number of mapped reads for this individual (i.e. number of reads aligned. Not assembly), determined from bowtie.bam's
    
    for gene in gene_allParents:
        num_isos = len(gene_allParents[gene])
        fpkm_dict = {}
        gene_line = 'gene' + '\t' + gene + '\t' + str(num_isos) + '\t' + ( ('NA'+'\t')*(10-3) ) 
        gene_line = '\t'.join([ gene_line[:-1], gene_allParents[gene][0][1], gene_allParents[gene][1][1] ]) # the last two being added on here are the wild and dom isoform versions
        print gene_line                                       
        fpkm_matrix = []

        for ind in range(6): # num indiviuals
            new_row = []
            for iso in range(12): # num isoforms + 2 for the two versions, wild v dom
                new_row.append('NA') # place holder is NA
            fpkm_matrix.append(new_row)

        for iso in gene_allParents[gene]:
            iso_id = '_'.join(iso[0:2])
            fpkms = iso[2:]
            fpkm_dict[iso_id] = fpkms

        for i in range(num_isos):
            iso = gene_allParents[gene][i]
            for ind in range(6):
                ind_total = 0
                # loop to get individual total for this gene
                fasta_iso = '_'.join([iso[0], iso[1]])
                L =  int(assembly_dict[fasta_iso] ) # isoform length
                C = 1 # read read
                N = mappedReads[ind] # number of mapped reads for this individual (i.e. number of reads aligned. Not assembly), determined from bowtie.bam file
                one_reads_fpkm = ( ( float(10**9) * C) / float(N * L) )
                one_reads_tpm = (one_reads_fpkm / total_fpkms[ind]) * 1000000
                original_tpm = float(fpkm_dict['_'.join(iso[0:2])][ind])
                ind_tpm = (original_tpm + one_reads_tpm)
                fpkm_matrix[ind][i] = ind_tpm


        # calculating proportion of total fpkm for each isoform (for each individual)
        for ind in fpkm_matrix:
            ind_total = 0
            for iso in ind:
                if iso != 'NA':
                    ind_total += iso
            new_line = ''
            for iso in ind:
                if iso != 'NA':
                    new_line += str(float(iso)/float(ind_total) )+ '\t'
                else:
                    new_line += 'NA' + '\t'
            print new_line[:-1]


    infile.close()
    #outfile.close()

main()
