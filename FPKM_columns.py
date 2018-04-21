
# this script should produce columns of FPKM phenotype data for rQTL analysis

# (1) list of genes we're insterested in (sig q values)
# (2) the current RIL's ILR transformed fpkm values (after adding 1 reads FPKM like for the quantitative analysis)
# (3) the filtered FPKM file, so that we grab the two isoforms that survive filtering, because there will be some genes with 4 isoforms or whatever but we are only interested in a specific 2 of them
# (4) single line trinity assembly, so that we can calculate 1 reads FPKM
# (5) number of mapped reads file, first call ril ID, second call number of mapped reads
# (6) output type: 1 for header, 2 for data

def main():

    import sys
    gene_path = sys.argv[1]
    gene_file = open(gene_path, 'r')
    FPKM_path = sys.argv[2]
    FPKM_file = open(FPKM_path, 'r')
    allParents_path = sys.argv[3]
    allParents_file = open(allParents_path, 'r')
    assemblypath = sys.argv[4]
    assembly = open(assemblypath, 'r') 
    readsPath = sys.argv[5]
    mappedReads = open(readsPath, 'r')
    outputType = sys.argv[6]

    # making dictionary of genes we're interested in
    gene_file.readline()
    gene_dict = {}
    for line in gene_file:
        gene = line.strip().split()[0]
        gene_dict[gene] = [ [], [] ] # first member the iso numbers, second the info

    # looking for the two isoforms that survive in the filtered FPKM file
    # (or else we'll end up with more than 2 isoforms sometimes)
    allParents_file.readline()
    for line in allParents_file:
        current = line.strip().split()
        gene = current[0]
        if gene in gene_dict:
            iso = current[1]
            gene_dict[gene][0].append(iso)

    # finally grabbing FPKM for each of the genes, and the SPECIFIC isoforms we want
    total_fpkm = 0
    FPKM_file.readline()
    for line in FPKM_file:
        current =  line.strip().split()
        gene_col = current[0].split('_')
        gene = '_'.join(gene_col[0:2])
        fpkm = current[5]
        total_fpkm += float(fpkm)
        if gene in gene_dict:
            iso = gene_col[2]
            if iso in gene_dict[gene][0]:
                gene_dict[gene][1].append(current)

    # get a dictionary of isoforms lengths from the assembly                        
    assembly_dict = {}
    for line in assembly:
        if line[0] == '>':
            fasta_iso = line.split(' ')[0][1:]
            length = line.split(' ')[1].split('=')[1]
            assembly_dict[fasta_iso] = length

    # get mapped read counts
    readCounts = {}
    for line in mappedReads:
        newline = line.strip().split()
        rilName, reads = newline
        readCounts[rilName] = reads
        
    # getting near output
    outlist = []
    gene_file.close()    
    gene_file = open(gene_path, 'r')
    gene_file.readline()
    for row in gene_file:
        gene = row.strip().split()[0]
        value = list(gene_dict[gene]) # use list() to avoid linking variables

    # recalculating fpkm values to get rid of zeros by adding 1 reads TPM (have to calculate fpkm first)
    for gene in gene_dict: 
        num_isos = 2
        fpkm_dict = {}
        gene_dict[gene].append([]) # this is a new list inside the dict value: for iso1fpkm, iso2fpkm, totalfpkm 
        rilID = '_'.join(FPKM_path.split('/')[-1].split('_')[:-1])
        for i in range(2): # 2 isoforms
            iso = gene_dict[gene][1][i] 
            old_tpm = float(iso[5])
            fasta_iso = iso[0]            
            L =  int(assembly_dict[fasta_iso] ) # isoform length
            C = 1 # read read
            N = float(readCounts[rilID]) # number of mapped reads for this individual (i.e. number of reads aligned. Not assembly), determined from bowtie.bam file
            one_reads_fpkm = ( ( float(10**9) * C) / (N * L) )
            one_reads_tpm = (one_reads_fpkm / total_fpkm) * 1000000
            ind_tpm = (old_tpm + one_reads_tpm)
            gene_dict[gene][2].append(ind_tpm)
            if i == 1:
                gene_dict[gene][2].append(gene_dict[gene][2][0] + gene_dict[gene][2][1]) # total tpm

    # now, convert the new fpkms to proportions of the total fpkm
    for gene in gene_dict:
        gene_dict[gene].append([]) # new list for prop iso1 and iso2
        for iso in range(2):
            prop = float(gene_dict[gene][2][iso]) / float(gene_dict[gene][2][2]) # iso / total fpkm
            gene_dict[gene][3].append(prop)
        
        # output, finally
        if outputType == '1': # header
            iso_id1 = str(gene_dict[gene][1][0][0])
            iso_id2 = str(gene_dict[gene][1][1][0])
            outlist.append(iso_id1)                                                                                 
            outlist.append(iso_id2) 
        elif outputType == '2': # data
            fpkm1 = str( gene_dict[gene][3][0] )
            fpkm2 = str( gene_dict[gene][3][1] )
            outlist.append(fpkm1)
            outlist.append(fpkm2)

    print '\t'.join(outlist)

    gene_file.close()
    FPKM_file.close()
    allParents_file.close()
    assembly.close()
    mappedReads.close()
            
main()
