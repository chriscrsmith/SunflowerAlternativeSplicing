
# current version 01052016
# this version of the script is looking for clear examples of intron retention/deletion in alternative splicing
# takes as input 1 the filtered fpkm file, 2 the Trinity.fasta, and 3 the output type (1 for list, 2 for visualization)

# messy. Instead, want to:
# go through filtered_fpkms, once, make collection of genes and isoforms
#     so, here, want to collect infomation in a certain form: c100002_g1_i1. 
#     so, gene c100002_g2, with 2 indicies: [0] isoform info c100002_g1_i1, and [1] node info
# then go through the dom version file to get dom version
# THEN go through the fasta to get the node info

def main():

    import sys

    #open files
    filtered_fpkm_path = sys.argv[1]
    fasta_file_name = sys.argv[2]
    output_type = int(sys.argv[3]) # e.g. 1 or 2: 1 is a list of good looking genes. 2 attemps to visualize those node alignments
    filtered_fpkms = open(filtered_fpkm_path, 'r')
    fasta_file = open(fasta_file_name, "r")

    # making list of fpkm filtered genes, to help distinguish wild and dom versions
    # this little bit of code is only necessary for the subsequent section of code, 
    # but is a little redundant with the gene_dict I think
    # also need this for determining which isoforms to analyze (e.g. some isoforms don't survive filtering)
    filtered_gene_list = {}
    for line in filtered_fpkms:
        gene = line.strip().split()[0]
        iso = '_'.join([ gene, line.strip().split()[1] ])
        if gene not in filtered_gene_list:
            filtered_gene_list[gene] = [[iso],[]] # first index of each dict value is a list of isos, second node info, third the dom version
        else:
            filtered_gene_list[gene][0].append(iso)

    # dictionary with isoform node info
    for line in fasta_file:
        if line[0] == ">":
            gene = "_".join([line[1:].split("_")[0],line[1:].split("_")[1]]) # does not include isoform
            if gene in filtered_gene_list:
                iso = line.split()[0][1:]
                if iso in filtered_gene_list[gene][0]:
                    filtered_gene_list[gene][1].append(line.strip())

    # check that they share first and last nodes
    good_isoforms = {}
    first_and_last_count = 0
    for gene in filtered_gene_list:
        num_isos = len( filtered_gene_list[gene][0] )
        if num_isos == 2: 
            info1 = filtered_gene_list[gene][1][0].split("[")[1][0:-1].split(" ")
            info2 = filtered_gene_list[gene][1][1].split("[")[1][0:-1].split(" ")
            nodes1 = []
            nodes2 = []
            for node in info1:
                new_node = node.split(":")[0]
                nodes1.append(new_node)
            for node in info2:
                new_node = node.split(":")[0]
                nodes2.append(new_node)
            if nodes1[0] == nodes2[0]:
                if nodes1[-1] == nodes2[-1]:
                    first_and_last_count += 1

                    # aligning the introns and looking for intron retetion/deletion
                    iso1_dict = {} # I think this part is counting the number of times each node occurs
                    for node in nodes1:
                        if node in iso1_dict:
                            iso1_dict[node] += 1
                        else:
                            iso1_dict[node] = 1
                    iso2_dict = {}
                    for node in nodes2:
                        if node in iso2_dict:
                            iso2_dict[node] += 1
                        else:
                            iso2_dict[node] = 1
                    align1 = []
                    align2 = []
                    ind_iso1 = 0
                    ind_iso2 = 0
                    still_aligning = True
                    mismatch = False
                    while still_aligning == True:
                        node1 = nodes1[ind_iso1]
                        node2 = nodes2[ind_iso2]
                        shared1 = False
                        shared2 = False
                        if node1 == node2:
                            align1.append(node1)
                            align2.append(node2)
                            ind_iso1 += 1
                            ind_iso2 += 1
                            dict_num = iso1_dict[node1]
                            iso1_dict[node1] = dict_num - 1
                            dict_num = iso2_dict[node2]
                            iso2_dict[node2] = dict_num - 1
                        else:
                            if node1 in nodes2:
                                shared1 = True
                            if node2 in nodes1:
                                shared2 = True
                            if (shared1 == True) and (shared2 == False):
                                align1.append('-')
                                align2.append(node2)
                                ind_iso2 += 1
                                dict_num = iso2_dict[node2]
                                iso2_dict[node2] = dict_num - 1
                            if (shared1 == False) and (shared2 == True):
                                align1.append(node1)
                                align2.append('-')
                                ind_iso1 += 1
                                dict_num = iso1_dict[node1]
                                iso1_dict[node1] = dict_num - 1
                            if (shared1 == False) and (shared2 == False):
                                align1.append(node1)
                                align2.append(node2)
                                ind_iso1 += 1
                                ind_iso2 += 1
                                mismatch = True
                                still_aligning = False # stopping there, we're not interested in this gene
                            if (shared1 == True) and (shared2 == True):
                                mismatch = True                            
                                still_aligning = False # not interested in this gene
                        if (ind_iso1 == len(nodes1)) and (ind_iso2 == len(nodes2)):
                            still_aligning = False
                    if mismatch == False:
                        good_isoforms[gene] = 1
                        if output_type == 2:
                            print 
                            print gene
                            print '\t'.join(align1)
                            print '\t'.join(align2)

    # output
    if output_type == 1:
        for gene in good_isoforms:
            print gene
    
            
main()
