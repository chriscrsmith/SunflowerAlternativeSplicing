

# this script will take as input an rQTL table, an
#   output a table after adding 1 to each fpkm score 
#   and taking ratio of i1 to i2 for each gene


def main():

    import sys
    rqtl_path = sys.argv[1]
    rqtl_file = open(rqtl_path, 'r')
    
    #saving header, and getting nunber of isos
    header = rqtl_file.readline()
    header_split = header.strip().split()
    iso_list = []
    for iso in header_split:
        if iso[0] == 'c':
            iso_list.append(iso)

    # output header
    second_part = header_split[len(iso_list):]
    gene_list = []
    for iso in iso_list:
        gene = '_'.join(iso.split('_')[0:2])
        if gene not in gene_list:
            gene_list.append(gene)
    print '\t'.join(gene_list)

    # getting past the chromosome and position rows (blank for fpkm columns)
    rqtl_file.readline()
    rqtl_file.readline()

    # going through each row
    for line in rqtl_file:
        current_isos = line.strip().split()[:len(iso_list)]
        
        # adding 1 to all isos
        new_isos = []
        for iso in current_isos:
            new_iso = float(iso) + 1
            new_isos.append(new_iso)            

        # taking ratio of i1 and i2
        final_isos = []
        first_or_second = 1
        iso1 = 0
        iso2 = 0
        newline_list = []
        for iso in new_isos:
            if first_or_second ==1:
                iso1 = iso
                first_or_second = 2
            else:
                iso2 = iso
                first_or_second = 1
                ratio =  iso1 / iso2
                newline_list.append(str(ratio))
        print '\t'.join(newline_list)


main()
