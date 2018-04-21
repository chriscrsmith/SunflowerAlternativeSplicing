

# this quick n dirty script is simply taking a list of isoforms, and outputting the longest of every pair of isoforms


def main():

    import sys
    inpath = sys.argv[1]
    infile = open(inpath, 'r')
    length_path = sys.argv[2]
    lengths = open(length_path, 'r')

    gene_dict = {}
    infile.readline()
    for line in infile:
        newline = line.strip().split()
        gene = newline[0]
        gene_dict[gene] = newline[4:6]

    length_dict = {}
    for line in lengths:
        newline = line.strip().split()
        iso = newline[0].split('>')[1]
        gene = '_'.join(newline[0].split('_')[0:2]).split('>')[1]
        if iso in gene_dict[gene]:
            length = int(newline[1].split('=')[1])
            gene_dict[gene].remove(iso)
            if len(gene_dict[gene]) == 1: # if this is the first iso for this gene, set the longest length to this iso
                length_dict[gene] = [iso, length]
            if len(gene_dict[gene]) == 0: # if this is the second iso, then see which lenght is longer
                if length > int(length_dict[gene][1]):
                    length_dict[gene] = [iso, length]
        
    # print
    for gene in length_dict:
        print '\t'.join([ length_dict[gene][0], str(length_dict[gene][1]) ])
            



    infile.close()
    lengths.close()

main()
