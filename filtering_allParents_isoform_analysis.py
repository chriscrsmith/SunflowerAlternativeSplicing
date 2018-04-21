
# this script is for filtering based on fpkm cutoff
# specifically, it's totalling the fpkm across all 6 parents, and if that number is less than the cutoff it's not included

def main():

    import sys
    #open file
    in_path = sys.argv[1]
    infile = open(in_path, 'r')
    threshold = float(sys.argv[2])

    header = infile.readline()
    print header.strip()

    #fill dictionary with isoforms and FPKM values
    for line in infile:
        new_line = line.strip().split("\t")
        total_FPKM = 0
        for i in range(6):
            total_FPKM += float(new_line[i+2])
        if total_FPKM >= threshold:
            print line.strip()
            
    infile.close()

main()
