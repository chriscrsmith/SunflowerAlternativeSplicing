

# this script will take a fasta file, and a list of contigs to look for, and simply grab those sequences.

def main():

    import sys

    fastaPath = sys.argv[1]
    listPath = sys.argv[2]
    fasta = open(fastaPath, 'r')
    contigList = open(listPath, 'r')

    contigDict = {}
    for line in contigList:
        contigDict[line.strip()] = 0

    newHeader = False
    for line in fasta:
        newline = line.strip()
        if newHeader == True:
            newHeader = False
            if contig in contigDict:
                seq = line.strip()
                print '>'+contig
                print seq
        else:
            newHeader = True
            contig = newline[1:].split()[0]
main()


