# get the reverse-complement DNA sequence


def ReverseComplement1(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])


def main():

    import sys 
    infile = open(sys.argv[1], 'r')
    
    for line in infile:
        print ReverseComplement1(line.strip())

main()
