
# this script is finding snps in the parent vcf that are fixed between populations

def main():

    import sys
    vcfpath = sys.argv[1]
    vcf = open(vcfpath, 'r')

    # header
    print ' '.join([ 'contig', 'pos', 'HA89version', '1238version'])

    nuclist = ['A', 'C', 'T', 'G']
    dp = -1 # index of DP field

    for line in vcf:
        if line[0] != '#':
            newline = line.strip().split()
            info_field = newline[7]
            A1238info = newline[9].split(':')
            b1238info = newline[10].split(':')
            e1238info = newline[11].split(':')
            HA89Ainfo = newline[12].split(':')
            HA89binfo = newline[13].split(':')
            HA89einfo = newline[14].split(':')

            A1238 = A1238info[0]
            b1238 = b1238info[0]
            e1238 = e1238info[0]
            HA89A = HA89Ainfo[0]
            HA89b = HA89binfo[0]
            HA89e = HA89einfo[0]

            if (A1238 == '0/0') and (b1238 == '0/0') and (e1238 == '0/0') and (HA89A == '1/1') and (HA89b == '1/1') and (HA89e == '1/1'):
                if (int(A1238info[dp]) >= 50) and (int(b1238info[dp]) >= 50) and (int(e1238info[dp]) >= 50) and (int(HA89Ainfo[dp]) >= 50) and (int(HA89binfo[dp]) >= 50) and (int(HA89einfo[dp]) >= 50):
                    if (newline[3] in nuclist) and (newline[4] in nuclist):
                        HA89version = newline[4]
                        version1289 = newline[3]
                        print newline[0], newline[1], HA89version, version1289

            if (A1238 == '1/1') and (b1238 == '1/1') and (e1238 == '1/1') and (HA89A == '0/0') and (HA89b == '0/0') and (HA89e == '0/0'):
                if (int(A1238info[dp]) >= 50) and (int(b1238info[dp]) >= 50) and (int(e1238info[dp]) >= 50) and (int(HA89Ainfo[dp]) >= 50) and (int(HA89binfo[dp]) >= 50) and (int(HA89einfo[dp]) >= 50):
                    if (newline[3] in nuclist) and (newline[4] in nuclist):
                        HA89version = newline[3]
                        version1289 = newline[4]
                        print newline[0], newline[1], HA89version, version1289



    vcf.close()

main()
