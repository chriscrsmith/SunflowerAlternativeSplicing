

# sister functions for interpolating bp positions given cM






# these lines are required to import a python function:
# import sys
# sys.path.insert(0, '/data5/AltSplicing09182015/Scripts/')
# import interpolateBpPos




def fillGeneticMap(filePath):
    import decimal
    decimal.getcontext().prec = 22
    geneticMap = {}
    with open(filePath, "r") as infile:
        infile.readline() # header
        for line in infile:
            newline = line.strip().split()
            if newline[2] != "NA":
                chrom, bp, cm = newline[0:3]
                if chrom[0] == "0":
                    chrom = chrom[1:]
                if chrom not in geneticMap:
                    geneticMap[chrom] = []
                geneticMap[chrom].append([float(bp), float(cm)])
    return geneticMap

def interpolateBP(chrom, cM, myMap):
    import decimal
    decimal.getcontext().prec = 22
    cM = float(cM)
    
    data = myMap[chrom]
    countPos = 0
    while data[countPos][1] < cM:
         countPos += 1
    # I want to linearly interpolate, e.g. if the cM position is 1/5 the way between two map positions,
    # then make the bp position 1/5 the way between the two bp positions
    cMdiff = data[countPos][1] - data[countPos-1][1]
    bpDiff = data[countPos][0] - data[countPos-1][0]
    cMdist = cM - data[countPos-1][1]
    prop = cMdist / cMdiff
    newPos = data[countPos-1][0] + (prop*bpDiff)
    return newPos






                  # test lines
#import sys
#geneticMap = fillGeneticMap(sys.argv[1])
#print interpolateBP(sys.argv[2], sys.argv[3], geneticMap)








