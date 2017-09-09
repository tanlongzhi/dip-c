import sys
import copy

# input:
# output:

reverseComplements = {"A":"T", "C":"G", "G":"C", "T":"A"}

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
minGenotype = int(sys.argv[2]) # 1=regular, 2=strand

# parameters
minMapQ = 40
minReads = 5
minReadsPerStrand = 2
maxErrorReads = 0
minDistance = 100
maxBulkErrorReads = 1
minBulkDepth = 15
minBulkAlleleDepth = 5
minBulkBalance = 0.3

# read file and analyze SNVs
potentialData = []
numOfBulkData = 0
numOfDetectedBulkData = 0
numOfHomBulkData = 0
numOfDetectedHomBulkData = 0
for inputFileLine in inputFile:
    if inputFileLine.startswith("#"):
        continue # skip header
    inputFileLineData = inputFileLine.strip().split()
    refAllele = inputFileLineData[3]
    altAllele = inputFileLineData[4]
    currentChr = inputFileLineData[0]
    currentPos = int(inputFileLineData[1])
    
    # keep only biallelic SNVs
    if refAllele not in ["A","C","G","T"]:
        continue
    if altAllele not in ["A","C","G","T"]:
        continue
    
    # skip low mapping quality
    mapQData = inputFileLineData[7][4:].split(",")
    if mapQData[0] != ".":
        if int(mapQData[0]) < minMapQ:
            continue
    if mapQData[1] != ".":
        if int(mapQData[1]) < minMapQ:
            continue
            
    # find genotypes
    # col (starting from 0) 9 is bulk, 10 is single cell
    isBulkSnp = 0 # 0=no, 1=het, 2=homo
    bulkCounts = [0, 0]
    for i in range(len(inputFileLineData)-9):
        countData = inputFileLineData[i+9].split(":")[1:]
        FCountData = countData[0].split(",")
        RCountData = countData[1].split(",")
        refCountF = int(FCountData[0])
        refCountR = int(RCountData[0])
        altCountF = int(FCountData[1])
        altCountR = int(RCountData[1])
        if i == 0: # bulk
            # determine bulk SNPs
            bulkCounts = [refCountF + refCountR, altCountF + altCountR]
            if bulkCounts[1] >= minBulkAlleleDepth and bulkCounts[0] + bulkCounts[1] >= minBulkDepth and bulkCounts[0] <= maxBulkErrorReads:
                isBulkSnp = 2
                numOfHomBulkData += 1
            elif bulkCounts[1] >= minBulkAlleleDepth and bulkCounts[0] >= minBulkAlleleDepth and bulkCounts[0] + bulkCounts[1] >= minBulkDepth and float(bulkCounts[0])/(bulkCounts[0] + bulkCounts[1]) >= minBulkBalance and float(bulkCounts[1])/(bulkCounts[0] + bulkCounts[1]) >= minBulkBalance:
                isBulkSnp = 1
                numOfBulkData += 1
        else:
            genotype = 0 # initial call: not called
            if altCountF + altCountR >= minReads and refCountF + refCountR <= maxErrorReads:
                genotype = 1 # regular call
                if altCountF >= minReadsPerStrand and altCountR >= minReadsPerStrand:
                    genotype = 2 # stranded call
            # detect bulk SNP
            if genotype >= minGenotype:
                if isBulkSnp == 1:
                    numOfDetectedBulkData += 1
                elif isBulkSnp == 2:
                    numOfDetectedHomBulkData += 1
            # call potential de novo
            if genotype >= minGenotype and bulkCounts[1] == 0 and bulkCounts[0] + bulkCounts[1] >= minBulkDepth:
                potentialData.append([currentChr, currentPos, refAllele, altAllele])

# calculate FNR (detected * 2 because sperm inherits half of the genome)
sys.stdout.write('HET\t'+str(numOfBulkData)+'\t'+str(numOfDetectedBulkData)+'\t'+str(1-float(numOfDetectedBulkData*2)/numOfBulkData)+'\n')
sys.stdout.write('HOM\t'+str(numOfHomBulkData)+'\t'+str(numOfDetectedHomBulkData)+'\t'+str(1-float(numOfDetectedHomBulkData)/numOfHomBulkData)+'\n')

# second pass to remove clustered SNVs
for i in range(len(potentialData)):
    if i > 0:
        if potentialData[i][0] == potentialData[i-1][0] and potentialData[i][1] - potentialData[i-1][1] < minDistance:
            continue
    if i < len(potentialData) - 1:
        if potentialData[i+1][1] - potentialData[i][1] < minDistance:
            continue
    # determine type
    if potentialData[i][2] > reverseComplements[potentialData[i][2]]:
        mutationType = reverseComplements[potentialData[i][2]] + reverseComplements[potentialData[i][3]]
    else:
        mutationType = potentialData[i][2] + potentialData[i][3]
    sys.stdout.write(potentialData[i][0]+"\t"+str(potentialData[i][1])+"\t"+potentialData[i][2]+"\t"+potentialData[i][3]+"\t"+mutationType+"\n")