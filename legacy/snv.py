import sys
import copy

# input:
# output:

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
sampleIndex=int(sys.argv[2])

# parameters
minMapQ = 40
minReads = 5
minDistance = 100
minBulkDepth = 15
minBulkAltDepth = 5
minBulkBalance = 0.3
minBulkOccurrence = 2

# read file and analyze SNVs
potentialData = []
numOfBulkData = 0
numOfDetectedBulkData = 0
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
    genotypes = []
    bulkCounts = [0,0]
    bulkOccurrences = [0,0]
    for i in range(len(inputFileLineData)-9):
        countData = inputFileLineData[i+9].split(":")[1:]
        FCountData = countData[0].split(",")
        RCountData = countData[1].split(",")
        refCountF = int(FCountData[0])
        refCountR = int(RCountData[0])
        altCountF = int(FCountData[1])
        altCountR = int(RCountData[1])
        genotype = 0 # ref only
        if altCountF + altCountR > 0:
            genotype = 1 # alt present
            if altCountF + altCountR >= minReads:# and altCountF > 0 and altCountR > 0:
                genotype = 2 # alt called
        genotypes.append(genotype)
        # sum all cells as the bulk
        bulkCounts[0] += refCountF + refCountR
        bulkCounts[1] += altCountF + altCountR
        if refCountF + refCountR > 0:
            bulkOccurrences[0] += 1
        if altCountF + altCountR > 0:
            bulkOccurrences[1] += 1
    #print currentChr, currentPos, genotypes, bulkCounts, bulkOccurrences
    # call potential de novo
    if genotypes[sampleIndex] == 2 and genotypes.count(2) == 1 and genotypes.count(1) == 0:
        potentialData.append([currentChr, currentPos, refAllele, altAllele])
    # determine bulk SNPs
    if bulkOccurrences[0] >= minBulkOccurrence and bulkOccurrences[1] >= minBulkOccurrence and bulkCounts[1] >= minBulkAltDepth and bulkCounts[0] + bulkCounts[1] >= minBulkDepth and float(bulkCounts[0])/(bulkCounts[0] + bulkCounts[1]) >= minBulkBalance and float(bulkCounts[1])/(bulkCounts[0] + bulkCounts[1]) >= minBulkBalance:
        numOfBulkData += 1
        if genotypes[sampleIndex] == 2:
            numOfDetectedBulkData += 1

# calculate FNR
sys.stdout.write(str(numOfBulkData)+'\t'+str(numOfDetectedBulkData)+'\t'+str(1-float(numOfDetectedBulkData)/numOfBulkData)+'\n')

# second pass to remove clustered SNVs
for i in range(len(potentialData)):
    if i > 0:
        if potentialData[i][1] - potentialData[i-1][1] < minDistance:
            continue
    if i < len(potentialData) - 1:
        if potentialData[i+1][1] - potentialData[i][1] < minDistance:
            continue        
    sys.stdout.write(potentialData[i][0]+"\t"+str(potentialData[i][1])+"\t"+potentialData[i][2]+"\t"+potentialData[i][3]+"\n")