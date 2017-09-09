import sys
import copy
import rmsd
import numpy as np
import bisect

# parameters
regionChr = 47
regionStart = 2699520
regionEnd = 154931043

#regionChr = 4
#regionStart = 0
#regionEnd = 243199373

removeQuantile = 0.05
halfWindow = 5e5

# load structure data
inputData = []
numOfStructures = len(sys.argv[2:])
for inputFileName in sys.argv[2:]:
    sys.stderr.write('reading file: '+inputFileName+'\n')
    inputFile = open(inputFileName,"r")
    inputData.append({})
    for inputFileLine in inputFile:
        inputFileLineData = inputFileLine.strip().split()
        if int(inputFileLineData[0]) == regionChr and int(inputFileLineData[1]) >= regionStart and int(inputFileLineData[1]) <= regionEnd:
            inputData[-1][(int(inputFileLineData[0]),int(inputFileLineData[1]))] = [float(inputFileLineData[2]),float(inputFileLineData[3]),float(inputFileLineData[4])]

# find common particles
commonLoci = set(inputData[0])
for inputStructure in inputData[1:]:
    commonLoci = commonLoci.intersection(set(inputStructure))
numOfLoci = len(commonLoci)
commonLoci = list(commonLoci)
commonData = []
for inputStructure in inputData:
    commonData.append([])
    for commonLocus in commonLoci:
        commonData[-1].append(inputStructure[commonLocus])
sys.stderr.write('found common particles: '+str(numOfLoci)+'\n')

# load contact data
sys.stderr.write('reading contacts: '+sys.argv[1]+'\n')
inputContactFile = open(sys.argv[1],"r")
contactData = {}
for inputLine in inputContactFile:
    inputLineData = inputLine.strip().split("\t")
    for i in [2,5]: # phase data
        inputLineData[i] = int(inputLineData[i])
    for i in [1,4]: # position data, for future averaging
        inputLineData[i] = float(inputLineData[i])
    leftChr = int(inputLineData[0])*2 + inputLineData[2]
    rightChr = int(inputLineData[3])*2 + inputLineData[5]
    if leftChr != regionChr and rightChr != regionChr:
        continue
    # add to contact list
    if not leftChr in contactData:
        contactData[leftChr] = []
    if not rightChr in contactData:
        contactData[rightChr] = []
    contactData[leftChr].append(inputLineData[1])
    contactData[rightChr].append(inputLineData[4])

# sort contact data
for contactChr in contactData:
    contactData[contactChr] = sorted(contactData[contactChr])
    
# count number of legs to determine the min number of legs
countData = []
for commonLocus in commonLoci:
    inputPdbChr = commonLocus[0]
    inputPdbPos = commonLocus[1]
    numOfContacts = bisect.bisect_right(contactData[inputPdbChr][bisect.bisect_left(contactData[inputPdbChr],inputPdbPos - halfWindow):],inputPdbPos + halfWindow)
    countData.append(numOfContacts)
    #sys.stdout.write(str(numOfContacts)+'\n')
minContacts = np.percentile(countData, removeQuantile*100.0)
sys.stderr.write('min number of legs: '+str(minContacts)+'\n')

# output clean data
for i in range(len(countData)):
    if countData[i] < minContacts:
        continue
    sys.stdout.write(str(commonLoci[i][0])+'\t'+str(commonLoci[i][1]))
    for j in range(numOfStructures):
        sys.stdout.write('\t'+str(inputData[j][commonLoci[i]][0])+'\t'+str(inputData[j][commonLoci[i]][1])+'\t'+str(inputData[j][commonLoci[i]][2]))
    sys.stdout.write('\n')