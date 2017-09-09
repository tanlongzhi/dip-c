import sys
import copy
import math

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")

# parameters
cleanDistance = 10e6
sqrtCleanDistance = math.sqrt(cleanDistance)
minAgreeCounts = 2

# load contact data
contactData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in [1,2,4,5]:
        inputLineData[i] = int(inputLineData[i])
    # add to contact list
    if not inputLineData[0] in contactData:
        contactData[inputLineData[0]] = {}
    if not inputLineData[3] in contactData[inputLineData[0]]:
        contactData[inputLineData[0]][inputLineData[3]] = []
    contactData[inputLineData[0]][inputLineData[3]].append(inputLineData)

# clean data
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        print leftChr+', '+rightChr
        numOfContacts = len(contactData[leftChr][rightChr])
        contactData[leftChr][rightChr] = sorted(contactData[leftChr][rightChr])
        for i in range(numOfContacts):
            agreeCounts = 0
            disagreeCounts = 0
            for j in range(i+1, numOfContacts):
                if agreeCounts >= minAgreeCounts:
                    break
                if contactData[leftChr][rightChr][j][1] - contactData[leftChr][rightChr][i][1] > cleanDistance:
                    break
                if math.sqrt(contactData[leftChr][rightChr][j][1] - contactData[leftChr][rightChr][i][1]) + math.sqrt(abs(contactData[leftChr][rightChr][j][4] - contactData[leftChr][rightChr][i][4])) <= sqrtCleanDistance:
                    if contactData[leftChr][rightChr][j][2] == contactData[leftChr][rightChr][i][2] and contactData[leftChr][rightChr][j][5] == contactData[leftChr][rightChr][i][5]:
                        agreeCounts += 1
                    else:
                        disagreeCounts += 1
            for j in range(0, i):
                if agreeCounts >= minAgreeCounts:
                    break
                if contactData[leftChr][rightChr][i][1] - contactData[leftChr][rightChr][j][1] > cleanDistance:
                    break
                if math.sqrt(contactData[leftChr][rightChr][i][1] - contactData[leftChr][rightChr][j][1]) + math.sqrt(abs(contactData[leftChr][rightChr][j][4] - contactData[leftChr][rightChr][i][4])) <= sqrtCleanDistance:
                    if contactData[leftChr][rightChr][j][2] == contactData[leftChr][rightChr][i][2] and contactData[leftChr][rightChr][j][5] == contactData[leftChr][rightChr][i][5]:
                        agreeCounts += 1
                    else:
                        disagreeCounts += 1
            #print nearbyCounts
            if agreeCounts >= minAgreeCounts:           
                outputFile.write(contactData[leftChr][rightChr][i][0]+'\t'+str(contactData[leftChr][rightChr][i][1])+'\t'+str(contactData[leftChr][rightChr][i][2])+'\t'+contactData[leftChr][rightChr][i][3]+'\t'+str(contactData[leftChr][rightChr][i][4])+'\t'+str(contactData[leftChr][rightChr][i][5])+'\n')