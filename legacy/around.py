import sys
import copy
import numpy
import time

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")

# parameters
maxDistance = 10e6 # 10mb
minIntraDistance = 50e6 # 50mb

# read contacts
contactData = {}
legData = {}
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
    
# sort data
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        contactData[leftChr][rightChr] = sorted(contactData[leftChr][rightChr])

# find nearby contacts
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        
        if leftChr != rightChr:
            continue # skip intra or skip inter
            
        sys.stderr.write('analyzing '+leftChr+', '+rightChr+'\n')
        
        numOfContacts = len(contactData[leftChr][rightChr])
        for i in range(numOfContacts):
            if leftChr == rightChr and abs(contactData[leftChr][rightChr][i][1] - contactData[leftChr][rightChr][i][4]) < minIntraDistance:
                continue # skip short intra
                
            for j in range(i+1, numOfContacts):
                if contactData[leftChr][rightChr][j][1] - contactData[leftChr][rightChr][i][1] > maxDistance:
                    break
                if abs(contactData[leftChr][rightChr][j][4] - contactData[leftChr][rightChr][i][4]) <= maxDistance:
                    relativePos = [abs(contactData[leftChr][rightChr][j][1]-contactData[leftChr][rightChr][i][1]), abs(contactData[leftChr][rightChr][j][4]-contactData[leftChr][rightChr][i][4])]
                    sys.stdout.write(str(relativePos[0])+'\t'+str(relativePos[1])+'\n')
            for j in range(0, i):
                if contactData[leftChr][rightChr][i][1] - contactData[leftChr][rightChr][j][1] > maxDistance:
                    break
                if abs(contactData[leftChr][rightChr][j][4] - contactData[leftChr][rightChr][i][4]) <= maxDistance:
                    relativePos = [abs(contactData[leftChr][rightChr][j][1]-contactData[leftChr][rightChr][i][1]), abs(contactData[leftChr][rightChr][j][4]-contactData[leftChr][rightChr][i][4])]
                    sys.stdout.write(str(relativePos[0])+'\t'+str(relativePos[1])+'\n')
