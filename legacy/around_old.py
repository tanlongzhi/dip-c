import sys
import copy
import numpy
import time

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")

# parameters
maxDistance = 10000000 # 10mb
maxGrids = maxDistance//gridSize # = 1e3
maxEdgeDistance = maxDistance + gridSize/2
totalGrids = maxGrids+1
aroundCounts = [x[:] for x in [[0] * totalGrids] * totalGrids] 

# read contacts
contactData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in [2,5]: # phase data
        inputLineData[i] = int(inputLineData[i])
    for i in [1,4]: # position data, for future averaging
        inputLineData[i] = float(inputLineData[i])
    # swap two partners of a contact if out of order
    if inputLineData[0] > inputLineData[3] or (inputLineData[0] == inputLineData[3] and inputLineData[1] > inputLineData[4]):
        newContact = [inputLineData[4], inputLineData[5], inputLineData[1], inputLineData[2], 1]
    else:
        newContact = [inputLineData[1], inputLineData[2], inputLineData[4], inputLineData[5], 1]
    # add to contact list and merge
    if not inputLineData[0] in contactData:
        contactData[inputLineData[0]] = {}
    if not inputLineData[3] in contactData[inputLineData[0]]:
        contactData[inputLineData[0]][inputLineData[3]] = []
    contactData[inputLineData[0]][inputLineData[3]].append(newContact)

# sort data
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        contactData[leftChr][rightChr] = sorted(contactData[leftChr][rightChr])

# find nearby contacts
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        if leftChr == rightChr:
            continue # skip intra or skip inter
        sys.stderr.write('cleaning '+leftChr+', '+rightChr+'\n')
        
        numOfContacts = len(contactData[leftChr][rightChr])
        for i in range(numOfContacts):
            for j in range(i+1, numOfContacts):
                if contactData[leftChr][rightChr][j][1] - contactData[leftChr][rightChr][i][1] > maxEdgeDistance:
                    break
                if abs(contactData[leftChr][rightChr][j][4] - contactData[leftChr][rightChr][i][4])) <= maxEdgeDistance:
                    relativePos = [contact2[0] - contact1[0], contact2[2] - contact1[2]]
                    roundedPos = [int((relativePos[0]+maxDistance)//gridSize),int((relativePos[1]+maxDistance)//gridSize)]
                    if roundedPos[0] >= 0 and roundedPos[0] < totalGrids and roundedPos[1] >= 0 and roundedPos[1] < totalGrids:
                        aroundCounts[roundedPos[0]][roundedPos[1]] += 1
                        aroundCounts[roundedPos[1]][roundedPos[0]] += 1
            for j in range(0, i):
                if contactData[leftChr][rightChr][i][1] - contactData[leftChr][rightChr][j][1] > maxEdgeDistance:
                    break
                if abs(contactData[leftChr][rightChr][j][4] - contactData[leftChr][rightChr][i][4])) <= maxEdgeDistance:




# sort and merge contacts within maxDistance
for chr1 in contactData:
    for chr2 in contactData[chr1]:
        if chr1 == chr2:
            continue # skip intra or skip inter
        print chr1+', '+chr2
        for contact1 in contactData[chr1][chr2]:
            for contact2 in contactData[chr1][chr2]:
                relativePos = [contact2[0] - contact1[0], contact2[2] - contact1[2]]
                roundedPos = [int((relativePos[0]+maxDistance)//gridSize),int((relativePos[1]+maxDistance)//gridSize)]
                if roundedPos[0] >= 0 and roundedPos[0] < totalGrids and roundedPos[1] >= 0 and roundedPos[1] < totalGrids:
                    aroundCounts[roundedPos[0]][roundedPos[1]] += 1
                    aroundCounts[roundedPos[1]][roundedPos[0]] += 1

# output results
for i in range(totalGrids):
    outputFile.write("\t".join(map(str,aroundCounts[i]))+"\n")