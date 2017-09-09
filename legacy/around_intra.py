import sys
import copy
import numpy
import time

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")

# parameters
gridSize = 10000 # 10kb
maxDistance = 10000000 # 10mb
maxGrids = maxDistance//gridSize # =1e3
totalGrids = 2*maxGrids+1

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

# sort and merge contacts within maxDistance
aroundCounts = [x[:] for x in [[0] * totalGrids] * totalGrids] 
for chr1 in contactData:
    for chr2 in contactData[chr1]:
        contactData[chr1][chr2] = sorted(contactData[chr1][chr2])
        if chr1 != chr2:
            continue # skip inter
        print chr1+', '+chr2
        for i in range(len(contactData[chr1][chr2])):
            contact1 = contactData[chr1][chr2][i]
            for j in range(i,len(contactData[chr1][chr2])):
                contact2 = contactData[chr1][chr2][j]
                relativePos = [contact2[0] - contact1[0], contact2[2] - contact1[2]]
                if relativePos[0] > maxDistance+gridSize//2:
                    break
                roundedPos = [int((relativePos[0]+maxDistance)//gridSize),int((relativePos[1]+maxDistance)//gridSize)]
                if roundedPos[0] >= 0 and roundedPos[0] < totalGrids and roundedPos[1] >= 0 and roundedPos[1] < totalGrids:
                    aroundCounts[roundedPos[0]][roundedPos[1]] += 1
                    aroundCounts[roundedPos[1]][roundedPos[0]] += 1
            for j in range(i-1,-1,-1):
                contact2 = contactData[chr1][chr2][j]
                relativePos = [contact2[0] - contact1[0], contact2[2] - contact1[2]]
                if -relativePos[0] > maxDistance+gridSize//2:
                    break
                roundedPos = [int((relativePos[0]+maxDistance)//gridSize),int((relativePos[1]+maxDistance)//gridSize)]
                if roundedPos[0] >= 0 and roundedPos[0] < totalGrids and roundedPos[1] >= 0 and roundedPos[1] < totalGrids:
                    aroundCounts[roundedPos[0]][roundedPos[1]] += 1
                    aroundCounts[roundedPos[1]][roundedPos[0]] += 1

# output results
for i in range(totalGrids):
    outputFile.write("\t".join(map(str,aroundCounts[i]))+"\n")