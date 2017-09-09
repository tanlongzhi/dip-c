import sys
import copy
import math

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")

# parameters
imputeDistance = 10e6
sqrtImputeDistance = math.sqrt(imputeDistance)
minVotes = 3
minVoteFraction = 0.9

# load contact data
contactData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in [1,2,4,5]:
        inputLineData[i] = int(inputLineData[i])
    if inputLineData[2] < 0 and inputLineData[5] < 0: # purely unphased
        continue
    # add to contact list
    if not inputLineData[0] in contactData:
        contactData[inputLineData[0]] = {}
    if not inputLineData[3] in contactData[inputLineData[0]]:
        contactData[inputLineData[0]][inputLineData[3]] = []
    contactData[inputLineData[0]][inputLineData[3]].append(inputLineData)

# impute data
imputeData = []
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        if leftChr == rightChr: # intra
            for contact in contactData[leftChr][rightChr]:
                if (contact[2] == 1 and contact[5] == 0) or (contact[2] == 0 and contact[5] == 1): # not the same phase
                    imputeData.append([contact[0],contact[1],contact[2],contact[3],contact[4],contact[5]])
                else:
                    imputeData.append([contact[0],contact[1],max(contact[2],contact[5]),contact[3],contact[4],max(contact[2],contact[5])])
        else: # inter
            for contact in contactData[leftChr][rightChr]:
                if contact[2] >= 0 and contact[5] >= 0: # fully phased
                    imputeData.append([contact[0],contact[1],contact[2],contact[3],contact[4],contact[5]])
                elif contact[2] >= 0: # left side phased, impute right side
                    votes = [0] * 2
                    for nearbyContact in contactData[leftChr][rightChr]:
                        if math.sqrt(abs(nearbyContact[1]-contact[1])) + math.sqrt(abs(nearbyContact[4]-contact[4])) > sqrtImputeDistance:
                            continue
                        if (nearbyContact[2] == contact[2] or nearbyContact[2] < 0) and nearbyContact[5] >= 0:
                            votes[nearbyContact[5]] += 1
                    if votes[0] >= minVoteFraction*(votes[0]+votes[1]) and votes[0] >= minVotes:
                        imputeData.append([contact[0],contact[1],contact[2],contact[3],contact[4],0])
                    elif votes[1] >= minVoteFraction*(votes[0]+votes[1]) and votes[1] >= minVotes:
                        imputeData.append([contact[0],contact[1],contact[2],contact[3],contact[4],1])
                else: # right side phased, impute left side
                    votes = [0] * 2
                    for nearbyContact in contactData[leftChr][rightChr]:
                        if math.sqrt(abs(nearbyContact[1]-contact[1])) + math.sqrt(abs(nearbyContact[4]-contact[4])) > sqrtImputeDistance:
                            continue
                        if (nearbyContact[5] == contact[5] or nearbyContact[5] < 0) and nearbyContact[2] >= 0:
                            votes[nearbyContact[2]] += 1
                    if votes[0] >= minVoteFraction*(votes[0]+votes[1]) and votes[0] >= minVotes:
                        imputeData.append([contact[0],contact[1],0,contact[3],contact[4],contact[5]])
                    elif votes[1] >= minVoteFraction*(votes[0]+votes[1]) and votes[1] >= minVotes:
                        imputeData.append([contact[0],contact[1],1,contact[3],contact[4],contact[5]])
                    
                
# print results
for imputeContact in imputeData:
    outputFile.write(imputeContact[0]+'\t'+str(imputeContact[1])+'\t'+str(imputeContact[2])+'\t'+imputeContact[3]+'\t'+str(imputeContact[4])+'\t'+str(imputeContact[5])+'\n')