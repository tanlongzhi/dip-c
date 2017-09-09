import sys
import copy

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")

# parameters
imputeDistance = 2e6
minVotes = 3
minVoteFraction = 0.9
minScore = 1e-5

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
            numOfContacts = len(contactData[leftChr][rightChr])
            weightMatrix = [[0.0 for i in range(numOfContacts)] for j in range(numOfContacts)]
            for i in range(numOfContacts):
                for j in range(numOfContacts):
                    if i == j:
                        continue
                    weightMatrix[i][j] = abs(contactData[leftChr][rightChr][i][1] - contactData[leftChr][rightChr][j][1]) + abs(contactData[leftChr][rightChr][i][4] - contactData[leftChr][rightChr][j][4]) + 1.0 # add one to avoid zeros
                    weightMatrix[i][j] = 1.0/weightMatrix[i][j]
        
            for i in range(numOfContacts):
                contact = contactData[leftChr][rightChr][i]
                if contact[2] >= 0 and contact[5] >= 0: # fully phased
                    imputeData.append([contact[0],contact[1],contact[2],contact[3],contact[4],contact[5]])
                elif contact[2] >= 0: # left side phased, impute right side
                    scoreValue = 0.0 # for paternal is positive, for maternal is negative
                    for j in range(numOfContacts):
                        nearbyContact = contactData[leftChr][rightChr][j]
                        if (nearbyContact[2] == contact[2] or nearbyContact[2] < 0) and nearbyContact[5] >= 0:
                            scoreValue += weightMatrix[i][j] * (1 - 2 * nearbyContact[5])
                    if scoreValue > minScore:
                        imputeData.append([contact[0],contact[1],contact[2],contact[3],contact[4],0])
                    elif scoreValue < -minScore:
                        imputeData.append([contact[0],contact[1],contact[2],contact[3],contact[4],1])
                else: # right side phased, impute left side
                    scoreValue = 0.0 # for paternal is positive, for maternal is negative
                    for j in range(numOfContacts):
                        nearbyContact = contactData[leftChr][rightChr][j]
                        if (nearbyContact[5] == contact[5] or nearbyContact[5] < 0) and nearbyContact[2] >= 0:
                            scoreValue += weightMatrix[i][j] * (1 - 2 * nearbyContact[2])
                    if scoreValue > minScore:
                        imputeData.append([contact[0],contact[1],0,contact[3],contact[4],contact[5]])
                    elif scoreValue < -minScore:
                        imputeData.append([contact[0],contact[1],1,contact[3],contact[4],contact[5]])
                
# print results
for imputeContact in imputeData:
    outputFile.write(imputeContact[0]+'\t'+str(imputeContact[1])+'\t'+str(imputeContact[2])+'\t'+imputeContact[3]+'\t'+str(imputeContact[4])+'\t'+str(imputeContact[5])+'\n')