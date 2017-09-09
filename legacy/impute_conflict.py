import sys
import copy
import random

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")

# parameters


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
imputeData = {}
for leftChr in contactData:
    if not leftChr in imputeData:
        imputeData[leftChr] = {}
    for rightChr in contactData[leftChr]:
        if not rightChr in imputeData[leftChr]:
            imputeData[leftChr][rightChr] = []
            
        if leftChr == rightChr: # intra
            continue
            
        else: # inter
        
            # initialization (use negative numbers to indicate imputed phases)
            for contact in contactData[leftChr][rightChr]:
                if contact[2] < 0:
                    contact[2] = random.randrange(1) - 2
                if contact[5] < 0:
                    contact[5] = random.randrange(1) - 2
                imputeData[leftChr][rightChr].append(contact)
                
            # calculate pairwise weights
            weightMatrix = [[0.0 for i in range(len(imputeData[leftChr][rightChr]))] for j in range(len(imputeData[leftChr][rightChr]))]
            for i in range(len(imputeData[leftChr][rightChr])):
                for j in range(len(imputeData[leftChr][rightChr])):
                    if i == j:
                        continue
                    weightMatrix[i][j] = abs(imputeData[leftChr][rightChr][i][1] - imputeData[leftChr][rightChr][j][1]) + abs(imputeData[leftChr][rightChr][i][4] - imputeData[leftChr][rightChr][j][4])
                    weightMatrix[i][j] = 1.0/weightMatrix[i][j]
                
            # start optimization
            while True:
                changedPhase = False
                choiceList = []
                for i in range(len(imputeData[leftChr][rightChr])):
                    # possible choices
                    if contact[2] < 0:
                        leftChoices = [0-2, 1-2]
                    else:
                        leftChoices = [imputeData[leftChr][rightChr][i][2]]
                    if contact[5] < 0:
                        rightChoices = [0-2, 1-2]
                    else:
                        rightChoices = [imputeData[leftChr][rightChr][i][2]]
                    # evaluate each choice
                    scoreValue = 0.0
                    for leftChoice in leftChoices:
                        for rightChoice in rightChoices:
                            for j in range(len(imputeData[leftChr][rightChr])):
                                scoreSign = 0
                                if (leftChoice == imputeData[leftChr][rightChr][j][2] or leftChoice + 2 == imputeData[leftChr][rightChr][j][2]) and (rightChoice == imputeData[leftChr][rightChr][j][5] or rightChoice + 2== imputeData[leftChr][rightChr][j][5]):
                                    scoreSign = 1
                                scoreValue += scoreSign * weightMatrix[i][j]
                            choiceList.append([leftChoice, rightChoice, scoreValue])
                    # make a choice
                    bestChoise = max(choiceList, key=lambda x:x[2])
                    if bestChoise[0] != imputeData[leftChr][rightChr][i][2] or bestChoise[1] != imputeData[leftChr][rightChr][i][5]:
                        changedPhase = True
                        imputeData[leftChr][rightChr][i][2] = bestChoise[0]
                        imputeData[leftChr][rightChr][i][5] = bestChoise[1]
                if not changedPhase:
                    break

            print imputeData[leftChr][rightChr]