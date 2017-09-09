import sys
import copy
import math

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")
sex=sys.argv[3] # m (male) or f (female)

# parameters
imputeDistance = 10e6 # for GM12878
#imputeDistance = 20e6 # for blood
minVotes = 3
minVoteFraction = 0.8
intraThreshold = 10e6
interThreshold = 100e6
cleanDistance = imputeDistance
minAgreeCounts = 2

# positions (0-based) of PARs in hg19 (parPositions[PAR ID = 1 or 2][chr = X or Y][start or end])
xName = "23"
yName = "24"
parNames = ["PAR1","PAR2"]
parPositions = [[[60000,2699520],[10000,2649520]],[[154931043,155260560],[59034049,59363566]]]
# if female: disregard
# if male: remove PAR, and set phase=1 for non-PAR X, phase=2 for non-PAR Y

# function to impute
def imputeContact (contactData, evidenceData, imputeDistance, minVotes, minVoteFraction, intraThreshold, interThreshold):
    sqrtImputeDistance = math.sqrt(imputeDistance)
    imputeData = {}
    for leftChr in contactData:
        imputeData[leftChr] = {}
        for rightChr in contactData[leftChr]:
            imputeData[leftChr][rightChr] = []
            print 'imputing '+leftChr+', '+rightChr
            isIntra = (leftChr == rightChr)
            
            # impute each contact
            for contact in contactData[leftChr][rightChr]:
                imputedContact = []
                if contact[2] >= 0 and contact[5] >= 0: # fully phased
                    imputedContact = [contact[0],contact[1],contact[2],contact[3],contact[4],contact[5]]
                elif isIntra and abs(contact[1] - contact[4]) <= intraThreshold: # assume intra-homolog
                    imputedContact = [contact[0],contact[1],max(contact[2],contact[5]),contact[3],contact[4],max(contact[2],contact[5])]
                elif contact[2] >= 0: # left side phased, impute right side
                    votes = [0, 0]
                    for nearbyContact in evidenceData[leftChr][rightChr]:
                        if isIntra and math.sqrt(abs(nearbyContact[1]-contact[1])) + math.sqrt(abs(nearbyContact[4]-contact[4])) > sqrtImputeDistance:
                            continue
                        if ~isIntra and abs(nearbyContact[1]-contact[1]) > imputeDistance and abs(nearbyContact[4]-contact[4]) > imputeDistance:
                            continue
                        if (nearbyContact[2] == contact[2] or nearbyContact[2] < 0) and nearbyContact[5] >= 0:
                            votes[nearbyContact[5]] += 1
                    if votes[0] >= minVoteFraction*(votes[0]+votes[1]) and votes[0] >= minVotes:
                        imputedContact = [contact[0],contact[1],contact[2],contact[3],contact[4],0]
                    elif votes[1] >= minVoteFraction*(votes[0]+votes[1]) and votes[1] >= minVotes:
                        imputedContact = [contact[0],contact[1],contact[2],contact[3],contact[4],1]
                else: # right side phased, impute left side
                    votes = [0, 0]
                    for nearbyContact in evidenceData[leftChr][rightChr]:
                        if isIntra and math.sqrt(abs(nearbyContact[1]-contact[1])) + math.sqrt(abs(nearbyContact[4]-contact[4])) > sqrtImputeDistance:
                            continue
                        if ~isIntra and abs(nearbyContact[1]-contact[1]) > imputeDistance and abs(nearbyContact[4]-contact[4]) > imputeDistance:
                            continue
                        if (nearbyContact[5] == contact[5] or nearbyContact[5] < 0) and nearbyContact[2] >= 0:
                            votes[nearbyContact[2]] += 1
                    if votes[0] >= minVoteFraction*(votes[0]+votes[1]) and votes[0] >= minVotes:
                        imputedContact = [contact[0],contact[1],0,contact[3],contact[4],contact[5]]
                    elif votes[1] >= minVoteFraction*(votes[0]+votes[1]) and votes[1] >= minVotes:
                        imputedContact = [contact[0],contact[1],1,contact[3],contact[4],contact[5]]
                            
                # discard if imputation failed
                if imputedContact == []:
                    continue
                
                # discard inter-homolog if too close
                if isIntra and imputedContact[2] != imputedContact[5] and abs(imputedContact[1] - imputedContact[4]) < interThreshold: # discard inter-homolog if too close
                    continue
                
                # store results
                imputeData[leftChr][rightChr].append(imputedContact)

    return imputeData
    
# function to clean imputed data
def cleanImputedContact (contactData, cleanDistance, minAgreeCounts):
    sqrtCleanDistance = math.sqrt(cleanDistance)
    cleanData = {}
    for leftChr in contactData:
        cleanData[leftChr] = {}
        for rightChr in contactData[leftChr]:
            cleanData[leftChr][rightChr] = []
            print 'cleaning '+leftChr+', '+rightChr
            numOfContacts = len(contactData[leftChr][rightChr])
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
                if agreeCounts >= minAgreeCounts:           
                    cleanData[leftChr][rightChr].append(contactData[leftChr][rightChr][i])
    return cleanData
    
# function to impute completely unphased data
def imputeEmptyContact (contactData, evidenceData, imputeDistance, minVotes, minVoteFraction):
    sqrtImputeDistance = math.sqrt(imputeDistance)
    imputeData = {}
    for leftChr in contactData:
        if leftChr not in evidenceData:
            continue
        imputeData[leftChr] = {}
        for rightChr in contactData[leftChr]:
            if rightChr not in evidenceData[leftChr]:
                continue
            imputeData[leftChr][rightChr] = []
            if leftChr == rightChr: # skip intra
                continue
            print 'imputing empty contacts '+leftChr+', '+rightChr
            
            for contact in contactData[leftChr][rightChr]:
                # count votes
                votes = [[0, 0], [0, 0]]
                for nearbyContact in evidenceData[leftChr][rightChr]:
                    if math.sqrt(abs(nearbyContact[1]-contact[1])) + math.sqrt(abs(nearbyContact[4]-contact[4])) > sqrtImputeDistance:
                        continue
                    votes[nearbyContact[2]][nearbyContact[5]] += 1
                # analyze votes
                haplotypeVotes = []
                totalVotes = 0
                for hap1 in range(2):
                    for hap2 in range(2):
                        totalVotes += votes[hap1][hap2]
                        haplotypeVotes.append((hap1,hap2,votes[hap1][hap2]))
                bestHaplotypeVotes = max(haplotypeVotes, key=lambda x:x[2])
                if bestHaplotypeVotes[2] < minVotes or float(bestHaplotypeVotes[2])/totalVotes < minVoteFraction:
                    continue
                imputeData[leftChr][rightChr].append([contact[0],contact[1],bestHaplotypeVotes[0],contact[3],contact[4],bestHaplotypeVotes[1]])
    return imputeData

# function to find which chromosome or PAR a leg belongs to (0, 1, or -1=none), and phase
def convertMale (chrName, pos, phase):
    maleChrName = chrName
    malePhase = phase
    sexChrId = -1
    if chrName == xName:
        sexChrId = 0
        malePhase = 1 # maternal unless in PAR
    elif chrName == yName:
        sexChrId = 1
        malePhase = 0 # paternal unless in PAR
    if sexChrId >= 0:
        for i in range(len(parNames)):
            if pos > parPositions[i][sexChrId][0] and pos <= parPositions[i][sexChrId][1]:
                maleChrName = parNames[i]
                malePhase = phase # in PAR
                break
    return maleChrName, malePhase
                    
# load contact data
# for male: denote each PAR as a separate chromosome ("PAR1" and "PAR2")
contactData = {}
emptyData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in [1,2,4,5]:
        inputLineData[i] = int(inputLineData[i])
        
    # for male
    if sex == "m":
        inputLineData[0], inputLineData[2] = convertMale(inputLineData[0], inputLineData[1], inputLineData[2])        
        inputLineData[3], inputLineData[5] = convertMale(inputLineData[3], inputLineData[4], inputLineData[5])
        if inputLineData[0] in parNames or inputLineData[3] in parNames: # skip PARs
            continue
    
    # add to contact list
    if not inputLineData[0] in contactData:
        contactData[inputLineData[0]] = {}
        emptyData[inputLineData[0]] = {}
    if not inputLineData[3] in contactData[inputLineData[0]]:
        contactData[inputLineData[0]][inputLineData[3]] = []
        emptyData[inputLineData[0]][inputLineData[3]] = []
    if inputLineData[2] < 0 and inputLineData[5] < 0: # purely unphased
        emptyData[inputLineData[0]][inputLineData[3]].append(inputLineData)
    else: # at least one leg phased
        contactData[inputLineData[0]][inputLineData[3]].append(inputLineData)

# sort data
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        contactData[leftChr][rightChr] = sorted(contactData[leftChr][rightChr])

# impute data
imputeData = imputeContact(contactData, contactData, imputeDistance, minVotes, minVoteFraction, intraThreshold, interThreshold)
imputeData = imputeContact(contactData, imputeData, imputeDistance, minVotes, minVoteFraction, intraThreshold, interThreshold)
imputeData = cleanImputedContact(imputeData, cleanDistance, minAgreeCounts)
imputeEmptyData = imputeEmptyContact(emptyData, imputeData, imputeDistance, minVotes, minVoteFraction)

# print results
for leftChr in imputeData:
    for rightChr in imputeData[leftChr]:
        for imputeContact in imputeData[leftChr][rightChr]:
            outputFile.write(imputeContact[0]+'\t'+str(imputeContact[1])+'\t'+str(imputeContact[2])+'\t'+imputeContact[3]+'\t'+str(imputeContact[4])+'\t'+str(imputeContact[5])+'\n')
        for imputeContact in imputeEmptyData[leftChr][rightChr]:
            outputFile.write(imputeContact[0]+'\t'+str(imputeContact[1])+'\t'+str(imputeContact[2])+'\t'+imputeContact[3]+'\t'+str(imputeContact[4])+'\t'+str(imputeContact[5])+'\n')