import sys
import copy
from scipy import interpolate
import math
import numpy as np

# input:
# output:

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
inputPdbFile=open(sys.argv[2],"r")
outputFile=open(sys.argv[3],"w")
sex=sys.argv[4]

# parameters
#maxDistance = 10 # default
#minFraction = 0.8 # default
maxDistance = 20
minFraction = 0.0 # useless
minRatio = 0.5 # min ratio of shortest and 2nd shortest distances
minSeparationFactor = 1.0 # times resolution (e.g. 100 kb) for the smallest separation allowed for purely unphased contacts
cleanDistance = 10e6 # for GM12878
#cleanDistance = 20e6 # for blood
minAgreeCounts = 2

# PAR (as in hg19)
xName = 23
yName = 24
sexNames = [24,23] # pat=Y, mat=X
parNames = [0,-1]
# format: parPositions[PAR ID = 1 or 2][chr = X or Y][start or end]
parPositions = [[[60000,2699520],[10000,2649520]],[[154931043,155260560],[59034049,59363566]]]
# format: parShifts[PAR ID = 1 or 2][chr = X or Y]
parShifts = [] # shift of real pos relative to PAR pos (zero for X)
for i in range(len(parNames)):
    parShifts.append([0,parPositions[i][1][0]-parPositions[i][0][0]])

# function to clean imputed data
def cleanImputedContact (contactData, cleanDistance, minAgreeCounts):
    sqrtCleanDistance = math.sqrt(cleanDistance)
    cleanData = {}
    for leftChr in contactData:
        cleanData[leftChr] = {}
        for rightChr in contactData[leftChr]:
            cleanData[leftChr][rightChr] = []
            #print 'cleaning '+str(leftChr)+', '+str(rightChr)
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
    
    
# function to find interpolated 3D position of a locus
def findPos (chrName, locus, phase, interpolatedStructure, minLocus, maxLocus):
    outLocus = locus
    if sex == "m" and chrName in parNames:
        parId = parNames.index(chrName)
        outChr = sexNames[phase]*2+phase
        outLocus += parShifts[parId][1-phase]
        #sys.stderr.write("original: "+str(chrName)+"/"+str(phase)+":"+str(locus)+"\tmale: "+str(outChr)+":"+str(outLocus)+'\n')
    else:
        outChr = chrName*2+phase
    outLocus = max(outLocus, minLocus[outChr])
    outLocus = min(outLocus, maxLocus[outChr])
    return interpolatedStructure[outChr](outLocus)

# determine resolution from structure data
incrementCounts = {}
previousChr = -1
previousPosValue = -1
for inputPdbLine in inputPdbFile:
    inputPdbLineData = inputPdbLine.strip().split()
    inputChr = int(inputPdbLineData[0])
    inputPosValue = int(inputPdbLineData[1])
    if inputChr == previousChr:
        increment = inputPosValue - previousPosValue
        if increment not in incrementCounts:
            incrementCounts[increment] = 0
        incrementCounts[increment] += 1
    previousChr = inputChr
    previousPosValue = inputPosValue
resolution = max(incrementCounts.iterkeys(), key=(lambda key: incrementCounts[key]))
sys.stderr.write("resolution: "+str(resolution)+"\n")
minSeparation = minSeparationFactor * resolution

# load structure data
inputPdbFile.seek(0)
inputPdbData = {}
for inputPdbLine in inputPdbFile:
    inputPdbLineData = inputPdbLine.strip().split()
    inputPdbChr = int(inputPdbLineData[0])
    inputPdbLocus = int(inputPdbLineData[1])
    inputPdbPos = [float(inputPdbLineData[2]),float(inputPdbLineData[3]),float(inputPdbLineData[4])]
    if not inputPdbChr in inputPdbData:
        inputPdbData[inputPdbChr]=[]
    inputPdbData[inputPdbChr].append((inputPdbLocus,inputPdbPos))
    
# linear interpolate
inputPdbLinearInterpolate = {}
inputPdbMinLocus = {}
inputPdbMaxLocus = {}
for inputPdbChr in inputPdbData:
    inputPdbLinearInterpolate[inputPdbChr] = interpolate.interp1d(np.array([x[0] for x in inputPdbData[inputPdbChr]]),np.array([x[1] for x in inputPdbData[inputPdbChr]]),axis=0)
    inputPdbMinLocus[inputPdbChr] = min([x[0] for x in inputPdbData[inputPdbChr]])
    inputPdbMaxLocus[inputPdbChr] = max([x[0] for x in inputPdbData[inputPdbChr]])

# assign phase to each contact
contactData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in range(6):
        inputLineData[i] = int(inputLineData[i])
    isInter = 0
    if inputLineData[0] != inputLineData[3]:
        isInter = 1
        
    # for male
    if sex == "m":
        inputLineData[0], inputLineData[2] = convertMale(inputLineData[0], inputLineData[1], inputLineData[2])        
        inputLineData[3], inputLineData[5] = convertMale(inputLineData[3], inputLineData[4], inputLineData[5])
    
    # for female
    if sex == "f" and (inputLineData[0] == yName or inputLineData[3] == yName):
        continue
    
    # find candidate pairs
    candidateSquaredDistance = []
    for hap1 in range(2) if inputLineData[2] < 0 else [inputLineData[2]]:
        pos1 = findPos(inputLineData[0], inputLineData[1], hap1, inputPdbLinearInterpolate, inputPdbMinLocus, inputPdbMaxLocus)
        for hap2 in range(2) if inputLineData[5] < 0 else [inputLineData[5]]:
            pos2 = findPos(inputLineData[3], inputLineData[4], hap2, inputPdbLinearInterpolate, inputPdbMinLocus, inputPdbMaxLocus)
            candidateSquaredDistance.append((hap1,hap2,(pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2))
    #print candidateSquaredDistance
    
    # find shortest distance
    candidateSquaredDistance = sorted(candidateSquaredDistance, key=lambda x:x[2])
    bestCandidateSquaredDistance = candidateSquaredDistance[0]
    if bestCandidateSquaredDistance[2] == 0:
        #print "too close" (a problem with dedup that leads to two same legs in a contact)
        continue
    if len(candidateSquaredDistance) >= 2 and bestCandidateSquaredDistance[2]/candidateSquaredDistance[1][2] > minRatio:
        #print "too uncertain"
        continue
    
    # output to stderr violation status
    sys.stdout.write(inputLine.strip()+'\t'+str(isInter)+'\t'+str(len(candidateSquaredDistance))+'\t'+str(math.sqrt(bestCandidateSquaredDistance[2]))+'\t'+str(1/bestCandidateSquaredDistance[2]/sum(1/x[2] for x in candidateSquaredDistance))+'\n')
    
    # remove intra contacts that are unphased and too local
    if inputLineData[0] == inputLineData[3] and abs(inputLineData[1] - inputLineData[4]) < minSeparation and max(inputLineData[2],inputLineData[5]) < 0:
        continue    

    # assign potential contacts
    if bestCandidateSquaredDistance[2] > maxDistance**2:
        #print "too far away"
        continue
    if 1/bestCandidateSquaredDistance[2]/sum(1/x[2] for x in candidateSquaredDistance) < minFraction:
        #print "too uncertain"
        continue
        
    # shift if contact is in male PAR
    if sex == "m" and inputLineData[0] in parNames:
        parId = parNames.index(inputLineData[0])
        inputLineData[0] = sexNames[bestCandidateSquaredDistance[0]]
        inputLineData[1] += parShifts[parId][1-bestCandidateSquaredDistance[0]]
    if sex == "m" and inputLineData[3] in parNames:
        parId = parNames.index(inputLineData[3])
        inputLineData[3] = sexNames[bestCandidateSquaredDistance[1]]
        inputLineData[4] += parShifts[parId][1-bestCandidateSquaredDistance[1]]
    # add to contact list
    if not inputLineData[0] in contactData:
        contactData[inputLineData[0]] = {}
    if not inputLineData[3] in contactData[inputLineData[0]]:
        contactData[inputLineData[0]][inputLineData[3]] = []
    contactData[inputLineData[0]][inputLineData[3]].append([inputLineData[0],inputLineData[1],bestCandidateSquaredDistance[0],inputLineData[3],inputLineData[4],bestCandidateSquaredDistance[1]])

# sort data
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        contactData[leftChr][rightChr] = sorted(contactData[leftChr][rightChr])

# clean data
contactData = cleanImputedContact(contactData, cleanDistance, minAgreeCounts)

# print results
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        for contact in contactData[leftChr][rightChr]:
            outputFile.write(str(contact[0])+'\t'+str(contact[1])+'\t'+str(contact[2])+'\t'+str(contact[3])+'\t'+str(contact[4])+'\t'+str(contact[5])+'\n')
