import sys
import copy
import random

# input:
# output:

# read IO locations from arguments
inputAssignFile=open(sys.argv[1],"r")
inputRawFile=open(sys.argv[2],"r")
outputFile=open(sys.argv[3],"w")

# parameters
intraDistance = 10e6
phasedProb = 0.0923

# load fully phased (assigned) contacts
assignedContacts = {}
for inputAssignLine in inputAssignFile:
    inputAssignLineData = inputAssignLine.strip().split()
    assignedContacts[(inputAssignLineData[0],inputAssignLineData[1],inputAssignLineData[3],inputAssignLineData[4])] = [inputAssignLineData[2], inputAssignLineData[5], False]
    # the last tag indicates whether this contact has already been written

# load raw contacts and simulate phasing
for inputRawLine in inputRawFile:
    inputRawLineData = inputRawLine.strip().split()
    inputRawLocus = (inputRawLineData[0],inputRawLineData[1],inputRawLineData[3],inputRawLineData[4])
    trueHap = []
    if inputRawLocus in assignedContacts: # assigned
        if assignedContacts[inputRawLocus][2]: # already written
            continue
        trueHap = [assignedContacts[inputRawLocus][0],assignedContacts[inputRawLocus][1]]
        assignedContacts[inputRawLocus][2] = True
    elif inputRawLineData[0] ==  inputRawLineData[3] and abs(int(inputRawLineData[1]) - int(inputRawLineData[4])) <= intraDistance: # assume intra
        randomHap = random.getrandbits(1)
        trueHap = [randomHap, randomHap]
    else:
        trueHap = [random.getrandbits(1), random.getrandbits(1)]
        
    # randomly discard phase info
    for i in range(2):
        if random.random() > phasedProb:
            trueHap[i] = -1
    
    # output
    outputFile.write(inputRawLineData[0]+'\t'+inputRawLineData[1]+'\t'+str(trueHap[0])+'\t'+inputRawLineData[3]+'\t'+inputRawLineData[4]+'\t'+str(trueHap[1])+'\n')
    
    
