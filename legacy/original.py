import sys
import copy
import random

# input:
# output:

# read IO locations from arguments
inputAssignFile=open(sys.argv[1],"r")
inputRawFile=open(sys.argv[2],"r")
outputFile=open(sys.argv[3],"w")

# load fully phased (assigned) contacts
assignedContacts = {}
for inputAssignLine in inputAssignFile:
    inputAssignLineData = inputAssignLine.strip().split()
    assignedContacts[(inputAssignLineData[0],inputAssignLineData[1],inputAssignLineData[3],inputAssignLineData[4])] = [inputAssignLineData[2], inputAssignLineData[5]]

# load raw contacts and find ones corresponding to assigned contacts
for inputRawLine in inputRawFile:
    inputRawLineData = inputRawLine.strip().split()
    inputRawLocus = (inputRawLineData[0],inputRawLineData[1],inputRawLineData[3],inputRawLineData[4])
    if inputRawLocus in assignedContacts: # assigned
        rawHap = [int(inputRawLineData[2]),int(inputRawLineData[5])]
        trueHap = [int(assignedContacts[inputRawLocus][0]),int(assignedContacts[inputRawLocus][1])]
        if (rawHap[0] < 0 or rawHap[0] == trueHap[0]) and (rawHap[1] < 0 or rawHap[1] == trueHap[1]): # phases agree
            outputFile.write(inputRawLineData[0]+'\t'+inputRawLineData[1]+'\t'+inputRawLineData[2]+'\t'+inputRawLineData[3]+'\t'+inputRawLineData[4]+'\t'+inputRawLineData[5]+'\n')
    
    
