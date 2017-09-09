import sys
import copy
from scipy import interpolate
import numpy as np

# input:
# output:

# read IO locations from arguments
inputPdbFile=open(sys.argv[1],"r")
inputFile=open(sys.argv[2],"r")
outputFile=open(sys.argv[3],"w")

# parameters
expandFactor = 0.0
halfWindow = 5e5
minContacts = 15
pdb_format = '%-6.6s%5.1d %4.4s%s%3.3s %s%4.1d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  %10d\n'

# load structure data
inputPdbData = {}
for inputPdbLine in inputPdbFile:
    if not inputPdbLine.startswith("HETATM"):
        continue
    inputPdbLineData = inputPdbLine[6:].strip().split()
    inputPdbChr = int(inputPdbLineData[2].strip("C_"))
    inputPdbLocus = int(inputPdbLineData[-1])
    inputPdbPos = [float(inputPdbLineData[5]),float(inputPdbLineData[6]),float(inputPdbLineData[7])]
    if not inputPdbChr in inputPdbData:
        inputPdbData[inputPdbChr]=[]
    inputPdbData[inputPdbChr].append((inputPdbLocus,inputPdbPos))
    
# calculate mean position of each chromosomes
meanChrPos = {}
for inputPdbChr in inputPdbData:
    numOfMolecules = len(inputPdbData[inputPdbChr])
    meanChrPos[inputPdbChr] = [0.0, 0.0, 0.0]
    for inputPdbMolecule in inputPdbData[inputPdbChr]:
        for i in range(3):
            meanChrPos[inputPdbChr][i] += inputPdbMolecule[1][i]/numOfMolecules
            
# load contact data
contactData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in [2,5]: # phase data
        inputLineData[i] = int(inputLineData[i])
    for i in [1,4]: # position data, for future averaging
        inputLineData[i] = float(inputLineData[i])
    leftChr = int(inputLineData[0])*2 + inputLineData[2]
    rightChr = int(inputLineData[3])*2 + inputLineData[5]
    # swap two partners of a contact if out of order
    if leftChr > rightChr or (leftChr == rightChr and inputLineData[1] > inputLineData[4]):
        newContact = [inputLineData[4], inputLineData[5], inputLineData[1], inputLineData[2], 1]
        temp = leftChr
        leftChr = rightChr
        rightChr = temp
    else:
        newContact = [inputLineData[1], inputLineData[2], inputLineData[4], inputLineData[5], 1]
    # add to contact list and merge
    if not leftChr in contactData:
        contactData[leftChr] = {}
    if not rightChr in contactData[leftChr]:
        contactData[leftChr][rightChr] = []
    contactData[leftChr][rightChr].append(newContact)

# print shifted PDB
inputPdbFile.seek(0)
for inputPdbLine in inputPdbFile:
    if not inputPdbLine.startswith("HETATM"):
        outputFile.write(inputPdbLine)
        continue
    inputPdbLineData = inputPdbLine[6:].strip().split()
    inputPdbChr = int(inputPdbLineData[2].strip("C_"))
    #inputPdbChrDip = inputPdbChr//2
    #inputPdbChrHap = inputPdbChr - inputPdbChrDip*2
    # remove if too few contacts
    numOfContacts = 0
    for chr1 in contactData:
        for chr2 in contactData[chr1]:
            if chr1 != inputPdbChr and chr2 != inputPdbChr:
                continue
            for contact in contactData[chr1][chr2]:
                if (chr1 == inputPdbChr and abs(contact[0]-int(inputPdbLineData[11])) <= halfWindow) or (chr2 == inputPdbChr and abs(contact[2]-int(inputPdbLineData[11])) <= halfWindow):
                    numOfContacts += 1
    if numOfContacts < minContacts:
        continue
    
    outputFile.write(pdb_format % ("HETATM",int(inputPdbLineData[0]),'%-3s' % inputPdbLineData[1]," ",inputPdbLineData[2],inputPdbLineData[3],int(inputPdbLineData[4])," ",float(inputPdbLineData[5])+meanChrPos[inputPdbChr][0]*expandFactor,float(inputPdbLineData[6])+meanChrPos[inputPdbChr][1]*expandFactor,float(inputPdbLineData[7])+meanChrPos[inputPdbChr][2]*expandFactor,0.0,0.0,inputPdbLineData[10],int(inputPdbLineData[11])))

    
    
