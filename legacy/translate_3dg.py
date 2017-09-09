import sys
import copy
from scipy import interpolate
import numpy as np

# input:
# output:

# read IO locations from arguments
inputPdbFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")

# parameters
expandFactor = 3.0

# load structure data
inputPdbData = {}
for inputPdbLine in inputPdbFile:
    inputPdbLineData = inputPdbLine.strip().split()
    inputPdbChr = int(inputPdbLineData[0])
    inputPdbLocus = int(inputPdbLineData[1])
    inputPdbPos = [float(inputPdbLineData[2]),float(inputPdbLineData[3]),float(inputPdbLineData[4])]
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

# print shifted PDB
inputPdbFile.seek(0)
for inputPdbLine in inputPdbFile:
    inputPdbLineData = inputPdbLine.strip().split()
    inputPdbChr=int(inputPdbLineData[0])
    outputFile.write(inputPdbLineData[0]+'\t'+inputPdbLineData[1]+'\t'+str(float(inputPdbLineData[2])+meanChrPos[inputPdbChr][0]*expandFactor)+'\t'+str(float(inputPdbLineData[3])+meanChrPos[inputPdbChr][1]*expandFactor)+'\t'+str(float(inputPdbLineData[4])+meanChrPos[inputPdbChr][2]*expandFactor)+'\n')
    
 