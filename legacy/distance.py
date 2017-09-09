import sys
import copy
from scipy import interpolate
import numpy as np
import math

# input:
# output:

# read IO locations from arguments
inputPdbFile=open(sys.argv[1],"r")

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
            
# calculate radius of gyration
inputPdbFile.seek(0)
radiusOfGyration = {}
for inputPdbLine in inputPdbFile:
    if not inputPdbLine.startswith("HETATM"):
        continue
    inputPdbLineData = inputPdbLine[6:].strip().split()
    inputPdbChr = int(inputPdbLineData[2].strip("C_"))
    if not inputPdbChr in radiusOfGyration:
        radiusOfGyration[inputPdbChr] = 0.0
    radiusOfGyration[inputPdbChr] += (float(inputPdbLineData[5]) - meanChrPos[inputPdbChr][0])**2 + (float(inputPdbLineData[6]) - meanChrPos[inputPdbChr][1])**2 + (float(inputPdbLineData[7]) - meanChrPos[inputPdbChr][2])**2

# print
for inputPdbChr in radiusOfGyration:
    print str(inputPdbChr)+'\t'+str(math.sqrt(radiusOfGyration[inputPdbChr]/len(inputPdbData[inputPdbChr])))