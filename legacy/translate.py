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
expandFactor = 2.0
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

# print shifted PDB
inputPdbFile.seek(0)
for inputPdbLine in inputPdbFile:
    if not inputPdbLine.startswith("HETATM"):
        outputFile.write(inputPdbLine)
        continue
    inputPdbLineData = inputPdbLine[6:].strip().split()
    inputPdbChr = int(inputPdbLineData[2].strip("C_"))
    outputFile.write(pdb_format % ("HETATM",int(inputPdbLineData[0]),'%-3s' % inputPdbLineData[1]," ",inputPdbLineData[2],inputPdbLineData[3],int(inputPdbLineData[4])," ",float(inputPdbLineData[5])+meanChrPos[inputPdbChr][0]*expandFactor,float(inputPdbLineData[6])+meanChrPos[inputPdbChr][1]*expandFactor,float(inputPdbLineData[7])+meanChrPos[inputPdbChr][2]*expandFactor,0.0,0.0,inputPdbLineData[10],int(inputPdbLineData[11])))
    
    #outputFile.write("HETATM\t"+inputPdbLineData[0]+'\t'+inputPdbLineData[1]+'\t'+inputPdbLineData[2]+'\t'+inputPdbLineData[3]+'\t'+inputPdbLineData[4]+'\t'+str(float(inputPdbLineData[5])+meanChrPos[inputPdbChr][0]*expandFactor)+'\t'+str(float(inputPdbLineData[6])+meanChrPos[inputPdbChr][1]*expandFactor)+'\t'+str(float(inputPdbLineData[7])+meanChrPos[inputPdbChr][2]*expandFactor)+'\t'+inputPdbLineData[8]+'\t'+inputPdbLineData[9]+'\n')
    
    