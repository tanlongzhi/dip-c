import sys
import copy
from scipy import interpolate
import math
import numpy as np
import bisect

# input:
# output:

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
inputPdbFile=open(sys.argv[2],"r")
resolution = 20e3

# load structure data
inputPdbData = {}
inputPdbLoci = {}
for inputPdbLine in inputPdbFile:
    inputPdbLineData = inputPdbLine.strip().split()
    inputPdbChr = int(inputPdbLineData[0])
    inputPdbLocus = int(inputPdbLineData[1])
    inputPdbPos = [float(inputPdbLineData[2]),float(inputPdbLineData[3]),float(inputPdbLineData[4])]
    if not inputPdbChr in inputPdbData:
        inputPdbData[inputPdbChr]=[]
        inputPdbLoci[inputPdbChr]=[]
    inputPdbData[inputPdbChr].append((inputPdbLocus,inputPdbPos))
    inputPdbLoci[inputPdbChr].append(inputPdbLocus)
    
# sort
for inputPdbChr in inputPdbLoci:
    inputPdbLoci[inputPdbChr] = sorted(inputPdbLoci[inputPdbChr])

# linear interpolate
inputPdbLinearInterpolate = {}
for inputPdbChr in inputPdbData:
    inputPdbLinearInterpolate[inputPdbChr] = interpolate.interp1d(np.array([x[0] for x in inputPdbData[inputPdbChr]]),np.array([x[1] for x in inputPdbData[inputPdbChr]]),axis=0)

# assign phase to each dedup-ed contact
inputFile.readline() # skip header
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    if inputLineData[0] == 'X':
        inputChr = 23
    elif inputLineData[0] == 'Y':
        inputChr = 24
    else:
        inputChr = int(inputLineData[0])
    inputPos1 = (int(inputLineData[1])+int(inputLineData[2]))/2
    inputPos2 = (int(inputLineData[4])+int(inputLineData[5]))/2
    outputString = str(inputChr)+'\t'+str(inputPos1)+'\t'+str(inputPos2)
    isOk = True
    for hap in range(2):
        # check if the loci is in the structure
        #if inputChr*2+hap not in inputPdbLoci:
            #isOk = False
            #break
        #index1 = bisect.bisect(inputPdbLoci[inputChr*2+hap],inputPos1)
        #index2 = bisect.bisect(inputPdbLoci[inputChr*2+hap],inputPos2)
        #if max(index1,index2) + 1 >= len(inputPdbLoci[inputChr*2+hap]):
            #isOk = False
            #break
        #print inputPdbLoci[inputChr*2+hap][index1+1] - inputPdbLoci[inputChr*2+hap][index1], inputPdbLoci[inputChr*2+hap][index2+1] - inputPdbLoci[inputChr*2+hap][index2]
        #if inputPdbLoci[inputChr*2+hap][index1+1] - inputPdbLoci[inputChr*2+hap][index1] > resolution or inputPdbLoci[inputChr*2+hap][index2+1] - inputPdbLoci[inputChr*2+hap][index2] > resolution:
            #isOk = False
            #break        
        # calculate distance
        try:
            loopDistance = np.linalg.norm(inputPdbLinearInterpolate[inputChr*2+hap](inputPos1) - inputPdbLinearInterpolate[inputChr*2+hap](inputPos2))
        except (KeyError, ValueError) as error:
            isOk = False
            break
        outputString += '\t'+str(loopDistance)
    if isOk:
        sys.stdout.write(outputString+'\n')
