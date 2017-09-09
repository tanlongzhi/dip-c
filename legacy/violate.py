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
    
# linear interpolate
inputPdbLinearInterpolate = {}
for inputPdbChr in inputPdbData:
    inputPdbLinearInterpolate[inputPdbChr] = interpolate.interp1d(np.array([x[0] for x in inputPdbData[inputPdbChr]]),np.array([x[1] for x in inputPdbData[inputPdbChr]]),axis=0)

# assign phase to each dedup-ed contact
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in range(6):
        inputLineData[i] = int(inputLineData[i])
    isInter = 0
    if inputLineData[0] != inputLineData[3]:
        isInter = 1
    candidateSquaredDistance = []
    outOfBound = False
    noChromosome = False
    for hap1 in range(2) if inputLineData[2] < 0 else [inputLineData[2]]:
        try: 
            pos1 = inputPdbLinearInterpolate[inputLineData[0]*2+hap1](inputLineData[1])
        except KeyError:
            noChromosome = True
            continue
        except ValueError:
            outOfBound = True
            continue
        for hap2 in range(2) if inputLineData[5] < 0 else [inputLineData[5]]:
            try: 
                pos2 = inputPdbLinearInterpolate[inputLineData[3]*2+hap2](inputLineData[4])
            except KeyError:
                noChromosome = True
                continue
            except ValueError:
                outOfBound = True
                continue
            candidateSquaredDistance.append((hap1,hap2,(pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2))
    if noChromosome:
        #print "no chromosome"
        continue
    if outOfBound:
        #print "out of bound"
        continue
    bestCandidateSquaredDistance = min(candidateSquaredDistance, key=lambda x:x[2])
    if bestCandidateSquaredDistance[2] == 0:
        #print "too close" (a problem with dedup that leads to two same legs in a contact)
        continue
    print inputLine.strip()+'\t'+str(isInter)+'\t'+str(len(candidateSquaredDistance))+'\t'+str(math.sqrt(bestCandidateSquaredDistance[2]))+'\t'+str(1/bestCandidateSquaredDistance[2]/sum(1/x[2] for x in candidateSquaredDistance))