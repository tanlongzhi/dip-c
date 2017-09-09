import sys
import copy
import numpy as np
from scipy.spatial import distance

inputFile=open(sys.argv[1],"r")

# load structure data
inputLoci = []
inputData = []
for inputFileLine in inputFile:
    inputFileLineData = inputFileLine.strip().split()
    inputLoci.append(int(inputFileLineData[1]))
    inputData.append(map(float, inputFileLineData[2:]))
numOfLoci = len(inputLoci)
numOfStructures = len(inputData[0])/3
numOfPairs = numOfLoci*(numOfLoci-1)/2
sys.stderr.write('read '+str(numOfLoci)+' loci from '+str(numOfStructures)+' structures\n')
inputLociNumpy = np.array(inputLoci, dtype=int)
inputDataNumpy = np.array(inputData, dtype=float)

# find resolution and max separation
lociPdist = distance.pdist(np.reshape(inputLociNumpy, [numOfLoci, 1])).astype(int)
resolution = min(lociPdist)
maxSeparation = max(lociPdist)
sys.stderr.write('resolution = '+str(resolution)+' max separation = '+str(maxSeparation)+'\n')

# calculate sum x matrix and sum x^2 matrix
for i in range(numOfStructures):
    sumNumMatrix = np.zeros([numOfLoci, numOfLoci, 3], dtype=float)
    sumMatrix = np.zeros([numOfLoci, numOfLoci, 3], dtype=float)
    sumSqMatrix = np.zeros([numOfLoci, numOfLoci, 3], dtype=float)
    for j in range(numOfLoci):
        sys.stderr.write('processing structure '+str(i+1)+', locus '+str(j+1)+'\n')
        sumNumMatrix[:(j+1), j:, :] += 1
        sumNumMatrix[j:, :(j+1), :] += 1
        sumNumMatrix[j, j, :] -= 1        
        position = inputDataNumpy[j, (3*i):(3*i+3)]
        sumMatrix[:(j+1), j:, :] += position
        sumMatrix[j:, :(j+1), :] += position
        sumMatrix[j, j, :] -= position
        sqPosition = position**2
        sumSqMatrix[:(j+1), j:, :] += sqPosition
        sumSqMatrix[j:, :(j+1), :] += sqPosition
        sumSqMatrix[j, j, :] -= sqPosition
    radiusMatrix = np.sqrt(np.sum(sumSqMatrix / sumNumMatrix - (sumMatrix / sumNumMatrix)**2, axis = 2))
    np.savetxt('radius_'+str(i+1)+'.txt', radiusMatrix, delimiter='\t')
    