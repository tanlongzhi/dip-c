import sys
import copy
import rmsd
import numpy as np


# load structure data
inputData = []
numOfStructures = len(sys.argv[1:])
for inputFileName in sys.argv[1:]:
    sys.stderr.write('reading file: '+inputFileName+'\n')
    inputFile = open(inputFileName,"r")
    inputData.append({})
    for inputFileLine in inputFile:
        inputFileLineData = inputFileLine.strip().split()
        inputData[-1][(inputFileLineData[0],int(inputFileLineData[1]))] = [float(inputFileLineData[2]),float(inputFileLineData[3]),float(inputFileLineData[4])]

# find common particles
commonLoci = set(inputData[0])
for inputStructure in inputData[1:]:
    commonLoci = commonLoci.intersection(set(inputStructure))
numOfLoci = len(commonLoci)
commonLoci = list(commonLoci)
commonData = []
for inputStructure in inputData:
    commonData.append([])
    for commonLocus in commonLoci:
        commonData[-1].append(inputStructure[commonLocus])
sys.stderr.write('found common particles: '+str(numOfLoci)+'\n')
        
# subtract centroid
commonData = np.array(commonData)
for commonStructure in commonData:
    commonStructure = np.array(commonStructure)
    commonStructure -= rmsd.centroid(commonStructure)
    
# calculate pairwise deviation
deviations = np.empty((numOfLoci, 0), float)
for i in range(numOfStructures):
    for j in range(i+1, numOfStructures):
        sys.stderr.write('calculating pair: '+str(i)+', '+str(j)+'\n')
        # mirror image if needed
        if rmsd.kabsch_rmsd(commonData[i], commonData[j]) > rmsd.kabsch_rmsd(commonData[i], -1 * commonData[j]):
            commonData[j] *= -1
        # calculate deviation
        deviations = np.c_[deviations, np.linalg.norm(rmsd.kabsch_rotate(commonData[j], commonData[i]) - commonData[i], axis = 1).T]

# summarize rmsd and print
rmsds = np.sqrt((deviations ** 2).mean(axis = 1))
totalrmsd = np.sqrt((rmsds ** 2).mean(axis = 0))
sys.stderr.write('RMS RMSD: '+str(totalrmsd)+'\n')
sys.stderr.write('median RMSD: '+str(np.median(rmsds,axis = 0))+'\n')
for i in range(numOfLoci):
    sys.stdout.write(str(commonLoci[i][0])+'\t'+str(commonLoci[i][1])+'\t'+str(rmsds[i])+'\n')