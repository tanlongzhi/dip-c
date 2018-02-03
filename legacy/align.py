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
centroidData = []
for commonStructure in commonData:
    commonStructure = np.array(commonStructure)
    centroidPos = rmsd.centroid(commonStructure)
    commonStructure -= centroidPos
    centroidData.append(centroidPos)
    
# calculate pairwise deviation and rotate
deviations = np.empty((numOfLoci, 0), float)
for i in range(numOfStructures):
    for j in range(numOfStructures):
        if j == i:
            continue
        # mirror image if needed
        isMirrored = False
        if rmsd.kabsch_rmsd(commonData[i], commonData[j]) > rmsd.kabsch_rmsd(commonData[i], -1 * commonData[j]):
            commonData[j] *= -1
            isMirrored = True
        # calculate deviation
        rotationMatrix = rmsd.kabsch(commonData[j], commonData[i])
        if j > i:
            deviation = np.linalg.norm(np.dot(commonData[j], rotationMatrix) - commonData[i], axis = 1).T
            deviations = np.c_[deviations, deviation]
            sys.stderr.write('median deviation between '+str(i)+' and '+str(j)+': '+str(np.median(deviation))+'\n')
        # rotate j to align with i
        sys.stderr.write('aligning '+str(i)+' to '+str(j)+'\n')
        alignedFilename = "aligned_" + str(j) + "_to_" + str(i) + ".3dg"
        alignedFile = open(alignedFilename, "wb")
        for inputLocus in inputData[j]:
            #sys.stderr.write("locus: "+str(inputLocus)+"\n")
            #sys.stderr.write("old: "+str(np.array(inputData[j][inputLocus]))+"\n")
            #sys.stderr.write("centroid: "+str(centroidData[j])+"\n")
            alignedPos = np.dot((np.array(inputData[j][inputLocus]) - centroidData[j]) * (-1.0 if isMirrored else 1.0), rotationMatrix) + centroidData[i]
            #sys.stderr.write("new: "+str(alignedPos)+"\n")
            #if inputLocus in commonLoci:
                #sys.stderr.write("reference: "+str(commonData[i][commonLoci.index(inputLocus), :])+"\n")
            alignedFile.write("\t".join([inputLocus[0], str(inputLocus[1]), str(alignedPos[0]), str(alignedPos[1]), str(alignedPos[2])]) + "\n")
        alignedFile.close()

# summarize rmsd and print
rmsds = np.sqrt((deviations ** 2).mean(axis = 1))
totalrmsd = np.sqrt((rmsds ** 2).mean(axis = 0))
sys.stderr.write('RMS RMSD: '+str(totalrmsd)+'\n')
sys.stderr.write('median RMSD: '+str(np.median(rmsds,axis = 0))+'\n')
for i in range(numOfLoci):
    sys.stdout.write(str(commonLoci[i][0])+'\t'+str(commonLoci[i][1])+'\t'+str(rmsds[i])+'\n')