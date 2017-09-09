import sys
import copy
import numpy as np
from scipy.spatial import distance

radiusMatrix = np.loadtxt(sys.argv[1])
numOfLoci = radiusMatrix.shape[0]
sys.stderr.write('read '+str(numOfLoci)+' loci\n')

# initialize
clusterData = [] # list of clusters [from, to, id] (inclusive), negative if merged
for i in range(numOfLoci):
    clusterData.append([i,i,i])
linkageMatrix = np.zeros([numOfLoci-1, 3], dtype = float)
    
# hierachical clustering by finding the lowest value radius[from, to] of adjacent merges
for clusteringStep in range(numOfLoci - 1):
    sys.stderr.write('clustering step '+str(clusteringStep+1)+'\n')
    # find all potential merges
    potentialMerge = []
    for i in range(numOfLoci - 1):
        if clusterData[i][0] < 0:
            continue # skip merged clusters
        for j in range(i+1, numOfLoci):
            if clusterData[j][0] < 0:
                continue # skip merged clusters
            potentialMerge.append([i, j, radiusMatrix[clusterData[i][0], clusterData[j][1]]])
            break
    
    # find the merge with lowest value
    bestMerge = sorted(potentialMerge, key=lambda x:x[2])[0]
    linkageMatrix[clusteringStep, :] = [clusterData[bestMerge[0]][2]+1, clusterData[bestMerge[1]][2]+1, bestMerge[2]]
    clusterData[bestMerge[0]][1] = clusterData[bestMerge[1]][1] # update the left one
    clusterData[bestMerge[0]][2] = numOfLoci + clusteringStep
    clusterData[bestMerge[1]][0] = -1 # remove the other one

np.savetxt('link.txt', linkageMatrix, delimiter='\t')
