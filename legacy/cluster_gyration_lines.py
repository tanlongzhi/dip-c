import sys
import copy
import numpy as np
from scipy.spatial import distance

radiusMatrix = np.loadtxt(sys.argv[1])
lociData = np.loadtxt(sys.argv[2])
numOfLoci = radiusMatrix.shape[0]
sys.stderr.write('read '+str(numOfLoci)+' loci\n')

# initialize
clusterData = [] # list of clusters [from, to, id] (inclusive), negative if merged
for i in range(numOfLoci):
    clusterData.append([i,i,0.0])
finalMatrix = np.zeros([numOfLoci-1, 6], dtype = float)
    
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
    #print bestMerge
    #print (clusterData[bestMerge[0]][0]+clusterData[bestMerge[0]][1])/2.0, (clusterData[bestMerge[1]][0]+clusterData[bestMerge[1]][1])/2.0
    #finalMatrix[clusteringStep, :] = [(clusterData[bestMerge[0]][0]+clusterData[bestMerge[0]][1])/2.0, (clusterData[bestMerge[0]][0]+clusterData[bestMerge[1]][1])/2.0, (clusterData[bestMerge[1]][0]+clusterData[bestMerge[1]][1])/2.0, clusterData[bestMerge[0]][2], bestMerge[2], clusterData[bestMerge[1]][2]]
    finalMatrix[clusteringStep, :] = [lociData[clusterData[bestMerge[0]][0]], (lociData[clusterData[bestMerge[0]][0]]+lociData[clusterData[bestMerge[1]][1]])/2.0, lociData[clusterData[bestMerge[1]][1]], 0, (-lociData[clusterData[bestMerge[0]][0]]+lociData[clusterData[bestMerge[1]][1]])/2.0, 0]
    clusterData[bestMerge[0]][1] = clusterData[bestMerge[1]][1] # update the left one
    clusterData[bestMerge[0]][2] = bestMerge[2]
    clusterData[bestMerge[1]][0] = -1 # remove the other one

np.savetxt(sys.stdout, finalMatrix, delimiter='\t')
