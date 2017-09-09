import sys
import copy
import math
import numpy as np

# input:
# output:

# parameters
minSeparation=1e4
interSeparation=1e9
minCount=3
resolution=1e6
numOfQuantiles=20
percentileValues=[]
for i in range(numOfQuantiles):
    percentileValues.append(float(i)/(numOfQuantiles-1)*100)

# read IO locations from arguments
inputContactFile=open(sys.argv[1],"r")

# read contacts to calculate new scores
contactScoreData = {}
for inputContactLine in inputContactFile:
    inputContactLineData = inputContactLine.strip().split()
    locus1 = (int(inputContactLineData[0]), int(round(int(inputContactLineData[1])/resolution)*resolution))
    locus2 = (int(inputContactLineData[3]), int(round(int(inputContactLineData[4])/resolution)*resolution))
    inputSeparation = abs(int(inputContactLineData[1]) - int(inputContactLineData[4]))
    if locus1[0] == locus2[0] and inputSeparation < minSeparation:
        continue # skip very close intra
    elif locus1[0] != locus2[0]:
        inputSeparation = interSeparation
    scoreValue = math.log10(inputSeparation)
    if locus1 not in contactScoreData:
        contactScoreData[locus1] = []
    if locus2 not in contactScoreData:
        contactScoreData[locus2] = []
    contactScoreData[locus2].append(scoreValue)
    contactScoreData[locus1].append(scoreValue)
        
# calculate final new scores
for locus in contactScoreData:
    if len(contactScoreData[locus]) < minCount:
        continue
    sys.stdout.write(str(locus[0])+'\t'+str(locus[1])+'\t'+'\t'.join(map(str,np.percentile(contactScoreData[locus],percentileValues)))+'\n')