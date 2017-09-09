import sys
import copy
from statsmodels.distributions.empirical_distribution import ECDF

# input:
# output:

# parameters
minSeparation=2e6
minCount=3

# read IO locations from arguments
inputBedgraphFile=open(sys.argv[1],"r") # initial scores
inputContactFile=open(sys.argv[2],"r")

# load bedgraph for initial scores
scoreData = {}
for inputBedgraphFileLine in inputBedgraphFile:
    inputBedgraphFileLineData = inputBedgraphFileLine.strip().split()
    if inputBedgraphFileLineData[0] == 'X':
        inputChr = 23
    elif inputBedgraphFileLineData[0] == 'Y':
        inputChr = 24
    else:
        inputChr = int(inputBedgraphFileLineData[0])
    if float(inputBedgraphFileLineData[3]) < 0:
        continue
    scoreData[(inputChr, (int(inputBedgraphFileLineData[1])+int(inputBedgraphFileLineData[2]))/2)] = float(inputBedgraphFileLineData[3])
resolution = int(inputBedgraphFileLineData[2])-int(inputBedgraphFileLineData[1])
sys.stderr.write("resolution: "+str(resolution)+'\n')

# convert to ECDF (so that it's uniform from 0 to 1)
scoreEcdf = ECDF(scoreData.values())
for locus in scoreData:
    scoreData[locus] = scoreEcdf(scoreData[locus])

# read contacts to calculate new scores
contactScoreData = {}
for inputContactLine in inputContactFile:
    inputContactLineData = inputContactLine.strip().split()
    locus1 = (int(inputContactLineData[0]), int(round(int(inputContactLineData[1])/resolution)*resolution))
    locus2 = (int(inputContactLineData[3]), int(round(int(inputContactLineData[4])/resolution)*resolution))
    if locus1[0] == locus2[0] and abs(locus1[1] - locus2[1]) < minSeparation:
        continue # skip very close intra
    if locus1 not in contactScoreData:
        contactScoreData[locus1] = []
    if locus2 not in contactScoreData:
        contactScoreData[locus2] = []
    if locus1 in scoreData:
        contactScoreData[locus2].append(scoreData[locus1])
    if locus2 in scoreData:
        contactScoreData[locus1].append(scoreData[locus2])
        
# calculate final new scores
for locus in contactScoreData:
    initialScore = -1
    if locus in scoreData:
        initialScore = scoreData[locus]
    if len(contactScoreData[locus]) < minCount:
        continue
    sys.stdout.write(str(locus[0])+'\t'+str(locus[1])+'\t'+str(sum(contactScoreData[locus])/len(contactScoreData[locus]))+'\t'+str(len(contactScoreData[locus]))+'\t'+str(initialScore)+'\n')