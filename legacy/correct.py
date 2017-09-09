import sys
import random

# input:
# output:

# read IO locations from arguments
inputFile1=open(sys.argv[1],"r")
inputFile2=open(sys.argv[2],"r")

# parameters
sampleProbability = 0.05

# read two sets of points
inputData1 = []
for inputLine in inputFile1:
    inputLineData = inputLine.strip().split('\t')
    inputChr = int(inputLineData[0][:-1])
    if inputLineData[0][-1] == 'p':
        inputPhase = 0
    elif inputLineData[0][-1] == 'm':
        inputPhase = 1
    else:
        continue
    inputLocus = int(inputLineData[1])
    inputData1.append([inputChr, inputLocus, inputPhase])
    
inputData2 = []
for inputLine in inputFile2:
    inputLineData = inputLine.strip().split('\t')
    inputChr = int(inputLineData[0][:-1])
    if inputLineData[0][-1] == 'p':
        inputPhase = 0
    elif inputLineData[0][-1] == 'm':
        inputPhase = 1
    else:
        continue
    inputLocus = int(inputLineData[1])
    inputData2.append([inputChr, inputLocus, inputPhase])

# random sampling of edges connecting the two sets
numOfEdges = int(len(inputData1)*len(inputData2)*sampleProbability)
sys.stderr.write('adding '+str(numOfEdges)+' contacts\n')
for i in range(numOfEdges):
    locus1 = random.choice(inputData1)
    locus2 = random.choice(inputData2)
    sys.stdout.write(str(locus1[0])+'\t'+str(locus1[1])+'\t'+str(locus1[2])+'\t'+str(locus2[0])+'\t'+str(locus2[1])+'\t'+str(locus2[2])+'\n')