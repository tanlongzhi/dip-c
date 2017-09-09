import sys
import copy

# load structure data
inputData = []
numOfStructures = len(sys.argv[1:])
for inputFileName in sys.argv[1:]:
    sys.stderr.write('reading file: '+inputFileName+'\n')
    inputFile = open(inputFileName,"r")
    inputData.append({})
    for inputFileLine in inputFile:
        inputFileLineData = inputFileLine.strip().split()
        inputData[-1][(int(inputFileLineData[0]),int(inputFileLineData[1]))] = float(inputFileLineData[2])

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

# output merged data
for i in range(len(commonLoci)):
    sys.stdout.write(str(commonLoci[i][0])+'\t'+str(commonLoci[i][1]))
    for j in range(numOfStructures):
        sys.stdout.write('\t'+str(inputData[j][commonLoci[i]]))
    sys.stdout.write('\n')