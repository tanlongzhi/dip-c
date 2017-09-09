import sys
import copy
import numpy as np

# load single-cell data
inputData = []
numOfCells = len(sys.argv[2:])
cellNames = []
for inputFileName in sys.argv[2:]:
    sys.stderr.write('reading file: '+inputFileName+'\n')
    cellNames.append(inputFileName.split('/')[1])
    inputFile = open(inputFileName,"r")
    inputData.append({})
    for inputFileLine in inputFile:
        inputData[-1][inputFileLine] = []
        
# read CHi-C data line by line
inputChicFile = open(sys.argv[1],'r')
sys.stdout.write('\t'.join(cellNames)+'\n')
for inputChicLine in inputChicFile:
    presenceData = []
    for i in range(numOfCells):
        if inputChicLine in inputData[i]:
            presenceData.append('1')
        else:
            presenceData.append('0')
    sys.stdout.write('\t'.join(presenceData)+'\n')
