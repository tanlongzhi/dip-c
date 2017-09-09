import sys
import copy
import numpy as np

numOfStructures = len(sys.argv[1:])
outputData = []
for inputFileName in sys.argv[1:]:
    sys.stderr.write('reading file: '+inputFileName+'\n')
    inputData = np.loadtxt(inputFileName)
    if outputData == []:
        outputData = inputData
    else:
        outputData += inputData
outputData /= numOfStructures
np.savetxt(sys.stdout, outputData, delimiter='\t')