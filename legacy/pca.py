import sys
import copy
import numpy as np
from sklearn.decomposition import PCA
import math

# input:
# output:

# read IO locations from arguments
inputPdbFile=open(sys.argv[1],"r")

# load structure data
inputPdbData = {}
for inputPdbLine in inputPdbFile:
    inputPdbLineData = inputPdbLine.strip().split()
    inputPdbChr = int(inputPdbLineData[0])
    inputPdbLocus = int(inputPdbLineData[1])
    inputPdbPos = [float(inputPdbLineData[2]),float(inputPdbLineData[3]),float(inputPdbLineData[4])]
    if not inputPdbChr in inputPdbData:
        inputPdbData[inputPdbChr]=[]
    inputPdbData[inputPdbChr].append(inputPdbPos)

# conduct PCA for each chromosome
chrData = {}
for inputPdbChr in inputPdbData:
    pca = PCA(n_components=3)
    pca.fit(np.array(inputPdbData[inputPdbChr]))
    chrData[inputPdbChr] = pca.explained_variance_

# print
for inputPdbChr in chrData:
    chr_values = chrData[inputPdbChr]
    print(f"{inputPdbChr}\t{chr_values[0]}\t{chr_values[1]}\t{chr_values[2]}")
