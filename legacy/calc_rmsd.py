import sys
import copy
import rmsd
import numpy as np

# input:
# output:

# read IO locations from arguments
inputFile1=open(sys.argv[1],"r")
inputFile2=open(sys.argv[2],"r")

# load structure data
inputData1 = {}
for inputFile1Line in inputFile1:
    inputFile1LineData = inputFile1Line.strip().split()
    inputData1[(int(inputFile1LineData[0]),int(inputFile1LineData[1]))] = [float(inputFile1LineData[2]),float(inputFile1LineData[3]),float(inputFile1LineData[4])]
inputData2 = {}
for inputFile2Line in inputFile2:
    inputFile2LineData = inputFile2Line.strip().split()
    inputData2[(int(inputFile2LineData[0]),int(inputFile2LineData[1]))] = [float(inputFile2LineData[2]),float(inputFile2LineData[3]),float(inputFile2LineData[4])]

# find common particles
commonData1 = []
commonData2 = []
for commonLocus in set(inputData1).intersection(set(inputData2)):
    commonData1.append(inputData1[commonLocus])
    commonData2.append(inputData2[commonLocus])

# align and calculate RMSD
commonData1Array = np.array(commonData1)
commonData2Array = np.array(commonData2)
commonData1Array -= rmsd.centroid(commonData1Array)
commonData2Array -= rmsd.centroid(commonData2Array)
print min(rmsd.kabsch_rmsd(commonData1Array, commonData2Array), rmsd.kabsch_rmsd(commonData1Array, -1 * commonData2Array))
