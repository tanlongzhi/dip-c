import sys
import copy
import operator
import math
from pdbx.reader.PdbxReader  import PdbxReader
from pdbx.writer.PdbxWriter  import PdbxWriter
from pdbx.reader.PdbxContainers import *

# input:
# output:

# parameters
phases = ["p", "m"]
chrLengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]

# read IO locations from arguments
inputPdbFile=open(sys.argv[1],"r")
inputFile=open(sys.argv[3],"r")

# write
ofh = open(sys.argv[2], "w")

myDataList=[]
curContainer=DataContainer("myblock")
aCat=DataCategory("atom_site")
aCat.appendAttribute("group_PDB")
aCat.appendAttribute("id")
aCat.appendAttribute("label_asym_id")
aCat.appendAttribute("label_comp_id")
aCat.appendAttribute("label_seq_id")  
aCat.appendAttribute("label_atom_id")  
aCat.appendAttribute("Cartn_x")  
aCat.appendAttribute("Cartn_y")  
aCat.appendAttribute("Cartn_z")
aCat.appendAttribute("B_iso_or_equiv")

sCat=DataCategory("struct_conn")
sCat.appendAttribute("conn_type_id")
sCat.appendAttribute("ptnr1_label_asym_id")
sCat.appendAttribute("ptnr1_label_comp_id")
sCat.appendAttribute("ptnr1_label_seq_id")
sCat.appendAttribute("ptnr1_label_atom_id")
sCat.appendAttribute("ptnr2_label_asym_id")
sCat.appendAttribute("ptnr2_label_comp_id")
sCat.appendAttribute("ptnr2_label_seq_id")
sCat.appendAttribute("ptnr2_label_atom_id")

# determine resolution by finding the most frequent increment
incrementCounts = {}
previousChr = -1
previousPosValue = -1
for inputPdbLine in inputPdbFile:
    inputPdbLineData = inputPdbLine.strip().split()
    inputChr = int(inputPdbLineData[0])
    inputPosValue = int(inputPdbLineData[1])
    if inputChr == previousChr:
        increment = inputPosValue - previousPosValue
        if increment not in incrementCounts:
            incrementCounts[increment] = 0
        incrementCounts[increment] += 1
    previousChr = inputChr
    previousPosValue = inputPosValue
resolution = max(incrementCounts.iterkeys(), key=(lambda key: incrementCounts[key]))
print "resolution: "+str(resolution)
    
# read lines
atomId = 1
previousChr = -1
previousPosValue = -1
previousPosList = []
inputPdbFile.seek(0)
for inputPdbLine in inputPdbFile:
    inputPdbLineData = inputPdbLine.strip().split()
    inputChr = int(inputPdbLineData[0])
    inputDipChr = inputChr // 2
    inputDipHap = inputChr % 2
    inputChrName = str(inputDipChr)+phases[inputDipHap]
    inputPosString = inputPdbLineData[1].rjust(9,'0')
    inputPosValue = int(inputPdbLineData[1])
    inputPosList = [inputPosString[0:3],1,inputPosString[3:6]]
    
    # for using position as B factor
    inputBFactor = float(inputPdbLineData[1])/chrLengths[inputDipChr - 1]
    
    aCat.append(('HETATM',atomId,inputChrName,inputPosList[0],inputPosList[1],inputPosList[2],inputPdbLineData[2],inputPdbLineData[3],inputPdbLineData[4],inputBFactor))
    atomId += 1
    #if inputChr == previousChr and inputPosValue - previousPosValue == resolution:
        #sCat.append(('covale',inputChrName,previousPosList[0],previousPosList[1],previousPosList[2],inputChrName,inputPosList[0],inputPosList[1],inputPosList[2]))
    previousChr = inputChr
    previousPosValue = inputPosValue
    previousPosList = copy.copy(inputPosList)
    
# add contacts
bondData = {}
for inputFileLine in inputFile:
    inputFileLineData = inputFileLine.strip().split()
    for i in [1,2,4,5]:
        inputFileLineData[i] = int(inputFileLineData[i])
    inputChrName1 = inputFileLineData[0]+phases[inputFileLineData[2]]
    inputChrName2 = inputFileLineData[3]+phases[inputFileLineData[5]]
    inputPosValue1 = int(round(inputFileLineData[1]/resolution))*resolution
    inputPosValue2 = int(round(inputFileLineData[4]/resolution))*resolution
    if inputChrName1 == inputChrName2 and inputPosValue1 == inputPosValue2:
        continue
    if (inputChrName1, inputPosValue1, inputChrName2, inputPosValue2) in bondData:
        continue
    bondData[(inputChrName1, inputPosValue1, inputChrName2, inputPosValue2)] = []
    inputPosString1 = str(inputPosValue1).rjust(9,'0')
    inputPosString2 = str(inputPosValue2).rjust(9,'0')
    inputPosList1 = [inputPosString1[0:3],1,inputPosString1[3:6]]
    inputPosList2 = [inputPosString2[0:3],1,inputPosString2[3:6]]
    sCat.append(('covale',inputChrName1,inputPosList1[0],inputPosList1[1],inputPosList1[2],inputChrName2,inputPosList2[0],inputPosList2[1],inputPosList2[2]))

curContainer.append(sCat)
curContainer.append(aCat)
myDataList.append(curContainer)
pdbxW=PdbxWriter(ofh)
pdbxW.write(myDataList)
ofh.close()
