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
inputInfoFile=open(sys.argv[3],"r")
inputExpressionFile=open(sys.argv[4],"r")

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

# load gene info data
infoData = {}
for inputInfoLine in inputInfoFile:
    inputInfoFileData = inputInfoLine.strip().split()
    if inputInfoFileData[5] != "protein_coding": # keep only coding genes
        continue
    if inputInfoFileData[1] == 'X':
        inputChr = 23
    elif inputInfoFileData[1] == 'Y':
        inputChr = 24
    else:
        try:
            inputChr = int(inputInfoFileData[1])
        except ValueError:
            continue
    inputPosValue = int(round((int(inputInfoFileData[2])+int(inputInfoFileData[3]))/2/resolution))*resolution
    infoData[inputInfoFileData[0]] = (inputChr, inputPosValue)

# load RNA-Seq data for B factors
bFactorData = {}
for inputExpressionLine in inputExpressionFile:
    inputExpressionLineData = inputExpressionLine.strip().split()
    inputGeneId = inputExpressionLineData[0].split(".")[0]
    if inputGeneId not in infoData:
        continue
    inputLocus = infoData[inputGeneId]
    inputExpression = math.log10(float(inputExpressionLineData[5]) + 1)
    if inputLocus in bFactorData and bFactorData[inputLocus] >= inputExpression:
        continue
    bFactorData[inputLocus] = inputExpression
    
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
    #inputBFactor = float(inputPdbLineData[1])/chrLengths[inputDipChr - 1]
    
    # for using external bedgraph file as B factor
    inputBFactor = 0
    if (inputDipChr, inputPosValue) in bFactorData:
        inputBFactor = bFactorData[(inputDipChr, inputPosValue)]
    
    aCat.append(('HETATM',atomId,inputChrName,inputPosList[0],inputPosList[1],inputPosList[2],inputPdbLineData[2],inputPdbLineData[3],inputPdbLineData[4],inputBFactor))
    atomId += 1
    if inputChr == previousChr and inputPosValue - previousPosValue == resolution:
        sCat.append(('covale',inputChrName,previousPosList[0],previousPosList[1],previousPosList[2],inputChrName,inputPosList[0],inputPosList[1],inputPosList[2]))
    previousChr = inputChr
    previousPosValue = inputPosValue
    previousPosList = copy.copy(inputPosList)

curContainer.append(sCat)
curContainer.append(aCat)
myDataList.append(curContainer)
pdbxW=PdbxWriter(ofh)
pdbxW.write(myDataList)
ofh.close()
