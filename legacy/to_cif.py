import sys
import copy
import operator
from pdbx.reader.PdbxReader  import PdbxReader
from pdbx.writer.PdbxWriter  import PdbxWriter
from pdbx.reader.PdbxContainers import *

# input:
# output:

# parameters
phases = ["p", "m"]
# human
chrLengths = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]

# mouse
#chrLengths = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698]

# read IO locations from arguments
inputPdbFile=open(sys.argv[1],"r")

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

cCat=DataCategory("chem_comp")
cCat.appendAttribute("id")
cCat.appendAttribute("name")


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
    if inputChr == previousChr and inputPosValue - previousPosValue == resolution:
        sCat.append(('covale',inputChrName,previousPosList[0],previousPosList[1],previousPosList[2],inputChrName,inputPosList[0],inputPosList[1],inputPosList[2]))
    previousChr = inputChr
    previousPosValue = inputPosValue
    previousPosList = copy.copy(inputPosList)

# chemical component dictionary
for i in range(1000):
    cCat.append((str(i).rjust(3,'0'),str(i).rjust(3,'0')))

curContainer.append(sCat)
curContainer.append(aCat)
#curContainer.append(cCat)
myDataList.append(curContainer)
pdbxW=PdbxWriter(ofh)
pdbxW.write(myDataList)
ofh.close()
