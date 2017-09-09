import sys
import bisect

# input:
# output:

# read IO locations from arguments
inputPdbFile=open(sys.argv[1],"r")
inputFile=open(sys.argv[2],"r")
outputFile=open(sys.argv[3],"w")

# parameters
halfWindow = 5e5
minContacts = 30 # for Meta-C
#minContacts = 10 # for Nagano
#minContacts = 3 # for Meta-C blood shallow
#minContacts = 10 # for Meta-C sperm shallow


# load contact data
contactData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in [2,5]: # phase data
        inputLineData[i] = int(inputLineData[i])
    for i in [1,4]: # position data, for future averaging
        inputLineData[i] = float(inputLineData[i])
    leftChr = int(inputLineData[0])*2 + inputLineData[2]
    rightChr = int(inputLineData[3])*2 + inputLineData[5]
    # add to contact list
    if not leftChr in contactData:
        contactData[leftChr] = []
    if not rightChr in contactData:
        contactData[rightChr] = []
    contactData[leftChr].append(inputLineData[1])
    contactData[rightChr].append(inputLineData[4])

# sort contact data
for contactChr in contactData:
    contactData[contactChr] = sorted(contactData[contactChr])
    
# clean 3dg file
for inputPdbLine in inputPdbFile:
    inputPdbLineData = inputPdbLine.strip().split()
    inputPdbChr = int(inputPdbLineData[0])
    inputPdbPos = int(inputPdbLineData[1])
    # remove if too few contacts
    numOfContacts = bisect.bisect_right(contactData[inputPdbChr][bisect.bisect_left(contactData[inputPdbChr],inputPdbPos - halfWindow):],inputPdbPos + halfWindow)
    #print numOfContacts
    if numOfContacts < minContacts:
        continue
    outputFile.write(inputPdbLine)

    
    
