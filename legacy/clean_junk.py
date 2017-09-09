import sys
import copy
import math
import bisect

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")

# parameters
cleanDistance = 10e6
sqrtCleanDistance = math.sqrt(cleanDistance)
#cleanCounts = 3 # for Nagano
#cleanCounts = 10 # for Meta-C
cleanCounts = 5 # for Meta-C blood
#cleanCounts = 2 # for Meta-C sperm shallow
legDistance = 500
legCounts = 10

# load contact data
contactData = {}
legData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in [1,2,4,5]:
        inputLineData[i] = int(inputLineData[i])
    # add to contact list
    if not inputLineData[0] in contactData:
        contactData[inputLineData[0]] = {}
    if not inputLineData[3] in contactData[inputLineData[0]]:
        contactData[inputLineData[0]][inputLineData[3]] = []
    contactData[inputLineData[0]][inputLineData[3]].append(inputLineData)
    # add to leg list
    if not inputLineData[0] in legData:
        legData[inputLineData[0]] = []
    if not inputLineData[3] in legData:
        legData[inputLineData[3]] = []
    legData[inputLineData[0]].append(inputLineData[1])
    legData[inputLineData[3]].append(inputLineData[4])
        
# sort data
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        contactData[leftChr][rightChr] = sorted(contactData[leftChr][rightChr])
for legChr in legData:
    legData[legChr] = sorted(legData[legChr])


# remove bad legs by counting identical legs
tempContactData = {}
for leftChr in contactData:
    tempContactData[leftChr] = {}
    for rightChr in contactData[leftChr]:
        tempContactData[leftChr][rightChr] = []
        print 'removing '+leftChr+', '+rightChr
        
        for contact in contactData[leftChr][rightChr]:
            numOfLeftLegs = bisect.bisect_right(legData[leftChr][bisect.bisect_left(legData[leftChr],contact[1] - legDistance):],contact[1] + legDistance)
            numOfRightLegs = bisect.bisect_right(legData[rightChr][bisect.bisect_left(legData[rightChr],contact[4] - legDistance):],contact[4] + legDistance)
            if numOfLeftLegs <= legCounts and numOfRightLegs <= legCounts:
                tempContactData[leftChr][rightChr].append(contact)
contactData = tempContactData

# clean data
for leftChr in contactData:
    for rightChr in contactData[leftChr]:
        print 'cleaning '+leftChr+', '+rightChr
        
        numOfContacts = len(contactData[leftChr][rightChr])
        for i in range(numOfContacts):
            nearbyCounts = 0
            for j in range(i+1, numOfContacts):
                if nearbyCounts >= cleanCounts:
                    break
                if contactData[leftChr][rightChr][j][1] - contactData[leftChr][rightChr][i][1] > cleanDistance:
                    break
                if math.sqrt(contactData[leftChr][rightChr][j][1] - contactData[leftChr][rightChr][i][1]) + math.sqrt(abs(contactData[leftChr][rightChr][j][4] - contactData[leftChr][rightChr][i][4])) <= sqrtCleanDistance:
                    nearbyCounts += 1
            for j in range(0, i):
                if nearbyCounts >= cleanCounts:
                    break
                if contactData[leftChr][rightChr][i][1] - contactData[leftChr][rightChr][j][1] > cleanDistance:
                    break
                if math.sqrt(contactData[leftChr][rightChr][i][1] - contactData[leftChr][rightChr][j][1]) + math.sqrt(abs(contactData[leftChr][rightChr][j][4] - contactData[leftChr][rightChr][i][4])) <= sqrtCleanDistance:
                    nearbyCounts += 1
            if nearbyCounts < cleanCounts:
                outputFile.write(contactData[leftChr][rightChr][i][0]+'\t'+str(contactData[leftChr][rightChr][i][1])+'\t'+str(contactData[leftChr][rightChr][i][2])+'\t'+contactData[leftChr][rightChr][i][3]+'\t'+str(contactData[leftChr][rightChr][i][4])+'\t'+str(contactData[leftChr][rightChr][i][5])+'\n')