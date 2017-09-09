import sys
import copy
import pysam
import time

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")

# parameters
maxDistance = 500

# read contacts
contactData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    for i in [2,5]: # phase data
        inputLineData[i] = int(inputLineData[i])
    for i in [1,4]: # position data, for future averaging
        inputLineData[i] = float(inputLineData[i])
    # swap two partners of a contact if out of order
    if inputLineData[0] > inputLineData[3] or (inputLineData[0] == inputLineData[3] and inputLineData[1] > inputLineData[4]):
        newContact = [inputLineData[4], inputLineData[5], inputLineData[1], inputLineData[2], 1]
    else:
        newContact = [inputLineData[1], inputLineData[2], inputLineData[4], inputLineData[5], 1]
    # add to contact list and merge
    if not inputLineData[0] in contactData:
        contactData[inputLineData[0]] = {}
    if not inputLineData[3] in contactData[inputLineData[0]]:
        contactData[inputLineData[0]][inputLineData[3]] = []
    contactData[inputLineData[0]][inputLineData[3]].append(newContact)

# sort and merge contacts within maxDistance
for chr1 in contactData:
    for chr2 in contactData[chr1]:
        print chr1+', '+chr2
        while True:
            merged = False
            contactData[chr1][chr2] = sorted(contactData[chr1][chr2])
            for i in range(len(contactData[chr1][chr2])):
                if i >= len(contactData[chr1][chr2]):
                    break
                for j in range(i+1, len(contactData[chr1][chr2])):
                    if contactData[chr1][chr2][j][0] - contactData[chr1][chr2][i][0] > maxDistance:
                        break
                    elif contactData[chr1][chr2][j][1] == contactData[chr1][chr2][i][1] and abs(contactData[chr1][chr2][j][2] - contactData[chr1][chr2][i][2]) <= maxDistance and contactData[chr1][chr2][j][3] == contactData[chr1][chr2][i][3]:
                        # merge two contacts
                        merged = True
                        mergedContact = [(contactData[chr1][chr2][i][0]*contactData[chr1][chr2][i][4]+contactData[chr1][chr2][j][0]*contactData[chr1][chr2][j][4])/(contactData[chr1][chr2][i][4]+contactData[chr1][chr2][j][4]),contactData[chr1][chr2][i][1],(contactData[chr1][chr2][i][2]*contactData[chr1][chr2][i][4]+contactData[chr1][chr2][j][2]*contactData[chr1][chr2][j][4])/(contactData[chr1][chr2][i][4]+contactData[chr1][chr2][j][4]),contactData[chr1][chr2][i][3],contactData[chr1][chr2][i][4]+contactData[chr1][chr2][j][4]]
                        contactData[chr1][chr2].pop(j)
                        contactData[chr1][chr2].pop(i)
                        contactData[chr1][chr2].insert(i,mergedContact)
                        break
            if merged == False:
                break

# output results
for chr1 in contactData:
    for chr2 in contactData[chr1]:
        for contact in contactData[chr1][chr2]:
            outputFile.write(chr1+'\t'+str(int(round(contact[0])))+'\t'+str(contact[1])+'\t'+chr2+'\t'+str(int(round(contact[2])))+'\t'+str(contact[3])+'\t'+str(contact[4])+'\n')