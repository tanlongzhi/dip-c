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
maxDistance = 50

# read contacts
contactData = {}
counter = 0
t = time.time()
for inputLine in inputFile:
    counter += 1
    if counter % 10000 == 0:
        print str(counter)+': '+str(time.time()-t)
    inputLineData = inputLine.strip().split("\t")
    for i in [1,2,4,5]:
        inputLineData[i] = int(inputLineData[i])
    # swap two partners of a contact if out of order
    if inputLineData[0] > inputLineData[3] or (inputLineData[0] == inputLineData[3] and inputLineData[1] > inputLineData[4]):
        inputLineData = [inputLineData[3], inputLineData[4], inputLineData[5], inputLineData[0], inputLineData[1], inputLineData[2], inputLineData[6]]
    
    # add to contact list and merge
    if not inputLineData[0] in contactData:
        contactData[inputLineData[0]] = {}
    if not inputLineData[3] in contactData[inputLineData[0]]:
        contactData[inputLineData[0]][inputLineData[3]] = []
    # each contact cluster: min pos1, max pos1, min pos2, max pos2, num of reads, mean pos1, mean pos2, phase1, phase2
    newCluster=[inputLineData[1],inputLineData[1],inputLineData[4],inputLineData[4],1,float(inputLineData[1]),float(inputLineData[4]),inputLineData[2],inputLineData[5]]
    overlapIndices=[]
    for clusterIndex in range(len(contactData[inputLineData[0]][inputLineData[3]])):
        cluster=contactData[inputLineData[0]][inputLineData[3]][clusterIndex]
        if inputLineData[1]<=cluster[1]+maxDistance and inputLineData[1]>=cluster[0]-maxDistance and inputLineData[4]<=cluster[3]+maxDistance and inputLineData[4]>=cluster[2]-maxDistance:
            newHaplo1 = max(newCluster[7],cluster[7])
            if newCluster[7] == -2 or cluster[7] == -2 or (newCluster[7] == 1 and cluster[7] == 2) or (newCluster[7] == 2 and cluster[7] == 1):
                newHaplo1 = -2                   
            newHaplo2 = max(newCluster[8],cluster[8])
            if newCluster[8] == -2 or cluster[8] == -2 or (newCluster[8] == 1 and cluster[8] == 2) or (newCluster[8] == 2 and cluster[8] == 1):
                newHaplo2 = -2
            newCluster=[min(newCluster[0],cluster[0]),max(newCluster[1],cluster[1]),min(newCluster[2],cluster[2]),max(newCluster[3],cluster[3]),newCluster[4]+cluster[4],(newCluster[5]*newCluster[4]+cluster[5]*cluster[4])/(newCluster[4]+cluster[4]),(newCluster[6]*newCluster[4]+cluster[6]*cluster[4])/(newCluster[4]+cluster[4]),newHaplo1,newHaplo2]
            overlapIndices.append(clusterIndex)
    for clusterIndex in reversed(overlapIndices):
        del contactData[inputLineData[0]][inputLineData[3]][clusterIndex]
    contactData[inputLineData[0]][inputLineData[3]].append(copy.copy(newCluster))
    
# print
for chr1 in contactData:
    for chr2 in contactData[chr1]:
        for cluster in contactData[chr1][chr2]:
            outputFile.write(chr1+'\t'+str(int(round(cluster[5])))+'\t'+str(cluster[7])+'\t'+chr2+'\t'+str(int(round(cluster[6])))+'\t'+str(cluster[8])+'\t'+str(cluster[4])+'\n')
        