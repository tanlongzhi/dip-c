import sys
import copy

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")
outputFile=open(sys.argv[2],"w")

# parameters
maxMergeDistance = 1000

inputData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    inputData[inputLineData[0]] = []
    # order the alignments by query start position, and put breakpoints in left-to-right order
    for inputAlignment in inputLineData[1:]:
        inputAlignmentData = inputAlignment.split(",")
        # re-order breakpoints based on read 1/2 and strand
        if (inputAlignmentData[0] == "0" and inputAlignmentData[6] == "-") or (inputAlignmentData[0] == "1" and inputAlignmentData[6] == "+"):
            tempString = inputAlignmentData[4]
            inputAlignmentData[4] = inputAlignmentData[5]
            inputAlignmentData[5] = tempString
        # re-order query positions based on read 1/2
        if inputAlignmentData[0] == "1":
            tempString = inputAlignmentData[1]
            inputAlignmentData[1] = inputAlignmentData[2]
            inputAlignmentData[2] = tempString
                    
        # insertion sort of alignments by query start position
        i = len(inputData[inputLineData[0]])
        for i in range(len(inputData[inputLineData[0]])-1,-1,-1):
            currentData = inputData[inputLineData[0]][i]
            if (currentData[0] == "0" and inputAlignmentData[0] == "0" and int(currentData[1]) < int(inputAlignmentData[1])) or (currentData[0] == "1" and inputAlignmentData[0] == "1" and int(currentData[2]) > int(inputAlignmentData[2])) or (currentData[0] == "0" and inputAlignmentData[0] == "1"):
                i += 1
                break
        inputData[inputLineData[0]].insert(i,inputAlignmentData)
        
    # merge the two adjacent alignments between mates, to utilize phasing information
    for i in range(len(inputData[inputLineData[0]])-1):
        leftData = inputData[inputLineData[0]][i]
        rightData = inputData[inputLineData[0]][i+1]
        if leftData[0] == "0" and rightData[0] == "1" and leftData[3] == rightData[3] and leftData[6] != rightData[6] and abs(int(leftData[5])-int(rightData[4])) <= maxMergeDistance:
            # merge phasing information
            mergedHaplotype = max(int(leftData[7]),int(rightData[7]))
            if leftData[7] == "-2" or rightData[7] == "-2" or (leftData[7] == "1" and rightData[7] == "2") or (leftData[7] == "2" and rightData[7] == "1"):
                mergedHaplotype = -2                   
            mergedData = ["2",leftData[1],rightData[2],leftData[3],leftData[4],rightData[5],leftData[6],str(mergedHaplotype)]
            # replace original with merged alignment
            inputData[inputLineData[0]].pop(i)
            inputData[inputLineData[0]].pop(i)
            inputData[inputLineData[0]].insert(i,mergedData)
            break
    
    # read out contact information
    for i in range(len(inputData[inputLineData[0]])-1):
        leftData = inputData[inputLineData[0]][i]
        rightData = inputData[inputLineData[0]][i+1]
        # flip if chromosomes and positions are out of order
        if leftData[3] > rightData[3] or (leftData[3] == rightData[3] and int(leftData[5]) > int(rightData[4])):
            outputFile.write(rightData[3]+'\t'+rightData[4]+'\t'+rightData[7]+'\t'+leftData[3]+'\t'+leftData[5]+'\t'+leftData[7])
        else:
            outputFile.write(leftData[3]+'\t'+leftData[5]+'\t'+leftData[7]+'\t'+rightData[3]+'\t'+rightData[4]+'\t'+rightData[7])
        if leftData[0] == "0" and rightData[0] == "1": # paired-end contact
            outputFile.write('\tPE\n')
        else:
            outputFile.write('\tSR\n')

