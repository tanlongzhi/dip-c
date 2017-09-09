import sys

inputFile=open(sys.argv[1],"r")

# parameters
snpsToBin = 100

# read input file line by line
counter = 0
for inputFileLine in inputFile:
    counter += 1
    inputFileLineData = inputFileLine.strip().split()
    currentChr = inputFileLineData[0]
    currentPos = int(inputFileLineData[1])
    if counter % snpsToBin == 1: # starting a bin
        binFirstChr = currentChr
        binFirstPos = currentPos
        leftCounts = 0
        rightCounts = 0
    leftCounts += int(inputFileLineData[2])
    rightCounts += int(inputFileLineData[3])
    if counter % snpsToBin == 0: # finished a bin
        if currentChr != binFirstChr: # if chromosome changed
            continue
        sys.stdout.write(currentChr+'\t'+str((binFirstPos+currentPos)/2.0)+'\t'+str(leftCounts)+'\t'+str(rightCounts)+'\n')