import sys

# input: 
# output: 

# read IO locations from arguments
inputFile=open(sys.argv[1],"r")

# positions (0-based) of PARs in hg19 (parPositions[PAR ID = 1 or 2][chr = X or Y][start or end])
xName = "23"
yName = "24"
parNames = ["PAR1","PAR2"]
parPositions = [[[60000,2699520],[10000,2649520]],[[154931043,155260560],[59034049,59363566]]]

for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    if inputLineData[0] != xName or inputLineData[3] != xName:
        continue
    inPar = False
    for parPosition in parPositions:
        if int(inputLineData[1]) >= parPosition[0][0] and int(inputLineData[1]) <= parPosition[0][1]:
            inPar = True
        if int(inputLineData[4]) >= parPosition[0][0] and int(inputLineData[4]) <= parPosition[0][1]:
            inPar = True
    if inPar:
        continue
    sys.stdout.write(inputLineData[0]+'\t'+inputLineData[1]+'\t1\t'+inputLineData[3]+'\t'+inputLineData[4]+'\t1\n')