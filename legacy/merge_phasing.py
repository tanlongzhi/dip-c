import sys

inputSpermFile=open(sys.argv[1],"r")
inputTrioFile=open(sys.argv[2],"r")

# read and write sperm data
spermLoci = set()
for inputSpermLine in inputSpermFile:
    sys.stdout.write(inputSpermLine)
    inputSpermLocus = inputSpermLine.strip().split()[0:2]
    spermLoci.add(tuple(inputSpermLocus))

# read and write trio data that not in sperm data
for inputTrioLine in inputTrioFile:
    inputTrioLocus = inputTrioLine.strip().split()[0:2]
    if tuple(inputTrioLocus) in spermLoci:
        continue
    sys.stdout.write(inputTrioLine)
