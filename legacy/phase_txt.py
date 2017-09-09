import sys
import copy
import pysam

# input: sorted BAM + phased tab-delimited genotypes (chr, pos, pat, mat) + precontact.py output
# output: read name + a list of segments (is read 2, query start, query end, ref name, ref start, ref end, strand, haplotype: 0=left, 1=right, -1=unknown, -2=disagree)

# read IO locations from arguments
inputBamFile=sys.argv[1]
inputGenoFile=open(sys.argv[2],"r")
inputFile=open(sys.argv[3],"r")
outputFile=open(sys.argv[4],"w")

# open BAM files
samfile = pysam.AlignmentFile(inputBamFile, "rb")

# load Meta-C data into memory
inputData = {}
for inputLine in inputFile:
    inputLineData = inputLine.strip().split("\t")
    inputData[inputLineData[0]] = []
    for inputAlignment in inputLineData[1:]:
        inputAlignmentData = inputAlignment.split(",")
        inputAlignmentData.append("-1")
        inputData[inputLineData[0]].append(inputAlignmentData)

# open genotype file
counter = 0
for inputGenoFileLine in inputGenoFile:
    inputGenoFileLineData = inputGenoFileLine.strip().split("\t")
    inputChr = inputGenoFileLineData[0]
    inputPos = int(inputGenoFileLineData[1])
    leftNucleotide = inputGenoFileLineData[2]
    rightNucleotide = inputGenoFileLineData[3]
    for pileupcolumn in samfile.pileup(inputChr, inputPos-1, inputPos):
        if pileupcolumn.pos == inputPos - 1: # find the position of the SNP
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    currentName = pileupread.alignment.query_name
                    currentNucleotide = pileupread.alignment.query_sequence[pileupread.query_position]
                    currentHaplotype = "-1"
                    if currentNucleotide == leftNucleotide:
                        currentHaplotype = "0"
                    elif currentNucleotide == rightNucleotide:
                        currentHaplotype = "1"
                    if currentHaplotype == "-1":
                        continue
                    if not currentName in inputData:
                        continue
                    inputReadData = inputData[currentName]
                    for inputAlignmentData in inputReadData:
                        if inputAlignmentData[3] == inputChr and inputPos - 1 >= int(inputAlignmentData[4]) and inputPos <= int(inputAlignmentData[5]):
                            if inputAlignmentData[-1] == "-1":
                                inputAlignmentData[-1] = currentHaplotype
                            elif inputAlignmentData[-1] != currentHaplotype:
                                inputAlignmentData[-1] = "-2"
    counter += 1
    if counter % 10000 == 0:
        sys.stderr.write("processed "+str(counter)+" SNPs, last one at "+inputChr+" "+str(inputPos)+"\n")
    
# print output
for inputRead in inputData:
    outputFile.write(inputRead)
    for inputAlignmentData in inputData[inputRead]:
        outputFile.write('\t'+",".join(inputAlignmentDataField for inputAlignmentDataField in inputAlignmentData))
    outputFile.write("\n")
