import sys
import copy
import pysam
import bisect

# input: sorted BAM + phased tab-delimited genotypes (chr, pos, pat, mat) + precontact.py output
# output: read name + a list of segments (is read 2, query start, query end, ref name, ref start, ref end, strand, haplotype: 0=left, 1=right, -1=unknown, -2=disagree)


# read IO locations from arguments
inputBamFile=sys.argv[1]
inputGenoFile=open(sys.argv[2],"r")
chrName=sys.argv[3]

# parameters
minQual = 20

# open BAM files
samfile = pysam.AlignmentFile(inputBamFile, "rb")

# initialize data for sorted fragments (each fragment: [start, end], and corresponding phase: phase)
fragmentData = []
phaseData = []

# find phased reads from the genotype file
counter = 0
for inputGenoFileLine in inputGenoFile:
    inputGenoFileLineData = inputGenoFileLine.strip().split("\t")
    inputChr = inputGenoFileLineData[0]
    if inputChr != chrName:
        continue
    inputPos = int(inputGenoFileLineData[1])
    leftNucleotide = inputGenoFileLineData[2]
    rightNucleotide = inputGenoFileLineData[3]
    leftCount = 0
    rightCount = 0
    for pileupcolumn in samfile.pileup(inputChr, inputPos-1, inputPos):
        if pileupcolumn.pos == inputPos - 1: # find the position of the SNP
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    read = pileupread.alignment
                    # skip bad reads
                    if read.is_duplicate or read.is_qcfail or read.is_secondary:
                        continue
                    if read.mapping_quality < minQual:
                        continue
                    
                    # find phase
                    currentNucleotide = pileupread.alignment.query_sequence[pileupread.query_position]
                    if currentNucleotide == leftNucleotide:
                        leftCount += 1
                    elif currentNucleotide == rightNucleotide:
                        rightCount += 1
    sys.stdout.write(inputChr+'\t'+str(inputPos)+'\t'+str(leftCount)+'\t'+str(rightCount)+'\n')
        