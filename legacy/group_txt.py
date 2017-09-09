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
expectedOverlap = 9
fuzzyRange = 2
minQual = 20

# open BAM files
samfile = pysam.AlignmentFile(inputBamFile, "rb")

# initialize data for sorted fragments (each fragment: [chr, start, end], and corresponding phase: phase)
fragmentData = []
phaseData = []

def mergePhase (phase1, phase2):
    if phase1 == phase2:
        return phase1
    elif phase1 == -2 or phase2 == -2:
        return -2
    elif phase1 == -1:
        return phase2
    elif phase2 == -1:
        return phase1
    return -2

# function to add reads
def addRead (fragmentData, phaseData, newFragment, newPhase):
    insertionIndex = bisect.bisect_left(fragmentData, newFragment)
    for i in range(insertionIndex, len(fragmentData)):
        if fragmentData[i] == newFragment:
            phaseData[i] = mergePhase(phaseData[i], newPhase)
            return
        else:
            continue
    fragmentData.insert(insertionIndex, newFragment)
    phaseData.insert(insertionIndex, newPhase)
    return

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
    for pileupcolumn in samfile.pileup(inputChr, inputPos-1, inputPos):
        if pileupcolumn.pos == inputPos - 1: # find the position of the SNP
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    read = pileupread.alignment
                    # skip bad reads
                    if read.is_supplementary or read.is_duplicate or read.is_qcfail or read.is_secondary:
                        continue
                    if read.has_tag("SA"):
                        continue
                    if read.mapping_quality < minQual:
                        continue
                    if read.is_paired and not read.is_proper_pair:
                        continue
                    
                    # find phase
                    currentNucleotide = pileupread.alignment.query_sequence[pileupread.query_position]
                    currentHaplotype = -1
                    if currentNucleotide == leftNucleotide:
                        currentHaplotype = 0
                    elif currentNucleotide == rightNucleotide:
                        currentHaplotype = 1
                    
                    # add to fragments
                    addRead(fragmentData, phaseData, [read.reference_name, read.reference_start, read.reference_end], currentHaplotype)
    counter += 1
    if counter % 1000 == 0:
        sys.stderr.write("processed "+str(counter)+" SNPs, last one at "+inputChr+" "+str(inputPos)+"\n")
        
# add regular reads
counter = 0
for read in samfile.fetch(chrName):
    # skip bad reads
    if read.is_supplementary or read.is_duplicate or read.is_qcfail or read.is_secondary:
        continue
    if read.has_tag("SA"):
        continue
    if read.mapping_quality < minQual:
        continue
    if read.is_paired and not read.is_proper_pair:
        continue
    # add to fragments
    addRead(fragmentData, phaseData, [read.reference_name, read.reference_start, read.reference_end], -1)

    counter += 1
    if counter % 10000 == 0:
        sys.stderr.write("processed "+str(counter)+" reads, last one at "+read.reference_name+" "+str(read.reference_start)+"\n")

# pass 1: fuzzy starts and ends
for i in range(len(fragmentData)):
    if phaseData[i] == -3: # deleted fragment
        continue
    for j in range(i+1, len(fragmentData)):
        if phaseData[j] == -3:
            continue
        if fragmentData[j][1] >= fragmentData[i][2]: # no overlaps
            break
        if fragmentData[j][0] != fragmentData[i][0]: # next chromosome
            break
        if fragmentData[j][1] - fragmentData[i][1] <= fuzzyRange or abs(fragmentData[j][2] - fragmentData[i][2]) <= fuzzyRange:
            # merge by removing j and updating i
            phaseData[i] = mergePhase(phaseData[i], phaseData[j])      
            if fragmentData[j][2] > fragmentData[i][2]:
                fragmentData[i][2] = fragmentData[j][2]
            phaseData[j] = -3
    if i % 10000 == 0:
        sys.stderr.write("merged starts/ends for "+str(i)+"/"+str(len(fragmentData))+" fragments, last one at "+fragmentData[i][0]+" "+str(fragmentData[i][1])+"\n")

# pass 2: merge 9 bp overlaps
for i in range(len(fragmentData)):
    if phaseData[i] == -3: # deleted fragment
        continue
    for j in range(i+1, len(fragmentData)):
        if phaseData[j] == -3:
            continue
        if fragmentData[j][1] >= fragmentData[i][2]: # no overlaps
            break
        if fragmentData[j][0] != fragmentData[i][0]: # next chromosome
            break
        if fragmentData[i][2] - fragmentData[j][1] >= expectedOverlap - fuzzyRange and fragmentData[i][2] - fragmentData[j][1] <= expectedOverlap + fuzzyRange:
            # merge by removing j and updating i
            phaseData[i] = mergePhase(phaseData[i], phaseData[j])      
            if fragmentData[j][2] > fragmentData[i][2]:
                fragmentData[i][2] = fragmentData[j][2]
            phaseData[j] = -3
    if i % 10000 == 0:
        sys.stderr.write("merged 9 bp for "+str(i)+"/"+str(len(fragmentData))+" fragments, last one at "+fragmentData[i][0]+" "+str(fragmentData[i][1])+"\n")

for i in range(len(fragmentData)):
    if phaseData[i] == -3:
        continue
    sys.stdout.write(fragmentData[i][0]+'\t'+str(fragmentData[i][1])+'\t'+str(fragmentData[i][2])+'\t'+str(phaseData[i])+'\n')