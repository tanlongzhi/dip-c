import sys
import copy
import pysam
import bisect

# parameters
minNumOfSperm = 2
numOfPhased = 5 # on each side to search
snpFold = 4 # >80% of nearby phased SNPs of each sperm must agree
spermFold = 3 # >75% of voted sperm must agree

# read IO locations from arguments
inputVcfFile = sys.argv[1]
inputDraftFile = open(sys.argv[2],"r")
chrName = sys.argv[3]

numOfSperm=99
childId=numOfSperm+1
fatherId=numOfSperm+0
motherId=numOfSperm+2

# load already phased SNPs
sys.stderr.write('reading phased SNPs\n')
phasedData = {}
for inputDraftLine in inputDraftFile:
    inputDraftLineData = inputDraftLine.strip().split()
    if inputDraftLineData[0][3:] != chrName: # assuming phased file has "chr1" style naming
        continue
    phasedData[int(inputDraftLineData[1])] = [inputDraftLineData[2],inputDraftLineData[3]]

# open VCF file
bcf_in = pysam.VariantFile(inputVcfFile)
vcfSamples = bcf_in.header.samples
childSample = vcfSamples[childId]
fatherSample = vcfSamples[fatherId]
motherSample = vcfSamples[motherId]
spermSamples = []
for spermId in range(numOfSperm):
    spermSamples.append(vcfSamples[spermId])
    
# build draft haplotypes based on phased SNPs
sys.stderr.write('using phased SNPs to determine sperm haplotypes\n')
draftLoci = []
draftHaplotypes = []
for rec in bcf_in.fetch():
    if rec.pos not in phasedData:
        continue # only use phased SNPs
    phasedNucleotides = phasedData[rec.pos]
    phasedHaplotype = []
    for spermSample in spermSamples:
        spermHaplotype = -1 # unknown
        if rec.samples[spermSample]["GT"][0] != None:
            spermNucleotide = rec.alleles[rec.samples[spermSample]["GT"][0]]
            if spermNucleotide == phasedNucleotides[0]:
                spermHaplotype = 0
            elif spermNucleotide == phasedNucleotides[1]:
                spermHaplotype = 1
        phasedHaplotype.append(spermHaplotype)
    draftLoci.append(rec.pos)
    draftHaplotypes.append(phasedHaplotype)
    #print rec.pos, phasedHaplotype
draftLoci, draftHaplotypes = zip(*sorted(zip(draftLoci, draftHaplotypes)))

# start phasing
phasingCounter = 0
sys.stderr.write('using sperm haplotypes to phase unknown SNPs\n')
for rec in bcf_in.fetch():
    childGenotype = rec.samples[childSample]["GT"]
    #fatherGenotype = rec.samples[fatherSample]["GT"]
    #motherGenotype = rec.samples[motherSample]["GT"]
    vcfAlleles = rec.alleles
    spermGenotypes = []
    if childGenotype[0] == childGenotype[1]:
        continue # only use het
    if rec.pos in phasedData:
        continue # only study unknown
    leftGenotype = childGenotype[0]
    rightGenotype = childGenotype[1]
    spermVotes = [0, 0] # left=pat vs. left=mat
    for spermId in range(numOfSperm):
        spermGenotype = rec.samples[spermSamples[spermId]]["GT"][0]
        if spermGenotype == leftGenotype:
            spermSide = 0
        elif spermGenotype == rightGenotype:
            spermSide = 1
        else:
            continue
        siteVotes = [0, 0]
        posId = bisect.bisect_left(draftLoci, rec.pos)

        counter = 0
        for i in range(posId-1, 0, -1):
            if draftHaplotypes[i][spermId] < 0:
                continue
            if draftHaplotypes[i][spermId] == spermSide:
                siteVotes[0] += 1
            else:
                siteVotes[1] += 1
            counter += 1
            if counter >= numOfPhased:
                break

        counter = 0
        for i in range(posId, len(draftLoci), 1):
            if draftHaplotypes[i][spermId] < 0:
                continue
            if draftHaplotypes[i][spermId] == spermSide:
                siteVotes[0] += 1
            else:
                siteVotes[1] += 1
            counter += 1
            if counter >= numOfPhased:
                break

        if siteVotes[0] > snpFold * siteVotes[1]: # > 80% of SNPs agree on a haplotype
            spermVotes[0] += 1
        elif siteVotes[1] > snpFold * siteVotes[0]:
            spermVotes[1] += 1
    if spermVotes[0] + spermVotes[1] < minNumOfSperm: # needs at least two sperms cover
        continue
    if spermVotes[0] > spermFold * spermVotes[1]: # needs > 3x majority vote
        phasingCounter += 1
        sys.stdout.write(chrName+'\t'+str(rec.pos)+'\t'+vcfAlleles[leftGenotype]+'\t'+vcfAlleles[rightGenotype]+'\n')
        if phasingCounter % 1000 == 0:
            sys.stderr.write('phased '+str(phasingCounter)+' SNPs\n')
    elif spermVotes[1] > spermFold * spermVotes[0]: # needs > 3x majority vote
        phasingCounter += 1
        sys.stdout.write(chrName+'\t'+str(rec.pos)+'\t'+vcfAlleles[rightGenotype]+'\t'+vcfAlleles[leftGenotype]+'\n')
        if phasingCounter % 1000 == 0:
            sys.stderr.write('phased '+str(phasingCounter)+' SNPs\n')