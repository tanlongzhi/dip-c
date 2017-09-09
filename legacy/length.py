import sys
import copy
import pysam

# input: raw (unsorted) BAM
# output: template length

# read IO locations from arguments
inputFile=sys.argv[1]

# parameters
minQual=20
maxSize=2000

# initialize
sizeCounts = [0] * (maxSize + 1)

# open BAM file
samfile = pysam.AlignmentFile(inputFile, "rb")
for read in samfile.fetch(until_eof=True):
    # exclude bad reads
    if read.is_unmapped or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.is_secondary:
        continue
    if read.has_tag("SA"):
        continue
    # exclude bad alignments
    if read.mapping_quality < minQual:
        continue
        
    # calculate insert size
    if not read.is_paired:
        # for single-end
        insertSize = read.reference_end - read.reference_start
    elif read.is_proper_pair and read.template_length > 0:
        # for paired-end
        insertSize = read.template_length
    else:
        continue
        
    # count
    if insertSize <= maxSize:
        sizeCounts[insertSize] += 1

# print
for i in range(maxSize + 1):
    print str(i)+'\t'+str(sizeCounts[i])
        