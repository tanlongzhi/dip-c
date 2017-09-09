import sys
import copy
import pysam

# input: raw (unsorted) BAM
# output: read name + a list of segments (is read 2, query start, query end, ref name, ref start, ref end, strand)

# read IO locations from arguments
inputFile=sys.argv[1]
outputFile=open(sys.argv[2],"w")

# parameters
minQual=20
maxEditDistance=3

# open BAM file
samfile = pysam.AlignmentFile(inputFile, "rb")
prevQueryName = ""
prevNumAlignment = 0
prevIsProper = False
lineToWrite = ""
for read in samfile.fetch(until_eof=True):
    # exclude bad reads
    if read.is_unmapped or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.is_secondary:
        continue
    # exclude bad alignments
    if read.mapping_quality < minQual:
        continue
        
    if read.query_name != prevQueryName:
        # write previous read
        if prevQueryName != "" and (prevNumAlignment > 2 or (prevNumAlignment == 2 and not prevIsProper)):
            outputFile.write(lineToWrite+'\n')
        # start a new read
        prevNumAlignment = 0
        lineToWrite = read.query_name

    # analyze read
    if read.is_reverse:
        readStrand = "-"
    else:
        readStrand = "+"
    lineToWrite += '\t'+str(int(read.is_read2))+','+str(read.query_alignment_start)+','+str(read.query_alignment_end)+','+str(read.reference_name)+','+str(read.reference_start)+','+str(read.reference_end)+','+readStrand
    prevNumAlignment += 1
    # analyze all SAs
    if read.has_tag("SA"):
        saTags = read.get_tag("SA")
        for saTag in saTags.split(";")[:-1]:
            saRefName, saRefStart, saStrand, saCigar, saMapQ, saEditDistance = saTag.split(",")
            if int(saMapQ) >= minQual and int(saEditDistance) <= maxEditDistance:
                prevNumAlignment += 1
                # get new mapping
                #saRead = pysam.AlignedSegment()
                read.reference_id = samfile.get_tid(saRefName)
                read.reference_start = int(saRefStart)
                read.cigarstring = saCigar
                if read.is_reverse == (saStrand == "-"):
                    saQueryStart = read.query_alignment_start
                    saQueryEnd = read.query_alignment_end
                else:
                    saQueryStart = read.query_length - read.query_alignment_end
                    saQueryEnd = read.query_length - read.query_alignment_start
                # output
                lineToWrite += '\t'+str(int(read.is_read2))+','+str(saQueryStart)+','+str(saQueryEnd)+','+str(read.reference_name)+','+str(read.reference_start)+','+str(read.reference_end)+','+saStrand
                
    prevQueryName = read.query_name
    prevIsProper = read.is_proper_pair
if prevQueryName != "" and (prevNumAlignment > 2 or (prevNumAlignment == 2 and not prevIsProper)):
    outputFile.write(lineToWrite+'\n') 
