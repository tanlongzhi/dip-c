import sys
import copy
import pysam

# input: sorted BAM + phased VCF + precontact.py output
# output: read name + a list of segments (is read 2, query start, query end, ref name, ref start, ref end, strand, haplotype: 0=left, 1=right, -1=unknown, -2=disagree)

# read IO locations from arguments
inputBamFile=sys.argv[1]
inputVcfFile=sys.argv[2]
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
print "finished loading Meta-C data"

# open VCF file
bcf_in = pysam.VariantFile(inputVcfFile)
vcfSample = bcf_in.header.samples[0]
for rec in bcf_in.fetch():
    vcfGenotype = rec.samples[vcfSample]["GT"]
    vcfAlleles = rec.alleles
    # keep only het biallelic SNPs
    if len(vcfGenotype) != 2 or vcfGenotype[0] + vcfGenotype[1] != 1 or len(vcfAlleles[0]) != 1 or len(vcfAlleles[1]) != 1:
        continue
    # find nucleotide for each haplotype
    leftNucleotide = vcfAlleles[vcfGenotype[0]]
    rightNucleotide = vcfAlleles[vcfGenotype[1]]
    vcfChr = rec.chrom[3:]
    for pileupcolumn in samfile.pileup(vcfChr, rec.pos-1, rec.pos): # here we remove "chr" from VCF's chromosome name
        if pileupcolumn.pos == rec.pos - 1: # find the position of the SNP
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
                        if inputAlignmentData[3] == vcfChr and rec.pos - 1 >= int(inputAlignmentData[4]) and rec.pos <= int(inputAlignmentData[5]):
                            if inputAlignmentData[-1] == "-1":
                                inputAlignmentData[-1] = currentHaplotype
                            elif inputAlignmentData[-1] != currentHaplotype:
                                inputAlignmentData[-1] = "-2"
print "finished phasing Meta-C data"

# print output
for inputRead in inputData:
    outputFile.write(inputRead)
    for inputAlignmentData in inputData[inputRead]:
        outputFile.write('\t'+",".join(inputAlignmentDataField for inputAlignmentDataField in inputAlignmentData))
    outputFile.write("\n")
print "finished writing output"