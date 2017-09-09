import sys
import copy
import pysam

# read IO locations from arguments
inputVcfFile=sys.argv[1]
childId=1
fatherId=0
motherId=2

# open VCF file
bcf_in = pysam.VariantFile(inputVcfFile)
vcfSamples = bcf_in.header.samples
childSample = vcfSamples[childId]
fatherSample = vcfSamples[fatherId]
motherSample = vcfSamples[motherId]
for rec in bcf_in.fetch():
    childGenotype = rec.samples[childSample]["GT"]
    fatherGenotype = rec.samples[fatherSample]["GT"]
    motherGenotype = rec.samples[motherSample]["GT"]
    if childGenotype != (0,1):
        continue # only keep het    
    vcfAlleles = rec.alleles
    leftNucleotide = vcfAlleles[0]
    rightNucleotide = vcfAlleles[1]
    if fatherGenotype == (0,1) and motherGenotype == (0,1):
        leftNucleotide = "unknown"
        rightNucleotide = "unknown"
    elif (fatherGenotype == (0,0) and motherGenotype == (0,0)) or (fatherGenotype == (1,1) and motherGenotype == (1,1)):
        leftNucleotide = "conflict"
        rightNucleotide = "conflict"
    elif fatherGenotype == (1,1) or motherGenotype == (0,0):
        leftNucleotide = vcfAlleles[1]
        rightNucleotide = vcfAlleles[0]
    sys.stdout.write(rec.chrom+'\t'+str(rec.pos)+'\t'+leftNucleotide+'\t'+rightNucleotide+'\n')
    
    