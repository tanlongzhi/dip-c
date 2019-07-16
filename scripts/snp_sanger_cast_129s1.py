import sys
import copy
import pysam

# this script generates SNP file for CAST/EiJ (female) x 129S1/SvImJ (male) F1 hybrid

# read IO locations from arguments
inputVcfFile=sys.argv[1]
sampleList=["CAST_EiJ","129S1_SvImJ"] # mother, father

# open VCF file
bcf_in = pysam.VariantFile(inputVcfFile)
bcf_in.subset_samples(sampleList)
for rec in bcf_in.fetch():
    vcfAlleles = rec.alleles
    vcfNucleotides = []
    toSkip = False
    for vcfSample in sampleList:
        vcfGenotype = rec.samples[vcfSample]["GT"]
        if rec.samples[vcfSample]["FI"] != 1: # did not pass filter
            toSkip = True
            break
        if vcfGenotype[0] is None:
            toSkip = True
            break
        vcfNucleotides.append(vcfAlleles[vcfGenotype[0]])
    if toSkip or vcfNucleotides[0] == vcfNucleotides[1]: # skips bad or homo SNPs
        continue
    sys.stdout.write("\t".join([rec.chrom, str(rec.pos), vcfNucleotides[1], vcfNucleotides[0]]) + "\n")
