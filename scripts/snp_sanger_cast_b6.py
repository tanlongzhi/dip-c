import sys
import copy
import pysam

# this script generates SNP file for CAST/EiJ (female) x C57BL/6J (male) F1 hybrid

# read IO locations from arguments
inputVcfFile=sys.argv[1]
vcfSample="CAST_EiJ"

# open VCF file
bcf_in = pysam.VariantFile(inputVcfFile)
bcf_in.subset_samples([vcfSample])
for rec in bcf_in.fetch():
    vcfAlleles = rec.alleles
    vcfNucleotides = []
    vcfGenotype = rec.samples[vcfSample]["GT"]
    if rec.samples[vcfSample]["FI"] != 1: # did not pass filter
        continue
    if vcfGenotype[0] is None:
        continue
    if vcfGenotype[0] == 0:
        continue
    vcfNucleotide = vcfAlleles[vcfGenotype[0]]
    refNucleotide = vcfAlleles[0]
    sys.stdout.write("\t".join([rec.chrom, str(rec.pos), refNucleotide, vcfNucleotide]) + "\n")
