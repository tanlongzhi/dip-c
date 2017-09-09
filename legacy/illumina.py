import sys
import copy
import pysam

# read IO locations from arguments
inputVcfFile=sys.argv[1]

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
    sys.stdout.write(vcfChr+'\t'+str(rec.pos)+'\t'+leftNucleotide+'\t'+rightNucleotide+'\n')