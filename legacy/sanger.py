import sys
import copy
import pysam

# read IO locations from arguments
inputVcfFile=sys.argv[1]
sampleList=["CAST_EiJ","129S1_SvImJ"]

# open VCF file
bcf_in = pysam.VariantFile(inputVcfFile)
bcf_in.subset_samples(sampleList)
for rec in bcf_in.fetch():
    vcfAlleles = rec.alleles
    vcfNucleotides = []
    hasNone = False
    for vcfSample in sampleList:
        vcfGenotype = rec.samples[vcfSample]["GT"]
        if vcfGenotype[0] == None:
            hasNone = True
            break
        vcfNucleotides.append(vcfAlleles[vcfGenotype[0]])
    if hasNone or vcfNucleotides[0] == vcfNucleotides[1]: # skips bad or homo SNPs
        continue
    sys.stdout.write(rec.chrom+'\t'+str(rec.pos)+'\t'+'\t'.join(vcfNucleotides)+'\n')