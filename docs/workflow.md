# Workflow

## Overview

A typical Dip-C analysis has four stages:

1. **Read alignment** — align Hi-C reads with BWA
2. **Segment extraction & haplotype imputation** — extract contacts with hickit and/or dip-c
3. **3D genome reconstruction** — model 3D structure with hickit
4. **Visualization & analysis** — color, visualize, and analyze with dip-c

The core Dip-C algorithm and 3D modeling are carried out with [hickit](https://github.com/lh3/hickit), a faster and more careful implementation. Dip-c provides the visualization, coloring, and analysis tools.

## Typical Workflow

Below is a typical workflow combining [hickit](https://github.com/lh3/hickit) and dip-c (with mm10 as an example genome):

### 1. Align reads

```sh
bwa mem -5SP genome.fa R1.fq.gz R2.fq.gz | gzip > aln.sam.gz  # for Nextera
# seqtk mergepe R1.fq.gz R2.fq.gz | pre-meta - | bwa mem -5SP -p genome.fa - | gzip > aln.sam.gz  # for META
```

### 2. Extract segments and impute haplotypes

```sh
# extract segments
hickit.js sam2seg -v snp.txt.gz aln.sam.gz | hickit.js chronly -y - | gzip > contacts.seg.gz  # female
# hickit.js sam2seg -v snp.txt.gz aln.sam.gz | hickit.js chronly - | hickit.js bedflt par.bed - | gzip > contacts.seg.gz  # male

# resolve haplotypes via imputation
hickit -i contacts.seg.gz -o - | bgzip > contacts.pairs.gz
hickit -i contacts.pairs.gz -u -o - | bgzip > impute.pairs.gz
```

### 3. Generate 3D structures

```sh
# 3 replicates for robustness
for rep in $(seq 1 3); do
  hickit -s${rep} -M -i impute.pairs.gz \
    -Sr1m -c1 -r10m -c2 \
    -b4m -b1m -O 1m.${rep}.3dg \
    -b200k -O 200k.${rep}.3dg -D5 \
    -b50k -O 50k.${rep}.3dg -D5 \
    -b20k -O 20k.${rep}.3dg
done
```

### 4. Convert formats and clean

```sh
# convert from hickit to dip-c formats
scripts/hickit_pairs_to_con.sh contacts.pairs.gz
scripts/hickit_impute_pairs_to_con.sh impute.pairs.gz

# rescale and clean 3D structures
for rep in $(seq 1 3); do
  scripts/hickit_3dg_to_3dg_rescale_unit.sh 20k.${rep}.3dg
  dip-c clean3 -c impute.con.gz 20k.${rep}.dip-c.3dg > 20k.${rep}.clean.3dg
done
```

### 5. Align replicates

```sh
dip-c align -o aligned.20k. 20k.[1-3].clean.3dg 2> 20k.align.log > 20k.align.color
```

### 6. Visualize and analyze

```sh
# contact map in Juicebox
scripts/con_imputed_to_juicer_pre_short.sh impute.con.gz
java -Xmx2g -jar juicer_tools.jar pre -n impute.juicer.txt.gz impute.hic mm10.chr.hom.len

# chromatin compartment values
dip-c color2 -b1000000 -H -c mm10.cpg.1m.txt -s contacts.con.gz > cpg_b1m.color2  # contact-based
dip-c color -c mm10.cpg.20k.txt -s3 20k.1.clean.3dg > cpg_s3.color                # 3D-structure-based

# radial positioning
dip-c color -C 20k.1.clean.3dg > C.color

# 3D visualization (mmCIF for PyMOL)
dip-c color -n mm10.chr.txt 20k.1.clean.3dg | dip-c vis -c /dev/stdin 20k.1.clean.3dg > 20k.1.clean.n.cif
```

!!! note "Bundled data files"
    Commands like `color` and `color2` accept bundled data file names directly (e.g., `mm10.cpg.20k.txt` instead of the full path). Run `dip-c color -h` to see all available files, or `dip-c data-path` to find the install directory.

## Dip-C Standalone Workflow (Legacy)

Before [hickit](https://github.com/lh3/hickit), Dip-C used its own imputation and a modified version of [nuc_dynamics](https://github.com/TheLaueLab/nuc_dynamics) for 3D modeling. This workflow is documented in an archived [document](https://github.com/tanlongzhi/dip-c/blob/master/README_old.md).

In the standalone workflow, `dip-c seg` extracts segments directly from BAM files:

```sh
dip-c seg input.bam > output.seg                        # unphased
dip-c seg -v snps.txt input.bam > output.phased.seg     # with SNP phasing
dip-c con output.phased.seg > output.con                # contacts
dip-c con -a output.phased.seg > output.adjacent.con    # adjacent-only contacts
```

## Helper Scripts

The `scripts/` directory contains shell and Python utilities for format conversion:

| Script | Purpose |
|--------|---------|
| `hickit_pairs_to_con.sh` | Convert hickit pairs to CON format |
| `hickit_impute_pairs_to_con.sh` | Convert imputed pairs to CON |
| `hickit_3dg_to_3dg_rescale_unit.sh` | Rescale hickit 3DG to unit scale |
| `con_to_juicer_pre_short.sh` | Convert CON to Juicer pre format |
| `con_imputed_to_juicer_pre_short.sh` | Convert imputed CON to Juicer |
| `cpg.sh` | Calculate CpG frequency (requires bedtools) |
| `cg.sh` | Calculate CG frequency |
| `color_chr_to_hom.sh` | Convert chromosome colors to homolog colors |
