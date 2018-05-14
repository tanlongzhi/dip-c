# Dip-C
![logo](images/logo_small.png)

## Introduction
**Dip**loid **C**hromatin Conformation Capture (Dip-C) reconstructs 3D diploid genomes from single cells by imputing the two chromosome haplotypes linked by each genomic contact.

## Citation
Tan, L.; Xing, D.; Chang, C.H., Xie, X.S. "3D Genome Structures of Single Diploid Human Cells," *in revision* (2018).

## Requirements
### Basic Requirements
Dip-C was tested on Python v2.7.13 (macOS and CentOS), with the following basic requirements:

* NumPy (tested on v1.12.1)
* SciPy (tested on v0.13.3)

### Additional Requirements
Some Dip-C commands have additional requirements:

* Read preprocessing for META: [LIANTI](https://github.com/lh3/lianti) (patch needed), which requires [seqtk](https://github.com/lh3/seqtk) for paired-end reads
* Read alignment: [BWA](https://github.com/lh3/bwa) (tested on v0.7.15), [SAMtools](http://www.htslib.org/download/) (tested on v1.3), and [Sambamba](https://github.com/biod/sambamba/releases) (tested on v0.6.3)
* `seg`: pysam (tested on v0.11.1)
* 3D reconstruction: [nuc_dynamics](https://github.com/TheLaueLab/nuc_dynamics) (patch needed)
* `vis` and other mmCIF scripts: [PDBx Python Parser](http://mmcif.wwpdb.org/docs/sw-examples/python/html/index.html)
* mmCIF viewing: [PyMol](https://pymol.org/2/)

### Patching LIANTI
For META read preprocessing, [LIANTI](https://github.com/lh3/lianti) needs a patch to replace the LIANTI adapters with the META ones:

1. Download the [LIANTI](https://github.com/lh3/lianti) source code.
2. Replace LIANTI's `trim.c` with Dip-C's `patch/trim.c`.
3. Compile LIANTI.

### Patching nuc_dynamics
For 3D reconstruction, [nuc_dynamics](https://github.com/TheLaueLab/nuc_dynamics) needs a patch to (1) change the backbone energy function, (2) skip the removal of isolated contacts, and (3) output in the 3D Genome (3DG) format instead of the original PDB format (which has a 99,999-atom limit):

1. Download the [nuc_dynamics](https://github.com/TheLaueLab/nuc_dynamics) source code.
2. Replace nuc_dynamics' `nuc_dynamics.py` with Dip-C's `patch/nuc_dynamics.py`.
3. Compile nuc_dynamics.

The above was tested in [May 2017](https://github.com/TheLaueLab/nuc_dynamics/commit/55148e68b919a59c2354d9157131f729495424c9). In Jun 2017, nuc_dynamics also changed its output format, along with other improvements. We have yet to update our patch.

## List of Phased SNPs
To resolve the two haplotypes, a list of phased single-nucleotide polymorphisms (SNPs) is required in a tab-delimited format (chromosome, coordinate, paternal nucleotide, maternal nucleotide).

An example SNP file for the GM12878 cell line (in hg19 coordinates), based on the 2016-1.0 release of [Illumina Platinum Genomes](https://www.illumina.com/platinumgenomes.html), is provided as `snps/NA12878.txt.gz`. The first few lines are shown below:

```
1	534324	G	T
1	550515	T	C
1	565419	G	C
1	655642	G	A
1	666543	G	A
1	676118	T	C
1	705882	A	G
1	727529	A	G
1	732772	G	A
1	742825	A	G
```

## Typical Workflow
Below is a typical workflow starting from paired-end META data (FASTQ), part of which is also included in `dip-c.sh`:

```sh
# read preprocessing and alignment
seqtk mergepe R1.fq.gz R2.fq.gz | lianti trim - | bwa mem -Cp hs37m.fa - | samtools view -uS | sambamba sort -o aln.bam /dev/stdin

# identify genomic contacts
dip-c seg -v snp.txt.gz aln.bam | gzip -c > phased.seg.gz
dip-c con phased.seg.gz | gzip -c > raw.con.gz
dip-c dedup raw.con.gz | gzip -c > dedup.con.gz
dip-c reg -p hf dedup.con.gz | gzip -c > reg.con.gz
dip-c clean reg.con.gz | gzip -c > clean.con.gz

# initial imputation of haplotypes
dip-c impute clean.con.gz | gzip -c > impute.con.gz

# further imputation and 3d reconstruction 
con_to_ncc.sh impute.con.gz
nuc_dynamics.sh impute.ncc 0.1
dip-c impute3 -3 impute.3dg clean.con.gz | gzip -c > impute3.round1.con.gz
dip-c clean3 -c impute.con.gz impute.3dg > impute.clean.3dg

con_to_ncc.sh impute3.round1.con.gz
nuc_dynamics.sh impute3.round1.ncc 0.1
dip-c impute3 -3 impute3.round1.3dg clean.con.gz | gzip -c > impute3.round2.con.gz
dip-c clean3 -c impute3.round1.con.gz impute3.round1.3dg > impute3.round1.clean.3dg

con_to_ncc.sh impute3.round2.con.gz
nuc_dynamics.sh impute3.round2.ncc 0.1
dip-c impute3 -3 impute3.round2.3dg clean.con.gz | gzip -c > impute3.round3.con.gz
dip-c clean3 -c impute3.round2.con.gz impute3.round2.3dg > impute3.round2.clean.3dg

con_to_ncc.sh impute3.round3.con.gz
nuc_dynamics.sh impute3.round3.ncc 0.02
dip-c impute3 -3 impute3.round3.3dg clean.con.gz | gzip -c > impute3.round4.con.gz
dip-c clean3 -c impute3.round3.con.gz impute3.round3.3dg > impute3.round3.clean.3dg

con_to_ncc.sh impute3.round4.con.gz
nuc_dynamics.sh impute3.round4.ncc 0.02
dip-c clean3 -c impute3.round4.con.gz impute3.round4.3dg > impute3.round4.clean.3dg

# color by chromosome number and visualize as mmCIF
dip-c color -n color/hg19.chr.txt impute3.round4.clean.3dg | dip-c vis -c /dev/stdin impute3.round4.clean.3dg > impute3.round4.clean.n.cif
```
