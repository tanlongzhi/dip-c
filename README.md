# Dip-C

## Introduction
**Dip**loid **C**hromatin Conformation Capture (Dip-C) reconstructs 3D diploid genomes from single cells by imputing the two chromosome haplotypes linked by each genomic contact.

## Citation
Tan, L.; Xing, D.; Chang, C.H., Xie, X.S. "3D Genome Structures of Single Diploid Human Cells," *submitted* (2018).

## Requirements
Dip-C was tested on Python v2.7.13 (macOS and CentOS), with the following basic requirements:

* NumPy (tested on v1.12.1)
* SciPy (tested on v0.13.3)

Some Dip-C commands have additional requirements:

* Read preprocessing: [LIANTI](https://github.com/lh3/lianti)
* Read alignment: [BWA](https://github.com/lh3/bwa) (tested on v0.7.15)
* `seg`: pysam (tested on v0.11.1)
* `vis` and other mmCIF scripts: [PDBx Python Parser](http://mmcif.wwpdb.org/docs/sw-examples/python/html/index.html)
* mmCIF viewing: [PyMol](https://pymol.org/2/)