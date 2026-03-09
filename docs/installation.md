# Installation

## Requirements

Dip-C requires Python 3.9+ and was tested on macOS and Linux.

### Core dependencies

- NumPy (>=1.22)
- SciPy (>=1.7)
- pysam (>=0.20)
- rmsd (>=1.5)
- mmcif-pdbx (>=2.0)

All are installed automatically via pip.

## pip install

```bash
pip install run-dipc
```

!!! note "Old Linux systems (CentOS/RHEL 7)"
    If you see `NumPy requires GCC >= 9.3`, your system's default compiler (GCC 4.8) is too old to build NumPy from source. Either load a newer compiler (`module load gcc`) or install NumPy from conda first (`conda install numpy scipy`). Please contact your system administrator if this is a problem for you.

## Additional requirements

Some Dip-C commands have additional requirements:

- **Read pre-processing for META** (not required for Nextera or other whole-genome amplification methods): pre-meta from [pre-pe](https://github.com/lh3/pre-pe), which requires [seqtk](https://github.com/lh3/seqtk) for paired-end reads
- **Read alignment**: [BWA](https://github.com/lh3/bwa) (tested on v0.7.15), [SAMtools](http://www.htslib.org/download/) (tested on v1.3)
- **Contact pre-processing & 3D reconstruction**: [hickit](https://github.com/lh3/hickit) (tested on v0.1.1)
- **`vis` and other mmCIF scripts**: [PDBx Python Parser](http://mmcif.wwpdb.org/docs/sw-examples/python/html/index.html)
- **mmCIF viewing**: [PyMol](https://pymol.org/2/)
- **`align`**: [rmsd](https://pypi.org/project/rmsd/)
