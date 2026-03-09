# Installation

## pip (recommended)

```bash
pip install run-dipc
```

This installs the `dip-c` command-line tool along with all Python dependencies:

- NumPy (>=1.22)
- SciPy (>=1.7)
- pysam (>=0.20)
- rmsd (>=1.5)
- mmcif-pdbx (>=2.0)

Verify the installation:

```bash
dip-c --help
```

!!! note "Old Linux systems (CentOS/RHEL 7)"
    If you see `NumPy requires GCC >= 9.3`, your system's default compiler (GCC 4.8) is too old to build NumPy from source. Either load a newer compiler (`module load gcc`) or install NumPy from conda first (`conda install numpy scipy`).

## From source (git clone)

```bash
git clone https://github.com/tanlongzhi/dip-c
cd dip-c
pip install .
```

Or run directly without installing:

```bash
./dip-c seg input.bam > output.seg
```

The `./dip-c` shim at the repository root automatically adds `src/` to the Python path.

## Development install

```bash
git clone https://github.com/tanlongzhi/dip-c
cd dip-c
pip install -e ".[dev]"
pytest tests/
```

## External tools

Some workflows require tools outside the Python package:

| Tool | Used by | Purpose |
|------|---------|---------|
| [hickit](https://github.com/lh3/hickit) | Typical workflow | Fast haplotype imputation and 3D modeling |
| [BWA](https://github.com/lh3/bwa) | Read alignment | Align Hi-C reads to reference genome |
| [SAMtools](http://www.htslib.org/download/) | Read processing | BAM file manipulation |
| [PyMOL](https://pymol.org/2/) | `vis` output | View mmCIF 3D structures |
| [Juicebox](http://www.aidenlab.org/juicebox/) | Contact visualization | Interactive contact map viewer |
| [bedtools](https://bedtools.readthedocs.io/) | `scripts/cpg.sh` | CpG frequency calculation |

## Supported platforms

- **Python**: 3.9, 3.10, 3.11, 3.12, 3.13
- **OS**: Linux (Ubuntu), macOS (Apple Silicon and Intel)
- Tested via GitHub Actions CI on every commit
