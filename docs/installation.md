# Installation

| Method | When to use | Install time |
|--------|-------------|--------------|
| conda + pip (recommended) | Works everywhere; best for HPC clusters and older Linux | ~1 minute |
| pip only | Modern systems: macOS, Ubuntu 20.04+, RHEL 8+ | ~15 seconds |
| pip only (source build) | Older Linux (e.g. CentOS/RHEL 7, Stanford Sherlock); not recommended — use conda + pip instead | 15–30 minutes |

## Recommended installation (conda + pip)

The easiest way to ensure that the installation will work on your system, regardless of the age of its software, is to create a conda environment and pre-install the compiled dependencies, then install Dip-C with pip:

```bash
conda create -n dipc python=3.11
conda activate dipc
conda install -c conda-forge -c bioconda numpy scipy pysam
pip install run-dipc
```

!!! note "HPC clusters"
    If you see `Run 'conda init' before 'conda activate'`, run this first:

    ```bash
    source $(conda info --base)/etc/profile.d/conda.sh
    ```

Verify the installation:

```bash
dip-c --help
```

## Alternative: pip-only install

On modern systems (macOS, Ubuntu 20.04+, RHEL 8+), pip can install everything directly. This is worth trying, but if it fails for any reason we suggest using the above conda + pip installation method:

```bash
pip install run-dipc
```

!!! note "Why does pip fail on older Linux? (optional reading)"
    If pip failed, use the conda + pip method above. Older systems like
    CentOS/RHEL 7 have an old C library (glibc < 2.28), so prebuilt packages
    for NumPy, SciPy, and pysam are not available. Pip falls back to compiling
    them from source, which takes **15–30 minutes** and requires a C++17-capable
    compiler. If that build also fails with
    `C++ Compiler does not support -std=c++17`, the system's C++ compiler is too
    old. You can install a newer one with:

    ```bash
    conda install -c conda-forge gcc_linux-64 gxx_linux-64
    ```

    But the simplest path is to skip all of this and use the recommended
    conda + pip installation above.

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

Some workflows require external tools that are **not** included in the pip package:

| Tool | Purpose |
|------|---------|
| [hickit](https://github.com/lh3/hickit) | Contact pre-processing, haplotype imputation, and 3D modeling |
| [BWA](https://github.com/lh3/bwa) | Align Hi-C reads to reference genome |
| [SAMtools](http://www.htslib.org/download/) | BAM file manipulation |
| [pre-meta](https://github.com/lh3/pre-pe) (+ [seqtk](https://github.com/lh3/seqtk)) | Read pre-processing for META |
| [PyMOL](https://pymol.org/2/) | View mmCIF 3D structures |
| [Juicebox](http://www.aidenlab.org/juicebox/) | Interactive contact map viewer |
| [bedtools](https://bedtools.readthedocs.io/) | CpG frequency calculation (`scripts/cpg.sh`) |

Dip-C does not call these tools directly. BWA, SAMtools, hickit, and pre-meta are used in upstream steps to generate the input files (`.seg`, `.con`, `.3dg`) that Dip-C operates on. PyMOL and Juicebox are used downstream to view the files that Dip-C produces.

## Supported platforms

- **Python**: 3.9, 3.10, 3.11, 3.12, 3.13
- **OS**: Linux (Ubuntu), macOS (Apple Silicon and Intel)
- Tested via GitHub Actions CI on every commit
