# Dip-C

**Dip**loid **C**hromatin Conformation Capture — reconstruct 3D diploid genomes from single cells.

## Installation

| Method | When to use | Install time |
|--------|-------------|--------------|
| conda + pip (recommended) | Works everywhere; best for HPC clusters and older Linux | ~1 minute |
| pip only | Modern systems: macOS, Ubuntu 20.04+, RHEL 8+ | ~15 seconds |
| pip only (source build) | Older Linux (e.g. CentOS/RHEL 7, Stanford Sherlock); not recommended — use conda + pip instead | 15–30 minutes |

### Recommended installation (conda + pip)

The easiest way to ensure that the installation will work on your system, regardless of the age of its software, is to create a conda environment and pre-install the compiled dependencies, then install Dip-C with pip:

```bash
conda create -n dipc python=3.11
conda activate dipc
conda install -c conda-forge -c bioconda numpy scipy pysam
pip install run-dipc
```

> **Note:** If you see `Run 'conda init' before 'conda activate'` (common on
> HPC clusters), run this first:
>
> ```bash
> source $(conda info --base)/etc/profile.d/conda.sh
> ```

Verify the installation:

```bash
dip-c --help
```

### Alternative: pip-only install

On modern systems (macOS, Ubuntu 20.04+, RHEL 8+), pip can install everything directly. This is worth trying, but if it fails for any reason we suggest using the above conda + pip installation method:

```bash
pip install run-dipc
```

> **Why does pip fail on older Linux? (optional reading):** If pip failed, use the conda + pip
> method above. Older systems like CentOS/RHEL 7 have an old C library
> (glibc < 2.28), so prebuilt packages for NumPy, SciPy, and pysam are not
> available. Pip falls back to compiling them from source, which takes
> **15–30 minutes** and requires a C++17-capable compiler. If that build also
> fails with `C++ Compiler does not support -std=c++17`, the system's C++
> compiler is too old. You can install a newer one with:
>
> ```bash
> conda install -c conda-forge gcc_linux-64 gxx_linux-64
> ```
>
> But the simplest path is to skip all of this and use the recommended
> conda + pip installation above.

## Upgrading

If you already have Dip-C installed and want to update to the latest version:

```bash
pip install --upgrade run-dipc
```

## Usage

```bash
dip-c <command> [options]
```

Run `dip-c` with no arguments to see all available commands.

## Documentation

Full documentation, workflows, and file format specifications are available on GitHub:

**https://github.com/tanlongzhi/dip-c**

## Citations

Please cite the original Dip-C paper:

> Tan, Longzhi\*; Xing, Dong\*; Chang, Chi-Han; Li, Heng; Xie, X. Sunney "Three-dimensional genome structures of single diploid human cells," *Science* **361**, 924-928 (2018). [DOI:10.1126/science.aat5641](https://doi.org/10.1126/science.aat5641)

## License

MIT
