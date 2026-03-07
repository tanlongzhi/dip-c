# Dip-C

**Dip**loid **C**hromatin Conformation Capture — reconstruct 3D diploid genomes from single cells.

## Installation

```bash
pip install run-dipc
```

### Troubleshooting: old Linux systems (CentOS/RHEL 7)

If `pip install` fails with **`NumPy requires GCC >= 9.3`**, your system's
default compiler is too old to build NumPy from source. This typically happens
on CentOS/RHEL 7 (end-of-life June 2024), which ships GCC 4.8.

**Option A** — load a newer compiler (HPC clusters usually have one):

```bash
module avail gcc          # see what's available
module load gcc/11.2.0    # load any version >= 9.3
pip install run-dipc
```

**Option B** — install NumPy and SciPy from conda first (provides prebuilt
binaries, no compiler needed):

```bash
conda install numpy scipy
pip install run-dipc
```

You can check your GCC version with `gcc --version` and your OS with
`cat /etc/os-release`. Please contact your system administrator if
this is a problem for you.

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
