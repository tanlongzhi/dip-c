# Dip-C

**Dip**loid **C**hromatin Conformation Capture — reconstruct 3D diploid genomes from single cells.

## Installation

**Recommended (works everywhere, ~1 minute):**

```bash
conda create -n dipc -c conda-forge -c bioconda python=3.11 hictkpy matplotlib numpy scipy pysam
conda activate dipc
pip install run-dipc
```

Verify:

```bash
dip-c version
```

### Why conda first, then pip?

Dip-C depends on `hictkpy` (a C++ library for reading `.hic` files) plus `pysam`, `numpy`, `scipy`, and `matplotlib`. A plain `pip install run-dipc` works on recent laptops where every dependency has a matching prebuilt wheel. On older Linux distributions or HPC clusters where pip can't find compatible wheels, pip falls back to compiling `hictkpy`'s C++ stack (Abseil, Eigen3, fmt, hictk) from source. That build takes 15–30 minutes and frequently fails with cryptic CMake or Conan errors. Bioconda ships `hictkpy` as a precompiled binary that works everywhere, so pre-installing it through conda skips all of that pain.

### Fast path for modern systems

On macOS or recent Linux (Ubuntu 20.04+, RHEL 8+) with a current pip, a direct install usually works:

```bash
pip install run-dipc
```

If pip starts compiling `hictkpy` from source, cancel it (Ctrl-C) and use the recommended conda + pip path above instead.

### HPC gotcha: `Run 'conda init' before 'conda activate'`

Common on Sherlock and similar clusters. Source conda's shell hook first:

```bash
source $(conda info --base)/etc/profile.d/conda.sh
```

## Troubleshooting

### `libpython3.X.so.1.0: cannot open shared object file`

You have a **stale `dip-c` entry-point script** from a previous `pip install --user` that was run while a different Python (often a system module like `python/3.9.0`) was loaded. The shebang in that old script still points to the now-unavailable interpreter, and bash may still have it cached in its hash table. Fix:

```bash
# Check whether a stale dip-c exists outside your conda env
ls -la ~/.local/bin/dip-c 2>/dev/null
head -1 ~/.local/bin/dip-c 2>/dev/null

# Remove it (and any sibling in ~/bin)
rm -f ~/.local/bin/dip-c ~/bin/dip-c

# Clear bash's command cache
hash -r

# Verify the live dip-c is the one in your conda env
type dip-c
dip-c version
```

### `Building wheel for hictkpy ...` hangs or fails with Eigen3/CMake errors

Pip is trying to compile `hictkpy` from source because it can't find a compatible prebuilt wheel. Install `hictkpy` from bioconda instead:

```bash
conda install -c bioconda hictkpy
pip install run-dipc
```

Or start clean with the recommended install command at the top of this page.

### `dip-c: command not found` (after reinstall)

Bash caches command locations per-shell. After uninstalling and reinstalling, the cached path may be stale:

```bash
hash -r
dip-c version
```

Or simply open a new shell.

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
