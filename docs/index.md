# Dip-C

![logo](images/logo_small.png)

![Coverage](images/coverage-badge.svg)

**Dip**loid **C**hromatin Conformation Capture (Dip-C) reconstructs 3D diploid genomes from single cells by imputing the two chromosome haplotypes linked by each genomic contact.

An alternative (faster and more careful) implementation of the Dip-C algorithm is included in [hickit](https://github.com/lh3/hickit).

## Quick Start

```sh
pip install run-dipc
dip-c seg input.bam > output.seg
```

Or clone and run directly:

```sh
git clone https://github.com/tanlongzhi/dip-c
cd dip-c
./dip-c seg input.bam > output.seg
```

See [Installation](installation.md) for full details.

## Commands

Dip-C provides 30 subcommands covering the full single-cell 3D genome reconstruction pipeline. Run `dip-c <command>` without arguments to see usage for any command.

### Contact processing

| Command | Description |
|---------|-------------|
| `seg` | Extract read segments from a BAM file |
| `con` | Extract contacts from a SEG file |
| `dedup` | Merge duplicates in a CON file |
| `reg` | Exclude genomic regions from a CON file |
| `clean` | Remove artifacts from a CON file |
| `impute` | Impute missing haplotypes in a CON file |
| `mkcon` | Make artificial CON from contact legs |
| `bincon` | Bin a CON file into a contact matrix |

### 3D genome reconstruction

| Command | Description |
|---------|-------------|
| `force` | Generate a 3DG file from a CON file (force-directed layout) |
| `impute3` | Impute missing haplotypes using a 3DG file |
| `clean3` | Remove contact-poor particles from a 3DG file |
| `reg3` | Exclude genomic regions from a 3DG file |
| `align` | Align replicate 3DG files and compute RMSD |
| `exp` | Expand a 3DG file by translating each chromosome |

### Visualization and coloring

| Command | Description |
|---------|-------------|
| `vis` | Convert a 3DG file to mmCIF format for PyMOL |
| `color` | Color 3DG particles by genomic features |
| `color2` | Color particles based on a CON file |
| `mgcolor` | Merge color files |

### Analysis

| Command | Description |
|---------|-------------|
| `info` | Print basic statistics about CON file(s) |
| `dist` | Calculate 3D distances from a 3DG file |
| `pd` | Calculate pairwise distances from a 3DG file |
| `rg` | Calculate the radius of gyration (Rg) matrix |
| `tad` | Find the TAD tree from an Rg matrix |
| `con3` | Generate CON from a 3DG file |
| `pos` | Find 3D positions from a 3DG file |
| `ard` | Contacts around reference points |
| `cv` | Cross-validation of CON files |

### Hi-C and .pairs processing

| Command | Description |
|---------|-------------|
| `hicplot` | Plot Hi-C contact maps from .hic files |
| `merge` | Hierarchical sort-merge of .pairs.gz files (requires GNU sort + pigz) |
| `pairs2con` | Convert hickit .pairs.gz to .con.gz format |

### Utilities

| Command | Description |
|---------|-------------|
| `data-path` | Print path to installed data directory |

## Tab Completion

Enable bash/zsh completion for all commands and bundled data files:

```sh
eval "$(dip-c --completion)"
```

Add this line to your `~/.bashrc` or `~/.zshrc` to make it permanent.

## Bundled Data

Dip-C ships with reference data files for **hg19**, **hg38**, and **mm10** (chromosome lengths, centromere positions, CpG frequencies at various resolutions). These are automatically available to commands like `color` and `bincon`. Run `dip-c data-path` to see where they are installed, or use `dip-c color -h` to list available files.
