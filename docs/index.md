# Dip-C

![logo](images/logo_small.png)

![Coverage](images/coverage-badge.svg)

**Dip**loid **C**hromatin Conformation Capture (Dip-C) reconstructs 3D diploid genomes from single cells by imputing the two chromosome haplotypes linked by each genomic contact.

An alternative (faster and more careful) implementation of the Dip-C algorithm is included in [hickit](https://github.com/lh3/hickit).

## Quick start

```sh
pip install run-dipc
```

Or clone and run directly:

```sh
git clone https://github.com/tanlongzhi/dip-c
cd dip-c
./dip-c
```

See [Installation](installation.md) for details and troubleshooting.

## Overview

Dip-C provides commands for each step of the single-cell 3D genome reconstruction pipeline:

- **`seg`** — Extract read segments from BAM files
- **`con`** — Generate chromatin contacts from segments
- **`clean`** / **`dedup`** — Filter and deduplicate contacts
- **`impute`** — Impute haplotypes
- **`vis`** — Generate mmCIF files for 3D visualization
- **`color`** — Color 3D structures by genomic features
- **`align`** — Align replicate 3D structures

See the full [Workflow](workflow.md) for a typical analysis pipeline.
