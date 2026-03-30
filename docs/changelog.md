# Changelog

All notable changes to dip-c (`run-dipc` on PyPI) will be documented in this file.

## [Unreleased]

## [1.9.8] - 2026-03-29

### Added
- `hicplot` subcommand: plot Hi-C contact maps from .hic files (absolute, difference, per-chromosome, and all-in-one modes)
- `merge` subcommand: hierarchical sort-merge of .pairs.gz files using GNU sort and pigz
- Full test coverage for hicplot, merge, hicplot_utils, and merge_utils (774 tests total, 100% coverage)

### Changed
- Replaced `hic-straw` dependency with `hictkpy>=1.4`
- Minimum Python version bumped from 3.9 to 3.10
- Fixed `_version.py` to use correct PyPI package name (`run-dipc`)

## [1.9.6] - 2026-03-08

### Added
- 100% test coverage (574 tests, 0 skips)
- Seg command tests with real and synthetic BAM data
- macOS + Python 3.13 to CI test matrix
- pysam as a required dependency
- Robustness guards across all commands

### Fixed
- `seg -Q` option: was documented in help text but never parsed
- `mgcolor`: error on empty input files instead of silent wrong output
- Dead code removed from `classes.py` (duplicate method definitions, broken methods)

### Changed
- Removed unused files from pre-PyPI build system

### Known issues
- Coverage instrumentation (`--cov`) is disabled on Python 3.13 due to a known incompatibility between pysam's C extension, coverage.py's trace hook, and Python 3.13's garbage collector that causes a segfault during test cleanup. All 574 tests still run and pass on Python 3.13; only coverage measurement is affected.

## [1.9.5] - 2026-03-07

### Fixed
- `vis.py`: `appendAttribute` -> `append_attribute` for mmcif-pdbx 2.x compatibility

## [1.9.4] - 2026-03-07

### Added
- Bundled data file listing in `color` help text
- Install troubleshooting section in PyPI readme

## [1.9.3] - 2026-03-06

### Added
- GitHub Actions CI workflow with test matrix (Python 3.9-3.13, Ubuntu + macOS)
- GitHub Actions publish workflow for PyPI
- Test suite with 3-tier test data (chr22, chr21+22, full genome)
- Bundled data file resolution and tab completion
- CITATION.cff and py.typed marker
- Bioconda recipe template
- Manual trigger for full genome CI tests

### Fixed
- numpy 2.x compatibility: `con3` index and `G3dData.to_np_arrays()` array shape
- Relaxed dependency upper bounds for numpy/scipy
- Python 2->3 bugs: dict mutation during iteration and division by zero
- mmcif-pdbx conda recipe package name

## [1.9.2] - 2026-03-06

### Changed
- Moved legacy scripts to separate Tan_Lab_Publications_Code repository

## [1.9.1] - 2026-03-06

### Added
- Concise PyPI readme

## [1.9.0] - 2026-03-05

### Added
- PyPI packaging with src layout, pyproject.toml, and CLI entry point (`dip-c`)
- MIT LICENSE file

### Changed
- Package renamed to `run-dipc` on PyPI

## [1.0.0] - 2026-03-02

First formally versioned release, marking the Python 3 port of the original dip-c codebase.

### Added
- Python 3 support (ported from Python 2)
- PyMOL coloring scripts and legend PDFs

### Fixed
- `locus_string` parsing in `g3d_particle_to_atom_data`

### History
The project has ~250 commits dating back to 2017, covering:
- **2018**: Core algorithm (seg, con, clean, impute, vis, color), hickit integration
- **2019**: Bug fixes, new color options (`-D`, `-R`)
- **2020**: `ard` normalize counts feature
- **2022**: Community contributions (haplotype separation score, `--c-hom`, `--s-exc` color options via PRs #48, #50, #51)
- **2023**: Cerebellum data, multiome analysis scripts (PR #59)
- **2024**: hg38 color files (PR #66)
