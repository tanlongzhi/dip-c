# Developing Dip-C

## Setup

```bash
# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# With pysam (optional, for seg command)
pip install -e ".[dev,seg]"
```

## Running tests

```bash
# Fast subset (chr22 + chr21_22, ~1 min) — same as CI
pytest tests/ -k "not FullGenome"

# Full suite including full genome (~15 min)
pytest tests/

# Single test file or test
pytest tests/test_commands.py
pytest tests/test_commands.py::TestCon3Command::test_con3_generates_contacts
```

### Test tiers

| Tier | Filter | Data | Time |
|------|--------|------|------|
| chr22 | default | ~200KB | fast |
| chr21+22 | default | ~500KB | ~1 min |
| Full genome | `-k FullGenome` | ~19MB | ~15 min |

Tests requiring pysam are skipped automatically if it's not installed.

## CI (GitHub Actions)

**Automatic:** Every push to `pypi-packaging` or `master` runs the fast test suite (no FullGenome) across Python 3.9–3.12 on Ubuntu + 3.12 on macOS.

**Manual full-genome run:**

```bash
# From CLI
gh workflow run test.yml -f full_genome=true --ref pypi-packaging -R conchoecia/dip-c

# Or from the web (only works once test.yml is on master):
# Actions tab → Tests → Run workflow → check "Run full genome tests"
```

Workflow files:
- `.github/workflows/test.yml` — test matrix
- `.github/workflows/publish.yml` — PyPI publish on GitHub release

## Releasing a new version

1. Update `VERSION` file (currently `1.9.2`)
2. Commit, push, and make sure CI is green
3. Create a GitHub release (tag matching the version, e.g. `v1.9.3`)
4. The publish workflow automatically builds and uploads to PyPI
5. Check the workflow logs for the SHA256 hash — you'll need it for the conda recipe

**Requires:** `PYPI_API_TOKEN` secret in the repo settings.

## Project layout

```
src/dip_c/
├── cli.py                 # Entry point, argument dispatch
├── classes.py             # Core data structures (Con, Leg, G3dData, etc.)
├── completion.py          # Bash completion generator
├── py.typed               # PEP 561 marker
├── dip-c.bash-completion  # Static bash completion loader
├── commands/              # One module per subcommand
│   ├── ard.py, clean.py, color.py, con3.py, dist.py, ...
└── data/
    ├── __init__.py        # resolve_data_file(), list_bundled_data_files()
    └── color/             # Bundled data files (chr lengths, color maps)

tests/
├── test_classes.py        # Unit tests for core classes
├── test_commands.py       # Integration tests (3 tiers)
├── test_cli.py            # CLI interface tests
├── test_completion.py     # Completion script tests
├── test_data.py           # Data resolution tests
└── data/                  # Test data (chr22, chr21+22, full genome)

recipes/
├── dip-c/meta.yaml        # conda-forge recipe
└── bioconda/meta.yaml     # bioconda recipe template
```

## Key files

| File | What it controls |
|------|-----------------|
| `VERSION` | Package version (read by setuptools) |
| `pyproject.toml` | Dependencies, build config, extras |
| `CITATION.cff` | GitHub "Cite this repository" button |
| `recipes/*/meta.yaml` | Conda recipes (SHA256 is PLACEHOLDER until PyPI publish) |

## Adding a new command

1. Create `src/dip_c/commands/yourcommand.py` with a function `yourcommand(argv)`
2. Add it to the dispatch in `src/dip_c/cli.py`
3. Add tests in `tests/test_commands.py`
4. The completion script auto-discovers commands from `cli.py`

## Bundled data files

Commands that need data files (chromosome lengths, color maps) use `resolve_data_file()`:

```python
from dip_c.data import resolve_data_file

f = resolve_data_file("hg19.chr.len")  # Returns open file handle
```

This checks for a local file first, then falls back to the bundled package data.
