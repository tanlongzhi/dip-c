# TODO

## v1.9.9 (planned)

### Merge command portability
- Replace `pigz` with Python `gzip` module (or fall back gracefully)
- Handle macOS BSD `sort` (no `-V` flag) — either detect and use `gsort`, or implement Python-level sort-merge
- Goal: `dip-c merge` works on macOS, not just Linux/HPC

### .hic test files
- Create and check in 1-2 small real `.hic` files (converted from the `.cool` test fixtures)
- Add format-specific smoke tests in `test_hicplot.py`

### Custom color files for `dip-c color`
- Allow users to provide their own color data files (e.g. custom CpG frequencies, epigenetic tracks)
- When the user runs `dip-c color` with a color option:
  1. Check if a file with that name exists in the current working directory
  2. If it does, use the user-provided file instead of the bundled one
  3. If the user explicitly names a file that matches a bundled filename, raise an error to prevent silent clashes between user files and preinstalled data
  4. Only check for the clash when the user actually calls that specific filename. Don't scan the directory preemptively.
- Add a help/list mode: when the user types `dip-c color --list` (or similar), print all available color options (both bundled and any user-provided files in the current directory)
