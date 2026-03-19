#!/bin/bash
#SBATCH -J dipc-tests
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -p normal
#SBATCH --mem-per-cpu=8000
#SBATCH -o dipc-tests-%j.out
#
# Run dip-c test suite on Sherlock.
#
# Usage (from repo root):
#   sbatch scripts/run_tests.sh
#
# Results will be in dipc-tests-<jobid>.out
#
set -euo pipefail

module load python/3.9.0 2>/dev/null || true

# Install in editable mode if needed
pip install -e ".[dev]" --quiet 2>/dev/null || pip install -e . --quiet

python3 -m pytest tests/ --cov --cov-fail-under=100 -v
