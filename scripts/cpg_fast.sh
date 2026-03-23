#!/bin/bash
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -p tttt
#SBATCH --mem-per-cpu=8000
# SLURM copies batch scripts to a spool directory, so dirname($0)
# will not resolve to the original location.  The caller must set
# CPG_FAST_PY to the absolute path of cpg_fast.py.
if [ -z "$CPG_FAST_PY" ]; then
    echo "Error: CPG_FAST_PY is not set" >&2
    exit 2
fi
python3 "$CPG_FAST_PY" "$@"
