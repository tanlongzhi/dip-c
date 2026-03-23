#!/bin/bash
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -p tttt
#SBATCH --mem-per-cpu=8000
python3 "$(dirname "$0")/cpg_fast.py" "$@"
