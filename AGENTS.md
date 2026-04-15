# Agent guidance for dip-c

These notes help AI coding agents (Claude Code, OpenAI Codex, Cursor, etc.)
work effectively in this repo.

## SLURM job submission

- **Time limit:** at least 12 hours — long-running jobs aren't penalized on
  most research partitions, but a timeout wastes the whole run
- **Memory:** request 2× what you estimate the job needs (e.g. if you think
  8 GB is enough, request 16 GB)
- **SLURM spool caveat:** `dirname "$0"` in SBATCH scripts resolves to the
  SLURM spool directory, not the original script location. Use environment
  variables (e.g. `CPG_FAST_PY`) to pass paths to helper scripts. This was written
  with Stanford Sherlock in mind, so there may be caveats on other clusters.
