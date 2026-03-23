#!/bin/bash
#
# Generate missing color data files for hg38, mm10, and mm39.
# Designed to run on Sherlock. Submit from /scratch/users/darrints/.
#
# Usage:
#   bash generate_color_data.sh
#
# After all SBATCH jobs complete, run:
#   bash post_process.sh
#
set -euo pipefail

WORK=/scratch/users/darrints/color_gen
SCRIPTS="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$SCRIPTS/.." && pwd)"

mkdir -p "$WORK"
cd "$WORK"

# ---------------------------------------------------------------------------
# Step 1: Symlink references and generate missing .fai indexes
# ---------------------------------------------------------------------------

echo "=== Step 1: Setting up reference symlinks and indexes ==="

# GRCh38 (hg38)
if [ ! -f GRCh38.primary_assembly.genome.fa ]; then
    ln -s /home/users/tttt/data/references/GRCh38/raw_data/GRCh38.primary_assembly.genome.fa .
fi
if [ ! -f GRCh38.primary_assembly.genome.fa.fai ]; then
    echo "Indexing GRCh38..."
    samtools faidx GRCh38.primary_assembly.genome.fa
fi

# GRCm38 (mm10) — .fai already exists in source dir, symlink both
if [ ! -f GRCm38.primary_assembly.genome.fa ]; then
    ln -s /home/users/tttt/data/references/GRCm38/raw_data/GRCm38.primary_assembly.genome.fa .
    ln -s /home/users/tttt/data/references/GRCm38/raw_data/GRCm38.primary_assembly.genome.fa.fai .
fi

# GRCm39 (mm39)
if [ ! -f GRCm39.primary_assembly.genome.fa ]; then
    ln -s /home/users/darrints/data/references/GRCm39/GRCm39.primary_assembly.genome.fa .
fi
if [ ! -f GRCm39.primary_assembly.genome.fa.fai ]; then
    echo "Indexing GRCm39..."
    samtools faidx GRCm39.primary_assembly.genome.fa
fi

# ---------------------------------------------------------------------------
# Step 2: Generate mm39 basic reference files (chr.txt, chr.len, chr.cen)
# ---------------------------------------------------------------------------

echo "=== Step 2: Generating mm39 reference files ==="

# chr.txt — autosomes + X + Y only (chr1-chr19, chrX, chrY for mouse)
awk -F'\t' '$1 ~ /^chr[0-9]+$/ || $1 == "chrX" || $1 == "chrY"' \
    GRCm39.primary_assembly.genome.fa.fai \
    | cut -f1 | sort -V > mm39.chr.txt

# chr.len — tab-delimited chr + length
awk -F'\t' '$1 ~ /^chr[0-9]+$/ || $1 == "chrX" || $1 == "chrY" {print $1"\t"$2}' \
    GRCm39.primary_assembly.genome.fa.fai \
    | sort -V > mm39.chr.len

# chr.cen — mm39 is acrocentric like mm10; centromeres at position 1555000
awk -F'\t' '{print $1"\t"$2"\t1555000"}' mm39.chr.len > mm39.chr.cen

# chr.hom.len
bash "$SCRIPTS/color_chr_len_to_hom.sh" mm39.chr.len

echo "mm39 reference files:"
wc -l mm39.chr.txt mm39.chr.len mm39.chr.cen mm39.chr.hom.len

# ---------------------------------------------------------------------------
# Step 3: Submit SBATCH jobs for CpG/CG frequency computation
# ---------------------------------------------------------------------------

echo "=== Step 3: Submitting SBATCH jobs ==="

# hg38: cpg.500k and cg.20k
echo "Submitting hg38 jobs..."
SBATCH_OPTS="-p tttt -t 2:00:00 --mem-per-cpu=8000"
export CPG_FAST_PY="$SCRIPTS/cpg_fast.py"

sbatch $SBATCH_OPTS -J hg38_cpg500k -o hg38.cpg.500k.raw.txt \
    "$SCRIPTS/cpg_fast.sh" GRCh38.primary_assembly.genome.fa 500000
sbatch $SBATCH_OPTS -J hg38_cg20k -o hg38.cg.20k.raw.txt \
    "$SCRIPTS/cg.sh" GRCh38.primary_assembly.genome.fa 20000

# mm10: cpg.500k
echo "Submitting mm10 jobs..."
sbatch $SBATCH_OPTS -J mm10_cpg500k -o mm10.cpg.500k.raw.txt \
    "$SCRIPTS/cpg_fast.sh" GRCm38.primary_assembly.genome.fa 500000

# mm39: cpg.20k, cpg.100k, cpg.500k, cpg.1m, cg.20k
echo "Submitting mm39 jobs..."
sbatch $SBATCH_OPTS -J mm39_cpg20k -o mm39.cpg.20k.raw.txt \
    "$SCRIPTS/cpg_fast.sh" GRCm39.primary_assembly.genome.fa 20000
sbatch $SBATCH_OPTS -J mm39_cpg100k -o mm39.cpg.100k.raw.txt \
    "$SCRIPTS/cpg_fast.sh" GRCm39.primary_assembly.genome.fa 100000
sbatch $SBATCH_OPTS -J mm39_cpg500k -o mm39.cpg.500k.raw.txt \
    "$SCRIPTS/cpg_fast.sh" GRCm39.primary_assembly.genome.fa 500000
sbatch $SBATCH_OPTS -J mm39_cpg1m -o mm39.cpg.1m.raw.txt \
    "$SCRIPTS/cpg_fast.sh" GRCm39.primary_assembly.genome.fa 1000000
sbatch $SBATCH_OPTS -J mm39_cg20k -o mm39.cg.20k.raw.txt \
    "$SCRIPTS/cg.sh" GRCm39.primary_assembly.genome.fa 20000

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo "When all jobs complete, run: bash $SCRIPTS/post_process.sh"
