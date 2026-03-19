#!/bin/bash
#
# Post-process color data after SBATCH jobs from generate_color_data.sh complete.
# Run from Sherlock after all jobs finish.
#
# Usage:
#   bash post_process.sh
#
set -euo pipefail

WORK=/scratch/users/darrints/color_gen
SCRIPTS="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$SCRIPTS/.." && pwd)"
DATA="$REPO/src/dip_c/data/color"

cd "$WORK"

# ---------------------------------------------------------------------------
# Step 1: Verify all SBATCH outputs exist
# ---------------------------------------------------------------------------

echo "=== Checking SBATCH outputs ==="
EXPECTED_FILES=(
    hg38.cpg.500k.raw.txt
    hg38.cg.20k.raw.txt
    mm10.cpg.500k.raw.txt
    mm39.cpg.20k.raw.txt
    mm39.cpg.100k.raw.txt
    mm39.cpg.500k.raw.txt
    mm39.cpg.1m.raw.txt
    mm39.cg.20k.raw.txt
)
MISSING=0
for f in "${EXPECTED_FILES[@]}"; do
    if [ ! -s "$f" ]; then
        echo "MISSING or empty: $f"
        MISSING=1
    fi
done
if [ "$MISSING" -eq 1 ]; then
    echo "ERROR: Some SBATCH outputs are missing. Check job logs."
    exit 1
fi
echo "All SBATCH outputs present."

# ---------------------------------------------------------------------------
# Step 2: Filter scaffold chromosomes — keep only main chromosomes
# ---------------------------------------------------------------------------

echo "=== Filtering to main chromosomes ==="

# hg38: keep chr1-chr22, chrX, chrY
for f in hg38.cpg.500k.raw.txt hg38.cg.20k.raw.txt; do
    out="${f/.raw.txt/.txt}"
    grep -E '^chr([0-9]+|[XY])\b' "$f" > "$out"
    echo "  $f -> $out ($(wc -l < "$out") lines)"
done

# mm10: keep chr1-chr19, chrX, chrY
for f in mm10.cpg.500k.raw.txt; do
    out="${f/.raw.txt/.txt}"
    grep -E '^chr([0-9]+|[XY])\b' "$f" > "$out"
    echo "  $f -> $out ($(wc -l < "$out") lines)"
done

# mm39: keep chr1-chr19, chrX, chrY
for f in mm39.cpg.20k.raw.txt mm39.cpg.100k.raw.txt mm39.cpg.500k.raw.txt mm39.cpg.1m.raw.txt mm39.cg.20k.raw.txt; do
    out="${f/.raw.txt/.txt}"
    grep -E '^chr([0-9]+|[XY])\b' "$f" > "$out"
    echo "  $f -> $out ($(wc -l < "$out") lines)"
done

# ---------------------------------------------------------------------------
# Step 3: Generate .hom variants
# ---------------------------------------------------------------------------

echo "=== Generating .hom variants ==="

# hg38
bash "$SCRIPTS/color_chr_to_hom.sh" hg38.cpg.500k.txt
bash "$SCRIPTS/color_chr_to_hom.sh" hg38.cg.20k.txt

# mm10
bash "$SCRIPTS/color_chr_to_hom.sh" mm10.cpg.500k.txt

# mm39
bash "$SCRIPTS/color_chr_to_hom.sh" mm39.cpg.20k.txt
bash "$SCRIPTS/color_chr_to_hom.sh" mm39.cpg.100k.txt
bash "$SCRIPTS/color_chr_to_hom.sh" mm39.cpg.500k.txt
bash "$SCRIPTS/color_chr_to_hom.sh" mm39.cpg.1m.txt
bash "$SCRIPTS/color_chr_to_hom.sh" mm39.cg.20k.txt

echo "  Done."

# ---------------------------------------------------------------------------
# Step 4: Generate rank files
# ---------------------------------------------------------------------------

echo "=== Generating rank files ==="
python3 "$SCRIPTS/color_to_rank.py" mm39.cpg.20k.txt > mm39.cpg.20k.rank.txt
echo "  mm39.cpg.20k.rank.txt ($(wc -l < mm39.cpg.20k.rank.txt) lines)"

# ---------------------------------------------------------------------------
# Step 5: Gzip and copy all new files to repo
# ---------------------------------------------------------------------------

echo "=== Gzipping and copying to repo ==="

# mm39 reference files
for f in mm39.chr.txt mm39.chr.len mm39.chr.hom.len; do
    gzip -c "$f" > "$DATA/${f}.gz"
    echo "  -> ${f}.gz"
done
gzip -c mm39.chr.cen > "$DATA/mm39.chr.cen.gz"
echo "  -> mm39.chr.cen.gz"

# hg38 new cpg/cg files + hom variants
for f in hg38.cpg.500k.txt hg38.cpg.500k.hom.txt hg38.cg.20k.txt hg38.cg.20k.hom.txt; do
    gzip -c "$f" > "$DATA/${f}.gz"
    echo "  -> ${f}.gz"
done

# mm10 new cpg files + hom variant
for f in mm10.cpg.500k.txt mm10.cpg.500k.hom.txt; do
    gzip -c "$f" > "$DATA/${f}.gz"
    echo "  -> ${f}.gz"
done

# mm39 cpg/cg files + hom variants + rank
for f in mm39.cpg.20k.txt mm39.cpg.100k.txt mm39.cpg.500k.txt mm39.cpg.1m.txt mm39.cg.20k.txt \
         mm39.cpg.20k.hom.txt mm39.cpg.100k.hom.txt mm39.cpg.500k.hom.txt mm39.cpg.1m.hom.txt mm39.cg.20k.hom.txt \
         mm39.cpg.20k.rank.txt; do
    gzip -c "$f" > "$DATA/${f}.gz"
    echo "  -> ${f}.gz"
done

# ---------------------------------------------------------------------------
# Step 6: Verification
# ---------------------------------------------------------------------------

echo ""
echo "=== Verification ==="
echo "Total files in color dir:"
ls "$DATA" | wc -l

echo ""
echo "New file listing:"
ls -lhS "$DATA"

echo ""
echo "Done! Next steps:"
echo "  1. cd $REPO"
echo "  2. pytest tests/ --cov --cov-fail-under=100"
echo "  3. dip-c color  # verify new files listed"
