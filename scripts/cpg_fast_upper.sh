#!/bin/bash
# Fast reimplementation of dip-c cpg.sh -- UPPERCASE variant.
#
# Same as cpg_fast.sh but ALWAYS uppercases the sequence before counting, so
# soft-masked (lowercase) bases are treated as ordinary A/C/G/T and ARE counted.
# This processes ANY genome (e.g. hg38) the same way the mm10 CpG tracks were
# built -- repeats included -- rather than cpg.sh's default case-sensitive
# behaviour (which silently drops soft-masked bases, like N's).
#
# Use this when you want hg38 (or any soft-masked FASTA) handled like mm10.
# (Equivalent to: cpg_fast.sh <fa> <W> [sizes] upper)
#
# Computes per-bin CpG dinucleotide fraction = (#CG) / (#valid dinucleotides),
# where a dinucleotide is "valid" iff both bases are A/C/G/T (case-insensitive
# here, because we uppercase first). Bins are centered at multiples of
# window_size, width = window_size (matches makewindows -w W -s W/2, keep
# end%W==W/2).
#
# Speedups vs original:
#   * ONE bedtools getfasta call for all windows (was 1 call per window).
#   * CG / valid-dinucleotide counts via regex gsub (C engine) instead of a
#     per-character awk loop.  Identical math:
#       validBases = #[ACGT];  runs = #maximal [ACGT]+ runs
#       valid dinucleotides = validBases - runs   (each run of length L gives L-1)
#       cg = #CG  (CG never self-overlaps, so non-overlapping gsub is exact)
#
# Usage: cpg_fast_upper.sh <genome.fa> <window_size> [<genome_sizes_for_windows>]
#   genome_sizes default: <genome.fa>.fai
set -euo pipefail

genome_file="$1"
window_size="$2"
win_g="${3:-${genome_file}.fai}"

half=$(( window_size / 2 ))
tmpdir=$(mktemp -d --tmpdir="${TMPDIR:-.}")
trap 'rm -rf "$tmpdir"' EXIT

# centered, non-overlapping windows at multiples of window_size
bedtools makewindows -g "$win_g" -w "$window_size" -s "$half" \
  | awk -v w="$window_size" -v h="$half" '($3 % w == h)' > "$tmpdir/windows.bed"

# extend end by +1 so the boundary dinucleotide is counted (as in original),
# fetch all sequences in a single call, compute fraction per window.
# toupper() makes soft-masked (lowercase) bases count as ordinary A/C/G/T.
awk -F'\t' 'BEGIN{OFS=FS}{print $1,$2,$3+1}' "$tmpdir/windows.bed" \
  | bedtools getfasta -fi "$genome_file" -bed - -tab \
  | awk -F'\t' '{
        s=toupper($2);
        a=s; vb=gsub(/[ACGT]/,"",a);
        b=s; runs=gsub(/[ACGT]+/,"",b);
        c=s; cg=gsub(/CG/,"",c);
        sum=vb-runs;
        if(sum>0){ print cg/sum; } else { print -1; }
    }' > "$tmpdir/values.txt"

paste "$tmpdir/windows.bed" "$tmpdir/values.txt" \
  | awk -v OFS='\t' '($4>=0){print $1,($2+$3)/2,$4}'
