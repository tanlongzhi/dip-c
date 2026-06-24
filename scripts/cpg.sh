#!/bin/bash
#SBATCH -n 1
#SBATCH -t 3-00:00
#SBATCH -p general
#SBATCH --mem-per-cpu=8000
#
# cpg.sh -- VERSION 2 (faster reimplementation of the original cpg.sh)
#
# Computes, per genomic bin, the CpG dinucleotide fraction
#     (# CG) / (# valid dinucleotides)
# where a dinucleotide is "valid" iff both bases are A/C/G/T. Bins are centered
# at multiples of <window_size> (same windows and same output format as v1:
#     <chrom>  <bin_center>  <cpg_fraction>
# one line per bin, bins with no valid dinucleotide are dropped).
#
# Changes from version 1:
#   * Much faster: ONE bedtools getfasta call for all windows (v1 called
#     getfasta once per window in a shell loop), and CG / valid-dinucleotide
#     counts via regex gsub (C engine) instead of a per-character awk loop.
#     The math is identical:
#       validBases = #[ACGT];  runs = #maximal [ACGT]+ runs
#       valid dinucleotides = validBases - runs   (a run of length L has L-1)
#       cg = #CG  (CG never self-overlaps, so non-overlapping gsub is exact)
#   * Counts SOFT-MASKED (lowercase) bases: the sequence is upper-cased
#     (toupper) before counting, so soft-masked repeats ARE included. Version 1
#     was case-sensitive (only matched uppercase ACGT) and therefore silently
#     dropped soft-masked bases, treating them like Ns. Use this version on a
#     soft-masked genome (e.g. UCSC hg38) to count it the same way an unmasked
#     genome (e.g. GENCODE) would.
#
# Usage: cpg.sh <genome.fa> <window_size> [<genome_sizes_for_windows>]
#   <genome.fa> must be indexed (<genome.fa>.fai is required by getfasta).
#   <genome_sizes_for_windows> defaults to <genome.fa>.fai (i.e. all contigs);
#   pass a chrom-sizes file (e.g. <name>.chr.len) to restrict windows to
#   specific chromosomes.
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

# extend end by +1 so the boundary dinucleotide is counted (as in v1), fetch all
# sequences in a single call, compute the fraction per window. toupper() makes
# soft-masked (lowercase) bases count as ordinary A/C/G/T.
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
