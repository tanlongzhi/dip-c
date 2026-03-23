#!/usr/bin/env python3
"""Compute CpG dinucleotide frequency in sliding genomic windows.

Functionally equivalent to cpg.sh but runs orders of magnitude faster
by extracting all sequences in a single bedtools getfasta call instead
of forking bedtools once per window.
"""
import re
import subprocess
import sys
import tempfile
import os


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <genome.fa> <window_size>", file=sys.stderr)
        sys.exit(1)

    genome_file = sys.argv[1]
    window_size = int(sys.argv[2])
    half_window_size = window_size // 2

    # Step 1: Generate sliding windows (50% overlap)
    makewindows = subprocess.run(
        ["bedtools", "makewindows", "-g", f"{genome_file}.fai",
         "-w", str(window_size), "-s", str(half_window_size)],
        capture_output=True, text=True, check=True,
    )

    # Filter: keep only windows where end % window_size == half_window_size
    # (matches the awk filter in cpg.sh)
    windows = []
    for line in makewindows.stdout.splitlines():
        chrom, start, end = line.split("\t")
        start, end = int(start), int(end)
        if end % window_size == half_window_size:
            windows.append((chrom, start, end))

    if not windows:
        return

    # Step 2: Write windows to temp BED with end+1 (matches cpg.sh's $3+1)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False,
    ) as bed_fh:
        bed_path = bed_fh.name
        for chrom, start, end in windows:
            bed_fh.write(f"{chrom}\t{start}\t{end + 1}\n")

    try:
        # Step 3: Extract all sequences in one pass (streaming to save memory)
        getfasta = subprocess.Popen(
            ["bedtools", "getfasta", "-fi", genome_file,
             "-bed", bed_path, "-tab"],
            stdout=subprocess.PIPE, text=True,
        )

        # Step 4: Compute CpG frequency per window
        # Only uppercase ACGT are counted (matches original awk behavior —
        # soft-masked lowercase bases are treated as invalid).
        split_re = re.compile("[^ACGT]+")
        for (chrom, start, end), seq_line in zip(windows, getfasta.stdout):
            seq = seq_line.rstrip("\n").split("\t", 1)[1]

            segments = split_re.split(seq)
            total = sum(len(s) - 1 for s in segments if len(s) >= 2)
            if total == 0:
                continue

            cg = sum(s.count("CG") for s in segments)
            center = (start + end) // 2
            # %.6g matches awk's default OFMT
            print(f"{chrom}\t{center}\t{cg / total:.6g}")

        getfasta.wait()
        if getfasta.returncode != 0:
            sys.exit(getfasta.returncode)

    finally:
        os.unlink(bed_path)


if __name__ == "__main__":
    main()
