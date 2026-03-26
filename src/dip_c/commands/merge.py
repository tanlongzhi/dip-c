"""Hierarchical sort-merge of .pairs.gz files.

Merges all .pairs.gz files in a directory into a single sorted,
compressed .pairs.gz file using GNU sort (merge mode) and pigz.

Requires:  sort, gunzip, pigz  (standard on HPC systems)

Usage:
    dip-c merge <dir> -g <genome> [-o <output>] [-H <header>]
                [-j <jobs>] [-m <memory>]
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from dip_c.merge_utils import (
    SORT_KEYS,
    check_required_tools,
    chk_chrom_and_fields,
    chrom_index_map,
    choose_grouping_factors,
    compute_batch_size,
    count_lines,
    decompress_strip_header,
    detect_cpus,
    estimate_file_sizes,
    extract_header,
)


# ══════════════════════════════════════════════════════════════════════════
# Logging
# ══════════════════════════════════════════════════════════════════════════

def _log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    sys.stderr.write("[M::merge] [%s] %s\n" % (ts, msg))


def _elapsed(start):
    secs = int(time.time() - start)
    h, m, s = secs // 3600, (secs % 3600) // 60, secs % 60
    return "%02d:%02d:%02d" % (h, m, s)


# ══════════════════════════════════════════════════════════════════════════
# Sort-merge helper (called in worker processes)
# ══════════════════════════════════════════════════════════════════════════

def _sort_merge_batch(file_list, output_path, sort_mem_gb, tmpdir):
    """Run ``sort -m`` on a list of pre-sorted files."""
    cmd = (
        ["sort", "-m", "-S", "%dG" % sort_mem_gb, "-T", tmpdir]
        + SORT_KEYS
        + file_list
    )
    with open(output_path, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)
    return output_path


# ══════════════════════════════════════════════════════════════════════════
# Phase implementations
# ══════════════════════════════════════════════════════════════════════════

def _phase1_decompress(pairs_files, work_dir, max_jobs):
    """Decompress .pairs.gz files, stripping headers."""
    _log("Phase 1: Decompressing %d .pairs.gz files..." % len(pairs_files))

    noheader_files = []
    with ProcessPoolExecutor(max_workers=max_jobs) as pool:
        futures = {}
        for pf in pairs_files:
            base = os.path.basename(pf).replace(".pairs.gz", ".noheader.pairs")
            out = os.path.join(work_dir, base)
            noheader_files.append(out)
            futures[pool.submit(decompress_strip_header, pf, out)] = pf

        for fut in as_completed(futures):
            fut.result()  # raise on error

    _log("Phase 1: Decompressed %d files" % len(noheader_files))
    return noheader_files


def _phase2_estimate(noheader_files):
    """Estimate file sizes by sampling."""
    _log("Phase 2: Estimating file sizes...")
    median, mean, stdev, minv, maxv = estimate_file_sizes(noheader_files)

    if len(noheader_files) < 10:
        _log("All %d files counted. Median: %d contacts"
             % (len(noheader_files), median))
    else:
        _log("Sample statistics (10 random files):")
        sys.stderr.write("  Mean: %d contacts\n" % mean)
        sys.stderr.write("  StdDev: %d contacts\n" % stdev)
        sys.stderr.write("  Median: %d contacts\n" % median)
        sys.stderr.write("  Range: %d - %d contacts\n" % (minv, maxv))

    return median


def _phase3_batch_size(num_files, median, sort_mem_per_job, jobs):
    """Compute optimal batch size."""
    batch_size, mem_limit, cpu_opt = compute_batch_size(
        num_files, median, sort_mem_per_job, jobs,
    )
    _log("Phase 3: Batch size calculation:")
    sys.stderr.write("  Memory limit batch size: %d\n" % mem_limit)
    sys.stderr.write("  CPU optimal batch size: %d\n" % cpu_opt)
    sys.stderr.write("  Selected batch size: %d\n" % batch_size)
    sys.stderr.write("  Expected batches: %d\n"
                     % ((num_files + batch_size - 1) // batch_size))
    return batch_size


def _phase4_initial_sort(noheader_files, batch_size, sort_mem_gb,
                         max_jobs, work_dir):
    """Split into batches and parallel sort-merge each batch."""
    # Chunk the file list
    batches = []
    for i in range(0, len(noheader_files), batch_size):
        batches.append(noheader_files[i:i + batch_size])

    _log("Phase 4: Sorting %d files in %d batches..."
         % (len(noheader_files), len(batches)))

    tmpdir = os.path.join(work_dir, "sort_tmp")
    os.makedirs(tmpdir, exist_ok=True)

    sorted_paths = []
    with ProcessPoolExecutor(max_workers=max_jobs) as pool:
        futures = {}
        for idx, batch in enumerate(batches):
            out = os.path.join(work_dir, "batch_%04d_sorted" % idx)
            sorted_paths.append(out)
            futures[pool.submit(
                _sort_merge_batch, batch, out, sort_mem_gb, tmpdir,
            )] = idx

        for fut in as_completed(futures):
            fut.result()  # raise on error

    _log("Phase 4: Batch sort complete (%d sorted batches)" % len(sorted_paths))
    return sorted_paths


def _phase5_validate(sorted_batches, genome, max_jobs):
    """Chromosome validation on each sorted batch."""
    chrom_idx = chrom_index_map(genome)
    if chrom_idx is None:
        _log("Phase 5: genome=any, checking field count only")
    else:
        _log("Phase 5: Validating chromosome ordering (%s)..." % genome)

    total_removed = 0
    with ProcessPoolExecutor(max_workers=max_jobs) as pool:
        futures = {
            pool.submit(chk_chrom_and_fields, path, chrom_idx): path
            for path in sorted_batches
        }
        for fut in as_completed(futures):
            removed = fut.result()
            if removed > 0:
                path = futures[fut]
                sys.stderr.write("  %s: %d lines removed\n"
                                 % (os.path.basename(path), removed))
            total_removed += removed

    _log("Phase 5: Validation complete (%d total lines removed)"
         % total_removed)


def _phase6_hierarchical_merge(sorted_batches, median, sort_mem_gb,
                                max_jobs, work_dir):
    """Multi-level hierarchical merge using sort -m."""
    factors = choose_grouping_factors(len(sorted_batches), median)
    _log("Phase 6: Hierarchical merge (%d batches, factors: %s)"
         % (len(sorted_batches), factors))

    tmpdir = os.path.join(work_dir, "sort_tmp")
    os.makedirs(tmpdir, exist_ok=True)

    current_files = list(sorted_batches)
    level = 1

    while len(current_files) > 1:
        factor_idx = min(level - 1, len(factors) - 1)
        group_size = factors[factor_idx]
        group_size = min(group_size, len(current_files))

        _log("Level %d: %d files, group size %d"
             % (level, len(current_files), group_size))

        level_dir = os.path.join(work_dir, "merge_level%d" % level)
        os.makedirs(level_dir, exist_ok=True)

        # Build groups
        groups = []
        for i in range(0, len(current_files), group_size):
            groups.append(current_files[i:i + group_size])

        # Parallel merge
        next_files = []
        with ProcessPoolExecutor(max_workers=max_jobs) as pool:
            futures = {}
            for idx, group in enumerate(groups):
                out = os.path.join(level_dir, "merged_%04d" % idx)
                next_files.append(out)
                futures[pool.submit(
                    _sort_merge_batch, group, out, sort_mem_gb, tmpdir,
                )] = idx

            for fut in as_completed(futures):
                fut.result()

        # Cleanup previous level
        if level == 1:
            for f in sorted_batches:
                if os.path.exists(f):
                    os.remove(f)
        else:
            prev_dir = os.path.join(work_dir, "merge_level%d" % (level - 1))
            if os.path.isdir(prev_dir):
                shutil.rmtree(prev_dir)

        current_files = next_files
        level += 1

    return current_files[0]


def _phase7_compress(final_sorted, header_text, output_path, pigz_threads):
    """Prepend header and compress with pigz."""
    _log("Phase 7: Compressing with pigz (%d threads)..." % pigz_threads)

    with open(output_path, "wb") as out:
        pigz = subprocess.Popen(
            ["pigz", "-p", str(pigz_threads)],
            stdin=subprocess.PIPE, stdout=out, stderr=subprocess.DEVNULL,
        )
        # Write header
        pigz.stdin.write(header_text.encode("utf-8"))
        # Stream sorted data
        with open(final_sorted, "rb") as f:
            shutil.copyfileobj(f, pigz.stdin, length=1024 * 1024)
        pigz.stdin.close()
        pigz.wait()

    if pigz.returncode != 0:
        sys.stderr.write("[E::merge] pigz failed with exit code %d\n"
                         % pigz.returncode)
        raise SystemExit(1)


# ══════════════════════════════════════════════════════════════════════════
# Argument parser
# ══════════════════════════════════════════════════════════════════════════

VALID_GENOMES = list(sorted(
    ["mm10", "hg19", "hg38", "any"]
))

def _build_parser():
    p = argparse.ArgumentParser(
        prog="dip-c merge",
        description="Hierarchical sort-merge of .pairs.gz files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
        add_help=False,
        epilog="""\
Usage:
  dip-c merge -i <dir> -o <output> [-g <genome>]
              [-h <header.pairs>] [-j <max_jobs>] [-m <max_mem_gb>]

Required:
  -i <dir>            Input directory containing .pairs.gz files
  -o <output>         Output file name (.pairs.gz appended if not present)

Optional:
  -g <genome>         Genome ID: mm10, hg19, hg38
                      If omitted, chrom-order validation is skipped
  -h <header>         Header file (default: auto-extracted from first file)
  -j <max_jobs>       Max parallel jobs (default: auto-detect CPUs, cap 16)
  -m <max_mem_gb>     Max memory in GB (default: 32)

Resource limits (-j, -m) are upper bounds for the adaptive algorithm.
Sort memory per job = min(total_mem / jobs, 30 GB).
Batch sizes and grouping factors are computed dynamically from file sizes.

Examples:
  # Basic merge (auto-detect resources)
  dip-c merge -i /path/to/pairsgz_dir -g mm10 -o merged

  # Specify resource limits
  dip-c merge -i /path/to/pairsgz_dir -g mm10 -o merged.pairs.gz -j 8 -m 64

  # Custom header file
  dip-c merge -i /path/to/pairsgz_dir -g hg38 -o merged -h header.pairs

  # Skip chrom validation (omit -g)
  dip-c merge -i /path/to/pairsgz_dir -o merged
""",
    )

    p.add_argument(
        "-i", "--input", dest="dir", required=True,
        metavar="DIR",
        help="Input directory containing .pairs.gz files to merge.",
    )
    p.add_argument(
        "-g", "--genome", required=False, default=None,
        metavar="GENOME",
        help="Genome ID: mm10, hg19, hg38. "
             "If omitted, chromosome-order validation is skipped "
             "(only field count is checked).",
    )
    p.add_argument(
        "-o", "--output", required=True,
        metavar="FILE",
        help="Output file name. If it does not end with .pairs.gz, "
             "the suffix is appended automatically.",
    )
    p.add_argument(
        "-h", "--header", default=None,
        metavar="FILE",
        help="Header file to prepend. "
             "Default: auto-extracted from the first .pairs.gz file.",
    )
    p.add_argument(
        "-j", "--jobs", type=int, default=None,
        metavar="N",
        help="Max parallel jobs (default: auto-detect, capped at 16).",
    )
    p.add_argument(
        "-m", "--memory", type=int, default=32,
        metavar="GB",
        help="Max memory in GB (default: 32).",
    )
    p.add_argument(
        "--help", action="help",
        help="Show this help message and exit.",
    )

    return p


# ══════════════════════════════════════════════════════════════════════════
# CLI entry point  –  called from dip_c.cli as  ``dip-c merge …``
# ══════════════════════════════════════════════════════════════════════════

def merge(argv):
    """Parse flags and run the hierarchical sort-merge pipeline."""
    parser = _build_parser()
    args = parser.parse_args(argv[1:])

    # -- Validate inputs ---------------------------------------------------
    pairs_dir = os.path.abspath(args.dir)
    if not os.path.isdir(pairs_dir):
        parser.error("Directory does not exist: %s" % pairs_dir)

    pairs_files = sorted(glob.glob(os.path.join(pairs_dir, "*.pairs.gz")))
    if not pairs_files:
        parser.error("No .pairs.gz files found in %s" % pairs_dir)

    genome = args.genome
    if genome is None:
        genome = "any"
    if genome not in VALID_GENOMES:
        parser.error(
            "Unknown genome '%s'. Supported: %s (or omit -g to skip validation)"
            % (genome, ", ".join(g for g in VALID_GENOMES if g != "any"))
        )

    output = args.output
    if not output.endswith(".pairs.gz"):
        output = output + ".pairs.gz"
    max_jobs = args.jobs or detect_cpus()
    max_mem = args.memory
    sort_mem_per_job = min(max_mem // max_jobs, 30) if max_jobs > 0 else 30

    # -- Check tools -------------------------------------------------------
    check_required_tools()

    # -- Log setup ---------------------------------------------------------
    start_time = time.time()
    _log("Merge started")
    _log("Directory: %s" % pairs_dir)
    _log("Output: %s" % output)
    _log("Genome: %s" % genome)
    _log("Found %d .pairs.gz files" % len(pairs_files))
    _log("Resource limits: %d jobs, %d GB memory (%d GB/job for sort)"
         % (max_jobs, max_mem, sort_mem_per_job))

    # -- Single-file shortcut ----------------------------------------------
    if len(pairs_files) == 1:
        _log("Only 1 file found, copying to output")
        shutil.copy2(pairs_files[0], output)
        _log("Done in %s" % _elapsed(start_time))
        return 0

    # -- Header extraction -------------------------------------------------
    if args.header:
        with open(args.header, "r") as f:
            header_text = f.read()
        _log("Header loaded from %s" % args.header)
    else:
        header_text = extract_header(pairs_files[0])
        _log("Header auto-extracted from %s (%d lines)"
             % (os.path.basename(pairs_files[0]),
                header_text.count("\n")))

    # -- Create work directory ---------------------------------------------
    work_dir = tempfile.mkdtemp(prefix=".merge_work_", dir=pairs_dir)
    _log("Work directory: %s" % work_dir)

    try:
        # Phase 1: Decompress
        noheader_files = _phase1_decompress(pairs_files, work_dir, max_jobs)

        # Phase 2: Estimate sizes
        median = _phase2_estimate(noheader_files)

        # Phase 3: Batch size
        batch_size = _phase3_batch_size(
            len(noheader_files), median, sort_mem_per_job, max_jobs,
        )

        # Phase 4: Initial batch sort
        sorted_batches = _phase4_initial_sort(
            noheader_files, batch_size, sort_mem_per_job, max_jobs, work_dir,
        )

        # Clean up decompressed files (no longer needed)
        for f in noheader_files:
            if os.path.exists(f):
                os.remove(f)

        # Phase 5: Chromosome validation
        _phase5_validate(sorted_batches, genome, max_jobs)

        # Phase 6: Hierarchical merge
        final_sorted = _phase6_hierarchical_merge(
            sorted_batches, median, sort_mem_per_job, max_jobs, work_dir,
        )

        # Count contacts
        total_contacts = count_lines(final_sorted)
        _log("Total contacts: %d" % total_contacts)

        # Phase 7: Compress
        pigz_threads = min(max_jobs, 16)
        _phase7_compress(final_sorted, header_text, output, pigz_threads)

    finally:
        # Cleanup work directory
        if os.path.isdir(work_dir):
            shutil.rmtree(work_dir)
            _log("Work directory cleaned up")

    # -- Final stats -------------------------------------------------------
    if os.path.exists(output):
        size_mb = os.path.getsize(output) / (1024 * 1024)
        _log("Output: %s (%.1f MB, %d contacts)" % (output, size_mb, total_contacts))
    else:
        sys.stderr.write("[E::merge] Output file not created: %s\n" % output)
        return 1

    _log("Merge completed in %s" % _elapsed(start_time))
    return 0
