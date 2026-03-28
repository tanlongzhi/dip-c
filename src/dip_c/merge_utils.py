"""Shared utilities for hierarchical sort-merge of .pairs.gz files.

Provides resource detection, chromosome validation, batch sizing,
and file-estimation helpers used by the merge command.
"""

import os
import random
import shutil
import statistics
import subprocess
import sys


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BYTES_PER_CONTACT = 75

# Chromosome orders for .pairs files (matches what appears in actual data).
# NOTE: These differ from hicplot_utils._ORDERED_CHROMS which uses .hic
# file conventions (bare names for hg19/hg38).
PAIRS_CHROM_ORDER = {
    "mm10": [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chrX", "chrY",
    ],
    "hg19": [str(i) for i in range(1, 23)] + ["X", "Y"],
    "hg38": [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
        "chrX", "chrY",
    ],
}

# Sort key for .pairs files: chrom1 (version sort), chrom2 (version sort),
# pos1 (numeric), pos2 (numeric)
SORT_KEYS = ["-k2,2V", "-k4,4V", "-k3,3n", "-k5,5n"]


# ---------------------------------------------------------------------------
# Resource detection
# ---------------------------------------------------------------------------

def detect_cpus(cap=16):
    """Return usable CPU count, capped at *cap*."""
    try:
        n = len(os.sched_getaffinity(0))
    except AttributeError:
        n = os.cpu_count() or 4
    return min(n, cap)


# ---------------------------------------------------------------------------
# Tool checks
# ---------------------------------------------------------------------------

def check_required_tools():
    """Verify sort, pigz, gunzip are on PATH.  Exit if any missing."""
    missing = [t for t in ("sort", "pigz", "gunzip") if shutil.which(t) is None]
    if missing:
        sys.stderr.write(
            "[E::merge] Missing required tools: %s\n"
            "    Install them before running dip-c merge.\n"
            % ", ".join(missing)
        )
        raise SystemExit(1)


# ---------------------------------------------------------------------------
# Header extraction
# ---------------------------------------------------------------------------

def extract_header(pairs_gz_path):
    """Read header lines (starting with ``#``) from a .pairs.gz file.

    Returns the header as a string (with trailing newline).
    """
    proc = subprocess.run(
        ["gunzip", "-c", pairs_gz_path],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    lines = []
    for line in proc.stdout.decode("utf-8", errors="replace").splitlines():
        if line.startswith("#"):
            lines.append(line)
        else:
            break
    if not lines:
        sys.stderr.write(
            "[W::merge] No header lines found in %s\n" % pairs_gz_path
        )
    return "\n".join(lines) + "\n" if lines else ""


# ---------------------------------------------------------------------------
# Decompress helpers
# ---------------------------------------------------------------------------

def decompress_strip_header(pairs_gz, output_path):
    """``gunzip -c FILE | grep -v '^#' > output_path``."""
    with open(output_path, "w") as out:
        gunzip = subprocess.Popen(
            ["gunzip", "-c", pairs_gz],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        )
        grep = subprocess.Popen(
            ["grep", "-v", "^#"],
            stdin=gunzip.stdout, stdout=out, stderr=subprocess.DEVNULL,
        )
        gunzip.stdout.close()
        grep.wait()
        gunzip.wait()


# ---------------------------------------------------------------------------
# Line counting
# ---------------------------------------------------------------------------

def count_lines(path):
    """Fast line count via ``wc -l``."""
    result = subprocess.run(
        ["wc", "-l", path],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    return int(result.stdout.split()[0])


# ---------------------------------------------------------------------------
# File size estimation
# ---------------------------------------------------------------------------

def estimate_file_sizes(file_list):
    """Estimate median contacts from *file_list*.

    If fewer than 10 files, counts all.  Otherwise samples 10 random files.
    Returns ``(median, mean, stdev, min_val, max_val)``.
    """
    if len(file_list) < 10:
        counts = [count_lines(f) for f in file_list]
    else:
        sample = random.sample(file_list, 10)
        counts = [count_lines(f) for f in sample]

    if not counts:
        return (1000000, 1000000, 0, 1000000, 1000000)

    med = int(statistics.median(counts))
    mean = int(statistics.mean(counts))
    stdev = int(statistics.stdev(counts)) if len(counts) > 1 else 0
    return (med, mean, stdev, min(counts), max(counts))


# ---------------------------------------------------------------------------
# Batch size calculation
# ---------------------------------------------------------------------------

def compute_batch_size(num_files, median_contacts, sort_mem_per_job_gb,
                       parallel_jobs):
    """Compute optimal batch size balancing memory and CPU.

    Direct port of the bash script's adaptive algorithm.
    """
    # Memory-limited batch size
    mem_bytes = sort_mem_per_job_gb * 1024 * 1024 * 1024 * 0.8
    denom = median_contacts * BYTES_PER_CONTACT * 1.5
    if denom > 0:
        memory_limit = int(mem_bytes / denom)
    else:
        memory_limit = 1000
    if memory_limit <= 0:
        memory_limit = 1000

    # CPU-optimal batch size
    target_batches = parallel_jobs + parallel_jobs // 5
    cpu_optimal = num_files // target_batches if target_batches > 0 else num_files

    # Min/max clamps based on file count
    if num_files > 200:
        min_bs = 5
    elif num_files > 100:
        min_bs = 10
    else:
        min_bs = 15

    # Max clamp based on median contacts
    if median_contacts > 10000000:
        max_bs = 100
    elif median_contacts > 1000000:
        max_bs = 150
    else:
        max_bs = 500

    # Pick the smaller of memory and CPU limits
    batch_size = min(memory_limit, cpu_optimal)
    batch_size = max(batch_size, min_bs)
    batch_size = min(batch_size, max_bs)

    # Can't be larger than file count
    batch_size = min(batch_size, num_files)

    return batch_size, memory_limit, cpu_optimal


# ---------------------------------------------------------------------------
# Grouping factors for hierarchical merge
# ---------------------------------------------------------------------------

def choose_grouping_factors(num_batches, median_contacts):
    """Select merge grouping factors per level.

    Returns a list of ints — factor[i] is the group size at merge level i+1.
    """
    if num_batches > 500:
        large = [10, 15, 20, 25, 30]
        medium = [8, 12, 18, 25, 30]
        small = [5, 10, 15, 20, 25]
    elif num_batches > 200:
        large = [8, 12, 18, 25, 30]
        medium = [5, 10, 15, 20, 25]
        small = [5, 8, 12, 18, 25]
    elif num_batches > 50:
        large = [5, 10, 15, 20]
        medium = [5, 8, 12, 18]
        small = [5, 8, 12, 15]
    else:
        large = [5, 10, 15]
        medium = [5, 8, 12]
        small = [5, 8, 10]

    if median_contacts > 5000000:
        return small
    elif median_contacts > 1000000:
        return medium
    else:
        return large


# ---------------------------------------------------------------------------
# Chromosome validation
# ---------------------------------------------------------------------------

def chrom_index_map(genome):
    """Return ``{chrom_name: index}`` for *genome*.

    Returns ``None`` for ``genome='any'`` (skips chrom-order validation).
    """
    if genome == "any":
        return None
    order = PAIRS_CHROM_ORDER.get(genome)
    if order is None:
        sys.stderr.write(
            "[E::merge] Unknown genome '%s'. "
            "Supported: %s, or 'any' to skip chrom validation.\n"
            % (genome, ", ".join(sorted(PAIRS_CHROM_ORDER.keys())))
        )
        raise SystemExit(1)
    return {c: i for i, c in enumerate(order)}


def chk_chrom_and_fields(path, chrom_idx):
    """Filter *path* in-place: keep only valid .pairs lines.

    Checks:
    - Exactly 7 tab-separated fields
    - If *chrom_idx* is not None: fields[1] and fields[3] must be known
      chromosomes, and chrom_idx[fields[1]] <= chrom_idx[fields[3]]

    Returns the number of removed lines.
    """
    tmp = path + ".tmp"
    removed = 0
    with open(path, "r") as fin, open(tmp, "w") as fout:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 7:
                removed += 1
                continue
            if chrom_idx is not None:
                c1, c2 = fields[1], fields[3]
                if c1 not in chrom_idx or c2 not in chrom_idx:
                    removed += 1
                    continue
                if chrom_idx[c1] > chrom_idx[c2]:
                    removed += 1
                    continue
            fout.write(line)

    os.replace(tmp, path)
    return removed
