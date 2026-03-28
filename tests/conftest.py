"""Shared fixtures for dip-c test suite."""

import gzip
import os
import shutil

import h5py
import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Platform / tool markers
# ---------------------------------------------------------------------------

needs_merge_tools = pytest.mark.skipif(
    shutil.which("pigz") is None,
    reason="pigz not available (merge integration tests require GNU sort + pigz)",
)

needs_gunzip = pytest.mark.skipif(
    shutil.which("gunzip") is None,
    reason="gunzip not available",
)


# ---------------------------------------------------------------------------
# Matplotlib backend (session-wide)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session", autouse=True)
def _matplotlib_agg():
    """Force non-interactive Agg backend for all tests."""
    import matplotlib
    matplotlib.use("Agg")


# ---------------------------------------------------------------------------
# .cool file fixtures (read by hictkpy identically to .hic)
# ---------------------------------------------------------------------------

def _create_cool(path, counts, assembly="unknown"):
    """Create a minimal 2-chromosome .cool file.

    Chromosomes: chr1 (500 kb, 5 bins) and chr2 (300 kb, 3 bins)
    at 100 kb resolution.

    *counts* is a dict mapping ``(bin1_id, bin2_id)`` to an integer count.
    Only upper-triangle entries (bin1_id <= bin2_id) should be provided.
    """
    import cooler

    resolution = 100_000
    bins = pd.DataFrame({
        "chrom": ["chr1"] * 5 + ["chr2"] * 3,
        "start": list(range(0, 500_000, resolution))
                 + list(range(0, 300_000, resolution)),
        "end": list(range(resolution, 600_000, resolution))
               + list(range(resolution, 400_000, resolution)),
    })

    bin1_ids = sorted(counts.keys(), key=lambda k: (k[0], k[1]))
    pixels = pd.DataFrame({
        "bin1_id": [k[0] for k in bin1_ids],
        "bin2_id": [k[1] for k in bin1_ids],
        "count": [counts[k] for k in bin1_ids],
    })

    cooler.create_cooler(path, bins, pixels, assembly=assembly)
    return path


def _add_balancing_weights(path, weight_name, weights):
    """Add a normalization weight vector to a .cool file via h5py."""
    with h5py.File(path, "a") as f:
        f["bins"].create_dataset(weight_name, data=np.array(weights))


# Deterministic contact counts for the "A" fixture.
_COOL_A_COUNTS = {
    # chr1 intra (5x5 upper triangle)
    (0, 0): 100, (0, 1): 50, (1, 1): 200,
    (1, 2): 75,  (2, 2): 150, (3, 3): 60, (4, 4): 40,
    # chr2 intra (3x3 upper triangle)
    (5, 5): 80,  (5, 6): 40,  (6, 6): 120, (7, 7): 90,
    # chr1-chr2 inter
    (0, 5): 25,  (1, 6): 15,  (3, 7): 10,
}

# Different counts for the "B" fixture (diff tests).
_COOL_B_COUNTS = {
    (0, 0): 80,  (0, 1): 30,  (1, 1): 150,
    (1, 2): 60,  (2, 2): 100, (3, 3): 45, (4, 4): 55,
    (5, 5): 70,  (5, 6): 35,  (6, 6): 110, (7, 7): 85,
    (0, 5): 20,  (1, 6): 10,  (3, 7): 5,
}

# Simple balancing weights (8 bins total).
_SCALE_WEIGHTS = [1.0, 0.8, 1.2, 0.9, 1.1, 0.7, 1.3, 0.6]
_KR_WEIGHTS = [0.5, 0.6, 0.4, 0.7, 0.3, 0.8, 0.2, 0.9]


@pytest.fixture(scope="session")
def tiny_cool(tmp_path_factory):
    """A small .cool file with chr1 (500 kb) and chr2 (300 kb) at 100 kb.

    Has known contact counts plus KR and SCALE normalization weights.
    """
    path = str(tmp_path_factory.mktemp("cool") / "test_a.cool")
    _create_cool(path, _COOL_A_COUNTS, assembly="hg38")
    _add_balancing_weights(path, "SCALE", _SCALE_WEIGHTS)
    _add_balancing_weights(path, "KR", _KR_WEIGHTS)
    return path


@pytest.fixture(scope="session")
def tiny_cool_b(tmp_path_factory):
    """A second .cool file with different counts, for difference-map tests."""
    path = str(tmp_path_factory.mktemp("cool") / "test_b.cool")
    _create_cool(path, _COOL_B_COUNTS, assembly="hg38")
    _add_balancing_weights(path, "SCALE", _SCALE_WEIGHTS)
    _add_balancing_weights(path, "KR", _KR_WEIGHTS)
    return path


@pytest.fixture(scope="session")
def empty_cool(tmp_path_factory):
    """A .cool file where all contact counts are zero."""
    path = str(tmp_path_factory.mktemp("cool") / "test_empty.cool")
    _create_cool(path, {}, assembly="hg38")
    return path


@pytest.fixture(scope="session")
def partial_cool(tmp_path_factory):
    """A .cool file where chr1 has data but chr2 is all zeros.

    Used to test 'skipped chromosomes' warning branches.
    """
    counts = {
        (0, 0): 100, (0, 1): 50, (1, 1): 200,
        (1, 2): 75,  (2, 2): 150, (3, 3): 60, (4, 4): 40,
        # chr2 bins (5,6,7) have NO contacts
    }
    path = str(tmp_path_factory.mktemp("cool") / "test_partial.cool")
    _create_cool(path, counts, assembly="hg38")
    return path


# ---------------------------------------------------------------------------
# .pairs.gz fixtures
# ---------------------------------------------------------------------------

_DEFAULT_PAIRS_HEADER = [
    "## pairs format v1.0",
    "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2",
]

# Pre-sorted mm10 .pairs lines (chrom1 <= chrom2 in version-sort order).
_SAMPLE_PAIRS_LINES_MM10 = [
    "read1\tchr1\t1000\tchr1\t2000\t+\t-",
    "read2\tchr1\t3000\tchr1\t5000\t+\t+",
    "read3\tchr1\t8000\tchr2\t1000\t-\t+",
    "read4\tchr2\t2000\tchr2\t4000\t+\t-",
    "read5\tchr2\t6000\tchr3\t1000\t-\t-",
]


@pytest.fixture
def pairs_gz_factory(tmp_path):
    """Factory that creates .pairs.gz files with configurable content."""

    def _make(filename, lines=None, header=None):
        if header is None:
            header = list(_DEFAULT_PAIRS_HEADER)
        if lines is None:
            lines = list(_SAMPLE_PAIRS_LINES_MM10)
        path = tmp_path / filename
        with gzip.open(str(path), "wt") as f:
            for h in header:
                f.write(h + "\n")
            for line in lines:
                f.write(line + "\n")
        return str(path)

    return _make


@pytest.fixture
def pairs_dir_factory(tmp_path):
    """Factory that creates a directory with N .pairs.gz files."""

    def _make(n_files=3, lines_per_file=None):
        if lines_per_file is None:
            lines_per_file = list(_SAMPLE_PAIRS_LINES_MM10)
        d = tmp_path / "pairs_input"
        d.mkdir(exist_ok=True)
        paths = []
        for i in range(n_files):
            p = d / ("cell_%03d.pairs.gz" % i)
            with gzip.open(str(p), "wt") as f:
                for h in _DEFAULT_PAIRS_HEADER:
                    f.write(h + "\n")
                for line in lines_per_file:
                    f.write(line + "\n")
            paths.append(str(p))
        return str(d), paths

    return _make
