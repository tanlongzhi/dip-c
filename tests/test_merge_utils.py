"""Unit tests for dip_c.merge_utils."""

import os
import shutil
import sys

import pytest

from dip_c.merge_utils import (
    PAIRS_CHROM_ORDER,
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

needs_gunzip = pytest.mark.skipif(
    shutil.which("gunzip") is None,
    reason="gunzip not available",
)


# ===================================================================
# detect_cpus
# ===================================================================

class TestDetectCpus:
    def test_returns_positive_int(self):
        n = detect_cpus()
        assert isinstance(n, int)
        assert n > 0

    def test_respects_cap(self):
        n = detect_cpus(cap=1)
        assert n == 1

    def test_fallback_no_sched_getaffinity(self, monkeypatch):
        monkeypatch.delattr(os, "sched_getaffinity", raising=False)
        n = detect_cpus()
        assert n > 0

    def test_fallback_cpu_count_none(self, monkeypatch):
        monkeypatch.delattr(os, "sched_getaffinity", raising=False)
        monkeypatch.setattr(os, "cpu_count", lambda: None)
        n = detect_cpus()
        assert n == 4  # default fallback


# ===================================================================
# check_required_tools
# ===================================================================

class TestCheckRequiredTools:
    def test_all_present(self, monkeypatch):
        monkeypatch.setattr(shutil, "which",
                            lambda t: "/usr/bin/" + t)
        check_required_tools()  # should not raise

    def test_missing_pigz(self, monkeypatch):
        original = shutil.which

        def fake_which(t):
            if t == "pigz":
                return None
            return original(t)

        monkeypatch.setattr(shutil, "which", fake_which)
        with pytest.raises(SystemExit) as exc_info:
            check_required_tools()
        assert exc_info.value.code == 1

    def test_missing_all(self, monkeypatch):
        monkeypatch.setattr(shutil, "which", lambda t: None)
        with pytest.raises(SystemExit):
            check_required_tools()


# ===================================================================
# extract_header
# ===================================================================

class TestExtractHeader:
    @needs_gunzip
    def test_normal_header(self, pairs_gz_factory):
        path = pairs_gz_factory("test.pairs.gz")
        header = extract_header(path)
        assert header.startswith("##")
        assert "#columns:" in header

    @needs_gunzip
    def test_no_header_warns(self, pairs_gz_factory, capsys):
        path = pairs_gz_factory("noheader.pairs.gz", header=[])
        header = extract_header(path)
        assert header == ""
        captured = capsys.readouterr()
        assert "No header lines" in captured.err


# ===================================================================
# decompress_strip_header
# ===================================================================

class TestDecompressStripHeader:
    @needs_gunzip
    def test_strips_header(self, pairs_gz_factory, tmp_path):
        path = pairs_gz_factory("test.pairs.gz")
        out = str(tmp_path / "out.pairs")
        decompress_strip_header(path, out)

        with open(out) as f:
            lines = f.readlines()
        # No header lines should remain
        assert all(not line.startswith("#") for line in lines)
        assert len(lines) > 0


# ===================================================================
# count_lines
# ===================================================================

class TestCountLines:
    def test_known_count(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("a\nb\nc\n")
        assert count_lines(str(f)) == 3

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_text("")
        assert count_lines(str(f)) == 0


# ===================================================================
# estimate_file_sizes
# ===================================================================

class TestEstimateFileSizes:
    def test_small_set(self, tmp_path):
        files = []
        for i, n_lines in enumerate([10, 20, 30]):
            p = tmp_path / ("f%d.txt" % i)
            p.write_text("\n".join(["x"] * n_lines) + "\n")
            files.append(str(p))

        med, mean, stdev, minv, maxv = estimate_file_sizes(files)
        assert med == 20
        assert minv == 10
        assert maxv == 30

    def test_large_set_samples(self, tmp_path):
        """With >= 10 files, only 10 are sampled."""
        files = []
        for i in range(15):
            p = tmp_path / ("f%d.txt" % i)
            p.write_text("\n".join(["x"] * 100) + "\n")
            files.append(str(p))

        med, mean, stdev, minv, maxv = estimate_file_sizes(files)
        assert med == 100
        assert mean == 100

    def test_empty_list(self):
        med, mean, stdev, minv, maxv = estimate_file_sizes([])
        assert med == 1000000  # default

    def test_single_file(self, tmp_path):
        p = tmp_path / "f.txt"
        p.write_text("a\nb\nc\n")
        med, mean, stdev, minv, maxv = estimate_file_sizes([str(p)])
        assert med == 3
        assert stdev == 0


# ===================================================================
# compute_batch_size — parametrized
# ===================================================================

class TestComputeBatchSize:
    def test_normal(self):
        bs, ml, co = compute_batch_size(100, 500000, 4, 8)
        assert bs > 0
        assert bs <= 100

    def test_zero_median(self):
        """denom=0 -> memory_limit=1000."""
        bs, ml, co = compute_batch_size(50, 0, 4, 8)
        assert ml == 1000

    def test_many_files(self):
        bs, _, _ = compute_batch_size(300, 500000, 4, 8)
        assert bs >= 5  # min_bs for >200 files

    def test_medium_files(self):
        bs, _, _ = compute_batch_size(150, 500000, 4, 8)
        assert bs >= 10  # min_bs for >100 files

    def test_few_files(self):
        bs, _, _ = compute_batch_size(50, 500000, 4, 8)
        assert bs >= 15  # min_bs for <=100 files

    def test_large_median(self):
        bs, _, _ = compute_batch_size(100, 20000000, 4, 8)
        assert bs <= 100  # max_bs for >10M median

    def test_medium_median(self):
        bs, _, _ = compute_batch_size(100, 5000000, 4, 8)
        assert bs <= 150  # max_bs for >1M median

    def test_small_median(self):
        bs, _, _ = compute_batch_size(100, 100000, 4, 8)
        assert bs <= 500  # max_bs for <=1M median

    def test_zero_jobs(self):
        bs, ml, co = compute_batch_size(100, 500000, 4, 0)
        assert co == 100  # target_batches=0 -> cpu_optimal=num_files

    def test_batch_capped_at_file_count(self):
        bs, _, _ = compute_batch_size(3, 100, 4, 8)
        assert bs <= 3

    def test_negative_memory_limit(self):
        """Huge median + tiny memory -> memory_limit rounds to 0, clamped to 1000."""
        bs, ml, _ = compute_batch_size(100, 999999999, 0, 8)
        assert ml == 1000


# ===================================================================
# choose_grouping_factors — parametrized
# ===================================================================

class TestChooseGroupingFactors:
    @pytest.mark.parametrize("num_batches,median,expected_len", [
        # >500 batches
        (600, 10000000, 5),   # large median -> small factors
        (600, 3000000, 5),    # medium median
        (600, 500000, 5),     # small median -> large factors
        # 200-500 batches
        (300, 10000000, 5),
        (300, 3000000, 5),
        (300, 500000, 5),
        # 50-200 batches
        (100, 10000000, 4),
        (100, 3000000, 4),
        (100, 500000, 4),
        # <=50 batches
        (30, 10000000, 3),
        (30, 3000000, 3),
        (30, 500000, 3),
    ])
    def test_factor_selection(self, num_batches, median, expected_len):
        factors = choose_grouping_factors(num_batches, median)
        assert len(factors) == expected_len
        assert all(isinstance(f, int) and f > 0 for f in factors)


# ===================================================================
# chrom_index_map
# ===================================================================

class TestChromIndexMap:
    def test_any_returns_none(self):
        assert chrom_index_map("any") is None

    def test_mm10(self):
        idx = chrom_index_map("mm10")
        assert isinstance(idx, dict)
        assert "chr1" in idx
        assert "chrX" in idx
        assert idx["chr1"] < idx["chrX"]

    def test_hg19(self):
        idx = chrom_index_map("hg19")
        assert "1" in idx
        assert "X" in idx

    def test_hg38(self):
        idx = chrom_index_map("hg38")
        assert "chr1" in idx
        assert "chrY" in idx

    def test_unknown_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            chrom_index_map("hg99")
        assert exc_info.value.code == 1


# ===================================================================
# chk_chrom_and_fields
# ===================================================================

class TestChkChromAndFields:
    def test_valid_lines_kept(self, tmp_path):
        f = tmp_path / "test.pairs"
        f.write_text(
            "read1\tchr1\t100\tchr1\t200\t+\t-\n"
            "read2\tchr1\t300\tchr2\t400\t+\t+\n"
        )
        idx = chrom_index_map("mm10")
        removed = chk_chrom_and_fields(str(f), idx)
        assert removed == 0
        assert f.read_text().count("\n") == 2

    def test_wrong_field_count(self, tmp_path):
        f = tmp_path / "test.pairs"
        f.write_text(
            "read1\tchr1\t100\tchr1\t200\t+\t-\n"
            "bad_line_too_few_fields\n"
        )
        idx = chrom_index_map("mm10")
        removed = chk_chrom_and_fields(str(f), idx)
        assert removed == 1

    def test_unknown_chrom(self, tmp_path):
        f = tmp_path / "test.pairs"
        f.write_text("read1\tchrZ\t100\tchr1\t200\t+\t-\n")
        idx = chrom_index_map("mm10")
        removed = chk_chrom_and_fields(str(f), idx)
        assert removed == 1
        assert f.read_text() == ""

    def test_order_violation(self, tmp_path):
        f = tmp_path / "test.pairs"
        # chr2 before chr1 violates ordering
        f.write_text("read1\tchr2\t100\tchr1\t200\t+\t-\n")
        idx = chrom_index_map("mm10")
        removed = chk_chrom_and_fields(str(f), idx)
        assert removed == 1

    def test_no_chrom_validation(self, tmp_path):
        """chrom_idx=None -> only field count is checked."""
        f = tmp_path / "test.pairs"
        f.write_text(
            "read1\tchrANY\t100\tchrANY\t200\t+\t-\n"
            "bad_line\n"
        )
        removed = chk_chrom_and_fields(str(f), None)
        assert removed == 1
        lines = f.read_text().strip().split("\n")
        assert len(lines) == 1
        assert "chrANY" in lines[0]

    def test_inplace_replacement(self, tmp_path):
        """Verify the file is replaced in-place (no .tmp left behind)."""
        f = tmp_path / "test.pairs"
        f.write_text("read1\tchr1\t100\tchr1\t200\t+\t-\n")
        chk_chrom_and_fields(str(f), chrom_index_map("mm10"))
        assert not (tmp_path / "test.pairs.tmp").exists()
