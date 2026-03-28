"""Tests for dip_c.commands.merge.

Most tests mock _sort_merge_batch (the only function that calls GNU sort -m)
and pigz so they work on macOS and Linux without requiring these tools.
Integration tests marked @needs_merge_tools run with real tools on Sherlock.
"""

import gzip
import os
import shutil
import subprocess
import sys
import time

import pytest

from dip_c.commands.merge import (
    _build_parser,
    _elapsed,
    _log,
    _phase1_decompress,
    _phase2_estimate,
    _phase3_batch_size,
    _phase4_initial_sort,
    _phase5_validate,
    _phase6_hierarchical_merge,
    _phase7_compress,
    _sort_merge_batch,
    merge,
)

needs_merge_tools = pytest.mark.skipif(
    shutil.which("pigz") is None,
    reason="pigz not available (integration tests need GNU sort + pigz)",
)

needs_gunzip = pytest.mark.skipif(
    shutil.which("gunzip") is None,
    reason="gunzip not available",
)

# Pre-sorted .pairs lines (mm10 chromosome order, tab-separated, 7 fields).
SORTED_LINES = [
    "r1\tchr1\t1000\tchr1\t2000\t+\t-",
    "r2\tchr1\t3000\tchr1\t5000\t+\t+",
    "r3\tchr1\t8000\tchr2\t1000\t-\t+",
    "r4\tchr2\t2000\tchr2\t4000\t+\t-",
    "r5\tchr2\t6000\tchr3\t1000\t-\t-",
]

HEADER_LINES = [
    "## pairs format v1.0",
    "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2",
]


def _make_pairs_gz(path, lines=None, header=None):
    """Write a .pairs.gz file with given header and data lines."""
    if header is None:
        header = list(HEADER_LINES)
    if lines is None:
        lines = list(SORTED_LINES)
    with gzip.open(str(path), "wt") as f:
        for h in header:
            f.write(h + "\n")
        for line in lines:
            f.write(line + "\n")


def _make_pairs_dir(base_dir, n_files=3, lines=None):
    """Create a directory with N .pairs.gz files."""
    d = base_dir / "pairs_input"
    d.mkdir(exist_ok=True)
    for i in range(n_files):
        _make_pairs_gz(d / ("cell_%03d.pairs.gz" % i), lines=lines)
    return str(d)


def _fake_sort_merge_batch(file_list, output_path, sort_mem_gb, tmpdir):
    """Mock for _sort_merge_batch: concatenates input files (no real sort)."""
    with open(output_path, "w") as out:
        for f in file_list:
            with open(f) as inp:
                out.write(inp.read())
    return output_path


def _fake_phase7_compress(final_sorted, header_text, output_path,
                          pigz_threads):
    """Mock for _phase7_compress: gzip without pigz."""
    with gzip.open(output_path, "wt") as out:
        out.write(header_text)
        with open(final_sorted) as f:
            out.write(f.read())


# ===================================================================
# Logging helpers
# ===================================================================

class TestLogging:
    def test_log_writes_timestamp(self, capsys):
        _log("test message")
        captured = capsys.readouterr()
        assert "[M::merge]" in captured.err
        assert "test message" in captured.err

    def test_elapsed_format(self):
        assert _elapsed(time.time() - 3661) == "01:01:01"
        assert _elapsed(time.time()) == "00:00:00"


# ===================================================================
# _sort_merge_batch
# ===================================================================

class TestSortMergeBatch:
    @needs_merge_tools
    def test_merges_files_real(self, tmp_path):
        """Integration test with real sort -m."""
        f1 = tmp_path / "a.txt"
        f1.write_text("r1\tchr1\t100\tchr1\t200\t+\t-\n")
        f2 = tmp_path / "b.txt"
        f2.write_text("r2\tchr1\t300\tchr1\t400\t+\t+\n")
        out = str(tmp_path / "merged.txt")
        tmpdir = str(tmp_path / "tmp")
        os.makedirs(tmpdir)

        _sort_merge_batch([str(f1), str(f2)], out, 1, tmpdir)
        assert os.path.exists(out)
        with open(out) as f:
            lines = f.readlines()
        assert len(lines) == 2

    def test_subprocess_call(self, tmp_path, monkeypatch):
        """Cover _sort_merge_batch body by mocking subprocess.run."""
        f1 = tmp_path / "a.txt"
        f1.write_text("line1\n")
        out = str(tmp_path / "out.txt")
        tmpdir = str(tmp_path / "tmp")
        os.makedirs(tmpdir)

        def fake_run(cmd, **kw):
            # Simulate sort by copying input to stdout
            stdout = kw.get("stdout")
            if stdout:
                for fp in cmd[len(cmd) - 1:]:  # last arg is the file
                    with open(fp) as inp:
                        stdout.write(inp.read())

        monkeypatch.setattr(subprocess, "run", fake_run)
        result = _sort_merge_batch([str(f1)], out, 1, tmpdir)
        assert result == out
        assert os.path.exists(out)


# ===================================================================
# Phase tests
# ===================================================================

class TestPhase1:
    @needs_gunzip
    def test_decompresses_all(self, tmp_path):
        d = _make_pairs_dir(tmp_path, n_files=2)
        work = str(tmp_path / "work")
        os.makedirs(work)
        files = sorted(
            str(tmp_path / "pairs_input" / f)
            for f in os.listdir(d) if f.endswith(".pairs.gz")
        )
        result = _phase1_decompress(files, work, 2)
        assert len(result) == 2
        for f in result:
            assert os.path.exists(f)
            with open(f) as fh:
                content = fh.read()
            assert not content.startswith("#")


class TestPhase2:
    def test_few_files(self, tmp_path, capsys):
        files = []
        for i in range(3):
            p = tmp_path / ("f%d.txt" % i)
            p.write_text("\n".join(["x"] * 100) + "\n")
            files.append(str(p))

        median = _phase2_estimate(files)
        assert median == 100
        captured = capsys.readouterr()
        assert "All 3 files counted" in captured.err

    def test_many_files(self, tmp_path, capsys):
        files = []
        for i in range(12):
            p = tmp_path / ("f%d.txt" % i)
            p.write_text("\n".join(["x"] * 50) + "\n")
            files.append(str(p))

        median = _phase2_estimate(files)
        assert median == 50
        captured = capsys.readouterr()
        assert "Sample statistics" in captured.err


class TestPhase3:
    def test_returns_batch_size(self, capsys):
        bs = _phase3_batch_size(100, 500000, 4, 8)
        assert bs > 0
        captured = capsys.readouterr()
        assert "Batch size calculation" in captured.err


class TestPhase4:
    def test_creates_sorted_batches(self, tmp_path, monkeypatch):
        """Test phase4 orchestration with mocked sort."""
        monkeypatch.setattr(
            "dip_c.commands.merge._sort_merge_batch",
            _fake_sort_merge_batch)

        files = []
        for i in range(4):
            p = tmp_path / ("f%d.txt" % i)
            p.write_text("r%d\tchr1\t%d\tchr1\t%d\t+\t-\n"
                         % (i, i * 100, i * 200))
            files.append(str(p))

        work = str(tmp_path / "work")
        os.makedirs(work)
        result = _phase4_initial_sort(files, 2, 1, 2, work)
        assert len(result) == 2
        for f in result:
            assert os.path.exists(f)

    @needs_merge_tools
    def test_creates_sorted_batches_real(self, tmp_path):
        """Integration: phase4 with real sort -m."""
        files = []
        for i in range(4):
            p = tmp_path / ("f%d.txt" % i)
            p.write_text("r%d\tchr1\t%d\tchr1\t%d\t+\t-\n"
                         % (i, i * 100, i * 200))
            files.append(str(p))

        work = str(tmp_path / "work")
        os.makedirs(work)
        result = _phase4_initial_sort(files, 2, 1, 2, work)
        assert len(result) == 2


class TestPhase5:
    def test_any_genome(self, tmp_path, capsys):
        f = tmp_path / "batch.txt"
        f.write_text(
            "r1\tchr1\t100\tchr1\t200\t+\t-\n"
            "bad_line\n"
        )
        _phase5_validate([str(f)], "any", 1)
        captured = capsys.readouterr()
        assert "genome=any" in captured.err

    def test_removes_bad_lines(self, tmp_path, capsys):
        f = tmp_path / "batch.txt"
        f.write_text(
            "r1\tchr1\t100\tchr1\t200\t+\t-\n"
            "bad\n"
        )
        _phase5_validate([str(f)], "mm10", 1)
        captured = capsys.readouterr()
        assert "1 lines removed" in captured.err

    def test_all_valid(self, tmp_path, capsys):
        f = tmp_path / "batch.txt"
        f.write_text("r1\tchr1\t100\tchr1\t200\t+\t-\n")
        _phase5_validate([str(f)], "mm10", 1)
        captured = capsys.readouterr()
        assert "0 total lines removed" in captured.err


class TestPhase6:
    def test_single_batch(self, tmp_path, monkeypatch):
        """Single file: while loop not entered."""
        monkeypatch.setattr(
            "dip_c.commands.merge._sort_merge_batch",
            _fake_sort_merge_batch)

        f = tmp_path / "batch.txt"
        f.write_text("r1\tchr1\t100\tchr1\t200\t+\t-\n")
        work = str(tmp_path / "work")
        os.makedirs(work)
        result = _phase6_hierarchical_merge([str(f)], 500000, 1, 2, work)
        assert os.path.exists(result)

    def test_multi_level(self, tmp_path, monkeypatch):
        """Multiple files: hierarchical merge through multiple levels."""
        monkeypatch.setattr(
            "dip_c.commands.merge._sort_merge_batch",
            _fake_sort_merge_batch)

        files = []
        for i in range(6):
            p = tmp_path / ("batch_%d.txt" % i)
            p.write_text("r%d\tchr1\t%d\tchr1\t%d\t+\t-\n"
                         % (i, i * 100, i * 200))
            files.append(str(p))

        work = str(tmp_path / "work")
        os.makedirs(work)
        result = _phase6_hierarchical_merge(files, 500000, 1, 2, work)
        assert os.path.exists(result)
        with open(result) as fh:
            lines = fh.readlines()
        assert len(lines) == 6

    @needs_merge_tools
    def test_multi_level_real(self, tmp_path):
        """Integration: hierarchical merge with real sort -m."""
        files = []
        for i in range(6):
            p = tmp_path / ("batch_%d.txt" % i)
            p.write_text("r%d\tchr1\t%d\tchr1\t%d\t+\t-\n"
                         % (i, i * 100, i * 200))
            files.append(str(p))

        work = str(tmp_path / "work")
        os.makedirs(work)
        result = _phase6_hierarchical_merge(files, 500000, 1, 2, work)
        assert os.path.exists(result)


class TestPhase7:
    def test_creates_gz_mocked(self, tmp_path):
        """Test phase7 interface with gzip mock (no pigz)."""
        sorted_file = tmp_path / "sorted.txt"
        sorted_file.write_text("r1\tchr1\t100\tchr1\t200\t+\t-\n")
        out = str(tmp_path / "output.pairs.gz")
        _fake_phase7_compress(str(sorted_file), "## header\n", out, 1)
        assert os.path.exists(out)
        assert os.path.getsize(out) > 0

    @needs_merge_tools
    def test_creates_gz_real(self, tmp_path):
        """Integration: phase7 with real pigz."""
        sorted_file = tmp_path / "sorted.txt"
        sorted_file.write_text("r1\tchr1\t100\tchr1\t200\t+\t-\n")
        out = str(tmp_path / "output.pairs.gz")
        _phase7_compress(str(sorted_file), "## header\n", out, 1)
        assert os.path.exists(out)

    def test_pigz_failure(self, tmp_path, monkeypatch):
        """pigz returning non-zero exit code -> SystemExit(1)."""
        sorted_file = tmp_path / "sorted.txt"
        sorted_file.write_text("data\n")
        out = str(tmp_path / "output.pairs.gz")

        original_popen = subprocess.Popen

        class FakePigz:
            def __init__(self, *a, **kw):
                self.stdin = open(os.devnull, "wb")
                self.returncode = 1

            def wait(self):
                pass

        def fake_popen(cmd, **kw):
            if cmd[0] == "pigz":
                return FakePigz()
            return original_popen(cmd, **kw)

        monkeypatch.setattr(subprocess, "Popen", fake_popen)
        with pytest.raises(SystemExit) as exc_info:
            _phase7_compress(str(sorted_file), "## header\n", out, 1)
        assert exc_info.value.code == 1


# ===================================================================
# merge() main CLI entry point
# ===================================================================

class TestMergeCli:
    """Tests for the merge(argv) main function.

    Tests that don't need real tools mock check_required_tools and
    the sort/pigz-dependent phases.
    """

    def _mock_tools(self, monkeypatch):
        """Mock check_required_tools, sort, and pigz for local testing."""
        monkeypatch.setattr(
            "dip_c.commands.merge.check_required_tools", lambda: None)
        monkeypatch.setattr(
            "dip_c.commands.merge._sort_merge_batch",
            _fake_sort_merge_batch)
        monkeypatch.setattr(
            "dip_c.commands.merge._phase7_compress",
            _fake_phase7_compress)

    # -- Argument validation (no tools needed) -------------------------

    def test_bad_dir(self):
        with pytest.raises(SystemExit) as exc_info:
            merge(["merge", "-i", "/nonexistent/path", "-o", "out"])
        assert exc_info.value.code == 2

    def test_empty_dir(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        with pytest.raises(SystemExit) as exc_info:
            merge(["merge", "-i", str(d), "-o", "out"])
        assert exc_info.value.code == 2

    def test_unknown_genome(self, tmp_path):
        d = _make_pairs_dir(tmp_path)
        with pytest.raises(SystemExit) as exc_info:
            merge(["merge", "-i", d, "-o", "out", "-g", "hg99"])
        assert exc_info.value.code == 2

    # -- Single-file shortcut ------------------------------------------

    def test_auto_suffix(self, tmp_path, monkeypatch):
        """Output without .pairs.gz gets suffix appended."""
        self._mock_tools(monkeypatch)
        d = _make_pairs_dir(tmp_path, n_files=1)
        out = str(tmp_path / "result")
        merge(["merge", "-i", d, "-o", out])
        assert os.path.exists(out + ".pairs.gz")

    def test_single_file_copies(self, tmp_path, capsys, monkeypatch):
        self._mock_tools(monkeypatch)
        d = _make_pairs_dir(tmp_path, n_files=1)
        out = str(tmp_path / "single.pairs.gz")
        ret = merge(["merge", "-i", d, "-o", out])
        assert ret == 0
        assert os.path.exists(out)
        captured = capsys.readouterr()
        assert "Only 1 file found" in captured.err

    def test_no_genome_defaults_any(self, tmp_path, capsys, monkeypatch):
        self._mock_tools(monkeypatch)
        d = _make_pairs_dir(tmp_path, n_files=1)
        out = str(tmp_path / "out.pairs.gz")
        merge(["merge", "-i", d, "-o", out])
        captured = capsys.readouterr()
        assert "Genome: any" in captured.err

    def test_custom_header(self, tmp_path, monkeypatch):
        self._mock_tools(monkeypatch)
        d = _make_pairs_dir(tmp_path, n_files=1)
        hdr = tmp_path / "header.txt"
        hdr.write_text("## custom header\n#columns: a b c d e f g\n")
        out = str(tmp_path / "out.pairs.gz")
        ret = merge(["merge", "-i", d, "-o", out,
                     "-h", str(hdr)])
        assert ret == 0

    # -- Full pipeline (mocked sort/pigz) ------------------------------

    @needs_gunzip
    def test_full_pipeline_mocked(self, tmp_path, capsys, monkeypatch):
        """Full 7-phase pipeline with mocked sort and pigz."""
        self._mock_tools(monkeypatch)
        d = _make_pairs_dir(tmp_path, n_files=3)
        out = str(tmp_path / "merged.pairs.gz")
        ret = merge(["merge", "-i", d, "-o", out, "-g", "mm10",
                     "-j", "2", "-m", "4"])
        assert ret == 0
        assert os.path.exists(out)
        assert os.path.getsize(out) > 0
        captured = capsys.readouterr()
        assert "Merge completed" in captured.err

    @needs_gunzip
    def test_cleanup_on_error(self, tmp_path, monkeypatch):
        """Work directory is cleaned up even if a phase raises."""
        self._mock_tools(monkeypatch)
        d = _make_pairs_dir(tmp_path, n_files=2)
        out = str(tmp_path / "out.pairs.gz")

        def boom(*args, **kwargs):
            raise RuntimeError("boom")

        monkeypatch.setattr("dip_c.commands.merge._phase2_estimate", boom)
        with pytest.raises(RuntimeError):
            merge(["merge", "-i", d, "-o", out, "-j", "1", "-m", "4"])

        # Work dir should be cleaned up
        work_dirs = [f for f in os.listdir(d) if f.startswith(".merge_work_")]
        assert len(work_dirs) == 0

    @needs_gunzip
    def test_output_missing_returns_1(self, tmp_path, monkeypatch):
        """If output file is missing after pipeline, return 1."""
        self._mock_tools(monkeypatch)
        d = _make_pairs_dir(tmp_path, n_files=2)
        out = str(tmp_path / "out.pairs.gz")

        # Override _phase7_compress to be a true no-op (don't create file)
        monkeypatch.setattr(
            "dip_c.commands.merge._phase7_compress",
            lambda *a, **kw: None)
        ret = merge(["merge", "-i", d, "-o", out, "-j", "1", "-m", "4"])
        assert ret == 1

    @needs_gunzip
    def test_auto_extracted_header(self, tmp_path, capsys, monkeypatch):
        """Default path: header auto-extracted from first .pairs.gz file."""
        self._mock_tools(monkeypatch)
        d = _make_pairs_dir(tmp_path, n_files=2)
        out = str(tmp_path / "out.pairs.gz")
        ret = merge(["merge", "-i", d, "-o", out, "-j", "1", "-m", "4"])
        assert ret == 0
        captured = capsys.readouterr()
        assert "Header auto-extracted" in captured.err

    @needs_gunzip
    def test_custom_header_multi_file(self, tmp_path, capsys, monkeypatch):
        """Custom header with multi-file pipeline."""
        self._mock_tools(monkeypatch)
        d = _make_pairs_dir(tmp_path, n_files=2)
        hdr = tmp_path / "header.txt"
        hdr.write_text("## custom\n")
        out = str(tmp_path / "out.pairs.gz")
        ret = merge(["merge", "-i", d, "-o", out,
                     "-h", str(hdr), "-j", "1", "-m", "4"])
        assert ret == 0
        captured = capsys.readouterr()
        assert "Header loaded from" in captured.err

    # -- Integration tests (real tools) --------------------------------

    @needs_merge_tools
    @needs_gunzip
    def test_full_pipeline_real(self, tmp_path, capsys):
        """Full pipeline with real sort -m and pigz."""
        d = _make_pairs_dir(tmp_path, n_files=3)
        out = str(tmp_path / "merged.pairs.gz")
        ret = merge(["merge", "-i", d, "-o", out, "-g", "mm10",
                     "-j", "2", "-m", "4"])
        assert ret == 0
        assert os.path.exists(out)
        captured = capsys.readouterr()
        assert "Merge completed" in captured.err
