"""Tests for the pairs2con command."""

import gzip
import os
import tempfile

import pytest

from dip_c.commands.pairs2con import pairs2con


def _write_pairs_gz(path, lines):
    """Write lines to a gzipped file."""
    with gzip.open(path, "wt") as f:
        for line in lines:
            f.write(line + "\n")


def _read_con_gz(path):
    """Read lines from a gzipped con file."""
    with gzip.open(path, "rt") as f:
        return [line.rstrip("\n") for line in f]


class TestPairs2con:
    def test_usage_with_no_args(self, capsys):
        ret = pairs2con(["pairs2con"])
        assert ret == 1
        err = capsys.readouterr().err
        assert "Usage:" in err

    def test_rejects_non_pairs_gz(self, capsys):
        ret = pairs2con(["pairs2con", "foo.txt"])
        assert ret == 1
        err = capsys.readouterr().err
        assert ".pairs.gz" in err

    def test_basic_conversion(self, tmp_path):
        input_file = str(tmp_path / "test.pairs.gz")
        _write_pairs_gz(input_file, [
            "#comment",
            "read1\tchr1\t100\tchr2\t200\t+\t-\t0\t1",
        ])
        ret = pairs2con(["pairs2con", input_file])
        assert ret == 0
        output_file = str(tmp_path / "test.con.gz")
        assert os.path.exists(output_file)
        lines = _read_con_gz(output_file)
        assert lines == ["chr1,100,0\tchr2,200,1"]

    def test_skips_comments(self, tmp_path):
        input_file = str(tmp_path / "test.pairs.gz")
        _write_pairs_gz(input_file, [
            "#header1",
            "#header2",
            "read1\tchr1\t100\tchr2\t200\t+\t-\t0\t1",
        ])
        ret = pairs2con(["pairs2con", input_file])
        assert ret == 0
        lines = _read_con_gz(str(tmp_path / "test.con.gz"))
        assert len(lines) == 1

    def test_missing_phases_become_dots(self, tmp_path):
        input_file = str(tmp_path / "test.pairs.gz")
        _write_pairs_gz(input_file, [
            "read1\tchr3\t300\tchr4\t400\t+\t+\t\t",
        ])
        ret = pairs2con(["pairs2con", input_file])
        assert ret == 0
        lines = _read_con_gz(str(tmp_path / "test.con.gz"))
        assert lines == ["chr3,300,.\tchr4,400,."]

    def test_absent_phase_columns(self, tmp_path):
        input_file = str(tmp_path / "test.pairs.gz")
        _write_pairs_gz(input_file, [
            "read1\tchr5\t500\tchr6\t600\t-\t+",
        ])
        ret = pairs2con(["pairs2con", input_file])
        assert ret == 0
        lines = _read_con_gz(str(tmp_path / "test.con.gz"))
        assert lines == ["chr5,500,.\tchr6,600,."]

    def test_multiple_records(self, tmp_path):
        input_file = str(tmp_path / "test.pairs.gz")
        _write_pairs_gz(input_file, [
            "#header",
            "r1\tchr1\t10\tchr2\t20\t+\t-\t0\t1",
            "r2\tchr3\t30\tchr4\t40\t+\t+\t.\t.",
            "r3\tchr5\t50\tchr6\t60\t-\t+",
        ])
        ret = pairs2con(["pairs2con", input_file])
        assert ret == 0
        lines = _read_con_gz(str(tmp_path / "test.con.gz"))
        assert len(lines) == 3
        assert lines[0] == "chr1,10,0\tchr2,20,1"
        assert lines[1] == "chr3,30,.\tchr4,40,."
        assert lines[2] == "chr5,50,.\tchr6,60,."

    def test_output_filename(self, tmp_path):
        input_file = str(tmp_path / "contacts_unisex.pairs.gz")
        _write_pairs_gz(input_file, [
            "r1\tchr1\t10\tchr2\t20\t+\t-\t0\t1",
        ])
        pairs2con(["pairs2con", input_file])
        expected = str(tmp_path / "contacts_unisex.con.gz")
        assert os.path.exists(expected)
