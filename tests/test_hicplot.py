"""Tests for dip_c.commands.hicplot."""

import argparse
import os
import sys

import numpy as np
import pytest

from dip_c.commands.hicplot import (
    _build_parser,
    _diff,
    _parse_file_spec,
    _parse_region,
    _plot_all,
    _plot_allinone,
    _plot_diff_all,
    _plot_diff_allinone,
    _plot_diff_map,
    _plot_map,
    hicplot,
)
from dip_c.hicplot_utils import make_colormap, REDMAP_SPEC, BWRMAP_SPEC


# ===================================================================
# _parse_file_spec
# ===================================================================

class TestParseFileSpec:
    def test_int_norm(self):
        assert _parse_file_spec("abc.hic:123456789") == ("abc.hic", 123456789)

    def test_float_norm(self):
        assert _parse_file_spec("abc.hic:1.5") == ("abc.hic", 1.5)

    def test_no_norm(self):
        assert _parse_file_spec("abc.hic") == ("abc.hic", 1)

    def test_bad_norm_treated_as_path(self):
        assert _parse_file_spec("abc.hic:notanumber") == ("abc.hic:notanumber", 1)

    def test_colon_at_start(self):
        # idx == 0 => not a valid split point
        assert _parse_file_spec(":123") == (":123", 1)

    def test_path_with_directory(self):
        assert _parse_file_spec("/tmp/data/file.hic:100") == (
            "/tmp/data/file.hic", 100)


# ===================================================================
# _parse_region
# ===================================================================

class TestParseRegion:
    def test_valid(self):
        assert _parse_region("chr1:1-10000000") == ("chr1", 1, 10000000)

    def test_no_colon(self):
        with pytest.raises(argparse.ArgumentTypeError, match="CHR:START-END"):
            _parse_region("chr1_1_10000000")

    def test_multi_dash(self):
        with pytest.raises(argparse.ArgumentTypeError, match="START-END"):
            _parse_region("chr1:1-2-3")

    def test_non_integer(self):
        with pytest.raises(argparse.ArgumentTypeError, match="integers"):
            _parse_region("chr1:abc-def")


# ===================================================================
# _diff
# ===================================================================

class TestDiff:
    def test_none_normalization(self):
        m1 = np.array([[10.0, 5.0], [5.0, 20.0]])
        m2 = np.array([[8.0, 3.0], [3.0, 15.0]])
        result = _diff(m1, 2, m2, 4, "NONE")
        expected = m1 / 2 - m2 / 4
        np.testing.assert_array_almost_equal(result, expected)

    def test_non_none_normalization(self):
        m1 = np.array([[10.0, 5.0]])
        m2 = np.array([[8.0, 3.0]])
        result = _diff(m1, 999, m2, 999, "KR")
        expected = m1 - m2
        np.testing.assert_array_equal(result, expected)


# ===================================================================
# Plotting functions — integration tests using tiny_cool fixtures
# ===================================================================

class TestPlotMap:
    def test_normal_order(self, tiny_cool, tmp_path):
        out = str(tmp_path / "map.png")
        _plot_map(tiny_cool, 1000, out,
                  ["chr1", "chr1"], [0, 500000, 0, 500000],
                  100_000, 200, "NONE")
        assert os.path.exists(out)
        assert os.path.getsize(out) > 0

    def test_reversed_order(self, tiny_cool, tmp_path):
        out = str(tmp_path / "map_rev.png")
        _plot_map(tiny_cool, 1000, out,
                  ["chr2", "chr1"], [0, 300000, 0, 500000],
                  100_000, 200, "NONE")
        assert os.path.exists(out)

    def test_with_normalization(self, tiny_cool, tmp_path):
        out = str(tmp_path / "map_scale.png")
        _plot_map(tiny_cool, 1, out,
                  ["chr1", "chr1"], [0, 500000, 0, 500000],
                  100_000, 500, "SCALE")
        assert os.path.exists(out)

    def test_custom_cmap(self, tiny_cool, tmp_path):
        out = str(tmp_path / "map_cmap.png")
        cmap = make_colormap(("custom", [(0, 0, 0), (1, 1, 0)]))
        _plot_map(tiny_cool, 1000, out,
                  ["chr1", "chr1"], [0, 500000, 0, 500000],
                  100_000, 200, "NONE", cmap=cmap)
        assert os.path.exists(out)


class TestPlotAll:
    def test_creates_per_chrom_pngs(self, tiny_cool, tmp_path):
        out = str(tmp_path / "all.png")
        _plot_all(tiny_cool, 1000, out, 100_000, 200, "NONE")
        assert os.path.exists(str(tmp_path / "all_chr1.png"))
        assert os.path.exists(str(tmp_path / "all_chr2.png"))

    def test_with_normalization(self, tiny_cool, tmp_path):
        out = str(tmp_path / "all_scale.png")
        _plot_all(tiny_cool, 1, out, 100_000, 500, "SCALE")
        assert os.path.exists(str(tmp_path / "all_scale_chr1.png"))

    def test_skipped_chroms_warning(self, partial_cool, tmp_path, capsys):
        out = str(tmp_path / "all_partial.png")
        _plot_all(partial_cool, 1000, out, 100_000, 200, "NONE")
        captured = capsys.readouterr()
        assert "Skipped chromosomes" in captured.err


class TestPlotAllinone:
    def test_with_gridlines(self, tiny_cool, tmp_path):
        out = str(tmp_path / "aio.png")
        _plot_allinone(tiny_cool, 1000, out, 100_000, 200, "NONE",
                       gridlines=True)
        assert os.path.exists(out)
        assert os.path.getsize(out) > 0

    def test_no_gridlines(self, tiny_cool, tmp_path):
        out = str(tmp_path / "aio_ng.png")
        _plot_allinone(tiny_cool, 1000, out, 100_000, 200, "NONE",
                       gridlines=False)
        assert os.path.exists(out)

    def test_skipped_chroms_warning(self, partial_cool, tmp_path, capsys):
        out = str(tmp_path / "aio_partial.png")
        _plot_allinone(partial_cool, 1000, out, 100_000, 200, "NONE")
        captured = capsys.readouterr()
        assert "Skipped chromosomes" in captured.err


class TestPlotDiffMap:
    def test_normal(self, tiny_cool, tiny_cool_b, tmp_path):
        out = str(tmp_path / "diff.png")
        _plot_diff_map(tiny_cool, 1000, tiny_cool_b, 1000, out,
                       ["chr1", "chr1"], [0, 500000, 0, 500000],
                       100_000, 200, "NONE")
        assert os.path.exists(out)

    def test_reversed_order(self, tiny_cool, tiny_cool_b, tmp_path):
        out = str(tmp_path / "diff_rev.png")
        _plot_diff_map(tiny_cool, 1000, tiny_cool_b, 1000, out,
                       ["chr2", "chr1"], [0, 300000, 0, 500000],
                       100_000, 200, "NONE")
        assert os.path.exists(out)


class TestPlotDiffAll:
    def test_creates_per_chrom_pngs(self, tiny_cool, tiny_cool_b, tmp_path):
        out = str(tmp_path / "diffall.png")
        _plot_diff_all(tiny_cool, 1000, tiny_cool_b, 1000, out,
                       100_000, 200, "NONE")
        assert os.path.exists(str(tmp_path / "diffall_chr1.png"))
        assert os.path.exists(str(tmp_path / "diffall_chr2.png"))

    def test_skipped_chroms_warning(self, partial_cool, tmp_path, capsys):
        out = str(tmp_path / "diffall_partial.png")
        _plot_diff_all(partial_cool, 1000, partial_cool, 1000, out,
                       100_000, 200, "NONE")
        captured = capsys.readouterr()
        assert "Skipped chromosomes" in captured.err


class TestPlotDiffAllinone:
    def test_with_gridlines(self, tiny_cool, tiny_cool_b, tmp_path):
        out = str(tmp_path / "daio.png")
        _plot_diff_allinone(tiny_cool, 1000, tiny_cool_b, 1000, out,
                            100_000, 200, "NONE", gridlines=True)
        assert os.path.exists(out)

    def test_no_gridlines(self, tiny_cool, tiny_cool_b, tmp_path):
        out = str(tmp_path / "daio_ng.png")
        _plot_diff_allinone(tiny_cool, 1000, tiny_cool_b, 1000, out,
                            100_000, 200, "NONE", gridlines=False)
        assert os.path.exists(out)

    def test_skipped_chroms_warning(self, partial_cool, tmp_path, capsys):
        out = str(tmp_path / "daio_partial.png")
        _plot_diff_allinone(partial_cool, 1000, partial_cool, 1000, out,
                            100_000, 200, "NONE")
        captured = capsys.readouterr()
        assert "Skipped chromosomes" in captured.err


# ===================================================================
# hicplot() main CLI entry point
# ===================================================================

class TestHicplotCli:
    """Test the main hicplot(argv) dispatcher — all argument branches."""

    def test_absolute_map_mode(self, tiny_cool, tmp_path):
        out = str(tmp_path / "out.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200",
                       "-c1", "chr1:0-500000"])
        assert ret == 0
        assert os.path.exists(out)

    def test_absolute_all_mode(self, tiny_cool, tmp_path):
        out = str(tmp_path / "out.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200"])
        assert ret == 0
        assert os.path.exists(str(tmp_path / "out_chr1.png"))

    def test_absolute_allinone_mode(self, tiny_cool, tmp_path):
        out = str(tmp_path / "out.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200",
                       "--allinone"])
        assert ret == 0
        assert os.path.exists(out)

    def test_absolute_allinone_no_gridlines(self, tiny_cool, tmp_path, capsys):
        out = str(tmp_path / "out.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200",
                       "--allinone", "-G"])
        assert ret == 0
        captured = capsys.readouterr()
        assert "Gridlines: off" in captured.err

    def test_asymmetric_region(self, tiny_cool, tmp_path):
        out = str(tmp_path / "out.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200",
                       "-c1", "chr1:0-500000",
                       "-c2", "chr2:0-300000"])
        assert ret == 0

    def test_diff_map_mode(self, tiny_cool, tiny_cool_b, tmp_path):
        out = str(tmp_path / "diff.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-2", "%s:1000" % tiny_cool_b,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200",
                       "-c1", "chr1:0-500000"])
        assert ret == 0
        assert os.path.exists(out)

    def test_diff_all_mode(self, tiny_cool, tiny_cool_b, tmp_path):
        out = str(tmp_path / "diff.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-2", "%s:1000" % tiny_cool_b,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200"])
        assert ret == 0

    def test_diff_allinone_mode(self, tiny_cool, tiny_cool_b, tmp_path):
        out = str(tmp_path / "diff.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-2", "%s:1000" % tiny_cool_b,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200",
                       "--allinone"])
        assert ret == 0
        assert os.path.exists(out)

    def test_custom_colormap(self, tiny_cool, tmp_path):
        out = str(tmp_path / "out.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200",
                       "-c1", "chr1:0-500000",
                       "-C", "viridis"])
        assert ret == 0

    def test_with_scale_normalization(self, tiny_cool, tmp_path):
        out = str(tmp_path / "out.png")
        ret = hicplot(["hicplot",
                       "-1", tiny_cool,
                       "-o", out,
                       "-r", "100000",
                       "-s", "500",
                       "-n", "SCALE",
                       "-c1", "chr1:0-500000"])
        assert ret == 0

    # -- Error branches ------------------------------------------------

    def test_none_norm_no_factor_errors(self, tiny_cool, tmp_path):
        out = str(tmp_path / "out.png")
        with pytest.raises(SystemExit) as exc_info:
            hicplot(["hicplot",
                     "-1", tiny_cool,      # no :NORM
                     "-o", out,
                     "-r", "100000",
                     "-s", "200"])
        assert exc_info.value.code == 2

    def test_diff_none_norm_no_factor2_errors(self, tiny_cool, tmp_path):
        out = str(tmp_path / "out.png")
        with pytest.raises(SystemExit) as exc_info:
            hicplot(["hicplot",
                     "-1", "%s:1000" % tiny_cool,
                     "-2", tiny_cool,      # no :NORM on -2
                     "-o", out,
                     "-r", "100000",
                     "-s", "200"])
        assert exc_info.value.code == 2

    def test_c2_without_c1_errors(self, tiny_cool, tmp_path):
        out = str(tmp_path / "out.png")
        with pytest.raises(SystemExit) as exc_info:
            hicplot(["hicplot",
                     "-1", "%s:1000" % tiny_cool,
                     "-o", out,
                     "-r", "100000",
                     "-s", "200",
                     "-c2", "chr1:0-500000"])
        assert exc_info.value.code == 2

    def test_allinone_with_region_warns(self, tiny_cool, tmp_path, capsys):
        out = str(tmp_path / "out.png")
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200",
                       "-c1", "chr1:0-500000",
                       "--allinone"])
        assert ret == 0
        captured = capsys.readouterr()
        assert "--allinone ignored" in captured.err

    def test_memory_error(self, tiny_cool, tmp_path, monkeypatch):
        out = str(tmp_path / "out.png")

        def raise_memory(*args, **kwargs):
            raise MemoryError("out of memory")

        monkeypatch.setattr(
            "dip_c.commands.hicplot._plot_map", raise_memory)
        ret = hicplot(["hicplot",
                       "-1", "%s:1000" % tiny_cool,
                       "-o", out,
                       "-r", "100000",
                       "-s", "200",
                       "-c1", "chr1:0-500000"])
        assert ret == 1

    def test_allinone_gridlines_on_message(self, tiny_cool, tmp_path, capsys):
        out = str(tmp_path / "out.png")
        hicplot(["hicplot",
                 "-1", "%s:1000" % tiny_cool,
                 "-o", out,
                 "-r", "100000",
                 "-s", "200",
                 "--allinone"])
        captured = capsys.readouterr()
        assert "Gridlines: on" in captured.err
