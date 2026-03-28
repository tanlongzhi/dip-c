"""Unit tests for dip_c.hicplot_utils."""

import importlib
import os
import sys

import numpy as np
import pytest

import dip_c.hicplot_utils as hu
from dip_c.hicplot_utils import (
    BWRMAP_SPEC,
    DEFAULT_SKIP_CHROMS,
    REDMAP_SPEC,
    _ensure_hicplot_deps,
    _suppress_native_stderr,
    chroms_per_row_for_genome,
    filter_chroms,
    get_chrom_names,
    get_positions_and_matrix_size,
    make_colormap,
    ordering_check,
    plot_matrix,
    resolve_colormap,
    strip_png_ext,
    validate_chroms,
)


# ===================================================================
# _ensure_hicplot_deps
# ===================================================================

class TestEnsureHicplotDeps:
    def test_imports_successfully(self):
        _ensure_hicplot_deps()
        assert hu._HICTKPY is not None
        assert hu._PLT is not None
        assert hu._LSCM is not None

    def test_noop_on_second_call(self):
        _ensure_hicplot_deps()
        hictkpy_ref = hu._HICTKPY
        _ensure_hicplot_deps()
        assert hu._HICTKPY is hictkpy_ref

    def test_missing_deps_raises_system_exit(self, monkeypatch):
        import builtins
        real_import = builtins.__import__

        # Save and reset the global so the guard triggers
        saved = hu._HICTKPY
        hu._HICTKPY = None

        def fake_import(name, *args, **kwargs):
            if name == "hictkpy":
                raise ImportError("no hictkpy")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        try:
            with pytest.raises(SystemExit) as exc_info:
                _ensure_hicplot_deps()
            assert exc_info.value.code == 1
        finally:
            hu._HICTKPY = saved


# ===================================================================
# get_chrom_names
# ===================================================================

class TestGetChromNames:
    def test_returns_chrom_list(self, tiny_cool):
        import hictkpy
        hic = hictkpy.File(tiny_cool, 100_000)
        names = get_chrom_names(hic)
        assert names == ["chr1", "chr2"]


# ===================================================================
# ordering_check
# ===================================================================

class TestOrderingCheck:
    def test_same_chrom(self):
        assert ordering_check("chr1", "chr1", ["chr1", "chr2"]) is True

    def test_less(self):
        assert ordering_check("chr1", "chr2", ["chr1", "chr2"]) is True

    def test_greater(self):
        assert ordering_check("chr2", "chr1", ["chr1", "chr2"]) is False

    def test_unknown_chrom_returns_true(self):
        assert ordering_check("chrZ", "chr1", ["chr1", "chr2"]) is True


# ===================================================================
# make_colormap
# ===================================================================

class TestMakeColormap:
    def test_redmap(self):
        from matplotlib.colors import LinearSegmentedColormap
        cmap = make_colormap(REDMAP_SPEC)
        assert isinstance(cmap, LinearSegmentedColormap)

    def test_bwrmap(self):
        from matplotlib.colors import LinearSegmentedColormap
        cmap = make_colormap(BWRMAP_SPEC)
        assert isinstance(cmap, LinearSegmentedColormap)


# ===================================================================
# resolve_colormap
# ===================================================================

class TestResolveColormap:
    def test_named_matplotlib_cmap(self):
        cmap = resolve_colormap("viridis")
        assert cmap is not None

    def test_comma_separated_colors(self):
        cmap = resolve_colormap("white,red")
        assert cmap is not None

    def test_three_color_stops(self):
        cmap = resolve_colormap("blue,white,red")
        assert cmap is not None

    def test_hex_colors(self):
        cmap = resolve_colormap("#0000ff,#ffffff,#ff0000")
        assert cmap is not None

    def test_single_stop_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            resolve_colormap("justred")
        assert exc_info.value.code == 1

    def test_bad_color_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            resolve_colormap("white,notacolor")
        assert exc_info.value.code == 1


# ===================================================================
# plot_matrix
# ===================================================================

class TestPlotMatrix:
    def test_small_matrix(self, tmp_path):
        mat = np.array([[100, 50], [50, 200]], dtype=float)
        out = str(tmp_path / "small.png")
        cmap = make_colormap(REDMAP_SPEC)
        plot_matrix(mat, out, cmap, 0, 200)

        assert os.path.exists(out)
        assert os.path.getsize(out) > 0

        import matplotlib.image as mpimg
        data = mpimg.imread(out)
        assert data.ndim == 3
        assert data.shape[2] == 4  # RGBA

    def test_large_matrix_auto_1px(self, tmp_path):
        mat = np.zeros((1000, 1000))
        out = str(tmp_path / "large.png")
        cmap = make_colormap(REDMAP_SPEC)
        plot_matrix(mat, out, cmap, 0, 1)

        assert os.path.exists(out)
        assert os.path.getsize(out) > 0

    def test_with_gridlines(self, tmp_path):
        mat = np.random.default_rng(42).random((10, 10))
        out = str(tmp_path / "grid.png")
        cmap = make_colormap(REDMAP_SPEC)
        plot_matrix(mat, out, cmap, 0, 1, hlines=[3, 6], vlines=[3, 6])

        assert os.path.exists(out)
        assert os.path.getsize(out) > 0

    def test_pixel_spot_check_red_white(self, tmp_path):
        """Verify high-count bins are red, zero bins are white."""
        mat = np.array([[200, 0], [0, 0]], dtype=float)
        out = str(tmp_path / "spot.png")
        cmap = make_colormap(REDMAP_SPEC)
        plot_matrix(mat, out, cmap, 0, 200)

        import matplotlib.image as mpimg
        data = mpimg.imread(out)
        h, w = data.shape[:2]
        # Top-left quadrant should be reddish (R high, G low, B low)
        tl = data[h // 4, w // 4]
        assert tl[0] > 0.7, "Expected red channel high for max-value bin"
        assert tl[1] < 0.4, "Expected green channel low for max-value bin"
        # Bottom-right quadrant should be whitish (all channels high)
        br = data[3 * h // 4, 3 * w // 4]
        assert br[0] > 0.8, "Expected R high for zero-value bin"
        assert br[1] > 0.8, "Expected G high for zero-value bin"


# ===================================================================
# get_positions_and_matrix_size
# ===================================================================

class TestGetPositionsAndMatrixSize:
    def test_single_matrix(self):
        matrices = [np.zeros((5, 5))]
        positions, size = get_positions_and_matrix_size(matrices, 1)
        assert positions == [(0, 0)]
        assert size == (5, 5)

    def test_2x2_grid(self):
        matrices = [np.zeros((5, 5)), np.zeros((5, 3)),
                     np.zeros((3, 5)), np.zeros((3, 3))]
        positions, size = get_positions_and_matrix_size(matrices, 2)
        assert len(positions) == 4
        assert positions[0] == (0, 0)
        assert positions[1] == (0, 5)
        assert positions[2] == (5, 0)
        assert positions[3] == (5, 5)
        assert size == (8, 8)

    def test_uneven_grid(self):
        matrices = [np.zeros((4, 4)), np.zeros((4, 3)),
                     np.zeros((3, 4))]
        positions, size = get_positions_and_matrix_size(matrices, 2)
        assert len(positions) == 3
        # First row: (0,0) and (0,4)
        assert positions[0] == (0, 0)
        assert positions[1] == (0, 4)
        # Second row: (4,0)
        assert positions[2] == (4, 0)


# ===================================================================
# chroms_per_row_for_genome
# ===================================================================

class TestChromsPerRowForGenome:
    def test_returns_input(self):
        assert chroms_per_row_for_genome(5) == 5
        assert chroms_per_row_for_genome(22) == 22


# ===================================================================
# filter_chroms
# ===================================================================

class TestFilterChroms:
    def test_removes_all(self):
        result = filter_chroms([("All", 1000), ("chr1", 500000)])
        assert result == [("chr1", 500000)]

    def test_removes_ALL(self):
        result = filter_chroms([("ALL", 1000), ("chr1", 500000)])
        assert result == [("chr1", 500000)]

    def test_passes_through(self):
        chrom_list = [("chr1", 500000), ("chr2", 300000)]
        assert filter_chroms(chrom_list) == chrom_list


# ===================================================================
# _suppress_native_stderr
# ===================================================================

class TestSuppressNativeStderr:
    def test_restores_fd(self, capsys):
        sys.stderr.write("before\n")
        with _suppress_native_stderr():
            # Writing inside context — should be suppressed at fd level
            os.write(2, b"suppressed\n")
        sys.stderr.write("after\n")
        captured = capsys.readouterr()
        assert "before" in captured.err
        assert "after" in captured.err
        assert "suppressed" not in captured.err

    def test_restores_on_error(self):
        try:
            with _suppress_native_stderr():
                raise ValueError("test")
        except ValueError:
            pass
        # stderr should still work after exception
        sys.stderr.write("works\n")
        sys.stderr.flush()


# ===================================================================
# validate_chroms
# ===================================================================

class TestValidateChroms:
    def test_keeps_valid(self, tiny_cool):
        import hictkpy
        hic = hictkpy.File(tiny_cool, 100_000)
        chroms = [("chr1", 500000), ("chr2", 300000)]
        valid, skipped = validate_chroms(hic, chroms, "NONE", 100_000)
        assert len(valid) == 2
        assert len(skipped) == 0

    def test_skips_empty(self, empty_cool):
        import hictkpy
        hic = hictkpy.File(empty_cool, 100_000)
        chroms = [("chr1", 500000), ("chr2", 300000)]
        valid, skipped = validate_chroms(hic, chroms, "NONE", 100_000)
        assert len(valid) == 0
        assert len(skipped) == 2

    def test_with_normalization(self, tiny_cool):
        import hictkpy
        hic = hictkpy.File(tiny_cool, 100_000)
        chroms = [("chr1", 500000)]
        valid, skipped = validate_chroms(hic, chroms, "SCALE", 100_000)
        assert len(valid) == 1

    def test_exception_skips(self):
        """A chromosome whose fetch raises is skipped."""
        class FakeHic:
            def fetch(self, *args, **kwargs):
                raise RuntimeError("test error")

        chroms = [("chr1", 500000)]
        valid, skipped = validate_chroms(FakeHic(), chroms, "NONE", 100_000)
        assert len(valid) == 0
        assert skipped == ["chr1"]


# ===================================================================
# strip_png_ext
# ===================================================================

class TestStripPngExt:
    def test_removes_extension(self):
        assert strip_png_ext("foo.png") == "foo"

    def test_noop_on_other_ext(self):
        assert strip_png_ext("foo.jpg") == "foo.jpg"

    def test_noop_on_no_ext(self):
        assert strip_png_ext("foo") == "foo"

    def test_only_removes_trailing(self):
        assert strip_png_ext("foo.png.bak") == "foo.png.bak"
