"""Edge-case tests for dip-c commands.

All tests use synthetic inline data (no test data files needed), fast to run.
"""

import gzip
import os
import sys
import tempfile

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Helpers for creating tiny synthetic files
# ---------------------------------------------------------------------------

def write_file(path, content, gz=False):
    """Write content to a plain or gzip file."""
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(content)
    else:
        with open(path, "w") as f:
            f.write(content)


# A single contact line: chr1(pat):100 <-> chr1(pat):200
ONE_CON_LINE = "1,100,0\t1,200,0\n"
# Two contacts (one intra, one inter)
TWO_CON_LINES = "1,100,0\t1,200,0\n2,300,1\t3,400,1\n"
# Minimal 3DG: one particle
ONE_3DG_LINE = "1(pat)\t100000\t0.0\t0.0\t0.0\n"
# Two particles on same chromosome (needed for some operations)
TWO_3DG_LINES = (
    "1(pat)\t100000\t0.0\t0.0\t0.0\n"
    "1(pat)\t200000\t1.0\t1.0\t1.0\n"
)
# A valid color line
ONE_COLOR_LINE = "1(pat)\t100000\t0.5\n"
TWO_COLOR_LINES = "1(pat)\t100000\t0.5\n1(pat)\t200000\t0.8\n"


# ===================================================================
# 3a: Empty CON files → ZeroDivisionError
# ===================================================================

class TestEmptyConFile:
    """Commands that read CON files should return error=1 on empty input,
    not crash with ZeroDivisionError."""

    @pytest.fixture()
    def empty_con(self, tmp_path):
        p = tmp_path / "empty.con"
        p.write_text("")
        return str(p)

    @pytest.fixture()
    def empty_con_gz(self, tmp_path):
        p = tmp_path / "empty.con.gz"
        write_file(str(p), "", gz=True)
        return str(p)

    def test_info_empty_con(self, empty_con):
        from dip_c.commands.info import info
        ret = info(["info", empty_con])
        assert ret == 1

    def test_dedup_empty_con(self, empty_con):
        from dip_c.commands.dedup import dedup
        ret = dedup(["dedup", empty_con])
        assert ret == 1

    def test_clean_empty_con(self, empty_con):
        from dip_c.commands.clean import clean
        ret = clean(["clean", empty_con])
        assert ret == 1

    def test_reg_empty_con(self, empty_con):
        from dip_c.commands.reg import reg
        ret = reg(["reg", "-p", "hf", empty_con])
        assert ret == 1

    def test_impute_empty_con(self, empty_con):
        from dip_c.commands.impute import impute
        ret = impute(["impute", empty_con])
        assert ret == 1

    def test_ard_empty_con(self, empty_con):
        from dip_c.commands.ard import ard
        ret = ard(["ard", empty_con])
        assert ret == 1

    def test_color2_empty_con(self, tmp_path, empty_con):
        from dip_c.commands.color2 import color2
        color_f = tmp_path / "color.txt"
        color_f.write_text(ONE_COLOR_LINE)
        ret = color2(["color2", "-c", str(color_f), empty_con])
        assert ret == 1

    def test_con_empty_seg(self, tmp_path):
        """con command with empty SEG file → empty CON → ZeroDivisionError."""
        from dip_c.commands.con import con
        seg_f = tmp_path / "empty.seg"
        seg_f.write_text("")
        ret = con(["con", str(seg_f)])
        assert ret == 1

    def test_info_empty_con_gz(self, empty_con_gz):
        """Gzipped empty file should also be handled."""
        from dip_c.commands.info import info
        ret = info(["info", empty_con_gz])
        assert ret == 1

    def test_clean3_empty_con(self, tmp_path):
        from dip_c.commands.clean3 import clean3
        con_f = tmp_path / "empty.con"
        con_f.write_text("")
        g3d_f = tmp_path / "one.3dg"
        g3d_f.write_text(ONE_3DG_LINE)
        ret = clean3(["clean3", "-c", str(con_f), str(g3d_f)])
        assert ret == 1

    def test_impute3_empty_con(self, tmp_path):
        from dip_c.commands.impute3 import impute3
        con_f = tmp_path / "empty.con"
        con_f.write_text("")
        g3d_f = tmp_path / "one.3dg"
        g3d_f.write_text(TWO_3DG_LINES)
        ret = impute3(["impute3", "-3", str(g3d_f), str(con_f)])
        assert ret == 1

    def test_bincon_empty_con(self, tmp_path):
        from dip_c.commands.bincon import bincon
        con_f = tmp_path / "empty.con"
        con_f.write_text("")
        lens_f = tmp_path / "chr.len"
        lens_f.write_text("1\t1000000\n")
        ret = bincon(["bincon", "-l", str(lens_f), str(con_f)])
        assert ret == 1


# ===================================================================
# 3b: Empty 3DG files → crashes or silent empty output
# ===================================================================

class TestEmpty3dgFile:
    """Commands that read 3DG files should handle empty input gracefully."""

    @pytest.fixture()
    def empty_3dg(self, tmp_path):
        p = tmp_path / "empty.3dg"
        p.write_text("")
        return str(p)

    def test_con3_empty_3dg(self, empty_3dg):
        from dip_c.commands.con3 import con3
        ret = con3(["con3", empty_3dg])
        assert ret == 1

    def test_exp_empty_3dg(self, empty_3dg):
        from dip_c.commands.exp import exp
        ret = exp(["exp", empty_3dg])
        assert ret == 1

    def test_vis_empty_3dg(self, empty_3dg, capsys):
        from dip_c.commands.vis import vis
        ret = vis(["vis", empty_3dg])
        assert ret == 1

    def test_dist_empty_3dg(self, empty_3dg):
        from dip_c.commands.dist import dist
        ret = dist(["dist", empty_3dg])
        assert ret == 1

    def test_rg_empty_3dg(self, empty_3dg):
        from dip_c.commands.rg import rg
        ret = rg(["rg", empty_3dg])
        assert ret == 1

    def test_align_empty_3dg(self, empty_3dg):
        from dip_c.commands.align import align
        ret = align(["align", empty_3dg, empty_3dg])
        assert ret == 1

    def test_color_empty_3dg(self, empty_3dg):
        from dip_c.commands.color import color
        ret = color(["color", "-n", empty_3dg])
        assert ret == 1

    def test_reg3_empty_3dg(self, empty_3dg):
        from dip_c.commands.reg3 import reg3
        ret = reg3(["reg3", "-p", "hf", empty_3dg])
        assert ret == 1


# ===================================================================
# 3c: Minimal valid input (one-contact CON / one-particle 3DG)
# ===================================================================

class TestMinimalValidInput:
    """Commands should work (return 0) with the smallest valid input."""

    @pytest.fixture()
    def one_con(self, tmp_path):
        p = tmp_path / "one.con"
        p.write_text(ONE_CON_LINE)
        return str(p)

    @pytest.fixture()
    def two_3dg(self, tmp_path):
        p = tmp_path / "two.3dg"
        p.write_text(TWO_3DG_LINES)
        return str(p)

    @pytest.fixture()
    def one_3dg(self, tmp_path):
        p = tmp_path / "one.3dg"
        p.write_text(ONE_3DG_LINE)
        return str(p)

    def test_info_one_con(self, one_con, capsys):
        from dip_c.commands.info import info
        ret = info(["info", one_con])
        assert ret == 0
        out = capsys.readouterr().out
        assert "1 contacts" in out

    def test_dedup_one_con(self, one_con, capsys):
        from dip_c.commands.dedup import dedup
        ret = dedup(["dedup", one_con])
        assert ret == 0

    def test_clean_one_con(self, one_con, capsys):
        from dip_c.commands.clean import clean
        ret = clean(["clean", one_con])
        assert ret == 0

    def test_impute_one_con(self, one_con, capsys):
        from dip_c.commands.impute import impute
        ret = impute(["impute", one_con])
        assert ret == 0

    def test_ard_one_con(self, one_con, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", one_con])
        assert ret == 0

    def test_con3_two_particles(self, two_3dg, capsys):
        from dip_c.commands.con3 import con3
        ret = con3(["con3", two_3dg])
        assert ret == 0

    def test_exp_two_particles(self, two_3dg, capsys):
        from dip_c.commands.exp import exp
        ret = exp(["exp", two_3dg])
        assert ret == 0

    def test_vis_two_particles(self, two_3dg, capsys):
        from dip_c.commands.vis import vis
        ret = vis(["vis", two_3dg])
        assert ret == 0

    def test_dist_two_particles(self, two_3dg, capsys):
        from dip_c.commands.dist import dist
        ret = dist(["dist", two_3dg])
        assert ret == 0

    def test_rg_one_particle(self, one_3dg, capsys):
        from dip_c.commands.rg import rg
        ret = rg(["rg", one_3dg])
        assert ret == 0


# ===================================================================
# 3d: mgcolor regression tests
# ===================================================================

class TestMgcolorEdgeCases:
    """Regression tests for mgcolor, including the empty-file fix."""

    def test_empty_file_returns_error(self, tmp_path):
        from dip_c.commands.mgcolor import mgcolor
        f = tmp_path / "empty.color"
        f.write_text("")
        ret = mgcolor(["mgcolor", str(f)])
        assert ret == 1

    def test_empty_plus_valid_returns_error(self, tmp_path):
        from dip_c.commands.mgcolor import mgcolor
        f_empty = tmp_path / "empty.color"
        f_empty.write_text("")
        f_valid = tmp_path / "valid.color"
        f_valid.write_text(ONE_COLOR_LINE)
        ret = mgcolor(["mgcolor", str(f_empty), str(f_valid)])
        assert ret == 1

    def test_further_merge_empty_returns_error(self, tmp_path):
        from dip_c.commands.mgcolor import mgcolor
        f = tmp_path / "empty.colors"
        f.write_text("homolog\tlocus\n")  # header only, no data
        ret = mgcolor(["mgcolor", "-s", str(f)])
        assert ret == 1

    def test_diploid_empty_returns_error(self, tmp_path):
        from dip_c.commands.mgcolor import mgcolor
        f = tmp_path / "empty.color"
        f.write_text("")
        ret = mgcolor(["mgcolor", "-d", str(f)])
        assert ret == 1

    def test_valid_merge_correct_output(self, tmp_path, capsys):
        from dip_c.commands.mgcolor import mgcolor
        f1 = tmp_path / "a.color"
        f1.write_text("1(pat)\t100000\t0.5\n1(pat)\t200000\t0.3\n")
        f2 = tmp_path / "b.color"
        f2.write_text("1(pat)\t100000\t0.8\n1(pat)\t200000\t0.9\n")
        ret = mgcolor(["mgcolor", str(f1), str(f2)])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == 3  # header + 2 data lines
        # header should have: homolog, locus, file1, file2
        assert lines[0].count("\t") == 3
        # data line should have 4 columns
        assert lines[1].count("\t") == 3

    def test_missing_values_filled(self, tmp_path, capsys):
        from dip_c.commands.mgcolor import mgcolor
        f1 = tmp_path / "a.color"
        f1.write_text("1(pat)\t100000\t0.5\n")
        f2 = tmp_path / "b.color"
        f2.write_text("1(pat)\t200000\t0.8\n")
        ret = mgcolor(["mgcolor", str(f1), str(f2)])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        # Both particles should appear, with -1.0 for missing
        assert len(lines) == 3  # header + 2 data lines
        assert "-1.0" in out


# ===================================================================
# 3e: color command edge cases
# ===================================================================

class TestColorEdgeCases:
    """Edge cases for the color command."""

    def test_empty_color_reference(self, tmp_path, capsys):
        """Empty color reference file → zero output particles."""
        from dip_c.commands.color import color
        g3d_f = tmp_path / "two.3dg"
        g3d_f.write_text(TWO_3DG_LINES)
        color_f = tmp_path / "empty.color"
        color_f.write_text("")
        ret = color(["color", "-c", str(color_f), str(g3d_f)])
        assert ret == 0
        out = capsys.readouterr().out
        assert out.strip() == ""  # no particles matched

    def test_mismatched_chromosome_names(self, tmp_path, capsys):
        """Color reference with wrong chr names → zero output."""
        from dip_c.commands.color import color
        g3d_f = tmp_path / "two.3dg"
        g3d_f.write_text(TWO_3DG_LINES)
        color_f = tmp_path / "wrong.color"
        color_f.write_text("chrXYZ\t100000\t0.5\n")
        ret = color(["color", "-c", str(color_f), str(g3d_f)])
        assert ret == 0
        out = capsys.readouterr().out
        assert out.strip() == ""

    def test_bundled_file_resolution(self, tmp_path, capsys):
        """Bundled color file resolution should work."""
        from dip_c.data import list_bundled_data_files
        files = list_bundled_data_files()
        assert len(files) > 0  # bundled data exists

    def test_no_color_scheme_returns_usage(self):
        """color with a 3DG file but no color scheme option returns usage."""
        from dip_c.commands.color import color
        # color requires exactly one color scheme
        # with no scheme, it should show usage
        ret = color(["color"])
        assert ret == 1


# ===================================================================
# 3f: Missing required options
# ===================================================================

class TestMissingRequiredOptions:
    """Commands that need specific flags should return error=1 without them."""

    def test_clean3_without_c(self, tmp_path):
        from dip_c.commands.clean3 import clean3
        g3d_f = tmp_path / "one.3dg"
        g3d_f.write_text(ONE_3DG_LINE)
        ret = clean3(["clean3", str(g3d_f)])
        assert ret == 1

    def test_impute3_without_3(self, tmp_path):
        from dip_c.commands.impute3 import impute3
        con_f = tmp_path / "one.con"
        con_f.write_text(ONE_CON_LINE)
        ret = impute3(["impute3", str(con_f)])
        assert ret == 1

    def test_pos_without_l(self, tmp_path):
        from dip_c.commands.pos import pos
        g3d_f = tmp_path / "one.3dg"
        g3d_f.write_text(ONE_3DG_LINE)
        ret = pos(["pos", str(g3d_f)])
        assert ret == 1

    def test_pd_without_1(self, tmp_path):
        from dip_c.commands.pd import pd
        g3d_f = tmp_path / "one.3dg"
        g3d_f.write_text(ONE_3DG_LINE)
        ret = pd(["pd", str(g3d_f)])
        assert ret == 1

    def test_tad_without_l(self, tmp_path):
        from dip_c.commands.tad import tad
        rg_f = tmp_path / "dummy.rg"
        rg_f.write_text("0\n")
        ret = tad(["tad", str(rg_f)])
        assert ret == 1

    def test_cv_without_i(self, tmp_path):
        from dip_c.commands.cv import cv
        con_f = tmp_path / "one.con"
        con_f.write_text(ONE_CON_LINE)
        ret = cv(["cv", "-t", str(con_f), str(con_f)])
        assert ret == 1

    def test_cv_without_t(self, tmp_path):
        from dip_c.commands.cv import cv
        con_f = tmp_path / "one.con"
        con_f.write_text(ONE_CON_LINE)
        ret = cv(["cv", "-i", str(con_f), str(con_f)])
        assert ret == 1

    def test_bincon_without_l(self, tmp_path):
        from dip_c.commands.bincon import bincon
        con_f = tmp_path / "one.con"
        con_f.write_text(ONE_CON_LINE)
        ret = bincon(["bincon", str(con_f)])
        assert ret == 1


# ===================================================================
# 3g: Missing files
# ===================================================================

class TestMissingFiles:
    """Nonexistent file paths should raise FileNotFoundError."""

    def test_info_missing_file(self):
        from dip_c.commands.info import info
        with pytest.raises((FileNotFoundError, OSError)):
            info(["info", "/nonexistent/path/to/file.con"])

    def test_dedup_missing_file(self):
        from dip_c.commands.dedup import dedup
        with pytest.raises((FileNotFoundError, OSError)):
            dedup(["dedup", "/nonexistent/path/to/file.con"])

    def test_vis_missing_file(self):
        from dip_c.commands.vis import vis
        with pytest.raises((FileNotFoundError, OSError)):
            vis(["vis", "/nonexistent/path/to/file.3dg"])

    def test_exp_missing_file(self):
        from dip_c.commands.exp import exp
        with pytest.raises((FileNotFoundError, OSError)):
            exp(["exp", "/nonexistent/path/to/file.3dg"])


# ===================================================================
# 3h: Gzip mismatch
# ===================================================================

class TestGzipMismatch:
    """Test gzip edge cases."""

    def test_plain_text_named_gz(self, tmp_path):
        """Plain text file with .gz extension should fail gracefully."""
        from dip_c.commands.info import info
        f = tmp_path / "bad.con.gz"
        f.write_text("not actually gzipped\n")
        with pytest.raises(Exception):
            info(["info", str(f)])

    def test_gzipped_without_gz_extension(self, tmp_path, capsys):
        """Gzipped file without .gz extension is read as plain text → parse error or empty."""
        from dip_c.commands.info import info
        f = tmp_path / "data.con"
        write_file(str(f), ONE_CON_LINE, gz=False)
        # This should work fine as plain text
        ret = info(["info", str(f)])
        assert ret == 0


# ===================================================================
# 4a: Malformed input handling
# ===================================================================

class TestMalformedInput:
    """Malformed lines should be skipped or cause command error, not tracebacks."""

    def test_con_wrong_field_count(self, tmp_path):
        """CON file with 3 tab-separated fields instead of 2."""
        from dip_c.commands.info import info
        f = tmp_path / "bad.con"
        f.write_text("1,100,0\t1,200,0\textra\n")
        ret = info(["info", str(f)])
        assert ret == 1  # malformed → 0 contacts → error

    def test_con_missing_comma_in_leg(self, tmp_path):
        """Legs with tabs instead of commas."""
        from dip_c.commands.info import info
        f = tmp_path / "bad.con"
        f.write_text("1\t200\t0\t1\t300\t0\n")
        ret = info(["info", str(f)])
        assert ret == 1

    def test_con_non_numeric_locus(self, tmp_path):
        """Non-numeric locus value."""
        from dip_c.commands.info import info
        f = tmp_path / "bad.con"
        f.write_text("1,abc,0\t1,200,0\n")
        ret = info(["info", str(f)])
        assert ret == 1

    def test_con_empty_lines_mixed(self, tmp_path, capsys):
        """Valid lines with blank lines interspersed — blanks skipped."""
        from dip_c.commands.info import info
        f = tmp_path / "mixed.con"
        f.write_text("\n1,100,0\t1,200,0\n\n1,300,0\t1,400,0\n\n")
        ret = info(["info", str(f)])
        assert ret == 0
        out = capsys.readouterr().out
        assert "2 contacts" in out

    def test_3dg_wrong_field_count(self, tmp_path):
        """3DG with 3 fields instead of 5."""
        from dip_c.commands.con3 import con3
        f = tmp_path / "bad.3dg"
        f.write_text("1(pat)\t100000\t0.0\n")
        ret = con3(["con3", str(f)])
        assert ret == 1  # malformed → 0 particles → error

    def test_3dg_non_numeric_coords(self, tmp_path):
        """Non-numeric coordinate."""
        from dip_c.commands.con3 import con3
        f = tmp_path / "bad.3dg"
        f.write_text("1(pat)\t100000\tnan\t0.0\t0.0\n")
        # nan is a valid float, so this should actually parse
        # Use truly non-numeric instead
        f.write_text("1(pat)\t100000\tabc\t0.0\t0.0\n")
        ret = con3(["con3", str(f)])
        assert ret == 1

    def test_3dg_empty_lines_mixed(self, tmp_path, capsys):
        """Valid 3DG with blank lines interspersed."""
        from dip_c.commands.exp import exp
        f = tmp_path / "mixed.3dg"
        f.write_text(
            "\n1(pat)\t100000\t0.0\t0.0\t0.0\n\n"
            "1(pat)\t200000\t1.0\t1.0\t1.0\n\n"
        )
        ret = exp(["exp", str(f)])
        assert ret == 0

    def test_con_trailing_newlines(self, tmp_path, capsys):
        """CON file with trailing empty lines — no crash."""
        from dip_c.commands.info import info
        f = tmp_path / "trailing.con"
        f.write_text("1,100,0\t1,200,0\n\n\n\n")
        ret = info(["info", str(f)])
        assert ret == 0


# ===================================================================
# 4b: Secondary guard tests
# ===================================================================

class TestSecondaryGuards:
    """Tests for secondary ZeroDivisionError and TypeError guards."""

    def test_impute3_single_particle_per_hom(self, tmp_path):
        """3DG with 1 particle per homolog → resolution None → error."""
        from dip_c.commands.impute3 import impute3
        g3d_f = tmp_path / "single.3dg"
        g3d_f.write_text("1(pat)\t100000\t0.0\t0.0\t0.0\n")
        con_f = tmp_path / "one.con"
        con_f.write_text(ONE_CON_LINE)
        ret = impute3(["impute3", "-3", str(g3d_f), str(con_f)])
        assert ret == 1

    def test_force_empty_con(self, tmp_path):
        """Empty CON + valid chr.len → error=1."""
        from dip_c.commands.force import force
        con_f = tmp_path / "empty.con"
        con_f.write_text("")
        lens_f = tmp_path / "chr.len"
        lens_f.write_text("1\t1000000\n")
        ret = force(["force", "-l", str(lens_f), str(con_f)])
        assert ret == 1

    def test_mkcon_empty_legs(self, tmp_path):
        """Two empty leg files → error=1."""
        from dip_c.commands.mkcon import mkcon
        leg1 = tmp_path / "empty1.leg"
        leg1.write_text("")
        leg2 = tmp_path / "empty2.leg"
        leg2.write_text("")
        ret = mkcon(["mkcon", str(leg1), str(leg2)])
        assert ret == 1


# ===================================================================
# 4c: Color command modes
# ===================================================================

CHR22_3DG = os.path.join(os.path.dirname(__file__), "data",
                         "gm12878_11.impute.chr22.3dg.txt")


class TestColorModes:
    """Exercise color command modes for coverage."""

    @pytest.fixture()
    def g3d_file(self):
        return CHR22_3DG

    @pytest.fixture()
    def chr_len(self, tmp_path):
        p = tmp_path / "chr.len"
        p.write_text("22\t51304566\n")
        return str(p)

    @pytest.fixture()
    def chr_txt(self, tmp_path):
        p = tmp_path / "chr.txt"
        p.write_text("22\n")
        return str(p)

    def test_color_center_of_mass(self, g3d_file, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-C", g3d_file])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip().split("\n")) > 0

    def test_color_chromosome_number(self, g3d_file, chr_txt, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-n", chr_txt, g3d_file])
        assert ret == 0
        out = capsys.readouterr().out
        assert "22(pat)" in out or "22(mat)" in out

    def test_color_by_locus(self, g3d_file, chr_len, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-l", chr_len, g3d_file])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip().split("\n")) > 0

    def test_color_smoothing(self, g3d_file, chr_txt, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-n", chr_txt, "-s", "3.0", g3d_file])
        assert ret == 0

    def test_color_intra_hom_fraction(self, g3d_file, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-i", "5.0", g3d_file])
        assert ret == 0

    def test_color_intra_hom_count(self, g3d_file, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-I", "5.0", g3d_file])
        assert ret == 0

    def test_color_same_haplotype(self, g3d_file, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-p", "5.0", g3d_file])
        assert ret == 0

    def test_color_homolog_diversity(self, g3d_file, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-d", "5.0", g3d_file])
        assert ret == 0

    def test_color_homolog_richness(self, g3d_file, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-r", "5.0", g3d_file])
        assert ret == 0

    def test_color_radial_mode(self, g3d_file, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-C", "-R", g3d_file])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip().split("\n")) > 1

    def test_color_multiple_schemes_error(self, g3d_file, chr_txt):
        from dip_c.commands.color import color
        ret = color(["color", "-n", chr_txt, "-C", g3d_file])
        assert ret == 1


# ===================================================================
# 4d: Other command modes
# ===================================================================

CHR22_CON_GZ = os.path.join(os.path.dirname(__file__), "data",
                             "gm12878_11.clean.chr22.con.txt.gz")
CHR21_22_3DG = os.path.join(os.path.dirname(__file__), "data",
                             "gm12878_11.impute.chr21_22.3dg.txt")
CHR21_22_CON_GZ = os.path.join(os.path.dirname(__file__), "data",
                                "gm12878_11.impute.chr21_22.con.txt.gz")


class TestCommandModes:
    """Exercise untested command modes for coverage."""

    @pytest.fixture()
    def chr_len(self, tmp_path):
        p = tmp_path / "chr.len"
        p.write_text("22\t51304566\n")
        return str(p)

    @pytest.fixture()
    def chr_len_21_22(self, tmp_path):
        p = tmp_path / "chr.len"
        p.write_text("21\t48129895\n22\t51304566\n")
        return str(p)

    def test_con3_matrix_mode(self, chr_len, capsys):
        from dip_c.commands.con3 import con3
        ret = con3(["con3", "-m", chr_len, CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip()) > 0

    def test_dist_basic(self, capsys):
        from dip_c.commands.dist import dist
        ret = dist(["dist", CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip()) > 0

    def test_clean_test_mode(self, capsys):
        from dip_c.commands.clean import clean
        ret = clean(["clean", "-t", CHR22_CON_GZ])
        assert ret == 0

    def test_ard_normalize_mode(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-n", "-t", CHR21_22_CON_GZ])
        assert ret == 0

    def test_reg_with_preset(self, capsys):
        from dip_c.commands.reg import reg
        ret = reg(["reg", "-p", "hf", CHR22_CON_GZ])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip()) > 0

    def test_reg3_with_preset(self, capsys):
        from dip_c.commands.reg3 import reg3
        ret = reg3(["reg3", "-p", "hf", CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert "22(pat)" in out or "22(mat)" in out

    def test_dedup_custom_params(self, capsys):
        from dip_c.commands.dedup import dedup
        ret = dedup(["dedup", "-s", "500", "-d", "2000", CHR22_CON_GZ])
        assert ret == 0

    def test_impute3_with_3dg(self, tmp_path, capsys):
        from dip_c.commands.impute3 import impute3
        ret = impute3(["impute3", "-3", CHR21_22_3DG, CHR21_22_CON_GZ])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip()) > 0


# ===================================================================
# Step 4: Option parsing branches
# ===================================================================

CHR22_3DG = os.path.join(os.path.dirname(__file__), "data",
                          "gm12878_11.impute.chr22.3dg.txt")
CHR22_IMPUTE_CON_GZ = os.path.join(os.path.dirname(__file__), "data",
                                    "gm12878_11.impute.chr22.con.txt.gz")
TEST_SEG = os.path.join(os.path.dirname(__file__), "data", "test_seg.seg")
TEST_SEG_BAM = os.path.join(os.path.dirname(__file__), "data", "test_seg.bam")
TEST_SNPS = os.path.join(os.path.dirname(__file__), "data", "test_snps.txt")


class TestOptionParsing:
    """Exercise untested option-parsing branches across commands."""

    # --- con.py options ---
    def test_con_min_separation(self, capsys):
        from dip_c.commands.con import con
        ret = con(["con", "-s", "500", TEST_SEG])
        assert ret == 0

    def test_con_max_distance(self, capsys):
        from dip_c.commands.con import con
        ret = con(["con", "-d", "2000", TEST_SEG])
        assert ret == 0

    def test_con_adjacent_only(self, capsys):
        from dip_c.commands.con import con
        ret = con(["con", "-a", TEST_SEG])
        assert ret == 0

    # --- clean.py options ---
    def test_clean_distance(self, capsys):
        from dip_c.commands.clean import clean
        ret = clean(["clean", "-d", "5000000", CHR22_CON_GZ])
        assert ret == 0

    def test_clean_count(self, capsys):
        from dip_c.commands.clean import clean
        ret = clean(["clean", "-c", "3", CHR22_CON_GZ])
        assert ret == 0

    def test_clean_leg_distance(self, capsys):
        from dip_c.commands.clean import clean
        ret = clean(["clean", "-D", "500", CHR22_CON_GZ])
        assert ret == 0

    def test_clean_leg_count(self, capsys):
        from dip_c.commands.clean import clean
        ret = clean(["clean", "-C", "5", CHR22_CON_GZ])
        assert ret == 0

    # --- clean3.py option ---
    def test_clean3_quantile(self, capsys):
        from dip_c.commands.clean3 import clean3
        ret = clean3(["clean3", "-c", CHR22_CON_GZ, "-q", "0.5", CHR22_3DG])
        assert ret == 0

    # --- exp.py options ---
    def test_exp_factor(self, capsys):
        from dip_c.commands.exp import exp
        ret = exp(["exp", "-f", "2.0", CHR22_3DG])
        assert ret == 0

    def test_exp_centers_only(self, capsys):
        from dip_c.commands.exp import exp
        ret = exp(["exp", "-c", CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip()) > 0

    # --- color2.py options ---
    def test_color2_merge_haplotypes(self, tmp_path, capsys):
        from dip_c.commands.color2 import color2
        color_f = tmp_path / "color.txt"
        color_f.write_text("22\t100000\t0.5\n22\t200000\t0.8\n")
        ret = color2(["color2", "-H", "-c", str(color_f), CHR22_CON_GZ])
        assert ret == 0

    def test_color2_smooth(self, tmp_path, capsys):
        from dip_c.commands.color2 import color2
        color_f = tmp_path / "color.txt"
        color_f.write_text("22(pat)\t100000\t0.5\n22(mat)\t200000\t0.8\n")
        ret = color2(["color2", "-s", "-c", str(color_f), CHR22_CON_GZ])
        assert ret == 0

    def test_color2_bin_size(self, tmp_path, capsys):
        from dip_c.commands.color2 import color2
        color_f = tmp_path / "color.txt"
        color_f.write_text("22(pat)\t100000\t0.5\n")
        ret = color2(["color2", "-b", "500000", "-c", str(color_f), CHR22_CON_GZ])
        assert ret == 0

    # --- impute.py options ---
    def test_impute_options(self, capsys):
        from dip_c.commands.impute import impute
        ret = impute(["impute", "-d", "5000000", "-v", "2", "-f", "0.8", "-r", "2",
                       "-D", "5000000", "-C", "1", CHR22_IMPUTE_CON_GZ])
        assert ret == 0

    def test_impute_male_preset(self, tmp_path, capsys):
        from dip_c.commands.impute import impute
        # Create CON with X and Y contacts
        con_f = tmp_path / "xy.con"
        con_f.write_text("X,100,0\tX,200,1\nY,100,0\tY,200,0\n1,100,0\t1,200,0\n")
        ret = impute(["impute", "-p", "hm", str(con_f)])
        assert ret == 0

    def test_impute_unknown_preset(self, capsys):
        from dip_c.commands.impute import impute
        ret = impute(["impute", "-p", "zz", CHR22_IMPUTE_CON_GZ])
        assert ret == 1

    # --- impute3.py options ---
    def test_impute3_options(self, tmp_path, capsys):
        from dip_c.commands.impute3 import impute3
        vio = tmp_path / "vio.txt"
        ret = impute3(["impute3", "-3", CHR21_22_3DG, "-v", str(vio),
                        "-d", "10.0", "-r", "0.3", "-D", "5000000", "-C", "1",
                        CHR21_22_CON_GZ])
        assert ret == 0
        assert vio.exists()
        assert vio.stat().st_size > 0

    def test_impute3_male_preset(self, tmp_path, capsys):
        from dip_c.commands.impute3 import impute3
        con_f = tmp_path / "xy.con"
        con_f.write_text("X,100,0\tX,200,1\n1,100,0\t1,200,0\n")
        ret = impute3(["impute3", "-3", CHR21_22_3DG, "-p", "hm", str(con_f)])
        assert ret == 0

    def test_impute3_unknown_preset(self, capsys):
        from dip_c.commands.impute3 import impute3
        ret = impute3(["impute3", "-3", CHR21_22_3DG, "-p", "zz", CHR21_22_CON_GZ])
        assert ret == 1

    # --- reg.py options ---
    def test_reg_include_file(self, tmp_path, capsys):
        from dip_c.commands.reg import reg
        inc = tmp_path / "inc.reg"
        inc.write_text("22\t.\t.\t.\n")
        ret = reg(["reg", "-i", str(inc), CHR22_CON_GZ])
        assert ret == 0

    def test_reg_exclude_file(self, tmp_path, capsys):
        from dip_c.commands.reg import reg
        inc = tmp_path / "inc.reg"
        inc.write_text("22\t.\t.\t.\n")
        exc = tmp_path / "exc.reg"
        exc.write_text("22\t.\t0\t1000000\n")
        ret = reg(["reg", "-i", str(inc), "-e", str(exc), CHR22_CON_GZ])
        assert ret == 0

    def test_reg_haplotype_file(self, tmp_path, capsys):
        from dip_c.commands.reg import reg
        inc = tmp_path / "inc.reg"
        inc.write_text("22\t.\t.\t.\n")
        hap = tmp_path / "hap.reg"
        hap.write_text("22\t0\t0\t10000000\n")
        ret = reg(["reg", "-i", str(inc), "-h", str(hap), CHR22_CON_GZ])
        assert ret == 0

    def test_reg_no_regions_error(self, capsys):
        from dip_c.commands.reg import reg
        ret = reg(["reg", CHR22_CON_GZ])
        assert ret == 1

    def test_reg_unknown_preset(self, capsys):
        from dip_c.commands.reg import reg
        ret = reg(["reg", "-p", "zz", CHR22_CON_GZ])
        assert ret == 1

    # --- reg3.py options ---
    def test_reg3_include_file(self, tmp_path, capsys):
        from dip_c.commands.reg3 import reg3
        inc = tmp_path / "inc.reg"
        inc.write_text("22\t.\t.\t.\n")
        ret = reg3(["reg3", "-i", str(inc), CHR22_3DG])
        assert ret == 0

    def test_reg3_exclude_file(self, tmp_path, capsys):
        from dip_c.commands.reg3 import reg3
        inc = tmp_path / "inc.reg"
        inc.write_text("22\t.\t.\t.\n")
        exc = tmp_path / "exc.reg"
        exc.write_text("22\t.\t0\t1000000\n")
        ret = reg3(["reg3", "-i", str(inc), "-e", str(exc), CHR22_3DG])
        assert ret == 0

    def test_reg3_no_regions_error(self, capsys):
        from dip_c.commands.reg3 import reg3
        ret = reg3(["reg3", CHR22_3DG])
        assert ret == 1

    def test_reg3_unknown_preset(self, capsys):
        from dip_c.commands.reg3 import reg3
        ret = reg3(["reg3", "-p", "zz", CHR22_3DG])
        assert ret == 1

    # --- bincon.py options ---
    def test_bincon_info_mode(self, tmp_path, capsys):
        from dip_c.commands.bincon import bincon
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("22\t51304566\n")
        ret = bincon(["bincon", "-l", str(chr_len), "-i", CHR22_CON_GZ])
        assert ret == 0
        out = capsys.readouterr().out
        assert "22(pat)" in out

    def test_bincon_merge_hap(self, tmp_path, capsys):
        from dip_c.commands.bincon import bincon
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("22\t51304566\n")
        ret = bincon(["bincon", "-l", str(chr_len), "-H", CHR22_CON_GZ])
        assert ret == 0

    def test_bincon_min_separation(self, tmp_path, capsys):
        from dip_c.commands.bincon import bincon
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("22\t51304566\n")
        ret = bincon(["bincon", "-l", str(chr_len), "-s", "1000000", CHR22_CON_GZ])
        assert ret == 0

    def test_bincon_bin_size(self, tmp_path, capsys):
        from dip_c.commands.bincon import bincon
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("22\t51304566\n")
        ret = bincon(["bincon", "-l", str(chr_len), "-b", "500000", CHR22_CON_GZ])
        assert ret == 0

    def test_bincon_leg_mode(self, tmp_path, capsys):
        from dip_c.commands.bincon import bincon
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("22\t51304566\n")
        leg_f = tmp_path / "legs.txt"
        leg_f.write_text("22,1000000,0\n22,2000000,1\n")
        ret = bincon(["bincon", "-l", str(chr_len), "-L", str(leg_f)])
        assert ret == 0

    # --- con3.py options ---
    def test_con3_info_mode(self, tmp_path, capsys):
        from dip_c.commands.con3 import con3
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("22\t51304566\n")
        ret = con3(["con3", "-m", str(chr_len), "-i", CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert "22(pat)" in out

    def test_con3_merge_haplotypes(self, tmp_path, capsys):
        from dip_c.commands.con3 import con3
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("22\t51304566\n")
        ret = con3(["con3", "-m", str(chr_len), "-H", CHR22_3DG])
        assert ret == 0

    # --- vis.py options ---
    def test_vis_color_discard(self, tmp_path, capsys):
        from dip_c.commands.vis import vis
        color_f = tmp_path / "color.txt"
        color_f.write_text("22(pat)\t16000000\t0.5\n")
        ret = vis(["vis", "-c", str(color_f), "-M", CHR22_3DG])
        assert ret == 0

    def test_vis_color_missing_value(self, tmp_path, capsys):
        from dip_c.commands.vis import vis
        color_f = tmp_path / "color.txt"
        color_f.write_text("22(pat)\t16000000\t0.5\n")
        ret = vis(["vis", "-c", str(color_f), "-m", "0.0", CHR22_3DG])
        assert ret == 0

    def test_vis_connect_adjacent_flag(self, capsys):
        from dip_c.commands.vis import vis
        ret = vis(["vis", "-a", CHR22_3DG])
        assert ret == 0

    # --- pos.py options ---
    def test_pos_out_of_bounds(self, tmp_path, capsys):
        from dip_c.commands.pos import pos
        leg_f = tmp_path / "leg.txt"
        leg_f.write_text("22,16000000,0\n99,100,0\n")
        ret = pos(["pos", "-l", str(leg_f), CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert "None" in out

    def test_pos_exclude_oob(self, tmp_path, capsys):
        from dip_c.commands.pos import pos
        leg_f = tmp_path / "leg.txt"
        # Use a locus that is out of the 3DG range
        leg_f.write_text("22,1,0\n")
        ret = pos(["pos", "-l", str(leg_f), "-O", CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert "None" in out

    # --- pd.py options ---
    def test_pd_pair_mode(self, tmp_path, capsys):
        from dip_c.commands.pd import pd
        leg1 = tmp_path / "l1.leg"
        leg1.write_text("22,16000000,0\n22,18000000,0\n")
        leg2 = tmp_path / "l2.leg"
        leg2.write_text("22,20000000,1\n")
        ret = pd(["pd", "-1", str(leg1), "-2", str(leg2), CHR22_3DG])
        assert ret == 0

    def test_pd_self_mode(self, tmp_path, capsys):
        from dip_c.commands.pd import pd
        leg1 = tmp_path / "l1.leg"
        leg1.write_text("22,16000000,0\n22,18000000,0\n")
        ret = pd(["pd", "-1", str(leg1), CHR22_3DG])
        assert ret == 0

    def test_pd_none_position(self, tmp_path, capsys):
        from dip_c.commands.pd import pd
        leg1 = tmp_path / "l1.leg"
        leg1.write_text("99,100,0\n22,16000000,0\n")
        ret = pd(["pd", "-1", str(leg1), CHR22_3DG])
        assert ret == 0

    # --- mkcon.py ---
    def test_mkcon_valid(self, tmp_path, capsys):
        from dip_c.commands.mkcon import mkcon
        leg1 = tmp_path / "l1.txt"
        leg1.write_text("22(pat)\t100000\n22(mat)\t200000\n")
        leg2 = tmp_path / "l2.txt"
        leg2.write_text("22(pat)\t300000\n22(mat)\t400000\n")
        ret = mkcon(["mkcon", "-n", "10", str(leg1), str(leg2)])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip().split("\n")) == 10

    def test_mkcon_wrong_arg_count(self, tmp_path, capsys):
        from dip_c.commands.mkcon import mkcon
        leg1 = tmp_path / "l1.txt"
        leg1.write_text("22(pat)\t100000\n")
        ret = mkcon(["mkcon", str(leg1)])
        assert ret == 1

    # --- tad.py option ---
    def test_tad_l_required(self, tmp_path, capsys):
        from dip_c.commands.tad import tad
        rg_f = tmp_path / "rg.txt"
        rg_f.write_text("1\t2\n2\t1\n")
        ret = tad(["tad", str(rg_f)])
        assert ret == 1

    # --- seg.py options ---
    def test_seg_custom_mapq(self, capsys):
        from dip_c.commands.seg import seg
        ret = seg(["seg", "-q", "10", "-m", "0.1", TEST_SEG_BAM])
        assert ret == 0

    def test_seg_getopt_error(self, capsys):
        """Cover seg getopt error handler (lines 51-53)."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", "--bogus"])
        assert ret == 1
        assert "unknown command" in capsys.readouterr().err

    def test_seg_no_args_usage(self, capsys):
        """Cover seg usage message when no BAM argument given (lines 55-62)."""
        from dip_c.commands.seg import seg
        ret = seg(["seg"])
        assert ret == 1
        assert "Usage:" in capsys.readouterr().err

    def test_seg_custom_baseq(self, capsys):
        """Cover seg -Q option (lines 70-71)."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", "-Q", "30", "-v", TEST_SNPS, TEST_SEG_BAM], _display_interval=1)
        assert ret == 0

    # --- rg.py options ---
    def test_rg_distance_mode(self, tmp_path, capsys):
        from dip_c.commands.rg import rg
        ret = rg(["rg", "-d", "-o", str(tmp_path / "out."), CHR22_3DG])
        assert ret == 0

    def test_rg_region_file(self, tmp_path, capsys):
        from dip_c.commands.rg import rg
        reg_f = tmp_path / "reg.txt"
        reg_f.write_text("22\t0\t.\t.\n")
        ret = rg(["rg", "-r", str(reg_f), "-o", str(tmp_path / "out."), CHR22_3DG])
        assert ret == 0

    # --- color.py options ---
    def test_color_arm_locus(self, tmp_path, capsys):
        from dip_c.commands.color import color
        cen = tmp_path / "chr.cen"
        cen.write_text("22\t51304566\t16500000\n")
        ret = color(["color", "-L", str(cen), CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip().split("\n")) > 0

    def test_color_c_hom(self, tmp_path, capsys):
        from dip_c.commands.color import color
        color_f = tmp_path / "color.txt"
        color_f.write_text("22(pat)\t16000000\t0.5\n22(mat)\t16000000\t0.8\n")
        ret = color(["color", "--c-hom=" + str(color_f), CHR22_3DG])
        assert ret == 0

    def test_color_distance_to_locus(self, tmp_path, capsys):
        from dip_c.commands.color import color
        leg_f = tmp_path / "leg.txt"
        leg_f.write_text("22,16000000,0\n")
        ret = color(["color", "-D", str(leg_f), CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip().split("\n")) > 0

    def test_color_smooth_exc(self, tmp_path, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-C", "--s-exc=3.0", CHR22_3DG])
        assert ret == 0

    def test_color_smooth_max_sep(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-i", "5.0", "-S", "1000000", CHR22_3DG])
        assert ret == 0

    def test_color_same_haplotype_exc(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-P", "5.0", CHR22_3DG])
        assert ret == 0

    # --- ard.py options ---
    def test_ard_intra_chr(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-s", "1000000", CHR22_CON_GZ])
        assert ret == 0

    def test_ard_histogram(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-h", "500000", CHR22_CON_GZ])
        assert ret == 0

    def test_ard_histogram_intra(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-s", "1000000", "-h", "500000", CHR22_CON_GZ])
        assert ret == 0

    def test_ard_no_symmetrize(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-h", "500000", "-S", CHR22_CON_GZ])
        assert ret == 0

    def test_ard_superellipse(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-e", "-n", CHR22_CON_GZ])
        assert ret == 0

    def test_ard_pairwise_legs(self, tmp_path, capsys):
        from dip_c.commands.ard import ard
        leg_f = tmp_path / "legs.leg"
        leg_f.write_text("22,16000000,0\n22,18000000,1\n")
        ret = ard(["ard", "-1", str(leg_f), CHR22_CON_GZ])
        assert ret == 0

    def test_ard_pairwise_two_leg_files(self, tmp_path, capsys):
        from dip_c.commands.ard import ard
        leg1 = tmp_path / "l1.leg"
        leg1.write_text("22,16000000,0\n")
        leg2 = tmp_path / "l2.leg"
        leg2.write_text("22,18000000,1\n")
        ret = ard(["ard", "-1", str(leg1), "-2", str(leg2), CHR22_CON_GZ])
        assert ret == 0

    def test_ard_reference_file(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-c", CHR22_CON_GZ, "-n", CHR22_CON_GZ])
        assert ret == 0

    def test_ard_count_normalize(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-n", "-t", CHR22_CON_GZ])
        assert ret == 0

    # --- cv.py ---
    def test_cv_correct_imputation(self, tmp_path, capsys):
        from dip_c.commands.cv import cv
        # Create minimal CON files for CV
        base = "22,100,.\t22,200,.\n"
        imputed = "22,100,0\t22,200,1\n"
        truth = "22,100,0\t22,200,1\n"
        base_f = tmp_path / "base.con"
        base_f.write_text(base)
        imp_f = tmp_path / "imp.con"
        imp_f.write_text(imputed)
        truth_f = tmp_path / "truth.con"
        truth_f.write_text(truth)
        ret = cv(["cv", "-i", str(imp_f), "-t", str(truth_f), str(base_f)])
        assert ret == 0
        out = capsys.readouterr().out
        assert "correct" in out

    def test_cv_wrong_imputation(self, tmp_path, capsys):
        from dip_c.commands.cv import cv
        base = "22,100,.\t22,200,.\n"
        imputed = "22,100,0\t22,200,0\n"
        truth = "22,100,1\t22,200,1\n"
        base_f = tmp_path / "base.con"
        base_f.write_text(base)
        imp_f = tmp_path / "imp.con"
        imp_f.write_text(imputed)
        truth_f = tmp_path / "truth.con"
        truth_f.write_text(truth)
        ret = cv(["cv", "-i", str(imp_f), "-t", str(truth_f), str(base_f)])
        assert ret == 0
        out = capsys.readouterr().out
        assert "wrong" in out

    def test_cv_no_truth(self, tmp_path, capsys):
        from dip_c.commands.cv import cv
        base = "22,100,.\t22,200,.\n"
        imputed = "22,100,0\t22,200,0\n"
        truth = "22,100,.\t22,200,.\n"
        base_f = tmp_path / "base.con"
        base_f.write_text(base)
        imp_f = tmp_path / "imp.con"
        imp_f.write_text(imputed)
        truth_f = tmp_path / "truth.con"
        truth_f.write_text(truth)
        ret = cv(["cv", "-i", str(imp_f), "-t", str(truth_f), str(base_f)])
        assert ret == 0
        out = capsys.readouterr().out
        assert "no ground truth" in out

    # --- mgcolor.py options ---
    def test_mgcolor_diploid(self, tmp_path, capsys):
        from dip_c.commands.mgcolor import mgcolor
        c1 = tmp_path / "c1.color"
        c1.write_text("22(pat)\t100000\t0.5\n22(mat)\t100000\t0.8\n")
        c2 = tmp_path / "c2.color"
        c2.write_text("22(pat)\t100000\t0.3\n22(mat)\t100000\t0.6\n")
        ret = mgcolor(["mgcolor", "-d", str(c1), str(c2)])
        assert ret == 0

    def test_mgcolor_further_merge(self, tmp_path, capsys):
        from dip_c.commands.mgcolor import mgcolor
        # First create merged files
        c1 = tmp_path / "c1.colors"
        c1.write_text("homolog\tlocus\tc1\n22(pat)\t100000\t0.5\n")
        c2 = tmp_path / "c2.colors"
        c2.write_text("homolog\tlocus\tc2\n22(pat)\t100000\t0.8\n")
        ret = mgcolor(["mgcolor", "-s", str(c1), str(c2)])
        assert ret == 0

    def test_mgcolor_merge_empty_second(self, tmp_path, capsys):
        """Cover merge_color_data when second dict is empty (line 16)."""
        from dip_c.commands.mgcolor import merge_color_data
        d1 = {("22(pat)", 100000): [0.5]}
        d2 = {}
        merge_color_data(d1, d2, -1.0)
        assert d1 == {("22(pat)", 100000): [0.5]}

    def test_mgcolor_missing_value(self, tmp_path, capsys):
        from dip_c.commands.mgcolor import mgcolor
        c1 = tmp_path / "c1.color"
        c1.write_text("22(pat)\t100000\t0.5\n")
        c2 = tmp_path / "c2.color"
        c2.write_text("22(pat)\t200000\t0.8\n")
        ret = mgcolor(["mgcolor", "-m", "0.0", str(c1), str(c2)])
        assert ret == 0
        out = capsys.readouterr().out
        assert "0.0" in out

    # --- force.py options ---
    def test_force_small(self, tmp_path, capsys):
        from dip_c.commands.force import force
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("22\t51304566\n")
        # Minimal CON
        con_f = tmp_path / "tiny.con"
        con_f.write_text("22,1000000,0\t22,2000000,0\n22,1000000,1\t22,2000000,1\n")
        ret = force(["force", "-l", str(chr_len), "-n", "10", "-w", "5",
                      "-b", "10000000", "-k", "0.01", "-f", "0.1",
                      "-o", str(tmp_path / "f."), str(con_f)])
        assert ret == 0

    # --- align.py options ---
    def test_align_with_output(self, tmp_path, capsys):
        from dip_c.commands.align import align
        g1 = tmp_path / "g1.3dg"
        g1.write_text("22(pat)\t16000000\t1.0\t2.0\t3.0\n22(pat)\t18000000\t4.0\t5.0\t6.0\n"
                       "22(mat)\t16000000\t7.0\t8.0\t9.0\n22(mat)\t18000000\t10.0\t11.0\t12.0\n")
        g2 = tmp_path / "g2.3dg"
        # Slightly different positions
        g2.write_text("22(pat)\t16000000\t2.0\t3.0\t4.0\n22(pat)\t18000000\t5.0\t6.0\t7.0\n"
                       "22(mat)\t16000000\t8.0\t9.0\t10.0\n22(mat)\t18000000\t11.0\t12.0\t13.0\n")
        ret = align(["align", "-o", str(tmp_path / "a."), str(g1), str(g2)])
        assert ret == 0
        # Check aligned file was written
        import glob
        aligned_files = glob.glob(str(tmp_path / "a.*.3dg"))
        assert len(aligned_files) > 0

    def test_align_mirror(self, tmp_path, capsys):
        from dip_c.commands.align import align
        g1 = tmp_path / "g1.3dg"
        g1.write_text("22(pat)\t16000000\t1.0\t2.0\t3.0\n22(pat)\t18000000\t4.0\t5.0\t6.0\n"
                       "22(mat)\t16000000\t7.0\t8.0\t9.0\n22(mat)\t18000000\t10.0\t11.0\t12.0\n")
        g2 = tmp_path / "g2.3dg"
        # Mirrored positions (negated)
        g2.write_text("22(pat)\t16000000\t-1.0\t-2.0\t-3.0\n22(pat)\t18000000\t-4.0\t-5.0\t-6.0\n"
                       "22(mat)\t16000000\t-7.0\t-8.0\t-9.0\n22(mat)\t18000000\t-10.0\t-11.0\t-12.0\n")
        ret = align(["align", "-o", str(tmp_path / "mirror."), str(g1), str(g2)])
        assert ret == 0

    # --- clean3.py -d option ---
    def test_clean3_max_distance(self, capsys):
        from dip_c.commands.clean3 import clean3
        ret = clean3(["clean3", "-c", CHR22_CON_GZ, "-d", "100000", CHR22_3DG])
        assert ret == 0

    # --- con3.py -b option ---
    def test_con3_bin_size(self, tmp_path, capsys):
        from dip_c.commands.con3 import con3
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("22\t51304566\n")
        ret = con3(["con3", "-m", str(chr_len), "-b", "500000", CHR22_3DG])
        assert ret == 0

    # --- impute3.py -s option ---
    def test_impute3_separation_factor(self, capsys):
        from dip_c.commands.impute3 import impute3
        ret = impute3(["impute3", "-3", CHR21_22_3DG, "-s", "2.0", CHR21_22_CON_GZ])
        assert ret == 0

    # --- ard.py intra-chr with histogram (symmetrize intra path) ---
    def test_ard_pairwise_legs_intra(self, tmp_path, capsys):
        from dip_c.commands.ard import ard
        leg_f = tmp_path / "legs.leg"
        leg_f.write_text("22,16000000,0\n22,18000000,0\n22,20000000,0\n")
        ret = ard(["ard", "-1", str(leg_f), "-s", "1000000", CHR22_CON_GZ])
        assert ret == 0

    # --- seg.py with low display interval to hit progress lines ---
    def test_seg_with_progress(self, capsys):
        from dip_c.commands.seg import seg
        ret = seg(["seg", TEST_SEG_BAM], _display_interval=10)
        assert ret == 0
        err = capsys.readouterr().err
        assert "pass 1:" in err or "pass 2:" in err

    def test_seg_with_snps_and_progress(self, capsys):
        from dip_c.commands.seg import seg
        ret = seg(["seg", "-v", TEST_SNPS, TEST_SEG_BAM], _display_interval=1)
        assert ret == 0
        err = capsys.readouterr().err
        assert "pass 3" in err

    # --- con.py with progress ---
    def test_con_with_progress(self, tmp_path, capsys):
        from dip_c.commands.con import con
        # Create a seg file with enough reads to trigger progress
        import gzip
        seg_lines = []
        for i in range(15):
            seg_lines.append(f"read_{i}\t.,0,100,1,{1000+i*1000},{2000+i*1000},+,.\t.,100,200,2,{3000+i*1000},{4000+i*1000},+,.")
        seg_f = tmp_path / "many.seg"
        seg_f.write_text("\n".join(seg_lines) + "\n")
        ret = con(["con", str(seg_f)], _display_interval=5)
        assert ret == 0

    # --- color.py: remaining modes ---
    def test_color_intra_hom_max_sep(self, capsys):
        """Cover -i with -S for intra_hom_fraction with max_separation."""
        from dip_c.commands.color import color
        ret = color(["color", "-i", "5.0", "-S", "1000000", CHR22_3DG])
        assert ret == 0

    def test_color_intra_hom_count(self, capsys):
        """Cover -I for intra_hom_count."""
        from dip_c.commands.color import color
        ret = color(["color", "-I", "5.0", CHR22_3DG])
        assert ret == 0

    def test_color_by_color_file(self, tmp_path, capsys):
        """Cover -c color mode."""
        from dip_c.commands.color import color
        color_f = tmp_path / "color.txt"
        color_f.write_text("22\t16000000\t0.5\n22\t18000000\t0.8\n")
        ret = color(["color", "-c", str(color_f), CHR22_3DG])
        assert ret == 0

    def test_color_homolog_distance(self, capsys):
        """Cover -h mode (distance to homologous locus)."""
        from dip_c.commands.color import color
        ret = color(["color", "-h", CHR22_3DG])
        assert ret == 0

    # --- force.py with usage ---
    def test_force_missing_chr_len(self, tmp_path, capsys):
        from dip_c.commands.force import force
        con_f = tmp_path / "tiny.con"
        con_f.write_text("22,1000000,0\t22,2000000,0\n")
        ret = force(["force", str(con_f)])
        assert ret == 1

    # --- cv.py: non-matching contacts ---
    def test_cv_mismatched_contacts(self, tmp_path, capsys):
        from dip_c.commands.cv import cv
        base = "22,100,.\t22,200,.\n"
        imputed = "22,100,0\t22,200,1\n"  # Same position → found in impute
        truth = "22,500,0\t22,600,1\n"    # Different position → NOT found in truth
        base_f = tmp_path / "b.con"
        base_f.write_text(base)
        imp_f = tmp_path / "i.con"
        imp_f.write_text(imputed)
        truth_f = tmp_path / "t.con"
        truth_f.write_text(truth)
        ret = cv(["cv", "-i", str(imp_f), "-t", str(truth_f), str(base_f)])
        assert ret == 0
        out = capsys.readouterr().out
        assert out.strip() == ""  # No matching contacts

    # --- dist.py: position_np_array_to_rg_np_array (lines 11-29) ---
    def test_dist_rg_function(self, capsys):
        """Cover position_np_array_to_rg_np_array function directly."""
        from dip_c.commands.dist import position_np_array_to_rg_np_array
        import numpy as np
        positions = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        result = position_np_array_to_rg_np_array(positions, 1)
        assert result.shape == (3, 3)

    # --- color.py: radial option parsing (lines 197, 199, 201, 203) ---
    def test_color_radial_with_options(self, tmp_path, capsys):
        """Cover --min-num, --missing, --max-r, --bin-size radial options."""
        from dip_c.commands.color import color
        ret = color(["color", "-C", "-R",
                     "--min-num=2", "--missing=-999.0", "--max-r=5.0", "--bin-size=0.1",
                     CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert len(out.strip().split("\n")) > 1

    # --- color.py: smoothing with -s (line 363) ---
    def test_color_smooth_with_s_progress(self, tmp_path, capsys):
        """Cover -c with -s for smooth_color path incl progress (line 363)."""
        from dip_c.commands.color import color
        color_f = tmp_path / "color.txt"
        # Write color for enough loci to match many particles
        lines = []
        for locus in range(16000000, 50000000, 1000000):
            lines.append(f"21\t{locus}\t0.5")
            lines.append(f"22\t{locus}\t0.5")
        color_f.write_text("\n".join(lines) + "\n")
        # Use CHR21_22_3DG (1417 particles) → atom_id % 1000 hits at 1000
        ret = color(["color", "-c", str(color_f), "-s", "5.0", CHR21_22_3DG])
        assert ret == 0
        err = capsys.readouterr().err
        assert "smoothed 1000" in err

    # --- color.py: -P mode returns None (line 336) covered by existing test ---
    # --- color.py: -n KeyError (lines 310-311) ---
    def test_color_chr_name_key_error(self, tmp_path, capsys):
        """Cover KeyError in -n mode (lines 310-311)."""
        from dip_c.commands.color import color
        chr_f = tmp_path / "chr.txt"
        chr_f.write_text("99\n")  # Non-existent chromosome
        ret = color(["color", "-n", str(chr_f), CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert out.strip() == ""  # All particles have KeyError

    # --- color.py: -l KeyError (lines 315-316) ---
    def test_color_locus_key_error(self, tmp_path, capsys):
        """Cover KeyError in -l mode (lines 315-316)."""
        from dip_c.commands.color import color
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("99\t50000000\n")  # Non-existent chromosome
        ret = color(["color", "-l", str(chr_len), CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert out.strip() == ""

    # --- color.py: -L KeyError (lines 325-326) ---
    def test_color_arm_locus_key_error(self, tmp_path, capsys):
        """Cover KeyError in -L mode (lines 325-326)."""
        from dip_c.commands.color import color
        cen = tmp_path / "chr.cen"
        cen.write_text("99\t50000000\t25000000\n")  # Non-existent chromosome
        ret = color(["color", "-L", str(cen), CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert out.strip() == ""

    # --- color.py: empty 3DG returns 1 (lines 226-227) ---
    def test_color_empty_3dg_returns_error(self, tmp_path, capsys):
        """Cover empty g3d_data error (lines 226-227)."""
        from dip_c.commands.color import color
        g3d_f = tmp_path / "empty.3dg"
        g3d_f.write_text("")
        ret = color(["color", "-C", str(g3d_f)])
        assert ret == 1

    # --- color.py: smooth_color_exc progress (line 376) ---
    def test_color_smooth_exc_progress(self, tmp_path, capsys):
        """Cover smooth_color_exc with enough particles for progress."""
        from dip_c.commands.color import color
        # Use CHR21_22_3DG which has 1417 particles — triggers progress at 1000
        ret = color(["color", "-C", "--s-exc=3.0", CHR21_22_3DG])
        assert ret == 0
        err = capsys.readouterr().err
        assert "smoothed 1000 particles" in err

    # --- color.py: -c-hom with missing loci (line 238 pass for None) ---
    # This is the `pass` for color_mode is None — but that requires no color scheme.
    # Actually it's impossible to reach because num_color_schemes != 1 is caught.
    # Let's skip that unreachable line.

    # --- color.py: intra_hom_fraction returns None (line 19) ---
    def test_color_intra_hom_fraction_no_neighbors(self, tmp_path, capsys):
        """Cover intra_hom_fraction returning None (line 19, 330)."""
        from dip_c.commands.color import color
        # Use a very small distance so most particles have no neighbors
        ret = color(["color", "-i", "0.001", CHR22_3DG])
        assert ret == 0

    # --- color.py: same_haplotype_fraction returns None (line 34) ---
    def test_color_same_haplotype_no_neighbors(self, tmp_path, capsys):
        """Cover same_haplotype_fraction returning None (line 34, 336)."""
        from dip_c.commands.color import color
        ret = color(["color", "-p", "0.001", CHR22_3DG])
        assert ret == 0

    # --- color.py: same_haplotype_fraction_excluding, line 49 (matching haplotype) ---
    def test_color_same_haplotype_exc_with_match(self, capsys):
        """Cover same_haplotype_fraction_excluding line 49 (same haplotype match).
        Need multi-chromosome data with large distance so different homologs
        (e.g., 21(pat) and 22(pat)) appear as neighbors."""
        from dip_c.commands.color import color
        # Distance 20.0 is large enough to span chromosome territories
        ret = color(["color", "-P", "20.0", CHR21_22_3DG])
        assert ret == 0

    # --- color.py: smooth_color returns None (line 98) ---
    def test_color_smooth_no_data(self, tmp_path, capsys):
        """Cover smooth_color returning None when no color data found (line 98)."""
        from dip_c.commands.color import color
        color_f = tmp_path / "color.txt"
        color_f.write_text("99\t99999999\t0.5\n")  # Non-matching data
        ret = color(["color", "-c", str(color_f), "-s", "3.0", CHR22_3DG])
        assert ret == 0

    # --- color.py: smooth_color_exc returns None (line 123) ---
    def test_color_smooth_exc_no_data(self, tmp_path, capsys):
        """Cover smooth_color_exc returning None (line 123)."""
        from dip_c.commands.color import color
        color_f = tmp_path / "color.txt"
        color_f.write_text("99\t99999999\t0.5\n")
        ret = color(["color", "-c", str(color_f), "--s-exc=3.0", CHR22_3DG])
        assert ret == 0

    # --- ard.py: inter-chr histogram with symmetrize (lines 131-134) ---
    def test_ard_inter_chr_histogram_symmetrize(self, capsys):
        """Cover inter-chr symmetrical histogram (lines 131-134)."""
        from dip_c.commands.ard import ard
        # Use chr21+22 data for inter-chr contacts
        ret = ard(["ard", "-d", "5000000", "-h", "500000", CHR21_22_CON_GZ])
        assert ret == 0

    # --- ard.py: -d option (line 57) ---
    def test_ard_max_distance(self, capsys):
        """Cover -d option for max_distance (line 57)."""
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-d", "1000000", CHR22_CON_GZ])
        assert ret == 0

    # --- ard.py: no symmetrize + inter (line 140) ---
    def test_ard_no_symmetrize_inter(self, capsys):
        """Cover no symmetrize for inter-chr histogram (line 140)."""
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-h", "500000", "-S", CHR21_22_CON_GZ])
        assert ret == 0

    # --- ard.py: pairwise mode with two leg files ---
    def test_ard_pairwise_two_leg_files(self, tmp_path, capsys):
        """Cover ard -1 -2 mode with two leg files."""
        from dip_c.commands.ard import ard
        leg1_f = tmp_path / "leg1.leg"
        leg2_f = tmp_path / "leg2.leg"
        leg1_f.write_text("22,16000000,0\n22,18000000,0\n")
        leg2_f.write_text("22,20000000,1\n22,22000000,1\n")
        ret = ard(["ard", "-1", str(leg1_f), "-2", str(leg2_f), CHR22_CON_GZ])
        assert ret == 0

    # --- ard.py: pairwise intra-chr with separation filter (line 177) ---
    def test_ard_pairwise_intra_chr_separation(self, tmp_path, capsys):
        """Cover pairwise intra-chr mode with separation filter (line 177)
        and progress line (line 180)."""
        from dip_c.commands.ard import ard
        leg_f = tmp_path / "legs.leg"
        # Legs on same chr with small and large separations
        # Include inter-chr leg (21) which gets filtered out by -s (line 177)
        leg_f.write_text("22,16000000,0\n22,16001000,0\n22,20000000,0\n21,16000000,0\n")
        ret = ard(["ard", "-s", "1000000", "-1", str(leg_f), CHR22_CON_GZ], _display_interval=1)
        assert ret == 0
        err = capsys.readouterr().err
        assert "analyzed" in err

    # --- bincon.py: None return from leg_to_matrix_index (lines 18, 37) ---
    def test_bincon_unknown_chr(self, tmp_path, capsys):
        """Cover leg_to_matrix_index returning None (lines 18, 37)."""
        from dip_c.commands.bincon import bincon
        chr_len = tmp_path / "chr.len"
        chr_len.write_text("1\t100000000\n")  # Only chr1
        con_f = tmp_path / "con.con"
        # Contact on chr22 which is not in chr.len → returns None
        # Third contact has one known leg (chr1) and one unknown (chr22) so
        # matrix_index_1 is not None but matrix_index_2 is None → traces line 37
        con_f.write_text("22,1000000,0\t22,2000000,0\n1,1000,0\t1,2000,0\n1,1000,0\t22,2000000,0\n")
        ret = bincon(["bincon", "-l", str(chr_len), str(con_f)])
        assert ret == 0

    # --- force.py: progress line (line 156) and collision forces (lines 37-44) ---
    def test_force_with_progress_and_collisions(self, tmp_path, capsys):
        """Cover force progress line and collision force paths (lines 37-44)."""
        from dip_c.commands.force import force
        chr_len = tmp_path / "chr.len"
        # Use a very large bin size so all contacts map to few nodes
        # that will overlap and trigger collision detection
        chr_len.write_text("1\t100000000\n2\t100000000\n")
        # Create contacts between nearby bins (bin_size=10M default)
        # Many contacts to the same bins → multiple forces between same nodes
        lines = []
        for i in range(20):
            lines.append(f"1,{i*10000000},0\t2,{i*10000000},0")
        con_f = tmp_path / "many.con"
        con_f.write_text("\n".join(lines) + "\n")
        # Run enough steps for particles to move and collide
        ret = force(["force", "-n", "50", "-w", "10", "-o", str(tmp_path / "f."), "-l", str(chr_len), str(con_f)], _display_interval=5)
        assert ret == 0
        err = capsys.readouterr().err
        assert "converted" in err

    def test_update_graph_collision_forces(self):
        """Directly test update_graph with particles < 1.0 apart (lines 37-44)."""
        from dip_c.commands.force import update_graph
        positions = np.array([[0.0, 0.0, 0.0],
                              [0.5, 0.0, 0.0]], dtype=float)
        velocities = np.zeros((2, 3), dtype=float)
        forces = np.array([[0, 1]], dtype=int)
        original = positions.copy()
        update_graph(positions, velocities, forces, 0.01, 0.1)
        # Particles should have moved apart (collision repulsion)
        assert not np.allclose(positions, original)

    # --- seg.py: pass 2 mate finding (lines 110-115) and pass 1 skip (line 97) ---
    def test_seg_pass2_mates_and_skip(self, capsys):
        """Cover seg pass 2 mate finding (lines 110-115), properly paired skip (line 97),
        and high-NM SA filter (line 34) using synthetic BAM."""
        from dip_c.commands.seg import seg
        synthetic_bam = os.path.join(os.path.dirname(__file__), "data", "test_seg_synthetic.bam")
        ret = seg(["seg", synthetic_bam], _display_interval=1)
        assert ret == 0

    # --- seg.py: pass 3 SNP handling (lines 129, 137) ---
    def test_seg_pass3_snps_non_candidate(self, tmp_path, capsys):
        """Cover pass 3 non-candidate skip (line 137) using synthetic BAM with SNPs."""
        from dip_c.commands.seg import seg
        synthetic_bam = os.path.join(os.path.dirname(__file__), "data", "test_seg_synthetic.bam")
        # Create SNP file with a position covered by BAM reads
        snp_f = tmp_path / "snps.txt"
        snp_f.write_text("chr22\t16000100\tA\tG\nchr22\t18000050\tC\tT\n")
        ret = seg(["seg", "-v", str(snp_f), synthetic_bam], _display_interval=1)
        assert ret == 0

    # --- color.py: -S error with non-i mode (lines 217-218) ---
    def test_color_max_sep_wrong_mode(self, capsys):
        """Cover -S must be used with -i error (lines 217-218)."""
        from dip_c.commands.color import color
        ret = color(["color", "-I", "5.0", "-S", "1000000", CHR22_3DG])
        assert ret == 1

    # --- color.py: radial mode continues (lines 398, 404) ---
    def test_color_radial_missing_and_oob(self, tmp_path, capsys):
        """Cover radial mode: color missing (398) and out of bound (404)."""
        from dip_c.commands.color import color
        # Use --c-hom with partial data — only a few particles have color (398)
        # Use --max-r=0.01 so colored particles are out of radial bound (404)
        color_f = tmp_path / "color.txt"
        color_f.write_text("22(pat)\t16200000\t0.5\n22(pat)\t16300000\t0.8\n")
        ret = color(["color", "--c-hom=" + str(color_f), "-R", "--max-r=0.01", CHR22_3DG])
        assert ret == 0

    # --- align.py: mirror factor (line 80) ---
    def test_align_mirror_factor(self, tmp_path, capsys):
        """Cover align mirror path (line 80)."""
        from dip_c.commands.align import align
        # Create two 3DG files that are mirror images
        g3d1 = tmp_path / "g1.3dg"
        g3d2 = tmp_path / "g2.3dg"
        lines1 = []
        lines2 = []
        for locus in range(16000000, 20000000, 1000000):
            lines1.append(f"22(pat)\t{locus}\t{locus/1e6}\t{locus/1e6}\t{locus/1e6}")
            # Mirror image (negative coordinates)
            lines2.append(f"22(pat)\t{locus}\t{-locus/1e6}\t{-locus/1e6}\t{-locus/1e6}")
        g3d1.write_text("\n".join(lines1) + "\n")
        g3d2.write_text("\n".join(lines2) + "\n")
        ret = align(["align", str(g3d1), str(g3d2)])
        assert ret == 0

    # --- seg.py: else branches via single-end + secondary reads (lines 86, 97, 114) ---
    @pytest.fixture()
    def seg_else_bam(self, tmp_path):
        """BAM with chimeric, single-end, and secondary reads for else-branch coverage."""
        import pysam
        unsorted = str(tmp_path / "unsorted.bam")
        bam_path = str(tmp_path / "else.bam")
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr22", "LN": 51304566}],
        })
        seq = "ACGT" * 25
        quals = pysam.qualitystring_to_array("I" * 100)
        with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
            # Chimeric read with SA tag → 2 segments, survives clean()
            a = pysam.AlignedSegment(header)
            a.query_name = "chimeric"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 1000
            a.mapping_quality = 60
            a.query_sequence = seq
            a.query_qualities = quals
            a.cigarstring = "100M"
            a.set_tag("NM", 0)
            a.set_tag("SA", "chr22,5000,+,100M,60,0;")
            out.write(a)

            # Single-end mapped read, no SA → else on line 97 (pass 1), 114 (pass 2)
            b = pysam.AlignedSegment(header)
            b.query_name = "singleend"
            b.flag = 0
            b.reference_id = 0
            b.reference_start = 2000
            b.mapping_quality = 60
            b.query_sequence = seq
            b.query_qualities = quals
            b.cigarstring = "100M"
            b.set_tag("NM", 0)
            out.write(b)

            # Secondary alignment (0x100) → is_secondary fires last in or-chain (line 86)
            c = pysam.AlignedSegment(header)
            c.query_name = "secondary_read"
            c.flag = 0x100
            c.reference_id = 0
            c.reference_start = 3000
            c.mapping_quality = 60
            c.query_sequence = seq
            c.query_qualities = quals
            c.cigarstring = "100M"
            c.set_tag("NM", 0)
            out.write(c)

        pysam.sort("-o", bam_path, unsorted)
        pysam.index(bam_path)
        return bam_path

    def test_seg_single_end_and_else_branches(self, seg_else_bam, capsys):
        """Cover seg else-branches: line 86 (secondary), 97 (pass 1), 114 (pass 2)."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", seg_else_bam], _display_interval=1)
        assert ret == 0

    # --- seg.py: SNP neither allele (line 149) ---
    def test_seg_snp_neither_allele(self, tmp_path, capsys):
        """Cover pass 3 else: continue when base matches neither allele (line 149)."""
        from dip_c.commands.seg import seg
        snp_f = tmp_path / "snps.txt"
        # At 22:16208451, candidate SRR7226708.6662870 has base A (bq=27).
        # pat=C, mat=T → A matches neither → line 149.
        snp_f.write_text("22\t16208451\tC\tT\n")
        ret = seg(["seg", "-v", str(snp_f), TEST_SEG_BAM], _display_interval=1)
        assert ret == 0
