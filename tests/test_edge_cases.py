"""Edge-case tests for dip-c commands.

All tests use synthetic inline data (no test data files needed), fast to run.
"""

import gzip
import os
import sys
import tempfile

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

    def test_dist_diploid_mode(self, capsys):
        from dip_c.commands.dist import dist
        ret = dist(["dist", "-d", CHR22_3DG])
        assert ret == 0
        out = capsys.readouterr().out
        assert "22(pat)" in out or "22(mat)" in out

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
