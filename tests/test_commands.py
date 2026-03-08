"""Integration tests for dip-c command modules.

Three test data tiers:
  1. Chr22 subset (fast, code coverage): 693 particles, 2,143/11,114 contacts
  2. Chr21+22 subset (CI, algorithm validation): 1,417 particles, 11,382/33,849 contacts
  3. Full genome (local only): 57,864 particles, 249,817/1,126,475 contacts

Test data: GM12878 Cell 11 (GSM3271358) from Tan et al. Science 2018.

Run tiers:
  pytest tests/test_commands.py -k "not FullGenome"    # CI (tiers 1+2, ~30s)
  pytest tests/test_commands.py                         # local (all tiers, ~15 min)
"""

import gzip
import math
import os
import sys
import tempfile

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Paths to fixture data
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# Chr22 subset (fast tests)
CHR22_G3D = os.path.join(DATA_DIR, "gm12878_11.impute.chr22.3dg.txt")
CHR22_IMPUTE_CON = os.path.join(DATA_DIR, "gm12878_11.impute.chr22.con.txt.gz")
CHR22_CLEAN_CON = os.path.join(DATA_DIR, "gm12878_11.clean.chr22.con.txt.gz")

# Chr21+22 subset (CI algorithm tests — has inter-chr contacts for impute)
CHR21_22_G3D = os.path.join(DATA_DIR, "gm12878_11.impute.chr21_22.3dg.txt")
CHR21_22_IMPUTE_CON = os.path.join(DATA_DIR, "gm12878_11.impute.chr21_22.con.txt.gz")
CHR21_22_CLEAN_CON = os.path.join(DATA_DIR, "gm12878_11.clean.chr21_22.con.txt.gz")

# Full genome (local-only algorithm tests)
FULL_G3D = os.path.join(DATA_DIR, "gm12878_11.impute.3dg.txt")
FULL_IMPUTE_CON = os.path.join(DATA_DIR, "gm12878_11.impute.con.txt.gz")
FULL_CLEAN_CON = os.path.join(DATA_DIR, "gm12878_11.clean.con.txt.gz")

NUM_CHR22_PARTICLES = 693
NUM_CHR21_22_PARTICLES = 1417
NUM_FULL_PARTICLES = 57864
NUM_CHR22_IMPUTE_CONS = 2143
NUM_CHR21_22_IMPUTE_CONS = 11382
NUM_CHR21_22_CLEAN_CONS = 33849
NUM_FULL_IMPUTE_CONS = 249817
NUM_FULL_CLEAN_CONS = 1126475


# ===================================================================
# FAST TESTS — chr22 subset (code coverage)
# ===================================================================

# ===== info =====

class TestInfoCommand:
    def test_info_reports_contact_count(self, capsys):
        from dip_c.commands.info import info
        ret = info(["info", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        assert "2143 contacts" in out

    def test_info_reports_intra_chr_and_phased(self, capsys):
        from dip_c.commands.info import info
        info(["info", CHR22_IMPUTE_CON])
        out = capsys.readouterr().out
        assert "% intra-chromosomal" in out
        assert "% legs phased" in out

    def test_info_multiple_files(self, capsys):
        from dip_c.commands.info import info
        ret = info(["info", CHR22_IMPUTE_CON, CHR22_CLEAN_CON])
        assert ret == 0
        out = capsys.readouterr().out
        assert "2143" in out
        assert "11114" in out

    def test_info_no_args_returns_error(self, capsys):
        from dip_c.commands.info import info
        assert info(["info"]) == 1


# ===== dedup =====

class TestDedupCommand:
    def test_dedup_produces_output(self, capsys):
        from dip_c.commands.dedup import dedup
        ret = dedup(["dedup", CHR22_CLEAN_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        assert "\t" in lines[0]

    def test_dedup_reduces_contacts(self, capsys):
        from dip_c.commands.dedup import dedup
        dedup(["dedup", CHR22_CLEAN_CON])
        out = capsys.readouterr().out
        output_count = len(out.strip().split("\n"))
        assert output_count <= 11114

    def test_dedup_custom_params(self, capsys):
        from dip_c.commands.dedup import dedup
        ret = dedup(["dedup", "-s", "500", "-d", "500", CHR22_CLEAN_CON])
        assert ret == 0

    def test_dedup_no_args_returns_error(self, capsys):
        from dip_c.commands.dedup import dedup
        assert dedup(["dedup"]) == 1


# ===== clean =====

class TestCleanCommand:
    def test_clean_produces_output(self, capsys):
        from dip_c.commands.clean import clean
        ret = clean(["clean", CHR22_CLEAN_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0

    def test_clean_output_count(self, capsys):
        from dip_c.commands.clean import clean
        clean(["clean", CHR22_CLEAN_CON])
        out = capsys.readouterr().out
        output_count = len(out.strip().split("\n"))
        assert output_count > 0
        assert output_count <= 11114

    def test_clean_test_mode(self, capsys):
        from dip_c.commands.clean import clean
        ret = clean(["clean", "-t", CHR22_CLEAN_CON])
        assert ret == 0
        out = capsys.readouterr().out
        assert out.strip() == ""

    def test_clean_no_args_returns_error(self, capsys):
        from dip_c.commands.clean import clean
        assert clean(["clean"]) == 1


# ===== reg =====

class TestRegCommand:
    def test_reg_preset_hf(self, capsys):
        from dip_c.commands.reg import reg
        ret = reg(["reg", "-p", "hf", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0

    def test_reg_preset_hm(self, capsys):
        from dip_c.commands.reg import reg
        ret = reg(["reg", "-p", "hm", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_CHR22_IMPUTE_CONS

    def test_reg_unknown_preset_returns_error(self, capsys):
        from dip_c.commands.reg import reg
        assert reg(["reg", "-p", "xx", CHR22_IMPUTE_CON]) == 1

    def test_reg_no_args_returns_error(self, capsys):
        from dip_c.commands.reg import reg
        assert reg(["reg"]) == 1


# ===== impute (chr22) =====

class TestImputeCommandChr22:
    def test_impute_single_chr_produces_output(self, capsys):
        """Impute on single-chromosome data (all contacts are intra-chr/phased)."""
        from dip_c.commands.impute import impute
        ret = impute(["impute", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0

    def test_impute_no_args_returns_error(self, capsys):
        from dip_c.commands.impute import impute
        assert impute(["impute"]) == 1


# ===== reg3 =====

class TestReg3Command:
    def test_reg3_preset_hf(self, capsys):
        from dip_c.commands.reg3 import reg3
        ret = reg3(["reg3", "-p", "hf", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_CHR22_PARTICLES

    def test_reg3_preset_hm(self, capsys):
        from dip_c.commands.reg3 import reg3
        ret = reg3(["reg3", "-p", "hm", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_CHR22_PARTICLES

    def test_reg3_unknown_preset_returns_error(self, capsys):
        from dip_c.commands.reg3 import reg3
        assert reg3(["reg3", "-p", "zz", CHR22_G3D]) == 1

    def test_reg3_no_regions_returns_error(self, capsys):
        from dip_c.commands.reg3 import reg3
        assert reg3(["reg3"]) == 1


# ===== exp =====

class TestExpCommand:
    def test_exp_default_expansion(self, capsys):
        from dip_c.commands.exp import exp
        ret = exp(["exp", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_CHR22_PARTICLES
        parts = lines[0].split("\t")
        assert len(parts) == 5

    def test_exp_centers_only(self, capsys):
        from dip_c.commands.exp import exp
        ret = exp(["exp", "-c", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == 2  # 22(pat) and 22(mat)

    def test_exp_custom_factor(self, capsys):
        from dip_c.commands.exp import exp
        ret = exp(["exp", "-f", "5.0", CHR22_G3D])
        assert ret == 0

    def test_exp_no_args_returns_error(self, capsys):
        from dip_c.commands.exp import exp
        assert exp(["exp"]) == 1


# ===== clean3 =====

class TestClean3Command:
    def test_clean3_with_contacts(self, capsys):
        from dip_c.commands.clean3 import clean3
        ret = clean3(["clean3", "-c", CHR22_IMPUTE_CON, CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert 0 < len(lines) < NUM_CHR22_PARTICLES

    def test_clean3_custom_quantile(self, capsys):
        from dip_c.commands.clean3 import clean3
        ret = clean3(["clean3", "-c", CHR22_IMPUTE_CON, "-q", "0.2", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) < NUM_CHR22_PARTICLES

    def test_clean3_missing_con_returns_error(self, capsys):
        from dip_c.commands.clean3 import clean3
        assert clean3(["clean3", CHR22_G3D]) == 1

    def test_clean3_no_args_returns_error(self, capsys):
        from dip_c.commands.clean3 import clean3
        assert clean3(["clean3"]) == 1


# ===== con3 =====

class TestCon3Command:
    def test_con3_generates_contacts(self, capsys):
        from dip_c.commands.con3 import con3
        ret = con3(["con3", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        assert "\t" in lines[0]
        parts = lines[0].split("\t")
        assert len(parts) == 2

    def test_con3_custom_distance(self, capsys):
        from dip_c.commands.con3 import con3
        ret = con3(["con3", "-d", "1.0", CHR22_G3D])
        assert ret == 0
        out_small = capsys.readouterr().out

        ret = con3(["con3", "-d", "5.0", CHR22_G3D])
        assert ret == 0
        out_large = capsys.readouterr().out

        n_small = len(out_small.strip().split("\n"))
        n_large = len(out_large.strip().split("\n"))
        assert n_small < n_large

    def test_con3_no_args_returns_error(self, capsys):
        from dip_c.commands.con3 import con3
        assert con3(["con3"]) == 1


# ===== dist =====

class TestDistCommand:
    def test_dist_produces_output(self, capsys):
        from dip_c.commands.dist import dist
        ret = dist(["dist", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        parts = lines[0].split("\t")
        assert len(parts) == 5

    def test_dist_output_is_deterministic(self, capsys):
        from dip_c.commands.dist import dist
        dist(["dist", CHR22_G3D])
        out1 = capsys.readouterr().out
        dist(["dist", CHR22_G3D])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_dist_values_make_sense(self, capsys):
        from dip_c.commands.dist import dist
        dist(["dist", CHR22_G3D])
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        parts = lines[0].split("\t")
        hom_name = parts[0]
        separation = float(parts[1])
        num_pairs = int(parts[2])
        mean_dist = float(parts[3])
        rms_dist = float(parts[4])
        assert hom_name in ("22(pat)", "22(mat)")
        assert separation >= 0
        assert num_pairs > 0
        assert mean_dist > 0
        assert rms_dist >= mean_dist

    def test_dist_no_args_returns_error(self, capsys):
        from dip_c.commands.dist import dist
        assert dist(["dist"]) == 1


# ===== color =====

class TestColorCommand:
    def test_color_by_chromosome_name(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-n", "hg19.chr.txt", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        parts = lines[0].split("\t")
        assert len(parts) == 3

    def test_color_by_locus_fraction(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-l", "hg19.chr.len", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        for line in lines[:10]:
            parts = line.split("\t")
            val = float(parts[2])
            assert 0.0 <= val <= 1.0

    def test_color_by_homolog_distance(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-h", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        for line in lines[:10]:
            val = float(line.split("\t")[2])
            assert val >= 0.0

    def test_color_by_center_of_mass_distance(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-C", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_CHR22_PARTICLES

    def test_color_is_deterministic(self, capsys):
        from dip_c.commands.color import color
        color(["color", "-l", "hg19.chr.len", CHR22_G3D])
        out1 = capsys.readouterr().out
        color(["color", "-l", "hg19.chr.len", CHR22_G3D])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_color_no_scheme_returns_error(self, capsys):
        from dip_c.commands.color import color
        assert color(["color", CHR22_G3D]) == 1

    def test_color_no_args_returns_error(self, capsys):
        from dip_c.commands.color import color
        assert color(["color"]) == 1


# ===== bincon =====

class TestBinconCommand:
    def test_bincon_info_mode(self, capsys):
        from dip_c.commands.bincon import bincon
        ret = bincon(["bincon", "-l", "hg19.chr.len", "-i", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        parts = lines[0].split("\t")
        assert len(parts) == 2

    def test_bincon_matrix_mode(self, capsys):
        from dip_c.commands.bincon import bincon
        ret = bincon(["bincon", "-l", "hg19.chr.len", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        ncols = len(lines[0].split("\t"))
        assert ncols == len(lines)

    def test_bincon_merge_haplotypes(self, capsys):
        from dip_c.commands.bincon import bincon
        ret = bincon(["bincon", "-l", "hg19.chr.len", "-H", "-i", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines_merged = out.strip().split("\n")

        bincon(["bincon", "-l", "hg19.chr.len", "-i", CHR22_IMPUTE_CON])
        out = capsys.readouterr().out
        lines_unmerged = out.strip().split("\n")
        assert len(lines_merged) < len(lines_unmerged)

    def test_bincon_requires_chr_len(self, capsys):
        from dip_c.commands.bincon import bincon
        assert bincon(["bincon", CHR22_IMPUTE_CON]) == 1

    def test_bincon_no_args_returns_error(self, capsys):
        from dip_c.commands.bincon import bincon
        assert bincon(["bincon"]) == 1


# ===== vis =====

class TestVisCommand:
    """vis uses pdbx which may have API changes (appendAttribute vs append_attribute)."""

    def _has_pdbx_compat(self):
        try:
            from pdbx import DataCategory
            dc = DataCategory("test")
            dc.appendAttribute("x")
            return True
        except AttributeError:
            return False

    def test_vis_produces_mmcif(self, capsys):
        if not self._has_pdbx_compat():
            pytest.skip("pdbx API incompatible (appendAttribute renamed to append_attribute)")
        from dip_c.commands.vis import vis
        ret = vis(["vis", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        assert "data_myblock" in out
        assert "atom_site" in out
        assert "HETATM" in out

    def test_vis_connect_adjacent(self, capsys):
        if not self._has_pdbx_compat():
            pytest.skip("pdbx API incompatible")
        from dip_c.commands.vis import vis
        ret = vis(["vis", "-a", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        assert "struct_conn" in out

    def test_vis_no_args_returns_error(self, capsys):
        from dip_c.commands.vis import vis
        assert vis(["vis"]) == 1


# ===== impute3 =====

class TestImpute3Command:
    def test_impute3_produces_output(self, capsys):
        from dip_c.commands.impute3 import impute3
        ret = impute3(["impute3", "-3", CHR22_G3D, CHR22_CLEAN_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0

    def test_impute3_missing_3dg_returns_error(self, capsys):
        from dip_c.commands.impute3 import impute3
        assert impute3(["impute3", CHR22_CLEAN_CON]) == 1

    def test_impute3_no_args_returns_error(self, capsys):
        from dip_c.commands.impute3 import impute3
        assert impute3(["impute3"]) == 1


# ===== color2 =====

class TestColor2Command:
    def test_color2_with_generated_color(self, capsys, tmp_path):
        from dip_c.commands.color import color
        color(["color", "-C", CHR22_G3D])
        color_output = capsys.readouterr().out

        color_file = tmp_path / "colors.txt"
        color_file.write_text(color_output)

        from dip_c.commands.color2 import color2
        ret = color2(["color2", "-c", str(color_file), CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0

    def test_color2_no_args_returns_error(self, capsys):
        from dip_c.commands.color2 import color2
        assert color2(["color2"]) == 1


# ===== ard =====

class TestArdCommand:
    def test_ard_count_mode_intra_chr(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-n", "-s", "1000000", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        for line in lines[:5]:
            int(line.strip())

    def test_ard_count_mode_normalized(self, capsys):
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-n", "-t", "-s", "1000000", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        for line in lines[:5]:
            val = float(line.strip())
            assert 0.0 <= val

    def test_ard_inter_chr_mode(self, capsys):
        """Default inter-chr mode (no -s flag)."""
        from dip_c.commands.ard import ard
        ret = ard(["ard", "-n", CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0

    def test_ard_no_args_returns_error(self, capsys):
        from dip_c.commands.ard import ard
        assert ard(["ard"]) == 1


# ===== align =====

class TestAlignCommand:
    def test_align_two_identical_structures(self, capsys):
        from dip_c.commands.align import align
        ret = align(["align", CHR22_G3D, CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_CHR22_PARTICLES
        for line in lines[:10]:
            parts = line.split("\t")
            rmsd_val = float(parts[2])
            assert rmsd_val == pytest.approx(0.0, abs=1e-6)

    def test_align_needs_two_files(self, capsys):
        from dip_c.commands.align import align
        assert align(["align", CHR22_G3D]) == 1

    def test_align_no_args_returns_error(self, capsys):
        from dip_c.commands.align import align
        assert align(["align"]) == 1


# ===== rg =====

class TestRgCommand:
    def test_rg_writes_output_files(self, tmp_path):
        from dip_c.commands.rg import rg
        prefix = str(tmp_path / "test.")
        ret = rg(["rg", "-o", prefix, "-d", CHR22_G3D])
        assert ret == 0

        for hom in ["22(pat)", "22(mat)"]:
            rg_file = prefix + hom + ".rg"
            loc_file = prefix + hom + ".loc"
            assert os.path.isfile(rg_file), f"Missing {rg_file}"
            assert os.path.isfile(loc_file), f"Missing {loc_file}"

    def test_rg_distance_mode(self, tmp_path):
        from dip_c.commands.rg import rg
        prefix = str(tmp_path / "test.")
        ret = rg(["rg", "-o", prefix, "-d", CHR22_G3D])
        assert ret == 0

        rg_file = prefix + "22(pat).rg"
        matrix = np.loadtxt(rg_file, delimiter="\t")
        assert matrix.shape[0] == matrix.shape[1]
        np.testing.assert_array_almost_equal(matrix, matrix.T, decimal=6)
        np.testing.assert_array_almost_equal(np.diag(matrix), 0.0, decimal=6)

    def test_rg_no_args_returns_error(self, capsys):
        from dip_c.commands.rg import rg
        assert rg(["rg"]) == 1


# ===== tad =====

class TestTadCommand:
    def test_tad_from_rg_output(self, tmp_path, capsys):
        from dip_c.commands.rg import rg
        prefix = str(tmp_path / "test.")
        rg(["rg", "-o", prefix, CHR22_G3D])
        capsys.readouterr()

        rg_file = prefix + "22(pat).rg"
        loc_file = prefix + "22(pat).loc"

        from dip_c.commands.tad import tad
        ret = tad(["tad", "-l", loc_file, rg_file])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        loci = np.loadtxt(loc_file, dtype=int, delimiter="\t")
        assert len(lines) == len(loci) - 1
        parts = lines[0].split("\t")
        assert len(parts) == 2

    def test_tad_requires_loc_file(self, capsys, tmp_path):
        rg_file = tmp_path / "dummy.rg"
        np.savetxt(str(rg_file), np.eye(3), delimiter="\t")
        from dip_c.commands.tad import tad
        assert tad(["tad", str(rg_file)]) == 1

    def test_tad_no_args_returns_error(self, capsys):
        from dip_c.commands.tad import tad
        assert tad(["tad"]) == 1


# ===== force (chr22) =====

class TestForceCommandChr22:
    def test_force_runs_and_produces_output(self, tmp_path, capsys):
        from dip_c.commands.force import force
        prefix = str(tmp_path / "force.")
        ret = force(["force", "-l", "hg19.chr.len", "-o", prefix,
                      "-n", "10", "-w", "10", "-b", "50000000",
                      CHR22_IMPUTE_CON])
        assert ret == 0
        output_file = prefix + "10.3dg"
        assert os.path.isfile(output_file)
        with open(output_file) as f:
            lines = f.readlines()
        assert len(lines) > 0
        parts = lines[0].strip().split("\t")
        assert len(parts) == 5

    def test_force_requires_chr_len(self, capsys):
        from dip_c.commands.force import force
        assert force(["force", CHR22_IMPUTE_CON]) == 1

    def test_force_no_args_returns_error(self, capsys):
        from dip_c.commands.force import force
        assert force(["force"]) == 1


# ===== mgcolor =====

class TestMgcolorCommand:
    def test_mgcolor_merges_two_files(self, capsys, tmp_path):
        from dip_c.commands.color import color
        color(["color", "-C", CHR22_G3D])
        color1_text = capsys.readouterr().out
        color(["color", "-h", CHR22_G3D])
        color2_text = capsys.readouterr().out

        f1 = tmp_path / "c1.txt"
        f2 = tmp_path / "c2.txt"
        f1.write_text(color1_text)
        f2.write_text(color2_text)

        from dip_c.commands.mgcolor import mgcolor
        ret = mgcolor(["mgcolor", str(f1), str(f2)])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert "homolog" in lines[0] or "locus" in lines[0]
        parts = lines[1].split("\t")
        assert len(parts) == 4

    def test_mgcolor_no_args_returns_error(self, capsys):
        from dip_c.commands.mgcolor import mgcolor
        assert mgcolor(["mgcolor"]) == 1


# ===== pos =====

class TestPosCommand:
    def test_pos_converts_legs_to_positions(self, capsys, tmp_path):
        leg_file = tmp_path / "test.leg"
        leg_file.write_text("22,20000000,0\n22,30000000,1\n")

        from dip_c.commands.pos import pos
        ret = pos(["pos", "-l", str(leg_file), CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == 2
        for line in lines:
            if line != "None":
                parts = line.split("\t")
                assert len(parts) == 3

    def test_pos_requires_leg_file(self, capsys):
        from dip_c.commands.pos import pos
        assert pos(["pos", CHR22_G3D]) == 1

    def test_pos_no_args_returns_error(self, capsys):
        from dip_c.commands.pos import pos
        assert pos(["pos"]) == 1


# ===== pd =====

class TestPdCommand:
    def test_pd_pairwise_distances(self, capsys, tmp_path):
        leg_file = tmp_path / "legs.leg"
        leg_file.write_text("22,20000000,0\n22,30000000,1\n22,40000000,0\n")

        from dip_c.commands.pd import pd
        ret = pd(["pd", "-1", str(leg_file), CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == 3
        for line in lines:
            parts = line.split("\t")
            assert len(parts) == 3

    def test_pd_requires_leg_file(self, capsys):
        from dip_c.commands.pd import pd
        assert pd(["pd", CHR22_G3D]) == 1

    def test_pd_no_args_returns_error(self, capsys):
        from dip_c.commands.pd import pd
        assert pd(["pd"]) == 1


# ===== cv =====

class TestCvCommand:
    def test_cv_with_same_files(self, capsys):
        from dip_c.commands.cv import cv
        ret = cv(["cv", "-i", CHR22_IMPUTE_CON, "-t", CHR22_IMPUTE_CON,
                   CHR22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        for line in lines[:10]:
            assert "wrong" not in line

    def test_cv_requires_both_files(self, capsys):
        from dip_c.commands.cv import cv
        assert cv(["cv", "-i", CHR22_IMPUTE_CON, CHR22_IMPUTE_CON]) == 1
        assert cv(["cv", "-t", CHR22_IMPUTE_CON, CHR22_IMPUTE_CON]) == 1

    def test_cv_no_args_returns_error(self, capsys):
        from dip_c.commands.cv import cv
        assert cv(["cv"]) == 1


# ===== con3 matrix mode =====

class TestCon3MatrixMode:
    def test_con3_matrix_mode(self, capsys):
        from dip_c.commands.con3 import con3
        ret = con3(["con3", "-m", "hg19.chr.len", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        ncols = len(lines[0].split("\t"))
        assert ncols == len(lines)

    def test_con3_matrix_info_mode(self, capsys):
        from dip_c.commands.con3 import con3
        ret = con3(["con3", "-m", "hg19.chr.len", "-i", CHR22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0


# ===== Deterministic computation tests (chr22) =====

class TestDeterministicComputations:
    """Verify that commands produce identical output on repeated runs."""

    def test_info_deterministic(self, capsys):
        from dip_c.commands.info import info
        info(["info", CHR22_IMPUTE_CON])
        out1 = capsys.readouterr().out
        info(["info", CHR22_IMPUTE_CON])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_dedup_deterministic(self, capsys):
        from dip_c.commands.dedup import dedup
        dedup(["dedup", CHR22_CLEAN_CON])
        out1 = capsys.readouterr().out
        dedup(["dedup", CHR22_CLEAN_CON])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_clean_deterministic(self, capsys):
        from dip_c.commands.clean import clean
        clean(["clean", CHR22_CLEAN_CON])
        out1 = capsys.readouterr().out
        clean(["clean", CHR22_CLEAN_CON])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_reg_deterministic(self, capsys):
        from dip_c.commands.reg import reg
        reg(["reg", "-p", "hm", CHR22_IMPUTE_CON])
        out1 = capsys.readouterr().out
        reg(["reg", "-p", "hm", CHR22_IMPUTE_CON])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_exp_deterministic(self, capsys):
        from dip_c.commands.exp import exp
        exp(["exp", CHR22_G3D])
        out1 = capsys.readouterr().out
        exp(["exp", CHR22_G3D])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_con3_deterministic(self, capsys):
        from dip_c.commands.con3 import con3
        con3(["con3", CHR22_G3D])
        out1 = capsys.readouterr().out
        con3(["con3", CHR22_G3D])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_dist_deterministic(self, capsys):
        from dip_c.commands.dist import dist
        dist(["dist", CHR22_G3D])
        out1 = capsys.readouterr().out
        dist(["dist", CHR22_G3D])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_bincon_deterministic(self, capsys):
        from dip_c.commands.bincon import bincon
        bincon(["bincon", "-l", "hg19.chr.len", CHR22_IMPUTE_CON])
        out1 = capsys.readouterr().out
        bincon(["bincon", "-l", "hg19.chr.len", CHR22_IMPUTE_CON])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_vis_deterministic(self, capsys):
        try:
            from pdbx import DataCategory
            DataCategory("test").appendAttribute("x")
        except AttributeError:
            pytest.skip("pdbx API incompatible")
        from dip_c.commands.vis import vis
        vis(["vis", CHR22_G3D])
        out1 = capsys.readouterr().out
        vis(["vis", CHR22_G3D])
        out2 = capsys.readouterr().out
        assert out1 == out2

    def test_align_deterministic(self, capsys):
        from dip_c.commands.align import align
        align(["align", CHR22_G3D, CHR22_G3D])
        out1 = capsys.readouterr().out
        align(["align", CHR22_G3D, CHR22_G3D])
        out2 = capsys.readouterr().out
        assert out1 == out2


# ===== Seg command =====

SEG_BAM = os.path.join(DATA_DIR, "test_seg.bam")
SEG_SNPS = os.path.join(DATA_DIR, "test_snps.txt")


class TestSegCommand:
    def test_seg_no_args_returns_usage(self, capsys):
        """seg with no arguments should print usage and return 1."""
        from dip_c.commands.seg import seg
        ret = seg(["seg"])
        assert ret == 1
        err = capsys.readouterr().err
        assert "Usage:" in err

    def test_seg_basic_extraction(self, capsys):
        """seg extracts 12 reads with 44 segments from test BAM."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", SEG_BAM])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == 12
        # Count total segments (tab-separated, first field is read name)
        total_segs = sum(len(line.split("\t")) - 1 for line in lines)
        assert total_segs == 44

    def test_seg_unphased_output_matches(self, capsys):
        """seg output matches expected unphased reference."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", SEG_BAM])
        assert ret == 0
        out = capsys.readouterr().out
        with open(os.path.join(DATA_DIR, "test_seg.seg")) as f:
            expected = f.read()
        assert out == expected

    def test_seg_phased_output(self, capsys):
        """seg -v phases segments using SNP file."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", "-v", SEG_SNPS, SEG_BAM])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == 12
        # Count phased segments (last field is 0 or 1, not .)
        phased = 0
        total = 0
        for line in lines:
            for seg_field in line.split("\t")[1:]:
                total += 1
                hap = seg_field.split(",")[-1]
                if hap in ("0", "1"):
                    phased += 1
        assert phased == 14
        assert total == 44

    def test_seg_phased_output_matches(self, capsys):
        """seg -v output matches expected phased reference."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", "-v", SEG_SNPS, SEG_BAM])
        assert ret == 0
        out = capsys.readouterr().out
        with open(os.path.join(DATA_DIR, "test_seg_phased.seg")) as f:
            expected = f.read()
        assert out == expected

    def test_seg_chimeric_reads_present(self, capsys):
        """Chimeric reads with SA tags produce multi-segment output."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", SEG_BAM])
        assert ret == 0
        out = capsys.readouterr().out
        # SRR7226708.9938350 has multi-SA — should have 5 segments
        for line in out.strip().split("\n"):
            if "9938350" in line:
                segs = line.split("\t")[1:]
                assert len(segs) == 5
                break
        else:
            pytest.fail("Multi-SA read SRR7226708.9938350 not found in output")

    def test_seg_filtered_reads_excluded(self, capsys):
        """Low MAPQ and high NM reads are filtered out."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", SEG_BAM])
        assert ret == 0
        out = capsys.readouterr().out
        # These reads should NOT appear in output
        assert "2363794" not in out   # low MAPQ
        assert "11247069" not in out  # low MAPQ
        assert "30346024" not in out  # high NM
        assert "20582799" not in out  # high NM

    def test_seg_duplicate_skipped(self, capsys):
        """Duplicate-flagged reads are skipped."""
        from dip_c.commands.seg import seg
        ret = seg(["seg", SEG_BAM])
        assert ret == 0
        err = capsys.readouterr().err
        # The duplicate read should not appear in candidate reads
        # 12 reads survive out of 26 unique read names
        assert "12 candidate reads" in err

    def test_seg_custom_baseq(self, capsys):
        """seg -Q lowers min base quality, phasing more segments."""
        from dip_c.commands.seg import seg
        # Default -Q 20 gives 14 phased segments
        ret = seg(["seg", "-v", SEG_SNPS, SEG_BAM])
        assert ret == 0
        default_out = capsys.readouterr().out
        default_phased = sum(
            1 for line in default_out.strip().split("\n")
            for seg_field in line.split("\t")[1:]
            if seg_field.split(",")[-1] in ("0", "1")
        )
        assert default_phased == 14
        # -Q 10 should phase at least one more (SRR7226708.10101251, bq=16)
        ret = seg(["seg", "-Q", "10", "-v", SEG_SNPS, SEG_BAM])
        assert ret == 0
        lowq_out = capsys.readouterr().out
        lowq_phased = sum(
            1 for line in lowq_out.strip().split("\n")
            for seg_field in line.split("\t")[1:]
            if seg_field.split(",")[-1] in ("0", "1")
        )
        assert lowq_phased > default_phased

    def test_seg_to_con_pipeline(self, capsys):
        """Full seg → con pipeline produces expected contacts."""
        from dip_c.commands.seg import seg
        from dip_c.commands.con import con
        # Run seg with phasing
        ret = seg(["seg", "-v", SEG_SNPS, SEG_BAM])
        assert ret == 0
        seg_out = capsys.readouterr().out
        # Write seg output to temp file, then run con
        seg_file = os.path.join(tempfile.mkdtemp(), "test.seg")
        with open(seg_file, "w") as f:
            f.write(seg_out)
        ret = con(["con", seg_file])
        assert ret == 0
        con_out = capsys.readouterr().out
        lines = con_out.strip().split("\n")
        assert len(lines) == 18
        # Check that phased contacts exist
        assert ",0\t" in con_out or ",1\t" in con_out or ",0\n" in con_out or ",1\n" in con_out


# ===================================================================
# CI TESTS — chr21+22 subset (algorithm validation, ~30s)
# ===================================================================

class TestChr2122Info:
    """Validate info on multi-chromosome subset."""

    def test_info_impute_con(self, capsys):
        from dip_c.commands.info import info
        ret = info(["info", CHR21_22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        assert "11382 contacts" in out

    def test_info_clean_con(self, capsys):
        from dip_c.commands.info import info
        ret = info(["info", CHR21_22_CLEAN_CON])
        assert ret == 0
        out = capsys.readouterr().out
        assert "33849 contacts" in out
        # Should have inter-chromosomal contacts
        assert "intra-chromosomal" in out


class TestChr2122Impute:
    """Test impute with multi-chromosome data (needs inter-chr contacts).

    The chr21+22 clean CON has ~10% legs phased and 10,048 unphased
    inter-chromosomal contacts — enough to exercise all 3 impute passes.
    """

    def test_impute_with_clean_con(self, capsys):
        """Impute unphased contacts from the chr21+22 clean CON."""
        from dip_c.commands.impute import impute
        ret = impute(["impute", CHR21_22_CLEAN_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        # All output contacts should be fully phased
        for line in lines[:20]:
            parts = line.split("\t")
            assert len(parts) == 2
            for part in parts:
                fields = part.split(",")
                assert fields[2] in ("0", "1")

    def test_impute_phased_input(self, capsys):
        """Already-phased contacts have 0 unphased inter-chr — should still succeed."""
        from dip_c.commands.impute import impute
        ret = impute(["impute", CHR21_22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0

    def test_impute_deterministic(self, capsys):
        """Impute should produce identical output on repeated runs."""
        from dip_c.commands.impute import impute
        impute(["impute", CHR21_22_CLEAN_CON])
        out1 = capsys.readouterr().out
        impute(["impute", CHR21_22_CLEAN_CON])
        out2 = capsys.readouterr().out
        assert out1 == out2


class TestChr2122Structure:
    """Test commands on the chr21+22 3D structure (4 homologs)."""

    def test_exp_4_homologs(self, capsys):
        from dip_c.commands.exp import exp
        ret = exp(["exp", "-c", CHR21_22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        # 21(pat), 21(mat), 22(pat), 22(mat) = 4
        assert len(lines) == 4

    def test_color_center_of_mass(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-C", CHR21_22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_CHR21_22_PARTICLES

    def test_color_homolog_distance(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-h", CHR21_22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        for line in lines[:10]:
            val = float(line.split("\t")[2])
            assert val >= 0.0

    def test_dist_4_homologs(self, capsys):
        from dip_c.commands.dist import dist
        ret = dist(["dist", CHR21_22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        hom_names = set(line.split("\t")[0] for line in lines)
        assert len(hom_names) == 4

    def test_align_self(self, capsys):
        from dip_c.commands.align import align
        ret = align(["align", CHR21_22_G3D, CHR21_22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_CHR21_22_PARTICLES
        for line in lines[:10]:
            rmsd_val = float(line.split("\t")[2])
            assert rmsd_val == pytest.approx(0.0, abs=1e-6)

    def test_con3(self, capsys):
        from dip_c.commands.con3 import con3
        ret = con3(["con3", CHR21_22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        # Should produce contacts across the 4 homologs
        assert len(lines) > 100

    def test_clean3(self, capsys):
        from dip_c.commands.clean3 import clean3
        ret = clean3(["clean3", "-c", CHR21_22_IMPUTE_CON, CHR21_22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert 0 < len(lines) < NUM_CHR21_22_PARTICLES

    def test_pos_cross_chr(self, capsys, tmp_path):
        """Position interpolation for loci on different chromosomes."""
        leg_file = tmp_path / "test.leg"
        leg_file.write_text("21,20000000,0\n22,30000000,1\n")

        from dip_c.commands.pos import pos
        ret = pos(["pos", "-l", str(leg_file), CHR21_22_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == 2
        for line in lines:
            assert line != "None"
            parts = line.split("\t")
            assert len(parts) == 3


class TestChr2122Bincon:
    """Test contact binning with multi-chromosome data."""

    def test_bincon_matrix(self, capsys):
        from dip_c.commands.bincon import bincon
        ret = bincon(["bincon", "-l", "hg19.chr.len", CHR21_22_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        ncols = len(lines[0].split("\t"))
        assert ncols == len(lines)
        # Matrix should be larger than chr22-only
        assert ncols > 10


class TestChr2122Force:
    """Test force-directed layout with multi-chromosome contacts."""

    def test_force(self, tmp_path, capsys):
        from dip_c.commands.force import force
        prefix = str(tmp_path / "force.")
        ret = force(["force", "-l", "hg19.chr.len", "-o", prefix,
                      "-n", "10", "-w", "10", "-b", "50000000",
                      CHR21_22_IMPUTE_CON])
        assert ret == 0
        output_file = prefix + "10.3dg"
        assert os.path.isfile(output_file)
        with open(output_file) as f:
            lines = f.readlines()
        assert len(lines) > 0
        parts = lines[0].strip().split("\t")
        assert len(parts) == 5


class TestChr2122Rg:
    """Test radius of gyration with multi-chromosome data."""

    def test_rg_both_chromosomes(self, tmp_path):
        from dip_c.commands.rg import rg
        prefix = str(tmp_path / "test.")
        ret = rg(["rg", "-o", prefix, "-d", CHR21_22_G3D])
        assert ret == 0

        for hom in ["21(pat)", "21(mat)", "22(pat)", "22(mat)"]:
            rg_file = prefix + hom + ".rg"
            loc_file = prefix + hom + ".loc"
            assert os.path.isfile(rg_file), f"Missing {rg_file}"
            assert os.path.isfile(loc_file), f"Missing {loc_file}"

        # Verify distance matrix properties
        matrix = np.loadtxt(prefix + "21(pat).rg", delimiter="\t")
        assert matrix.shape[0] == matrix.shape[1]
        np.testing.assert_array_almost_equal(matrix, matrix.T, decimal=6)
        np.testing.assert_array_almost_equal(np.diag(matrix), 0.0, decimal=6)
        # Chr21 should have more particles than chr22
        assert matrix.shape[0] > 300


# ===================================================================
# FULL-GENOME TESTS — local-only deep validation (skip with -k "not FullGenome")
# ===================================================================

class TestFullGenomeInfo:
    """Validate info on full genome data."""

    def test_info_full_genome(self, capsys):
        from dip_c.commands.info import info
        ret = info(["info", FULL_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        assert "249817 contacts" in out
        assert "intra-chromosomal" in out

    def test_info_clean_con(self, capsys):
        from dip_c.commands.info import info
        ret = info(["info", FULL_CLEAN_CON])
        assert ret == 0
        out = capsys.readouterr().out
        assert "1126475 contacts" in out


class TestFullGenomeImpute:
    """Test impute with full genome (needs inter-chr contacts).

    The clean CON has 9.05% legs phased — real imputation input.
    Using the impute CON (already 100% phased) hits a division-by-zero
    because there are 0 unphased inter-chr contacts.
    """

    def test_impute_full_genome_phased_input(self, capsys):
        """Already-phased contacts have 0 unphased inter-chr — should still succeed."""
        from dip_c.commands.impute import impute
        ret = impute(["impute", FULL_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0

    def test_impute_full_genome_with_clean_con(self, capsys):
        """Impute unphased contacts from the full-genome clean CON file."""
        from dip_c.commands.impute import impute
        ret = impute(["impute", FULL_CLEAN_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        # All output contacts should be fully phased
        for line in lines[:20]:
            parts = line.split("\t")
            assert len(parts) == 2
            for part in parts:
                fields = part.split(",")
                assert fields[2] in ("0", "1")

    def test_impute_deterministic(self, capsys):
        from dip_c.commands.impute import impute
        impute(["impute", FULL_CLEAN_CON])
        out1 = capsys.readouterr().out
        impute(["impute", FULL_CLEAN_CON])
        out2 = capsys.readouterr().out
        assert out1 == out2


class TestFullGenomeStructure:
    """Test commands that operate on the full 3D genome structure."""

    def test_exp_full_genome_46_homologs(self, capsys):
        from dip_c.commands.exp import exp
        ret = exp(["exp", "-c", FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        # GM12878 is female: 22 autosomes + X = 23 chromosomes x 2 haplotypes = 46
        assert len(lines) == 46

    def test_reg3_full_genome_hm(self, capsys):
        from dip_c.commands.reg3 import reg3
        ret = reg3(["reg3", "-p", "hm", FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_FULL_PARTICLES

    def test_color_center_of_mass_full_genome(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-C", FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_FULL_PARTICLES

    def test_color_homolog_distance_full_genome(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-h", FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0
        # All distances should be non-negative
        for line in lines[:20]:
            val = float(line.split("\t")[2])
            assert val >= 0.0

    def test_color_by_locus_fraction_full_genome(self, capsys):
        from dip_c.commands.color import color
        ret = color(["color", "-l", "hg19.chr.len", FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_FULL_PARTICLES
        for line in lines[:20]:
            val = float(line.split("\t")[2])
            assert 0.0 <= val <= 1.0

    def test_dist_full_genome_46_homologs(self, capsys):
        from dip_c.commands.dist import dist
        ret = dist(["dist", FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        hom_names = set(line.split("\t")[0] for line in lines)
        # GM12878 is female: 22 autosomes + X = 23 chromosomes x 2 = 46
        assert len(hom_names) == 46

    def test_align_full_genome_self(self, capsys):
        from dip_c.commands.align import align
        ret = align(["align", FULL_G3D, FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == NUM_FULL_PARTICLES
        # Self-alignment RMSD should be 0
        for line in lines[:10]:
            rmsd_val = float(line.split("\t")[2])
            assert rmsd_val == pytest.approx(0.0, abs=1e-6)

    def test_con3_full_genome(self, capsys):
        from dip_c.commands.con3 import con3
        ret = con3(["con3", FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        # Full genome should produce many inter-chromosomal contacts
        assert len(lines) > NUM_CHR22_PARTICLES  # more contacts than chr22 particles

    def test_clean3_full_genome(self, capsys):
        from dip_c.commands.clean3 import clean3
        ret = clean3(["clean3", "-c", FULL_IMPUTE_CON, FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert 0 < len(lines) < NUM_FULL_PARTICLES

    def test_impute3_full_genome(self, capsys):
        from dip_c.commands.impute3 import impute3
        ret = impute3(["impute3", "-3", FULL_G3D, FULL_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) > 0

    def test_pos_full_genome_cross_chr(self, capsys, tmp_path):
        """Test position interpolation for loci on different chromosomes."""
        leg_file = tmp_path / "test.leg"
        leg_file.write_text("1,50000000,0\n10,30000000,1\n22,20000000,0\n")

        from dip_c.commands.pos import pos
        ret = pos(["pos", "-l", str(leg_file), FULL_G3D])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        assert len(lines) == 3
        # All loci should resolve (not None) in full genome
        for line in lines:
            assert line != "None"
            parts = line.split("\t")
            assert len(parts) == 3


class TestFullGenomeForce:
    """Test force-directed 3D structure generation (needs full genome)."""

    def test_force_full_genome(self, tmp_path, capsys):
        from dip_c.commands.force import force
        prefix = str(tmp_path / "force.")
        ret = force(["force", "-l", "hg19.chr.len", "-o", prefix,
                      "-n", "10", "-w", "10", "-b", "50000000",
                      FULL_IMPUTE_CON])
        assert ret == 0
        output_file = prefix + "10.3dg"
        assert os.path.isfile(output_file)
        with open(output_file) as f:
            lines = f.readlines()
        assert len(lines) > 0
        parts = lines[0].strip().split("\t")
        assert len(parts) == 5
        # Force uses hg19.chr.len which includes Y: 24 chromosomes x 2 = 48
        hom_names = set(line.strip().split("\t")[0] for line in lines)
        assert len(hom_names) == 48


class TestFullGenomeBincon:
    """Test contact binning with full genome data."""

    def test_bincon_full_genome_matrix(self, capsys):
        from dip_c.commands.bincon import bincon
        ret = bincon(["bincon", "-l", "hg19.chr.len", FULL_IMPUTE_CON])
        assert ret == 0
        out = capsys.readouterr().out
        lines = out.strip().split("\n")
        ncols = len(lines[0].split("\t"))
        assert ncols == len(lines)
        # Matrix should be larger than chr22-only
        assert ncols > 10


class TestFullGenomeRg:
    """Test radius of gyration with full genome (restricted to one homolog)."""

    def test_rg_full_genome_single_chr(self, tmp_path):
        from dip_c.commands.rg import rg
        prefix = str(tmp_path / "test.")
        reg_file = tmp_path / "chr1.reg"
        reg_file.write_text("1\t.\t.\t.\n")
        ret = rg(["rg", "-o", prefix, "-r", str(reg_file), "-d", FULL_G3D])
        assert ret == 0

        for hom in ["1(pat)", "1(mat)"]:
            rg_file = prefix + hom + ".rg"
            loc_file = prefix + hom + ".loc"
            assert os.path.isfile(rg_file)
            assert os.path.isfile(loc_file)

        # Verify chr1 distance matrix
        matrix = np.loadtxt(prefix + "1(pat).rg", delimiter="\t")
        assert matrix.shape[0] == matrix.shape[1]
        np.testing.assert_array_almost_equal(matrix, matrix.T, decimal=6)
        np.testing.assert_array_almost_equal(np.diag(matrix), 0.0, decimal=6)
        # Chr1 should have more particles than chr22
        assert matrix.shape[0] > 100
