import pytest
import numpy as np
from io import StringIO

from dip_c.classes import (
    Haplotypes, is_known_haplotype, haplotype_to_string, string_to_haplotype,
    known_haplotypes, homologous_haplotype, update_haplotype,
    homologous_hom_name, homologous_ref_name_haplotype,
    ref_name_haplotype_to_hom_name, hom_name_to_ref_name_haplotype,
    ref_name_tuple_to_string,
    Leg, string_to_leg, Con, string_to_con, ConData, file_to_con_data,
    LegData, LegList,
    G3dParticle, string_to_g3d_particle, G3dData, file_to_g3d_data,
    Reg, string_to_reg, file_to_reg_list, get_phased_regs,
    Par, ParData,
    DupLeg, DupCon, DupConList, DupConData,
    vote_from_hap_tuples, winning_vote,
    counts_to_hist_num_with_zero, hist_num_to_string, hist_num_to_string_with_zero,
)


# ---------------------------------------------------------------------------
# Haplotype utilities
# ---------------------------------------------------------------------------

class TestHaplotypes:
    def test_haplotype_values(self):
        assert Haplotypes.paternal == 0
        assert Haplotypes.maternal == 1
        assert Haplotypes.unknown == -1
        assert Haplotypes.conflict == -2
        assert Haplotypes.minus_infinity == -3
        assert Haplotypes.infinity == 2

    def test_is_known_haplotype(self):
        assert is_known_haplotype(Haplotypes.paternal) is True
        assert is_known_haplotype(Haplotypes.maternal) is True
        assert is_known_haplotype(Haplotypes.unknown) is False
        assert is_known_haplotype(Haplotypes.conflict) is False

    def test_haplotype_to_string(self):
        assert haplotype_to_string(Haplotypes.paternal) == "0"
        assert haplotype_to_string(Haplotypes.maternal) == "1"
        assert haplotype_to_string(Haplotypes.unknown) == "."
        assert haplotype_to_string(Haplotypes.conflict) == "."

    def test_string_to_haplotype(self):
        assert string_to_haplotype("0") == Haplotypes.paternal
        assert string_to_haplotype("1") == Haplotypes.maternal
        assert string_to_haplotype(".") == Haplotypes.unknown

    def test_known_haplotypes_generator(self):
        haps = list(known_haplotypes())
        assert set(haps) == {Haplotypes.paternal, Haplotypes.maternal}

    def test_homologous_haplotype(self):
        assert homologous_haplotype(Haplotypes.paternal) == Haplotypes.maternal
        assert homologous_haplotype(Haplotypes.maternal) == Haplotypes.paternal
        assert homologous_haplotype(Haplotypes.unknown) is None

    def test_update_haplotype_same(self):
        assert update_haplotype(Haplotypes.paternal, Haplotypes.paternal) == Haplotypes.paternal

    def test_update_haplotype_with_unknown(self):
        assert update_haplotype(Haplotypes.unknown, Haplotypes.maternal) == Haplotypes.maternal
        assert update_haplotype(Haplotypes.paternal, Haplotypes.unknown) == Haplotypes.paternal

    def test_update_haplotype_conflict(self):
        assert update_haplotype(Haplotypes.paternal, Haplotypes.maternal) == Haplotypes.conflict


# ---------------------------------------------------------------------------
# Homolog name conversions
# ---------------------------------------------------------------------------

class TestHaploidNames:
    def test_ref_name_haplotype_to_hom_name(self):
        assert ref_name_haplotype_to_hom_name(("chr1", Haplotypes.paternal)) == "chr1(pat)"
        assert ref_name_haplotype_to_hom_name(("22", Haplotypes.maternal)) == "22(mat)"

    def test_hom_name_to_ref_name_haplotype(self):
        assert hom_name_to_ref_name_haplotype("chr1(pat)") == ("chr1", Haplotypes.paternal)
        assert hom_name_to_ref_name_haplotype("22(mat)") == ("22", Haplotypes.maternal)

    def test_hom_name_roundtrip(self):
        t = ("chr1", Haplotypes.paternal)
        assert hom_name_to_ref_name_haplotype(ref_name_haplotype_to_hom_name(t)) == t

    def test_homologous_ref_name_haplotype(self):
        assert homologous_ref_name_haplotype(("chr1", Haplotypes.paternal)) == ("chr1", Haplotypes.maternal)

    def test_homologous_hom_name(self):
        assert homologous_hom_name("chr1(pat)") == "chr1(mat)"
        assert homologous_hom_name("22(mat)") == "22(pat)"

    def test_ref_name_tuple_to_string(self):
        assert ref_name_tuple_to_string(("chr1", "chr2")) == "chr1,chr2"


# ---------------------------------------------------------------------------
# Leg
# ---------------------------------------------------------------------------

class TestLeg:
    def test_leg_creation(self):
        leg = Leg("chr1", 1000000, Haplotypes.paternal)
        assert leg.get_ref_name() == "chr1"
        assert leg.get_ref_locus() == 1000000
        assert leg.get_haplotype() == Haplotypes.paternal

    def test_leg_equality_and_hash(self):
        leg1 = Leg("chr1", 1000000, Haplotypes.paternal)
        leg2 = Leg("chr1", 1000000, Haplotypes.paternal)
        assert leg1 == leg2
        assert hash(leg1) == hash(leg2)
        assert leg1 != Leg("chr1", 1000000, Haplotypes.maternal)

    def test_leg_ordering(self):
        leg1 = Leg("chr1", 1000000, Haplotypes.paternal)
        leg2 = Leg("chr1", 2000000, Haplotypes.paternal)
        leg3 = Leg("chr2", 500000, Haplotypes.paternal)
        assert leg1 < leg2
        assert leg1 < leg3

    def test_string_to_leg(self):
        leg = string_to_leg("chr1,1000000,0")
        assert leg.get_ref_name() == "chr1"
        assert leg.get_ref_locus() == 1000000
        assert leg.get_haplotype() == Haplotypes.paternal

    def test_leg_to_string(self):
        assert Leg("chr1", 1000000, Haplotypes.paternal).to_string() == "chr1,1000000,0"
        assert Leg("chr2", 2000000, Haplotypes.unknown).to_string() == "chr2,2000000,."

    def test_leg_roundtrip(self):
        s = "chr1,1000000,0"
        assert string_to_leg(s).to_string() == s

    def test_leg_repr(self):
        leg = Leg("chr1", 100, Haplotypes.paternal)
        assert repr(leg) == "chr1,100,0"

    def test_is_phased(self):
        assert Leg("chr1", 0, Haplotypes.paternal).is_phased() is True
        assert Leg("chr1", 0, Haplotypes.unknown).is_phased() is False

    def test_is_conflict(self):
        assert Leg("chr1", 0, Haplotypes.conflict).is_conflict() is True
        assert Leg("chr1", 0, Haplotypes.paternal).is_conflict() is False

    def test_ref_name_haplotype(self):
        leg = Leg("chr1", 100, Haplotypes.maternal)
        assert leg.ref_name_haplotype() == ("chr1", Haplotypes.maternal)

    def test_same_chr_with(self):
        a = Leg("chr1", 100, Haplotypes.paternal)
        b = Leg("chr1", 200, Haplotypes.maternal)
        c = Leg("chr2", 100, Haplotypes.paternal)
        assert a.same_chr_with(b) is True
        assert a.same_chr_with(c) is False

    def test_separation_with(self):
        a = Leg("chr1", 1000, Haplotypes.paternal)
        b = Leg("chr1", 3000, Haplotypes.paternal)
        assert a.separation_with(b) == 2000
        assert b.separation_with(a) == 2000

    def test_signed_separation_with(self):
        a = Leg("chr1", 1000, Haplotypes.paternal)
        b = Leg("chr1", 3000, Haplotypes.paternal)
        assert a.signed_separation_with(b) == -2000
        assert b.signed_separation_with(a) == 2000

    def test_merge_with(self):
        a = Leg("chr1", 1000, Haplotypes.unknown)
        b = Leg("chr1", 3000, Haplotypes.paternal)
        a.merge_with(b)
        assert a.get_ref_locus() == 2000
        assert a.get_haplotype() == Haplotypes.paternal

    def test_compatible_haps_phased(self):
        leg = Leg("chr1", 100, Haplotypes.paternal)
        assert leg.compatible_haps() == [Haplotypes.paternal]

    def test_compatible_haps_unphased(self):
        leg = Leg("chr1", 100, Haplotypes.unknown)
        assert set(leg.compatible_haps()) == {Haplotypes.paternal, Haplotypes.maternal}

    def test_compatible_legs_female(self):
        leg = Leg("chr1", 100, Haplotypes.unknown)
        legs = leg.compatible_legs_female()
        assert len(legs) == 2

    def test_in_reg(self):
        reg = Reg("chr1")
        reg.add_start(1000)
        reg.add_end(5000)
        assert Leg("chr1", 2000, Haplotypes.paternal).in_reg(reg) is True
        assert Leg("chr1", 500, Haplotypes.paternal).in_reg(reg) is False
        assert Leg("chr2", 2000, Haplotypes.paternal).in_reg(reg) is False

    def test_in_regs(self):
        reg1 = Reg("chr1")
        reg1.add_start(0)
        reg1.add_end(1000)
        reg2 = Reg("chr2")
        reg2.add_start(0)
        reg2.add_end(1000)
        assert Leg("chr1", 500, Haplotypes.paternal).in_regs([reg1, reg2]) is True
        assert Leg("chr3", 500, Haplotypes.paternal).in_regs([reg1, reg2]) is False

    def test_satisfy_regs(self):
        inc = Reg("chr1")
        exc = Reg("chr1")
        exc.add_start(500)
        exc.add_end(600)
        assert Leg("chr1", 100, Haplotypes.paternal).satisfy_regs([inc], [exc]) is True
        assert Leg("chr1", 550, Haplotypes.paternal).satisfy_regs([inc], [exc]) is False


# ---------------------------------------------------------------------------
# Con
# ---------------------------------------------------------------------------

class TestCon:
    def test_con_creation(self):
        leg1 = Leg("chr1", 1000000, Haplotypes.paternal)
        leg2 = Leg("chr1", 2000000, Haplotypes.maternal)
        con = Con(leg1, leg2)
        assert con.leg_1() == leg1
        assert con.leg_2() == leg2

    def test_con_ordering(self):
        con = Con(Leg("chr1", 2000000, Haplotypes.maternal), Leg("chr1", 1000000, Haplotypes.paternal))
        assert con.leg_1() < con.leg_2()

    def test_string_to_con(self):
        con = string_to_con("chr1,1000000,0\tchr1,2000000,1")
        assert con.leg_1().get_ref_locus() == 1000000
        assert con.leg_2().get_haplotype() == Haplotypes.maternal

    def test_con_to_string(self):
        con = Con(Leg("chr1", 1000000, Haplotypes.paternal), Leg("chr1", 2000000, Haplotypes.maternal))
        assert con.to_string() == "chr1,1000000,0\tchr1,2000000,1"

    def test_con_roundtrip(self):
        s = "chr1,1000000,0\tchr1,2000000,1"
        assert string_to_con(s).to_string() == s

    def test_is_intra_chr(self):
        intra = string_to_con("chr1,100,0\tchr1,200,1")
        inter = string_to_con("chr1,100,0\tchr2,200,1")
        assert intra.is_intra_chr() is True
        assert inter.is_intra_chr() is False

    def test_is_inter_hom(self):
        inter_hom = string_to_con("chr1,100,0\tchr1,200,1")
        intra_hom = string_to_con("chr1,100,0\tchr1,200,0")
        assert inter_hom.is_inter_hom() is True
        assert intra_hom.is_inter_hom() is False

    def test_separation(self):
        con = string_to_con("chr1,1000,0\tchr1,3000,1")
        assert con.separation() == 2000

    def test_num_phased_legs(self):
        assert string_to_con("chr1,1,0\tchr1,2,1").num_phased_legs() == 2
        assert string_to_con("chr1,1,.\tchr1,2,1").num_phased_legs() == 1
        assert string_to_con("chr1,1,.\tchr1,2,.").num_phased_legs() == 0

    def test_num_conflict_legs(self):
        leg_c = Leg("chr1", 100, Haplotypes.conflict)
        leg_p = Leg("chr1", 200, Haplotypes.paternal)
        assert Con(leg_c, leg_p).num_conflict_legs() == 1

    def test_ref_name_tuple(self):
        con = string_to_con("chr1,100,0\tchr2,200,1")
        assert con.ref_name_tuple() == ("chr1", "chr2")

    def test_hap_tuple(self):
        con = string_to_con("chr1,100,0\tchr1,200,1")
        assert con.hap_tuple() == (Haplotypes.paternal, Haplotypes.maternal)

    def test_compatible_hap_tuples(self):
        # Fully phased: 1 compatible
        con = string_to_con("chr1,100,0\tchr1,200,1")
        assert con.compatible_hap_tuples() == [(0, 1)]
        # Fully unphased: 4 compatible
        con2 = string_to_con("chr1,100,.\tchr1,200,.")
        assert len(con2.compatible_hap_tuples()) == 4

    def test_set_hap_tuple(self):
        con = string_to_con("chr1,100,.\tchr1,200,.")
        con.set_hap_tuple((Haplotypes.paternal, Haplotypes.maternal))
        assert con.leg_1().get_haplotype() == Haplotypes.paternal
        assert con.leg_2().get_haplotype() == Haplotypes.maternal

    def test_distance_functions(self):
        c1 = string_to_con("chr1,1000,0\tchr1,5000,0")
        c2 = string_to_con("chr1,2000,0\tchr1,4000,0")
        assert c1.distance_leg_1_with(c2) == 1000
        assert c1.distance_leg_2_with(c2) == 1000
        assert c1.distance_inf_with(c2) == 1000

    def test_to_rel_locus_around(self):
        c1 = string_to_con("chr1,3000,0\tchr1,5000,0")
        c2 = string_to_con("chr1,1000,0\tchr1,2000,0")
        assert c1.to_rel_locus_around(c2) == (2000, 3000)

    def test_to_string_around(self):
        c1 = string_to_con("chr1,3000,0\tchr1,5000,0")
        c2 = string_to_con("chr1,1000,0\tchr1,2000,0")
        assert c1.to_string_around(c2) == "2000\t3000"

    def test_merge_with(self):
        c1 = string_to_con("chr1,1000,.\tchr1,3000,.")
        c2 = string_to_con("chr1,1000,0\tchr1,3000,1")
        c1.merge_with(c2)
        assert c1.leg_1().get_haplotype() == Haplotypes.paternal
        assert c1.leg_2().get_haplotype() == Haplotypes.maternal

    def test_equality_and_hash(self):
        c1 = string_to_con("chr1,100,0\tchr1,200,1")
        c2 = string_to_con("chr1,100,0\tchr1,200,1")
        assert c1 == c2
        assert hash(c1) == hash(c2)

    def test_lt(self):
        c1 = string_to_con("chr1,100,0\tchr1,200,0")
        c2 = string_to_con("chr1,100,0\tchr1,300,0")
        assert c1 < c2


# ---------------------------------------------------------------------------
# ConData
# ---------------------------------------------------------------------------

class TestConData:
    def test_condata_creation(self):
        assert ConData().num_cons() == 0

    def test_condata_add_con(self):
        cd = ConData()
        cd.add_con(string_to_con("chr1,1000,0\tchr1,2000,1"))
        assert cd.num_cons() == 1

    def test_file_to_con_data(self):
        text = "chr1,1000,0\tchr1,2000,1\nchr1,3000,.\tchr1,4000,.\nchr2,500,0\tchr2,1500,1"
        cd = file_to_con_data(StringIO(text))
        assert cd.num_cons() == 3

    def test_num_phased_legs(self):
        cd = ConData()
        cd.add_con(string_to_con("chr1,1,0\tchr1,2,1"))
        cd.add_con(string_to_con("chr1,3,.\tchr1,4,."))
        assert cd.num_phased_legs() == 2

    def test_num_intra_chr(self):
        cd = ConData()
        cd.add_con(string_to_con("chr1,1,0\tchr1,2,1"))
        cd.add_con(string_to_con("chr1,3,0\tchr2,4,0"))
        assert cd.num_intra_chr() == 1

    def test_num_phased_cons(self):
        cd = ConData()
        cd.add_con(string_to_con("chr1,1,0\tchr1,2,1"))
        cd.add_con(string_to_con("chr1,3,.\tchr1,4,1"))
        assert cd.num_phased_cons() == 1

    def test_num_conflict_legs(self):
        cd = ConData()
        con = Con(Leg("chr1", 100, Haplotypes.conflict), Leg("chr1", 200, Haplotypes.paternal))
        cd.add_con(con)
        assert cd.num_conflict_legs() == 1

    def test_sort_cons(self):
        cd = ConData()
        cd.add_con(string_to_con("chr1,3000,0\tchr1,4000,0"))
        cd.add_con(string_to_con("chr1,1000,0\tchr1,2000,0"))
        cd.sort_cons()
        cons = list(cd.get_cons())
        assert cons[0].leg_1().get_ref_locus() == 1000

    def test_to_string(self):
        cd = ConData()
        cd.add_con(string_to_con("chr1,100,0\tchr1,200,1"))
        assert cd.to_string() == "chr1,100,0\tchr1,200,1"

    def test_to_string_roundtrip(self):
        text = "chr1,100,0\tchr1,200,1\nchr1,300,0\tchr1,400,0"
        cd = file_to_con_data(StringIO(text))
        cd.sort_cons()
        assert cd.to_string() == text

    def test_merge_with(self):
        cd1 = ConData()
        cd1.add_con(string_to_con("chr1,100,0\tchr1,200,1"))
        cd2 = ConData()
        cd2.add_con(string_to_con("chr2,100,0\tchr2,200,1"))
        cd1.merge_with(cd2)
        assert cd1.num_cons() == 2

    def test_get_cons(self):
        cd = ConData()
        cd.add_con(string_to_con("chr1,100,0\tchr1,200,1"))
        cd.add_con(string_to_con("chr2,100,0\tchr2,200,1"))
        cons = list(cd.get_cons())
        assert len(cons) == 2


# ---------------------------------------------------------------------------
# Voting / imputation helpers
# ---------------------------------------------------------------------------

class TestVoting:
    def test_vote_from_hap_tuples_one_match(self):
        self_haps = [(0, 0), (0, 1)]
        other_haps = [(0, 1), (1, 1)]
        assert vote_from_hap_tuples(self_haps, other_haps) == (0, 1)

    def test_vote_from_hap_tuples_no_match(self):
        self_haps = [(0, 0)]
        other_haps = [(1, 1)]
        assert vote_from_hap_tuples(self_haps, other_haps) is None

    def test_vote_from_hap_tuples_multiple_matches(self):
        self_haps = [(0, 0), (0, 1)]
        other_haps = [(0, 0), (0, 1)]
        assert vote_from_hap_tuples(self_haps, other_haps) is None

    def test_winning_vote_meets_criteria(self):
        votes = [(0, 1), (0, 1), (0, 1), (1, 0)]
        assert winning_vote(votes, 2, 0.5) == (0, 1)

    def test_winning_vote_below_min_votes(self):
        votes = [(0, 1)]
        assert winning_vote(votes, 2, 0.5) is None

    def test_winning_vote_below_fraction(self):
        votes = [(0, 1), (1, 0), (0, 0)]
        assert winning_vote(votes, 1, 0.5) is None

    def test_winning_vote_empty(self):
        assert winning_vote([], 1, 0.5) is None


# ---------------------------------------------------------------------------
# Histogram helpers
# ---------------------------------------------------------------------------

class TestHistogramHelpers:
    def test_counts_to_hist_num_with_zero(self):
        hist = counts_to_hist_num_with_zero([0, 1, 1, 2, 3])
        assert hist == [1, 2, 1, 1]

    def test_hist_num_to_string(self):
        s = hist_num_to_string([5, 3, 2])
        assert "1\t5" in s
        assert ">=3\t2" in s

    def test_hist_num_to_string_with_zero(self):
        s = hist_num_to_string_with_zero([5, 3, 2])
        assert "0\t5" in s
        assert ">=2\t2" in s


# ---------------------------------------------------------------------------
# Reg
# ---------------------------------------------------------------------------

class TestReg:
    def test_reg_creation(self):
        reg = Reg("chr1")
        reg.add_haplotype(Haplotypes.paternal)
        reg.add_start(1000000)
        reg.add_end(2000000)
        assert reg.ref_name == "chr1"
        assert reg.haplotype == Haplotypes.paternal
        assert reg.start == 1000000
        assert reg.end == 2000000

    def test_reg_defaults(self):
        reg = Reg("chr1")
        assert reg.has_haplotype is False
        assert reg.has_start is False
        assert reg.has_end is False

    def test_remove_haplotype(self):
        reg = Reg("chr1")
        reg.add_haplotype(Haplotypes.paternal)
        assert reg.has_haplotype is True
        reg.remove_haplotype()
        assert reg.has_haplotype is False

    def test_string_to_reg(self):
        reg = string_to_reg("chr1\t0\t1000000\t2000000")
        assert reg.ref_name == "chr1"
        assert reg.haplotype == Haplotypes.paternal
        assert reg.start == 1000000
        assert reg.end == 2000000

    def test_string_to_reg_dots(self):
        reg = string_to_reg("chr1\t.\t.\t.")
        assert reg.has_haplotype is False
        assert reg.has_start is False
        assert reg.has_end is False

    def test_file_to_reg_list(self):
        text = "chr1\t.\t.\t.\nchr2\t1\t232800000\t."
        assert len(file_to_reg_list(StringIO(text))) == 2

    def test_to_string(self):
        reg = string_to_reg("chr1\t0\t1000\t2000")
        assert reg.to_string() == "chr1\t0\t1000\t2000"

    def test_to_string_roundtrip(self):
        s = "chr1\t.\t.\t."
        assert string_to_reg(s).to_string() == s

    def test_to_name_string(self):
        reg = string_to_reg("chr1\t0\t1000\t2000")
        assert reg.to_name_string() == "chr1(pat)_1000-2000"

    def test_to_name_string_no_hap(self):
        reg = string_to_reg("chr1\t.\t.\t.")
        assert reg.to_name_string() == "chr1"

    def test_get_phased_with_haplotype(self):
        reg = string_to_reg("chr1\t0\t100\t200")
        phased = list(reg.get_phased())
        assert len(phased) == 1
        assert phased[0].haplotype == Haplotypes.paternal

    def test_get_phased_without_haplotype(self):
        reg = string_to_reg("chr1\t.\t100\t200")
        phased = list(reg.get_phased())
        assert len(phased) == 2

    def test_get_unphased(self):
        reg = string_to_reg("chr1\t0\t100\t200")
        u = reg.get_unphased()
        assert u.has_haplotype is False

    def test_get_phased_regs(self):
        regs = file_to_reg_list(StringIO("chr1\t.\t.\t.\nchr2\t0\t100\t200"))
        phased = list(get_phased_regs(regs))
        # chr1 unphased → 2 phased copies, chr2 already phased → 1
        assert len(phased) == 3


# ---------------------------------------------------------------------------
# LegData / LegList
# ---------------------------------------------------------------------------

class TestLegData:
    def test_legdata_creation(self):
        assert LegData().num_legs() == 0

    def test_legdata_add_leg(self):
        ld = LegData()
        ld.add_leg(Leg("chr1", 1000, Haplotypes.paternal))
        assert ld.num_legs() == 1

    def test_legdata_to_string(self):
        ld = LegData()
        ld.add_leg(Leg("chr1", 1000, Haplotypes.paternal))
        s = ld.to_string()
        assert "chr1,1000,0" in s

    def test_legdata_add_con(self):
        ld = LegData()
        con = string_to_con("chr1,100,0\tchr1,200,1")
        ld.add_con(con)
        assert ld.num_legs() == 2


class TestLegList:
    def test_sort_legs(self):
        ll = LegList()
        ll.add_leg(Leg("chr1", 300, Haplotypes.paternal))
        ll.add_leg(Leg("chr1", 100, Haplotypes.paternal))
        ll.sort_legs()
        assert ll.legs[0].get_ref_locus() == 100

    def test_num_legs(self):
        ll = LegList()
        ll.add_leg(Leg("chr1", 100, Haplotypes.paternal))
        assert ll.num_legs() == 1


# ---------------------------------------------------------------------------
# DupLeg / DupCon / DupConData
# ---------------------------------------------------------------------------

class TestDuplicateClasses:
    def test_dup_leg(self):
        leg = Leg("chr1", 100, Haplotypes.paternal)
        dl = DupLeg(leg)
        assert dl.num_dups() == 1
        assert dl.dups[Haplotypes.paternal] == 1

    def test_dup_leg_merge(self):
        dl1 = DupLeg(Leg("chr1", 100, Haplotypes.paternal))
        dl2 = DupLeg(Leg("chr1", 200, Haplotypes.maternal))
        dl1.merge_with(dl2)
        assert dl1.num_dups() == 2
        assert dl1.dups[Haplotypes.paternal] == 1
        assert dl1.dups[Haplotypes.maternal] == 1

    def test_dup_con(self):
        con = string_to_con("chr1,100,0\tchr1,200,1")
        dc = DupCon(con)
        assert dc.num_dups() == 1

    def test_dup_con_data_stats(self):
        cd = ConData()
        cd.add_con(string_to_con("chr1,100,0\tchr1,200,0"))
        cd.add_con(string_to_con("chr1,300,0\tchr1,400,0"))
        cd.sort_cons()
        dcd = DupConData(cd)
        hist = dcd.dup_stats(3)
        assert sum(hist) == 2
        assert hist[0] == 2  # all have 1 dup


# ---------------------------------------------------------------------------
# Par / ParData
# ---------------------------------------------------------------------------

class TestPar:
    def test_par_contain_leg(self):
        par = Par("X", 60000, 2700000, "Y", 10000)
        assert par.contain_leg(Leg("X", 100000, Haplotypes.paternal)) is True
        assert par.contain_leg(Leg("X", 3000000, Haplotypes.paternal)) is False

    def test_par_data_set_non_par(self):
        pd = ParData("X", "Y")
        pd.add_par(Par("X", 60000, 2700000, "Y", 10000))
        leg = Leg("X", 5000000, Haplotypes.unknown)
        pd.set_non_par_leg_haplotype_male(leg)
        assert leg.get_haplotype() == Haplotypes.maternal

        leg_y = Leg("Y", 5000000, Haplotypes.unknown)
        pd.set_non_par_leg_haplotype_male(leg_y)
        assert leg_y.get_haplotype() == Haplotypes.paternal


# ---------------------------------------------------------------------------
# G3dParticle
# ---------------------------------------------------------------------------

class TestG3dParticle:
    def test_g3d_particle_creation(self):
        p = G3dParticle("chr1(pat)", 1000000, [10.5, 20.3, -5.7])
        assert p.get_hom_name() == "chr1(pat)"
        assert p.get_ref_locus() == 1000000
        assert p.get_x() == 10.5
        assert p.get_y() == 20.3
        assert p.get_z() == -5.7

    def test_string_to_g3d_particle(self):
        p = string_to_g3d_particle("chr1(pat)\t1000000\t10.5\t20.3\t-5.7")
        assert p.get_hom_name() == "chr1(pat)"

    def test_g3d_particle_to_string_roundtrip(self):
        s = "chr1(pat)\t1000000\t10.5\t20.3\t-5.7"
        assert string_to_g3d_particle(s).to_string() == s

    def test_g3d_particle_ordering(self):
        p1 = G3dParticle("chr1(pat)", 1000000, [0, 0, 0])
        p2 = G3dParticle("chr1(pat)", 2000000, [0, 0, 0])
        assert p1 < p2

    def test_get_ref_name(self):
        p = G3dParticle("chr1(pat)", 100, [0, 0, 0])
        assert p.get_ref_name() == "chr1"

    def test_get_haplotype(self):
        p = G3dParticle("chr1(mat)", 100, [0, 0, 0])
        assert p.get_haplotype() == Haplotypes.maternal

    def test_to_leg(self):
        p = G3dParticle("chr1(pat)", 100, [0, 0, 0])
        leg = p.to_leg()
        assert leg.get_ref_name() == "chr1"
        assert leg.get_ref_locus() == 100
        assert leg.get_haplotype() == Haplotypes.paternal

    def test_get_set_position(self):
        p = G3dParticle("chr1(pat)", 100, [1.0, 2.0, 3.0])
        assert p.get_position() == [1.0, 2.0, 3.0]
        p.set_position([4.0, 5.0, 6.0])
        assert p.get_x() == 4.0

    def test_ref_name_haplotype(self):
        p = G3dParticle("22(mat)", 100, [0, 0, 0])
        assert p.ref_name_haplotype() == ("22", Haplotypes.maternal)

    def test_in_reg(self):
        reg = Reg("22")
        reg.add_start(100)
        reg.add_end(500)
        assert G3dParticle("22(pat)", 200, [0, 0, 0]).in_reg(reg) is True
        assert G3dParticle("22(pat)", 50, [0, 0, 0]).in_reg(reg) is False
        assert G3dParticle("1(pat)", 200, [0, 0, 0]).in_reg(reg) is False

    def test_satisfy_regs(self):
        inc = Reg("22")
        exc = Reg("22")
        exc.add_start(400)
        exc.add_end(600)
        assert G3dParticle("22(pat)", 200, [0, 0, 0]).satisfy_regs([inc], [exc]) is True
        assert G3dParticle("22(pat)", 500, [0, 0, 0]).satisfy_regs([inc], [exc]) is False


# ---------------------------------------------------------------------------
# G3dData
# ---------------------------------------------------------------------------

class TestG3dData:
    def _make_g3d_data(self):
        text = (
            "chr1(pat)\t1000000\t10.0\t20.0\t-5.0\n"
            "chr1(pat)\t2000000\t11.0\t21.0\t-6.0\n"
            "chr1(pat)\t3000000\t12.0\t22.0\t-7.0\n"
            "chr1(mat)\t1000000\t-10.0\t-20.0\t5.0\n"
            "chr1(mat)\t2000000\t-11.0\t-21.0\t6.0"
        )
        return file_to_g3d_data(StringIO(text))

    def test_g3d_data_creation(self):
        assert G3dData().num_g3d_particles() == 0

    def test_g3d_data_add_particle(self):
        gd = G3dData()
        gd.add_g3d_particle(G3dParticle("chr1(pat)", 100, [0, 0, 0]))
        assert gd.num_g3d_particles() == 1

    def test_file_to_g3d_data(self):
        gd = self._make_g3d_data()
        assert gd.num_g3d_particles() == 5

    def test_get_hom_names(self):
        gd = self._make_g3d_data()
        names = list(gd.get_hom_names())
        assert "chr1(mat)" in names
        assert "chr1(pat)" in names

    def test_get_g3d_list_from_hom_name(self):
        gd = self._make_g3d_data()
        assert gd.get_g3d_list_from_hom_name("chr1(pat)").num_g3d_particles() == 3
        assert gd.get_g3d_list_from_hom_name("nonexistent") is None

    def test_get_g3d_particle_from_hom_name_ref_locus(self):
        gd = self._make_g3d_data()
        p = gd.get_g3d_particle_from_hom_name_ref_locus("chr1(pat)", 1000000)
        assert p is not None
        assert p.get_x() == 10.0
        assert gd.get_g3d_particle_from_hom_name_ref_locus("chr1(pat)", 999) is None

    def test_sort_g3d_particles(self):
        gd = G3dData()
        gd.add_g3d_particle(G3dParticle("chr1(pat)", 300, [0, 0, 0]))
        gd.add_g3d_particle(G3dParticle("chr1(pat)", 100, [0, 0, 0]))
        gd.sort_g3d_particles()
        particles = list(gd.get_g3d_particles())
        assert particles[0].get_ref_locus() == 100

    def test_resolution(self):
        gd = self._make_g3d_data()
        gd.sort_g3d_particles()
        assert gd.resolution() == 1000000

    def test_resolution_empty(self):
        assert G3dData().resolution() is None

    def test_to_string(self):
        gd = self._make_g3d_data()
        s = gd.to_string()
        assert "chr1(pat)" in s
        assert "chr1(mat)" in s

    def test_to_string_roundtrip(self):
        text = "chr1(mat)\t100\t-1.0\t-2.0\t3.0\nchr1(pat)\t100\t1.0\t2.0\t-3.0"
        gd = file_to_g3d_data(StringIO(text))
        assert gd.to_string() == text

    def test_to_np_arrays(self):
        gd = self._make_g3d_data()
        hom_names, loci, positions = gd.to_np_arrays()
        assert len(hom_names) == 5
        assert loci.shape == (5, 1)
        assert positions.shape == (5, 3)

    def test_get_g3d_particles(self):
        gd = self._make_g3d_data()
        assert len(list(gd.get_g3d_particles())) == 5

    def test_get_adjacent_g3d_particle_tuples(self):
        gd = self._make_g3d_data()
        gd.sort_g3d_particles()
        tuples = list(gd.get_adjacent_g3d_particle_tuples())
        # chr1(mat): 2 particles → 1 tuple, chr1(pat): 3 → 2 tuples
        assert len(tuples) == 3

    def test_get_adjacent_g3d_particle_tuples_given_separation(self):
        gd = self._make_g3d_data()
        gd.sort_g3d_particles()
        tuples = list(gd.get_adjacent_g3d_particle_tuples_given_separation(1000000))
        assert len(tuples) == 3

    def test_prepare_interpolate_and_interpolate_leg(self):
        gd = self._make_g3d_data()
        gd.sort_g3d_particles()
        gd.prepare_interpolate()
        # Interpolate at a known point
        is_out, pos = gd.interpolate_leg(Leg("chr1", 1000000, Haplotypes.paternal))
        assert is_out is False
        assert abs(pos[0] - 10.0) < 0.01
        # Interpolate midpoint
        is_out, pos = gd.interpolate_leg(Leg("chr1", 1500000, Haplotypes.paternal))
        assert is_out is False
        assert abs(pos[0] - 10.5) < 0.01
        # Unknown homolog
        is_out, pos = gd.interpolate_leg(Leg("chr99", 100, Haplotypes.paternal))
        assert pos is None

    def test_apply_regs(self):
        gd = self._make_g3d_data()
        inc = Reg("chr1")
        inc.add_haplotype(Haplotypes.paternal)
        gd.apply_regs([inc], [])
        assert gd.num_g3d_particles() == 3


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
