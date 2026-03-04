import pytest
from io import StringIO

from dip_c.classes import (
    Haplotypes, is_known_haplotype, haplotype_to_string, string_to_haplotype,
    known_haplotypes, Leg, string_to_leg, Con, string_to_con,
    ConData, file_to_con_data, LegData,
    G3dParticle, string_to_g3d_particle,
    G3dData, file_to_g3d_data, Reg, string_to_reg, file_to_reg_list,
    ref_name_haplotype_to_hom_name, hom_name_to_ref_name_haplotype
)


class TestHaplotypes:
    """Test haplotype enumeration and utility functions."""
    
    def test_haplotype_values(self):
        """Test that haplotype enum values are correct."""
        assert Haplotypes.paternal == 0
        assert Haplotypes.maternal == 1
        assert Haplotypes.unknown == -1
        assert Haplotypes.conflict == -2
        assert Haplotypes.minus_infinity == -3
        assert Haplotypes.infinity == 2
    
    def test_is_known_haplotype(self):
        """Test haplotype identification."""
        assert is_known_haplotype(Haplotypes.paternal) == True
        assert is_known_haplotype(Haplotypes.maternal) == True
        assert is_known_haplotype(0) == True
        assert is_known_haplotype(1) == True
        assert is_known_haplotype(Haplotypes.unknown) == False
        assert is_known_haplotype(Haplotypes.conflict) == False
        assert is_known_haplotype(-1) == False
    
    def test_haplotype_to_string(self):
        """Test haplotype to string conversion."""
        assert haplotype_to_string(Haplotypes.paternal) == "0"
        assert haplotype_to_string(Haplotypes.maternal) == "1"
        assert haplotype_to_string(Haplotypes.unknown) == "."
        assert haplotype_to_string(Haplotypes.conflict) == "."
        assert haplotype_to_string(-1) == "."
    
    def test_string_to_haplotype(self):
        """Test string to haplotype conversion."""
        assert string_to_haplotype("0") == Haplotypes.paternal
        assert string_to_haplotype("1") == Haplotypes.maternal
        assert string_to_haplotype(".") == Haplotypes.unknown
    
    def test_known_haplotypes_generator(self):
        """Test that known_haplotypes yields correct values."""
        haps = list(known_haplotypes())
        assert len(haps) == 2
        assert Haplotypes.paternal in haps
        assert Haplotypes.maternal in haps


class TestLeg:
    """Test Leg class and related functions."""
    
    def test_leg_creation(self):
        """Test creating a Leg object."""
        leg = Leg("chr1", 1000000, Haplotypes.paternal)
        assert leg.get_ref_name() == "chr1"
        assert leg.get_ref_locus() == 1000000
        assert leg.get_haplotype() == Haplotypes.paternal
    
    def test_leg_equality(self):
        """Test Leg equality comparison."""
        leg1 = Leg("chr1", 1000000, Haplotypes.paternal)
        leg2 = Leg("chr1", 1000000, Haplotypes.paternal)
        leg3 = Leg("chr1", 1000000, Haplotypes.maternal)
        leg4 = Leg("chr1", 2000000, Haplotypes.paternal)
        
        assert leg1 == leg2
        assert leg1 != leg3  # different haplotype
        assert leg1 != leg4  # different locus
    
    def test_leg_ordering(self):
        """Test Leg comparison operators."""
        leg1 = Leg("chr1", 1000000, Haplotypes.paternal)
        leg2 = Leg("chr1", 2000000, Haplotypes.paternal)
        leg3 = Leg("chr2", 500000, Haplotypes.paternal)
        
        assert leg1 < leg2
        assert leg2 > leg1
        assert leg1 < leg3  # chr1 < chr2
    
    def test_string_to_leg(self):
        """Test parsing leg from string."""
        leg = string_to_leg("chr1,1000000,0")
        assert leg.get_ref_name() == "chr1"
        assert leg.get_ref_locus() == 1000000
        assert leg.get_haplotype() == Haplotypes.paternal
        
        leg_unknown = string_to_leg("chr2,2000000,.")
        assert leg_unknown.get_ref_name() == "chr2"
        assert leg_unknown.get_ref_locus() == 2000000
        assert leg_unknown.get_haplotype() == Haplotypes.unknown
    
    def test_leg_to_string(self):
        """Test converting leg to string."""
        leg = Leg("chr1", 1000000, Haplotypes.paternal)
        assert leg.to_string() == "chr1,1000000,0"
        
        leg_unknown = Leg("chr2", 2000000, Haplotypes.unknown)
        assert leg_unknown.to_string() == "chr2,2000000,."
    
    def test_leg_roundtrip(self):
        """Test that string_to_leg and to_string are inverses."""
        original = "chr1,1000000,0"
        leg = string_to_leg(original)
        back_to_string = leg.to_string()
        assert original == back_to_string


class TestCon:
    """Test Con (contact) class and related functions."""
    
    def test_con_creation(self):
        """Test creating a Con object."""
        leg1 = Leg("chr1", 1000000, Haplotypes.paternal)
        leg2 = Leg("chr1", 2000000, Haplotypes.maternal)
        con = Con(leg1, leg2)
        
        assert con.leg_1() == leg1
        assert con.leg_2() == leg2
    
    def test_con_ordering(self):
        """Test that legs are automatically ordered."""
        leg1 = Leg("chr1", 2000000, Haplotypes.maternal)
        leg2 = Leg("chr1", 1000000, Haplotypes.paternal)
        con = Con(leg1, leg2)
        
        # Legs should be reordered so leg_1() < leg_2()
        assert con.leg_1() < con.leg_2()
    
    def test_string_to_con(self):
        """Test parsing contact from string."""
        con = string_to_con("chr1,1000000,0\tchr1,2000000,1")
        
        assert con.leg_1().get_ref_name() == "chr1"
        assert con.leg_1().get_ref_locus() == 1000000
        assert con.leg_1().get_haplotype() == Haplotypes.paternal
        assert con.leg_2().get_ref_name() == "chr1"
        assert con.leg_2().get_ref_locus() == 2000000
        assert con.leg_2().get_haplotype() == Haplotypes.maternal
    
    def test_con_to_string(self):
        """Test converting contact to string."""
        leg1 = Leg("chr1", 1000000, Haplotypes.paternal)
        leg2 = Leg("chr1", 2000000, Haplotypes.maternal)
        con = Con(leg1, leg2)
        
        assert con.to_string() == "chr1,1000000,0\tchr1,2000000,1"
    
    def test_con_roundtrip(self):
        """Test that string_to_con and to_string are inverses."""
        original = "chr1,1000000,0\tchr1,2000000,1"
        con = string_to_con(original)
        back_to_string = con.to_string()
        assert original == back_to_string


class TestConData:
    """Test ConData class for managing multiple contacts."""
    
    def test_condata_creation(self):
        """Test creating an empty ConData object."""
        con_data = ConData()
        assert con_data.num_cons() == 0
    
    def test_condata_add_con(self):
        """Test adding contacts to ConData."""
        con_data = ConData()
        leg1 = Leg("chr1", 1000000, Haplotypes.paternal)
        leg2 = Leg("chr1", 2000000, Haplotypes.maternal)
        con = Con(leg1, leg2)
        
        con_data.add_con(con)
        assert con_data.num_cons() == 1
    
    def test_file_to_con_data(self):
        """Test reading contacts from file-like object."""
        con_text = """chr1,1000000,0\tchr1,2000000,1
chr1,3000000,.\tchr1,4000000,.
chr2,500000,0\tchr2,1500000,1"""
        
        con_file = StringIO(con_text)
        con_data = file_to_con_data(con_file)
        
        assert con_data.num_cons() == 3
    
    def test_condata_phased_legs(self):
        """Test counting phased legs."""
        con_data = ConData()
        
        # Add contact with both legs phased
        con1 = string_to_con("chr1,1000000,0\tchr1,2000000,1")
        con_data.add_con(con1)
        
        # Add contact with one leg phased
        con2 = string_to_con("chr1,3000000,.\tchr1,4000000,1")
        con_data.add_con(con2)
        
        # Add contact with no legs phased
        con3 = string_to_con("chr1,5000000,.\tchr1,6000000,.")
        con_data.add_con(con3)
        
        # Total: 3 phased legs out of 6 total
        assert con_data.num_phased_legs() == 3


class TestG3dParticle:
    """Test G3dParticle class for 3D genome particles."""
    
    def test_g3d_particle_creation(self):
        """Test creating a G3dParticle."""
        particle = G3dParticle("chr1(pat)", 1000000, [10.5, 20.3, -5.7])
        
        assert particle.get_hom_name() == "chr1(pat)"
        assert particle.get_ref_locus() == 1000000
        assert particle.get_x() == 10.5
        assert particle.get_y() == 20.3
        assert particle.get_z() == -5.7
    
    def test_string_to_g3d_particle(self):
        """Test parsing G3dParticle from string."""
        particle = string_to_g3d_particle("chr1(pat)\t1000000\t10.5\t20.3\t-5.7")
        
        assert particle.get_hom_name() == "chr1(pat)"
        assert particle.get_ref_locus() == 1000000
        assert abs(particle.get_x() - 10.5) < 0.001
        assert abs(particle.get_y() - 20.3) < 0.001
        assert abs(particle.get_z() - (-5.7)) < 0.001
    
    def test_g3d_particle_to_string(self):
        """Test converting G3dParticle to string."""
        particle = G3dParticle("chr1(pat)", 1000000, [10.5, 20.3, -5.7])
        result = particle.to_string()
        
        assert "chr1(pat)" in result
        assert "1000000" in result
    
    def test_g3d_particle_ordering(self):
        """Test G3dParticle comparison."""
        p1 = G3dParticle("chr1(pat)", 1000000, [0, 0, 0])
        p2 = G3dParticle("chr1(pat)", 2000000, [0, 0, 0])
        p3 = G3dParticle("chr2(pat)", 1000000, [0, 0, 0])
        
        assert p1 < p2  # same chromosome, different locus
        assert p1 < p3  # different chromosome


class TestG3dData:
    """Test G3dData class for managing 3D genome structures."""
    
    def test_g3d_data_creation(self):
        """Test creating an empty G3dData object."""
        g3d_data = G3dData()
        assert g3d_data.num_g3d_particles() == 0
    
    def test_g3d_data_add_particle(self):
        """Test adding particles to G3dData."""
        g3d_data = G3dData()
        particle = G3dParticle("chr1(pat)", 1000000, [10.5, 20.3, -5.7])
        
        g3d_data.add_g3d_particle(particle)
        assert g3d_data.num_g3d_particles() == 1
    
    def test_file_to_g3d_data(self):
        """Test reading 3D structure from file-like object."""
        g3d_text = """chr1(pat)\t1000000\t10.5\t20.3\t-5.7
chr1(pat)\t2000000\t11.2\t21.1\t-6.3
chr1(mat)\t1000000\t-10.5\t-20.3\t5.7"""
        
        g3d_file = StringIO(g3d_text)
        g3d_data = file_to_g3d_data(g3d_file)
        
        assert g3d_data.num_g3d_particles() == 3


class TestReg:
    """Test Reg (region) class."""
    
    def test_reg_creation(self):
        """Test creating a Reg object."""
        reg = Reg("chr1", Haplotypes.paternal, 1000000, 2000000)
        
        assert reg.ref_name == "chr1"
        assert reg.haplotype == Haplotypes.paternal
        assert reg.start_locus() == 1000000
        assert reg.end_locus() == 2000000
    
    def test_string_to_reg(self):
        """Test parsing region from string."""
        reg = string_to_reg("chr1\t0\t1000000\t2000000")
        
        assert reg.ref_name() == "chr1"
        assert reg.haplotype() == Haplotypes.paternal
        assert reg.start_locus() == 1000000
        assert reg.end_locus() == 2000000
    
    def test_string_to_reg_unknown_haplotype(self):
        """Test parsing region with unknown haplotype."""
        reg = string_to_reg("chr1\t.\t1000000\t2000000")
        
        assert reg.haplotype() == Haplotypes.unknown
    
    def test_file_to_reg_list(self):
        """Test reading regions from file-like object."""
        reg_text = """chr1\t.\t.\t.
chr2\t1\t232800000\t.
chr2\t.\t238500000\t."""
        
        reg_file = StringIO(reg_text)
        reg_list = file_to_reg_list(reg_file)
        
        assert len(reg_list) == 3


class TestHaploidNames:
    """Test haploid name conversion functions."""
    
    def test_ref_name_haplotype_to_hom_name(self):
        """Test converting ref name + haplotype tuple to homolog name."""
        assert ref_name_haplotype_to_hom_name(("chr1", Haplotypes.paternal)) == "chr1(pat)"
        assert ref_name_haplotype_to_hom_name(("chr1", Haplotypes.maternal)) == "chr1(mat)"
    
    def test_hom_name_to_ref_name_haplotype(self):
        """Test parsing homolog name into ref name + haplotype tuple."""
        ref_name, haplotype = hom_name_to_ref_name_haplotype("chr1(pat)")
        assert ref_name == "chr1"
        assert haplotype == Haplotypes.paternal
        
        ref_name, haplotype = hom_name_to_ref_name_haplotype("chr1(mat)")
        assert ref_name == "chr1"
        assert haplotype == Haplotypes.maternal
    
    def test_hom_name_roundtrip(self):
        """Test that conversion functions are inverses."""
        ref_name_haplotype = ("chr1", Haplotypes.paternal)
        
        hom_name = ref_name_haplotype_to_hom_name(ref_name_haplotype)
        back_tuple = hom_name_to_ref_name_haplotype(hom_name)
        
        assert back_tuple == ref_name_haplotype


class TestReg:
    """Test Reg (region) class."""
    
    def test_reg_creation(self):
        """Test creating a Reg object."""
        reg = Reg("chr1")
        reg.add_haplotype(Haplotypes.paternal)
        reg.add_start(1000000)
        reg.add_end(2000000)
        
        assert reg.ref_name == "chr1"
        assert reg.haplotype == Haplotypes.paternal
        assert reg.start == 1000000
        assert reg.end == 2000000
    
    def test_string_to_reg(self):
        """Test parsing region from string."""
        reg = string_to_reg("chr1\t0\t1000000\t2000000")
        
        assert reg.ref_name == "chr1"
        assert reg.haplotype == Haplotypes.paternal
        assert reg.start == 1000000
        assert reg.end == 2000000
    
    def test_string_to_reg_unknown_haplotype(self):
        """Test parsing region with unknown haplotype."""
        reg = string_to_reg("chr1\t.\t1000000\t2000000")
        
        assert reg.haplotype == Haplotypes.unknown
    
    def test_file_to_reg_list(self):
        """Test reading regions from file-like object."""
        reg_text = """chr1\t.\t.\t.
chr2\t1\t232800000\t.
chr2\t.\t238500000\t."""
        
        reg_file = StringIO(reg_text)
        reg_list = file_to_reg_list(reg_file)
        
        assert len(reg_list) == 3


class TestLegData:
    """Test LegData class for managing leg collections."""
    
    def test_legdata_creation(self):
        """Test creating an empty LegData object."""
        leg_data = LegData()
        assert leg_data.num_legs() == 0
    
    def test_legdata_add_leg(self):
        """Test adding legs to LegData."""
        leg_data = LegData()
        leg = Leg("chr1", 1000000, Haplotypes.paternal)
        
        leg_data.add_leg(leg)
        assert leg_data.num_legs() == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
