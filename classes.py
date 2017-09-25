import bisect
import sys
from functools import total_ordering
import math
from collections import Counter
from scipy import interpolate
from scipy import spatial
import numpy as np
import copy


# enum for haplotypes
class Haplotypes:
    minus_infinity = -3
    conflict = -2
    unknown = -1
    paternal = 0
    maternal = 1
    infinity = 2
def is_known_haplotype(haplotype):
    return haplotype >= 0
# the operations below will set "." for all unphased haplotypes
def haplotype_to_string(haplotype):
    return str(haplotype) if is_known_haplotype(haplotype) else "."
def string_to_haplotype(haplotype_string):
    return Haplotypes.unknown if haplotype_string == "." else int(haplotype_string) 
def known_haplotypes():
    yield Haplotypes.paternal
    yield Haplotypes.maternal
def ref_name_haplotype_to_hom_name(ref_name, haplotype):
    return ref_name + "(" + ("pat" if haplotype == Haplotypes.paternal else "mat") + ")"
def hom_name_to_ref_name_haplotype(hom_name):
    ref_name, haplotype = hom_name.split("(")
    haplotype = haplotype.strip(")")
    haplotype = Haplotypes.paternal if haplotype == "pat" else Haplotypes.maternal
    return ref_name, haplotype

# rules for updating one haplotype with another (merging)
def update_haplotype(haplotype_1, haplotype_2):
    if haplotype_1 == haplotype_2:
        return haplotype_1
    elif haplotype_1 == Haplotypes.unknown:
        return haplotype_2
    elif haplotype_2 == Haplotypes.unknown:
        return haplotype_1
    return Haplotypes.conflict

# a read segment (alignment)
@total_ordering
class Seg:
    
    def __init__(self, is_read2, query_start, query_end, ref_name, ref_start, ref_end, is_reverse, haplotype = Haplotypes.unknown):
        self.is_read2 = is_read2
        self.query_start = query_start
        self.query_end = query_end
        self.ref_name = ref_name
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.is_reverse = is_reverse
        self.haplotype = haplotype # if not given, default is unknown
    
    # order: from left to right on the fragment
    # namely, read 1 before read 2
    # on read 1, order by left side
    # on read 2, order by minus right side
    def __lt__(self, other):
        return (self.is_read2, (-self.query_end if self.is_read2 else self.query_start)) < (other.is_read2, (-other.query_end if other.is_read2 else other.query_start))
    def __eq__(self, other):
        return (self.is_read2, (-self.query_end if self.is_read2 else self.query_start)) == (other.is_read2, (-other.query_end if other.is_read2 else other.query_start))

                
    def update_haplotype(self, is_read2, ref_name, ref_locus, haplotype):
        if self.is_read2 == is_read2 and self.ref_name == ref_name and ref_locus - 1 >= self.ref_start and ref_locus <= self.ref_end:
            self.haplotype = update_haplotype(self.haplotype, haplotype)
    
    def is_phased(self):
        return is_known_haplotype(self.haplotype)
    
    # output mapped loci of the left and right ends (ordered based on fragment)
    def ref_left(self):
        if self.is_read2 == self.is_reverse:
            return self.ref_start
        return self.ref_end

    def ref_right(self):
        if self.is_read2 == self.is_reverse:
            return self.ref_end
        return self.ref_start 
        
    def set_ref_left(self, ref_left):
        if self.is_read2 == self.is_reverse:
            self.ref_start = ref_left
        else:
            self.ref_end = ref_left

    def set_ref_right(self, ref_right):
        if self.is_read2 == self.is_reverse:
            self.ref_end = ref_right
        else:
            self.ref_start = ref_right
            
    def to_con_with(self, other):
        return Con(Leg(self.ref_name, self.ref_right(), self.haplotype), Leg(other.ref_name, other.ref_left(), other.haplotype))
    
    def to_string(self): # "m" is for mate
        return ",".join(["m" if self.is_read2 else ".", str(self.query_start), str(self.query_end), self.ref_name, str(self.ref_start), str(self.ref_end), "-" if self.is_reverse else "+", haplotype_to_string(self.haplotype)])

# create a segment from a string ("." will be set to unknown)
def string_to_seg(seg_string):
    is_read2, query_start, query_end, ref_name, ref_start, ref_end, is_reverse, haplotype = seg_string.split(",")
    is_read2 = True if is_read2 == "m" else False
    query_start = int(query_start)
    query_end = int(query_end)
    ref_start = int(ref_start)
    ref_end = int(ref_end)
    is_reverse = True if is_reverse == "-" else False
    haplotype = string_to_haplotype(haplotype)
    return Seg(is_read2, query_start, query_end, ref_name, ref_start, ref_end, is_reverse, haplotype)

# a read, containing all its segments
class Read:
    
    def __init__(self, name):
        self.name = name
        self.segs = []

    def add_seg(self, seg):
        self.segs.append(seg)
    
    def add_segs_from_read(self, read):
        self.segs += read.segs
    
    def num_segs(self):
        return len(self.segs)

    def num_phased_segs(self):
        num_phased_segs = 0
        for seg in self.segs:
            if seg.is_phased():
                num_phased_segs += 1
        return num_phased_segs
            
    def update_haplotype(self, is_read2, ref_name, ref_locus, haplotype):
        for seg in self.segs:
            seg.update_haplotype(is_read2, ref_name, ref_locus, haplotype)
            
    def sort_segs(self):
        self.segs.sort()
        
    def to_con_data(self, adjacent_only):
        self.sort_segs()
        con_data = ConData()
        for i in range(self.num_segs() - 1):
            for j in range(i + 1, self.num_segs()):
                if adjacent_only and j > i + 1:
                    break
                con_data.add_con(self.segs[i].to_con_with(self.segs[j]))
        return con_data
    
    def to_string(self):
        return self.name + "\t" + "\t".join([seg.to_string() for seg in self.segs])
        
# create a read from a string
def string_to_read(read_string):
    read_string_data = read_string.split("\t")
    read = Read(read_string_data[0])
    for seg_string in read_string_data[1:]:
        read.add_seg(string_to_seg(seg_string))
    return read

# a hash map of reads (a SEG file)
class SegData:
    
    def __init__(self):
        self.reads = {}

    def contains_read_name(self, name):
        return name in self.reads

    # add a read; merge if exists; ignore if empty
    def add_read(self, read):
        if read.num_segs() == 0:
            return # ignore empty reads
        if read.name not in self.reads: # add a new read
            self.reads[read.name] = read
        else: # add segments to an existing read
            self.reads[read.name].add_segs_from_read(read)
            
    # discard reads with a single segments
    def clean(self):
        for name in self.reads.keys():
            if self.reads[name].num_segs() < 2:
                del self.reads[name]
                
    # update haplotype for a specific read, if exists
    def update_haplotype(self, name, is_read2, ref_name, ref_locus, haplotype):
        if name in self.reads:
            self.reads[name].update_haplotype(is_read2, ref_name, ref_locus, haplotype)
            
    def num_reads(self):
        return len(self.reads)
        
    def num_segs(self):
        num_segs = 0
        for read in self.reads.values():
            num_segs += read.num_segs()
        return num_segs
    
    def num_phased_segs(self):
        num_phased_segs = 0
        for read in self.reads.values():
            num_phased_segs += read.num_phased_segs()
        return num_phased_segs
          
    def to_string(self): # no tailing new line
        return "\n".join(read.to_string() for read in self.reads.values())

# a leg
@total_ordering
class Leg:
    
    def __init__(self, ref_name, ref_locus, haplotype):
        self.ref_name = ref_name
        self.ref_locus = ref_locus
        self.haplotype = haplotype
    def __repr__(self):
        return self.to_string()
                
    def __lt__(self, other):
        return (self.ref_name, self.ref_locus, self.haplotype) < (other.ref_name, other.ref_locus, other.haplotype)
    def __eq__(self, other):
        return (self.ref_name, self.ref_locus, self.haplotype) == (other.ref_name, other.ref_locus, other.haplotype)
    def __hash__(self):
        return hash((self.ref_name, self.ref_locus, self.haplotype))
        
    def get_ref_name(self):
        return self.ref_name
    def get_ref_locus(self):
        return self.ref_locus
    def get_haplotype(self):
        return self.haplotype
    def set_haplotype(self, haplotype):
        self.haplotype = haplotype
    def in_x(self, par_data):
        return self.ref_name == par_data.get_x_name()
    def in_y(self, par_data):
        return self.ref_name == par_data.get_y_name()
    def ref_name_haplotype(self):
        return (self.ref_name, self.haplotype)
            
    # a list of compatible haplotypes for imputation
    def compatible_haps(self):
        return [haplotype for haplotype in ([self.haplotype] if self.is_phased() else known_haplotypes())]
            
    # a list of compatible phased legs   
    def compatible_legs_female(self):
        return [Leg(self.ref_name, self.ref_locus, haplotype) for haplotype in ([self.haplotype] if self.is_phased() else known_haplotypes())]
    def compatible_legs_male(self, par_data):
        return par_data.compatible_legs_male(self)
    def compatible_legs(self, is_male, par_data):
        return self.compatible_legs_male(par_data) if is_male else self.compatible_legs_female()        
           
    def is_phased(self):
        return is_known_haplotype(self.haplotype)
    def is_conflict(self):
        return self.haplotype == Haplotypes.conflict  
        
    def merge_with(self, other):
        self.ref_locus = (self.ref_locus + other.ref_locus)/2
        self.haplotype = update_haplotype(self.haplotype, other.haplotype)

    def same_chr_with(self, other):
        return self.ref_name == other.ref_name
        
    def separation_with(self, other):
        return abs(self.ref_locus - other.ref_locus)
    def signed_separation_with(self, other):
        return self.ref_locus - other.ref_locus 
         
    # check if a leg is in a region
    def in_reg(self, reg):
        if self.ref_name != reg.ref_name:
            return False
        if reg.has_haplotype and self.haplotype != reg.haplotype:
            return False
        if reg.has_start and self.ref_locus < reg.start:
            return False
        if reg.has_end and self.ref_locus > reg.end:
            return False
        return True
    # check if a leg is in a list of regions
    def in_regs(self, regs):
        for reg in regs:
            if self.in_reg(reg):
                return True
        return False
    # check if a leg is in a list of included regions, but not in a list of excluded regions
    def satisfy_regs(self, inc_regs, exc_regs):
        return self.in_regs(inc_regs) and not self.in_regs(exc_regs)
                        
    def to_string(self):
        return ",".join([self.ref_name, str(self.ref_locus), haplotype_to_string(self.haplotype)])

def string_to_leg(leg_string):
    ref_name, ref_locus, haplotype = leg_string.split(",")
    ref_locus = int(ref_locus)
    haplotype = string_to_haplotype(haplotype)
    return Leg(ref_name, ref_locus, haplotype)

# a list of legs, can be sorted to query number of legs in a region
class LegList:
    def __init__(self):
        self.legs = []
    def num_legs(self):
        return(len(self.legs))
    def sort_legs(self):
        self.legs.sort()
    def add_leg(self, leg):
        self.legs.append(leg)
    def add_con_data(self, con_data):
        for con in con_data.get_cons():
            self.add_con(con)
            
    # query a leg, regardless of haplotypes, assume sorted and that the list includes the leg itself
    def is_leg_promiscuous(self, leg, max_leg_distance, max_leg_count):
        index_to_check = bisect.bisect_left(self.legs, Leg(leg.get_ref_name(), leg.get_ref_locus() - max_leg_distance, Haplotypes.minus_infinity)) + max_leg_count
        if index_to_check >= len(self.legs):
            return False
        return self.legs[index_to_check] <= Leg(leg.get_ref_name(), leg.get_ref_locus() + max_leg_distance, Haplotypes.infinity)
        
    ## clean3
    # sort fully phased legs
    def sort_phased_legs(self):
        self.legs.sort(key = lambda x:(x.get_ref_name(), x.get_haplotype(), x.get_ref_locus()))
    # return num of legs near a G3dParticle, assume fully phased and sorted
    def num_legs_near_g3d_particle(self, g3d_particle, max_distance):
        # modified from bisect to support key
        lo = 0
        hi = len(self.legs)
        while lo < hi:
            mid = (lo + hi) // 2
            if (self.legs[mid].get_ref_name(), self.legs[mid].get_haplotype(), self.legs[mid].get_ref_locus()) < (g3d_particle.get_ref_name(), g3d_particle.get_haplotype(), g3d_particle.get_ref_locus() - max_distance):
                lo = mid + 1
            else:
                hi = mid
        left_index = lo
        hi = len(self.legs)
        while lo < hi:
            mid = (lo + hi) // 2
            if (g3d_particle.get_ref_name(), g3d_particle.get_haplotype(), g3d_particle.get_ref_locus() + max_distance) < (self.legs[mid].get_ref_name(), self.legs[mid].get_haplotype(), self.legs[mid].get_ref_locus()):
                hi = mid
            else:
                lo = mid + 1
        return lo - left_index
    def is_g3d_particle_leg_poor(self, g3d_particle, max_distance, min_leg_count):
        # modified from bisect to support key
        lo = 0
        hi = len(self.legs)
        while lo < hi:
            mid = (lo + hi) // 2
            if (self.legs[mid].get_ref_name(), self.legs[mid].get_haplotype(), self.legs[mid].get_ref_locus()) < (g3d_particle.get_ref_name(), g3d_particle.get_haplotype(), g3d_particle.get_ref_locus() - max_distance):
                lo = mid + 1
            else:
                hi = mid
        index_to_check = lo + min_leg_count - 1
        if index_to_check >= len(self.legs):
            return True
        return (g3d_particle.get_ref_name(), g3d_particle.get_haplotype(), g3d_particle.get_ref_locus() + max_distance) < (self.legs[index_to_check].get_ref_name(), self.legs[index_to_check].get_haplotype(), self.legs[index_to_check].get_ref_locus())
        
            
    def to_string(self):
        return "\n".join([leg.to_string() for leg in self.legs])

class LegData:
    def __init__(self):
        self.leg_lists = {}
    def num_legs(self):
        num_legs = 0
        for leg_list in self.leg_lists.values():
            num_legs += leg_list.num_legs()
        return num_legs
    def add_empty_leg_list(self, ref_name):
        self.leg_lists[ref_name] = LegList()
    def add_leg(self, leg):
        if leg.get_ref_name() not in self.leg_lists:
            self.add_empty_leg_list(leg.get_ref_name())
        self.leg_lists[leg.get_ref_name()].add_leg(leg)
    def add_con(self, con):
        self.add_leg(con.leg_1())
        self.add_leg(con.leg_2())   
    def add_con_data(self, con_data):
        for con in con_data.get_cons():
            self.add_con(con)
    def sort_legs(self):
        for leg_list in self.leg_lists.values():
            leg_list.sort_legs()
    def sort_phased_legs(self):
        for leg_list in self.leg_lists.values():
            leg_list.sort_phased_legs()
                
    def is_leg_promiscuous(self, leg, max_leg_distance, max_leg_count):
        return self.leg_lists[leg.get_ref_name()].is_leg_promiscuous(leg, max_leg_distance, max_leg_count)
    
    def num_legs_near_g3d_particle(self, g3d_particle, max_distance):
        return self.leg_lists[g3d_particle.get_ref_name()].num_legs_near_g3d_particle(g3d_particle, max_distance)
    def is_g3d_particle_leg_poor(self, g3d_particle, max_distance, min_leg_count):
        return self.leg_lists[g3d_particle.get_ref_name()].is_g3d_particle_leg_poor(g3d_particle, max_distance, min_leg_count)

    def to_string(self):
        return "\n".join([self.leg_lists[ref_name].to_string() for ref_name in sorted(self.leg_lists.keys())])
        
# a contact (legs always sorted)
@total_ordering
class Con:
    def __init__(self, leg_1, leg_2):
        self.legs = sorted([leg_1, leg_2])
    def __repr__(self):
        return self.to_string()    
    def __eq__(self, other):
        return (self.legs[0], self.legs[1]) == (other.legs[0], other.legs[1])
    def __lt__(self, other):
        return (self.legs[0], self.legs[1]) < (other.legs[0], other.legs[1])
    def __hash__(self):
        return hash((self.legs[0], self.legs[1]))
    def deep_copy_from_con(self, other):
        self.legs = copy.deepcopy(other.legs)
    
    def leg_1(self):
        return self.legs[0]
    def leg_2(self):
        return self.legs[1]
    def num_phased_legs(self):
        num_phased_legs = 0
        for i in range(2):
            if self.legs[i].is_phased():
                num_phased_legs += 1
        return num_phased_legs
    def num_conflict_legs(self):
        num_conflict_legs = 0
        for i in range(2):
            if self.legs[i].is_conflict():
                num_conflict_legs += 1
        return num_conflict_legs
    def ref_name_tuple(self):
        return tuple([leg.get_ref_name() for leg in self.legs])
    def hap_tuple(self):
        return tuple([leg.get_haplotype() for leg in self.legs])

    def compatible_hap_tuples(self):
        return [(haplotype_1, haplotype_2) for haplotype_1 in self.legs[0].compatible_haps() for haplotype_2 in self.legs[1].compatible_haps()]

    def compatible_cons(self, is_male, par_data):
        return [Con(leg_1, leg_2) for leg_1 in self.legs[0].compatible_legs(is_male, par_data) for leg_2 in self.legs[1].compatible_legs(is_male, par_data)]
        
    ## impute based on nearby contacts in a ConData
    def votes_from_con_data(self, con_data, max_impute_distance):
        compatible_hap_tuples = self.compatible_hap_tuples()
        voted_hap_tuples = []
        for con in con_data.get_cons_near(self, max_impute_distance):
            voted_hap_tuple = vote_from_hap_tuples(compatible_hap_tuples, con.compatible_hap_tuples())
            if not voted_hap_tuple is None:
                voted_hap_tuples.append(voted_hap_tuple)
        return voted_hap_tuples
    def set_hap_tuple_from_votes(self, voted_hap_tuples, min_impute_votes, min_impute_vote_fraction):
        hap_tuple = winning_vote(voted_hap_tuples, min_impute_votes, min_impute_vote_fraction)
        if not hap_tuple is None:
            self.set_hap_tuple(hap_tuple)
    def impute_from_con_data(self, con_data, max_impute_distance, min_impute_votes, min_impute_vote_fraction, max_intra_hom_separation):
        if self.num_phased_legs() == 2:
            # already phased
            return
        if self.num_phased_legs() == 1 and self.is_intra_chr() and self.separation() <= max_intra_hom_separation:
            # assume intra-homologous
            if self.legs[0].is_phased():
                self.legs[1].set_haplotype(self.legs[0].get_haplotype())
            else:
                self.legs[0].set_haplotype(self.legs[1].get_haplotype())
            return
        self.set_hap_tuple_from_votes(self.votes_from_con_data(con_data, max_impute_distance), min_impute_votes, min_impute_vote_fraction)
        
    ## impute based on 3D structure in a G3dData
    # for a fully phased contact, return whether both legs are out of bounds, and the 3D distance between the two legs
    def distance_in_g3d_data(self, g3d_data):
        is_leg_1_out, leg_1_position = g3d_data.interpolate_leg(self.legs[0])
        is_leg_2_out, leg_2_position = g3d_data.interpolate_leg(self.legs[1])
        is_both_out = is_leg_1_out and is_leg_2_out
        if leg_1_position is None or leg_2_position is None:
            return True, None
        distance = spatial.distance.euclidean(leg_1_position, leg_2_position)
        return is_both_out, distance    
    
    def impute_from_g3d_data(self, g3d_data, max_impute3_distance, max_impute3_ratio, min_impute3_separation, is_male, par_data):
        compatible_cons = self.compatible_cons(is_male, par_data)
        num_compatible_cons = len(compatible_cons)
        if num_compatible_cons == 1:
            # already phased (copy in case of X and Y)
            self.deep_copy_from_con(compatible_cons[0])
            return
        if num_compatible_cons == 4 and self.is_intra_chr() and self.separation() < min_impute3_separation:
            # too close and completely unphased, cannot impute
            return
        con_distance_tuples = []
        any_is_both_out = False
        for con in compatible_cons:
            is_both_out, distance = con.distance_in_g3d_data(g3d_data)
            if is_both_out:
                any_is_both_out = True
            con_distance_tuples.append((con, distance))
        
        if self.is_intra_chr() and any_is_both_out:
            # intra-chromosomal, and both legs are out of range (in any of the compatible), cannot impute (this does not work for CNV loss)
            return
            
        # core imputation
        con_distance_tuples.sort(key = lambda x:x[1])
        impute3_con, impute3_distance = con_distance_tuples[0]
        impute3_ratio = impute3_distance / con_distance_tuples[1][1]
        #if not vio_file is None:
            #vio_file.write("\t".join([self.to_string(), str(num_compatible_cons), str(impute3_distance), str(impute3_ratio]) + "\n")
        if impute3_distance <= max_impute3_distance and impute3_ratio <= max_impute3_ratio:
            self.deep_copy_from_con(impute3_con)
        
        
    def set_hap_tuple(self, hap_tuple):
        for i in range(2):
            self.legs[i].set_haplotype(hap_tuple[i])
    def set_non_par_hap_tuple_male(self, par_data):
        for i in range(2):
            par_data.set_non_par_leg_haplotype_male(self.legs[i])
    def in_par(self, par_data):
        return par_data.contain_leg(self.legs[0]) or par_data.contain_leg(self.legs[1])
    
    def sort_legs(self):
        self.legs.sort()
       
    def is_intra_chr(self):
        return self.leg_1().same_chr_with(self.leg_2())
    def is_inter_hom(self):
        return self.is_intra_chr() and self.num_phased_legs() == 2 and self.leg_1().get_haplotype() != self.leg_2().get_haplotype()
    
    def separation(self):
        return self.leg_2().get_ref_locus() - self.leg_1().get_ref_locus()
        
    def merge_with(self, other):
        for i in range(2):
            self.legs[i].merge_with(other.legs[i])
        self.sort_legs()
    
    # different distance functions w. r. t. another contact, assuming the same chromosome
    def distance_leg_1_with(self, other):
        return self.leg_1().separation_with(other.leg_1())
    def distance_leg_2_with(self, other):
        return self.leg_2().separation_with(other.leg_2())
    def distance_inf_with(self, other): # L-inf norm
        return max(self.distance_leg_1_with(other), self.distance_leg_2_with(other))
    def distance_half_with(self, other): # L-1/2 norm
        return (math.sqrt(self.distance_leg_1_with(other)) + math.sqrt(self.distance_leg_2_with(other))) ** 2
            
    def satisfy_regs(self, inc_regs, exc_regs):
        return self.leg_1().satisfy_regs(inc_regs, exc_regs) and self.leg_2().satisfy_regs(inc_regs, exc_regs)
    
    def is_promiscuous(self, leg_data, max_leg_distance, max_leg_count):
        return leg_data.is_leg_promiscuous(self.leg_1(), max_leg_distance, max_leg_count) or leg_data.is_leg_promiscuous(self.leg_2(), max_leg_distance, max_leg_count)
    def is_isolated(self, con_data, max_clean_distance, min_clean_count):
        num_neighbors = 0
        for con in con_data.get_cons_near(self, max_clean_distance):
            num_neighbors += 1
            if num_neighbors >= min_clean_count + 1: # assume con is included in ConData, thus plus 1
                return False
        return True
    def test_isolated(self, con_data, max_clean_distance, min_clean_count):
        num_neighbors = 0
        for con in con_data.get_cons_near(self, max_clean_distance):
            num_neighbors += 1
            if num_neighbors >= min_clean_count + 1: # assume con is included in ConData, thus plus 1
                break
        return num_neighbors - 1
    def is_isolated_phased(self, con_data, max_clean_distance, min_clean_count):
        num_neighbors = 0
        for con in con_data.get_cons_near(self, max_clean_distance):
            if self.hap_tuple() != con.hap_tuple():
                continue # count only if hap_tuples match
            num_neighbors += 1
            if num_neighbors >= min_clean_count + 1: # assume con is included in ConData, thus plus 1
                return False
        return True
            
    def to_string(self):
        return "\t".join([leg.to_string() for leg in self.legs])
        
    # ard: print relative locus with respect to a reference point
    def to_rel_locus_around(self, other):
        return (self.leg_1().signed_separation_with(other.leg_1())), self.leg_2().signed_separation_with(other.leg_2())
    def to_string_around(self, other):
        return str(self.leg_1().signed_separation_with(other.leg_1())) + "\t" + str(self.leg_2().signed_separation_with(other.leg_2()))
        
def ref_name_tuple_to_string(ref_name_tuple):
    return ",".join(ref_name_tuple)

def string_to_con(con_string):
    leg_1, leg_2 = con_string.split("\t")
    return Con(string_to_leg(leg_1), string_to_leg(leg_2))

## voting during imputation
# a list (self) being voted by another list (other), return None if not compatible or more than one compatible
def vote_from_hap_tuples(self_hap_tuples, other_hap_tuples):
    compatible_hap_tuples = set(self_hap_tuples) & set(other_hap_tuples)
    if len(compatible_hap_tuples) == 1:
        return compatible_hap_tuples.pop()
    return None
# return the voted winner, None if not meeting criteria
def winning_vote(hap_tuples, min_impute_votes, min_impute_vote_fraction):
    if len(hap_tuples) == 0:
        return None
    hap_tuple_counter = Counter(hap_tuples)
    hap_tuple, vote = hap_tuple_counter.most_common(1)[0]
    if vote >= min_impute_votes and float(vote) / len(hap_tuples) >= min_impute_vote_fraction:
        return hap_tuple
    return None

# a sorted list of contacts
class ConList:
    def __init__(self):
        self.cons = []
    
    # generator for all its contacts
    def get_cons(self):
        for con in self.cons:
            yield con
    
    # generator for contacts near a given contact (L-1/2 norm), assume sorted
    def get_cons_near(self, con, max_distance):
        start_index = bisect.bisect_left(self.cons, con)
        for i in range(start_index, len(self.cons)):
            if self.cons[i].distance_leg_1_with(con) > max_distance:
                break
            if self.cons[i].distance_half_with(con) <= max_distance:
                yield self.cons[i]
        for i in range(start_index - 1, -1, -1):
            if self.cons[i].distance_leg_1_with(con) > max_distance:
                break
            if self.cons[i].distance_half_with(con) <= max_distance:
                yield self.cons[i]            
    # generator for contacts near a given contact (L-inf norm), assume sorted
    def get_cons_near_inf(self, con, max_distance):
        start_index = bisect.bisect_left(self.cons, con)
        for i in range(start_index, len(self.cons)):
            if self.cons[i].distance_leg_1_with(con) > max_distance:
                break
            if self.cons[i].distance_leg_2_with(con) <= max_distance:
                yield self.cons[i]
        for i in range(start_index - 1, -1, -1):
            if self.cons[i].distance_leg_1_with(con) > max_distance:
                break
            if self.cons[i].distance_leg_2_with(con) <= max_distance:
                yield self.cons[i]         
                        
    def num_cons(self):
        return(len(self.cons))
    def num_phased_legs(self):
        num_phased_legs = 0
        for con in self.cons:
            num_phased_legs += con.num_phased_legs()
        return num_phased_legs
    def num_phased_cons(self):
        num_phased_con = 0
        for con in self.cons:
            if con.num_phased_legs() == 2:
                num_phased_con += 1
        return num_phased_con
    def num_conflict_legs(self):
        num_conflict_legs = 0
        for con in self.cons:
            num_conflict_legs += con.num_conflict_legs()
        return num_conflict_legs
    def num_intra_chr(self):
        num_intra_chr = 0
        for con in self.cons:
            if con.is_intra_chr():
                num_intra_chr += 1
        return num_intra_chr
                        
    def sort_cons(self):
        self.cons.sort()
    
    def add_con(self, con):
        self.cons.append(con)
    
    def merge_with(self, other):
        self.cons += other.cons
        
    # remove intra-chromosomal contacts with small separations
    def clean_separation(self, min_separation):
        self.cons[:] = [con for con in self.cons if not con.is_intra_chr() or con.separation() >= min_separation]
    # remove contacts containing promiscuous legs
    def clean_promiscuous(self, leg_data, max_leg_distance, max_leg_count):
        self.cons[:] = [con for con in self.cons if not con.is_promiscuous(leg_data, max_leg_distance, max_leg_count)]
    # remove isolated contacts
    def clean_isolated(self, con_data, max_clean_distance, min_clean_count):
        self.cons[:] = [con for con in self.cons if not con.is_isolated(con_data, max_clean_distance, min_clean_count)]
    def test_isolated(self, con_data, max_clean_distance, min_clean_count):
        return [con.test_isolated(con_data, max_clean_distance, min_clean_count) for con in self.cons]
    # remove contacts in PARs
    def clean_in_par(self, par_data):
        self.cons[:] = [con for con in self.cons if not con.in_par(par_data)]
    # set haplotypes in haploid regions of male
    def set_non_par_hap_tuple_male(self, par_data):
        for con in self.cons:
            con.set_non_par_hap_tuple_male(par_data)
    # remove inter-homologous contacts with small separations
    def clean_separation_hom(self, min_separation):
        self.cons[:] = [con for con in self.cons if not con.is_inter_hom() or con.separation() >= min_separation]
    # impute
    def impute_from_con_data(self, con_data, max_impute_distance, min_impute_votes, min_impute_vote_fraction, max_intra_hom_separation, min_inter_hom_separation):
        for con in self.cons:
            con.impute_from_con_data(con_data, max_impute_distance, min_impute_votes, min_impute_vote_fraction, max_intra_hom_separation)
        self.clean_separation_hom(min_inter_hom_separation)
    # remove everything but fully phased contacts
    def clean_unphased(self):
        self.cons[:] = [con for con in self.cons if con.num_phased_legs() == 2]
    # remove isolated contacts for fully phased contacts
    def clean_isolated_phased(self, con_data, max_clean_distance, min_clean_count):
        self.cons[:] = [con for con in self.cons if not con.is_isolated_phased(con_data, max_clean_distance, min_clean_count)]
    # impute3
    def impute_from_g3d_data(self, g3d_data, max_impute3_distance, max_impute3_ratio, min_impute3_separation, is_male, par_data):
        for con in self.cons:
            con.impute_from_g3d_data(g3d_data, max_impute3_distance, max_impute3_ratio, min_impute3_separation, is_male, par_data)

        
    # simple dedup within a read (no binary search), assuming the same chromosome
    def dedup_within_read(self, max_distance):
        while True:
            merged = False
            for i in range(len(self.cons)):
                for j in range(i + 1, len(self.cons)):
                    if self.cons[i].distance_inf_with(self.cons[j]) <= max_distance:
                        self.cons[i].merge_with(self.cons[j])
                        self.cons.pop(j)
                        merged = True
                        break
            if merged == False:
                break
    
    # faster dedup, assuming the same chromosome
    def dedup(self, max_distance):
        self.cons.sort()
        while True:
            merged = False
            for i in range(len(self.cons)):
                for j in range(i + 1, len(self.cons)):
                    if self.cons[i].distance_leg_1_with(self.cons[j]) > max_distance:
                        break
                    if self.cons[i].distance_leg_2_with(self.cons[j]) <= max_distance:
                        self.cons[i].merge_with(self.cons[j])
                        self.cons.pop(j)
                        merged = True
                        break
            if merged == False:
                break
            self.cons[i:j] = sorted(self.cons[i:j])
    
    def apply_regs(self, inc_regs, exc_regs):
        self.cons[:] = [con for con in self.cons if con.satisfy_regs(inc_regs, exc_regs)]
        
    def to_string(self):
        return "\n".join([con.to_string() for con in self.cons])

        
# a hashmap (tuples of two sorted chromosome names) of lists of contacts (a CON file)
class ConData:
    def __init__(self):
        self.con_lists = {}
        
    # generator for all its cons, with ref_name_tuple sorted
    def get_cons(self):
        for ref_name_tuple in sorted(self.con_lists.keys()):
            for con in self.con_lists[ref_name_tuple].get_cons():
                yield con
    
    # wrappers for generators
    def get_cons_near(self, con, max_distance):
        if con.ref_name_tuple() in self.con_lists:
            for con in self.con_lists[con.ref_name_tuple()].get_cons_near(con, max_distance):
                yield con
    def get_cons_near_inf(self, con, max_distance):
        if con.ref_name_tuple() in self.con_lists:
            for con in self.con_lists[con.ref_name_tuple()].get_cons_near_inf(con, max_distance):
                yield con
                    
    def add_empty_con_list(self, ref_name_tuple):
        self.con_lists[ref_name_tuple] = ConList()
    
    def add_con(self, con):
        if con.ref_name_tuple() not in self.con_lists:
            self.add_empty_con_list(con.ref_name_tuple())
        self.con_lists[con.ref_name_tuple()].add_con(con)
    
    def merge_with(self, other):
        for ref_name_tuple in other.con_lists.keys():
            if ref_name_tuple in self.con_lists:
                self.con_lists[ref_name_tuple].merge_with(other.con_lists[ref_name_tuple])
            else:
                self.con_lists[ref_name_tuple] = other.con_lists[ref_name_tuple]
    
    # remove all intra-chromosomal contacts
    def clean_intra_chr(self):
        for ref_name_tuple in self.con_lists.keys():
            if ref_name_tuple[0] == ref_name_tuple[1]:
                del self.con_lists[ref_name_tuple]

    # remove all inter-chromosomal contacts
    def clean_inter_chr(self):
        for ref_name_tuple in self.con_lists.keys():
            if ref_name_tuple[0] != ref_name_tuple[1]:
                del self.con_lists[ref_name_tuple]
        
    # wrappers for all ConList operations
    def sort_cons(self):
        for con_list in self.con_lists.values():
            con_list.sort_cons()
    def clean_separation(self, min_separation):
        for ref_name_tuple in self.con_lists.keys():
            self.con_lists[ref_name_tuple].clean_separation(min_separation)
            if self.con_lists[ref_name_tuple].num_cons() == 0:
                del self.con_lists[ref_name_tuple]
    def clean_promiscuous(self, leg_data, max_leg_distance, max_leg_count):
        for ref_name_tuple in self.con_lists.keys():
            self.con_lists[ref_name_tuple].clean_promiscuous(leg_data, max_leg_distance, max_leg_count)
            if self.con_lists[ref_name_tuple].num_cons() == 0:
                del self.con_lists[ref_name_tuple]
    def clean_isolated(self, con_data, max_clean_distance, min_clean_count):
        for ref_name_tuple in self.con_lists.keys():
            original_num_cons = self.con_lists[ref_name_tuple].num_cons()
            self.con_lists[ref_name_tuple].clean_isolated(con_data, max_clean_distance, min_clean_count)
            sys.stderr.write("[M::" + __name__ + "] cleaned isolated contacts for chromosome pair (" + ref_name_tuple_to_string(ref_name_tuple) + "): " + str(original_num_cons) + " -> " + str(self.con_lists[ref_name_tuple].num_cons()) + " contacts\n")
            if self.con_lists[ref_name_tuple].num_cons() == 0:
                del self.con_lists[ref_name_tuple]
    def test_isolated(self, con_data, max_clean_distance, min_clean_count):
        neighbor_counts = []
        for ref_name_tuple in self.con_lists.keys():
            original_num_cons = self.con_lists[ref_name_tuple].num_cons()
            neighbor_counts.extend(self.con_lists[ref_name_tuple].test_isolated(con_data, max_clean_distance, min_clean_count))
            sys.stderr.write("[M::" + __name__ + "] tested isolated contacts for chromosome pair (" + ref_name_tuple_to_string(ref_name_tuple) + "): " + str(original_num_cons) + " contacts\n")
        return neighbor_counts
    def clean_in_par(self, par_data):
        for ref_name_tuple in self.con_lists.keys():
            if par_data.get_x_name() not in ref_name_tuple and par_data.get_y_name() not in ref_name_tuple:
                continue
            self.con_lists[ref_name_tuple].clean_in_par(par_data)
            if self.con_lists[ref_name_tuple].num_cons() == 0:
                del self.con_lists[ref_name_tuple]
    def clean_unphased(self):
        for ref_name_tuple in self.con_lists.keys():
            self.con_lists[ref_name_tuple].clean_unphased()
            if self.con_lists[ref_name_tuple].num_cons() == 0:
                del self.con_lists[ref_name_tuple]
    def clean_isolated_phased(self, con_data, max_clean_distance, min_clean_count):
        for ref_name_tuple in self.con_lists.keys():
            original_num_cons = self.con_lists[ref_name_tuple].num_cons()
            self.con_lists[ref_name_tuple].clean_isolated_phased(con_data, max_clean_distance, min_clean_count)
            #sys.stderr.write("[M::" + __name__ + "] cleaned isolated contacts for chromosome pair (" + ref_name_tuple_to_string(ref_name_tuple) + "): " + str(original_num_cons) + " -> " + str(self.con_lists[ref_name_tuple].num_cons()) + " contacts\n")
            if self.con_lists[ref_name_tuple].num_cons() == 0:
                del self.con_lists[ref_name_tuple]
    def set_non_par_hap_tuple_male(self, par_data):
        for ref_name_tuple in self.con_lists.keys():
            if par_data.get_x_name() not in ref_name_tuple and par_data.get_y_name() not in ref_name_tuple:
                continue
            self.con_lists[ref_name_tuple].set_non_par_hap_tuple_male(par_data)
    def impute_from_con_data(self, con_data, max_impute_distance, min_impute_votes, min_impute_vote_fraction, max_intra_hom_separation, min_inter_hom_separation):
        for ref_name_tuple in self.con_lists.keys():
            self.con_lists[ref_name_tuple].impute_from_con_data(con_data, max_impute_distance, min_impute_votes, min_impute_vote_fraction, max_intra_hom_separation, min_inter_hom_separation)
            #sys.stderr.write("[M::" + __name__ + "] imputed haplotypes for chromosome pair (" + ref_name_tuple_to_string(ref_name_tuple) + "): " + str(self.con_lists[ref_name_tuple].num_cons()) + " contacts\n")
    def impute_from_g3d_data(self, g3d_data, max_impute3_distance, max_impute3_ratio, min_impute3_separation, is_male, par_data):
        for ref_name_tuple in self.con_lists.keys():
            self.con_lists[ref_name_tuple].impute_from_g3d_data(g3d_data, max_impute3_distance, max_impute3_ratio, min_impute3_separation, is_male, par_data)
            sys.stderr.write("[M::" + __name__ + "] imputed haplotypes for chromosome pair (" + ref_name_tuple_to_string(ref_name_tuple) + "): " + str(self.con_lists[ref_name_tuple].num_cons()) + " contacts (" + str(round(100.0 * self.con_lists[ref_name_tuple].num_phased_cons() / self.con_lists[ref_name_tuple].num_cons(), 2)) + "% phased)\n")
    def dedup_within_read(self, max_distance):
        for con_list in self.con_lists.values():
            con_list.dedup_within_read(max_distance)
    def dedup(self, max_distance):
        for ref_name_tuple in self.con_lists.keys():
            sys.stderr.write("[M::" + __name__ + "] merging duplicates for chromosome pair (" + ref_name_tuple_to_string(ref_name_tuple) + "): " + str(self.con_lists[ref_name_tuple].num_cons()) + " putative contacts\n")
            self.con_lists[ref_name_tuple].dedup(max_distance)
    def apply_regs(self, inc_regs, exc_regs):
        for ref_name_tuple in self.con_lists.keys():
            self.con_lists[ref_name_tuple].apply_regs(inc_regs, exc_regs)
            if self.con_lists[ref_name_tuple].num_cons() == 0:
                del self.con_lists[ref_name_tuple]
    def num_cons(self):
        num_cons = 0
        for con_list in self.con_lists.values():
            num_cons += con_list.num_cons()
        return num_cons
    def num_phased_cons(self):
        num_cons = 0
        for con_list in self.con_lists.values():
            num_cons += con_list.num_phased_cons()
        return num_cons
    def num_phased_legs(self):
        num_phased_legs = 0
        for con_list in self.con_lists.values():
            num_phased_legs += con_list.num_phased_legs()
        return num_phased_legs
    def num_conflict_legs(self):
        num_conflict_legs = 0
        for con_list in self.con_lists.values():
            num_conflict_legs += con_list.num_conflict_legs()
        return num_conflict_legs
    def num_intra_chr(self):
        num_intra_chr = 0
        for con_list in self.con_lists.values():
            num_intra_chr += con_list.num_intra_chr()
        return num_intra_chr
    
 
                         
    def to_string(self): # no tailing new line
        return "\n".join([self.con_lists[ref_name_tuple].to_string() for ref_name_tuple in sorted(self.con_lists.keys())])

def file_to_con_data(con_file):
    con_data = ConData()
    for con_file_line in con_file:
        con_data.add_con(string_to_con(con_file_line.strip()))
    return con_data

# augmented data for dedup: each leg records haplotypes of all duplicates
class DupLeg(Leg):
    def __init__(self, leg):
        Leg.__init__(self, leg.ref_name, leg.ref_locus, leg.haplotype)
        self.dups = {Haplotypes.unknown: 0, Haplotypes.paternal: 0, Haplotypes.maternal: 0}
        self.dups[leg.haplotype] += 1
    def num_dups(self):
        return sum(self.dups.values())
    def merge_with(self, other):
        Leg.merge_with(self, other)
        for haplotype in self.dups.keys():
            self.dups[haplotype] += other.dups[haplotype]
    #def to_string(self):
        #if self.haplotype != Haplotypes.conflict:
            #return ""
        #return Leg.to_string(self) + "(" + ",".join([str(self.dups[Haplotypes.unknown]), str(self.dups[Haplotypes.paternal]), str(self.dups[Haplotypes.maternal])]) + ")"

class DupCon(Con):
    def __init__(self, con):
        Con.__init__(self, DupLeg(con.legs[0]), DupLeg(con.legs[1]))
    def num_dups(self):
        return self.legs[0].num_dups() # should be the same for both legs

class DupConList(ConList):
    def __init__(self, con_list):
        ConList.__init__(self)
        for con in con_list.cons:
            self.add_con(DupCon(con))
    def dup_stats(self, display_max_num_dups):
        hist_num_dups = [0] * display_max_num_dups
        for dup_con in self.cons:
            hist_num_dups[min(dup_con.num_dups(), display_max_num_dups) - 1] += 1
        return hist_num_dups
        
class DupConData(ConData):
    def __init__(self, con_data):
        ConData.__init__(self)
        for ref_name_tuple in con_data.con_lists.keys():
            self.con_lists[ref_name_tuple] = DupConList(con_data.con_lists[ref_name_tuple])
            
    def dup_stats(self, display_max_num_dups):
        hist_num_dups = [0] * display_max_num_dups
        for con_list in self.con_lists.values():
            list_hist_num_dups = con_list.dup_stats(display_max_num_dups)
            for i in range(display_max_num_dups):
                hist_num_dups[i] += list_hist_num_dups[i]
        return hist_num_dups
        
    def add_empty_con_list(self, ref_name_tuple):
        self.con_lists[ref_name_tuple] = DupConList()
        
# print a histogram of counts to a string
def counts_to_hist_num_with_zero(counts):
    hist_num = [0] * (max(counts) + 1)
    for count in counts:
        hist_num[count] += 1
    return hist_num
def hist_num_to_string(hist_num):
    return "\n".join([(">=" if i == len(hist_num) - 1 else "") + str(i + 1) + "\t" + str(hist_num[i]) + " (" + str(round(100.0 * hist_num[i] / sum(hist_num), 2))+ "%)" for i in range(len(hist_num))])        
def hist_num_to_string_with_zero(hist_num):
    return "\n".join([(">=" if i == len(hist_num) - 1 else "") + str(i) + "\t" + str(hist_num[i]) + " (" + str(round(100.0 * hist_num[i] / sum(hist_num), 2))+ "%)" for i in range(len(hist_num))])        


# structures for included and excluded regions
class Reg:
    def __init__(self, ref_name):
        self.ref_name = ref_name
        self.has_haplotype = False
        self.has_start = False
        self.has_end = False
        self.haplotype = Haplotypes.unknown
        self.start = -1
        self.end = -1
    def add_haplotype(self, haplotype):
        if is_known_haplotype(haplotype):
            self.has_haplotype = True
            self.haplotype = haplotype
    def add_start(self, start):
        self.has_start = True
        self.start = start
    def add_end(self, end):
        self.has_end = True
        self.end = end
    def to_string(self):
        return "\t".join([self.ref_name, (haplotype_to_string(self.haplotype) if self.has_haplotype else "."), (str(self.start) if self.has_start else "."), (str(self.end) if self.has_end else ".")])


# structures for PARs
class Par:
    def __init__(self, x_name, x_start, x_end, y_name, y_start):
        self.x_name = x_name
        self.y_name = y_name
        self.x_start = x_start
        self.x_end = x_end
        self.y_start = y_start
        
        # derived quantities
        self.shift = y_start - x_start # Y relative to X
        
        self.x_reg = Reg(x_name)
        self.x_reg.add_start(self.x_start)
        self.x_reg.add_end(self.x_end)
        
        self.y_reg = Reg(x_name)
        self.y_reg.add_start(self.y_start)
        self.y_reg.add_end(self.y_start + self.x_end - self.x_start)
    def contain_leg(self, leg):
        return leg.in_reg(self.x_reg) or leg.in_reg(self.y_reg)
    def compatible_par_legs_male(self, leg):
        if leg.in_reg(self.x_reg):
            if leg.get_haplotype() == Haplotypes.paternal:
                # shift to Y
                return [Leg(self.y_name, leg.get_ref_locus() + self.shift, Haplotypes.paternal)]
            if leg.get_haplotype() == Haplotypes.maternal:
                # keep on X
                return [Leg(self.x_name, leg.get_ref_locus(), Haplotypes.maternal)]
            # both possible
            return [Leg(self.x_name, leg.get_ref_locus(), Haplotypes.maternal), Leg(self.y_name, leg.get_ref_locus() + self.shift, Haplotypes.paternal)]

        if leg.in_reg(self.y_reg):
            if leg.get_haplotype() == Haplotypes.paternal:
                # keep on Y
                return [Leg(self.y_name, leg.get_ref_locus(), Haplotypes.paternal)]
            if leg.get_haplotype() == Haplotypes.maternal:
                # shift to X
                return [Leg(self.x_name, leg.get_ref_locus() - self.shift, Haplotypes.maternal)]
            # both possible
            return [Leg(self.y_name, leg.get_ref_locus(), Haplotypes.paternal), Leg(self.x_name, leg.get_ref_locus() - self.shift, Haplotypes.maternal)]
        return None # leg not in this PAR

class ParData:
    def __init__(self, x_name, y_name):
        self.x_name = x_name
        self.y_name = y_name
        self.pars = []
    def add_par(self, par):
        self.pars.append(par)
    def get_pars(self):
        for par in self.pars:
            yield par
    def get_x_name(self):
        return self.x_name
    def get_y_name(self):
        return self.y_name
    def contain_leg(self, leg):
        for par in self.pars:
            if par.contain_leg(leg):
                return True
        return False
    def set_non_par_leg_haplotype_male(self, leg):
        if leg.get_ref_name() == self.x_name:
            leg.set_haplotype(Haplotypes.maternal)
        elif leg.get_ref_name() == self.y_name:
            leg.set_haplotype(Haplotypes.paternal)
    def compatible_legs_male(self, leg):
        for par in self.pars:
            legs = par.compatible_par_legs_male(leg)
            if not legs is None:
                return legs
        if leg.get_ref_name() == self.x_name:
            return [Leg(self.x_name, leg.get_ref_locus(), Haplotypes.maternal)]
        if leg.get_ref_name() == self.y_name:
            return [Leg(self.y_name, leg.get_ref_locus(), Haplotypes.paternal)]
        return leg.compatible_legs_female()

# structures for 3D genome (3DG) files
@total_ordering
class G3dParticle:
    def __init__(self, ref_name, ref_locus, haplotype, position):
        self.ref_name = ref_name
        self.ref_locus = ref_locus
        self.haplotype = haplotype
        self.position = position
    def __eq__(self, other):
        return (self.ref_name, self.haplotype, self.ref_locus) == (other.ref_name, other.haplotype, other.ref_locus)
    def __lt__(self, other):
        return (self.ref_name, self.haplotype, self.ref_locus) < (other.ref_name, other.haplotype, other.ref_locus)
    def __repr__(self):
        return self.to_string()
    def to_string(self):
        return "\t".join([ref_name_haplotype_to_hom_name(self.ref_name, self.haplotype), str(self.ref_locus)] + map(str, self.position))
    def get_ref_name(self):
        return self.ref_name
    def get_haplotype(self):
        return self.haplotype
    def get_ref_locus(self):
        return self.ref_locus
    def get_position(self):
        return self.position
    def get_x(self):
        return self.position[0]
    def get_y(self):
        return self.position[1]
    def get_z(self):
        return self.position[2]
    def ref_name_haplotype(self):
        return (self.ref_name, self.haplotype)
    def hom_name(self):
        return ref_name_haplotype_to_hom_name(self.ref_name, self.haplotype)
        
def string_to_g3d_particle(g3d_particle_string):
    hom_name, ref_locus, x, y, z = g3d_particle_string.split("\t")
    ref_name, haplotype = hom_name_to_ref_name_haplotype(hom_name)
    ref_locus = int(ref_locus)
    position = [float(x), float(y), float(z)]
    return G3dParticle(ref_name, ref_locus, haplotype, position)
    
class G3dList:
    def __init__(self):
        self.g3d_particles = []
        
        # for linear interpolation
        self.min_ref_locus = None
        self.max_ref_locus = None
        self.interpolate = None # scipy function for locus -> position
    def num_g3d_particles(self):
        return len(self.g3d_particles)
    def add_g3d_particle(self, g3d_particle):
        self.g3d_particles.append(g3d_particle)
    def sort_g3d_particles(self):
        self.g3d_particles.sort()
    def prepare_interpolate(self):
        self.min_ref_locus = min([g3d_particle.get_ref_locus() for g3d_particle in self.g3d_particles])
        self.max_ref_locus = max([g3d_particle.get_ref_locus() for g3d_particle in self.g3d_particles])
        self.interpolate = interpolate.interp1d(np.array([g3d_particle.get_ref_locus() for g3d_particle in self.g3d_particles]), np.array([g3d_particle.get_position() for g3d_particle in self.g3d_particles]),axis=0)
    def interpolate_ref_locus(self, ref_locus):
        is_out = ref_locus < self.min_ref_locus or ref_locus > self.max_ref_locus # if out of bounds
        position = self.interpolate(min(max(ref_locus, self.min_ref_locus), self.max_ref_locus))
        return is_out, position
    def to_string(self):
        return "\n".join([g3d_particle.to_string() for g3d_particle in self.g3d_particles])
    
    # return a list of locus increments, for inferring resolution, must be sorted
    def ref_locus_increments(self):
        return [self.g3d_particles[i + 1].get_ref_locus() - self.g3d_particles[i].get_ref_locus() for i in range(len(self.g3d_particles) - 1)]
    
    # return a list of leg counts
    def leg_counts(self, leg_data, max_distance):
        return [leg_data.num_legs_near_g3d_particle(g3d_particle, max_distance) for g3d_particle in self.g3d_particles]
    
    def clean_leg_poor(self, leg_data, max_distance, min_leg_count):
        self.g3d_particles[:] = [g3d_particle for g3d_particle in self.g3d_particles if not leg_data.is_g3d_particle_leg_poor(g3d_particle, max_distance, min_leg_count)]

    def get_g3d_particles(self):
        for g3d_particle in self.g3d_particles:
            yield g3d_particle
    def get_adjacent_g3d_particle_tuples(self, resolution): # assume sorted
        for i in range(len(self.g3d_particles) - 1):
            g3d_particle_tuple = (self.g3d_particles[i], self.g3d_particles[i + 1])
            if g3d_particle_tuple[1].get_ref_locus() - g3d_particle_tuple[0].get_ref_locus() == resolution:
                yield g3d_particle_tuple

class G3dData:
    def __init__(self):
        self.g3d_lists = {}
    def add_empty_g3d_list(self, ref_name_haplotype):
        self.g3d_lists[ref_name_haplotype] = G3dList()
    def add_g3d_particle(self, g3d_particle):
        if g3d_particle.ref_name_haplotype() not in self.g3d_lists:
            self.add_empty_g3d_list(g3d_particle.ref_name_haplotype())
        self.g3d_lists[g3d_particle.ref_name_haplotype()].add_g3d_particle(g3d_particle)
    
    # infer resolution, must be sorted
    def ref_locus_increments(self):
        ref_locus_increments = []
        for g3d_list in self.g3d_lists.values():
            ref_locus_increments.extend(g3d_list.ref_locus_increments())
        return ref_locus_increments
    def resolution(self):
        ref_locus_increments_counter = Counter(self.ref_locus_increments())
        ref_locus_increment = ref_locus_increments_counter.most_common(1)[0][0]
        return ref_locus_increment
        
    # linear interpolation, return whether out of bounds, and position (None if no homologs are found)
    def interpolate_leg(self, leg):
        if leg.ref_name_haplotype() in self.g3d_lists:
            return self.g3d_lists[leg.ref_name_haplotype()].interpolate_ref_locus(leg.get_ref_locus())
        return True, None

        
    # wrappers for G3dList functions:
    def num_g3d_particles(self):
        num_g3d_particles = 0
        for g3d_list in self.g3d_lists.values():
            num_g3d_particles += g3d_list.num_g3d_particles()
        return num_g3d_particles
    def sort_g3d_particles(self):
        for g3d_list in self.g3d_lists.values():
            g3d_list.sort_g3d_particles()
    def prepare_interpolate(self):
        for g3d_list in self.g3d_lists.values():
            g3d_list.prepare_interpolate()
    def leg_counts(self, leg_data, max_distance):
        leg_counts = []
        for g3d_list in self.g3d_lists.values():
            leg_counts.extend(g3d_list.leg_counts(leg_data, max_distance))
        return leg_counts
    def clean_leg_poor(self, leg_data, max_distance, min_leg_count):
        for g3d_list in self.g3d_lists.values():
            g3d_list.clean_leg_poor(leg_data, max_distance, min_leg_count)
    def get_g3d_particles(self):
        for ref_name_haplotype in sorted(self.g3d_lists.keys()):
            for g3d_particle in self.g3d_lists[ref_name_haplotype].get_g3d_particles():
                yield g3d_particle
    def get_adjacent_g3d_particle_tuples(self, resolution): # assume sorted
        for ref_name_haplotype in sorted(self.g3d_lists.keys()):
            for g3d_particle_tuple in self.g3d_lists[ref_name_haplotype].get_adjacent_g3d_particle_tuples(resolution):
                yield g3d_particle_tuple
    def to_string(self):
        return "\n".join([self.g3d_lists[ref_name_haplotype].to_string() for ref_name_haplotype in sorted(self.g3d_lists.keys())])

def file_to_g3d_data(g3d_file):
    g3d_data = G3dData()
    for g3d_file_line in g3d_file:
        g3d_data.add_g3d_particle(string_to_g3d_particle(g3d_file_line.strip()))
    return g3d_data
