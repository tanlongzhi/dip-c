import bisect
import sys

# enum for haplotypes
class Haplotypes:
    conflict = -2
    unknown = -1
    paternal = 0
    maternal = 1
def is_known_haplotype(haplotype):
    return haplotype >= 0
# the operations below will set "." for all unphased haplotypes
def haplotype_to_string(haplotype):
    return str(haplotype) if is_known_haplotype(haplotype) else "."
def string_to_haplotype(haplotype_string):
    return Haplotypes.unknown if haplotype_string == "." else int(haplotype_string) 

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
        if not self.is_read2 and other.is_read2: # read 1 < read 2
            return True
        if not self.is_read2 and not other.is_read2: # both read 1
            if self.query_start < other.query_start:
                return True
        if self.is_read2 and other.is_read2: # both read 2
            if self.query_end > other.query_end:
                return True
        return False
        
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
class Leg:
    
    def __init__(self, ref_name, ref_locus, haplotype):
        self.ref_name = ref_name
        self.ref_locus = ref_locus
        self.haplotype = haplotype
        
    def __lt__(self, other):
        if self.ref_name < other.ref_name:
            return True
        if self.ref_name == other.ref_name and self.ref_locus < other.ref_locus:
            return True
        if self.ref_name == other.ref_name and self.ref_locus == other.ref_locus and self.haplotype < other.haplotype:
            return True
        return False
        
    def get_ref_name(self):
        return self.ref_name

    def get_ref_locus(self):
        return self.ref_locus
    
    def is_phased(self):
        return is_known_haplotype(self.haplotype)
        
    def merge_with(self, other):
        self.ref_locus = (self.ref_locus + other.ref_locus)/2
        self.haplotype = update_haplotype(self.haplotype, other.haplotype)

    def same_chr_with(self, other):
        return self.ref_name == other.ref_name
        
    def separation_with(self, other):
        return abs(self.ref_locus - other.ref_locus)
    
    def to_string(self):
        return ",".join([self.ref_name, str(self.ref_locus), haplotype_to_string(self.haplotype)])

def string_to_leg(leg_string):
    ref_name, ref_locus, haplotype = leg_string.split(",")
    ref_locus = int(ref_locus)
    haplotype = string_to_haplotype(haplotype)
    return Leg(ref_name, ref_locus, haplotype)


# a contact (legs always sorted)
class Con:
    def __init__(self, leg_1, leg_2):
        self.legs = sorted([leg_1, leg_2])
    
    def __lt__(self, other):
        return self.legs[0] < other.legs[0] or self.legs[0] == other.legs[0] and self.legs[1] < other.legs[1]
    
    def leg_1(self):
        return self.legs[0]
    def leg_2(self):
        return self.legs[1]
    def num_phased_legs(self):
        num_phased_legs = 0
        for i in range(2):
            num_phased_legs += 1 if self.legs[i].is_phased() else 0
        return num_phased_legs
    def ref_names(self):
        return tuple([leg.get_ref_name() for leg in self.legs])
    
    def sort_legs(self):
        self.legs.sort()
       
    def is_intra_chr(self):
        return self.leg_1().same_chr_with(self.leg_2())
    
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
            
    def to_string(self):
        return "\t".join([leg.to_string() for leg in self.legs])
def ref_names_to_string(ref_names):
    return ",".join(ref_names)

def string_to_con(con_string):
    leg_1, leg_2 = con_string.split("\t")
    return Con(string_to_leg(leg_1), string_to_leg(leg_2))

# a sorted list of contacts
class ConList:
    def __init__(self):
        self.cons = []
        self.is_sorted = True
        
    def num_cons(self):
        return(len(self.cons))
        
    def num_phased_legs(self):
        num_phased_legs = 0
        for con in self.cons:
            num_phased_legs += con.num_phased_legs()
        return num_phased_legs
        
    def sort_cons(self):
        self.cons.sort()
        self.is_sorted = True
    
    def add_con(self, con):
        self.cons.append(con)
        self.is_sorted = False
    
    def merge_with(self, other):
        self.cons += other.cons
        if other.num_cons() > 0:
            self.is_sorted = False
        
    # remove intra-chromosomal contacts with small separations, no sorting needed
    def clean_separation(self, min_separation):
        self.cons[:] = [con for con in self.cons if not con.is_intra_chr() or con.separation() > min_separation]

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
        self.is_sorted = False
    
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
        self.is_sorted = True
        
    def to_string(self):
        return "\n".join([con.to_string() for con in self.cons])
        
# a hashmap (tuples of two sorted chromosome names) of lists of contacts (a CON file)
class ConData:
    def __init__(self):
        self.con_lists = {}
        self.is_sorted = True
    
    def add_empty_con_list(self, ref_names):
        self.con_lists[ref_names] = ConList()
    
    def add_con(self, con):
        if con.ref_names() not in self.con_lists:
            self.add_empty_con_list(con.ref_names())
        self.con_lists[con.ref_names()].add_con(con)
        self.is_sorted = False
    
    def merge_with(self, other):
        for ref_names in other.con_lists.keys():
            if ref_names in self.con_lists:
                self.con_lists[ref_names].merge_with(other.con_lists[ref_names])
                if not self.con_lists[ref_names].is_sorted:
                    self.is_sorted = False
            else:
                self.con_lists[ref_names] = other.con_lists[ref_names]
        
    # wrappers for all ConList operations
    def sort_cons(self):
        for con_list in self.con_lists.values():
            con_list.sort_cons()
        self.is_sorted = True
    def clean_separation(self, min_separation):
        for ref_names in self.con_lists.keys():
            self.con_lists[ref_names].clean_separation(min_separation)
            if self.con_lists[ref_names].num_cons() == 0:
                del self.con_lists[ref_names]
    def dedup_within_read(self, max_distance):
        for con_list in self.con_lists.values():
            con_list.dedup_within_read(max_distance)
        self.is_sorted = False
    def dedup(self, max_distance):
        for ref_names in self.con_lists.keys():
            sys.stderr.write("[M::" + __name__ + "] merging duplicates for chromosome pair (" + ref_names_to_string(ref_names) + "): " + str(self.con_lists[ref_names].num_cons()) + " putative contacts\n")
            self.con_lists[ref_names].dedup(max_distance)
        self.is_sorted = True
    def num_cons(self):
        num_cons = 0
        for con_list in self.con_lists.values():
            num_cons += con_list.num_cons()
        return num_cons
    def num_phased_legs(self):
        num_phased_legs = 0
        for con_list in self.con_lists.values():
            num_phased_legs += con_list.num_phased_legs()
        return num_phased_legs
    
    def to_string(self): # no tailing new line
        return "\n".join([self.con_lists[ref_names].to_string() for ref_names in sorted(ref_names for ref_names in self.con_lists.keys())])

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
    def to_string(self):
        return Leg.to_string(self) + "(" + ",".join([str(self.dups[Haplotypes.unknown]), str(self.dups[Haplotypes.paternal]), str(self.dups[Haplotypes.maternal])]) + ")"

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
    def dup_stats(self, max_num_dups):
        hist_num_dups = [0] * max_num_dups
        for dup_con in self.cons:
            hist_num_dups[min(dup_con.num_dups(), max_num_dups) - 1] += 1
        return hist_num_dups
        
class DupConData(ConData):
    def __init__(self, con_data):
        ConData.__init__(self)
        for ref_names in con_data.con_lists.keys():
            self.con_lists[ref_names] = DupConList(con_data.con_lists[ref_names])
            
    def dup_stats(self, max_num_dups):
        hist_num_dups = [0] * max_num_dups
        for con_list in self.con_lists.values():
            list_hist_num_dups = con_list.dup_stats(max_num_dups)
            for i in range(max_num_dups):
                hist_num_dups[i] += list_hist_num_dups[i]
        return hist_num_dups
        
    def add_empty_con_list(self, ref_names):
        self.con_lists[ref_names] = DupConList()
