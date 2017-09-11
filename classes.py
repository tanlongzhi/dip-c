import bisect

# enum for haplotypes
class Haplotypes:
    conflict = -2
    unknown = -1
    paternal = 0
    maternal = 1
def is_known_haplotype(haplotype):
    return haplotype >= 0

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
    
    def to_string(self): # "m" is for mate, "." for all unphased haplotypes
        return ",".join(["m" if self.is_read2 else ".", str(self.query_start), str(self.query_end), self.ref_name, str(self.ref_start), str(self.ref_end), "-" if self.is_reverse else "+", str(self.haplotype) if is_known_haplotype(self.haplotype) else "."])

# create a segment from a string ("." will be set to unknown)
def string_to_seg(seg_string):
    is_read2, query_start, query_end, ref_name, ref_start, ref_end, is_reverse, haplotype = seg_string.split(",")
    is_read2 = True if is_read2 == "m" else False
    query_start = int(query_start)
    query_end = int(query_end)
    ref_start = int(ref_start)
    ref_end = int(ref_end)
    is_reverse = True if is_reverse == "-" else False
    haplotype = Haplotypes.unknown if haplotype == "." else int(haplotype)
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
        print self.to_string()
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
        
    def __eq__(self, other):
        return self.ref_name == other.ref_name and self.ref_locus == other.ref_locus and self.haplotype == other.haplotype
    
    def get_ref_name(self):
        return self.ref_name

    def get_ref_locus(self):
        return self.ref_locus
    
    def to_string(self):
        return ",".join([self.ref_name, str(self.ref_locus), "." if not is_known_haplotype(self.haplotype) else str(self.haplotype)])

# a contact (legs sorted)
class Con:
    def __init__(self, leg_1, leg_2):
        self.legs = [leg_1, leg_2] if leg_1 < leg_2 else [leg_2, leg_1]
    
    def __lt__(self, other):
        return self.leg_1() < other.leg_1() or self.leg_1() == other.leg_1() and self.leg_2() < other.leg_2()
    
    def leg_1(self):
        return self.legs[0]

    def leg_2(self):
        return self.legs[1]
   
    def sort_legs(self):
        self.legs.sort()
    
    def ref_names(self):
        return tuple([leg.get_ref_name() for leg in self.legs])
        
    def is_intra_chr(self):
        return self.leg_1().get_ref_name() == self.leg_2().get_ref_name()
    
    def separation(self):
        return self.leg_2().get_ref_locus() - self.leg_1().get_ref_locus()
            
    def to_string(self):
        return "\t".join([leg.to_string() for leg in self.legs])
        
# a sorted list of contacts
class ConList:
    def __init__(self):
        self.cons = []
    
    def add_con(self, con):
        bisect.insort(self.cons, con)
        
    # remove intra-chromosomal contacts with small separations
    def clean_separation(self, min_separation):
        self.cons[:] = [con for con in self.cons if not con.is_intra_chr() or con.separation() > min_separation]
                
    def to_string(self):
        return "\n".join([con.to_string() for con in self.cons])
        
# a hashmap (tuples of two sorted chromosome names) of lists of contacts (a CON file)
class ConData:
    def __init__(self):
        self.con_lists = {}
    
    def add_con(self, con):
        if con.ref_names() not in self.con_lists:
            self.con_lists[con.ref_names()] = ConList()
        self.con_lists[con.ref_names()].add_con(con)
    
    def clean_separation(self, min_separation):
        for con_list in self.con_lists.values():
            con_list.clean_separation(min_separation)
    
    def to_string(self): # no tailing new line
        return "\n".join([con_list.to_string() for con_list in self.con_lists.values()])
    