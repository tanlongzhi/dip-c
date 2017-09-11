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
    is_read2 = False
    query_start = -1
    query_end = -1
    ref_name = ""
    ref_start = -1
    ref_end = -1
    is_reverse = False
    haplotype = Haplotypes.unknown
    
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
                
    # merge a read 1 with its adjacent read 2 (now query_end becomes undefined, set to -1), return True if worked
    def merge_mate(self, other, max_merge_distance):
        # only merge read 1 with read 2
        if self.is_read2 or not other.is_read2:
            return False
        # skip bad pairs
        if self.ref_name != other.ref_name or self.is_reverse == other.is_reverse or abs(self.ref_right() - other.ref_left()) > max_merge_distance:
            return False
        self.query_end = -1
        self.set_ref_right(other.ref_right())
        self.haplotype = update_haplotype(self.haplotype, other.haplotype) 
        return True
    
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
    name = ""
    segs = []
    
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
        self.segs = sorted(self.segs)
            
    # try to merge its segs (must be sorted)
    def merge_mate(self, max_merge_distance):
        for i in range(len(self.segs) - 1):
            if self.segs[i].merge_mate(self.segs[i + 1], max_merge_distance):
                self.segs.pop(i + 1)
                return True
                break
        return False
        
    # keep merging until not possible (must be sorted)
    def merge_mates(self, max_merge_distance):
        while True:
            if not self.merge_mate(max_merge_distance):
                break
    
    def to_string(self):
        return self.name + "\t" + "\t".join([seg.to_string() for seg in self.segs])
        
# create a read from a string
def string_to_read(read_string):
    read_string_data = read_string.split("\t")
    read = Read(read_string_data[0])
    for seg_string in read_string_data[1:]:
        read.add_seg(string_to_seg(seg_string))
    return read

# a hash map of reads (a .seg file)
class SegData:
    reads = {}
    
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

# a contact
#class Con:
