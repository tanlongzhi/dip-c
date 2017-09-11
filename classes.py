# enum for haplotypes
class Haplotypes:
    conflict = -2
    unknown = -1
    paternal = 0
    maternal = 1
def is_known_haplotype(haplotype):
    return haplotype >= 0

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
    
    def to_string(self): # "m" is for mate, "." for all unknown haplotypes
        return ",".join(["m" if self.is_read2 else ".", str(self.query_start), str(self.query_end), self.ref_name, str(self.ref_start), str(self.ref_end), "-" if self.is_reverse else "+", str(self.haplotype) if is_known_haplotype(self.haplotype) else "."])
    
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
    
    def to_string(self):
        return self.name + "\t" + "\t".join([seg.to_string() for seg in self.segs])
        
# a hash map of reads (a .seg file)
class SegData:
    reads = {}
    
    def __init__(self):
        self.reads = {}

    def contains_read_name(self, name):
        return name in self.reads

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
            
    def num_reads(self):
        return len(self.reads)
            
    def to_string(self): # no tailing new line
        return "\n".join(read.to_string() for read in self.reads.values())
            