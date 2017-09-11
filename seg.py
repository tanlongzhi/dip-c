import sys
import getopt
import pysam
from classes import Haplotypes, is_known_haplotype, Seg, Read, SegData

def add_primary_seg(read, bam_read, min_mapq, max_nm_per_bp):
    if bam_read.mapping_quality < min_mapq or float(bam_read.get_tag("NM")) / (bam_read.query_alignment_end - bam_read.query_alignment_start) > max_nm_per_bp:
        return
    read.add_seg(Seg(bam_read.is_read2, bam_read.query_alignment_start, bam_read.query_alignment_end, bam_read.reference_name, bam_read.reference_start, bam_read.reference_end, bam_read.is_reverse))

def add_sa_segs(read, bam_read, bam_file, min_mapq, max_nm_per_bp):
    sa_tags = bam_read.get_tag("SA")
    # read each part of the SA tag
    for sa_tag in sa_tags.split(";")[:-1]:
        ref_name, ref_start, strand, cigar, mapping_quality, edit_distance = sa_tag.split(",")
        mapping_quality = int(mapping_quality)
        edit_distance = int(edit_distance)
        ref_start = int(ref_start)
        if mapping_quality < min_mapq:
            continue
        # hack pysam's current alignment (primary) to get SA alignment 
        bam_read.reference_id = bam_file.get_tid(ref_name)
        bam_read.reference_start = int(ref_start)
        bam_read.cigarstring = cigar
        if bam_read.is_reverse == (strand == "-"):
            query_start = bam_read.query_alignment_start
            query_end = bam_read.query_alignment_end
        else:
            query_start = bam_read.query_length - bam_read.query_alignment_end
            query_end = bam_read.query_length - bam_read.query_alignment_start
        if float(edit_distance) / (query_end - query_start) > max_nm_per_bp:
            continue
        read.add_seg(Seg(bam_read.is_read2, query_start, query_end, ref_name, ref_start, bam_read.reference_end, True if strand == "-" else False))

def seg(argv):
    # default parameters
    min_mapq = 20
    max_nm_per_bp = 0.05
    
    # progress display parameters
    display_num_bam_reads = 1e5
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "q:m:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac seg [options] <in.bam>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -q INT     min mapping quality [" + str(min_mapq) + "]\n")
        sys.stderr.write("  -m FLOAT   max edit distance per bp of aligned read [" + str(max_nm_per_bp) + "]\n")
        return 1
    for o, a in opts:
        if o == "-q":
            min_mapq = int(a)
        elif o == "-m":
            max_nm = int(a)
            
    # initialize data structure
    seg_data = SegData()
    
    # pass 1: find reads with SA tags and/or discordant reads
    bam_file = pysam.AlignmentFile(args[0], "rb")
    num_bam_reads = 0
    for bam_read in bam_file.fetch(until_eof=True):
        num_bam_reads += 1
        if num_bam_reads % display_num_bam_reads == 0:
            sys.stderr.write("[M::" + __name__ + "] pass 1: read " + str(num_bam_reads) + " alignments, last at " + (bam_read.reference_name + ":" + str(bam_read.reference_start) if not bam_read.is_unmapped else "*") + "; added " + str(seg_data.num_reads()) + " candidate reads\n")

        # exclude bad reads
        if bam_read.is_unmapped or bam_read.is_supplementary or bam_read.is_duplicate or bam_read.is_qcfail or bam_read.is_secondary:
            continue
        
        # get segments from reads
        if bam_read.has_tag("SA"):
            read = Read(bam_read.query_name)
            add_primary_seg(read, bam_read, min_mapq, max_nm_per_bp)
            add_sa_segs(read, bam_read, bam_file, min_mapq, max_nm_per_bp)
        elif bam_read.is_paired and not bam_read.is_proper_pair:
            read = Read(bam_read.query_name)
            add_primary_seg(read, bam_read, min_mapq, max_nm_per_bp)
        else:
            continue
        seg_data.add_read(read)

    # pass 2: find mates of all reads from pass 1
    num_bam_reads = 0
    bam_file.reset()
    for bam_read in bam_file.fetch(until_eof=True):
        num_bam_reads += 1
        if num_bam_reads % display_num_bam_reads == 0:
            sys.stderr.write("[M::" + __name__ + "] pass 2: read " + str(num_bam_reads) + " alignments, last at " + (bam_read.reference_name + ":" + str(bam_read.reference_start) if not bam_read.is_unmapped else "*") + "\n")
        if bam_read.has_tag("SA") or bam_read.is_paired and not bam_read.is_proper_pair:
            continue # would have been added in the first round
        if seg_data.contains_read_name(bam_read.query_name):
            read = Read(bam_read.query_name)
            add_primary_seg(read, bam_read, min_mapq, max_nm_per_bp)
        else:
            continue
        seg_data.add_read(read)
            
    sys.stderr.write("[M::" + __name__ + "] cleaning " + str(seg_data.num_reads()) + " candidate reads\n")
    seg_data.clean()
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(seg_data.num_reads()) + " candidate reads \n")
    sys.stdout.write(seg_data.to_string()+"\n")
    
    # pass 3: find haplotype information
        
    return 0

            
