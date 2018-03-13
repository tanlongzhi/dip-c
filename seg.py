import sys
import getopt
import pysam
import gzip
from classes import Haplotypes, Seg, Read, SegData

# add the primary alignment as a segment to a read, if satisfying conditions
def add_primary_seg(read, bam_read, min_mapq, max_nm_per_bp):
    if bam_read.mapping_quality < min_mapq or float(bam_read.get_tag("NM")) / (bam_read.query_alignment_end - bam_read.query_alignment_start) > max_nm_per_bp:
        return
    read.add_seg(Seg(bam_read.is_read2, bam_read.query_alignment_start, bam_read.query_alignment_end, bam_read.reference_name, bam_read.reference_start, bam_read.reference_end, bam_read.is_reverse))

# add all SA alignments as segments to a read, if satisfying conditions
def add_sa_segs(read, bam_read, min_mapq, max_nm_per_bp):
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
    min_baseq = 20
    snp_file_name = None
    
    # progress display parameters
    display_num_bam_reads = 1e5
    display_num_snps = 1e4
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "q:m:v:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c seg [options] <in.bam>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -q INT          min mapping quality [" + str(min_mapq) + "]\n")
        sys.stderr.write("  -m FLOAT        max edit distance per bp of aligned read [" + str(max_nm_per_bp) + "]\n")
        sys.stderr.write("  -v <snps.txt>   SNP file:\n")
        sys.stderr.write("                    tab-delimited: chr, pos, paternal nucleotide, maternal nucleotide\n")
        sys.stderr.write("  -Q              min base quality (if -v) [" + str(min_baseq) + "]\n")
        return 1
    for o, a in opts:
        if o == "-q":
            min_mapq = int(a)
        elif o == "-m":
            max_nm_per_bp = float(a)
        elif o == "-v":
            snp_file_name = a
            
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
            add_sa_segs(read, bam_read, min_mapq, max_nm_per_bp)
        elif bam_read.is_paired and not bam_read.is_proper_pair:
            read = Read(bam_read.query_name)
            add_primary_seg(read, bam_read, min_mapq, max_nm_per_bp)
        else:
            continue
        seg_data.add_read(read)
    sys.stderr.write("[M::" + __name__ + "] pass 1 done: read " + str(num_bam_reads) + " alignments; added " + str(seg_data.num_reads()) + " candidate reads (" + str(round(100.0 * seg_data.num_reads() / num_bam_reads, 2)) + "% of alignments)\n")

    # pass 2: find mates of all reads from pass 1, and then removes reads with only one segment
    num_bam_reads = 0
    bam_file.reset() # not a documented API, but works
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
    sys.stderr.write("[M::" + __name__ + "] pass 2: cleaning " + str(seg_data.num_reads()) + " candidate reads\n")
    seg_data.clean()
    sys.stderr.write("[M::" + __name__ + "] pass 2 done: read " + str(num_bam_reads) + " alignments; kept " + str(seg_data.num_reads()) + " candidate reads (" + str(round(100.0 * seg_data.num_reads() / num_bam_reads, 2)) + "% of alignments)\n")
    
    # pass 3: (if -v) find haplotype information based on the SNP file
    if not snp_file_name is None:
        num_snps = 0
        snp_file = gzip.open(snp_file_name, "rb") if snp_file_name.endswith(".gz") else open(snp_file_name, "rb")
        for snp_file_line in snp_file:
            num_snps += 1
            snp_chr, snp_locus, snp_pat, snp_mat = snp_file_line.strip().split()
            snp_locus = int(snp_locus)
            if num_snps % display_num_snps == 0:
                sys.stderr.write("[M::" + __name__ + "] pass 3: read " + str(num_snps) + " SNPs, last at " + snp_chr + ":" + str(snp_locus) + "\n")
            
            for pileup_column in bam_file.pileup(snp_chr, snp_locus - 1, snp_locus):
                if pileup_column.pos == snp_locus - 1:
                    for pileup_read in pileup_column.pileups:
                        if not pileup_read.is_del and not pileup_read.is_refskip:
                            # skip read if not a candidate
                            if not seg_data.contains_read_name(pileup_read.alignment.query_name):
                                continue
                            # discard if low base quality
                            if pileup_read.alignment.query_qualities[pileup_read.query_position] < min_baseq:
                                continue
                            
                            # phase the corresponding segment
                            seg_haplotype = -1
                            if pileup_read.alignment.query_sequence[pileup_read.query_position] == snp_pat:
                                seg_haplotype = 0
                            elif pileup_read.alignment.query_sequence[pileup_read.query_position] == snp_mat:
                                seg_haplotype = 1
                            else:
                                continue
                            seg_data.update_haplotype(pileup_read.alignment.query_name, pileup_read.alignment.is_read2, snp_chr, snp_locus, seg_haplotype)
        sys.stderr.write("[M::" + __name__ + "] pass 3 done: read " + str(num_snps) + " SNPs\n")
    
    # write output
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(seg_data.num_reads()) + " candidate reads, containing " + str(seg_data.num_segs()) + " segments (" + str(round(100.0 * seg_data.num_phased_segs() / seg_data.num_segs(), 2)) + "% phased)\n")
    sys.stdout.write(seg_data.to_string()+"\n")
  
    return 0
    
