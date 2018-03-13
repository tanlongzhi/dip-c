import sys
import getopt
import gzip
from classes import Read, string_to_read, ConData

def con(argv):
    # default parameters
    min_separation = 1000
    max_distance = 1000
    adjacent_only = False
    
    # progress display parameters
    display_num_reads = 1e4
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "s:d:a")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c con [options] <in.seg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -s INT     min separation (bp) for intra-chromosomal contacts [" + str(min_separation) + "]\n")
        sys.stderr.write("  -d INT     max distance (bp, L-inf norm) for merging contacts within a read [" + str(max_distance) + "]\n")
        sys.stderr.write("  -a         contact between adjacent segments only\n")
        return 1
    for o, a in opts:
        if o == "-s":
            min_separation = int(a)
        elif o == "-d":
            max_distance = int(a)
        elif o == "-a":
            adjacent_only = True
               
    # read SEG file
    seg_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    num_reads = 0   
    con_data = ConData()
    for seg_file_line in seg_file:
        num_reads += 1
        if num_reads % display_num_reads == 0:
            sys.stderr.write("[M::" + __name__ + "] read " + str(num_reads) + " candidate reads; added " + str(con_data.num_cons()) + " putative contacts\n")
        read = string_to_read(seg_file_line.strip())
        read_con_data = read.to_con_data(adjacent_only)
        
        # dedup and clean within each read
        read_con_data.dedup_within_read(max_distance)
        read_con_data.clean_separation(min_separation)
        
        # add to full data
        con_data.merge_with(read_con_data)
    
    # sort and output
    sys.stderr.write("[M::" + __name__ + "] read " + str(num_reads) + " candidate reads; sorting " + str(con_data.num_cons()) + " putative contacts\n")
    con_data.sort_cons()
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(con_data.num_cons()) + " putative contacts (" +  str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    sys.stdout.write(con_data.to_string()+"\n")
    
    return 0
    
