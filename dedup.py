import sys
import getopt
import gzip
from classes import ConData, file_to_con_data, DupConData

def dup_stats_to_string(hist_num_dups):
    return "\n".join(["#dups\t#cons"] + [(">=" if i == len(hist_num_dups) - 1 else "") + str(i + 1) + "\t" + str(hist_num_dups[i]) + " (" + str(round(100.0 * hist_num_dups[i] / sum(hist_num_dups), 2))+ "%)" for i in range(len(hist_num_dups))])        

def dedup(argv):
    # default parameters
    min_separation = 1000
    max_distance = 1000
    
    # progress display parameters
    max_num_dups = 10
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "s:d:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac dedup [options] <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -s INT     min separation (bp) for intra-chromosomal contacts [" + str(min_separation) + "]\n")
        sys.stderr.write("  -d INT     max distance (bp, L-inf norm) for merging duplicates [" + str(max_distance) + "]\n")
        return 1
    for o, a in opts:
        if o == "-s":
            min_separation = int(a)
        elif o == "-d":
            max_distance = int(a)
               
    # read CON file
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " putative contacts (" +  str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")

    # dedup and clean
    dup_con_data = DupConData(con_data)
    dup_con_data.dedup(max_distance)
    con_data.clean_separation(min_separation)
    
    # show stats and output
    sys.stderr.write("[M::" + __name__ + "] merged duplicates into " + str(dup_con_data.num_cons()) + " putative contacts (legs: " + str(round(100.0 * dup_con_data.num_phased_legs() / dup_con_data.num_cons() / 2, 2)) + "% phased, " + str(round(100.0 * dup_con_data.num_conflict_legs() / dup_con_data.num_cons() / 2, 2)) + "% conflict); statistics:\n")
    sys.stderr.write(dup_stats_to_string(dup_con_data.dup_stats(max_num_dups)) + "\n")
    sys.stdout.write(dup_con_data.to_string() + "\n")
    
    return 0
    