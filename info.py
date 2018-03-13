import sys
import getopt
import gzip
import copy
from classes import ConData, file_to_con_data, LegData, counts_to_hist_num_with_zero, hist_num_to_string_with_zero, Leg

def info(argv):
    # default parameters
    reference_file_name = None
    min_separation = None
    max_distance = 10000000
    grid_size = None
    
    # progress display parameters
    display_num_ref_cons = 1000
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c info <in1.con> <in2.con> ...\n")
        return 1
                                        
    # read each CON file
    for a in args:
        con_file = gzip.open(a, "rb") if a.endswith(".gz") else open(a, "rb")
        con_data = file_to_con_data(con_file)
        con_file.close()
        sys.stdout.write(a + ": " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    
    return 0
    
