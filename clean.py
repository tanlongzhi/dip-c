import sys
import getopt
import gzip
import copy
from classes import ConData, file_to_con_data, LegData, counts_to_hist_num_with_zero, hist_num_to_string_with_zero, Leg

def clean(argv):
    # default parameters
    max_clean_distance = 10000000
    min_clean_count = 5
    max_leg_distance = 1000
    max_leg_count = 10
    test_mode = False
    
    # progress display parameters
    display_max_num_legs = 20
    display_num_cons = 10000
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "d:c:D:C:t")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c clean [options] <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -d INT     max distance (bp, L-1/2 norm) for removing isolated contacts [" + str(max_clean_distance) + "]\n")
        sys.stderr.write("  -c INT     min neighbor count for an unisolated contact [" + str(min_clean_count) + "]\n")
        sys.stderr.write("  -D INT     max distance (bp) for removing promiscuous legs [" + str(max_leg_distance) + "]\n")
        sys.stderr.write("  -C INT     max neighbor count for a nonpromiscuous leg [" + str(max_leg_count) + "]\n")
        sys.stderr.write("  -t         test mode for isolated contacts: statistics only, no removal\n")
        sys.stderr.write("               (a larger -c, for example 20, is recommended for this mode)\n")
        return 1
    for o, a in opts:
        if o == "-d":
            max_clean_distance = int(a)
        elif o == "-c":
            min_clean_count = int(a)
        elif o == "-D":
            max_leg_distance = int(a)
        elif o == "-C":
            max_leg_count = int(a)
        elif o == "-t":
            test_mode = True     
                            
    # read CON file
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    original_num_cons = con_data.num_cons()
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " putative contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    #sys.stdout.write(con_data.to_string() + "\n")
    
    # pass 1: remove contacts containing promiscuous legs
    leg_data = LegData()
    leg_data.add_con_data(con_data)
    leg_data.sort_legs()
    sys.stderr.write("[M::" + __name__ + "] pass 1: sorted " + str(leg_data.num_legs()) + " legs\n")
    con_data.clean_promiscuous(leg_data, max_leg_distance, max_leg_count)
    pass_1_num_cons = con_data.num_cons()
    sys.stderr.write("[M::" + __name__ + "] pass 1 done: removed " + str(original_num_cons - pass_1_num_cons) + " contacts (" + str(round(100.0 - 100.0 * pass_1_num_cons / original_num_cons, 2)) + "%)\n")
    
    # pass 2: remove isolated contacts
    con_data.sort_cons()
    sys.stderr.write("[M::" + __name__ + "] pass 2: sorted " + str(pass_1_num_cons) + " contacts\n")
    if test_mode:
        neighbor_counts = con_data.test_isolated(copy.deepcopy(con_data), max_clean_distance, min_clean_count)
        sys.stderr.write("[M::" + __name__ + "] pass 2 test done. statistics:\n")
        sys.stderr.write("#neighbors\t#cons\n")
        sys.stderr.write(hist_num_to_string_with_zero(counts_to_hist_num_with_zero(neighbor_counts)) + "\n")
    else:
        con_data.clean_isolated(copy.deepcopy(con_data), max_clean_distance, min_clean_count)
        pass_2_num_cons = con_data.num_cons()
        sys.stderr.write("[M::" + __name__ + "] pass 2 done: removed " + str(pass_1_num_cons - pass_2_num_cons) + " contacts (" + str(round(100.0 * (pass_1_num_cons - pass_2_num_cons) / original_num_cons, 2)) + "%)\n")
        sys.stderr.write("[M::" + __name__ + "] writing output for " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
        sys.stdout.write(con_data.to_string()+"\n")
        
    
    return 0
    
