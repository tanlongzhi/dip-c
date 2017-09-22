import sys
import getopt
import gzip
import copy
from classes import ConData, file_to_con_data, LegData, counts_to_hist_num_with_zero, hist_num_to_string_with_zero, Leg

def ard(argv):
    # default parameters
    reference_file_name = None
    min_separation = None
    max_distance = 10000000
    
    # progress display parameters
    display_num_ref_cons = 10000
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "c:s:d:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac ard [options] <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -c <ref.con>    contact file for reference points [<in.con> itself]\n")
        sys.stderr.write("  -s INT          only use intra-chromosomal reference points, min separation (bp) [only use inter-chromosomal] \n")
        sys.stderr.write("  -d INT          max distance (bp, L-inf norm) around reference points [" + str(max_distance) + "]\n")
        return 1
    for o, a in opts:
        if o == "-c":
            reference_file_name = a
        elif o == "-s":
            min_separation = int(a)
        elif o == "-d":
            max_distance = int(a)
                            
    # read CON file
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    
    # read reference CON file
    if reference_file_name is None:
        # use itself
        ref_con_data = copy.deepcopy(con_data)
    else:
        # open another file
        ref_con_file = gzip.open(reference_file_name, "rb") if reference_file_name.endswith(".gz") else open(reference_file_name, "rb")
        ref_con_data = file_to_con_data(ref_con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(ref_con_data.num_cons()) + " reference points (" + str(round(100.0 * ref_con_data.num_intra_chr() / ref_con_data.num_cons(), 2)) + "% intra-chromosomal)\n")
    
    # keep only desired reference points
    if min_separation is None:
        # inter-chromosomal only
        ref_con_data.clean_intra_chr()
    else:
        # intra-chromosmal only, remove small separations
        ref_con_data.clean_inter_chr()
        ref_con_data.clean_separation(min_separation)
    sys.stderr.write("[M::" + __name__ + "] kept " + str(ref_con_data.num_cons()) + " reference points (" + str(round(100.0 * ref_con_data.num_intra_chr() / ref_con_data.num_cons(), 2)) + "% intra-chromosomal)\n")
    
    # find relation positions
    con_data.sort_cons()
    num_ref_cons = 0
    for ref_con in ref_con_data.get_cons():
        num_ref_cons += 1
        if num_ref_cons % display_num_ref_cons == 0:
            sys.stderr.write("[M::" + __name__ + "] analyzed " + str(num_ref_cons) + " reference points\n")
        for con in con_data.get_cons_near_inf(ref_con, max_distance):
            sys.stdout.write(con.to_string_around(ref_con) + "\n")
            
    
    return 0
    