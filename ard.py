import sys
import getopt
import gzip
import copy
from classes import ConData, file_to_con_data, LegData, counts_to_hist_num_with_zero, hist_num_to_string_with_zero, Leg, string_to_leg, Con
import numpy as np

def add_ref_locus_to_hist(around_hist, rel_locus, max_distance, grid_size):
    leg_1_bin = (rel_locus[0] + max_distance) / grid_size
    leg_2_bin = (rel_locus[1] + max_distance) / grid_size
    grid_num = 2 * max_distance / grid_size
    if leg_1_bin < grid_num and leg_1_bin >= 0 and leg_2_bin < grid_num and leg_2_bin >= 0:
        around_hist[leg_1_bin][leg_2_bin] += 1

def ard(argv):
    # default parameters
    reference_file_name = None
    min_separation = None
    max_distance = 10000000
    grid_size = None
    is_symmetrical = True
    superellipse_mode = False
    count_mode = False
    normalize_by_num_cons = False
    leg_file_1_name = None
    leg_file_2_name = None
        
    # progress display parameters
    display_num_ref_cons = 1000
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "c:s:d:h:Sent1:2:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c ard [options] <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -c <ref.con>    contact file for reference points [<in.con> itself]\n")
        sys.stderr.write("  -s INT          only use intra-chromosomal reference points, min separation (bp) [only use inter-chromosomal] \n")
        sys.stderr.write("  -d INT          max distance (bp, L-inf norm) around reference points [" + str(max_distance) + "]\n")
        sys.stderr.write("  -h INT          output 2D histogram, grid size (bp) (useful for too many contacts)\n")
        sys.stderr.write("  -e              use L-1/2 norm (superellipse) instead\n")
        sys.stderr.write("  -S              does not symmetrize for \"-h\"\n\n")
        sys.stderr.write("  -n              output the number of nearby contacts for each reference point\n")
        sys.stderr.write("  -t              normalize by the total number of contacts for \"-n\"\n\n")
        sys.stderr.write("  -1 <in1.leg>    generate a pairwise count matrix between reference legs\n")
        sys.stderr.write("  -2 <in2.leg>    generate a pairwise count matrix between two sets of reference legs [<in2.leg>]\n")
        return 1
    for o, a in opts:
        if o == "-c":
            reference_file_name = a
        elif o == "-s":
            min_separation = int(a)
        elif o == "-d":
            max_distance = int(a)
        elif o == "-h":
            grid_size = int(a)
        elif o == "-S":
            is_symmetrical = False
        elif o == "-e":
            superellipse_mode = True     
        elif o == "-n":
            count_mode = True  
        elif o == "-t":
            normalize_by_num_cons = True  
        elif o == "-1":
            leg_file_1_name = a
        elif o == "-2":
            leg_file_2_name = a

    # read CON file
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    
    if leg_file_1_name is None:
        # regular mode
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

        # initialize 2D histogram
        if not grid_size is None:
            grid_num = 2 * max_distance / grid_size
            around_hist = np.zeros((grid_num, grid_num), dtype = np.int)
    
        # find relation positions
        con_data.sort_cons()
        num_ref_cons = 0
        for ref_con in ref_con_data.get_cons():
            num_ref_cons += 1
            if num_ref_cons % display_num_ref_cons == 0:
                sys.stderr.write("[M::" + __name__ + "] analyzed " + str(num_ref_cons) + " reference points\n")
            num_nearby_cons = 0
            for con in (con_data.get_cons_near(ref_con, max_distance) if superellipse_mode else con_data.get_cons_near_inf(ref_con, max_distance)):
                num_nearby_cons += 1
                if count_mode:
                    continue
                if grid_size is None:
                    # output relative positions
                    sys.stdout.write(con.to_string_around(ref_con) + "\n")
                else:
                    # calculate histogram
                    rel_locus = con.to_rel_locus_around(ref_con)
                    if is_symmetrical:
                        # symmetrize
                        if min_separation is None:
                            # inter-chromosomal: 8 copies
                            for sign_1 in [-1, 1]:
                                for sign_2 in [-1, 1]:
                                    add_ref_locus_to_hist(around_hist, (sign_1 * rel_locus[0], sign_2 * rel_locus[1]), max_distance, grid_size)
                                    add_ref_locus_to_hist(around_hist, (sign_2 * rel_locus[1], sign_1 * rel_locus[0]), max_distance, grid_size)
                        else:
                            # intra-chromosomal: 2 copies
                            add_ref_locus_to_hist(around_hist, (rel_locus[0], rel_locus[1]), max_distance, grid_size)
                            add_ref_locus_to_hist(around_hist, (-1 * rel_locus[1], -1 * rel_locus[0]), max_distance, grid_size)
                    else:
                        add_ref_locus_to_hist(around_hist, (rel_locus[0], rel_locus[1]), max_distance, grid_size)
            if count_mode:
                if normalize_by_num_cons:
                    sys.stdout.write(str(float(num_nearby_cons) / con_data.num_cons()) + "\n")
                else:
                    sys.stdout.write(str(num_nearby_cons) + "\n")
        # output 2D histogram
        if not grid_size is None:
            sys.stderr.write("[M::" + __name__ + "] writing output for 2D histogram\n")
            np.savetxt(sys.stdout, around_hist, delimiter='\t')
    else:
        # pairwise leg mode
        # read legs
        legs_1 = [string_to_leg(leg_file_line.strip()) for leg_file_line in open(leg_file_1_name, "rb")]
        if leg_file_2_name is None:
            legs_2 = legs_1
        else:
            legs_2 = [string_to_leg(leg_file_line.strip()) for leg_file_line in open(leg_file_2_name, "rb")]
    
        # initilize pariwise count matrix
        num_legs_1 = len(legs_1)
        num_legs_2 = len(legs_2)
        count_matrix = np.empty([num_legs_1, num_legs_2], dtype=int)
        count_matrix[:] = -1
    
        # for each pair of legs
        num_ref_cons = 0
        for i in range(num_legs_1):
            for j in (range(i + 1, num_legs_2) if leg_file_2_name is None else range(num_legs_2)):
                ref_con = Con(legs_1[i], legs_2[j])
                if min_separation is None:
                    # inter-chromosomal only
                    if ref_con.is_intra_chr():
                        continue
                else:
                    # intra-chromosmal only, remove small separations
                    if not ref_con.is_intra_chr() or ref_con.separation() < min_separation:
                        continue
                num_ref_cons += 1
                if num_ref_cons % display_num_ref_cons == 0:
                    sys.stderr.write("[M::" + __name__ + "] analyzed " + str(num_ref_cons) + " reference points\n")
                
                # count
                num_nearby_cons = 0
                for con in (con_data.get_cons_near(ref_con, max_distance) if superellipse_mode else con_data.get_cons_near_inf(ref_con, max_distance)):
                    num_nearby_cons += 1
                count_matrix[i, j] = num_nearby_cons
                if leg_file_2_name is None:
                    count_matrix[j, i] = num_nearby_cons
                
        # write pariwise count matrix
        sys.stderr.write("[M::" + __name__ + "] writing output for pairwise count matrix\n")
        np.savetxt(sys.stdout, count_matrix, fmt='%i', delimiter='\t')         
    
    return 0
    
