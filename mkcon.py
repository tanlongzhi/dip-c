import sys
import getopt
import gzip
import copy
from classes import ConData, hom_name_to_ref_name_haplotype, Leg, LegList, Con

def phased_leg_list_from_file(input_file):
    leg_list = LegList()
    for input_file_line in input_file:
        hom_name, ref_locus = input_file_line.strip().split("\t")
        ref_locus = int(ref_locus)
        ref_name, haplotype = hom_name_to_ref_name_haplotype(hom_name)
        leg_list.add_leg(Leg(ref_name, ref_locus, haplotype))
    return leg_list

def mkcon(argv):
    # default parameters
    num_con = 500
    
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
        sys.stderr.write("Usage: dip-c mkcon [options] <leg1.txt> <leg2.txt>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -n INT     number of contacts to output [" + str(num_con) + "]\n")
        return 1
    elif len(args) != 2:
        sys.stderr.write("[E::" + __name__ + "] exactly two leg files are needed\n")
        return 1
    for o, a in opts:
        if o == "-n":
            num_con = int(a)
    
                            
    # read input file
    leg_list_1 = phased_leg_list_from_file(open(args[0], "rb"))
    sys.stderr.write("[M::" + __name__ + "] leg 1: read " + str(leg_list_1.num_legs()) + " legs\n")
    leg_list_2 = phased_leg_list_from_file(open(args[1], "rb"))
    sys.stderr.write("[M::" + __name__ + "] leg 2: read " + str(leg_list_2.num_legs()) + " legs\n")
    
    # make CON
    con_data = ConData()
    for i in range(num_con):
        con_data.add_con(Con(leg_list_1.get_random_leg(), leg_list_2.get_random_leg()))
    con_data.sort_cons()
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    sys.stdout.write(con_data.to_string()+"\n")
    
    
    return 0
    
