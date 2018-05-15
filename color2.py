import sys
import getopt
import gzip
from classes import Haplotypes, hom_name_to_ref_name_haplotype, ref_name_haplotype_to_hom_name, ConData, file_to_con_data, Con, Leg, G3dData, file_to_g3d_data
import numpy as np
from scipy import spatial
import math
    
def leg_to_bin_ref_locus(leg, bin_size):
    return int(round(float(leg.get_ref_locus()) / bin_size)) * bin_size

def color2(argv):
    # default parameters
    chr_len_file_name = None
    color_file_name = None
    bin_size = 1000000
    max_step = 1000
    force_constant = 1e-2
    friction = 1e-1
    disp_num_steps = 100
    disp_num_cons = 100000
    merge_haplotypes = False
    smooth_mode = False
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "b:Hc:s")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c color2 [options] <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -b INT            bin size (bp) (bins are centered around multiples of bin size) [" + str(bin_size) + "]\n")
        sys.stderr.write("  -H                merge the two haplotypes\n")
        sys.stderr.write("  -c <color.txt>    color by a list of locus-color pairs (tab-delimited: chr, locus, color)\n")
        sys.stderr.write("  -s                smooth color by averaging over all contacted bins\n")
        sys.stderr.write("Output:\n")
        sys.stderr.write("  tab-delimited: homolog or chr if \"-H\", locus, color\n")
        return 1
        
    num_color_schemes = 0
    for o, a in opts:
        if o == "-o":
            output_prefix = a
        elif o == "-b":
            bin_size = int(a)
        elif o == "-H":
            merge_haplotypes = True
        elif o == "-s":
            smooth_mode = True
        elif o == "-c":
            color_file_name = a

    # open color file
    if not color_file_name is None:
        color_file = open(color_file_name, "rb")
                                                        
    # read CON file
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")

    # calculate colors for each particle
    color_data = {}
    for color_file_line in color_file:
        hom_name, ref_locus, color = color_file_line.strip().split("\t")
        ref_locus = int(ref_locus)
        color = float(color)
        color_data[(hom_name, ref_locus)] = color
    
    # smoothing
    if not smooth_mode is None:
        smooth_color_raw_data = {}
        for con in con_data.get_cons():
            bin_ref_locus_1 = leg_to_bin_ref_locus(con.leg_1(), bin_size)
            bin_ref_locus_2 = leg_to_bin_ref_locus(con.leg_2(), bin_size)
            if merge_haplotypes:
                hom_name_1 = con.leg_1().get_ref_name()
                hom_name_2 = con.leg_2().get_ref_name()
            else:
                hom_name_1 = ref_name_haplotype_to_hom_name(con.leg_1().ref_name_haplotype())
                hom_name_2 = ref_name_haplotype_to_hom_name(con.leg_2().ref_name_haplotype())
            hom_name_bin_ref_locus_1 = (hom_name_1, bin_ref_locus_1)
            hom_name_bin_ref_locus_2 = (hom_name_2, bin_ref_locus_2)
            
            #sys.stderr.write(str(hom_name_bin_ref_locus_1)+" - " + str(hom_name_bin_ref_locus_2) + "\n")
            
            
            # add raw data
            if hom_name_bin_ref_locus_1 == hom_name_bin_ref_locus_2:
                continue
            if hom_name_bin_ref_locus_1 in color_data:
                if hom_name_bin_ref_locus_2 not in smooth_color_raw_data:
                    smooth_color_raw_data[hom_name_bin_ref_locus_2] = []
                smooth_color_raw_data[hom_name_bin_ref_locus_2].append(color_data[hom_name_bin_ref_locus_1])
            if hom_name_bin_ref_locus_2 in color_data:
                if hom_name_bin_ref_locus_1 not in smooth_color_raw_data:
                    smooth_color_raw_data[hom_name_bin_ref_locus_1] = []
                smooth_color_raw_data[hom_name_bin_ref_locus_1].append(color_data[hom_name_bin_ref_locus_2])
        
        # averaging
        for hom_name_bin_ref_locus in smooth_color_raw_data:
            smooth_color_raw_data[hom_name_bin_ref_locus] = np.mean(smooth_color_raw_data[hom_name_bin_ref_locus])
            
        color_data = smooth_color_raw_data
    
    # output
    sys.stderr.write("[M::" + __name__ + "] writing " + str(len(color_data)) + " colors\n")
    for hom_name, ref_locus in sorted(color_data.keys()):
        sys.stdout.write("\t".join([hom_name, str(ref_locus), str(color_data[(hom_name, ref_locus)])]) + "\n")
    
    return 0
    
