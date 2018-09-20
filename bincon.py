import sys
import gzip
import getopt
from classes import Haplotypes, string_to_leg, ref_name_haplotype_to_hom_name, ConData, file_to_con_data, Con, Leg
import numpy as np
from scipy import spatial
import math
    
def leg_to_matrix_index(leg, hom_offsets, matrix_bin_size, merge_haplotypes):
    ref_name, haplotype = leg.ref_name_haplotype()
    if merge_haplotypes:
        # merge both haplotypes into "paternal"
        haplotype = Haplotypes.paternal
    hom_name = ref_name_haplotype_to_hom_name((ref_name, haplotype))
    ref_locus = leg.get_ref_locus()
    if hom_name not in hom_offsets:
        return None
    return hom_offsets[hom_name] + int(round(float(ref_locus) / matrix_bin_size))

def con_to_matrix_index(con, hom_offsets, matrix_bin_size, merge_haplotypes):
    matrix_index_1 = leg_to_matrix_index(con.leg_1(), hom_offsets, matrix_bin_size, merge_haplotypes)
    matrix_index_2 = leg_to_matrix_index(con.leg_2(), hom_offsets, matrix_bin_size, merge_haplotypes)
    return matrix_index_1, matrix_index_2
    
def con_data_to_matrix(con_data, hom_offsets, matrix_bin_size, matrix_size, merge_haplotypes, display_num_cons):
    # convert to matrix
    sys.stderr.write("[M::" + __name__ + "] generating a " + str(matrix_size) + " x " + str(matrix_size) + " matrix\n")
    matrix_data = np.zeros((matrix_size, matrix_size), dtype = int)
    num_cons = 0
    for con in con_data.get_cons():
        num_cons += 1
        if num_cons % display_num_cons == 0:
            sys.stderr.write("[M::" + __name__ + "] processed " + str(num_cons) + " contacts\n")
        matrix_index_1, matrix_index_2 = con_to_matrix_index(con, hom_offsets, matrix_bin_size, merge_haplotypes)
        if matrix_index_1 is None or matrix_index_2 is None:
            continue
        matrix_data[matrix_index_1, matrix_index_2] += 1
        matrix_data[matrix_index_2, matrix_index_1] += 1
    return matrix_data

def bincon(argv):
    # default parameters
    chr_len_file_name = None
    matrix_bin_size = 1000000
    merge_haplotypes = False
    info_mode = False
    leg_mode = False
    min_separation = 0
    
    # progress display parameters
    display_num_cons = 1e4
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "l:b:HiLs:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c bincon [options] -l <chr.len> <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -l <chr.len>   file containing chromosome lengths (tab-delimited: chr, len)\n")
        sys.stderr.write("  -L             analyze LEG instead of CON\n")
        sys.stderr.write("  -b INT         bin size (bp) (bins are centered around multiples of bin size) [" + str(matrix_bin_size) + "]\n")
        sys.stderr.write("  -H             merge the two haplotypes\n")
        sys.stderr.write("  -s INT         min separation (bp) for intra-chromosomal contacts [" + str(min_separation) + "]\n")
        sys.stderr.write("  -i             output bin info (tab-delimited: homolog or chr if \"-H\", bin center) instead\n")
        return 1
        
    num_color_schemes = 0
    for o, a in opts:
        if o == "-l":
            matrix_mode = True
            chr_len_file_name = a
        elif o == "-s":
            min_separation = int(a)
        elif o == "-b":
            matrix_bin_size = int(a)
        elif o == "-H":
            merge_haplotypes = True
        elif o == "-i":
            info_mode = True
        elif o == "-L":
            leg_mode = True
    if chr_len_file_name is None:
        sys.stderr.write("[E::" + __name__ + "] -l is required\n")
        return 1

    # read chromosome lengths
    hom_lens = {}
    hom_bin_lens = {}
    hom_offsets = {}
    matrix_size = 0
    chr_len_file = open(chr_len_file_name, "rb")
    for chr_len_file_line in chr_len_file:
        ref_name, ref_len = chr_len_file_line.strip().split("\t")
        ref_len = int(ref_len)
        for haplotype in ([Haplotypes.paternal] if merge_haplotypes else [Haplotypes.paternal, Haplotypes.maternal]):
            hom_name = ref_name_haplotype_to_hom_name((ref_name, haplotype))
            hom_bin_len = int(round(float(ref_len) / matrix_bin_size)) + 1
            hom_lens[hom_name] = ref_len
            hom_bin_lens[hom_name] = hom_bin_len
            hom_offsets[hom_name] = matrix_size
            matrix_size += hom_bin_len
            
            if info_mode:
                for bin_id in range(hom_bin_len):
                    sys.stdout.write("\t".join([(ref_name if merge_haplotypes else hom_name), str(bin_id * matrix_bin_size)]) + "\n")
    
    # generate matrix
    if not info_mode:
        if leg_mode:
            matrix_data = np.zeros((matrix_size, 1), dtype = int)
            for leg_file_line in open(args[0], "rb"):
                leg = string_to_leg(leg_file_line.strip())
                matrix_data[leg_to_matrix_index(leg, hom_offsets, matrix_bin_size, merge_haplotypes)] += 1
        else:
            con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
            con_data = file_to_con_data(con_file)
            con_data.clean_separation(min_separation)
            sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " putative contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
            matrix_data = con_data_to_matrix(con_data, hom_offsets, matrix_bin_size, matrix_size, merge_haplotypes, display_num_cons)
        np.savetxt(sys.stdout, matrix_data, fmt='%i', delimiter='\t')
    
    return 0
    
