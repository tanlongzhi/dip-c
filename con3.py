import sys
import getopt
from classes import Haplotypes, hom_name_to_ref_name_haplotype, ref_name_haplotype_to_hom_name, ConData, file_to_con_data, Con, Leg, G3dData, file_to_g3d_data
import numpy as np
from scipy import spatial
import math

def g3d_np_arrays_to_leg(hom_names, loci_np_array, index):
    hom_name = hom_names[index]
    ref_locus = int(loci_np_array[index])
    ref_name, haplotype = hom_name_to_ref_name_haplotype(hom_name)
    return Leg(ref_name, ref_locus, haplotype)
    
def g3d_np_arrays_to_matrix_index(hom_names, loci_np_array, index, hom_offsets, matrix_bin_size, merge_haplotypes):
    hom_name = hom_names[index]
    if merge_haplotypes:
        # merge both haplotypes into "paternal"
        ref_name, haplotype = hom_name_to_ref_name_haplotype(hom_name)
        hom_name = ref_name_haplotype_to_hom_name((ref_name, Haplotypes.paternal))
    ref_locus = int(loci_np_array[index])
    return hom_offsets[hom_name] + int(round(float(ref_locus) / matrix_bin_size))
    
def g3d_data_to_con_data(g3d_data, max_distance):
    hom_names, loci_np_array, position_np_array = g3d_data.to_np_arrays()
    sys.stderr.write("[M::" + __name__ + "] preparing genome-wide KD tree\n")
    kdtree = spatial.KDTree(position_np_array)
    sys.stderr.write("[M::" + __name__ + "] finding particle pairs\n")
    index_pairs = kdtree.query_pairs(max_distance)
    sys.stderr.write("[M::" + __name__ + "] found " + str(len(index_pairs)) + " pairs of particles (" + str(round(float(len(index_pairs)) / g3d_data.num_g3d_particles(), 2)) + " pairs per particle)\n")
    
    # convert to contacts
    con_data = ConData()
    for index_pair in index_pairs:
        index_1, index_2 = index_pair
        con_data.add_con(Con(g3d_np_arrays_to_leg(hom_names, loci_np_array, index_1), g3d_np_arrays_to_leg(hom_names, loci_np_array, index_2)))
    return con_data
    
def g3d_data_to_matrix(g3d_data, max_distance, hom_offsets, matrix_bin_size, matrix_size, merge_haplotypes):
    hom_names, loci_np_array, position_np_array = g3d_data.to_np_arrays()
    sys.stderr.write("[M::" + __name__ + "] preparing genome-wide KD tree\n")
    kdtree = spatial.KDTree(position_np_array)
    sys.stderr.write("[M::" + __name__ + "] finding particle pairs\n")
    index_pairs = kdtree.query_pairs(max_distance)
    sys.stderr.write("[M::" + __name__ + "] found " + str(len(index_pairs)) + " pairs of particles (" + str(round(float(len(index_pairs)) / g3d_data.num_g3d_particles(), 2)) + " pairs per particle)\n")
    
    # convert to matrix
    sys.stderr.write("[M::" + __name__ + "] generating a " + str(matrix_size) + " x " + str(matrix_size) + " matrix\n")
    matrix_data = np.zeros((matrix_size, matrix_size), dtype = int)
    for index_pair in index_pairs:
        index_1, index_2 = index_pair
        matrix_index_1 = g3d_np_arrays_to_matrix_index(hom_names, loci_np_array, index_1, hom_offsets, matrix_bin_size, merge_haplotypes)
        matrix_index_2 = g3d_np_arrays_to_matrix_index(hom_names, loci_np_array, index_2, hom_offsets, matrix_bin_size, merge_haplotypes)
        matrix_data[matrix_index_1, matrix_index_2] += 1
        matrix_data[matrix_index_2, matrix_index_1] += 1
    return matrix_data

def con3(argv):
    # default parameters
    max_distance = 3.0
    matrix_mode = False
    chr_len_file_name = None
    matrix_bin_size = 1000000
    merge_haplotypes = False
    info_mode = False
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "d:m:b:Hi")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c con3 [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -d FLOAT       max distance for generating a contact [" + str(max_distance) + "]\n")
        sys.stderr.write("  -m <chr.len>   output a matrix of binned counts based on chromosome lengths (tab-delimited: chr, len)\n")
        sys.stderr.write("  -b INT         bin size (bp) for \"-m\" (bins are centered around multiples of bin size) [" + str(matrix_bin_size) + "]\n")
        sys.stderr.write("  -H             merge the two haplotypes (for \"-m\")\n")
        sys.stderr.write("  -i             output bin info (tab-delimited: homolog or chr if \"-H\", bin center) instead (for \"-m\")\n")
        return 1
        
    num_color_schemes = 0
    for o, a in opts:
        if o == "-d":
            max_distance = float(a)
        elif o == "-m":
            matrix_mode = True
            chr_len_file_name = a
        elif o == "-b":
            matrix_bin_size = int(a)
        elif o == "-H":
            merge_haplotypes = True
        elif o == "-i":
            info_mode = True
                                                        
    # read 3DG file
    if not info_mode:
        g3d_data = file_to_g3d_data(open(args[0], "rb"))
        g3d_data.sort_g3d_particles()
        g3d_resolution = g3d_data.resolution()
        sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + ("N.A." if g3d_resolution is None else str(g3d_resolution)) + " bp resolution\n")
    
    # matrix mode
    if matrix_mode:
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
            matrix_data = g3d_data_to_matrix(g3d_data, max_distance, hom_offsets, matrix_bin_size, matrix_size, merge_haplotypes)
            np.savetxt(sys.stdout, matrix_data, fmt='%i', delimiter='\t')
        
    else:
        if not info_mode:
            con_data = g3d_data_to_con_data(g3d_data, max_distance)
            con_data.sort_cons()
            sys.stderr.write("[M::" + __name__ + "] writing output for " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
            sys.stdout.write(con_data.to_string() + "\n")
    
    return 0
    
