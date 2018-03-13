import sys
import getopt
import gzip
import copy
from classes import Haplotypes, Reg, file_to_reg_list, get_phased_regs, G3dList, G3dParticle, G3dData, file_to_g3d_data, ref_name_haplotype_to_hom_name
import numpy as np
from scipy.spatial.distance import pdist
import math

def position_np_array_to_rg_np_array(position_np_array, disp_num_particles):
    num_particles = position_np_array.shape[0]
    sum_zero_np_array = np.zeros([num_particles, num_particles, 3], dtype=float)
    sum_one_np_array = np.zeros([num_particles, num_particles, 3], dtype=float)
    sum_two_np_array = np.zeros([num_particles, num_particles, 3], dtype=float)
    for i in range(num_particles):
        if i % disp_num_particles == 0:
            sys.stderr.write("[M::" + __name__ + "] analyzed " + str(i) + " particles (" + str(round(100.0 * i / num_particles, 2)) + "%)\n")
        sum_zero_np_array[:(i + 1), i:, :] += 1
        sum_zero_np_array[i:, :(i + 1), :] += 1
        sum_zero_np_array[i, i, :] -= 1        
        position = position_np_array[i, :]
        sum_one_np_array[:(i + 1), i:, :] += position
        sum_one_np_array[i:, :(i + 1), :] += position
        sum_one_np_array[i, i, :] -= position
        squared_position = position ** 2
        sum_two_np_array[:(i + 1), i:, :] += squared_position
        sum_two_np_array[i:, :(i + 1), :] += squared_position
        sum_two_np_array[i, i, :] -= squared_position
    return np.sqrt(np.sum(sum_two_np_array / sum_zero_np_array - (sum_one_np_array / sum_zero_np_array) ** 2, axis = 2))

def dist(argv):
    # default parameters

    
    # display parameters
    disp_num_particles = 100
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "d")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c dist [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -d            diploid mode\n")
        sys.stderr.write("Output:\n")
        sys.stderr.write("  tab-delimited: homolog (chr if \"-d\"), separation (in bp), #pairs, mean distance, r.m.s. distance\n")

        return 1
    for o, a in opts:
        if o == "-o":
            output_prefix = a
        elif o == "-r":
            reg_file_name = a

    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + str(g3d_resolution) + " bp resolution\n")
    
    # analyze each homolog
    for hom_name in g3d_data.get_hom_names():
        sys.stderr.write("analyzing " + hom_name + "\n")
        loci_np_array, position_np_array = g3d_data.get_g3d_list_from_hom_name(hom_name).to_np_arrays()
        sep_np_array = pdist(loci_np_array)
        dist_np_array = pdist(position_np_array)
        uniq_seps, uniq_indices = np.unique(sep_np_array, return_inverse = True)
        num_seps = len(uniq_seps)
        
        # calculate statistics
        nums_pairs = [0] * num_seps
        sums = [0.0] * num_seps
        sums_sq = [0.0] * num_seps
        for i in range(len(sep_np_array)):
            sep_index = uniq_indices[i]
            nums_pairs[sep_index] += 1
            sums[sep_index] += dist_np_array[i]
            sums_sq[sep_index] += dist_np_array[i] ** 2
        
        # print
        for i in range(num_seps):
            sys.stdout.write("\t".join([hom_name, str(int(uniq_seps[i])), str(nums_pairs[i]), str(sums[i] / nums_pairs[i]), str(math.sqrt(sums_sq[i] / nums_pairs[i]))]) + "\n")
        
    
    return 0
    
