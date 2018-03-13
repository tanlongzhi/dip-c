import sys
import getopt
import gzip
import copy
from classes import Haplotypes, Reg, file_to_reg_list, get_phased_regs, G3dList, G3dParticle, G3dData, file_to_g3d_data, ref_name_haplotype_to_hom_name
import numpy as np
import math

# output a N-1 by 2 matrix of TAD tree (each row: start locus, end locus)
def tad_tree(rg_np_array, loci_np_array, disp_num_steps):
    num_loci = rg_np_array.shape[0]
    
    # initialize output
    tad_np_array = np.empty((num_loci - 1, 2), dtype = int)
    
    # initialize clusters (each row: start locus, end locus, whether already merged)
    clustered_np_array = np.vstack((np.arange(num_loci), np.arange(num_loci), np.zeros(num_loci, dtype = int))).transpose()
    
    # hierarchical clustering
    for cluster_step in range(num_loci - 1):
        if cluster_step % disp_num_steps == 0:
            sys.stderr.write("[M::" + __name__ + "] step " + str(cluster_step) + "\n")
        # find the best merge
        best_rg = None
        best_start = None
        best_end = None
        for start in range(num_loci - 1):
            if clustered_np_array[start, 2] > 0:
                continue # skip merged
            for end in range(start + 1, num_loci):
                if clustered_np_array[end, 2] > 0:
                    continue # skip merged
                rg = rg_np_array[clustered_np_array[start, 0], clustered_np_array[end, 1]]
                if best_rg is None or rg < best_rg:
                    best_rg = rg
                    best_start = start
                    best_end = end
                break
        #sys.stderr.write(str(best_start)+"-"+str(best_end)+":"+str(best_rg)+"\n")
        best_start_locus_id = clustered_np_array[best_start, 0]
        best_end_locus_id = clustered_np_array[best_end, 1]
        # record the best merge
        tad_np_array[cluster_step, :] = [loci_np_array[best_start_locus_id], loci_np_array[best_end_locus_id]]
        # merge
        clustered_np_array[best_start, 1] = best_end_locus_id
        clustered_np_array[best_end, 2] = 1
    
    return tad_np_array
    

def tad(argv):
    # default parameters
    rg_file_name = None
    loci_file_name = None
    
    # display parameters
    disp_num_steps = 100
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "l:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c tad [options] -l <in.loc> <in.rg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -l <in.loc>   input locus file (one locus per line)\n")
        sys.stderr.write("Output:\n")
        sys.stderr.write("  tab-delimited: start locus, end locus\n")
        return 1
    for o, a in opts:
        if o == "-l":
            loci_file_name = a

    if loci_file_name is None:
        sys.stderr.write("[E::" + __name__ + "] -l is required\n")
        return 1
                
    # read files
    loci_np_array = np.loadtxt(loci_file_name, dtype = int, delimiter='\t')
    sys.stderr.write("[M::" + __name__ + "] read a list of " + str(loci_np_array.shape[0]) + " loci\n")
    rg_np_array = np.loadtxt(args[0], delimiter='\t')
    sys.stderr.write("[M::" + __name__ + "] read a matrix of " + str(rg_np_array.shape[0]) + " times " + str(rg_np_array.shape[1]) + " radii of gyration\n")
    
    # hierarchical clustering
    np.savetxt(sys.stdout, tad_tree(rg_np_array, loci_np_array, disp_num_steps), fmt='%i', delimiter='\t')
    
    return 0
    
