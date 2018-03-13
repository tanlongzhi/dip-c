import sys
import getopt
import gzip
import copy
from classes import Haplotypes, Reg, file_to_reg_list, get_phased_regs, G3dList, G3dParticle, G3dData, file_to_g3d_data, ref_name_haplotype_to_hom_name, hom_name_to_ref_name_haplotype
import numpy as np
from scipy.spatial.distance import pdist, squareform
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

def rg(argv):
    # default parameters
    output_prefix = None
    reg_file_name = None
    reg_list = []
    distance_mode = False
    
    # display parameters
    disp_num_particles = 100
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "o:r:d")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c rg [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -o STR        output prefix [<in.3dg>.]\n")
        sys.stderr.write("  -r <in.reg>   only analyze certain regions\n")
        sys.stderr.write("                  (will output two regions if haplotype is \".\")\n")
        sys.stderr.write("  -d            output pairwise distances instead\n")
        sys.stderr.write("Output:\n")
        sys.stderr.write("  <prefix><region_name>.rg    an N x N matrix of radii of gyration\n")
        sys.stderr.write("  <prefix><region_name>.loc   a list of N chromosomal loci\n")
        return 1
    for o, a in opts:
        if o == "-o":
            output_prefix = a
        elif o == "-r":
            reg_file_name = a
        elif o == "-d":
            distance_mode = True
            
    if output_prefix is None:
        output_prefix = args[0] + "."
        
    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + str(g3d_resolution) + " bp resolution\n")
    
    # prepare regions to analyze
    if reg_file_name is None:
        # analyze all homologs
        for hom_name in g3d_data.get_hom_names():
            ref_name, haplotype = hom_name_to_ref_name_haplotype(hom_name)
            reg = Reg(ref_name)
            reg.add_haplotype(haplotype)
            reg_list.append(reg)
    else:
        reg_file = gzip.open(reg_file_name, "rb") if reg_file_name.endswith(".gz") else open(reg_file_name, "rb")
        reg_list.extend(file_to_reg_list(reg_file))
        reg_file.close()
        reg_list = [reg for reg in get_phased_regs(reg_list)]
    
    sys.stderr.write("[M::" + __name__ + "] will analyze the following regions:\n")
    sys.stderr.write("name\tchr\thap\tstart\tend\n")
    for reg in reg_list:
        sys.stderr.write(reg.to_name_string() + "\t" + reg.to_string() + "\n")
    
    # calculate Rg matrix for each region
    for reg in reg_list:
        reg_name = reg.to_name_string()
        g3d_list = G3dList()
        for g3d_particle in g3d_data.get_g3d_particles_in_reg(reg):
            g3d_list.add_g3d_particle(g3d_particle)
        g3d_list.sort_g3d_particles()
        sys.stderr.write("[M::" + __name__ + "] processing region " + reg_name + ", with " + str(g3d_list.num_g3d_particles()) + " particles\n")

        loci_np_array, position_np_array = g3d_list.to_np_arrays()
        
        # write loci file
        loci_file_name = output_prefix + reg_name + ".loc"
        np.savetxt(loci_file_name, loci_np_array, fmt='%i', delimiter='\t')
        
        # calculate Rg
        rg_file_name = output_prefix + reg_name + ".rg"
        if distance_mode:
            output_matrix = squareform(pdist(position_np_array))
        else:
            output_matrix = position_np_array_to_rg_np_array(position_np_array, disp_num_particles)
        np.savetxt(rg_file_name, output_matrix, delimiter='\t')        
    
    return 0
    
