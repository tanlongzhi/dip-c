import sys
import getopt
import gzip
import copy
from classes import Haplotypes, G3dData, file_to_g3d_data, G3dParticle
import numpy as np
import math

# find center of mass
def center_g3d_particles(g3d_particles):
    sum_zero = 0
    sum_one = np.zeros([3], dtype=float)
    for g3d_particle in g3d_particles:
        sum_zero += 1
        position = np.array(g3d_particle.get_position())
        sum_one += position
    sum_one /= sum_zero
    return sum_one

def exp(argv):
    # default parameters
    expansion_factor = 3.0
    centers_only = False
     
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "f:c")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac exp [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -f FLOAT     expansion factor for translating away from nuclear center [" + str(expansion_factor) + "]\n")
        sys.stderr.write("  -c           output centers of mass\n")
        return 1
    for o, a in opts:
        if o == "-f":
            expansion_factor = float(a)
        if o == "-c":
            centers_only = True      
                          
    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + str(g3d_resolution) + " bp resolution\n")
    
    # center of nucleus
    nuc_center = center_g3d_particles(g3d_data.get_g3d_particles())
    
    # center of each homologs
    if centers_only:
        center_g3d_data = G3dData()
        for ref_name_haplotype in g3d_data.get_ref_name_haplotype():
            center_position = center_g3d_particles(g3d_data.get_g3d_particles_from_ref_name_haplotype(ref_name_haplotype))
            center_position += (center_position - nuc_center) * expansion_factor
            center_g3d_data.add_g3d_particle(G3dParticle(ref_name_haplotype[0], ref_name_haplotype[1], 0, center_position.tolist()))
        g3d_data = center_g3d_data
    else:
        ref_name_haplotype_centers = {}
        for ref_name_haplotype in g3d_data.get_ref_name_haplotype():
            ref_name_haplotype_centers[ref_name_haplotype] = center_g3d_particles(g3d_data.get_g3d_particles_from_ref_name_haplotype(ref_name_haplotype))
        # translate
        for ref_name_haplotype in g3d_data.get_ref_name_haplotype():
            translation_vector = (ref_name_haplotype_centers[ref_name_haplotype] - nuc_center) * expansion_factor
            for g3d_particle in g3d_data.get_g3d_particles_from_ref_name_haplotype(ref_name_haplotype):
                g3d_particle.set_position((np.array(g3d_particle.get_position()) + translation_vector).tolist())
    
    # output
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(g3d_data.num_g3d_particles()) + " particles\n")
    sys.stdout.write(g3d_data.to_string()+"\n")
    
    return 0
    