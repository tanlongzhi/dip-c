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

# for object name in pymol
def hom_name_to_object_name(hom_name):
    return hom_name.replace("(", "_").replace(")", "")

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
        sys.stderr.write("Usage: dip-c exp [options] <in.3dg>\n")
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
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + ("N.A." if g3d_resolution is None else str(g3d_resolution)) + " bp resolution\n")
    
    # center of nucleus
    nuc_center = center_g3d_particles(g3d_data.get_g3d_particles())
    
    # process data
    if centers_only:
        center_g3d_data = G3dData()
        for hom_name in g3d_data.get_hom_names():
            center_position = center_g3d_particles(g3d_data.get_g3d_particles_from_hom_name(hom_name))
            center_position += (center_position - nuc_center) * expansion_factor
            center_g3d_data.add_g3d_particle(G3dParticle(hom_name, 0, center_position.tolist()))
        g3d_data = center_g3d_data
    else:
        hom_centers = {}
        # center of each homologs
        for hom_name in g3d_data.get_hom_names():
            hom_centers[hom_name] = center_g3d_particles(g3d_data.get_g3d_particles_from_hom_name(hom_name))
            sys.stderr.write("extract " + hom_name_to_object_name(hom_name) + ", chain \"" + hom_name + "\"\n")
        # translate
        for hom_name in g3d_data.get_hom_names():
            translation_vector = (hom_centers[hom_name] - nuc_center) * expansion_factor
            for g3d_particle in g3d_data.get_g3d_particles_from_hom_name(hom_name):
                g3d_particle.set_position((np.array(g3d_particle.get_position()) + translation_vector).tolist())
            #sys.stderr.write("translate [" + ",".join(map(str, translation_vector)) + "], chain \"" + hom_name + "\"\n")
            sys.stderr.write("translate [" + ",".join(map(str, translation_vector)) + "], object=" + hom_name_to_object_name(hom_name) + ", camera=0\n")
        for hom_name in g3d_data.get_hom_names():
            sys.stderr.write("mview store, object=" + hom_name_to_object_name(hom_name)  + "\n")            
    
    # output
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(g3d_data.num_g3d_particles()) + " particles\n")
    sys.stdout.write(g3d_data.to_string()+"\n")
    
    return 0
    
