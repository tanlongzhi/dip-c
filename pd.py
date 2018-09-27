import sys
import getopt
import gzip
import copy
from classes import Haplotypes, Leg, string_to_leg, G3dData, file_to_g3d_data
import numpy as np
from scipy.spatial import distance

def pd(argv):
    # default parameters
    leg_file_1_name = None
    leg_file_2_name = None
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "1:2:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c pd [options] -1 <in1.leg> [-2 <in2.leg>] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -1 <in1.leg>    LEG file (required)\n")
        sys.stderr.write("  -2 <in2.leg>    LEG file [<in1.leg>]\n")
        return 1
    for o, a in opts:
        if o == "-1":
            leg_file_1_name = a
        elif o == "-2":
            leg_file_2_name = a
    if leg_file_1_name is None:
        sys.stderr.write("[E::" + __name__ + "] -1 is required\n")
        return 1
    if leg_file_2_name is None:
        leg_file_2_name = leg_file_1_name
                
    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + ("N.A." if g3d_resolution is None else str(g3d_resolution)) + " bp resolution\n")
    g3d_data.prepare_interpolate()
                            
    # convert LEG file to 3DG particles
    positions_1 = np.empty([0, 3])
    for leg_file_1_line in open(leg_file_1_name, "rb"):
        is_out, position = g3d_data.interpolate_leg(string_to_leg(leg_file_1_line.strip()))
        if position is None:
            position = np.array([np.nan, np.nan, np.nan])
        positions_1 = np.vstack([positions_1, position])
        
    positions_2 = np.empty([0, 3])
    for leg_file_2_line in open(leg_file_2_name, "rb"):
        is_out, position = g3d_data.interpolate_leg(string_to_leg(leg_file_2_line.strip()))
        if position is None:
            position = np.array([np.nan, np.nan, np.nan])
        positions_2 = np.vstack([positions_2, position])

    # calculate pairwise distances
    distances = distance.cdist(positions_1, positions_2)
    np.savetxt(sys.stdout, distances, delimiter='\t')
    
    return 0
    
