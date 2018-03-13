import sys
import getopt
import gzip
import copy
from classes import Haplotypes, Leg, string_to_leg, G3dData, file_to_g3d_data
import numpy as np

def pos(argv):
    # default parameters
    leg_file_name = None
    in_only = False
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "l:O")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c leg [options] -l <in.leg> <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -l <in.leg>    LEG file to convert to 3D positions (required)\n")
        sys.stderr.write("  -O             exclude out-of-bound legs\n")
        return 1
    for o, a in opts:
        if o == "-l":
            leg_file_name = a
        if o == "-O":
            in_only = True
    if leg_file_name is None:
        sys.stderr.write("[E::" + __name__ + "] -l is required\n")
        return 1
        
    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + ("N.A." if g3d_resolution is None else str(g3d_resolution)) + " bp resolution\n")
    g3d_data.prepare_interpolate()
                            
    # convert LEG file to 3DG particles
    for leg_file_line in open(leg_file_name, "rb"):
        is_out, position = g3d_data.interpolate_leg(string_to_leg(leg_file_line.strip()))
        if position is None or (is_out and in_only):
            sys.stdout.write("None\n")
        else:
            sys.stdout.write("\t".join(map(str, position)) + "\n")
    
    return 0
    
