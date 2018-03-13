import sys
import getopt
import gzip
import copy
from classes import Haplotypes, LegData, ConData, file_to_con_data, Leg, Par, ParData, G3dData, file_to_g3d_data
import numpy as np
import math

def clean3(argv):
    # default parameters
    con_file_name = None
    max_clean_distance = 500000
    clean_quantile = 0.06
    
    # display parameters
    display_quantiles = np.arange(0.0, 1.01, 0.01)
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "c:d:q:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c clean3 [options] -c <in.con> <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -c <in.con>    contact file for cleaning (required)\n")
        sys.stderr.write("  -d INT         max distance (bp) from a contact leg to a 3D genome particle [" + str(max_clean_distance) + "]\n")
        sys.stderr.write("  -q FLOAT       quantile of particles to remove [" + str(clean_quantile) + "]\n")
        return 1
    for o, a in opts:
        if o == "-c":
            con_file_name = a
        elif o == "-d":
            max_clean_distance = int(a)
        elif o == "-q":
            clean_quantile = float(a)

    if con_file_name is None:
        sys.stderr.write("[E::" + __name__ + "] -c is required\n")
        return 1
        
    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + ("N.A." if g3d_resolution is None else str(g3d_resolution)) + " bp resolution\n")
                            
    # read legs from CON file
    con_file = gzip.open(con_file_name, "rb") if con_file_name.endswith(".gz") else open(con_file_name, "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    leg_data = LegData()
    leg_data.add_con_data(con_data)
    leg_data.sort_phased_legs()
    sys.stderr.write("[M::" + __name__ + "] sorted " + str(leg_data.num_legs()) + " legs\n")
    
    # find cut-off
    leg_counts = g3d_data.leg_counts(leg_data, max_clean_distance)
    sys.stderr.write("[M::" + __name__ + "] statistics:\n")
    sys.stderr.write("quantile\t#legs\n")
    for display_quantile in display_quantiles:
        sys.stderr.write(str(display_quantile) + "\t" + str(int(round(np.percentile(leg_counts, display_quantile*100.0), 0))) + "\n")
    min_leg_count = int(math.ceil(np.percentile(leg_counts, clean_quantile*100.0)))
    sys.stderr.write("[M::" + __name__ + "] min leg count: " + str(min_leg_count) + "\n")
    
    # clean
    g3d_data.clean_leg_poor(leg_data, max_clean_distance, min_leg_count)
    leg_counts = g3d_data.leg_counts(leg_data, max_clean_distance)
    g3d_data.sort_g3d_particles()
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(g3d_data.num_g3d_particles()) + " particles\n")
    sys.stdout.write(g3d_data.to_string()+"\n")
    
    
    return 0
    
