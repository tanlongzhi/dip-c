import sys
import getopt
import gzip
from classes import Haplotypes, hom_name_to_ref_name_haplotype, ref_name_haplotype_to_hom_name, ConData, file_to_con_data, Con, Leg, G3dData, file_to_g3d_data
import numpy as np
from scipy import spatial
import math
    
def leg_to_index(leg, hom_offsets, bin_size):
    return hom_offsets[ref_name_haplotype_to_hom_name(leg.ref_name_haplotype())] + int(round(float(leg.get_ref_locus()) / bin_size))

def con_to_force(con, hom_offsets, bin_size):
    return np.array(sorted([leg_to_index(con.leg_1(), hom_offsets, bin_size), leg_to_index(con.leg_2(), hom_offsets, bin_size)]))
    
def init_graph(num_nodes):
    return (np.random.uniform(-100, 100, (num_nodes, 3)), np.zeros((num_nodes, 3), dtype = float))
    
def update_graph(positions, velocities, forces, force_constant, friction):
    num_nodes = positions.shape[0]
    force_values = np.zeros((num_nodes, 3), dtype = float)
    
    # calculate forces (spring with equilibrium distance 1)
    for i in range(forces.shape[0]):
        displacement = positions[forces[i, 1], ] - positions[forces[i, 0], ]
        distance = np.linalg.norm(displacement)
        if distance > 0.0:
            force_value = -force_constant * (distance - 1.0) * displacement / distance            
            #sys.stderr.write(str(force_value) + "\n")
            force_values[forces[i, 1], ] += force_value
            force_values[forces[i, 0], ] -= force_value
            
    # collision forces (same as normal forces, but only repulsion)
    kdtree = spatial.cKDTree(positions)
    index_pairs = kdtree.query_pairs(1.0)
    for index_pair in index_pairs:
        index_1, index_2 = index_pair
        displacement = positions[index_2, ] - positions[index_1, ]
        distance = np.linalg.norm(displacement)
        if distance > 0.0:
            force_value = -force_constant * (distance - 1.0) * displacement / distance            
            #sys.stderr.write(str(force_value) + "\n")
            force_values[index_2, ] += force_value
            force_values[index_1, ] -= force_value
                
    # apply forces
    velocities += force_values
    velocities *= 1 - friction
    positions += velocities
    
    return

def force(argv):
    # default parameters
    chr_len_file_name = None
    bin_size = 10000000
    max_step = 1000
    force_constant = 1e-2
    friction = 1e-1
    disp_num_steps = 100
    disp_num_cons = 100000
    output_prefix = "force."
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "l:b:k:f:o:n:w:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c force [options] -l <chr.len> <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -l <chr.len>   chromosome names and lengths (tab-delimited: chr, len) (required)\n")
        sys.stderr.write("  -o STR         output prefix [" + output_prefix + "]\n")
        sys.stderr.write("  -b INT         bin size (bp) (bins are centered around multiples of bin size) [" + str(bin_size) + "]\n")
        sys.stderr.write("  -k FLOAT       force constant [" + str(force_constant) + "]\n")
        sys.stderr.write("  -f FLOAT       friction [" + str(friction) + "]\n")
        sys.stderr.write("  -n INT         number of time steps [" + str(max_step) + "]\n")
        sys.stderr.write("  -w INT         number of time steps for each output [" + str(disp_num_steps) + "]\n")
        return 1
        
    num_color_schemes = 0
    for o, a in opts:
        if o == "-l":
            chr_len_file_name = a
        elif o == "-o":
            output_prefix = a
        elif o == "-b":
            bin_size = int(a)
        elif o == "-k":
            force_constant = float(a)
        elif o == "-f":
            friction = float(a)
        elif o == "-n":
            max_step = int(a)
        elif o == "-w":
            disp_num_steps = int(a)
            
    if chr_len_file_name is None:
        sys.stderr.write("[E::" + __name__ + "] -l is required\n")
        return 1
                                                        
    # read CON file
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")

    # initialize force matrix (id_1, id_2)
    forces = np.empty([0, 2], dtype = int)
    
    # read chromosome lengths
    # for position to node id
    hom_lens = {}
    hom_bin_lens = {}
    hom_offsets = {}
    num_nodes = 0
    
    # for node id to position
    hom_name_list = []
    ref_locus_list = []
    
    # read each chromosome and add backbone to force matrix (each row: id_1, id_2, with id_1 < id_2)
    chr_len_file = open(chr_len_file_name, "rb")
    for chr_len_file_line in chr_len_file:
        ref_name, ref_len = chr_len_file_line.strip().split("\t")
        ref_len = int(ref_len)
        for haplotype in [Haplotypes.paternal, Haplotypes.maternal]:
            hom_name = ref_name_haplotype_to_hom_name((ref_name, haplotype))
            hom_bin_len = int(round(float(ref_len) / bin_size)) + 1
            
            backbone_forces = np.arange(num_nodes, num_nodes + hom_bin_len - 1, 1).reshape((-1, 1))
            backbone_forces = np.hstack([backbone_forces, backbone_forces + 1])
            forces = np.vstack([forces, backbone_forces])
            
            hom_lens[hom_name] = ref_len
            hom_bin_lens[hom_name] = hom_bin_len
            hom_offsets[hom_name] = num_nodes
            num_nodes += hom_bin_len
            
            hom_name_list.extend([hom_name] * hom_bin_len)
            ref_locus_list.extend(range(0, bin_size * hom_bin_len, bin_size))
    sys.stderr.write("[M::" + __name__ + "] starting with " + str(num_nodes) + " nodes\n")
    sys.stderr.write("[M::" + __name__ + "] starting with " + str(forces.shape[0]) + " backbone forces\n")
    
    
    # convert CON to a force matrix
    con_forces = np.zeros((con_data.num_cons(), 2), dtype = int)
    counter = 0
    for con in con_data.get_cons():
        con_forces[counter, ] = con_to_force(con, hom_offsets, bin_size)            
        counter += 1
        if counter % disp_num_cons == 0:
            sys.stderr.write("[M::" + __name__ + "] converted " + str(counter) + " contacts\n")
            
    # remove redundant forces
    forces = np.vstack([forces, con_forces])
    forces = np.unique(forces, axis=0)
    sys.stderr.write("[M::" + __name__ + "] starting with a total of " + str(forces.shape[0]) + " forces\n")
    #sys.stderr.write(str(forces) + "\n")
    
    # initialize
    positions, velocities = init_graph(num_nodes)
    
    # run
    prev_positions = positions.copy()
    for step in range(max_step):
        update_graph(positions, velocities, forces, force_constant, friction)
        #sys.stderr.write(str(positions) + "\n")
        if (step + 1) % disp_num_steps == 0:
            sys.stderr.write("[M::" + __name__ + "] finished step " + str(step + 1) + ", average movement: " + str(np.mean(np.linalg.norm(positions - prev_positions, axis = 1))) + "\n")
            
            prev_positions = positions.copy()
            
            output_file = open(output_prefix + str(step + 1) + ".3dg", "wb")
            for i in range(num_nodes):
                output_file.write("\t".join(map(str, [hom_name_list[i], ref_locus_list[i], positions[i, 0], positions[i, 1], positions[i, 2]])) + "\n")
            output_file.close()
            prev_positions = positions.copy()
            
    
    return 0
    
