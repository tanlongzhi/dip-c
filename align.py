import sys
import getopt
import copy
import rmsd
import numpy as np

def align(argv):
    # default parameters
    output_prefix = None

    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "o:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c align [options] <in1.3dg> <in2.3dg> ...\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -o STR        output prefix [no output]\n")
        sys.stderr.write("Output:\n")
        sys.stderr.write("  tab-delimited: homolog, locus, RMSD\n")
        sys.stderr.write("  additionally with \"-o\": 3DG files aligned to each other\n")

        return 1
    for o, a in opts:
        if o == "-o":
            output_prefix = a            
            
    # load 3dg files
    input_data = []
    num_structures = len(args)
    if num_structures < 2:
        sys.stderr.write("[E::" + __name__ + "] at least 2 structures are required\n")
        return 1
    counter = 0
    for input_filename in args:
        sys.stderr.write("[M::" + __name__ + "] reading 3dg file " + str(counter) + ": " + input_filename + "\n")
        input_data.append({})
        for input_file_line in open(input_filename, "rb"):
            input_file_line_data = input_file_line.strip().split()
            input_data[-1][(input_file_line_data[0], int(input_file_line_data[1]))] = [float(input_file_line_data[2]),float(input_file_line_data[3]),float(input_file_line_data[4])]
        counter += 1
    
    # find common particles
    common_loci = set(input_data[0])
    for input_structure in input_data[1:]:
        common_loci = common_loci.intersection(set(input_structure))
    num_loci = len(common_loci)
    common_loci = list(common_loci)
    common_data = []
    for input_structure in input_data:
        common_data.append([])
        for common_locus in common_loci:
            common_data[-1].append(input_structure[common_locus])
    sys.stderr.write("[M::" + __name__ + "] found " + str(num_loci) + " common particles\n")
        
    # subtract centroid
    common_data = np.array(common_data)
    centroid_data = []
    for i in range(num_structures):
        common_data[i] = np.array(common_data[i])
        centroid_pos = rmsd.centroid(common_data[i])
        common_data[i] -= centroid_pos
        centroid_data.append(centroid_pos)
    sys.stderr.write("[M::" + __name__ + "] found centroids for " + str(num_structures) + " structures\n")
    
    # calculate pairwise deviation and rotate
    deviations = np.empty((num_loci, 0), float)
    for i in range(num_structures):
        for j in range(num_structures):
            if j == i:
                continue
            # mirror image if needed
            mirror_factor = 1.0
            if rmsd.kabsch_rmsd(common_data[i], common_data[j]) > rmsd.kabsch_rmsd(common_data[i], -1.0 * common_data[j]):
                mirror_factor = -1.0
            # calculate deviation
            rotation_matrix = rmsd.kabsch(mirror_factor * common_data[j], common_data[i])
            if j > i:
                deviation = np.linalg.norm(np.dot(mirror_factor * common_data[j], rotation_matrix) - common_data[i], axis = 1).T
                deviations = np.c_[deviations, deviation]
                sys.stderr.write("[M::" + __name__ + "] median deviation between file " + str(i) + " and file " + str(j) + ": " + str(np.median(deviation)) + "\n")
            
            # rotate
            if output_prefix is not None:
                # rotate j to align with i
                sys.stderr.write("[M::" + __name__ + "] aligning file " + str(j) + " to file " + str(i) + "\n")
                aligned_filename = output_prefix + str(j) + "_to_" + str(i) + ".3dg"
                aligned_file = open(aligned_filename, "wb")
                for input_locus in input_data[j]:
                    aligned_pos = np.dot((np.array(input_data[j][input_locus]) - centroid_data[j]) * mirror_factor, rotation_matrix) + centroid_data[i]
                    aligned_file.write("\t".join([input_locus[0], str(input_locus[1]), str(aligned_pos[0]), str(aligned_pos[1]), str(aligned_pos[2])]) + "\n")
                aligned_file.close()

    # summarize rmsd and print
    rmsds = np.sqrt((deviations ** 2).mean(axis = 1))
    totalrmsd = np.sqrt((rmsds ** 2).mean(axis = 0))
    sys.stderr.write("[M::" + __name__ + "] RMS RMSD: " + str(totalrmsd) + "\n")
    sys.stderr.write("[M::" + __name__ + "] median RMSD: " + str(np.median(rmsds,axis = 0)) + "\n")
    sys.stderr.write("[M::" + __name__ + "] writing output\n")
    for i in range(num_loci):
        sys.stdout.write("\t".join(map(str, [common_loci[i][0], common_loci[i][1], rmsds[i]])) + "\n")
        
    return 0
    