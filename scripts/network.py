import sys
import numpy as np
from scipy.sparse.csgraph import csgraph_from_dense, connected_components

distance_threshold = 3.0 # max distance to form an edge
separation_threshold = 20e6 # min separation (bp) for long-range intra

# read pairwise distance graph
pairwise_distances = np.loadtxt(sys.argv[1], dtype = float, delimiter='\t')

# read node info from a LEG file
node_hom_names = []
node_ref_loci = []
for input_leg_file_line in open(sys.argv[2], "rb"):
    ref_name, ref_locus, haplotype = input_leg_file_line.strip().split(",")
    ref_locus = int(ref_locus)
    hom_name = ref_name + ("(pat)" if haplotype == "0" else "(mat)")
    node_hom_names.append(hom_name)
    node_ref_loci.append(ref_locus)
    
# read node names
node_names = []
for node_name_line in open(sys.argv[3], "rb"):
    node_names.append(node_name_line.strip())

# form edges: 0=no, 1=short-range intra, 2=long-range intra, 3=inter
num_nodes = pairwise_distances.shape[0]
adjacency_matrix = np.zeros((num_nodes, num_nodes))
for i in range(num_nodes):
    for j in range(i + 1, num_nodes): # no self edges
        if pairwise_distances[i, j] <= distance_threshold:
            if node_hom_names[i] != node_hom_names[j]:
                adjacency_matrix[i,j] = 3
                adjacency_matrix[j,i] = 3
            elif abs(node_ref_loci[i] - node_ref_loci[j]) >= separation_threshold:
                adjacency_matrix[i,j] = 2
                adjacency_matrix[j,i] = 2
            else:
                adjacency_matrix[i,j] = 1
                adjacency_matrix[j,i] = 1

# get connected components
sparse_graph = csgraph_from_dense(adjacency_matrix)
num_comps, label_ids = connected_components(sparse_graph, directed=False)

# parse connected components and output
sys.stdout.write("#" + "\t".join(["n_node", "n_edge_inter", "n_edge_long_intra", "n_edge_short_intra", "n_hom", "homs", "nodes"]) + "\n")
for m in range(num_comps):
    comp_node_ids, = np.where(label_ids == m)
    comp_num_nodes = comp_node_ids.size
    comp_num_edge_inter = 0
    comp_num_edge_long_intra = 0
    comp_num_edge_short_intra = 0
    comp_names = []
    comp_hom_names = []
    for i in range(comp_num_nodes):
        index_1 = comp_node_ids[i]
        for j in range(i + 1, comp_num_nodes):
            index_2 = comp_node_ids[j]
            if adjacency_matrix[index_1, index_2] == 1:
                comp_num_edge_short_intra += 1
            elif adjacency_matrix[index_1, index_2] == 2:
                comp_num_edge_long_intra += 1
            elif adjacency_matrix[index_1, index_2] == 3:
                comp_num_edge_inter += 1
        comp_names.append(node_names[index_1])
        comp_hom_names.append(node_hom_names[index_1])
            
    sys.stdout.write("\t".join([str(comp_num_nodes), str(comp_num_edge_inter), str(comp_num_edge_long_intra), str(comp_num_edge_short_intra), str(len(set(comp_hom_names))), ",".join(comp_hom_names), ",".join(comp_names)]) + "\n")
    