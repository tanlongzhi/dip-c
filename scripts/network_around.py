import sys
import numpy as np

distance_threshold = float(sys.argv[1]) # max distance

# read pairwise distance graph
pairwise_distances = np.loadtxt(sys.argv[2], dtype = float, delimiter='\t')

# read row info from a LEG file
row_hom_names = []
row_ref_loci = []
for input_leg_file_line in open(sys.argv[3], "rb"):
    ref_name, ref_locus, haplotype = input_leg_file_line.strip().split(",")
    ref_locus = int(ref_locus)
    hom_name = ref_name + ("(pat)" if haplotype == "0" else "(mat)")
    row_hom_names.append(hom_name)
    row_ref_loci.append(ref_locus)
    
# read col info from a LEG file
col_hom_names = []
col_ref_loci = []
for input_leg_file_line in open(sys.argv[4], "rb"):
    ref_name, ref_locus, haplotype = input_leg_file_line.strip().split(",")
    ref_locus = int(ref_locus)
    hom_name = ref_name + ("(pat)" if haplotype == "0" else "(mat)")
    col_hom_names.append(hom_name)
    col_ref_loci.append(ref_locus)
    
# read row names
row_names = []
for name_line in open(sys.argv[5], "rb"):
    row_names.append(name_line.strip())

# read col names
col_names = []
for name_line in open(sys.argv[6], "rb"):
    col_names.append(name_line.strip())

# analyze
sys.stdout.write("#" + "\t".join(["name", "hom", "n_other_inter", "n_other_intra", "n_other_hom", "others", "other_homs"]) + "\n")
num_rows, num_cols = pairwise_distances.shape
for i in range(num_rows):
    near_col_names = []
    near_col_hom_names = []
    num_near_cols_inter = 0
    num_near_cols_intra = 0
    for j in range(num_cols): # no self edges
        if pairwise_distances[i, j] <= distance_threshold:
            if row_hom_names[i] == col_hom_names[j]:
                num_near_cols_intra += 1
            else:
                num_near_cols_inter += 1
            near_col_names.append(col_names[j])
            near_col_hom_names.append(col_hom_names[j])
                
    sys.stdout.write("\t".join([row_names[i], row_hom_names[i], str(num_near_cols_inter), str(num_near_cols_intra), str(len(set(near_col_hom_names))), ",".join(near_col_names), ",".join(near_col_hom_names)]) + "\n")
    