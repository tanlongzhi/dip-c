import sys
import numpy as np
from scipy.stats import rankdata

hom_names = []
ref_locus_strings = []
color_values = []
for file_line in open(sys.argv[1], "rb"):
    file_line_data = file_line.strip().split("\t")
    hom_names.append(file_line_data[0])
    ref_locus_strings.append(file_line_data[1])
    color_values.append(float(file_line_data[2]))

# rank
rank_color_values = rankdata(np.asarray(color_values))
# convert to 0-1
rank_color_values = (rank_color_values - 1) / (float(len(color_values) - 1))

sys.stderr.write("read " + str(len(color_values)) + " color values\n")

sys.stderr.write("writing output\n")
for i in range(len(hom_names)):
    sys.stdout.write("\t".join([hom_names[i], ref_locus_strings[i], str(rank_color_values[i])]) + "\n")
