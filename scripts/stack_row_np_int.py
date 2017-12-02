import sys
import numpy as np

num_files = 0
row_id = int(sys.argv[1])
sys.stderr.write("row to extract: " + int(row_id) + "\n")
for file_name in sys.argv[2:]:
    sys.stderr.write("reading file: " + file_name + "\n")
    file_data = np.loadtxt(file_name, dtype = int, delimiter='\t')
    file_data = file_data[row_id, :]
    if num_files == 0:
        sum_data = file_data
    else:
        sum_data = np.vstack((sum_data, file_data))
    num_files += 1
sys.stderr.write("writing output\n")
np.savetxt(sys.stdout, sum_data, fmt='%i', delimiter='\t')