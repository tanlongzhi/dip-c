import sys
import getopt

def mgcolor(argv):
    # default parameters
    missing_value = -1.0
        
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "m:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac mgcolor [options] <in1.color> <in2.color> ...\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -m FLOAT          color for particles that are missing from some files [" + str(missing_value) + "]\n\n")

        sys.stderr.write("Output:\n")
        sys.stderr.write("  tab-delimited with header: homolog, locus, in1 color, in2 color, ...\n")
        return 1
        
    for o, a in opts:
        if o == "-m":
            missing_value = float(a)

    # read color files
    merged_color_data = {}
    num_color_files = 0
    for a in args:
        color_file = gzip.open(a, "rb") if a.endswith(".gz") else open(a, "rb")
        num_colors = 0
        
        # read a color file
        color_data = {}
        for color_file_line in color_file:
            hom_name, ref_locus, color = color_file_line.strip().split("\t")
            ref_locus = int(ref_locus)
            color_data[(hom_name, ref_locus)] = color
        
        # merge with existing file
        for key in set().union(*[merged_color_data, color_data]):
            if key in merged_color_data:
                if key in color_data:
                    merged_color_data[key].append(color_data[key])
                else:
                    merged_color_data[key].append(missing_value)
            else:
                merged_color_data[key] = [0] * num_color_files
                merged_color_data[key].append(color_data[key])
        num_color_files += 1
        
        sys.stderr.write("[M::" + __name__ + "] read color file " + a + " with " + str(len(color_data)) + " particles\n")
        color_file.close()
        
    
    # output
    sys.stderr.write("[M::" + __name__ + "] writing " + str(len(merged_color_data)) + " particles from " + str(num_color_files) + " files\n")
    sys.stdout.write("\t".join(["homolog", "locus"] + args) + "\n")
    for hom_name, ref_locus in sorted(merged_color_data.keys()):
        sys.stdout.write("\t".join([hom_name, str(ref_locus)] + map(str, merged_color_data[(hom_name, ref_locus)])) + "\n")
    
    return 0
    