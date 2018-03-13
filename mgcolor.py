import sys
import getopt
from classes import Haplotypes, hom_name_to_ref_name_haplotype

def append_color_data(merged_color_data, color_data, missing_value):
    # find the number of existing files
    if merged_color_data:
        num_color_files = len(merged_color_data.itervalues().next())
    else:
        num_color_files = 0
    
    # append the new file
    for key in set().union(*[merged_color_data, color_data]):
        if key in merged_color_data:
            if key in color_data:
                merged_color_data[key].append(color_data[key])
            else:
                merged_color_data[key].append(missing_value)
        else:
            merged_color_data[key] = [missing_value] * num_color_files
            merged_color_data[key].append(color_data[key])

def mgcolor(argv):
    # default parameters
    missing_value = -1.0
    diploid_mode = False
        
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "m:d")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c mgcolor [options] <in1.color> <in2.color> ...\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -m FLOAT          color for particles that are missing from some files [" + str(missing_value) + "]\n")
        sys.stderr.write("  -d                diploid mode: two columns for each file (pat and mat)\n\n")

        sys.stderr.write("Output:\n")
        sys.stderr.write("  tab-delimited with header: homolog (chr if \"-d\"), locus, in1 color, in2 color, ...\n")
        return 1
        
    for o, a in opts:
        if o == "-m":
            missing_value = float(a)
        elif o == "-d":
            diploid_mode = True

    # read color files
    header_line = "\t".join([("chr" if diploid_mode else "homolog"), "locus"])
    merged_color_data = {}
    num_color_files = 0
    for a in args:
        color_file = gzip.open(a, "rb") if a.endswith(".gz") else open(a, "rb")
        num_colors = 0
        
        # read a color file
        if diploid_mode:
            pat_color_data = {}
            mat_color_data = {}
        else:
            color_data = {}
        for color_file_line in color_file:
            hom_name, ref_locus, color = color_file_line.strip().split("\t")
            ref_locus = int(ref_locus)
            if diploid_mode:
                ref_name, haplotype = hom_name_to_ref_name_haplotype(hom_name)
                if haplotype == Haplotypes.paternal:
                    pat_color_data[(ref_name, ref_locus)] = color
                elif haplotype == Haplotypes.maternal:
                    mat_color_data[(ref_name, ref_locus)] = color
            else:
                color_data[(hom_name, ref_locus)] = color
        
        # merge with existing file
        if diploid_mode:
            append_color_data(merged_color_data, pat_color_data, missing_value)
            append_color_data(merged_color_data, mat_color_data, missing_value)
            sys.stderr.write("[M::" + __name__ + "] read color file " + a + " with " + str(len(pat_color_data)) + " paternal and " + str(len(mat_color_data)) + " maternal particles\n")
        else:
            append_color_data(merged_color_data, color_data, missing_value)
            sys.stderr.write("[M::" + __name__ + "] read color file " + a + " with " + str(len(color_data)) + " particles\n")
        
        # append header
        if diploid_mode:
            header_line += "\t" + a + ".pat\t" + a + ".mat"
        else:
            header_line += "\t" + a
        
        color_file.close()
        
    
    # output
    sys.stderr.write("[M::" + __name__ + "] writing " + str(len(merged_color_data)) + " particles from " + str(len(args)) + " files\n")
    sys.stdout.write(header_line + "\n")
    for hom_name, ref_locus in sorted(merged_color_data.keys()):
        sys.stdout.write("\t".join([hom_name, str(ref_locus)] + map(str, merged_color_data[(hom_name, ref_locus)])) + "\n")
    
    return 0
    
