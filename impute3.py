import sys
import getopt
import gzip
import copy
from classes import Haplotypes, ConData, file_to_con_data, Leg, Par, ParData, G3dData, file_to_g3d_data

def impute3(argv):
    # default parameters
    g3d_file_name = None
    vio_file_name = None
    max_impute3_distance = 20
    max_impute3_ratio = 0.5
    min_impute3_separation_factor = 1.0
    max_clean_distance = 10000000
    min_clean_count = 2
    is_male = False
    par_data = None
    
    # presets
    h_par = ParData("X", "Y")
    h_par.add_par(Par("X", 60000, 2699520, "Y", 10000))
    h_par.add_par(Par("X", 154931043, 155260560, "Y", 59034049))

    m_par = ParData("chrX", "chrY")
    m_par.add_par(Par("chrX", 169969758, 170931299, "chrY", 90745844))
    
    presets = {
        "f":[False, None],
        "hf":[False, None],
        "mf":[False, None],
        "hm":[True, h_par],
        "mm":[True, m_par]}
    preset_descriptions = {
        "f":"female",
        "hf":"human female (same as f)",
        "mf":"mouse female (same as f)",
        "hm":"human male (hg19, no \"chr\" prefix)",
        "mm":"mouse male (mm10)"}
    
    # progress display parameters
    display_max_num_legs = 20
    display_num_cons = 10000
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "3:v:d:r:s:D:C:p:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac impute3 [options] -3 <in.3dg> [-v <out.vio>] <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -3 <in.3dg>    3D genome file for imputing haplotypes (required)\n")
        sys.stderr.write("  -v <out.vio>   output statistics to a contact violation file:\n")
        sys.stderr.write("                   tab-delimited: leg 1, leg 2, num of compatible haplotypes,\n")
        sys.stderr.write("                   shortest 3D distance, ratio between the shortest and the 2nd shortest distance\n\n")
        sys.stderr.write("  -d FLOAT       max 3D distance for imputing haplotypes [" + str(max_impute3_distance) + "]\n")
        sys.stderr.write("  -r FLOAT       max ratio between 3D distances for the best and 2nd best haplotypes [" + str(max_impute3_ratio) + "]\n")
        sys.stderr.write("  -s FLOAT       min separation (unit: 3D genome resolution) for imputing\n")
        sys.stderr.write("                   completely unphased, intra-chromosomal contacts [" + str(min_impute3_separation_factor) + "]\n\n")
        sys.stderr.write("  -D INT         max distance (bp, L-1/2 norm) for removing isolated contacts [" + str(max_clean_distance) + "]\n")
        sys.stderr.write("  -C INT         min neighbor count for an unisolated contact [" + str(min_clean_count) + "]\n\n")
        sys.stderr.write("  -p STR         presets for PARs and sex: [f]\n")
        for preset in sorted(presets.keys()):
            sys.stderr.write("                   " + preset + " = " + preset_descriptions[preset] + "\n")
        return 1
    for o, a in opts:
        if o == "-3":
            g3d_file_name = a
        elif o == "-v":
            vio_file_name = a
        elif o == "-d":
            max_impute3_distance = float(a)
        elif o == "-r":
            max_impute3_ratio = float(a)
        elif o == "-s":
            min_impute3_separation_factor = float(a)
        elif o == "-D":
            max_clean_distance = int(a)
        elif o == "-C":
            min_clean_count = int(a)  
        elif o == "-p":
            try:
                is_male, par_data = presets[a]
                sys.stderr.write("[M::" + __name__ + "] use preset " + a + " = " + preset_descriptions[a] + "\n")
            except KeyError:
                sys.stderr.write("[E::" + __name__ + "] unknown preset\n")
                return 1   
    if g3d_file_name is None:
        sys.stderr.write("[E::" + __name__ + "] -3 is required\n")
        return 1
        
    # read 3DG file
    g3d_data = file_to_g3d_data(open(g3d_file_name, "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + ("N.A." if g3d_resolution is None else str(g3d_resolution)) + " bp resolution\n")
    g3d_data.prepare_interpolate()
                            
    # read CON file
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    
    # impute3
    vio_file = None
    if not vio_file_name is None:
        vio_file = open(vio_file_name, "wb")
    con_data.impute_from_g3d_data(g3d_data, max_impute3_distance, max_impute3_ratio, max_impute3_ratio * g3d_resolution, is_male, par_data, vio_file)
    if not vio_file is None:
        vio_file.close()
    sys.stderr.write("[M::" + __name__ + "] imputed " + str(con_data.num_phased_cons()) + " contacts (" + str(round(100.0 * con_data.num_phased_cons() / con_data.num_cons(), 2)) + "%)\n")
    
    # clean imputed
    con_data.sort_cons()
    con_data.clean_unphased()
    before_clean_num_cons = con_data.num_cons()
    con_data.clean_isolated_phased(copy.deepcopy(con_data), max_clean_distance, min_clean_count)
    after_clean_num_cons = con_data.num_cons()
    sys.stderr.write("[M::" + __name__ + "] removed " + str(before_clean_num_cons - after_clean_num_cons) + " isolated contacts (" + str(round(100.0 * (before_clean_num_cons - after_clean_num_cons) / before_clean_num_cons, 2)) + "%)\n")
    
    # write output
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(con_data.num_cons()) + " contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    sys.stdout.write(con_data.to_string()+"\n")
    
    return 0
    
