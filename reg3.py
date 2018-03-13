import sys
import getopt
import gzip
from classes import Reg, file_to_reg_list, file_to_g3d_data, G3dData

def reg3(argv):
    # default parameters
    use_presets = False
    inc_file_name = None
    exc_file_name = None
    inc_regs = []
    exc_regs = []
    
    # presets
    presets = {
        "hm":[Reg(str(i + 1)) for i in range(22)] + [Reg("X"), Reg("Y")],
        "hf":[Reg(str(i + 1)) for i in range(22)] + [Reg("X")],
        "mm":[Reg("chr" + str(i + 1)) for i in range(19)] + [Reg("chrX"), Reg("chrY")],
        "mf":[Reg("chr" + str(i + 1)) for i in range(19)] + [Reg("chrX")]}
    preset_descriptions = {
        "hm":"human male (no \"chr\" prefix)",
        "hf":"human female (no \"chr\" prefix)",
        "mm":"mouse male",
        "mf":"mouse female"}
    
    # progress display parameters
    max_num_dups = 10
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "i:e:p:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c reg3 [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -i <inc.reg>   the only regions to include:\n")
        sys.stderr.write("                   tab-delimited: chr, haplotype, start, end\n")
        sys.stderr.write("                   haplotype: 0 = paternal, 1 = maternal, . = both\n")
        sys.stderr.write("                   start/end: . = denotes end of chromosome\n")
        sys.stderr.write("  -e <exc.reg>   regions to exclude (same format as above)\n")
        sys.stderr.write("  -p STR         presets for included regions:\n")
        for preset in sorted(presets.keys()):
            sys.stderr.write("                   " + preset + " = " + preset_descriptions[preset] + "\n")
        return 1
    for o, a in opts:
        if o == "-p":
            try:
                inc_regs = presets[a]
                sys.stderr.write("[M::" + __name__ + "] use preset " + a + " = " + preset_descriptions[a] + "\n")
            except KeyError:
                sys.stderr.write("[E::" + __name__ + "] unknown preset\n")
                return 1
            use_presets = True
        elif o == "-i":
            use_presets = False
            inc_file_name = a
        elif o == "-e":
            use_presets = False  
            exc_file_name = a

                
    # load regions from file
    if not inc_file_name is None:
        inc_file = gzip.open(inc_file_name, "rb") if inc_file_name.endswith(".gz") else open(inc_file_name, "rb")
        inc_regs.extend(file_to_reg_list(inc_file))
        inc_file.close()
    if not exc_file_name is None:
        exc_file = gzip.open(exc_file_name, "rb") if exc_file_name.endswith(".gz") else open(exc_file_name, "rb")
        exc_regs.extend(file_to_reg_list(exc_file))
        exc_file.close()        
        
    if not inc_regs:
        sys.stderr.write("[E::" + __name__ + "] must specify included regions\n")
        return 1
            
    # display regions
    sys.stderr.write("[M::" + __name__ + "] only include the following regions:\n")
    sys.stderr.write("chr\thap\tstart\tend\n")
    for reg in inc_regs:
        sys.stderr.write(reg.to_string() + "\n")
    sys.stderr.write("[M::" + __name__ + "] only exclude the following regions:\n")
    sys.stderr.write("chr\thap\tstart\tend\n")
    for reg in exc_regs:
        sys.stderr.write(reg.to_string() + "\n")
            
    # read CON file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + ("N.A." if g3d_resolution is None else str(g3d_resolution)) + " bp resolution\n")
    
    # apply regions
    g3d_data.apply_regs(inc_regs, exc_regs)
    g3d_data.sort_g3d_particles()
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(g3d_data.num_g3d_particles()) + " particles\n")
    sys.stdout.write(g3d_data.to_string()+"\n")
    
    return 0
    
