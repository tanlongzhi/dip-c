import sys
import getopt
import gzip
from classes import ConData, file_to_con_data, Reg, file_to_reg_list

def reg(argv):
    # default parameters
    use_presets = False
    inc_file_name = None
    exc_file_name = None
    hap_file_name = None
    inc_regs = []
    exc_regs = []
    hap_regs = []
    
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
        opts, args = getopt.getopt(argv[1:], "i:e:p:h:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c reg [options] <in.con>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -i <inc.reg>   the only regions to include:\n")
        sys.stderr.write("                   tab-delimited: chr, haplotype, start, end\n")
        sys.stderr.write("                   haplotype: 0 = paternal, 1 = maternal, . = both\n")
        sys.stderr.write("                   start/end: . = denotes end of chromosome\n")
        sys.stderr.write("  -e <exc.reg>   regions to exclude (same format as above)\n")
        sys.stderr.write("  -h <hap.reg>   regions to assume haploid (haplotype must be known)\n")
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
        elif o == "-h":
            use_presets = False  
            hap_file_name = a
                
    # load regions from file
    if not inc_file_name is None:
        inc_file = gzip.open(inc_file_name, "rb") if inc_file_name.endswith(".gz") else open(inc_file_name, "rb")
        inc_regs.extend(file_to_reg_list(inc_file))
        inc_file.close()
    if not exc_file_name is None:
        exc_file = gzip.open(exc_file_name, "rb") if exc_file_name.endswith(".gz") else open(exc_file_name, "rb")
        exc_regs.extend(file_to_reg_list(exc_file))
        exc_file.close()        
    if not hap_file_name is None:
        hap_file = gzip.open(hap_file_name, "rb") if hap_file_name.endswith(".gz") else open(hap_file_name, "rb")
        hap_regs.extend(file_to_reg_list(hap_file))
        hap_file.close() 
        
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
    sys.stderr.write("[M::" + __name__ + "] assume haploid in the following regions:\n")
    sys.stderr.write("chr\thap\tstart\tend\n")
    for reg in hap_regs:
        sys.stderr.write(reg.to_string() + "\n")
            
    # read CON file
    con_file = gzip.open(args[0], "rb") if args[0].endswith(".gz") else open(args[0], "rb")
    con_data = file_to_con_data(con_file)
    sys.stderr.write("[M::" + __name__ + "] read " + str(con_data.num_cons()) + " putative contacts (" +  str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    
    # apply regions
    con_data.apply_regs(inc_regs, exc_regs, hap_regs)
    sys.stderr.write("[M::" + __name__ + "] writing output for " + str(con_data.num_cons()) + " putative contacts (" + str(round(100.0 * con_data.num_intra_chr() / con_data.num_cons(), 2)) + "% intra-chromosomal, " + str(round(100.0 * con_data.num_phased_legs() / con_data.num_cons() / 2, 2)) + "% legs phased)\n")
    sys.stdout.write(con_data.to_string()+"\n")
    
    return 0
    
