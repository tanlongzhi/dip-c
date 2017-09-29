import sys
import getopt
from classes import Haplotypes, LegData, ConData, file_to_con_data, Leg, Par, ParData, G3dData, file_to_g3d_data
from pdbx.reader.PdbxReader  import PdbxReader
from pdbx.writer.PdbxWriter  import PdbxWriter
from pdbx.reader.PdbxContainers import *

def g3d_particle_to_atom_data(g3d_particle, atom_id, color):
    locus_string = str(g3d_particle.get_ref_locus()).rjust(9,'0')
    return ("HETATM", atom_id, g3d_particle.hom_name(), locus_string[0:3], 1, locus_string[3:6], g3d_particle.get_x(), g3d_particle.get_y(), g3d_particle.get_z(), color)

def g3d_particle_tuple_to_conn_data(g3d_particle_tuple):
    locus_1_string = str(g3d_particle_tuple[0].get_ref_locus()).rjust(9,'0')
    locus_2_string = str(g3d_particle_tuple[1].get_ref_locus()).rjust(9,'0')
    return ("covale", g3d_particle_tuple[0].hom_name(), locus_1_string[0:3], 1, locus_1_string[3:6], g3d_particle_tuple[1].hom_name(), locus_2_string[0:3], 1, locus_2_string[3:6])

def vis(argv):
    # default parameters
    color_file_name = None
    color_mode = None
    missing_value = -1.0
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "c:n:l:m:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac vis [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -c <color.txt>    color by a list of locus-color pairs (tab-delimited: chr, locus, color)\n")
        sys.stderr.write("  -n <chr.txt>      color by chromosome name (one chromosome per line)\n")
        sys.stderr.write("  -l <chr.len>      color by locus divided by chromosome length (tab-delimited: chr, len)\n\n")
        sys.stderr.write("  -m FLOAT          color for particles that are missing from the color scheme [" + str(missing_value) + "]\n\n")
        sys.stderr.write("Output mmCIF format:\n")
        sys.stderr.write("  label_asym_id     homolog name (e.g. \"1(mat)\")\n")
        sys.stderr.write("  label_comp_id     locus // 1 Mb, 3 digits with leading zeros\n")
        sys.stderr.write("  label_seq_id      1\n")
        sys.stderr.write("  label_atom_id     locus % 1 Mb // 1 kb, 3 digits with leading zeros\n")
        sys.stderr.write("  B_iso_or_equiv    scalar color\n")
        sys.stderr.write("  covale            backbone bond\n")
        return 1
        
    num_color_schemes = 0
    for o, a in opts:
        if o == "-m":
            missing_value = float(a)
        else:
            num_color_schemes += 1
            color_mode = o[1:]
            color_file_name = a
    if num_color_schemes > 1:
        sys.stderr.write("[E::" + __name__ + "] only one color scheme is allowed\n")
        return 1
                    
    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + str(g3d_resolution) + " bp resolution\n")

    # open color file
    if not color_file_name is None:
        color_file = open(color_file_name, "rb")
    
    # read color file
    if color_mode is None:
        pass
    elif color_mode == "c":
        ref_name_ref_locus_colors = {}
        for color_file_line in color_file:
            ref_name, ref_locus, color = color_file_line.strip().split("\t")
            ref_locus = int(ref_locus)
            color = float(color)
            ref_name_ref_locus_colors[(ref_name, ref_locus)] = color
    elif color_mode == "n":
        ref_name_colors = {}
        color_counter = 0
        for color_file_line in color_file:
            color_counter += 1
            ref_name = color_file_line.strip()
            ref_name_colors[ref_name] = color_counter
    elif color_mode == "l":
        ref_lens = {}
        for color_file_line in color_file:
            ref_name, ref_len = color_file_line.strip().split("\t")
            ref_len = int(ref_len)
            ref_lens[ref_name] = ref_len
        
    # open mmCIF file to write
    myDataList = []
    curContainer = DataContainer("myblock")
    aCat = DataCategory("atom_site")
    aCat.appendAttribute("group_PDB")
    aCat.appendAttribute("id")
    aCat.appendAttribute("label_asym_id")
    aCat.appendAttribute("label_comp_id")
    aCat.appendAttribute("label_seq_id")  
    aCat.appendAttribute("label_atom_id")  
    aCat.appendAttribute("Cartn_x")  
    aCat.appendAttribute("Cartn_y")  
    aCat.appendAttribute("Cartn_z")
    aCat.appendAttribute("B_iso_or_equiv")

    sCat = DataCategory("struct_conn")
    sCat.appendAttribute("conn_type_id")
    sCat.appendAttribute("ptnr1_label_asym_id")
    sCat.appendAttribute("ptnr1_label_comp_id")
    sCat.appendAttribute("ptnr1_label_seq_id")
    sCat.appendAttribute("ptnr1_label_atom_id")
    sCat.appendAttribute("ptnr2_label_asym_id")
    sCat.appendAttribute("ptnr2_label_comp_id")
    sCat.appendAttribute("ptnr2_label_seq_id")
    sCat.appendAttribute("ptnr2_label_atom_id")
    
    # write atoms
    atom_id = 0
    for g3d_particle in g3d_data.get_g3d_particles():
        atom_id += 1
        
        # color
        if color_mode is None:
            color = atom_id
        elif color_mode == "c":
            try:
                color = ref_name_ref_locus_colors[(g3d_particle.get_ref_name(), g3d_particle.get_ref_locus())]
            except KeyError:
                color = missing_value
        elif color_mode == "n":
            try:
                color = ref_name_colors[g3d_particle.get_ref_name()]
            except KeyError:
                color = missing_value
        elif color_mode == "l":
            try:
                color = float(g3d_particle.get_ref_locus()) / ref_lens[g3d_particle.get_ref_name()]
            except KeyError:
                color = missing_value        
        aCat.append(g3d_particle_to_atom_data(g3d_particle, atom_id, color))
    
    # write backbond bonds
    for g3d_particle_tuple in g3d_data.get_adjacent_g3d_particle_tuples(g3d_resolution):
        sCat.append(g3d_particle_tuple_to_conn_data(g3d_particle_tuple))
 
    # write output
    curContainer.append(sCat)
    curContainer.append(aCat)
    myDataList.append(curContainer)
    pdbxW = PdbxWriter(sys.stdout)
    pdbxW.write(myDataList)    
    
    return 0
    