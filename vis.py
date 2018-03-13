import sys
import getopt
from classes import Haplotypes, LegData, ConData, file_to_con_data, Leg, Par, ParData, G3dData, file_to_g3d_data
from pdbx.reader.PdbxReader  import PdbxReader
from pdbx.writer.PdbxWriter  import PdbxWriter
from pdbx.reader.PdbxContainers import *

def g3d_particle_to_atom_data(g3d_particle, atom_id, color):
    locus_string = str(g3d_particle.get_ref_locus()).rjust(9,'0')
    return ("HETATM", ".", atom_id, g3d_particle.get_hom_name(), locus_string[0:3], 1, locus_string[3:6], g3d_particle.get_x(), g3d_particle.get_y(), g3d_particle.get_z(), color)

def g3d_particle_tuple_to_conn_data(g3d_particle_tuple, conn_id):
    locus_1_string = str(g3d_particle_tuple[0].get_ref_locus()).rjust(9,'0')
    locus_2_string = str(g3d_particle_tuple[1].get_ref_locus()).rjust(9,'0')
    return (conn_id, "covale", g3d_particle_tuple[0].get_hom_name(), locus_1_string[0:3], 1, locus_1_string[3:6], g3d_particle_tuple[1].get_hom_name(), locus_2_string[0:3], 1, locus_2_string[3:6])

def vis(argv):
    # default parameters
    color_file_name = None
    missing_value = -1.0
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "c:m:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c vis [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -c <color.txt>    color by a list of locus-color pairs (tab-delimited: homolog, locus, color)\n")
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
        elif o == "-c":
            color_file_name = a
                    
    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + ("N.A." if g3d_resolution is None else str(g3d_resolution)) + " bp resolution\n")

    # read color file
    color_data = {}
    if not color_file_name is None:
        color_file = open(color_file_name, "rb")
        for color_file_line in color_file:
            hom_name, ref_locus, color = color_file_line.strip().split("\t")
            ref_locus = int(ref_locus)
            color = float(color)
            color_data[(hom_name, ref_locus)] = color
                    
    # open mmCIF file to write
    myDataList = []
    curContainer = DataContainer("myblock")
    aCat = DataCategory("atom_site")
    aCat.appendAttribute("group_PDB")
    aCat.appendAttribute("type_symbol")
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
    sCat.appendAttribute("id")
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
        try:
            color = color_data[(g3d_particle.get_hom_name(), g3d_particle.get_ref_locus())]
        except KeyError:
            color = missing_value    
        aCat.append(g3d_particle_to_atom_data(g3d_particle, atom_id, color))
    
    # write backbond bonds
    conn_id = 0
    for g3d_particle_tuple in g3d_data.get_adjacent_g3d_particle_tuples(g3d_resolution):
        conn_id += 1
        sCat.append(g3d_particle_tuple_to_conn_data(g3d_particle_tuple, conn_id))
 
    # write output
    curContainer.append(sCat)
    curContainer.append(aCat)
    myDataList.append(curContainer)
    pdbxW = PdbxWriter(sys.stdout)
    pdbxW.write(myDataList)    
    
    return 0
    
