import sys
import getopt
import gzip
import copy
from classes import Haplotypes, LegData, ConData, file_to_con_data, Leg, Par, ParData, G3dData, file_to_g3d_data
import numpy as np
import math
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
    
    # display parameters
    display_quantiles = np.arange(0.0, 1.01, 0.01)
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "c:d:q:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac vis [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -c <chr.txt>      color by chromosome name (one chromosome per line)\n")
        sys.stderr.write("  -l <chr.len>      color by locus along each chromosome (tab-delimited: chr, len)\n\n")
        sys.stderr.write("mmCIF format:\n")
        sys.stderr.write("  label_asym_id     homolog name (e.g. \"1(mat)\")\n")
        sys.stderr.write("  label_comp_id     locus // 1 Mb, 3 digits with leading zeros\n")
        sys.stderr.write("  label_seq_id      1\n")
        sys.stderr.write("  label_atom_id     locus % 1 Mb // 1 kb, 3 digits with leading zeros\n")
        sys.stderr.write("  B_iso_or_equiv    scalar color\n")
        sys.stderr.write("  covale            backbone bond\n")
        return 1
    for o, a in opts:
        if o == "-c":
            color_file_name = a
        elif o == "-l":
            max_clean_distance = int(a)
            
    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + str(g3d_resolution) + " bp resolution\n")

    # open color file
    color_file = open(color_file_name, "rb")
    
    # -c: read chromosome name file
    ref_name_colors = {}
    color_counter = 0
    for color_file_line in color_file:
        color_counter += 1
        ref_name = color_file_line.strip()
        ref_name_colors[ref_name] = color_counter
        
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
        color = ref_name_colors[g3d_particle.get_ref_name()]
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
    