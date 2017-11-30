import sys
import getopt
from classes import Haplotypes, homologous_hom_name, LegData, ConData, file_to_con_data, Leg, Par, ParData, G3dData, file_to_g3d_data
import math

def intra_hom_fraction(g3d_particle, nearby_g3d_particles, max_separation):
    hom_name = g3d_particle.get_hom_name()
    num_g3d_particles = -1 # to exclude self
    num_intra_hom = -1 # to exclude self
    for nearby_g3d_particle in nearby_g3d_particles:
        num_g3d_particles += 1
        if nearby_g3d_particle.get_hom_name() == hom_name:
            if max_separation is None or abs(nearby_g3d_particle.get_ref_locus() - g3d_particle.get_ref_locus()) <= max_separation:
                num_intra_hom += 1
    if num_g3d_particles == 0:
        return None
    return float(num_intra_hom) / num_g3d_particles

def smooth_color(g3d_particle, nearby_g3d_particles, color_data):
    hom_name = g3d_particle.get_hom_name()
    ref_locus = g3d_particle.get_ref_locus()
    num_g3d_particles = 0 # include self
    sum_color = 0
    for nearby_g3d_particle in nearby_g3d_particles:
        nearby_hom_name = nearby_g3d_particle.get_hom_name()
        nearby_ref_locus = nearby_g3d_particle.get_ref_locus()
        if (nearby_hom_name, nearby_ref_locus) in color_data:
            num_g3d_particles += 1
            sum_color += color_data[(nearby_hom_name, nearby_ref_locus)]
    
    if num_g3d_particles == 0:
        return None
    return sum_color / num_g3d_particles
    
def color(argv):
    # default parameters
    color_file_name = None
    color_mode = None
    intra_distance = None
    smooth_distance = None
    max_separation = None
    
    # display parameters
    disp_num_particles = 1000
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "c:n:l:m:L:i:s:S:h")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: metac color [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -c <color.txt>    color by a list of locus-color pairs (tab-delimited: chr, locus, color)\n")
        sys.stderr.write("  -n <chr.txt>      color by chromosome name (one chromosome per line)\n")
        sys.stderr.write("  -l <chr.len>      color by locus divided by chromosome length (tab-delimited: chr, len)\n")
        sys.stderr.write("  -L <chr.cen>      color by arm locus divided by arm length (tab-delimited: chr, len, center of centromere)\n")
        sys.stderr.write("  -h                color by distance to homologous locus\n\n")
        sys.stderr.write("  -i FLOAT          color by percentage of intra-homologous neighbors within a given distance\n")
        sys.stderr.write("  -S INT            (with \"-i\") max separation (bp) for intra-homologous neighbors\n\n")
        sys.stderr.write("  -s FLOAT          smooth color by averaging over a ball\n")
        sys.stderr.write("Output:\n")
        sys.stderr.write("  tab-delimited: homolog, locus, color\n")
        return 1
        
    num_color_schemes = 0
    for o, a in opts:
        if o == "-i":
            num_color_schemes += 1
            color_mode = "i"
            intra_distance = float(a)
        elif o == "-s":
            smooth_distance = float(a)
        elif o == "-S":
            max_separation = int(a)
        else:
            num_color_schemes += 1
            color_mode = o[1:]
            if a != "":
                color_file_name = a
    if not max_separation is None and color_mode != "i":
        sys.stderr.write("[E::" + __name__ + "] \"-S\" must be used with \"-i\"\n")
        return 1
    if num_color_schemes != 1:
        sys.stderr.write("[E::" + __name__ + "] exactly one color scheme is needed\n")
        return 1
                    
    # read 3DG file
    g3d_data = file_to_g3d_data(open(args[0], "rb"))
    g3d_data.sort_g3d_particles()
    g3d_resolution = g3d_data.resolution()
    sys.stderr.write("[M::" + __name__ + "] read a 3D structure with " + str(g3d_data.num_g3d_particles()) + " particles at " + ("N.A." if g3d_resolution is None else str(g3d_resolution)) + " bp resolution\n")

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
    elif color_mode == "L":
        ref_lens = {}
        ref_cens = {}
        for color_file_line in color_file:
            ref_name, ref_len, ref_cen = color_file_line.strip().split("\t")
            ref_len = int(ref_len)
            ref_cen = int(ref_cen)
            ref_lens[ref_name] = ref_len
            ref_cens[ref_name] = ref_cen
    elif color_mode == "i":
        sys.stderr.write("[M::" + __name__ + "] preparing KD tree\n")
        g3d_data.prepare_nearby()
                        
    # calculate colors
    color_data = {}
    atom_id = 0
    for g3d_particle in g3d_data.get_g3d_particles():
        atom_id += 1
        if atom_id % disp_num_particles == 0:
            sys.stderr.write("[M::" + __name__ + "] analyzed " + str(atom_id) + " particles (" + str(round(100.0 * atom_id / g3d_data.num_g3d_particles(), 2)) + "%)\n")
        
        # color
        if color_mode == "c":
            try:
                color = ref_name_ref_locus_colors[(g3d_particle.get_ref_name(), g3d_particle.get_ref_locus())]
            except KeyError:
                continue
        elif color_mode == "n":
            try:
                color = ref_name_colors[g3d_particle.get_ref_name()]
            except KeyError:
                continue
        elif color_mode == "l":
            try:
                color = float(g3d_particle.get_ref_locus()) / ref_lens[g3d_particle.get_ref_name()]
            except KeyError:
                continue       
        elif color_mode == "L":
            try:
                arm_locus = g3d_particle.get_ref_locus() - ref_cens[g3d_particle.get_ref_name()]
                if arm_locus > 0:
                    arm_len = ref_lens[g3d_particle.get_ref_name()] - ref_cens[g3d_particle.get_ref_name()]
                else:
                    arm_len = ref_cens[g3d_particle.get_ref_name()]
                color = float(abs(arm_locus)) / arm_len
            except KeyError:
                continue    
        elif color_mode == "i":
            color = intra_hom_fraction(g3d_particle, g3d_data.get_g3d_particles_near(g3d_particle.get_position(), intra_distance), max_separation)
            if color is None:
                continue
        elif color_mode == "h":
            homologous_g3d_particle = g3d_data.get_g3d_particle_from_hom_name_ref_locus(homologous_hom_name(g3d_particle.get_hom_name()), g3d_particle.get_ref_locus())
            if homologous_g3d_particle is None:
                continue
            color = math.sqrt((g3d_particle.get_x() - homologous_g3d_particle.get_x()) ** 2 + (g3d_particle.get_y() - homologous_g3d_particle.get_y()) ** 2 + (g3d_particle.get_z() - homologous_g3d_particle.get_z()) ** 2)
        
        color_data[g3d_particle.get_hom_name(), g3d_particle.get_ref_locus()] = color
        
    # smoothing
    if not smooth_distance is None:
        sys.stderr.write("[M::" + __name__ + "] preparing KD tree\n")
        g3d_data.prepare_nearby()
        smooth_color_data = {}
        atom_id = 0
        for g3d_particle in g3d_data.get_g3d_particles():
            atom_id += 1
            if atom_id % disp_num_particles == 0:
                sys.stderr.write("[M::" + __name__ + "] smoothed " + str(atom_id) + " particles (" + str(round(100.0 * atom_id / g3d_data.num_g3d_particles(), 2)) + "%)\n")
            color = smooth_color(g3d_particle, g3d_data.get_g3d_particles_near(g3d_particle.get_position(), smooth_distance), color_data)
            if not color is None:
                smooth_color_data[g3d_particle.get_hom_name(), g3d_particle.get_ref_locus()] = color
        color_data = smooth_color_data
            
    # output
    sys.stderr.write("[M::" + __name__ + "] writing " + str(len(color_data)) + " colors (" + str(round(100.0 * len(color_data) / g3d_data.num_g3d_particles(), 2)) + "%)\n")
    for hom_name, ref_locus in sorted(color_data.keys()):
        sys.stdout.write("\t".join([hom_name, str(ref_locus), str(color_data[(hom_name, ref_locus)])]) + "\n")
    
    return 0
    