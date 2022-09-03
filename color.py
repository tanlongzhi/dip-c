import sys
import getopt
from classes import Haplotypes, homologous_hom_name, LegData, ConData, file_to_con_data, Leg, Par, ParData, G3dData, file_to_g3d_data, string_to_leg
import math
import numpy as np


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

def same_haplotype_fraction(g3d_particle, nearby_g3d_particles):
    hom_name = g3d_particle.get_hom_name()
    haplotype = g3d_particle.get_haplotype()
    num_g3d_particles = -1 # to exclude self
    num_same_haplotype = -1 # to exclude self
    for nearby_g3d_particle in nearby_g3d_particles:
        num_g3d_particles += 1
        nearby_hom_name = nearby_g3d_particle.get_hom_name()
        nearby_haplotype = nearby_g3d_particle.get_haplotype()
        if nearby_haplotype == haplotype:
                num_same_haplotype += 1
    if num_g3d_particles == 0:
        return None
    return float(num_same_haplotype) / num_g3d_particles


def same_haplotype_fraction_excluding_same_homolog(g3d_particle, nearby_g3d_particles):
    hom_name = g3d_particle.get_hom_name()
    haplotype = g3d_particle.get_haplotype()
    num_g3d_particles = -1 # to exclude self
    num_same_haplotype = 0
    for nearby_g3d_particle in nearby_g3d_particles:
        nearby_hom_name = nearby_g3d_particle.get_hom_name()
        nearby_haplotype = nearby_g3d_particle.get_haplotype()
        if nearby_hom_name != hom_name:
            num_g3d_particles += 1
            if nearby_haplotype == haplotype:
                    num_same_haplotype += 1
    if num_g3d_particles == 0:
        return None
    return float(num_same_haplotype) / num_g3d_particles

def intra_hom_count(g3d_particle, nearby_g3d_particles, max_separation):
    hom_name = g3d_particle.get_hom_name()
    num_intra_hom = -1 # to exclude self
    for nearby_g3d_particle in nearby_g3d_particles:
        if nearby_g3d_particle.get_hom_name() == hom_name:
            if max_separation is None or abs(nearby_g3d_particle.get_ref_locus() - g3d_particle.get_ref_locus()) <= max_separation:
                num_intra_hom += 1
    return num_intra_hom

def hom_diversity(g3d_particles):
    particle_counts = {}
    total_count = 0
    for g3d_particle in g3d_particles:
        total_count += 1
        hom_name = g3d_particle.get_hom_name()
        if hom_name not in particle_counts:
            particle_counts[hom_name] = 0
        particle_counts[hom_name] += 1
    # calculate entropy
    entropy = 0.0
    for particle_count in particle_counts.values():
        hom_fraction = float(particle_count) / total_count
        entropy -= hom_fraction * math.log(hom_fraction)
    return entropy

def hom_richness(g3d_particles):
    hom_names = set()
    for g3d_particle in g3d_particles:
        hom_names.add(g3d_particle.get_hom_name())
    return len(hom_names)

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
    
def smooth_color_exc(g3d_particle, nearby_g3d_particles, color_data, g3d_resolution):
    hom_name = g3d_particle.get_hom_name()
    ref_locus = g3d_particle.get_ref_locus()
    ref_locus_up = ref_locus + g3d_resolution
    ref_locus_down = ref_locus - g3d_resolution
    num_g3d_particles = 0 
    sum_color = 0
    for nearby_g3d_particle in nearby_g3d_particles:
        nearby_hom_name = nearby_g3d_particle.get_hom_name()
        nearby_ref_locus = nearby_g3d_particle.get_ref_locus()

        if (nearby_hom_name == hom_name) and (nearby_ref_locus == ref_locus):
            continue
        elif (nearby_hom_name == hom_name) and (nearby_ref_locus == ref_locus_up):
            continue
        elif (nearby_hom_name == hom_name) and (nearby_ref_locus == ref_locus_down):
            continue
        elif (nearby_hom_name, nearby_ref_locus) in color_data:
            num_g3d_particles += 1
            sum_color += color_data[(nearby_hom_name, nearby_ref_locus)]
        
    if num_g3d_particles == 0:
        return None
    return sum_color / num_g3d_particles

def color(argv):
    # default parameters
    color_file_name = None
    color_mode = None
    max_distance = None
    smooth_distance = None
    smooth_distance_exc = None
    max_separation = None
    radial_mode = False
    radial_min_num_particles = 10
    radial_missing_value = -1.0
    radial_max_r = 3.0
    radial_bin_r = 0.05
    
    # display parameters
    disp_num_particles = 1000
    
    # read arguments
    try:
        opts, args = getopt.getopt(argv[1:], "c:n:l:m:L:i:s:S:hd:r:I:p:P:CD:R", ["min-num=", "missing=", "max-r=", "bin-size=", "c-hom=", "s-exc="])
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: dip-c color [options] <in.3dg>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -c <color.txt>    color by a list of locus-color pairs (tab-delimited: chr, locus, color)\n")
        sys.stderr.write("  --c-hom=<color.txt>    color by a list of homolog specific locus-color pairs (tab-delimited: homolog, locus, color)\n")
        sys.stderr.write("  -n <chr.txt>      color by chromosome name (one chromosome per line)\n")
        sys.stderr.write("  -l <chr.len>      color by locus divided by chromosome length (tab-delimited: chr, len)\n")
        sys.stderr.write("  -L <chr.cen>      color by arm locus divided by arm length (tab-delimited: chr, len, center of centromere)\n")
        sys.stderr.write("  -h                color by distance to homologous locus\n\n")
        sys.stderr.write("  -i FLOAT          color by percentage of intra-homologous neighbors within a given distance\n")
        sys.stderr.write("  -I FLOAT          color by number of intra-homologous neighbors within a given distance\n")
        sys.stderr.write("  -p FLOAT          color by percentage of same-haplotype neighbors within a given distance\n")
        sys.stderr.write("  -P FLOAT          color by percentage of same-haplotype neighbors (excluding the same homolog neighbors) within a given distance\n")
        sys.stderr.write("  -S INT            (with \"-i\" or \"-I\") max separation (bp) for intra-homologous neighbors\n\n")
        sys.stderr.write("  -d FLOAT          color by homolog diversity within a given distance\n")
        sys.stderr.write("  -r FLOAT          color by homolog richness within a given distance\n\n")
        sys.stderr.write("  -C                color by distance to the nuclear center of mass\n")
        sys.stderr.write("  -D <in.leg>       color by distance to a given locus (only the first line of the LEG file will be used)\n\n")
        sys.stderr.write("  -s FLOAT          smooth color by averaging over a ball\n\n")
        sys.stderr.write("  --s-exc=FLOAT     smooth color by averaging over a ball, excluding self and its two flanking particles\n")
        sys.stderr.write("  -R                special: output average color for different radial distances (normalized to 1.0)\n")
        sys.stderr.write("  --min-num=INT     (with \"-R\") min number of particles per bin [" + str(radial_min_num_particles) + "]\n")
        sys.stderr.write("  --missing=FLOAT   (with \"-R\") output value when \"--min-num\" is not met [" + str(radial_missing_value) + "]\n")
        sys.stderr.write("  --max-r=FLOAT     (with \"-R\") max radial distance [" + str(radial_max_r) + "]\n")
        sys.stderr.write("  --bin-size=FLOAT  (with \"-R\") bin size of radial distances [" + str(radial_bin_r) + "]\n\n")
        sys.stderr.write("Output:\n")
        sys.stderr.write("  tab-delimited: homolog, locus, color\n")
        sys.stderr.write("  (with \"-R\") tab-delimited: radial distance, average color, #particles\n")
        return 1
        
    num_color_schemes = 0
    for o, a in opts:
        if o == "-i" or o == "-I" or o == "-d" or o == "-r" or o == "-p" or o == "-P":
            num_color_schemes += 1
            color_mode = o[1:]
            max_distance = float(a)
        elif o == "-s":
            smooth_distance = float(a)
        elif o == "--s-exc":
            smooth_distance_exc = float(a)
        elif o == "-S":
            max_separation = int(a)
        elif o == "--min-num":
            radial_min_num_particles = int(a)
        elif o == "--missing":
            radial_missing_value = float(a)
        elif o == "--max-r":
            radial_max_r = float(a)
        elif o == "--bin-size":
            radial_bin_r = float(a)
        elif o == "-R":
            radial_mode = True
        elif o == "--c-hom":
            num_color_schemes += 1
            color_mode = o[2:]
            if a != "":
                color_file_name = a
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
    
    # prepare
    if color_mode is None:
        pass
    elif color_mode == "c":
        ref_name_ref_locus_colors = {}
        for color_file_line in color_file:
            ref_name, ref_locus, color = color_file_line.strip().split("\t")
            ref_locus = int(ref_locus)
            color = float(color)
            ref_name_ref_locus_colors[(ref_name, ref_locus)] = color
    elif color_mode == "c-hom":
        hom_name_ref_locus_colors = {}
        for color_file_line in color_file:
            hom_name, ref_locus, color = color_file_line.strip().split("\t")
            ref_locus = int(ref_locus)
            color = float(color)
            hom_name_ref_locus_colors[(hom_name, ref_locus)] = color
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
    elif color_mode == "i" or color_mode == "I" or color_mode == "d" or color_mode == "r" or color_mode == "p" or color_mode == "P":
        g3d_data.prepare_nearby()
    elif color_mode == "C":
        hom_names, loci_np_array, position_np_array = g3d_data.to_np_arrays()
        ref_pos = np.mean(position_np_array, axis = 0)
        sys.stderr.write("[M::" + __name__ + "] reference point (center of mass) is at (" + ", ".join(map(str, ref_pos)) + ")\n")
    elif color_mode == "D":
        # fine reference point position
        ref_leg = string_to_leg(color_file.readline().strip())
        g3d_data.prepare_interpolate()
        is_out, ref_pos = g3d_data.interpolate_leg(ref_leg)
        sys.stderr.write("[M::" + __name__ + "] reference point (" + ref_leg.to_string() + ") is at (" + ", ".join(map(str, ref_pos)) + ")\n")
                        
    # calculate colors for each particle
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
        elif color_mode == "c-hom":
            try:
                color = hom_name_ref_locus_colors[(g3d_particle.get_hom_name(), g3d_particle.get_ref_locus())]
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
            color = intra_hom_fraction(g3d_particle, g3d_data.get_g3d_particles_near(g3d_particle.get_position(), max_distance), max_separation)
            if color is None:
                continue
        elif color_mode == "I":
            color = intra_hom_count(g3d_particle, g3d_data.get_g3d_particles_near(g3d_particle.get_position(), max_distance), max_separation)
        elif color_mode == "p":
            color = same_haplotype_fraction(g3d_particle, g3d_data.get_g3d_particles_near(g3d_particle.get_position(), max_distance))
            if color is None:
                continue
        elif color_mode == "P":
            color = same_haplotype_fraction_excluding_same_homolog(g3d_particle, g3d_data.get_g3d_particles_near(g3d_particle.get_position(), max_distance))
            if color is None:
                continue
        elif color_mode == "h":
            homologous_g3d_particle = g3d_data.get_g3d_particle_from_hom_name_ref_locus(homologous_hom_name(g3d_particle.get_hom_name()), g3d_particle.get_ref_locus())
            if homologous_g3d_particle is None:
                continue
            color = math.sqrt((g3d_particle.get_x() - homologous_g3d_particle.get_x()) ** 2 + (g3d_particle.get_y() - homologous_g3d_particle.get_y()) ** 2 + (g3d_particle.get_z() - homologous_g3d_particle.get_z()) ** 2)
        elif color_mode == "d":
            color = hom_diversity(g3d_data.get_g3d_particles_near(g3d_particle.get_position(), max_distance))
        elif color_mode == "r":
            color = hom_richness(g3d_data.get_g3d_particles_near(g3d_particle.get_position(), max_distance))
        elif color_mode == "C" or color_mode == "D":
            color = math.sqrt((g3d_particle.get_x() - ref_pos[0]) ** 2 + (g3d_particle.get_y() - ref_pos[1]) ** 2 + (g3d_particle.get_z() - ref_pos[2]) ** 2)
        #sys.stderr.write(str(color) + "\n")
        color_data[g3d_particle.get_hom_name(), g3d_particle.get_ref_locus()] = color
        
    # smoothing
    if not smooth_distance is None:
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
        
    if not smooth_distance_exc is None:
        g3d_data.prepare_nearby()
        smooth_color_data = {}
        atom_id = 0
        for g3d_particle in g3d_data.get_g3d_particles():
            atom_id += 1
            if atom_id % disp_num_particles == 0:
                sys.stderr.write("[M::" + __name__ + "] smoothed " + str(atom_id) + " particles (" + str(round(100.0 * atom_id / g3d_data.num_g3d_particles(), 2)) + "%)\n")
            color = smooth_color_exc(g3d_particle, g3d_data.get_g3d_particles_near(g3d_particle.get_position(), smooth_distance_exc), color_data, g3d_resolution)
            if not color is None:
                smooth_color_data[g3d_particle.get_hom_name(), g3d_particle.get_ref_locus()] = color
        color_data = smooth_color_data
        
    # radial
    if radial_mode:
        num_radial_bins = int(radial_max_r / radial_bin_r) + 1
        radial_color_sums = [0.0] * num_radial_bins
        radial_color_nums = [0] * num_radial_bins
        
        # calculate center of mass, and normalization factor
        hom_names, loci_np_array, position_np_array = g3d_data.to_np_arrays()
        ref_pos = np.mean(position_np_array, axis = 0)
        mean_radial = np.mean(np.sum((position_np_array - ref_pos) ** 2, axis = -1) ** 0.5, axis = 0)
        sys.stderr.write("[M::" + __name__ + "] radial mode: average radial distance = " + str(mean_radial) + ", which will be normalize to 1.0\n")
        
        # examine each particle
        for g3d_particle in g3d_data.get_g3d_particles():
            atom_id += 1
            if atom_id % disp_num_particles == 0:
                sys.stderr.write("[M::" + __name__ + "] radial mode for " + str(atom_id) + " particles (" + str(round(100.0 * atom_id / g3d_data.num_g3d_particles(), 2)) + "%)\n")
            if (g3d_particle.get_hom_name(), g3d_particle.get_ref_locus()) not in color_data:
                continue
            color = color_data[g3d_particle.get_hom_name(), g3d_particle.get_ref_locus()]
            radial = math.sqrt((g3d_particle.get_x() - ref_pos[0]) ** 2 + (g3d_particle.get_y() - ref_pos[1]) ** 2 + (g3d_particle.get_z() - ref_pos[2]) ** 2) / mean_radial
            radial_bin_id = int(radial / radial_bin_r + 0.5)
            #sys.stderr.write(str(radial)+", " + str(radial_bin_id) + "=" + str(radial_bin_id*radial_bin_r)+ ", "+ str(color)+"\n")
            if radial_bin_id >= num_radial_bins:
                continue # out of bound, skip
            radial_color_sums[radial_bin_id] += color
            radial_color_nums[radial_bin_id] += 1
        
        # output
        sys.stderr.write("[M::" + __name__ + "] writing radial mode output\n")
        for radial_bin_id in range(num_radial_bins):
            if radial_color_nums[radial_bin_id] < radial_min_num_particles:
                output_value = radial_missing_value
            else:
                output_value = radial_color_sums[radial_bin_id] / radial_color_nums[radial_bin_id]
            sys.stdout.write("\t".join([str(radial_bin_id * radial_bin_r), str(output_value), str(radial_color_nums[radial_bin_id])]) + "\n")
    
        return 0        
            
    # output
    sys.stderr.write("[M::" + __name__ + "] writing " + str(len(color_data)) + " colors (" + str(round(100.0 * len(color_data) / g3d_data.num_g3d_particles(), 2)) + "%)\n")
    for hom_name, ref_locus in sorted(color_data.keys()):
        sys.stdout.write("\t".join([hom_name, str(ref_locus), str(color_data[(hom_name, ref_locus)])]) + "\n")
    
    return 0
