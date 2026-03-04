"""
PyMOL script: color chromosomes by homolog (mat/pat) using UCSC-inspired colors.

Maternal (mat) = darkened version of each chromosome's base color
Paternal (pat) = lightened (pastel) version of the same color

Usage (PyMOL console):
    run color_chromosomes.py
"""

from pymol import cmd

# UCSC genome browser-inspired base colors per chromosome (hex RGB)
CHR_COLORS = {
    "1":  "#9966FF",  # violet
    "2":  "#6666CC",  # blue-violet
    "3":  "#3399FF",  # cornflower blue
    "4":  "#00CCFF",  # sky blue
    "5":  "#00CC99",  # teal-green
    "6":  "#00CC00",  # green
    "7":  "#66CC00",  # yellow-green
    "8":  "#CCCC00",  # yellow
    "9":  "#CC9900",  # dark yellow
    "10": "#CC6600",  # amber
    "11": "#CC3300",  # burnt orange
    "12": "#FF0000",  # red
    "13": "#CC0066",  # crimson
    "14": "#990099",  # purple
    "15": "#CC3399",  # pink-purple
    "16": "#FF66CC",  # pink
    "17": "#FF99AA",  # light pink
    "18": "#AABB88",  # sage
    "19": "#558855",  # forest green
    "20": "#88AACC",  # steel blue
    "21": "#AACCDD",  # powder blue
    "22": "#BBBBBB",  # gray
    "X":  "#FF1493",  # deep pink
    "Y":  "#1E90FF",  # dodger blue
}


def hex_to_rgb(h):
    h = h.lstrip("#")
    return [int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4)]


def darken(rgb, factor=0.65):
    """Darken an RGB triplet by scaling toward black."""
    return [c * factor for c in rgb]


def lighten(rgb, factor=0.45):
    """Lighten an RGB triplet by blending toward white."""
    return [c + (1.0 - c) * factor for c in rgb]


def color_chromosomes():
    for chrom, hex_color in CHR_COLORS.items():
        base_rgb = hex_to_rgb(hex_color)
        mat_rgb = darken(base_rgb, factor=0.65)
        pat_rgb = lighten(base_rgb, factor=0.45)

        mat_name = f"chr{chrom}_mat"
        pat_name = f"chr{chrom}_pat"

        cmd.set_color(mat_name, mat_rgb)
        cmd.set_color(pat_name, pat_rgb)

        cmd.color(mat_name, f'chain "{chrom}(mat)"')
        cmd.color(pat_name, f'chain "{chrom}(pat)"')

    print("Chromosome colors applied.  Dark = maternal,  Light = paternal")


color_chromosomes()
