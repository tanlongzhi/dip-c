"""
dip_c_color.py — unified PyMOL coloring commands for dip-c CIF files.

Load once per session:
    run /path/to/dip-c/pymol/dip_c_color.py

Available commands:
    dipc_view
    color_chromosomes
    color_methylation cpg_file [palette [min_val [max_val]]]
    color_centromere cen_file [palette [min_val [max_val]]]
    color_locus [chr_len_file]
    color_radial [palette [min_val [max_val]]]

Power-user utility (loads any chr/locus/value file into B-factors without coloring):
    load_bfactor value_file [missing_val]

Examples:
    dipc_view
    color_chromosomes
    color_methylation /path/to/hg19.cpg.20k.txt
    color_methylation /path/to/hg19.cpg.20k.txt, cyan yellow, 0.001, 0.05
    color_centromere /path/to/hg19.chr.cen
    color_centromere /path/to/hg19.chr.cen, blue white red
    color_locus
    color_locus /path/to/hg19.chr.len
    color_radial
    color_radial blue white red, 0, 500
"""

import pymol
from pymol import cmd


# ---------------------------------------------------------------------------
# Chromosome colors (UCSC-inspired)
# ---------------------------------------------------------------------------

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


def _hex_to_rgb(h):
    h = h.lstrip("#")
    return [int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4)]


def _darken(rgb, factor=0.65):
    return [c * factor for c in rgb]


def _lighten(rgb, factor=0.45):
    return [c + (1.0 - c) * factor for c in rgb]


def _locus_to_resn_name(locus):
    """Convert an integer locus to (resn, name, seq_id) as encoded by vis.py."""
    locus_string = str(locus).rjust(9, "0")
    return locus_string[0:3], locus_string[3:6], int(locus_string[6:9])


def _push_bfactors(b_by_id, missing_val=-1.0):
    """Write a {atom_ID: float} dict into B-factors in one alter call."""
    pymol.stored.b_by_id = b_by_id
    cmd.alter("all", "b = pymol.stored.b_by_id.get(ID, -1.0)")
    cmd.rebuild()


def _collect_atom_info():
    """Return list of (ID, chain, resn, name, resi) for all atoms via one iterate call."""
    pymol.stored.atom_info = []
    cmd.iterate("all", "pymol.stored.atom_info.append((ID, chain, resn, name, resi))")
    return pymol.stored.atom_info


# ---------------------------------------------------------------------------
# User-facing commands
# ---------------------------------------------------------------------------

def _reset_view():
    """Zoom to fit all atoms and set slab to the molecule's bounding box extent."""
    cmd.zoom("all", complete=1)
    extent = cmd.get_extent("all")
    if extent:
        size = max(extent[1][i] - extent[0][i] for i in range(3))
        cmd.clip("slab", max(size, 100))  # at least 100 Å so structure isn't invisible


def dipc_view():
    """Apply standard dip-c visualization settings (sticks, lighting, slab)."""
    cmd.viewport(800, 800)
    cmd.set("ray_shadows", 0)
    cmd.hide("everything", "all")
    cmd.show("sticks", "all")
    cmd.set_bond("stick_radius", 0.5, "all")
    _reset_view()
    print("dip-c view settings applied.")


def color_chromosomes():
    """Color each chromosome chain by haplotype: dark=maternal, light=paternal."""
    for chrom, hex_color in CHR_COLORS.items():
        base_rgb = _hex_to_rgb(hex_color)
        cmd.set_color(f"chr{chrom}_mat", _darken(base_rgb, factor=0.65))
        cmd.set_color(f"chr{chrom}_pat", _lighten(base_rgb, factor=0.45))
        cmd.color(f"chr{chrom}_mat", f'chain "{chrom}(mat)"')
        cmd.color(f"chr{chrom}_pat", f'chain "{chrom}(pat)"')
    _reset_view()
    print("Chromosome colors applied.  Dark=maternal  Light=paternal")


def color_locus(chr_len_file=""):
    """Color by locus position along each chromosome using each chromosome's identity color.

    Each chromosome grades white (start) → base color (midpoint) → near-black (end),
    so you can trace the fiber while still identifying chromosomes by color.
    Chromosome lengths are estimated from the maximum locus in the structure if no
    file is provided.

    Args:
        chr_len_file:  path to tab-delimited chr.len file (chr, length)  [optional]

    Examples:
        color_locus
        color_locus /path/to/hg19.chr.len
    """
    # --- chromosome lengths ---
    chr_lens = {}
    if chr_len_file.strip():
        with open(chr_len_file.strip()) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                chr_lens[parts[0]] = int(parts[1])
        print(f"Loaded chromosome lengths from {chr_len_file.strip()}")
    else:
        # estimate from max locus observed in structure
        for atom_id, chain, resn, name, resi in _collect_atom_info():
            chrom = chain.split("(")[0]
            try:
                locus = int(resn) * 1_000_000 + int(name) * 1_000 + int(resi)
            except ValueError:
                continue
            if chrom not in chr_lens or locus > chr_lens[chrom]:
                chr_lens[chrom] = locus
        print("Chromosome lengths estimated from structure (pass a chr.len file for accuracy).")

    # --- b = locus / chr_len (0.0 = start, 1.0 = end) ---
    b_by_id = {}
    for atom_id, chain, resn, name, resi in _collect_atom_info():
        chrom = chain.split("(")[0]
        if chrom not in chr_lens or chr_lens[chrom] == 0:
            continue
        try:
            locus = int(resn) * 1_000_000 + int(name) * 1_000 + int(resi)
        except ValueError:
            continue
        b_by_id[atom_id] = locus / chr_lens[chrom]

    pymol.stored.b_by_id = b_by_id
    cmd.alter("all", "b = -1.0")
    cmd.alter("all", "b = pymol.stored.b_by_id[ID] if ID in pymol.stored.b_by_id else b")
    cmd.rebuild()

    # --- per-chromosome 3-stop gradient: white → base color → near-black ---
    for chrom, hex_color in CHR_COLORS.items():
        base_rgb = _hex_to_rgb(hex_color)
        start_color = _lighten(base_rgb, factor=0.85)   # nearly white
        end_color   = _darken(base_rgb,  factor=0.12)   # nearly black
        cmd.set_color(f"chr{chrom}_locus_start", start_color)
        cmd.set_color(f"chr{chrom}_locus_mid",   base_rgb)
        cmd.set_color(f"chr{chrom}_locus_end",   end_color)
        palette = f"chr{chrom}_locus_start chr{chrom}_locus_mid chr{chrom}_locus_end"
        for hap in ("mat", "pat"):
            sel = f'chain "{chrom}({hap})" and b > -0.5'
            cmd.spectrum("b", palette, sel, minimum=0.0, maximum=1.0)

    _reset_view()
    print("Locus coloring applied.  White=start  →  chromosome color  →  black=end")


def color_methylation(cpg_file, palette="magenta green", min_val=0.005, max_val=0.02):
    """Load a CpG/methylation file and color by value.

    The file format is tab-delimited: chr, locus, value  (e.g. hg19.cpg.20k.txt)
    Values are applied to both (mat) and (pat) chains.

    Args:
        cpg_file:  path to tab-delimited file (chr, locus, value)
        palette:   PyMOL color gradient ["magenta green"]
        min_val:   value mapped to first color [0.005]
        max_val:   value mapped to last color  [0.02]
    """
    min_val = float(min_val)
    max_val = float(max_val)
    load_bfactor(cpg_file)
    cmd.clip("slab", 10)
    cmd.spectrum("b", palette, "all", minimum=min_val, maximum=max_val)
    print(f"Methylation coloring: {palette}  [{min_val}, {max_val}]")


def color_centromere(cen_file, palette="red white blue", min_val=0.0, max_val=1.0):
    """Color by arm position (0=centromere, 1=telomere) using a .cen file.

    The file format is tab-delimited: chr, length, centromere_position

    Args:
        cen_file:  path to tab-delimited .cen file
        palette:   PyMOL color gradient ["red white blue"]
        min_val:   value mapped to first color [0.0]
        max_val:   value mapped to last color  [1.0]
    """
    min_val = float(min_val)
    max_val = float(max_val)

    ref_lens = {}
    ref_cens = {}
    with open(cen_file.strip()) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            chrom, length, cen = parts[0], int(parts[1]), int(parts[2])
            ref_lens[chrom] = length
            ref_cens[chrom] = cen

    def arm_position(locus, chrom):
        if chrom not in ref_cens:
            return None
        cen = ref_cens[chrom]
        arm_locus = locus - cen
        arm_len = (ref_lens[chrom] - cen) if arm_locus > 0 else cen
        if arm_len == 0:
            return None
        return abs(arm_locus) / arm_len

    b_by_id = {}
    for atom_id, chain, resn, name, resi in _collect_atom_info():
        try:
            locus = int(resn) * 1_000_000 + int(name) * 1_000 + int(resi)
        except ValueError:
            continue  # skip unparseable atoms (matches color.py -L behavior)
        chrom = chain.split("(")[0]
        val = arm_position(locus, chrom)
        if val is None:
            continue  # skip chromosomes not in .cen file (matches color.py -L behavior)
        b_by_id[atom_id] = val

    # set all to -1.0 first, then only overwrite matched particles
    pymol.stored.b_by_id = b_by_id
    cmd.alter("all", "b = -1.0")
    cmd.alter("all", "b = pymol.stored.b_by_id[ID] if ID in pymol.stored.b_by_id else b")
    cmd.rebuild()
    _reset_view()
    # only spectrum matched particles (b > -0.5), leaving unmatched at -1.0 uncolored
    # PyMOL selectors don't support >=, so use > -0.5 as equivalent
    cmd.spectrum("b", palette, "b > -0.5", minimum=min_val, maximum=max_val)
    print(f"Centromere coloring: {palette}  [0.0=centromere, 1.0=telomere]")


def color_radial(palette="blue white red", min_val="auto", max_val="auto"):
    """Color by distance from the nuclear center of mass (matching dip-c color -C).

    No input file required — purely geometry-based on the loaded structure.
    0 = center, max = periphery.

    Args:
        palette:   PyMOL color gradient ["blue white red"]
        min_val:   distance mapped to first color ["auto" = 0]
        max_val:   distance mapped to last color  ["auto" = mean radial distance * 2]
    """
    import math

    # collect positions
    pymol.stored.pos_info = []
    cmd.iterate_state(1, "all", "pymol.stored.pos_info.append((ID, x, y, z))")
    pos_info = pymol.stored.pos_info
    if not pos_info:
        print("No atoms found.")
        return

    # center of mass (matching color.py -C: simple mean of all particle positions)
    n = len(pos_info)
    cx = sum(p[1] for p in pos_info) / n
    cy = sum(p[2] for p in pos_info) / n
    cz = sum(p[3] for p in pos_info) / n

    # distance for each atom
    b_by_id = {}
    for atom_id, x, y, z in pos_info:
        b_by_id[atom_id] = math.sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)

    distances = list(b_by_id.values())
    mean_r = sum(distances) / len(distances)
    actual_min = min(distances)
    actual_max = max(distances)

    lo = float(min_val) if min_val != "auto" else 0.0
    hi = float(max_val) if max_val != "auto" else actual_max

    pymol.stored.b_by_id = b_by_id
    cmd.alter("all", "b = pymol.stored.b_by_id.get(ID, 0.0)")
    cmd.rebuild()
    _reset_view()
    cmd.spectrum("b", palette, "all", minimum=lo, maximum=hi)
    print(f"Radial coloring: {palette}  center=({cx:.1f}, {cy:.1f}, {cz:.1f})")
    print(f"  distances: min={actual_min:.1f}  mean={mean_r:.1f}  max={actual_max:.1f}  spectrum=[{lo:.1f}, {hi:.1f}]")


# ---------------------------------------------------------------------------
# Power-user utility
# ---------------------------------------------------------------------------

def load_bfactor(value_file, missing_val=-1.0):
    """Load a chr/locus/value file into B-factors without applying any coloring.

    File format is tab-delimited: chr, locus, value  (same as dip-c color -c input)
    Values are applied to both (mat) and (pat) chains.
    After loading, run any spectrum command to visualize.

    Args:
        value_file:   path to tab-delimited file (chr, locus, value)
        missing_val:  B-factor for particles not in the file [-1.0]
    """
    missing_val = float(missing_val)

    color_map = {}
    with open(value_file.strip()) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            chrom, locus, value = parts[0], int(parts[1]), float(parts[2])
            resn, name, seq_id = _locus_to_resn_name(locus)
            color_map[(f"{chrom}(mat)", resn, name, seq_id)] = value
            color_map[(f"{chrom}(pat)", resn, name, seq_id)] = value
            # Backward compat: old CIFs had resi=1 hardcoded; new CIFs store
            # sub-kb digits so >=1kb loci produce seq_id=0.  Register both so
            # old and new CIFs match without collision (a real locus ending in
            # ...001 would have seq_id=1 and be registered directly above).
            if seq_id == 0:
                color_map[(f"{chrom}(mat)", resn, name, 1)] = value
                color_map[(f"{chrom}(pat)", resn, name, 1)] = value

    b_by_id = {}
    for atom_id, chain, resn, name, resi in _collect_atom_info():
        b_by_id[atom_id] = color_map.get((chain, resn, name, int(resi)), missing_val)

    _push_bfactors(b_by_id, missing_val)

    matched = sum(1 for v in b_by_id.values() if v != missing_val)
    print(f"Loaded {matched} / {len(b_by_id)} particles from {value_file}")


# ---------------------------------------------------------------------------
# Register commands
# ---------------------------------------------------------------------------

cmd.extend("dipc_view",          dipc_view)
cmd.extend("color_chromosomes",  color_chromosomes)
cmd.extend("color_methylation",  color_methylation)
cmd.extend("color_centromere",   color_centromere)
cmd.extend("color_locus",        color_locus)
cmd.extend("color_radial",       color_radial)
cmd.extend("load_bfactor",       load_bfactor)

# common typo aliases
def _typo_centromeres(*args, **kwargs):
    print("Did you mean: color_centromere  (no trailing 's')?")
    color_centromere(*args, **kwargs)

def _typo_chromosomes(*args, **kwargs):
    print("Did you mean: color_chromosomes?")
    color_chromosomes(*args, **kwargs)

cmd.extend("color_centromeres",  _typo_centromeres)
cmd.extend("colour_centromere",  _typo_centromeres)
cmd.extend("colour_chromosomes", _typo_chromosomes)

print("dip-c color commands loaded:")
print("  dipc_view")
print("  color_chromosomes")
print("  color_methylation cpg_file [palette [min_val [max_val]]]")
print("  color_centromere cen_file [palette [min_val [max_val]]]")
print("  color_locus [chr_len_file]")
print("  color_radial [palette [min_val [max_val]]]")
print("  load_bfactor value_file [missing_val]  (power-user)")

