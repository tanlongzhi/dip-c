"""
PyMOL script: color a loaded dip-c CIF by centromere-arm position.

Replicates dip-c color -L logic entirely in PyMOL — no new CIF needed.
Each particle is colored by: abs(locus - centromere) / arm_length
  0.0 = centromere, 1.0 = telomere

Usage (PyMOL console):
    run color_centromere.py

Edit CEN_FILE to point to your .cen file.
"""

import sys
import pymol
from pymol import cmd

# Pass the .cen file path in one of two ways:
#
# 1. Command line (headless):
#      pymol -qc color_centromere.py -- color/hg19.chr.cen
#
# 2. PyMOL console (interactive):
#      cmd.stored.cen_file = "/Users/darrin/git/dip-c/color/hg19.chr.cen"
#      run /Users/darrin/git/dip-c/pymol/color_centromere.py

DEFAULT_CEN_FILE = "color/hg19.chr.cen"

args = [a for a in sys.argv if not a.startswith("-")]
if len(args) > 1:
    CEN_FILE = args[1]                              # command-line mode
else:
    CEN_FILE = getattr(pymol.stored, "cen_file", DEFAULT_CEN_FILE)  # console mode

# Color gradient: centromere -> telomere
PALETTE   = "red white blue"
MIN_VALUE = 0.0
MAX_VALUE = 1.0


def read_cen_file(path):
    ref_lens = {}
    ref_cens = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            chrom, length, cen = parts[0], int(parts[1]), int(parts[2])
            ref_lens[chrom] = length
            ref_cens[chrom] = cen
    return ref_lens, ref_cens


def arm_position(locus, chrom, ref_lens, ref_cens):
    """abs(locus - centromere) / arm_length, matching dip-c -L logic."""
    if chrom not in ref_cens:
        return None
    cen = ref_cens[chrom]
    arm_locus = locus - cen
    arm_len = (ref_lens[chrom] - cen) if arm_locus > 0 else cen
    if arm_len == 0:
        return None
    return abs(arm_locus) / arm_len


def color_centromere():
    ref_lens, ref_cens = read_cen_file(CEN_FILE)

    # Step 1: collect atom IDs and their chain/resn/name via a single iterate call
    pymol.stored.atom_info = []
    cmd.iterate("all", "pymol.stored.atom_info.append((ID, chain, resn, name))")

    # Step 2: compute b-value for each atom ID
    b_by_id = {}
    for atom_id, chain, resn, name in pymol.stored.atom_info:
        try:
            locus = int(resn) * 1_000_000 + int(name) * 1_000
        except ValueError:
            b_by_id[atom_id] = -1.0
            continue
        chrom = chain.split("(")[0]
        val = arm_position(locus, chrom, ref_lens, ref_cens)
        b_by_id[atom_id] = val if val is not None else -1.0

    # Step 3: push all b-values in a single alter call keyed by atom ID
    pymol.stored.b_by_id = b_by_id
    cmd.alter("all", "b = pymol.stored.b_by_id.get(ID, -1.0)")

    cmd.rebuild()
    cmd.spectrum("b", PALETTE, "all", minimum=MIN_VALUE, maximum=MAX_VALUE)
    print(f"Centromere coloring: {PALETTE}  [0.0=centromere, 1.0=telomere]")


color_centromere()
