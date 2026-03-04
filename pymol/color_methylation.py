"""
PyMOL script: color a dip-c CIF by methylation (CpG) stored in the B-factor.

Replicates:
    spectrumany b, magenta green, all, 0.005, 0.02

without requiring the spectrumany plugin, using PyMOL's built-in spectrum command.

Usage (PyMOL console):
    run color_methylation.py

Adjust PALETTE, MIN_VALUE, MAX_VALUE below as needed.
"""

from pymol import cmd

# Color gradient endpoints
PALETTE   = "magenta green"  # any two (or more) PyMOL color names
MIN_VALUE = 0.005             # B-factor mapped to the first color (magenta)
MAX_VALUE = 0.02              # B-factor mapped to the last color (green)


def color_methylation():
    cmd.spectrum("b", PALETTE, "all",
                 minimum=MIN_VALUE,
                 maximum=MAX_VALUE)
    print(f"Methylation coloring applied: {PALETTE}  [{MIN_VALUE}, {MAX_VALUE}]")


# visual settings
cmd.viewport(800, 800)
cmd.clip("slab", 10)
cmd.set("ray_shadows", 0)
cmd.set("ambient", 1)
cmd.set("specular", "off")
cmd.set("ray_opaque_background", "off")
cmd.show("sticks", "all")
cmd.set_bond("stick_radius", 0.5, "all")

color_methylation()
