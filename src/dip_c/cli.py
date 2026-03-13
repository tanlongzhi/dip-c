#!/usr/bin/env python3

import os
import sys
import time


def main():
    if len(sys.argv) == 1:
        sys.stderr.write("Usage: dip-c <command> <arguments>\n")
        sys.stderr.write("Commands:\n")
        sys.stderr.write("  seg      extract read segments (SEG) from a BAM file\n")
        sys.stderr.write("  con      extract contacts (CON) from a SEG file\n")
        sys.stderr.write("  dedup    merge duplicates in a CON file\n")
        sys.stderr.write("  reg      exclude genomic regions from a CON file\n")
        sys.stderr.write("  bincon   bin a CON file into a matrix\n")
        sys.stderr.write("  clean    remove artifacts from a CON file\n")
        sys.stderr.write("  impute   impute missing haplotypes in a CON file\n")
        sys.stderr.write("  impute3  impute missing haplotypes in a CON file based on a 3DG file\n")
        sys.stderr.write("  clean3   remove contact-poor particles (e.g. centromeres) from a 3DG file\n")
        sys.stderr.write("  reg3     exclude genomic regions from a 3DG file\n")
        sys.stderr.write("  mkcon    make artificial CON from contact legs\n")
        sys.stderr.write("\n")
        sys.stderr.write("  color    color particles in a 3DG file\n")
        sys.stderr.write("  color2   color particles based on a CON file\n")
        sys.stderr.write("  mgcolor  merge color files\n")
        sys.stderr.write("  vis      visualize a 3DG file into mmCIF format, and color particles\n")
        sys.stderr.write("\n")
        sys.stderr.write("  info     basic information about CON file(s)\n")
        sys.stderr.write("  ard      contacts around certain points (e.g. other contacts)\n")
        sys.stderr.write("  cv       cross validation of CON files\n")
        sys.stderr.write("  align    align 3DG files\n")
        sys.stderr.write("  dist     calculate 3D distances from a 3DG file\n")
        sys.stderr.write("  pd       calculate 3D pairwise distances from a 3DG file\n")
        sys.stderr.write("  rg       calculate the Rg matrix from a 3DG file\n")
        sys.stderr.write("  con3     generate CON from a 3DG file\n")
        sys.stderr.write("  pos      find 3D positions from a 3DG file\n")
        sys.stderr.write("  tad      find the TAD tree from a Rg matrix\n")
        sys.stderr.write("  exp      expand a 3DG file by translating each chromosome\n")
        sys.stderr.write("  force    generate a 3DG file from a CON file using force-directed layout\n")
        sys.stderr.write("\n")
        sys.stderr.write("  pairs2con  convert hickit .pairs.gz to .con.gz format\n")
        sys.stderr.write("\n")
        sys.stderr.write("  data-path  print path to installed data directory\n")
        return 1

    if sys.argv[1] == "--completion":
        from dip_c.completion import generate_completion_script
        print(generate_completion_script(), end="")
        return 0

    start_time = time.time()
    command = sys.argv[1]

    if command == "seg":
        from dip_c.commands.seg import seg
        return_value = seg(sys.argv[1:])
    elif command == "con":
        from dip_c.commands.con import con
        return_value = con(sys.argv[1:])
    elif command == "dedup":
        from dip_c.commands.dedup import dedup
        return_value = dedup(sys.argv[1:])
    elif command == "reg":
        from dip_c.commands.reg import reg
        return_value = reg(sys.argv[1:])
    elif command == "clean":
        from dip_c.commands.clean import clean
        return_value = clean(sys.argv[1:])
    elif command == "bincon":
        from dip_c.commands.bincon import bincon
        return_value = bincon(sys.argv[1:])
    elif command == "impute":
        from dip_c.commands.impute import impute
        return_value = impute(sys.argv[1:])
    elif command == "impute3":
        from dip_c.commands.impute3 import impute3
        return_value = impute3(sys.argv[1:])
    elif command == "clean3":
        from dip_c.commands.clean3 import clean3
        return_value = clean3(sys.argv[1:])
    elif command == "reg3":
        from dip_c.commands.reg3 import reg3
        return_value = reg3(sys.argv[1:])
    elif command == "mkcon":
        from dip_c.commands.mkcon import mkcon
        return_value = mkcon(sys.argv[1:])
    elif command == "color":
        from dip_c.commands.color import color
        return_value = color(sys.argv[1:])
    elif command == "color2":
        from dip_c.commands.color2 import color2
        return_value = color2(sys.argv[1:])
    elif command == "mgcolor":
        from dip_c.commands.mgcolor import mgcolor
        return_value = mgcolor(sys.argv[1:])
    elif command == "vis":
        from dip_c.commands.vis import vis
        return_value = vis(sys.argv[1:])
    elif command == "info":
        from dip_c.commands.info import info
        return_value = info(sys.argv[1:])
    elif command == "ard":
        from dip_c.commands.ard import ard
        return_value = ard(sys.argv[1:])
    elif command == "cv":
        from dip_c.commands.cv import cv
        return_value = cv(sys.argv[1:])
    elif command == "con3":
        from dip_c.commands.con3 import con3
        return_value = con3(sys.argv[1:])
    elif command == "dist":
        from dip_c.commands.dist import dist
        return_value = dist(sys.argv[1:])
    elif command == "align":
        from dip_c.commands.align import align
        return_value = align(sys.argv[1:])
    elif command == "pd":
        from dip_c.commands.pd import pd
        return_value = pd(sys.argv[1:])
    elif command == "rg":
        from dip_c.commands.rg import rg
        return_value = rg(sys.argv[1:])
    elif command == "pos":
        from dip_c.commands.pos import pos
        return_value = pos(sys.argv[1:])
    elif command == "tad":
        from dip_c.commands.tad import tad
        return_value = tad(sys.argv[1:])
    elif command == "exp":
        from dip_c.commands.exp import exp
        return_value = exp(sys.argv[1:])
    elif command == "force":
        from dip_c.commands.force import force
        return_value = force(sys.argv[1:])
    elif command == "pairs2con":
        from dip_c.commands.pairs2con import pairs2con
        return_value = pairs2con(sys.argv[1:])
    elif command == "data-path":
        print(os.path.join(os.path.dirname(__file__), "data"))
        return 0
    else:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1

    if return_value == 0:
        sys.stderr.write("[M::" + __name__ + "] command: " + " ".join(sys.argv) + "\n")
        sys.stderr.write("[M::" + __name__ + "] finished in " + str(round(time.time() - start_time, 1)) + " sec\n")
    return return_value


if __name__ == "__main__":
    sys.exit(main())
