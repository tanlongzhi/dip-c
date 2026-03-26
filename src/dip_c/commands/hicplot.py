"""Plot Hi-C contact maps from .hic files.

Unified flag-based interface — the mode is inferred from which flags are
present:

  -2 supplied?    → difference map  (file 1 minus file 2)
  -c1 supplied?   → single-region   (one chromosomal window)
  --allinone?     → whole genome composited into one PNG
  (none of above) → per-chromosome  (one PNG per chromosome)

Requires:  pip install run-dipc[hicplot]
"""

import argparse
import sys
import numpy as np

import dip_c.hicplot_utils as _hu

from dip_c.hicplot_utils import (
    _ensure_hicplot_deps,
    VALID_NORMALIZATIONS,
    validate_normalization,
    ordering_check,
    make_colormap,
    plot_matrix,
    get_positions_and_matrix_size,
    chroms_per_row_for_genome,
    filter_chroms,
    strip_png_ext,
    REDMAP_SPEC,
    BWRMAP_SPEC,
)


# ══════════════════════════════════════════════════════════════════════════
# Input-parsing helpers
# ══════════════════════════════════════════════════════════════════════════

def _parse_file_spec(spec):
    """Parse ``path.hic:NORM`` → ``(path, norm_float)``.

    *NORM* is the contact-count factor used when ``-n NONE``
    (raw counts are divided by it).  Defaults to **1** when the colon
    portion is omitted, which effectively means "no scaling".
    """
    idx = spec.rfind(":")
    if idx > 0:
        path, norm_str = spec[:idx], spec[idx + 1:]
        try:
            return path, int(norm_str)
        except ValueError:
            try:
                return path, float(norm_str)
            except ValueError:
                pass  # not a number – treat whole string as path
    return spec, 1


def _parse_region(spec):
    """Parse ``chr1:1-10000000`` → ``('chr1', 1, 10000000)``.

    Raises :class:`argparse.ArgumentTypeError` on bad input so that
    argparse can produce a clean error message.
    """
    colon = spec.rfind(":")
    if colon < 0:
        raise argparse.ArgumentTypeError(
            "Region must be CHR:START-END, e.g. chr1:1-10000000, "
            "got '%s'" % spec
        )
    chrom = spec[:colon]
    coords = spec[colon + 1:]
    parts = coords.split("-")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            "Region coordinates must be START-END, got '%s'" % coords
        )
    try:
        return chrom, int(parts[0]), int(parts[1])
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Region coordinates must be integers, got '%s'" % coords
        )


# ══════════════════════════════════════════════════════════════════════════
# Absolute contact-map helpers
# ══════════════════════════════════════════════════════════════════════════

def _plot_map(hic_path, norm, output, chroms, region, resolution, maxcolor,
              normalization):
    """Single chromosomal region from one .hic file."""
    _ensure_hicplot_deps()
    validate_normalization(normalization)
    cmap = make_colormap(REDMAP_SPEC)

    hic = _hu._HICSTRAW.HiCFile(hic_path)
    sys.stderr.write("[M::hicplot] .hic file loaded\n")
    sys.stderr.write("[M::hicplot] Normalization: %s\n" % normalization)
    sys.stderr.write("[M::hicplot] Use Ctrl+C to force quit if normalization not present in the file\n")

    mo = hic.getMatrixZoomData(
        chroms[0], chroms[1], "observed", normalization, "BP", resolution,
    )
    sys.stderr.write("[M::hicplot] Matrix zoom data retrieved\n")

    if ordering_check(chroms[0], chroms[1], hic.getGenomeID()):
        mat = mo.getRecordsAsMatrix(region[0], region[1], region[2], region[3])
    else:
        mat = mo.getRecordsAsMatrix(region[2], region[3], region[0], region[1])
        mat = mat.transpose()

    if normalization == "NONE":
        mat = mat / norm

    sys.stderr.write("[M::hicplot] Matrix shape: %s\n" % (mat.shape,))
    plot_matrix(mat, output, cmap, 0, maxcolor)


def _plot_all(hic_path, norm, output, resolution, maxcolor, normalization):
    """Each chromosome individually — one PNG per chromosome."""
    _ensure_hicplot_deps()
    validate_normalization(normalization)
    cmap = make_colormap(REDMAP_SPEC)

    hic = _hu._HICSTRAW.HiCFile(hic_path)
    sys.stderr.write("[M::hicplot] .hic file loaded (genome: %s)\n"
                     % hic.getGenomeID())
    sys.stderr.write("[M::hicplot] Normalization: %s\n" % normalization)
    sys.stderr.write("[M::hicplot] Use Ctrl+C to force quit if normalization not present in the file\n")

    chroms = [(c.name, c.length) for c in hic.getChromosomes()]
    base = strip_png_ext(output)

    for name, length in chroms:
        mo = hic.getMatrixZoomData(
            name, name, "observed", normalization, "BP", resolution,
        )
        mat = mo.getRecordsAsMatrix(1, int(length), 1, int(length))
        if normalization == "NONE":
            mat = mat / norm

        plot_matrix(mat, "%s_%s.png" % (base, name), cmap, 0, maxcolor)
        sys.stderr.write("[M::hicplot] Chromosome %s complete\n" % name)


def _plot_allinone(hic_path, norm, output, resolution, maxcolor,
                   normalization, gridlines=True):
    """Whole genome composited into one image."""
    _ensure_hicplot_deps()
    validate_normalization(normalization)
    cmap = make_colormap(REDMAP_SPEC)

    hic = _hu._HICSTRAW.HiCFile(hic_path)
    genome_id = hic.getGenomeID()
    sys.stderr.write("[M::hicplot] .hic file loaded (genome: %s)\n"
                     % genome_id)
    sys.stderr.write("[M::hicplot] Normalization: %s\n" % normalization)
    sys.stderr.write("[M::hicplot] Use Ctrl+C to force quit if normalization not present in the file\n")

    chroms = filter_chroms([(c.name, c.length) for c in hic.getChromosomes()])
    chrom_names = [c[0] for c in chroms]
    matrices = []

    for row in chroms:
        for col in chroms:
            mo = hic.getMatrixZoomData(
                row[0], col[0], "observed", normalization, "BP", resolution,
            )
            if ordering_check(row[0], col[0], genome_id, chrom_names):
                m = mo.getRecordsAsMatrix(1, int(row[1]), 1, int(col[1]))
            else:
                m = mo.getRecordsAsMatrix(1, int(col[1]), 1, int(row[1]))
                m = m.transpose()
            matrices.append(m)

    per_row = chroms_per_row_for_genome(genome_id, len(chroms))
    positions, final_size = get_positions_and_matrix_size(matrices, per_row)
    sys.stderr.write("[M::hicplot] Sub-matrices: %d, final size: %s\n"
                     % (len(positions), final_size))

    final = np.zeros(final_size)
    for mat, (r, c) in zip(matrices, positions):
        rows, cols = mat.shape
        final[r:r + rows, c:c + cols] = mat

    if normalization == "NONE":
        final = final / norm

    hlines = vlines = None
    if gridlines:
        hlines = sorted(set(r for _, (r, _) in zip(matrices, positions)) - {0})
        vlines = sorted(set(c for _, (_, c) in zip(matrices, positions)) - {0})

    plot_matrix(final, output, cmap, 0, maxcolor, hlines=hlines, vlines=vlines)


# ══════════════════════════════════════════════════════════════════════════
# Difference-map helpers
# ══════════════════════════════════════════════════════════════════════════

def _diff(m1, norm1, m2, norm2, normalization):
    """(m1/norm1 - m2/norm2) when NONE, else (m1 - m2)."""
    if normalization == "NONE":
        return m1 / norm1 - m2 / norm2
    return m1 - m2


def _plot_diff_map(path1, norm1, path2, norm2, output, chroms, region,
                   resolution, maxcolor, normalization):
    """Difference map for a single chromosomal region."""
    _ensure_hicplot_deps()
    validate_normalization(normalization)
    cmap = make_colormap(BWRMAP_SPEC)

    hic1 = _hu._HICSTRAW.HiCFile(path1)
    hic2 = _hu._HICSTRAW.HiCFile(path2)
    sys.stderr.write("[M::hicplot] Both .hic files loaded\n")
    sys.stderr.write("[M::hicplot] Normalization: %s\n" % normalization)
    sys.stderr.write("[M::hicplot] Use Ctrl+C to force quit if normalization not present in the file\n")

    mo1 = hic1.getMatrixZoomData(
        chroms[0], chroms[1], "observed", normalization, "BP", resolution,
    )
    mo2 = hic2.getMatrixZoomData(
        chroms[0], chroms[1], "observed", normalization, "BP", resolution,
    )

    if ordering_check(chroms[0], chroms[1], hic1.getGenomeID()):
        m1 = mo1.getRecordsAsMatrix(region[0], region[1], region[2], region[3])
        m2 = mo2.getRecordsAsMatrix(region[0], region[1], region[2], region[3])
    else:
        m1 = mo1.getRecordsAsMatrix(region[2], region[3], region[0], region[1])
        m2 = mo2.getRecordsAsMatrix(region[2], region[3], region[0], region[1])
        m1, m2 = m1.transpose(), m2.transpose()

    result = _diff(m1, norm1, m2, norm2, normalization)
    sys.stderr.write("[M::hicplot] Diff matrix shape: %s\n" % (result.shape,))
    plot_matrix(result, output, cmap, -maxcolor, maxcolor)


def _plot_diff_all(path1, norm1, path2, norm2, output, resolution, maxcolor,
                   normalization):
    """Per-chromosome difference maps."""
    _ensure_hicplot_deps()
    validate_normalization(normalization)
    cmap = make_colormap(BWRMAP_SPEC)

    hic1 = _hu._HICSTRAW.HiCFile(path1)
    hic2 = _hu._HICSTRAW.HiCFile(path2)
    sys.stderr.write("[M::hicplot] Both .hic files loaded\n")
    sys.stderr.write("[M::hicplot] Normalization: %s\n" % normalization)
    sys.stderr.write("[M::hicplot] Use Ctrl+C to force quit if normalization not present in the file\n")

    chroms1 = [(c.name, c.length) for c in hic1.getChromosomes()]
    chroms2 = [(c.name, c.length) for c in hic2.getChromosomes()]
    base = strip_png_ext(output)

    for i in range(len(chroms1)):
        mo1 = hic1.getMatrixZoomData(
            chroms1[i][0], chroms1[i][0], "observed", normalization,
            "BP", resolution,
        )
        mo2 = hic2.getMatrixZoomData(
            chroms2[i][0], chroms2[i][0], "observed", normalization,
            "BP", resolution,
        )
        m1 = mo1.getRecordsAsMatrix(1, int(chroms1[i][1]),
                                     1, int(chroms1[i][1]))
        m2 = mo2.getRecordsAsMatrix(1, int(chroms2[i][1]),
                                     1, int(chroms2[i][1]))

        result = _diff(m1, norm1, m2, norm2, normalization)
        plot_matrix(result, "%s_%s.png" % (base, chroms1[i][0]),
                    cmap, -maxcolor, maxcolor)
        sys.stderr.write("[M::hicplot] Chromosome %s complete\n"
                         % chroms1[i][0])


def _plot_diff_allinone(path1, norm1, path2, norm2, output, resolution,
                        maxcolor, normalization, gridlines=True):
    """Whole-genome difference map in one image."""
    _ensure_hicplot_deps()
    validate_normalization(normalization)
    cmap = make_colormap(BWRMAP_SPEC)

    hic1 = _hu._HICSTRAW.HiCFile(path1)
    hic2 = _hu._HICSTRAW.HiCFile(path2)
    genome_id = hic1.getGenomeID()
    sys.stderr.write("[M::hicplot] Both .hic files loaded (genome: %s)\n"
                     % genome_id)
    sys.stderr.write("[M::hicplot] Normalization: %s\n" % normalization)
    sys.stderr.write("[M::hicplot] Use Ctrl+C to force quit if normalization not present in the file\n")

    chroms1 = filter_chroms([(c.name, c.length) for c in hic1.getChromosomes()])
    chrom_names = [c[0] for c in chroms1]
    matrices = []

    for i in range(len(chroms1)):
        for j in range(len(chroms1)):
            cr, cc = chroms1[i], chroms1[j]
            mo1 = hic1.getMatrixZoomData(
                cr[0], cc[0], "observed", normalization, "BP", resolution,
            )
            mo2 = hic2.getMatrixZoomData(
                cr[0], cc[0], "observed", normalization, "BP", resolution,
            )
            if ordering_check(cr[0], cc[0], genome_id, chrom_names):
                m1 = mo1.getRecordsAsMatrix(1, int(cr[1]), 1, int(cc[1]))
                m2 = mo2.getRecordsAsMatrix(1, int(cr[1]), 1, int(cc[1]))
            else:
                m1 = mo1.getRecordsAsMatrix(1, int(cc[1]), 1, int(cr[1]))
                m2 = mo2.getRecordsAsMatrix(1, int(cc[1]), 1, int(cr[1]))
                m1, m2 = m1.transpose(), m2.transpose()

            matrices.append(_diff(m1, norm1, m2, norm2, normalization))

    per_row = chroms_per_row_for_genome(genome_id, len(chroms1))
    positions, final_size = get_positions_and_matrix_size(matrices, per_row)
    sys.stderr.write("[M::hicplot] Sub-matrices: %d, final size: %s\n"
                     % (len(positions), final_size))

    final = np.zeros(final_size)
    for mat, (r, c) in zip(matrices, positions):
        rows, cols = mat.shape
        final[r:r + rows, c:c + cols] = mat

    hlines = vlines = None
    if gridlines:
        hlines = sorted(set(r for _, (r, _) in zip(matrices, positions)) - {0})
        vlines = sorted(set(c for _, (_, c) in zip(matrices, positions)) - {0})

    plot_matrix(final, output, cmap, -maxcolor, maxcolor,
                hlines=hlines, vlines=vlines)


# ══════════════════════════════════════════════════════════════════════════
# Argument parser
# ══════════════════════════════════════════════════════════════════════════

def _build_parser():
    p = argparse.ArgumentParser(
        prog="dip-c hicplot",
        description="Plot Hi-C contact maps from .hic files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
        epilog="""\
Usage:
  dip-c hicplot -1 <file.hic>:<norm> -o <out.png> -r <bp> -s <val>
                [-2 <file2.hic>:<norm>] [-c1 <chr:start-end>] [-c2 <chr:start-end>]
                [-n <TYPE>] [--allinone] [-G]

Mode is inferred from flags:
  -2 present        → difference map (file 1 minus file 2)
  -c1 present       → single region
  -c1 absent        → all chromosomes (one PNG each)
  --allinone        → all chromosomes composited into one PNG

File spec format:
  FILE:NORM         e.g. abc.hic:123456789
  The NORM factor is only used with -n NONE (default).
  Omit the colon+number when using built-in normalizations (KR, VC, …).

Region format:
  CHR:START-END     e.g. chr1:1-10000000
  If only -c1 is given, -c2 defaults to the same region (symmetric).

Examples:
  # Absolute map, single symmetric region (SCALE normalization)
  dip-c hicplot -1 abc.hic -c1 chr1:1-10000000 \\
                -o out.png -r 1000000 -s 1.5e-6 -n SCALE

  # Absolute map, asymmetric region (chr1 vs chr2)
  dip-c hicplot -1 abc.hic -c1 chr1:1-10000000 \\
                -c2 chr2:1-10000000 -o out.png -r 1000000 -s 1.5e-6 -n SCALE

  # Absolute map, all chromosomes (one PNG each)
  dip-c hicplot -1 abc.hic -o out.png -r 1000000 -s 1.5e-6 -n SCALE

  # Absolute map, whole genome in one image (with gridlines)
  dip-c hicplot -1 abc.hic --allinone \\
                -o out.png -r 1000000 -s 1.5e-6 -n SCALE

  # Same but without gridlines
  dip-c hicplot -1 abc.hic --allinone -G \\
                -o out.png -r 1000000 -s 1.5e-6 -n SCALE

  # Difference map, single region, KR normalization
  dip-c hicplot -1 abc.hic -2 def.hic -c1 chr1:1-10000000 \\
                -o diff.png -r 1000000 -s 1.5e-6 -n KR

  # Difference map, all chromosomes (NONE normalization, :NORM required)
  dip-c hicplot -1 abc.hic:123456789 -2 def.hic:987654321 \\
                -o diff.png -r 1000000 -s 1.5e-6

  # Difference map, whole genome in one image
  dip-c hicplot -1 abc.hic:123456789 -2 def.hic:987654321 --allinone \\
                -o diff.png -r 1000000 -s 1.5e-6
""",
    )

    # -- Required --------------------------------------------------------
    p.add_argument(
        "-1", dest="file1", required=True,
        metavar="HIC:NORM",
        help="First .hic file with optional normalization factor "
             "(e.g. abc.hic:123456789).  The NORM part is the total "
             "contact count used to scale raw matrices when -n NONE.  "
             "Omit the :NORM when using built-in normalizations.",
    )
    p.add_argument(
        "-o", "--output", required=True,
        metavar="PNG",
        help="Output PNG path. For per-chromosome mode, used as the base name "
             "(e.g. out.png → out_chr1.png, out_chr2.png, …).",
    )
    p.add_argument(
        "-r", "--resolution", type=int, required=True,
        metavar="BP",
        help="Resolution in base pairs (e.g. 1000000 for 1 Mb).",
    )
    p.add_argument(
        "-s", "--scale", type=float, required=True,
        metavar="VAL",
        help="Maximum value on the colour scale. For difference maps the "
             "range is [-scale, +scale].",
    )

    # -- Optional: second file (triggers difference mode) ----------------
    p.add_argument(
        "-2", dest="file2", default=None,
        metavar="HIC:NORM",
        help="Second .hic file (with optional :NORM).  When supplied, "
             "produces a difference map (file 1 − file 2).",
    )

    # -- Optional: regions (triggers single-region mode) -----------------
    p.add_argument(
        "-c1", dest="region1", default=None,
        metavar="CHR:START-END",
        help="First genomic region, e.g. chr1:1-10000000.  "
             "Defines the row axis of the contact matrix.  "
             "If -c2 is omitted, the same region is used for both axes "
             "(symmetric view, which is the most common case).",
    )
    p.add_argument(
        "-c2", dest="region2", default=None,
        metavar="CHR:START-END",
        help="Second genomic region for the column axis.  "
             "Only needed for asymmetric views (e.g. chr1 vs chr2).  "
             "Defaults to the same region as -c1.",
    )

    # -- Optional: all-in-one --------------------------------------------
    p.add_argument(
        "--allinone", action="store_true", default=False,
        help="Composite all chromosomes into a single image instead of "
             "producing one PNG per chromosome.  Chromosome boundary "
             "gridlines are drawn by default; use -G to turn them off.  "
             "Ignored when -c1 is given.",
    )

    # -- Optional: gridlines off -----------------------------------------
    p.add_argument(
        "-G", "--no-gridlines", dest="no_gridlines",
        action="store_true", default=False,
        help="Turn off chromosome boundary gridlines in --allinone mode.",
    )

    # -- Optional: normalization -----------------------------------------
    p.add_argument(
        "-n", "--normalization", default="NONE",
        choices=VALID_NORMALIZATIONS,
        metavar="TYPE",
        help="Normalization method (default: NONE). "
             "Choices: %(choices)s. "
             "When not NONE, the :NORM values in -1/-2 are ignored and "
             "hic-straw's built-in balancing is used instead.",
    )

    return p


# ══════════════════════════════════════════════════════════════════════════
# CLI entry point  –  called from dip_c.cli as  ``dip-c hicplot …``
# ══════════════════════════════════════════════════════════════════════════

def hicplot(argv):
    """Parse flags and dispatch to the correct plotting function."""
    parser = _build_parser()

    # argv arrives as ["hicplot", "-1", …]; skip the command word itself
    args = parser.parse_args(argv[1:])

    # Unpack file specs
    path1, norm1 = _parse_file_spec(args.file1)
    output = args.output
    resolution = args.resolution
    maxcolor = args.scale
    normalization = args.normalization

    # Validate: -n NONE requires :NORM on file specs
    if normalization == "NONE":
        if norm1 == 1:
            parser.error(
                "-n NONE (default) requires a contact count on -1, "
                "e.g. -1 abc.hic:123456789. "
                "Use -n SCALE (or KR, VC, …) to skip the :NORM.")

    is_diff = args.file2 is not None
    is_region = args.region1 is not None
    is_allinone = args.allinone
    gridlines = is_allinone and not args.no_gridlines
    if is_allinone and not gridlines:
        sys.stderr.write("[M::hicplot] Gridlines: off (-G)\n")
    elif is_allinone:
        sys.stderr.write("[M::hicplot] Gridlines: on (use -G to turn off)\n")

    # -c2 without -c1 is an error
    if args.region2 is not None and args.region1 is None:
        parser.error("-c2 requires -c1")

    # Parse regions if provided
    chroms = region = None
    if is_region:
        r1_chrom, r1_start, r1_end = _parse_region(args.region1)
        if args.region2 is not None:
            r2_chrom, r2_start, r2_end = _parse_region(args.region2)
        else:
            r2_chrom, r2_start, r2_end = r1_chrom, r1_start, r1_end
        chroms = [r1_chrom, r2_chrom]
        region = [r1_start, r1_end, r2_start, r2_end]

    # --allinone is meaningless with -c1 (single region already)
    if is_region and is_allinone:
        sys.stderr.write(
            "[W::hicplot] --allinone ignored when -c1 is specified "
            "(single region already)\n"
        )
        is_allinone = False

    # ── Dispatch ───────────────────────────────────────────────────────
    if is_diff:
        path2, norm2 = _parse_file_spec(args.file2)
        if normalization == "NONE" and norm2 == 1:
            parser.error(
                "-n NONE (default) requires a contact count on -2, "
                "e.g. -2 def.hic:987654321. "
                "Use -n SCALE (or KR, VC, …) to skip the :NORM.")

        if is_region:
            _plot_diff_map(path1, norm1, path2, norm2, output,
                           chroms, region, resolution, maxcolor,
                           normalization)
        elif is_allinone:
            _plot_diff_allinone(path1, norm1, path2, norm2, output,
                                resolution, maxcolor, normalization,
                                gridlines=gridlines)
        else:
            _plot_diff_all(path1, norm1, path2, norm2, output,
                           resolution, maxcolor, normalization)
    else:
        if is_region:
            _plot_map(path1, norm1, output, chroms, region,
                      resolution, maxcolor, normalization)
        elif is_allinone:
            _plot_allinone(path1, norm1, output, resolution, maxcolor,
                           normalization, gridlines=gridlines)
        else:
            _plot_all(path1, norm1, output, resolution, maxcolor,
                      normalization)

    return 0
