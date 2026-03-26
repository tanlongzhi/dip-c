"""Shared utilities for Hi-C contact map plotting.

Provides matrix visualization, chromosome ordering, and genome-aware
layout helpers used by the hicplot* commands.

Requires optional dependencies: pip install run-dipc[hicplot]
"""

import sys
import re
import numpy as np

# Lazy imports for optional deps -- checked at call sites
_HICSTRAW = None
_PLT = None
_LSCM = None


def _ensure_hicplot_deps():
    """Import hic-straw, matplotlib; raise informative error if missing."""
    global _HICSTRAW, _PLT, _LSCM
    if _HICSTRAW is None:
        try:
            import hicstraw
            import matplotlib.pyplot as plt
            from matplotlib.colors import LinearSegmentedColormap
        except ImportError as exc:
            sys.stderr.write(
                "[E::hicplot] Missing dependencies. Install with:\n"
                "    pip install run-dipc[hicplot]\n"
            )
            raise SystemExit(1) from exc
        _HICSTRAW = hicstraw
        _PLT = plt
        _LSCM = LinearSegmentedColormap


VALID_NORMALIZATIONS = ["NONE", "KR", "VC", "VC_SQRT", "SCALE"]

# Colormaps
REDMAP_SPEC = ("bright_red", [(1, 1, 1), (1, 0, 0)])
BWRMAP_SPEC = ("blue_white_red", [(0, 0, 1), (1, 1, 1), (1, 0, 0)])

# Genome-specific chromosome counts per row for all-in-one layouts
_CHROMS_PER_ROW = {"mm10": 21, "hg19": 24, "hg38": 24}

# Genome-specific ordered chromosome lists
_ORDERED_CHROMS = {
    "mm10": [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chrX", "chrY",
    ],
    "hg19": [str(i) for i in range(1, 23)] + ["X", "Y"],
    "hg38": [str(i) for i in range(1, 23)] + ["X", "Y"],
}

# Chromosomes to skip in all-in-one plots
_SKIP_CHROMS = {"All", "ALL", "chrM", "M", "MT"}


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def validate_normalization(normalization):
    """Exit with error if *normalization* is not recognised."""
    if normalization not in VALID_NORMALIZATIONS:
        sys.stderr.write(
            "[E::hicplot] Invalid normalization: %s. "
            "Must be one of %s\n" % (normalization, VALID_NORMALIZATIONS)
        )
        raise SystemExit(1)




def ordering_check(chrom_1, chrom_2, ref_genome, ordered_list=None):
    """Return True if *chrom_1* should come before (or equal) *chrom_2*."""
    if ordered_list is None:
        ordered_list = _ORDERED_CHROMS.get(ref_genome, [])
    ordered_list = [c for c in ordered_list if c not in _SKIP_CHROMS]
    try:
        return ordered_list.index(chrom_1) <= ordered_list.index(chrom_2)
    except ValueError:
        return True  # unknown chroms: accept as-is


def make_colormap(spec):
    """Create a LinearSegmentedColormap from a ``(name, color_list)`` tuple."""
    _ensure_hicplot_deps()
    return _LSCM.from_list(spec[0], spec[1])


# ---------------------------------------------------------------------------
# Matrix plotting
# ---------------------------------------------------------------------------

def plot_matrix(input_matrix, output_png_path, cmap, vmin, vmax,
                multipl_factor=10, hlines=None, vlines=None):
    """Render *input_matrix* as a publication-quality PNG.

    Automatically switches to 1 px/bin for matrices ≥ 1000 bins on a side
    to avoid memory exhaustion.
    """
    _ensure_hicplot_deps()

    dpi = 300
    if input_matrix.shape[0] >= 1000:
        multipl_factor = 1

    figsize = (
        input_matrix.shape[1] * multipl_factor / dpi,
        input_matrix.shape[0] * multipl_factor / dpi,
    )
    sys.stderr.write("[M::hicplot] Figure size: %s\n" % (figsize,))

    fig = _PLT.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(input_matrix, cmap=cmap, vmin=vmin, vmax=vmax,
              aspect="equal", interpolation="nearest")

    if hlines:
        for y in hlines:
            ax.axhline(y=y, color="black", linewidth=0.5)
    if vlines:
        for x in vlines:
            ax.axvline(x=x, color="black", linewidth=0.5)

    _PLT.axis("off")
    _PLT.savefig(output_png_path, format="png", dpi=dpi,
                 bbox_inches="tight", pad_inches=0, transparent=True)
    _PLT.close()


# ---------------------------------------------------------------------------
# All-in-one layout helpers
# ---------------------------------------------------------------------------

def get_positions_and_matrix_size(matrices, matrices_per_row):
    """Compute tile positions for a grid of variable-size sub-matrices.

    Returns ``(positions, (total_rows, total_cols))`` where *positions* is a
    list of ``(row_offset, col_offset)`` tuples.
    """
    row_sizes, col_sizes = [], []

    for i in range(0, len(matrices), matrices_per_row):
        row_sizes.append(
            max(m.shape[0] for m in matrices[i:i + matrices_per_row])
        )
    for i in range(matrices_per_row):
        col_sizes.append(
            max(m.shape[1] for m in matrices[i::matrices_per_row])
        )

    positions = []
    current_row = 0
    for i, rs in enumerate(row_sizes):
        current_col = 0
        for j in range(matrices_per_row):
            idx = i * matrices_per_row + j
            if idx < len(matrices):
                positions.append((current_row, current_col))
                current_col += matrices[idx].shape[1]
        current_row += rs

    return positions, (sum(row_sizes), sum(col_sizes))


def chroms_per_row_for_genome(genome_id, num_chroms=None):
    """Return the number of chromosome tiles per row for *genome_id*.

    For unknown genomes, uses *num_chroms* (the actual chromosome count
    from the .hic file) so the all-in-one layout is always square.
    """
    if genome_id in _CHROMS_PER_ROW:
        return _CHROMS_PER_ROW[genome_id]
    fallback = num_chroms if num_chroms is not None else 24
    sys.stderr.write(
        "[W::hicplot] Unknown genome '%s'; using %d chroms/row\n"
        % (genome_id, fallback)
    )
    return fallback


def filter_chroms(chrom_list):
    """Remove unwanted chromosomes (All, chrM, MT, …) from *chrom_list*."""
    return [c for c in chrom_list if c[0] not in _SKIP_CHROMS]


def strip_png_ext(path):
    """Remove a trailing ``.png`` extension, if present."""
    return re.sub(r"\.png$", "", path)
