"""Shared utilities for Hi-C contact map plotting.

Provides matrix visualization, chromosome ordering, and genome-aware
layout helpers used by the hicplot* commands.

Requires optional dependencies: pip install run-dipc[hicplot]
"""

import os
import sys
import re
from contextlib import contextmanager

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


# Colormaps
REDMAP_SPEC = ("bright_red", [(1, 1, 1), (1, 0, 0)])
BWRMAP_SPEC = ("blue_white_red", [(0, 0, 1), (1, 1, 1), (1, 0, 0)])

# Chromosomes skipped by default in all-in-one layouts.
# "All"/"ALL" are Juicer pseudo-chromosomes representing the whole genome.
DEFAULT_SKIP_CHROMS = {"All", "ALL"}


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def get_chrom_names(hic_file):
    """Return the ordered list of chromosome names from a HiCFile object."""
    return [c.name for c in hic_file.getChromosomes()]



def ordering_check(chrom_1, chrom_2, chrom_order):
    """Return True if *chrom_1* comes before (or equals) *chrom_2*.

    *chrom_order* should be the list returned by :func:`get_chrom_names`.
    """
    try:
        return chrom_order.index(chrom_1) <= chrom_order.index(chrom_2)
    except ValueError:
        return True


def make_colormap(spec):
    """Create a LinearSegmentedColormap from a ``(name, color_list)`` tuple."""
    _ensure_hicplot_deps()
    return _LSCM.from_list(spec[0], spec[1])


def resolve_colormap(spec):
    """Resolve a user-supplied colormap specification.

    *spec* is either a **matplotlib named colormap** (e.g. ``viridis``,
    ``coolwarm``, ``RdBu_r``) or a **comma-separated list of colour stops**
    understood by :func:`matplotlib.colors.to_rgba` (named colours, hex
    codes, etc.).

    Examples::

        resolve_colormap("viridis")
        resolve_colormap("white,red")
        resolve_colormap("blue,white,red")
        resolve_colormap("#0000ff,#ffffff,#ff0000")
    """
    _ensure_hicplot_deps()
    import matplotlib
    from matplotlib.colors import to_rgba

    try:
        return matplotlib.colormaps[spec]
    except (KeyError, ValueError):
        pass

    stops = [s.strip() for s in spec.split(",")]
    if len(stops) < 2:
        sys.stderr.write(
            "[E::hicplot] Colormap spec must be a matplotlib name or "
            "at least two comma-separated colors, got '%s'\n" % spec
        )
        raise SystemExit(1)

    try:
        colors = [to_rgba(s) for s in stops]
    except ValueError as exc:
        sys.stderr.write(
            "[E::hicplot] Bad color in colormap spec '%s': %s\n" % (spec, exc)
        )
        raise SystemExit(1) from exc

    return _LSCM.from_list("custom", colors)


# ---------------------------------------------------------------------------
# Matrix plotting
# ---------------------------------------------------------------------------

def plot_matrix(input_matrix, output_png_path, cmap, vmin, vmax,
                multipl_factor=10, hlines=None, vlines=None):
    """Render *input_matrix* as a publication-quality PNG.

    Automatically switches to 1 px/bin for matrices >= 1000 bins on a side
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


def chroms_per_row_for_genome(num_chroms):
    """Return the number of chromosome tiles per row for all-in-one layout.

    Uses the actual chromosome count from the .hic file so the
    all-in-one layout is always a square grid.
    """
    return num_chroms


def filter_chroms(chrom_list):
    """Remove Juicer's ``All`` pseudo-chromosome from *chrom_list*."""
    return [c for c in chrom_list if c[0] not in DEFAULT_SKIP_CHROMS]


@contextmanager
def _suppress_native_stderr():
    """Temporarily redirect file-descriptor 2 to /dev/null.

    Suppresses warnings printed by C/C++ libraries (e.g. hic-straw)
    that write directly to fd 2 and cannot be caught from Python.
    Python-level sys.stderr is also silenced while inside the context,
    so keep usage narrow.
    """
    sys.stderr.flush()
    old_fd = os.dup(2)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull_fd, 2)
    os.close(devnull_fd)
    try:
        yield
    finally:
        sys.stderr.flush()
        os.dup2(old_fd, 2)
        os.close(old_fd)


def validate_chroms(hic, chroms, normalization, resolution):
    """Return only chromosomes that have actual contact data.

    Probes each chromosome with a full self-vs-self query at the given
    resolution.  Chromosomes whose matrix is empty or all-zero are
    dropped.  C-level hic-straw warnings are suppressed during the
    probe so they do not clutter stderr.

    Returns ``(valid_chroms, skipped_names)``.
    """
    _ensure_hicplot_deps()
    valid = []
    skipped = []
    with _suppress_native_stderr():
        for name, length in chroms:
            try:
                mo = hic.getMatrixZoomData(
                    name, name, "observed", normalization, "BP", resolution,
                )
                m = mo.getRecordsAsMatrix(1, int(length), 1, int(length))
                if m.size > 0 and m.sum() != 0:
                    valid.append((name, length))
                else:
                    skipped.append(name)
            except Exception:
                skipped.append(name)
    return valid, skipped


def strip_png_ext(path):
    """Remove a trailing ``.png`` extension, if present."""
    return re.sub(r"\.png$", "", path)
