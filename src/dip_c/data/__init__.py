"""Helpers for accessing dip-c package data files."""

import gzip
import os

_DATA_DIR = os.path.dirname(__file__)
COLOR_DIR = os.path.join(_DATA_DIR, "color")


def list_bundled_data_files():
    """Return a sorted list of bundled data filenames (without .gz suffix)."""
    files = []
    for name in os.listdir(COLOR_DIR):
        if name.endswith(".gz"):
            files.append(name[:-3])
        else:
            files.append(name)
    return sorted(set(files))


def resolve_data_file(path):
    """Open a data file with automatic fallback to bundled data.

    Resolution order:
      1. If *path* exists on disk, open it directly (backward compat).
      2. Strip to basename, try ``<basename>.gz`` in bundled COLOR_DIR.
      3. Try ``<basename>`` plain in bundled COLOR_DIR.
      4. Raise ``FileNotFoundError`` with a helpful message.

    Returns a text-mode file object.
    """
    # 1. Local path exists — use it directly
    if os.path.isfile(path):
        return open(path, "r")

    # 2–3. Try bundled data by basename
    basename = os.path.basename(path)

    gz_path = os.path.join(COLOR_DIR, basename + ".gz")
    if os.path.isfile(gz_path):
        return gzip.open(gz_path, "rt")

    plain_path = os.path.join(COLOR_DIR, basename)
    if os.path.isfile(plain_path):
        return open(plain_path, "r")

    # 4. Not found anywhere
    available = list_bundled_data_files()
    raise FileNotFoundError(
        f"Data file not found: {path}\n"
        f"Available bundled files: {', '.join(available)}"
    )
