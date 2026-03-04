"""Helpers for accessing dip-c package data files."""

import gzip
import os

_DATA_DIR = os.path.dirname(__file__)
COLOR_DIR = os.path.join(_DATA_DIR, "color")


def open_color_file(filename):
    """Open a color data file, transparently handling .gz compression.

    Tries ``<filename>.gz`` first (shipped in the wheel), then falls back
    to the plain ``<filename>`` (for development / legacy clone-and-run).

    Returns a text-mode file object.
    """
    gz_path = os.path.join(COLOR_DIR, filename + ".gz")
    if os.path.isfile(gz_path):
        return gzip.open(gz_path, "rt")

    plain_path = os.path.join(COLOR_DIR, filename)
    if os.path.isfile(plain_path):
        return open(plain_path, "r")

    raise FileNotFoundError(
        f"Color data file not found: {filename} (looked in {COLOR_DIR})"
    )
