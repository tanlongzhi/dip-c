#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Version information for dip-c.

This module provides version information in a format compatible with
both Python 2 and Python 3.
"""

import os

# Read version from VERSION file
_version_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'VERSION')
try:
    with open(_version_file, 'r') as f:
        __version__ = f.read().strip()
except IOError:
    __version__ = 'unknown'

# Version components
version_info = tuple(int(x) for x in __version__.split('.') if x.isdigit())

# Convenience function
def get_version():
    """Return the version string."""
    return __version__

if __name__ == '__main__':
    print(__version__)
