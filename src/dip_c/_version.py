"""Version information for dip-c."""

try:
    from importlib.metadata import version, PackageNotFoundError
    try:
        __version__ = version("run-dipc")
    except PackageNotFoundError:
        __version__ = "unknown"
except ImportError:
    __version__ = "unknown"
