from ._version import __version__
# Project packaging. Imported here so ``xslope.pack`` / ``xslope.unpack`` are reachable
# from a bare ``import xslope``; the module is pure stdlib, so this costs the package
# import nothing.
from .package import PACKAGE_EXT, pack, project_files, unpack
