"""XSlope Studio — a desktop GUI for the xslope slope-stability engine.

The GUI is a thin PySide6 layer over the existing ``xslope`` package: it opens
Excel input files, renders the geometry/inputs with the engine's own Matplotlib
plotting, and (in later phases) edits inputs and runs analyses. Importing this
package is cheap; PySide6 is only imported when the app is actually launched.
"""

__all__ = ["main"]


def main(argv=None):
    """Launch the XSlope Studio desktop application."""
    from .app import main as _main
    return _main(argv)
