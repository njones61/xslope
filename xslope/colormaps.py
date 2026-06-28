"""Named color ramps for ramp-based plots, including the custom ``BGYR`` ramp.

Registers ``BGYR`` (blue → green → yellow → red) with Matplotlib's global
colormap registry so it resolves by name (``cmap="BGYR"``) anywhere a Matplotlib
colormap name is accepted — the same as the built-in ramps. ``RAMP_CHOICES`` is
the curated list backing the Studio's colormap selector: the screenshot ramps
plus the ramps the plots currently default to (Spectral / coolwarm), so each
plot's existing default is always a selectable option.
"""

from __future__ import annotations

import matplotlib
from matplotlib.colors import LinearSegmentedColormap

BGYR_NAME = "BGYR"

# Curated ramps offered in the colormap selector: (matplotlib name, display label).
# Order matches the Studio selector top-to-bottom. The set mixes the classic
# high-contrast "rainbow" ramps engineers expect for reading contour bands
# (BGYR, Turbo, Spectral), the perceptually-uniform / colorblind-safe sequential
# ramps (Viridis, Plasma, Inferno, Cividis), a single-hue option (Blues), and a
# diverging ramp for signed/centered data (Coolwarm). Spectral_r and coolwarm
# are kept because they are the plots' existing defaults.
RAMP_CHOICES = [
    (BGYR_NAME, "BGYR"),
    ("turbo", "Turbo"),
    ("Spectral_r", "Spectral"),
    ("viridis", "Viridis"),
    ("plasma", "Plasma"),
    ("inferno", "Inferno"),
    ("cividis", "Cividis"),
    ("Blues", "Blues"),
    ("coolwarm", "Coolwarm"),
]


def _bgyr():
    """Blue → green → yellow → red ramp (low to high)."""
    return LinearSegmentedColormap.from_list(
        BGYR_NAME, ["#0000ff", "#00b050", "#ffff00", "#ff0000"], N=256)


def register_colormaps():
    """Register the custom ramps with Matplotlib (idempotent) so they resolve by name."""
    if BGYR_NAME not in matplotlib.colormaps:
        matplotlib.colormaps.register(_bgyr())


register_colormaps()
