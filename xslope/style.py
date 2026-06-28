"""Feature styles for the plot functions (plans/plan_studio.md §8a).

A *style* is the persistent per-feature visual identity — color, line style, line
width, fill hatch, alpha. ``default_style_sheet()`` holds the baked-in defaults
(the "factory", expressed in code); a project may override a sparse subset, which
``resolve_style()`` deep-merges over the defaults. The ``plot_*`` functions take an
optional ``style=`` kwarg and fall back to these defaults, so **passing nothing
reproduces the original hardcoded look** — notebooks/CLI are unaffected.

Two kinds of entry:

- ``materials`` — per-``mat_id`` overrides ``{color, hatch, alpha}`` (sparse; keyed
  by the 0-based index, matching ``mat_id`` on polygons/profile lines). A material's
  default color is the tab10 palette (``plot.get_material_color``); ``material`` holds
  the defaults shared by all materials (hatch/alpha).
- ``features`` — per feature type ``{color, linestyle, linewidth, alpha, ...}``.
  Profile lines and polygons inherit their *color* from the material and carry only
  their own stroke/alpha here.
"""

from __future__ import annotations

import copy


def default_style_sheet():
    """The complete default style — the values the ``plot_*`` functions hardcoded.

    Material *color* defaults are not listed (they come from the tab10 palette via
    ``plot.get_material_color``); per-material overrides live under ``materials``.
    """
    return {
        "version": 1,
        # Per-material sparse overrides, keyed by str(mat_id): {color, hatch, alpha}.
        "materials": {},
        # Defaults shared by every material's fill (overridden per-id above).
        "material": {"hatch": None, "alpha": 0.6},
        # Per-feature-type styling. Profile/polygon take color from the material;
        # every other feature carries its own color/linestyle/width/alpha. Defaults
        # match the literals the plot_* functions used, so no-style output is
        # unchanged.
        "features": {
            "profile_line": {"linewidth": 1.0, "linestyle": "-"},
            "polygon": {"linewidth": 1.0},
            "piezo_line": {"color": "b", "linewidth": 2.0},
            "piezo_line2": {"color": "skyblue", "linewidth": 2.0},
            "max_depth": {"color": "black", "linewidth": 1.5},
            "tcrack": {"color": "red", "linestyle": ":", "linewidth": 1.5},
            "circles": {"color": "red", "linestyle": "--", "linewidth": 1.5},
            "noncirc": {"color": "red", "linestyle": "--", "linewidth": 1.5},
            "reinforcement": {"color": "darkgray", "linewidth": 3.0,
                              "linestyle": "-", "alpha": 0.8},
            "mesh": {"color": "gray", "linewidth": 0.5, "alpha": 0.25},
        },
    }


def _deep_merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def resolve_style(style):
    """``None`` → the defaults; a (sparse) override dict → deep-merged over them.
    Idempotent: resolving an already-resolved sheet returns an equivalent sheet."""
    if style is None:
        return default_style_sheet()
    return _deep_merge(default_style_sheet(), style)


def feature_style(style, name):
    """The resolved ``features[name]`` sub-dict (already merged over defaults)."""
    return (style.get("features") or {}).get(name, {})


def material_style(style, mat_id):
    """Resolved ``{color, hatch, alpha, ...}`` for a material: the shared ``material``
    defaults, overlaid with any per-``mat_id`` override, with the color falling back
    to the tab10 palette when not overridden."""
    out = dict(style.get("material") or {})
    out.update((style.get("materials") or {}).get(str(mat_id)) or {})
    if not out.get("color"):
        from .plot import get_material_color
        out["color"] = get_material_color(mat_id)
    return out


def style_delta(style):
    """The sparse diff of a (resolved) style vs the defaults — what a project sidecar
    stores. Inverse of ``resolve_style`` for round-tripping; empty dict if unchanged."""
    return _diff(default_style_sheet(), style or {})


def _diff(base, cur):
    out = {}
    for k, v in (cur or {}).items():
        bv = base.get(k)
        if isinstance(v, dict) and isinstance(bv, dict):
            d = _diff(bv, v)
            if d:
                out[k] = d
        elif v != bv:
            out[k] = copy.deepcopy(v)
    return out
