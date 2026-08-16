"""Regenerate the measured figures on docs/reliability/monte_carlo.md.

Each figure is a real campaign on a committed input file, drawn by the same
plot functions Studio uses, so the page's images are reproducible with:

    python tools/make_reliability_figures.py

Covered here: the sample-15 convergence trace and the VP34 histogram.
(lhs_vs_random.png and rs_surface_fit.png predate this script and do not yet
have producers.)
"""

import contextlib
import io
import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(_HERE)
OUT = os.path.join(REPO, "docs", "reliability", "images")

KEY = os.path.join(REPO, "docs", "lem", "files", "xslope_prob_submerged_KEY.xlsx")
VP34 = os.path.join(REPO, "docs", "verification", "files", "rocscience", "vp034.xlsx")


def _save(name, fig):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(path, REPO))


def convergence_trace():
    """Sample-15 Monte Carlo, drawn in full with the ±5% stop target it satisfies."""
    from xslope.fileio import load_slope_data
    from xslope.reliability import reliability_mc
    from xslope.plot import plot_reliability_convergence

    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = reliability_mc(load_slope_data(KEY), method="spencer",
                                 n_samples=10000, search=True)
    assert ok, res
    fig = plt.figure(figsize=(9, 5.5))
    plot_reliability_convergence(res, fig=fig, converge_rel=0.05)
    _save("mc_convergence_trace.png", fig)


def vp034_histogram():
    """VP34's Phase I fill at COV(phi) = 124% — the case Monte Carlo was added for."""
    from xslope.fileio import load_slope_data
    from xslope.reliability import reliability_mc
    from xslope.plot import plot_reliability_histogram

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with contextlib.redirect_stdout(io.StringIO()):
            ok, res = reliability_mc(load_slope_data(VP34), method="spencer",
                                     circular=False, search=False,
                                     n_samples=10000, num_slices=40)
    assert ok, res
    fig = plt.figure(figsize=(9, 5.5))
    plot_reliability_histogram(res, fig=fig)
    _save("reliability_mc_vp034.png", fig)


if __name__ == "__main__":
    convergence_trace()
    vp034_histogram()
