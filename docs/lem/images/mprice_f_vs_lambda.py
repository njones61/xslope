"""Regenerate mprice_f_vs_lambda.png for docs/lem/mprice.md.

Traces the two Morgenstern-Price factor-of-safety curves on the Arai & Tagyo
benchmark (verification corpus LEM-2b / VP14): at each interslice scale factor
lambda, F_f is the factor of safety that zeroes the force residual and F_m the
one that zeroes the moment residual. The Morgenstern-Price solution is the
lambda where the two curves cross.

The surface is the critical circle from the benchmark's own circular search
(method mprice, 50 slices), so the crossing reproduces the published factor of
safety for that row.

    PYTHONPATH=. python3 docs/lem/images/mprice_f_vs_lambda.py
"""
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, newton

from xslope.fileio import load_slope_data
from xslope.search import circular_search
from xslope.solve import _mp_extract, _mp_f_vals, _mp_residuals, bishop, mprice

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
XLSX = os.path.join(HERE, "..", "files", "xslope_arai_tagyo.xlsx")
OUT = os.path.join(HERE, "mprice_f_vs_lambda.png")

NUM_SLICES = 50
F_TYPE = "half_sine"
LAMBDAS = np.linspace(-0.2, 0.9, 23)


def critical_slices():
    """Slice dataframe for the benchmark's critical circle (mprice, 50 slices)."""
    slope_data = load_slope_data(XLSX)
    fs_cache, _converged, _path, _cache = circular_search(
        slope_data, "mprice", num_slices=NUM_SLICES)
    if not fs_cache:
        raise RuntimeError("circular_search found no valid surface")
    return fs_cache[0]["slices"]


def curves(slice_df):
    """(F_f, F_m) at every lambda in LAMBDAS, plus the solved (lambda*, F*)."""
    right_facing = bool(slice_df.attrs.get(
        "right_facing", slice_df["y_cb"].values[0] > slice_df["y_cb"].values[-1]))
    A = _mp_extract(slice_df, right_facing)
    f_vals = _mp_f_vals(slice_df, F_TYPE)

    ok, res = bishop(slice_df.copy())
    seed = res["FS"] if ok else 1.5

    def force_res(lam, FS):
        return _mp_residuals(A, f_vals, lam, FS)[2]

    def moment_res(lam, FS):
        return _mp_residuals(A, f_vals, lam, FS)[3]

    F_f, F_m = [], []
    for lam in LAMBDAS:
        # The force residual is monotonic in FS — a secant from the Bishop seed
        # lands on its single root.
        ff = newton(lambda FS: force_res(lam, FS), seed, tol=1e-8, maxiter=60)
        # The moment residual is multivalued in FS and carries an asymptote, so
        # it is bracketed around F_f rather than seeded.
        lo, hi = 0.25 * ff, 2.0 * ff
        r_lo, r_hi = moment_res(lam, lo), moment_res(lam, hi)
        while r_lo * r_hi > 0 and lo > 0.05:
            lo *= 0.8
            r_lo = moment_res(lam, lo)
        fm = brentq(lambda FS: moment_res(lam, FS), lo, hi, xtol=1e-10)
        F_f.append(ff)
        F_m.append(fm)

    ok, res = mprice(slice_df.copy(), f_type=F_TYPE)
    if not ok:
        raise RuntimeError(res)
    return np.array(F_f), np.array(F_m), res["lambda"], res["FS"]


def main():
    F_f, F_m, lam_star, FS = curves(critical_slices())

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(LAMBDAS, F_f, "-o", ms=4, lw=1.6, color="tab:blue",
            label=r"$F_f$ (force equilibrium)")
    ax.plot(LAMBDAS, F_m, "-s", ms=4, lw=1.6, color="tab:green",
            label=r"$F_m$ (moment equilibrium)")
    ax.axvline(lam_star, color="0.5", ls="--", lw=1.0)
    ax.axhline(FS, color="0.5", ls="--", lw=1.0)
    ax.plot([lam_star], [FS], "*", ms=18, color="red", zorder=6,
            label=rf"M-P solution: $F$ = {FS:.3f}, $\lambda$ = {lam_star:.3f}")

    ax.set_xlabel(r"interslice scale factor  $\lambda$")
    ax.set_ylabel(r"factor of safety  $F$")
    ax.set_title(r"Morgenstern-Price: force and moment factors of safety vs $\lambda$")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT}  (F = {FS:.5f}, lambda = {lam_star:.5f})")


if __name__ == "__main__":
    main()
