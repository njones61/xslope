"""Regenerate the figures of the Reliability pages (docs/reliability/images/).

Every figure is either a computed curve (no model behind it) or drawn from a
real XSLOPE run on a committed input file, so the pages' images are
reproducible with:

    python tools/make_reliability_figures.py                 # every figure
    python tools/make_reliability_figures.py same_f beta_pf  # named figures
    python tools/make_reliability_figures.py --list          # the names

The runs (each made once per invocation, and only for the figures asked for):

- LEM-11 (docs/lem/files/xslope_reliability.xlsx, Spencer, 40 slices): the
  Taylor series with a search at each solve; Monte Carlo, 10,000 Latin
  hypercube realizations at the default seed, searched once at the means; and
  the response surface, ten million realizations of the fitted formula.
- VP35, Cannon Dam (docs/verification/files/rocscience/vp035.xlsx, Bishop): the
  critical circle of a search at the means (50 slices, as the verification tag
  runs it) and the Taylor series on that circle and on the file's own circle,
  each held fixed (40 slices).
- VP36 (vp036.xlsx, Bishop): the Taylor series with a search at each solve.
- VP34 (vp034.xlsx, Spencer, its noncircular surface): Monte Carlo, 10,000.
- LEM-11 again, for the scatter between runs: twenty Monte Carlo runs of 1,000
  realizations each way (random and Latin hypercube, seeds 1 to 20), on the
  critical circle of a search at the means.
- LEM-11 again, for the segments of a Latin hypercube: Monte Carlo at ten Latin
  hypercube realizations, the default seed, searched once at the means.

``--cache DIR`` keeps the runs' results in DIR between invocations (for
restyling without re-solving); the figures are the same either way.

Figures (page: name):

- index.md: same_f, r_pf, sigma_normal, sigma_rules, beta_on_f, beta_pf,
  vp035_surfaces
- taylor.md: taylor_gamma, taylor_shapes, reliability_taylor_vp036
- monte_carlo.md: mc_inputs, mc_histogram, mc_check_shape, mc_sampling_band,
  mc_realizations_needed, lhs_vs_random, lhs_segments, lhs_bins, lhs_repeat,
  mc_convergence_trace, reliability_mc_vp034, rs_surface_fit

The Studio captures on fem.md come from tools/capture_reliability_fem.py.
"""

import contextlib
import copy
import io
import math
import os
import pickle
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings("ignore")

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(_HERE)
OUT = os.path.join(REPO, "docs", "reliability", "images")

LEM11 = os.path.join(REPO, "docs", "lem", "files", "xslope_reliability.xlsx")
ROCSCIENCE = os.path.join(REPO, "docs", "verification", "files", "rocscience")
VP34 = os.path.join(ROCSCIENCE, "vp034.xlsx")
VP35 = os.path.join(ROCSCIENCE, "vp035.xlsx")
VP36 = os.path.join(ROCSCIENCE, "vp036.xlsx")

DPI = 200
BLUE = "#1f4e9c"
ORANGE = "#d9822b"
RED = "#b03a2e"
GREEN = "#2e7d5b"
GREY = "0.45"
LIGHT = "#c9d5e3"      # histogram bars
SKY = "#dce6f2"        # the area above F = 1
# the house style of the computed figures; the figures drawn by XSLOPE's own plot
# functions (Studio's) keep theirs
STYLE = {"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False,
         "legend.frameon": False}
# label placement, in points on the page: every label is placed from the
# geometry it labels and clears every line, curve, marker and label by GAP; a
# pointer line starts GAP from its label and stops POINTER_GAP short of its target
GAP = 3.0
POINTER_GAP = 4.0
POINTER_LEN = 12.0

CACHE = None           # --cache DIR
_RUNS = {}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def _save(fig, name, tight=True):
    _grid_backs(fig)
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    if tight:
        fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.06,
                    facecolor="white")
    else:
        fig.savefig(path, dpi=DPI, facecolor="white")
    plt.close(fig)
    print("wrote", os.path.relpath(path, REPO))


def _run(name, fn):
    """A run's result, made once per invocation (and kept in --cache DIR)."""
    if name in _RUNS:
        return _RUNS[name]
    path = os.path.join(CACHE, name + ".pkl") if CACHE else None
    if path and os.path.exists(path):
        with open(path, "rb") as f:
            _RUNS[name] = pickle.load(f)
        return _RUNS[name]
    print(f"  running {name} ...", flush=True)
    out = fn()
    _RUNS[name] = out
    if path:
        os.makedirs(CACHE, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump(out, f)
    return out


def norm_pdf(x, m, s):
    x = np.asarray(x, float)
    return np.exp(-((x - m) / s) ** 2 / 2) / (s * math.sqrt(2 * math.pi))


def lognorm_pdf(x, mean, sd):
    """The lognormal density with mean ``mean`` and standard deviation ``sd``."""
    z2 = math.log(1 + (sd / mean) ** 2)
    lam = math.log(mean) - z2 / 2
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    ok = x > 0
    out[ok] = (np.exp(-(np.log(x[ok]) - lam) ** 2 / (2 * z2))
               / (x[ok] * math.sqrt(2 * math.pi * z2)))
    return out


def beta_ln(F, cov):
    return math.log(F / math.sqrt(1 + cov * cov)) / math.sqrt(math.log(1 + cov * cov))


def Phi(z):
    return 0.5 * math.erfc(-z / math.sqrt(2))


def _no_y(ax):
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def _ink(s, size=None):
    """A label's ink extents (x0, y0, x1, y1), in points from the start of its
    baseline (y0 < 0 is ink below the baseline, a subscript's)."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath
    size = size or plt.rcParams["font.size"]
    return TextPath((0, 0), s, size=size, prop=FontProperties(size=size)).get_extents().extents


def _to_pt(ax):
    """Data -> points on the page, and back (the axes' limits must be set)."""
    k = 72 / ax.figure.dpi

    def fwd(xy):
        return ax.transData.transform(np.atleast_2d(xy)) * k

    def inv(p):
        return ax.transData.inverted().transform(np.atleast_2d(p) / k)
    return fwd, inv


def _top_labels(ax, top, labels):
    """Labels (x, text, color) centered over the tops of vertical lines that end
    at y = top, on one baseline, the lowest ink of any GAP above the lines."""
    base = GAP - min(_ink(s)[1] for _, s, _ in labels)
    for x, s, col in labels:
        ax.annotate(s, xy=(x, top), xytext=(0, base), textcoords="offset points",
                    ha="center", va="baseline", color=col)


def _on_grid(t):
    """Mark a label that sits over grid lines: at save, _grid_backs puts a white
    box GAP around its ink, over the grid and under every line, marker and label,
    so the grid stops GAP short of the label and no data is hidden (the label
    must not sit on a shaded area, which the box would also cover)."""
    t._on_grid = True
    return t


def _grid_backs(fig):
    marked = [t for ax in fig.axes for t in ax.texts
              if getattr(t, "_on_grid", False) and not getattr(t, "_backed", False)]
    if not marked:
        return
    from matplotlib.patches import Rectangle
    dpi = fig.dpi
    fig.set_dpi(DPI)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    g = GAP * DPI / 72
    for t in marked:
        arrow = getattr(t, "arrow_patch", None)
        if arrow is not None:
            arrow.set_visible(False)              # the label's ink only
        r.clear()
        t.draw(r)
        if arrow is not None:
            arrow.set_visible(True)
        a = np.asarray(r.buffer_rgba())[..., 3]
        ys, xs = np.nonzero(a)
        assert len(xs), f"{t.get_text()!r} is off the canvas"
        H = a.shape[0]
        (x0, y0), (x1, y1) = t.axes.transData.inverted().transform(
            [(xs.min() - g, H - ys.max() - 1 - g), (xs.max() + 1 + g, H - ys.min() + g)])
        t.axes.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="white", ec="none",
                                   zorder=1.6, gid="grid_back"))
        t._backed = True
    r.clear()
    fig.set_dpi(dpi)


def _f1(ax, top, label=True):
    """The line F = 1, labeled at its top."""
    ax.plot([1, 1], [0, top], color=RED, lw=1.5, zorder=4, solid_capstyle="butt")
    if label:
        _top_labels(ax, top, [(1, "F = 1", RED)])


# ---------------------------------------------------------------------------
# the runs
# ---------------------------------------------------------------------------

def run_lem11_taylor():
    from xslope.fileio import load_slope_data
    from xslope.reliability import reliability_taylor

    def go():
        ok, r = _quiet(reliability_taylor, load_slope_data(LEM11), "spencer",
                       num_slices=40)
        assert ok, r
        return {k: r[k] for k in ("F_MLV", "sigma_F", "COV_F", "beta_ln",
                                  "beta_normal", "prob_failure", "pf_normal",
                                  "param_info")}
    return _run("lem11_taylor", go)


def run_lem11_mc():
    from xslope.fileio import load_slope_data
    from xslope.reliability import reliability_mc

    def go():
        ok, r = _quiet(reliability_mc, load_slope_data(LEM11), "spencer",
                       n_samples=10000, search=True, num_slices=40)
        assert ok, r
        out = dict(r)
        out["fs"] = np.asarray(r["fs_samples"], float)
        out["x"] = np.asarray(r["param_samples"], float)
        return out
    return _run("lem11_mc", go)


def run_lem11_rs():
    from xslope.fileio import load_slope_data
    from xslope.reliability import reliability_rs, _ccd_design

    def go():
        ok, r = _quiet(reliability_rs, load_slope_data(LEM11), "spencer",
                       search=True, num_slices=40)
        assert ok, r
        mlv = np.array([p["mlv"] for p in r["param_info"]])
        std = np.array([p["std"] for p in r["param_info"]])
        design = mlv + _ccd_design(len(mlv), r["alpha"]) * std
        return {"names": [p["param"] for p in r["param_info"]], "mlv": mlv,
                "std": std, "coef": np.asarray(r["rs_coefficients"], float),
                "design": design, "x": np.asarray(r["param_samples"], float),
                "pf": r["pf_empirical"], "beta_ln": r["beta_ln"],
                "gate_r2": r["gate_r2"], "n_design": r["n_design"]}
    return _run("lem11_rs", go)


def run_vp035():
    """Cannon Dam: the circle of lowest F (a Bishop search at the means, 50
    slices, as the VP35-fs tag runs it) and the file's own circle (the minimum-beta
    circle of the verification page), each with the Taylor series held on it."""
    from xslope.fileio import load_slope_data
    from xslope.search import circular_search
    from xslope.reliability import reliability_taylor
    from xslope.slice import generate_slices

    def go():
        sd = load_slope_data(VP35)
        cache, converged, _, _ = _quiet(circular_search, sd, "bishop",
                                        num_slices=50)
        a = cache[0]
        circles = {"A": {"Xo": a["Xo"], "Yo": a["Yo"], "Depth": a["Depth"],
                         "R": a["Yo"] - a["Depth"]},
                   "B": dict(sd["circles"][0])}
        out = {"search_F": a["FS"], "circles": {}}
        for key, circ in circles.items():
            s = copy.deepcopy(sd)
            s["circles"] = [circ]
            ok, r = _quiet(reliability_taylor, s, "bishop", search=False)
            assert ok, r
            ok, res = generate_slices(s, circle=circ, num_slices=40)
            assert ok, res
            out["circles"][key] = {
                "circle": {k: float(v) for k, v in circ.items()},
                "arc": np.asarray(res[1].coords, float),
                "F": float(r["F_MLV"]), "beta_ln": float(r["beta_ln"])}
        names = [m["name"] for m in sd["materials"]]
        out["zones"] = [(names[p["mat_id"]], np.asarray(p["polygon"].exterior.coords))
                        for p in sd["polygons"]]
        out["ground"] = np.asarray(sd["ground_surface"].coords)
        return out
    return _run("vp035", go)


# ---------------------------------------------------------------------------
# index.md
# ---------------------------------------------------------------------------

def same_f():
    """Two slopes at F = 1.5 whose F is known to +/-0.15 and +/-0.45 (one
    standard deviation), F lognormal: the area below F = 1 is P_f."""
    fig, axs = plt.subplots(1, 2, figsize=(8.0, 2.7), sharey=True)
    x = np.linspace(0.3, 3.0, 1200)
    for ax, sd, name in zip(axs, (0.15, 0.45), ("Slope A", "Slope B")):
        y = lognorm_pdf(x, 1.5, sd)
        pf = Phi(-beta_ln(1.5, sd / 1.5))
        ax.fill_between(x, y, where=x >= 1, color=SKY, lw=0)
        ax.fill_between(x, y, where=x <= 1, color=ORANGE, alpha=0.75, lw=0)
        ax.plot(x, y, color=BLUE, lw=2)
        _f1(ax, 2.95)
        ax.set_xlim(0.3, 3.0)
        ax.set_ylim(0, 3.3)
        _no_y(ax)
        ax.set_xlabel("Factor of safety, F")
        p = f"{pf * 100:.3f}%" if pf < 0.001 else f"{pf * 100:.0f}%"
        ax.text(1.95, 2.55, name, fontweight="bold", va="center")
        ax.text(1.95, 2.1, f"F = 1.5 ± {sd:.2f}", va="center")
        ax.text(1.95, 1.65, f"$P_f$ = {p}", va="center", color=RED)
    _save(fig, "same_f.png")


def r_pf():
    """R = 1 - P_f on a distribution of F."""
    fig, ax = plt.subplots(figsize=(5.2, 2.9))
    x = np.linspace(0.2, 2.8, 1000)
    y = lognorm_pdf(x, 1.3, 0.4)
    ax.fill_between(x, y, where=x >= 1, color=SKY, lw=0)
    ax.fill_between(x, y, where=x <= 1, color=ORANGE, alpha=0.75, lw=0)
    ax.plot(x, y, color=BLUE, lw=2)
    _f1(ax, 1.3)
    ax.text(0.86, 0.18, "$P_f$", fontsize=13, ha="center", va="center")
    ax.text(1.36, 0.3, "$R = 1 - P_f$", fontsize=13, ha="center", va="center")
    ax.set_xlim(0.2, 2.8)
    ax.set_ylim(0, 1.45)
    _no_y(ax)
    ax.set_xlabel("Factor of safety, F")
    _save(fig, "r_pf.png")


def sigma_normal():
    """One standard deviation on a normal distribution: 68% of the values."""
    fig, ax = plt.subplots(figsize=(5.2, 2.8))
    z = np.linspace(-3.6, 3.6, 800)
    y = norm_pdf(z, 0, 1)
    ax.fill_between(z, y, where=np.abs(z) <= 1, color=SKY, lw=0)
    ax.plot(z, y, color=BLUE, lw=2)
    for v in (-1, 1):
        ax.plot([v, v], [0, norm_pdf(v, 0, 1)], color=GREY, lw=1, ls="--")
    ax.plot([0, 0], [0, 0.13], color=GREY, lw=1)
    ax.text(0, 0.21, "68% of\nthe values", ha="center", va="center")
    ax.annotate("", xy=(1, 0.07), xytext=(0, 0.07),
                arrowprops=dict(arrowstyle="<->", color="0.15", lw=1.1,
                                shrinkA=0, shrinkB=0))
    ax.text(0.5, 0.08, "σ", ha="center", va="bottom", fontsize=12)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_xticklabels([r"$\bar{x} - 2σ$", r"$\bar{x} - σ$", r"$\bar{x}$",
                        r"$\bar{x} + σ$", r"$\bar{x} + 2σ$"])
    ax.set_xlim(-3.6, 3.6)
    ax.set_ylim(0, 0.44)
    _no_y(ax)
    _save(fig, "sigma_normal.png")


def sigma_rules():
    """The 3-sigma and N-sigma rules on LEM-11's cohesion: highest conceivable
    600 psf, lowest 200 psf, most likely 400 psf."""
    lcv, mlv, hcv = 200.0, 400.0, 600.0
    s4, s6 = (hcv - lcv) / 4, (hcv - lcv) / 6
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    x = np.linspace(50, 750, 900)
    y4 = norm_pdf(x, mlv, s4)
    ax.fill_between(x, y4, where=(x >= lcv) & (x <= hcv), color=SKY, lw=0)
    ax.plot(x, norm_pdf(x, mlv, s6), color=GREY, lw=1.4, ls="--",
            label=f"3σ rule: σ = (600 − 200)/6 = {s6:.0f} psf")
    ax.plot(x, y4, color=BLUE, lw=2,
            label=f"$N_σ$ rule: σ = (600 − 200)/4 = {s4:.0f} psf")
    top = norm_pdf(mlv, mlv, s6) * 1.08
    for v in (lcv, hcv):
        ax.plot([v, v], [0, top], color="0.25", lw=1, ls=":")
    _top_labels(ax, top, [(lcv, "LCV = 200", "black"), (hcv, "HCV = 600", "black")])
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1,
              handlelength=2.4, borderaxespad=0)
    ax.set_xlim(50, 750)
    ax.set_ylim(0, top * 1.16)
    # the +/-2 sigma arrow, its heads (half-width 0.2 of the text size, the "<->"
    # style's) GAP below the curve where the curve is lowest over it, at LCV and HCV
    fwd, inv = _to_pt(ax)
    lw = 1.1
    head = 0.2 * plt.rcParams["font.size"] + lw / 2
    (xl, yc), = fwd([(lcv, norm_pdf(lcv, mlv, s4))])
    (_, ya), = inv([(xl, yc - 1.0 - GAP - head)])
    ax.annotate("", xy=(lcv, ya), xytext=(hcv, ya),
                arrowprops=dict(arrowstyle="<->", color="0.15", lw=lw,
                                shrinkA=0, shrinkB=0))
    s = "±2σ: 95%"
    ax.annotate(s, xy=(mlv, ya), xytext=(0, lw / 2 + GAP - _ink(s)[1]),
                textcoords="offset points", ha="center", va="baseline")
    _no_y(ax)
    ax.set_xticks([100, 200, 300, 400, 500, 600, 700])
    ax.set_xlabel("Cohesion, c (psf)")
    _save(fig, "sigma_rules.png")


def beta_on_f():
    """beta on a normal distribution of F: beta sigma_F from F = 1 to F_MLV."""
    F, sd = 1.35, 0.25
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    x = np.linspace(0.4, 2.3, 1000)
    y = norm_pdf(x, F, sd)
    pk = norm_pdf(F, F, sd)
    ax.fill_between(x, y, where=x >= 1, color=SKY, lw=0)
    ax.fill_between(x, y, where=x <= 1, color=ORANGE, alpha=0.75, lw=0)
    ax.plot(x, y, color=BLUE, lw=2)
    top = pk * 1.1
    _f1(ax, top, label=False)
    ax.plot([F, F], [0, top], color="0.15", lw=1.3, solid_capstyle="butt")
    _top_labels(ax, top, [(1, "F = 1", RED), (F, "$F_{MLV}$", "0.15")])
    ym = pk * 0.55
    yl = pk * 0.22
    ax.plot([F + sd, F + sd], [0, norm_pdf(F + sd, F, sd)], color=GREY, lw=1,
            ls="--")
    for a, b, ya in ((1, F, ym), (F, F + sd, yl)):
        ax.annotate("", xy=(a, ya), xytext=(b, ya),
                    arrowprops=dict(arrowstyle="<->", color="0.15", lw=1.1,
                                    shrinkA=0, shrinkB=0))
    ax.set_xlim(0.4, 2.3)
    ax.set_ylim(0, top * 1.15)
    _no_y(ax)
    ax.set_xlabel("Factor of safety, F")

    # the labels, placed from the geometry (in points on the page)
    from shapely.geometry import LineString, Point, Polygon, box
    fwd, inv = _to_pt(ax)
    size = 12
    curve = LineString(fwd(np.column_stack([x, y])))
    half = {"curve": 1.0, "F1": 0.75, "MLV": 0.65, "dash": 0.5, "arrow": 0.55}
    lines = {"curve": curve,
             "F1": LineString(fwd([(1, 0), (1, top)])),
             "MLV": LineString(fwd([(F, 0), (F, top)])),
             "dash": LineString(fwd([(F + sd, 0), (F + sd, norm_pdf(F + sd, F, sd))]))}

    def under(s, a, b, y):
        """A dimension label under its arrow (from a to b at height y): its ink
        top GAP below the arrow, centered under it, or moved along it just far
        enough to clear the curve and the lines by GAP."""
        x0, y0, x1, y1 = _ink(s, size)
        (pa, py), (pb, _) = fwd([(a, y), (b, y)])
        base = py - half["arrow"] - GAP - y1
        w = x1 - x0

        def ink(cx):
            return box(cx - w / 2, base + y0, cx + w / 2, base + y1)

        def clear(cx):
            return min(ink(cx).distance(ln) - half[k] for k, ln in lines.items())
        steps = np.linspace(0, (pb - pa) / 2, 400)
        cx = next((pa + pb) / 2 + d * sgn for d in steps for sgn in (1, -1)
                  if clear((pa + pb) / 2 + d * sgn) >= GAP)
        assert pa <= cx - w / 2 and cx + w / 2 <= pb, s
        (xd, yd), = inv([(cx - (x0 + x1) / 2, base)])
        ax.text(xd, yd, s, ha="left", va="baseline", fontsize=size)
        return ink(cx)
    labels = [under("$\\beta\\sigma_F$", 1, F, ym), under("$σ_F$", F, F + sd, yl)]

    # P_f outside the orange area, its pointer aimed at the area's centroid along
    # the curve's normal and stopping POINTER_GAP short of the curve
    left = x <= 1
    area = Polygon(np.vstack([fwd(np.column_stack([x[left], y[left]])),
                              fwd([(1, 0), (x[left][0], 0)])]))
    c = np.array(area.centroid.coords[0])
    branch = LineString(fwd(np.column_stack([x[left], y[left]])))
    q = np.array(branch.interpolate(branch.project(Point(c))).coords[0])
    u = (q - c) / np.linalg.norm(q - c)
    p1 = q + u * (half["curve"] + POINTER_GAP)
    p0 = p1 + u * POINTER_LEN
    x0, y0, x1, y1 = _ink("$P_f$", size)
    w, h = x1 - x0, y1 - y0

    def pf_ink(t):
        cx, cy = p0 + u * t
        return box(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
    lo, hi = 0.0, w + h
    for _ in range(60):                              # pf_ink(t) GAP from p0
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if pf_ink(mid).distance(Point(p0)) < GAP else (lo, mid)
    pf = pf_ink(hi)
    cx, cy = pf.centroid.coords[0]
    (xd, yd), = inv([(cx - (x0 + x1) / 2, cy - (y0 + y1) / 2)])
    ax.text(xd, yd, "$P_f$", ha="left", va="baseline", fontsize=size)
    (ax0, ay0), (ax1, ay1) = inv([p0, p1])
    ax.plot([ax0, ax1], [ay0, ay1], color="0.15", lw=0.9, label="pointer")
    seg = LineString([p0, p1])
    assert abs(seg.distance(curve) - half["curve"] - POINTER_GAP) < 0.01
    assert seg.distance(area) > POINTER_GAP
    for ink_box in labels + [pf]:
        assert min(ink_box.distance(ln) - half[k] for k, ln in lines.items()) >= GAP - 1e-6
    assert pf.distance(seg) >= GAP - 1e-6 and not pf.intersects(area)
    _save(fig, "beta_on_f.png")


def beta_pf():
    """P_f = Phi(-beta): beta = 1, 2, 3 marked, and P_f = 1% at beta = 2.33."""
    b = np.linspace(0, 4, 400)
    p = np.array([Phi(-v) for v in b])
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    ax.semilogy(b, p, color=BLUE, lw=2)
    ax.set_xlim(0, 4)
    ax.set_ylim(1e-5, 1)
    ax.set_yticks([1, 0.1, 0.01, 0.001, 1e-4, 1e-5])
    ax.set_yticklabels(["100%", "10%", "1%", "0.1%", "0.01%", "0.001%"])
    ax.minorticks_off()
    ax.grid(True, color="0.9")
    ax.set_xlabel("Reliability index, β")
    ax.set_ylabel("Probability of failure, $P_f$")
    for v, txt in ((1, "16%"), (2, "2.3%"), (3, "0.13%")):
        pv = Phi(-v)
        ax.plot(v, pv, "o", color=BLUE, ms=6, zorder=5)
        _on_grid(ax.text(v + 0.09, pv * 1.25, txt, va="bottom"))
    bb = 2.326
    ax.plot([0, bb, bb], [0.01, 0.01, 1e-5], color=ORANGE, lw=1.3, ls="--")
    ax.plot(bb, 0.01, "o", color=ORANGE, ms=6, zorder=5)
    _on_grid(ax.text(0.08, 0.0125, "$P_f$ = 1% at β = 2.33", color=ORANGE,
                     va="bottom"))
    _save(fig, "beta_pf.png")


# Cannon Dam's zones (vp035.xlsx material names), in soil tints
VP035_ZONES = {"Phase II clay fill": "#f4e8d2", "Phase I clay fill": "#e3cfa0",
               "Sand filter": "#fbf8f1", "Spoil fill": "#efe6d0",
               "Foundation sand": "#d8ccb2", "Limestone": "#d5dadf"}
# zone names (m), each placed where it clears every line; the sand filter, a
# 5 m strip, is named from the upstream fill with a pointer line to it
VP035_NAMES = {"Phase II clay fill": (318.9, 173.2),
               "Phase I clay fill": (262.0, 156.7),
               "Foundation sand": (423.2, 148.7), "Limestone": (338.0, 134.1)}
VP035_FILTER = (250.5, 172.5)


def vp035_surfaces():
    """Cannon Dam at true scale: the downstream half of the section, the circle
    of lowest F and the circle of lowest beta, keyed with each one's F and beta.
    Every name is checked to clear every line, and the filter's pointer line to
    stop 4 pt short of the filter's face and clear everything else."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon as MPoly
    from shapely.geometry import LineString, Polygon, box
    from shapely.ops import unary_union
    V = run_vp035()
    W, H = 7.5, 3.3
    x0, x1, yb = 214.0, 458.0, 128.0
    yt = yb + (x1 - x0) * H / W
    pt = (x1 - x0) / W / 72                    # meters per point
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0, x1)
    ax.set_ylim(yb, yt)
    ax.set_aspect("equal")
    ax.axis("off")
    for name, poly in V["zones"]:
        ax.add_patch(MPoly(poly, closed=True, fc=VP035_ZONES[name], ec="#a3967a",
                           lw=0.7, zorder=1))
    g = V["ground"]
    ax.plot(g[:, 0], g[:, 1], color="0.15", lw=1.3, zorder=5)
    lines = [LineString(p) for _, p in V["zones"]] + [LineString(g)]
    section = ax.add_patch(MPoly(
        list(unary_union([Polygon(p) for _, p in V["zones"]]).exterior.coords),
        closed=True, fc="none", ec="none"))
    handles = []
    for k, col, name in (("A", BLUE, "Lowest F"), ("B", ORANGE, "Lowest β")):
        c = V["circles"][k]
        ci, arc = c["circle"], c["arc"]
        x = np.linspace(arc[0, 0], arc[-1, 0], 400)
        y = ci["Yo"] - np.sqrt(np.maximum(ci["R"] ** 2 - (x - ci["Xo"]) ** 2, 0))
        ax.plot(x, y, color=col, lw=2.4, zorder=4, clip_path=section)
        lines.append(LineString(np.column_stack([x, y])))
        handles.append(Line2D([], [], color=col, lw=2.4,
                              label=f"{name}:  F = {c['F']:.2f},  β = {c['beta_ln']:.2f}"))
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.985, 0.97),
              fontsize=11, handlelength=1.8, borderaxespad=0.0)
    kw = dict(fontsize=9.5, style="italic", color="0.15", va="center", zorder=5)
    names = [ax.text(x, y, n, ha="center", **kw) for n, (x, y) in VP035_NAMES.items()]
    names.append(ax.text(*VP035_FILTER, "Sand filter", ha="right", **kw))
    fig.canvas.draw()
    r, inv = fig.canvas.get_renderer(), ax.transData.inverted()

    def bbox(artist):
        (a0, b0), (a1, b1) = inv.transform(artist.get_window_extent(r).get_points())
        return box(a0, b0, a1, b1)

    every = unary_union(lines)
    boxes = [bbox(t) for t in names] + [bbox(ax.get_legend())]
    for t, b in zip(names + [ax.get_legend()], boxes):
        assert b.distance(every) > 3 * pt, getattr(t, "get_text", lambda: "key")()
    f = next(p for n, p in V["zones"] if n == "Sand filter")
    fa, fb = sorted(sorted(map(tuple, f[:-1]), key=lambda q: q[0])[:2],
                    key=lambda q: q[1])            # the upstream face
    y = VP035_FILTER[1]
    xf = fa[0] + (fb[0] - fa[0]) * (y - fa[1]) / (fb[1] - fa[1])
    nx = abs(fb[1] - fa[1]) / math.hypot(fb[0] - fa[0], fb[1] - fa[1])
    p0, p1 = (boxes[-2].bounds[2] + 3 * pt, y), (xf - 4 * pt / nx, y)
    ax.plot(*zip(p0, p1), color="0.15", lw=0.9, zorder=5)
    seg, face = LineString([p0, p1]), LineString([fa, fb])
    assert abs(seg.distance(face) - 4 * pt) < 0.01 * pt
    others = unary_union([ln.difference(face.buffer(0.01)) for ln in lines])
    assert seg.distance(others) > 3 * pt, "the filter's pointer meets a line"
    _save(fig, "vp035_surfaces.png", tight=False)


# ---------------------------------------------------------------------------
# taylor.md
# ---------------------------------------------------------------------------

def _taylor_rows():
    T = run_lem11_taylor()
    p = {q["param"]: q for q in T["param_info"]}
    return T, p["gamma"], p["c"]


def taylor_gamma():
    """F at the unit weight's three Taylor-series values, and the straight
    line from F- to F+ the method assumes: it passes above F_MLV."""
    T, g, _ = _taylor_rows()
    x = [g["mlv"] - g["std"], g["mlv"], g["mlv"] + g["std"]]
    y = [g["F_minus"], T["F_MLV"], g["F_plus"]]
    mid = (y[0] + y[2]) / 2
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    ax.plot([x[0], x[2]], [y[0], y[2]], color=GREY, lw=1.2, ls="--", zorder=2)
    ms = 7
    ax.plot(x, y, "o", ms=ms, color=BLUE, zorder=4)
    # the gap between the line and F_MLV, at the most likely value
    xd = x[1] + 2.4
    ax.annotate("", xy=(xd, mid), xytext=(xd, y[1]),
                arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.1,
                                shrinkA=0, shrinkB=0))
    ax.plot([x[1] + 0.9, xd + 0.6], [y[1], y[1]], color=ORANGE, lw=0.8)
    ax.plot([x[1] + 1.0, xd + 0.6], [mid, mid], color=ORANGE, lw=0.8)
    ax.text(xd + 0.9, (mid + y[1]) / 2, f"{mid - y[1]:.3f}", color=ORANGE,
            ha="left", va="center")
    ax.set_xlim(98, 136)
    ax.set_ylim(1.12, 1.64)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{x[0]:g}\nMLV − σ", f"{x[1]:g}\nMLV", f"{x[2]:g}\nMLV + σ"])
    ax.set_xlabel("Unit weight, γ (pcf)")
    ax.set_ylabel("Factor of safety, F")
    ax.grid(True, color="0.9")
    # each point's label to its left, at its height, all three the same distance
    # from their points: the ink GAP from the marker, or farther if a label's
    # grid box (GAP around its ink) would otherwise reach the dashed line
    from shapely.geometry import LineString, box
    fwd, inv = _to_pt(ax)
    line = LineString(fwd([(x[0], y[0]), (x[2], y[2])]))
    labels = []
    for xi, yi, lab in zip(x, y, ("$F^-$", "$F_{MLV}$", "$F^+$")):
        s = f"{lab} = {yi:.3f}"
        (px, py), = fwd([(xi, yi)])
        labels.append((s, _ink(s), px, py))

    def grid_box(ink, px, py, off):
        x0, y0, x1, y1 = ink
        r = px - off
        return box(r - (x1 - x0) - GAP, py - (y1 - y0) / 2 - GAP, r + GAP,
                   py + (y1 - y0) / 2 + GAP)
    off = ms / 2 + GAP
    while any(grid_box(ink, px, py, off).distance(line) < 1.2 / 2
              for _, ink, px, py in labels):
        off += 0.1
    for s, (x0, y0, x1, y1), px, py in labels:
        (xt, yt), = inv([(px - off - x1, py - (y0 + y1) / 2)])
        _on_grid(ax.text(xt, yt, s, ha="left", va="baseline"))
    _save(fig, "taylor_gamma.png")


def taylor_shapes():
    """The two shapes the Taylor series can assume for F, each with LEM-11's
    F_MLV and sigma_F, and the area below F = 1 under each."""
    T = run_lem11_taylor()
    F, sd = T["F_MLV"], T["sigma_F"]
    x = np.linspace(-0.3, 3.0, 1400)
    fig, axs = plt.subplots(1, 2, figsize=(8.0, 2.8), sharey=True)
    for ax, kind in zip(axs, ("normal", "lognormal")):
        y = norm_pdf(x, F, sd) if kind == "normal" else lognorm_pdf(x, F, sd)
        b = T["beta_normal"] if kind == "normal" else T["beta_ln"]
        pf = T["pf_normal"] if kind == "normal" else T["prob_failure"]
        ax.fill_between(x, y, where=x >= 1, color=SKY, lw=0)
        ax.fill_between(x, y, where=x <= 1, color=ORANGE, alpha=0.75, lw=0)
        ax.plot(x, y, color=BLUE, lw=2)
        ax.axhline(0, color="0.15", lw=0.8)
        _f1(ax, 1.27)
        ax.set_xlim(-0.3, 3.0)
        ax.set_ylim(0, 1.42)
        _no_y(ax)
        ax.spines["bottom"].set_visible(False)
        ax.set_xlabel("Factor of safety, F")
        ax.text(2.0, 1.05, f"F {kind}", fontweight="bold", va="center")
        bname = "$β_N$" if kind == "normal" else "$β_{LN}$"
        ax.text(2.0, 0.84, f"{bname} = {b:.3f}", va="center")
        ax.text(2.0, 0.63, f"$P_f$ = {pf * 100:.1f}%", va="center", color=RED)
    _save(fig, "taylor_shapes.png")


def reliability_taylor_vp036():
    """VP36's Taylor series, searched at each solve, drawn by Studio's own plot."""
    from xslope.fileio import load_slope_data
    from xslope.reliability import reliability_taylor
    from xslope.plot import plot_reliability_results

    sd = load_slope_data(VP36)

    def go():
        ok, r = _quiet(reliability_taylor, sd, "bishop")
        assert ok, r
        return r
    r = _run("vp036_taylor", go)
    fig = plt.figure(figsize=(9, 5.2))
    _quiet(plot_reliability_results, sd, r, fig=fig, frame="content")
    _save(fig, "reliability_taylor_vp036.png")


# ---------------------------------------------------------------------------
# monte_carlo.md
# ---------------------------------------------------------------------------

def _mc_columns(M):
    names = [p["param"] for p in M["param_info"]]
    return M["x"][:, names.index("gamma")], M["x"][:, names.index("c")]


def mc_inputs():
    """LEM-11's two input distributions, and realization 1's draw from each."""
    M = run_lem11_mc()
    g, c = _mc_columns(M)
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 2.4))
    for ax, (m, s, lab, v, fmt) in zip(axs, (
            (120, 8, "Unit weight, γ (pcf)", g[0], "γ = {:.1f}"),
            (400, 100, "Cohesion, c (psf)", c[0], "c = {:.0f}"))):
        x = np.linspace(m - 3.6 * s, m + 3.6 * s, 400)
        y = norm_pdf(x, m, s)
        ax.fill_between(x, y, color=SKY, lw=0)
        ax.plot(x, y, color=BLUE, lw=2)
        yv = norm_pdf(v, m, s)
        ax.plot([v, v], [0, yv], color=ORANGE, lw=1.8)
        ax.plot([v], [yv], "o", ms=5, color=ORANGE)
        right = v >= m
        ax.text(v + (1 if right else -1) * 0.2 * s, yv * 1.04, fmt.format(v),
                ha="left" if right else "right", va="bottom", color=ORANGE)
        _no_y(ax)
        ax.set_xticks([m - 2 * s, m - s, m, m + s, m + 2 * s])
        ax.set_xlim(m - 3.6 * s, m + 3.6 * s)
        ax.set_ylim(0, 0.5 / s)
        ax.set_xlabel(lab)
    fig.tight_layout(w_pad=2.5)
    _save(fig, "mc_inputs.png")


def _bins(F, w=0.05):
    """Bin edges on multiples of w, so F = 1 is an edge."""
    return np.arange(math.floor(F.min() / w) * w, math.ceil(F.max() / w) * w + w / 2, w)


def _hist(ax, F, color=LIGHT, tail=True):
    h, e = np.histogram(F, bins=_bins(F), density=True)
    for i, v in enumerate(h):
        hot = tail and e[i + 1] <= 1.0 + 1e-9
        ax.bar(e[i], v, width=e[i + 1] - e[i], align="edge",
               color=ORANGE if hot else color, alpha=0.85 if hot else 1.0,
               ec="white", lw=0.5)
    return h.max()


def mc_histogram():
    """LEM-11's 10,000 Monte Carlo factors of safety: the realizations below
    F = 1 in orange, their mean beside F_MLV."""
    M = run_lem11_mc()
    T = run_lem11_taylor()
    F = M["fs"][np.isfinite(M["fs"])]
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    top = _hist(ax, F) * 1.08
    _f1(ax, top)
    n_fail = int((F < 1).sum())
    hi = top * 1.2
    ax.plot([M["mean_FS"]] * 2, [0, hi], color="0.15", lw=1.2, ls="--", zorder=5)
    ax.text(M["mean_FS"], hi, f" mean = {M['mean_FS']:.3f}", ha="left",
            va="bottom")
    ax.plot([T["F_MLV"]] * 2, [0, hi], color="0.15", lw=1.2, ls=":", zorder=5)
    ax.text(T["F_MLV"], hi, f"$F_{{MLV}}$ = {T['F_MLV']:.3f} ", ha="right",
            va="bottom")
    ax.text(0.45, top * 0.45, f"{n_fail:,} of {len(F):,}\nbelow F = 1",
            ha="center", va="bottom", color=RED)
    ax.set_xlim(0, 3.0)
    ax.set_ylim(0, hi * 1.1)
    _no_y(ax)
    ax.set_xlabel("Factor of safety, F")
    _save(fig, "mc_histogram.png")


def mc_check_shape():
    """The Monte Carlo histogram under the two curves the Taylor series assumes
    from its F_MLV and sigma_F."""
    M = run_lem11_mc()
    T = run_lem11_taylor()
    F = M["fs"][np.isfinite(M["fs"])]
    x = np.linspace(0, 3.0, 1200)
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    top = _hist(ax, F, tail=False) * 1.12
    ax.plot(x, norm_pdf(x, T["F_MLV"], T["sigma_F"]), color=GREEN, lw=1.8,
            ls="--", label=f"Taylor series, normal ($P_f$ = {T['pf_normal'] * 100:.1f}%)")
    ax.plot(x, lognorm_pdf(x, T["F_MLV"], T["sigma_F"]), color=BLUE, lw=2,
            label=f"Taylor series, lognormal ($P_f$ = {T['prob_failure'] * 100:.1f}%)")
    ax.bar([0], [0], color=LIGHT,
           label=f"Monte Carlo ($P_f$ = {M['pf_empirical'] * 100:.1f}%, counted)")
    _f1(ax, top)
    ax.set_xlim(0, 3.0)
    ax.set_ylim(0, top * 1.12)
    _no_y(ax)
    ax.set_xlabel("Factor of safety, F")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=9.5,
              borderaxespad=0)
    _save(fig, "mc_check_shape.png")


def mc_sampling_band():
    """The 95% band, P_f +/- 1.96 sqrt(p(1 - p)/n) around LEM-11's counted P_f,
    and the run's own running count as the realizations accumulate."""
    M = run_lem11_mc()
    F = M["fs"][np.isfinite(M["fs"])]
    p = float((F < 1).mean())
    n = np.linspace(100, len(F), 400)
    half = 1.96 * np.sqrt(p * (1 - p) / n) * 100
    k = np.arange(1, len(F) + 1)
    run = np.cumsum(F < 1) / k * 100
    fig, ax = plt.subplots(figsize=(6.6, 3.5))
    ax.fill_between(n, p * 100 - half, p * 100 + half, color=SKY, lw=0)
    ax.plot([0, len(F)], [p * 100] * 2, color=GREY, lw=1, ls="--")
    ax.plot(k[99:], run[99:], color=BLUE, lw=1.4)
    h = 1.96 * math.sqrt(p * (1 - p) / len(F)) * 100
    _on_grid(ax.annotate(f"±{h:.1f} points\nat n = {len(F):,}",
                xy=(len(F) * 0.97, p * 100 + h + 0.25),
                xytext=(len(F) * 0.97, p * 100 + h + 3.2), ha="right",
                va="bottom", arrowprops=dict(arrowstyle="-", color="0.3", lw=0.8,
                                             shrinkA=2, shrinkB=0)))
    ax.text(350, p * 100 + 6.0, "95% band", color=BLUE, va="center")
    _on_grid(ax.text(len(F) / 2, p * 100 - 1.6, f"$P_f$ = {p * 100:.2f}%",
                     ha="center", va="top"))
    ax.set_xlim(0, len(F))
    ax.set_ylim(5, 27)
    ax.set_xticks([0, 2500, 5000, 7500, 10000])
    ax.set_xticklabels(["0", "2,500", "5,000", "7,500", "10,000"])
    ax.set_yticks([5, 10, 15, 20, 25])
    ax.set_yticklabels(["5%", "10%", "15%", "20%", "25%"])
    ax.set_xlabel("Number of realizations, n")
    ax.set_ylabel("Running $P_f$")
    ax.grid(True, color="0.9")
    _save(fig, "mc_sampling_band.png")


def mc_realizations_needed():
    """Realizations to count a P_f to within +/-10% of itself (95%):
    n = 1.96^2 (1 - p) / (p 0.1^2)."""
    p = np.logspace(-4, -0.3, 300)
    n = 1.96 ** 2 * (1 - p) / (p * 0.1 ** 2)
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    ax.loglog(p * 100, n, color=BLUE, lw=2)
    ax.grid(True, which="major", color="0.9")
    for pv in (0.1, 0.01, 0.001):
        nv = 1.96 ** 2 * (1 - pv) / (pv * 0.01)
        nr = float(f"{nv:.2g}")
        ax.plot(pv * 100, nv, "o", color=ORANGE, ms=6, zorder=5)
        _on_grid(ax.text(pv * 100 * 1.25, nv * 1.15, f"{pv * 100:g}%: {nr:,.0f}",
                         va="bottom"))
    ax.set_xlim(0.01, 50)
    ax.set_ylim(100, 1e7)
    ax.set_xticks([0.01, 0.1, 1, 10])
    ax.set_xticklabels(["0.01%", "0.1%", "1%", "10%"])
    ax.minorticks_off()
    ax.set_xlabel("Probability of failure, $P_f$")
    ax.set_ylabel("Realizations, n")
    _save(fig, "mc_realizations_needed.png")


def lhs_vs_random():
    """250 draws of LEM-11's unit weight, random and Latin hypercube, through the
    one routine both sampling engines draw with (same seed), counted in 25
    segments of equal probability: each segment should hold 10."""
    from scipy.stats import norm
    from xslope.reliability import _draw_sample_matrix, MC_DEFAULT_SEED
    M = run_lem11_mc()
    info = M["param_info"]
    j = [p["param"] for p in info].index("gamma")
    n, nb = 250, 25
    fig, axs = plt.subplots(1, 2, figsize=(8.0, 2.8), sharey=True)
    for ax, mode, title in zip(axs, ("random", "lhs"),
                               ("Random", "Latin hypercube")):
        X, err = _draw_sample_matrix(info, n, np.random.default_rng(MC_DEFAULT_SEED),
                                     "normal", sampling=mode)
        assert err is None, err
        u = norm.cdf(X[:, j], loc=info[j]["mlv"], scale=info[j]["std"])
        counts = np.bincount(np.minimum((u * nb).astype(int), nb - 1), minlength=nb)
        ax.bar(np.arange(nb) + 0.5, counts, width=0.82, color=BLUE if mode == "lhs" else LIGHT,
               ec="none")
        ax.axhline(n / nb, color=ORANGE, lw=1.3, ls="--")
        ax.set_title(f"{title}: {n} draws", fontsize=10.5)
        ax.set_xlim(0, nb)
        ax.set_xticks([0, nb / 2, nb])
        ax.set_xticklabels(["lowest", "middle", "highest"])
        ax.set_xlabel(f"Unit weight, in {nb} segments of equal probability")
        ax.grid(True, axis="y", color="0.9")
        ax.set_axisbelow(True)                   # the grid behind the bars
    axs[0].set_ylabel("Draws in the segment")
    axs[1].annotate("10 per segment", xy=(nb - 0.3, n / nb),
                    xytext=(0, 1.3 / 2 + GAP - _ink("10 per segment")[1]),
                    textcoords="offset points", color=ORANGE, ha="right", va="baseline")
    fig.tight_layout(w_pad=2.0)
    _save(fig, "lhs_vs_random.png")


def run_lem11_lhs10():
    """LEM-11's Monte Carlo at ten Latin hypercube realizations, the default
    seed, searched once at the means: small enough to show every draw."""
    from xslope.fileio import load_slope_data
    from xslope.reliability import reliability_mc

    def go():
        ok, r = _quiet(reliability_mc, load_slope_data(LEM11), "spencer",
                       n_samples=10, search=True, num_slices=40, sampling="lhs")
        assert ok, r
        return {"names": [p["param"] for p in r["param_info"]],
                "param_info": r["param_info"],
                "fs": np.asarray(r["fs_samples"], float),
                "x": np.asarray(r["param_samples"], float),
                "rng_seed": r["rng_seed"]}
    return _run("lem11_lhs10", go)


def run_lem11_repeats():
    """Twenty 1,000-realization Monte Carlo runs each way, seeds 1 to 20, on the
    critical circle of a search at the means (held fixed, as Monte Carlo holds
    it)."""
    from xslope.fileio import load_slope_data
    from xslope.search import circular_search
    from xslope.reliability import reliability_mc

    def go():
        sd = load_slope_data(LEM11)
        cache, _, _, _ = _quiet(circular_search, sd, "spencer", num_slices=40)
        a = cache[0]
        sd["circles"] = [{"Xo": a["Xo"], "Yo": a["Yo"], "Depth": a["Depth"],
                          "R": a["Yo"] - a["Depth"]}]
        out = {}
        for mode in ("random", "lhs"):
            pf = []
            for seed in range(1, 21):
                ok, r = _quiet(reliability_mc, sd, "spencer", n_samples=1000,
                               search=False, num_slices=40, rng_seed=seed,
                               sampling=mode)
                assert ok, r
                pf.append(r["pf_empirical"])
            out[mode] = np.array(pf)
        return out
    return _run("lem11_repeats", go)


def _lhs_rows(L):
    """Each realization's segment and place in its segment, for each input,
    recovered from its draw: u = Phi((x - mean)/sd), segment = floor(10 u) + 1,
    position = 10 u - floor(10 u); checked to be a Latin hypercube (each input's
    ten segments used once)
    and to give the draws back."""
    from scipy.stats import norm
    n = len(L["fs"])
    out = {}
    for j, p in enumerate(L["param_info"]):
        xj = L["x"][:, j]
        u = norm.cdf((xj - p["mlv"]) / p["std"])
        b = np.floor(n * u).astype(int) + 1
        pos = n * u - np.floor(n * u)
        assert sorted(b) == list(range(1, n + 1)), (p["param"], b)
        back = p["mlv"] + p["std"] * norm.ppf((b - 1 + pos) / n)
        assert np.abs(back - xj).max() < 1e-9, (p["param"], np.abs(back - xj).max())
        out[p["param"]] = {"seg": b, "pos": pos, "x": xj, "mlv": p["mlv"], "std": p["std"]}
    return out


def lhs_segments():
    """Step 1 of a Latin hypercube, before anything is drawn: LEM-11's two inputs,
    each cut into the ten segments of equal probability of a ten-realization run,
    bounded by dashed red lines and numbered near the base."""
    _lhs_figure("lhs_segments.png", draws=False)


def lhs_bins():
    """The ten-realization run's draws on its segments (``_lhs_figure``)."""
    _lhs_figure("lhs_bins.png", draws=True)


def _lhs_figure(fname, draws=True):
    """How a Latin hypercube draws, in the style of the course slide: LEM-11's two
    inputs, each split into the ten segments of equal probability of a
    ten-realization run (LEM-11 at the default seed), bounded by dashed red lines
    and numbered near the base, and each realization's draw on the curve in its
    segment: every segment used once, at a random place inside it. Prints the
    realizations, for the page's table."""
    from scipy.stats import norm
    from shapely.geometry import LineString, Point, box
    L = run_lem11_lhs10()
    rows = _lhs_rows(L)
    n = len(L["fs"])
    if draws:
        print("  n | gamma: segment, position, value | c: segment, position, value | F")
    for i in range(n if draws else 0):
        g, c = rows["gamma"], rows["c"]
        print(f"  {i + 1} | {g['seg'][i]}, {float(g['pos'][i])!r}, {float(g['x'][i])!r} | "
              f"{c['seg'][i]}, {float(c['pos'][i])!r}, {float(c['x'][i])!r} | "
              f"{float(L['fs'][i])!r}")
    fig, axs = plt.subplots(2, 1, figsize=(7.2, 4.6))
    ms, lw_seg = 6, 1.1
    for ax, (name, xlabel) in zip(axs, (("gamma", "Unit weight, γ (pcf)"),
                                        ("c", "Cohesion, c (psf)"))):
        R = rows[name]
        m, s = R["mlv"], R["std"]
        pk = norm_pdf(m, m, s)
        top = pk * 1.1                           # the segment lines' top, over the curve
        x = np.linspace(m - 3.6 * s, m + 3.6 * s, 1200)
        y = norm_pdf(x, m, s)
        ax.fill_between(x, y, color=SKY, lw=0)
        for e in m + s * norm.ppf(np.arange(1, n) / n):
            ax.plot([e, e], [0, top], color=RED, lw=lw_seg, ls="--", zorder=3)
        ax.plot(x, y, color=BLUE, lw=2, zorder=4)
        if draws:
            ax.plot(R["x"], norm_pdf(R["x"], m, s), "o", ms=ms, color=ORANGE, mec="white",
                    mew=0.8, zorder=5)
        _no_y(ax)
        ax.set_xticks([m - 2 * s, m - s, m, m + s, m + 2 * s])
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(0, top)
        ax.set_xlabel(xlabel)
    fig.tight_layout(h_pad=1.0)

    # the segment numbers inside the shaded area near its base: on one baseline,
    # their lowest ink GAP above the axis, each centered in the room its segment
    # has at that height (between its lines, and inside the curve for the two
    # outer ones); each clears the lines, the curve, the draws and the others by GAP
    inks = [_ink(str(k)) for k in range(1, n + 1)]
    for ax, name in zip(axs, ("gamma", "c")):
        R = rows[name]
        m, s = R["mlv"], R["std"]
        fwd, inv = _to_pt(ax)
        lo, hi = ax.get_xlim()
        _, top = ax.get_ylim()
        xx = np.linspace(lo, hi, 1200)
        curve = LineString(fwd(np.column_stack([xx, norm_pdf(xx, m, s)])))
        edges = m + s * norm.ppf(np.arange(1, n) / n)
        segs = [LineString(fwd([(e, 0), (e, top)])) for e in edges]
        dots = ([Point(p) for p in fwd(np.column_stack([R["x"], norm_pdf(R["x"], m, s)]))]
                if draws else [])
        (x_lo, y_axis), (x_hi, _) = fwd([(lo, 0), (hi, 0)])
        base = y_axis + 0.8 / 2 + GAP - min(i[1] for i in inks)
        y_top = base + max(i[3] for i in inks)
        px = np.r_[x_lo, [sg.coords[0][0] for sg in segs], x_hi]
        # where the curve is GAP over the numbers' top, on each side
        cv = np.asarray(curve.coords)
        tall = cv[:, 1] >= y_top + GAP + 1.0
        px[0], px[-1] = cv[tall, 0].min(), cv[tall, 0].max()
        placed = []
        for k, (x0, y0, x1, y1) in enumerate(inks):
            cx = (px[k] + px[k + 1]) / 2
            b = box(cx - (x1 - x0) / 2, base + y0, cx + (x1 - x0) / 2, base + y1)
            gaps = ([b.distance(sg) - lw_seg / 2 for sg in segs] + [b.distance(curve) - 1.0]
                    + [b.distance(d) - ms / 2 - 0.8 / 2 for d in dots]
                    + [b.distance(o) for o in placed])
            assert min(gaps) >= GAP - 1e-6, (name, k + 1, min(gaps))
            placed.append(b)
            (xt, yt), = inv([(cx - (x0 + x1) / 2, base)])
            ax.text(xt, yt, str(k + 1), ha="left", va="baseline", color="0.15")
    _save(fig, fname)


def lhs_repeat():
    """P_f from twenty 1,000-realization runs with different seeds, random and
    Latin hypercube, beside the 10,000-realization answer."""
    R = run_lem11_repeats()
    M = run_lem11_mc()
    fig, ax = plt.subplots(figsize=(6.0, 2.4))
    for k, (mode, name, col) in enumerate((("random", "Random", LIGHT),
                                           ("lhs", "Latin hypercube", BLUE))):
        pf = R[mode] * 100
        y = 1 - k
        ax.plot(pf, np.full(len(pf), y), "o", ms=7, mfc=col, mec="0.25", mew=0.6,
                alpha=0.9, zorder=3)
        _on_grid(ax.text(pf.min() - 0.25, y, name, ha="right", va="center"))
        print(f"  {name}: P_f {pf.min():.2f}% to {pf.max():.2f}%, sd {pf.std(ddof=1):.3f} points")
    ax.plot([M["pf_empirical"] * 100] * 2, [-0.7, 1.55], color=GREY, lw=1, ls="--",
            zorder=1)
    s = f"{M['pf_empirical'] * 100:.2f}% from 10,000 realizations"
    _on_grid(ax.annotate(s, xy=(M["pf_empirical"] * 100, 1.55),
                         xytext=(0, GAP - _ink(s)[1]), textcoords="offset points",
                         color=GREY, ha="center", va="baseline"))
    ax.set_ylim(-0.7, 1.85)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    lo = min(R["random"].min(), R["lhs"].min()) * 100
    hi = max(R["random"].max(), R["lhs"].max()) * 100
    ax.set_xlim(lo - 3.2, hi + 0.6)
    ax.set_xlabel("$P_f$ from one run of 1,000 realizations (%)")
    ax.grid(True, axis="x", color="0.9")
    print(f"  scatter ratio (random / LHS): {R['random'].std(ddof=1) / R['lhs'].std(ddof=1):.2f}")
    _save(fig, "lhs_repeat.png")


def _surrogate(R, X):
    """The fitted quadratic at the physical points X (n x d)."""
    z = (np.asarray(X, float) - R["mlv"]) / R["std"]
    d = z.shape[1]
    cols = ([np.ones(len(z))] + [z[:, j] for j in range(d)]
            + [z[:, j] * z[:, k] for j in range(d) for k in range(j, d)])
    return np.column_stack(cols) @ R["coef"]


def rs_surface_fit():
    """LEM-11's response surface: the formula's F as contours over the two
    inputs, F = 1 in red, the real solves it was fitted to, and a sample of the
    draws; beside it, the formula against the real factor of safety of the
    first 1,000 Monte Carlo realizations, solved on the same circle."""
    R = run_lem11_rs()
    M = run_lem11_mc()
    ig, ic = R["names"].index("gamma"), R["names"].index("c")
    g = np.linspace(R["mlv"][ig] - 3.6 * R["std"][ig], R["mlv"][ig] + 3.6 * R["std"][ig], 241)
    c = np.linspace(max(0.0, R["mlv"][ic] - 3.6 * R["std"][ic]),
                    R["mlv"][ic] + 3.6 * R["std"][ic], 241)
    G, C = np.meshgrid(g, c)
    X = np.zeros((G.size, len(R["mlv"])))
    X[:, ig], X[:, ic] = G.ravel(), C.ravel()
    Fz = _surrogate(R, X).reshape(G.shape)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.6, 4.0),
                                 gridspec_kw={"width_ratios": [1.15, 1]})
    ax.contourf(G, C, Fz, levels=[-10, 1, 100], colors=[ORANGE, SKY], alpha=0.45)
    cs = ax.contour(G, C, Fz, levels=[0.5, 1.5, 2.0, 2.5, 3.0], colors=GREY,
                    linewidths=0.8)
    one = ax.contour(G, C, Fz, levels=[1.0], colors=RED, linewidths=1.8)
    D = R["x"][:: max(1, len(R["x"]) // 2000)]
    ax.plot(D[:, ig], D[:, ic], ".", ms=2.0, color="0.2", alpha=0.35,
            label="Draws from the inputs")
    P = R["design"]
    ax.plot(P[:, ig], P[:, ic], "D", ms=7, mfc="white", mec="0.1", mew=1.4,
            label=f"The {R['n_design']} real solves the formula is fitted to")
    ax.set_xlim(g[0], g[-1])
    ax.set_ylim(c[0], c[-1])
    ax.set_xlabel("Unit weight, γ (pcf)")
    ax.set_ylabel("Cohesion, c (psf)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=9.5,
              borderaxespad=0.3, markerscale=1.0)
    gm, cm = _mc_columns(M)
    n = 1000
    Xm = np.zeros((n, len(R["mlv"])))
    Xm[:, ig], Xm[:, ic] = gm[:n], cm[:n]
    real, pred = M["fs"][:n], _surrogate(R, Xm)
    ok = np.isfinite(real)
    lo, hi = 0.2, 3.0
    bx.plot([lo, hi], [lo, hi], color=GREY, lw=1, ls="--", label="Perfect agreement")
    bx.plot(real[ok], pred[ok], "o", ms=3, color=BLUE, alpha=0.6,
            label=f"{ok.sum():,} realizations, solved")
    bx.axvline(1, color=RED, lw=1)
    bx.axhline(1, color=RED, lw=1)
    bx.set_xlim(lo, hi)
    bx.set_ylim(lo, hi)
    bx.set_aspect("equal")
    bx.set_xlabel("F from the solver")
    bx.set_ylabel("F from the formula")
    bx.grid(True, color="0.9")
    bx.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=9.5,
              borderaxespad=0.3)
    fig.tight_layout(w_pad=3.0)

    # each contour's label on its own contour, inside the frame: as near as it can
    # be to where the contour crosses gamma = 97.5 pcf (0.5 and F = 1, lower down,
    # further along) with its ink GAP clear of every draw, solve, other contour,
    # other label and the frame
    def on(level, gx):
        k = int(np.argmin(np.abs(g - gx)))
        col = Fz[:, k]
        i = int(np.argmax(col >= level)) if level > col[0] else 0
        return (g[k], c[i])
    fwd, _ = _to_pt(ax)
    (fx0, fy0), (fx1, fy1) = fwd([(g[0], c[0]), (g[-1], c[-1])])
    marks = [(fwd(D[:, [ig, ic]]), 2.0 / 2),                  # the dots, ms 2
             (fwd(P[:, [ig, ic]]), 7 * math.sqrt(2) / 2 + 1.4 / 2)]   # the diamonds' tips
    sets = [(cs, lv, "%.1f" % lv, 8.5, 0.4) for lv in cs.levels] + [(one, 1.0, "F = 1", 9.5, 0.9)]
    def dense(seg, step=0.5):
        """A polyline resampled every `step` points along its length."""
        s_ = np.r_[0, np.cumsum(np.hypot(*np.diff(seg, axis=0).T))]
        t_ = np.linspace(0, s_[-1], max(2, int(s_[-1] / step) + 1))
        return np.column_stack([np.interp(t_, s_, seg[:, 0]), np.interp(t_, s_, seg[:, 1])])
    lines = {}
    for cset, lv, _, _, hw in sets:
        k = list(cset.levels).index(lv)
        lines[lv] = ([dense(fwd(seg)) for seg in
                      cset.get_paths()[k].to_polygons(closed_only=False) if len(seg) > 1], hw)
    placed = []

    def clearance(p, t, a, b, lv):
        n = np.array([-t[1], t[0]])

        def to_rect(q, r):
            d = q - p
            u, v = np.abs(d @ t) - a, np.abs(d @ n) - b
            return np.hypot(np.maximum(u, 0), np.maximum(v, 0)) - r
        out = min(to_rect(q, r).min() for q, r in marks)
        for olv, (segs, hw) in lines.items():
            if olv != lv:
                out = min(out, min(to_rect(sg, hw).min() for sg in segs))
        corners = p + np.array([[sa * a, sb * b] for sa in (-1, 1) for sb in (-1, 1)]) @ np.array([t, n])
        out = min(out, (corners[:, 0] - fx0).min(), (fx1 - corners[:, 0]).min(),
                  (corners[:, 1] - fy0).min(), (fy1 - corners[:, 1]).min())
        for q in placed:                                        # other labels' inks
            out = min(out, to_rect(q, 0).min())
        return out
    manual = {}
    order = [(1.0, 104)] + [(lv, 122 if lv == 0.5 else 97.5) for lv in cs.levels]
    for lv, gx in order:
        cset, _, text, size, _ = next(st for st in sets if st[1] == lv and
                                      (st[0] is one) == (lv == 1.0))
        x0, y0, x1, y1 = _ink(text, size)
        a, b = (x1 - x0) / 2, (y1 - y0) / 2
        want = fwd([on(lv, gx)])[0]
        segs = lines[lv][0]
        home = int(np.argmin([np.hypot(*(sg - want).T).min() for sg in segs]))
        best = None
        for i, sg in enumerate(segs):
            seglen = np.r_[0, np.cumsum(np.hypot(*np.diff(sg, axis=0).T))]
            j0 = int(np.argmin(np.hypot(*(sg - want).T)))
            for j in range(1, len(sg) - 1):
                t = sg[j + 1] - sg[j - 1]
                t = t / np.hypot(*t)
                if clearance(sg[j], t, a, b, lv) >= GAP:
                    # along the contour from the preferred place
                    key = (abs(seglen[j] - seglen[j0]) if i == home
                           else 1e6 + np.hypot(*(sg[j] - want)))
                    if best is None or key < best[0]:
                        best = (key, sg[j], t)
        assert best is not None, f"no clear place for the {text} label"
        _, p, t = best
        n = np.array([-t[1], t[0]])
        placed.append(p + np.array([[u, v] for u in np.linspace(-a, a, 9)
                                    for v in np.linspace(-b, b, 5)]) @ np.array([t, n]))
        manual[lv] = ax.transData.inverted().transform(p * ax.figure.dpi / 72)
    ax.clabel(cs, fmt="%.1f", fontsize=8.5, inline_spacing=3,
              manual=[manual[lv] for lv in cs.levels])
    ax.clabel(one, fmt="F = 1", fontsize=9.5, inline_spacing=4, manual=[manual[1.0]])
    _save(fig, "rs_surface_fit.png")


def mc_convergence_trace():
    """LEM-11's Monte Carlo, drawn in full by Studio's convergence plot with the
    +/-5% stop target."""
    from xslope.plot import plot_reliability_convergence
    res = run_lem11_mc()
    fig = plt.figure(figsize=(9, 5.5))
    plot_reliability_convergence(res, fig=fig, converge_rel=0.05)
    _save(fig, "mc_convergence_trace.png")


def reliability_mc_vp034():
    """VP34's Phase I fill at COV(phi) = 124%: Monte Carlo on the noncircular
    surface, drawn by Studio's histogram plot."""
    from xslope.fileio import load_slope_data
    from xslope.reliability import reliability_mc
    from xslope.plot import plot_reliability_histogram

    def go():
        ok, res = _quiet(reliability_mc, load_slope_data(VP34), method="spencer",
                         circular=False, search=False, n_samples=10000,
                         num_slices=40)
        assert ok, res
        return res
    res = _run("vp034_mc", go)
    fig = plt.figure(figsize=(9, 5.5))
    plot_reliability_histogram(res, fig=fig)
    _save(fig, "reliability_mc_vp034.png")


FIGURES = {f.__name__: f for f in (
    same_f, r_pf, sigma_normal, sigma_rules, beta_on_f, beta_pf, vp035_surfaces,
    taylor_gamma, taylor_shapes, reliability_taylor_vp036,
    mc_inputs, mc_histogram, mc_check_shape, mc_sampling_band,
    mc_realizations_needed, lhs_vs_random, lhs_segments, lhs_bins, lhs_repeat, rs_surface_fit,
    mc_convergence_trace,
    reliability_mc_vp034)}
# drawn by XSLOPE's own plot functions, in their own style
STUDIO_PLOTS = ("reliability_taylor_vp036", "mc_convergence_trace",
                "reliability_mc_vp034")


def main(argv):
    global CACHE
    if "--cache" in argv:
        i = argv.index("--cache")
        CACHE = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    if "--list" in argv:
        print("\n".join(FIGURES))
        return
    names = argv or list(FIGURES)
    unknown = [n for n in names if n not in FIGURES]
    if unknown:
        raise SystemExit(f"unknown figure(s): {', '.join(unknown)}; --list names them")
    for n in names:
        if n in STUDIO_PLOTS:
            FIGURES[n]()
        else:
            with plt.rc_context(STYLE):
                FIGURES[n]()


if __name__ == "__main__":
    main(sys.argv[1:])
