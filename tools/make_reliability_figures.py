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

``--cache DIR`` keeps the runs' results in DIR between invocations (for
restyling without re-solving); the figures are the same either way.

Figures (page: name):

- index.md: same_f, r_pf, sigma_normal, sigma_rules, beta_on_f, beta_pf,
  vp035_surfaces
- taylor.md: taylor_gamma, taylor_shapes, reliability_taylor_vp036
- monte_carlo.md: mc_inputs, mc_histogram, mc_check_shape, mc_sampling_band,
  mc_realizations_needed, lhs_vs_random, lhs_repeat, mc_convergence_trace,
  reliability_mc_vp034, rs_surface_fit

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

CACHE = None           # --cache DIR
_RUNS = {}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def _save(fig, name, tight=True):
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


def _f1(ax, top, label=True, ha="center"):
    """The line F = 1, labeled at its top."""
    ax.plot([1, 1], [0, top], color=RED, lw=1.5, zorder=4)
    if label:
        ax.text(1, top, "F = 1", color=RED, ha=ha, va="bottom")


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
    for v, lab in ((lcv, "LCV = 200"), (hcv, "HCV = 600")):
        ax.plot([v, v], [0, top], color="0.25", lw=1, ls=":")
        ax.text(v, top, lab, ha="center", va="bottom")
    ya = 0.00045
    ax.annotate("", xy=(lcv, ya), xytext=(hcv, ya),
                arrowprops=dict(arrowstyle="<->", color="0.15", lw=1.1,
                                shrinkA=0, shrinkB=0))
    ax.text(mlv, ya + 0.00006, "±2σ: 95%", ha="center", va="bottom")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1,
              handlelength=2.4, borderaxespad=0)
    ax.set_xlim(50, 750)
    ax.set_ylim(0, top * 1.16)
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
    _f1(ax, top, ha="right")
    ax.plot([F, F], [0, top], color="0.15", lw=1.3)
    ax.text(F, top, "  $F_{MLV}$", ha="left", va="bottom")
    ym = pk * 0.55
    ax.annotate("", xy=(1, ym), xytext=(F, ym),
                arrowprops=dict(arrowstyle="<->", color="0.15", lw=1.1,
                                shrinkA=0, shrinkB=0))
    ax.text((1 + F) / 2, ym + 0.04, "$β\\,σ_F$", ha="center", va="bottom",
            fontsize=12)
    yl = pk * 0.22
    ax.plot([F + sd, F + sd], [0, norm_pdf(F + sd, F, sd)], color=GREY, lw=1,
            ls="--")
    ax.annotate("", xy=(F, yl), xytext=(F + sd, yl),
                arrowprops=dict(arrowstyle="<->", color="0.15", lw=1.1,
                                shrinkA=0, shrinkB=0))
    ax.text(F + sd / 2, yl + 0.04, "$σ_F$", ha="center", va="bottom",
            fontsize=12)
    ax.text(0.93, 0.03, "$P_f$", ha="center", va="bottom", fontsize=12)
    ax.set_xlim(0.4, 2.3)
    ax.set_ylim(0, top * 1.15)
    _no_y(ax)
    ax.set_xlabel("Factor of safety, F")
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
        ax.text(v + 0.09, pv * 1.25, txt, va="bottom")
    bb = 2.326
    ax.plot([0, bb, bb], [0.01, 0.01, 1e-5], color=ORANGE, lw=1.3, ls="--")
    ax.plot(bb, 0.01, "o", color=ORANGE, ms=6, zorder=5)
    ax.text(0.08, 0.0125, "$P_f$ = 1% at β = 2.33", color=ORANGE, va="bottom")
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
    ax.plot(x, y, "o", ms=7, color=BLUE, zorder=4)
    for xi, yi, lab in zip(x, y, ("$F^-$", "$F_{MLV}$", "$F^+$")):
        ax.text(xi - 1.0, yi, f"{lab} = {yi:.3f}", ha="right", va="center")
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
    _f1(ax, top, ha="right")
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
    ax.annotate(f"±{h:.1f} points\nat n = {len(F):,}",
                xy=(len(F) * 0.97, p * 100 + h + 0.25),
                xytext=(len(F) * 0.97, p * 100 + h + 3.2), ha="right",
                va="bottom", arrowprops=dict(arrowstyle="-", color="0.3", lw=0.8,
                                             shrinkA=2, shrinkB=0))
    ax.text(350, p * 100 + 6.0, "95% band", color=BLUE, va="center")
    ax.text(len(F) / 2, p * 100 - 1.6, f"$P_f$ = {p * 100:.2f}%", ha="center",
            va="top")
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
        ax.text(pv * 100 * 1.25, nv * 1.15, f"{pv * 100:g}%: {nr:,.0f}",
                va="bottom")
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
    one routine both sampling engines draw with (same seed), counted in 25 bins
    of equal probability: each bin should hold 10."""
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
        ax.set_xlabel(f"Unit weight, in {nb} bins of equal probability")
        ax.grid(True, axis="y", color="0.9")
    axs[0].set_ylabel("Draws in the bin")
    axs[1].text(nb - 0.3, n / nb + 0.4, "10 per bin", color=ORANGE, ha="right",
                va="bottom")
    fig.tight_layout(w_pad=2.0)
    _save(fig, "lhs_vs_random.png")


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
        ax.text(pf.min() - 0.25, y, name, ha="right", va="center")
        print(f"  {name}: P_f {pf.min():.2f}% to {pf.max():.2f}%, sd {pf.std(ddof=1):.3f} points")
    ax.plot([M["pf_empirical"] * 100] * 2, [-0.7, 1.55], color=GREY, lw=1, ls="--",
            zorder=1)
    ax.text(M["pf_empirical"] * 100, 1.58,
            f"{M['pf_empirical'] * 100:.2f}% from 10,000 realizations", color=GREY,
            ha="center", va="bottom")
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
    # each label on its own contour, inside the frame: where the contour crosses
    # gamma = 97.5 pcf (0.5 and F = 1, lower down, further along)

    def on(level, gx):
        k = int(np.argmin(np.abs(g - gx)))
        col = Fz[:, k]
        i = int(np.argmax(col >= level)) if level > col[0] else 0
        return (g[k], c[i])
    ax.clabel(cs, fmt="%.1f", fontsize=8.5, inline_spacing=3,
              manual=[on(v, 97.5) for v in (1.5, 2.0, 2.5, 3.0)] + [on(0.5, 122)])
    one = ax.contour(G, C, Fz, levels=[1.0], colors=RED, linewidths=1.8)
    ax.clabel(one, fmt="F = 1", fontsize=9.5, inline_spacing=4,
              manual=[on(1.0, 104)])
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
    mc_realizations_needed, lhs_vs_random, lhs_repeat, rs_surface_fit,
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
