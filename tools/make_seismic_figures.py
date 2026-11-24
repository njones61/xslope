"""Render the run figures of the Seismic Analysis pages (docs/seismic/images/).

- ``seismic_ky_curve.png``: VP104 (the Slide2 verification problem the Newmark
  example comes from), F against k on the circle with the minimum yield
  coefficient, Spencer's method, 40 slices, with k_y where F = 1.
- ``seismic_newmark_vp104.png``: VP104's acceleration record integrated at that
  k_y, both directions, in the form of Duncan, Wright, and Brandon (2014),
  Fig. 10.3.
- ``seismic_ky_surfaces.png``: the slope of DWB Fig. 10.7 with soil 2's R
  envelope, dry: the static critical circle and the circle with the minimum k_y.
- ``seismic_envelopes.png``: DWB Table 10.3's four soils, each one's R envelope
  against its K_c = 1 envelope and its drained envelope, read from the
  verification workbook, over the range of stress in the Fig. 10.7 slope.
- ``seismic_single_two_stage.png``: DWB Table 10.4, the single-stage factor of
  safety against the two-stage one, XSLOPE's values (the verification page's
  test tags) beside the book's.

Run from the repo root:
    python tools/make_seismic_figures.py            # every figure
    python tools/make_seismic_figures.py envelopes  # one (ky_curve, newmark,
                                                     # ky_surfaces, envelopes,
                                                     # single_two_stage)
"""
import contextlib
import io
import math
import os
import re
import sys
import warnings

warnings.filterwarnings("ignore")
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from xslope.fileio import load_slope_data  # noqa: E402
from xslope.plot import plot_inputs, plot_newmark  # noqa: E402
from xslope.seismic import (_OneSurface, newmark_displacement,  # noqa: E402
                            yield_coefficient)

OUT = os.path.join(ROOT, "docs", "seismic", "images")
VER = os.path.join(ROOT, "docs", "verification")
VP104 = os.path.join(VER, "files", "rocscience", "vp104a.xlsx")
DWB_DRY = os.path.join(VER, "files", "published", "dwb_t10_4_dry.xlsx")
PUBLISHED = os.path.join(VER, "published.md")
NUM_SLICES = 40
DPI = 200

# Duncan, Wright, and Brandon (2014), Table 10.4: the book's factors of safety,
# (single-stage, two-stage) for soils 1 to 4.
DWB_T10_4 = {
    "dry": [(0.95, 1.06), (1.56, 1.77), (1.07, 1.19), (2.76, 3.42)],
    "submerged": [(0.83, 0.95), (1.59, 1.79), (1.10, 1.21), (2.83, 3.49)],
}

BLUE = "#1f4e9c"
ORANGE = "#d9822b"
GREY = "0.45"


def _quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def _save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(path, ROOT))


def _one_soil(path, name):
    """The workbook with every material row a copy of the row named ``name``, as
    the verification tags' ``material=`` runs it."""
    sd = load_slope_data(path)
    row = [m for m in sd["materials"] if m.get("name") == name]
    assert len(row) == 1, name
    sd["materials"] = [dict(row[0]) for _ in sd["materials"]]
    return sd


def _crest_to(sd, x_end):
    """The model with its crest run out to ``x = x_end`` and nothing else changed:
    the profile line, its polygon, the domain, and the ground surface."""
    from shapely.geometry import LineString, Polygon
    pts = list(sd["profile_lines"][0]["coords"])
    assert len(sd["profile_lines"]) == 1 and pts[-1][1] == pts[-2][1], pts
    pts[-1] = (x_end, pts[-1][1])
    base = sd["max_depth"]
    poly = Polygon(pts + [(x_end, base), (pts[0][0], base)])
    out = dict(sd)
    out["profile_lines"] = [dict(sd["profile_lines"][0], coords=pts)]
    out["polygons"] = [dict(sd["polygons"][0], polygon=poly)]
    out["domain_polygon"] = poly
    out["ground_surface"] = LineString(pts)
    return out


# ---------------------------------------------------------------------------
# VP104: k_y, F against k, and the Newmark histories
# ---------------------------------------------------------------------------

_VP104 = {}


def vp104_yield():
    if not _VP104:
        sd = load_slope_data(VP104)
        bundle = _quiet(yield_coefficient, sd, "spencer", analysis="auto_search",
                        surface="circular", num_slices=NUM_SLICES, announce=False)
        _VP104.update(sd=sd, info=bundle["yield"])
    return _VP104


def ky_curve():
    v = vp104_yield()
    info = v["info"]
    ky = info["k_y"]
    one = _OneSurface(v["sd"], "spencer", circle=info["surface"], non_circ=None,
                      num_slices=NUM_SLICES, composite=False, seismic2=False)
    ks = [i * 0.01 for i in range(0, 26)]
    pts = [(k, one(k)[0]) for k in ks]
    pts = [(k, f) for k, f in pts if f is not None]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot([k for k, _ in pts], [f for _, f in pts], color=BLUE, lw=2)
    ax.axhline(1.0, color=GREY, lw=1, ls="--")
    ax.plot([ky], [1.0], "o", color=ORANGE, ms=8, zorder=5)
    ax.annotate(f"k$_y$ = {ky:.4f}", xy=(ky, 1.0), xytext=(ky + 0.02, 1.12),
                fontsize=11, arrowprops=dict(arrowstyle="-", color=GREY, lw=0.8,
                                             shrinkA=0, shrinkB=6))
    ax.set_xlabel("Horizontal seismic coefficient, k$_h$")
    ax.set_ylabel("Factor of safety, F")
    ax.set_xlim(0, ks[-1])
    ax.grid(True, color="0.9")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    _save(fig, "seismic_ky_curve.png")
    print(f"  VP104 static critical F = {info['static']['FS']:.3f}, k_y on it "
          f"{info['static']['k_y']:.4f}; minimum k_y = {ky:.4f}; F at k = 0 on the "
          f"k_y circle = {info['FS_k0']:.3f}")


def newmark():
    ky = vp104_yield()["info"]["k_y"]
    # VP104's record is inside Rocscience's own tutorial model, which is not in
    # this repository: it is read from the local vendor archive.
    sys.path.insert(0, os.path.join(ROOT, "benchmarks", "rocscience"))
    from newmark_vp104 import MODEL_NAME, MODEL_URL, VENDOR_FILES, load_record
    rec = load_record()
    if rec is None:
        raise SystemExit(f"{MODEL_NAME} is not under {VENDOR_FILES}; it is in "
                         f"{MODEL_URL}.")
    res = newmark_displacement(rec, round(ky, 4), units="si")
    fig = plt.figure(figsize=(9.0, 7.5))
    plot_newmark(res, fig=fig)
    _save(fig, "seismic_newmark_vp104.png")
    d = res["directions"]
    print(f"  k_y {round(ky, 4)}: D = {res['D']:.3f} {res['D_unit']} "
          f"(+: {d[1]['D']:.3f}, {d[1]['episodes']} episodes; -: {d[-1]['D']:.3f}, "
          f"{d[-1]['episodes']} episodes); record {rec['description']!r}, "
          f"{rec['npts']} points, dt {rec['dt']}, PGA {rec['pga']:.3f} g")


def newmark_steps():
    """The sliding-block integration on a record simple enough to follow: three
    cycles of a 1 Hz sine wave with a 0.3 g peak, k_y = 0.15, positive
    accelerations driving the block downslope. The block slides while its
    velocity relative to the ground is above zero; those spans are shaded in all
    three panels."""
    import numpy as np
    ky, peak, cycles = 0.15, 0.30, 3
    time = np.linspace(0.0, float(cycles), 3001)
    rec = {"time": time, "accel": peak * np.sin(2.0 * np.pi * time)}
    res = newmark_displacement(rec, ky, units="si")
    run = res["directions"][+1]
    v, d, a = run["v"], run["d"], res["accel"]
    sliding = v > 1e-9
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 7.0), sharex=True)
    for ax in axes:
        ax.fill_between(time, 0, 1, where=sliding, transform=ax.get_xaxis_transform(),
                        color=ORANGE, alpha=0.12, lw=0)
        ax.grid(True, color="0.93")
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    ax = axes[0]
    ax.plot(time, a, color=BLUE, lw=1.6)
    ax.axhline(ky, color=ORANGE, lw=1.4, ls="--")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.fill_between(time, ky, a, where=a > ky, color=ORANGE, alpha=0.55, lw=0)
    ax.text(time[-1], ky + 0.01, "k$_y$ = 0.15", color=ORANGE, fontsize=9,
            ha="right", va="bottom")
    ax.set_ylabel("Ground acceleration\na(t) (g)")
    ax = axes[1]
    ax.plot(time, v, color=BLUE, lw=1.6)
    ax.set_ylabel("Block velocity relative\nto the ground (cm/s)")
    ax = axes[2]
    ax.plot(time, d, color=BLUE, lw=1.6)
    ax.set_ylabel("Displacement\n(cm)")
    ax.set_xlabel("Time (s)")
    ax.set_xlim(0, time[-1])
    fig.align_ylabels(axes)
    _save(fig, "seismic_newmark_steps.png")
    print(f"  sine record: peak {peak} g, k_y {ky}: D = {run['D']:.2f} cm in "
          f"{run['episodes']} episodes, v_max {run['v_max']:.1f} cm/s")


# ---------------------------------------------------------------------------
# DWB Fig. 10.7: the static critical circle and the minimum-k_y circle
# ---------------------------------------------------------------------------

def ky_surfaces():
    # The verification model's crest ends at x = 165, which stops soil 2's
    # lowest-k_y circle. Run out to x = 230, the crest leaves room past where that
    # circle comes out (x = 217); nothing else in the model changes.
    sd = _crest_to(_one_soil(DWB_DRY, "Soil 2 R envelope"), 230.0)
    bundle = _quiet(yield_coefficient, sd, "spencer", analysis="auto_search",
                    surface="circular", num_slices=NUM_SLICES, announce=False)
    info = bundle["yield"]
    st = info["static"]
    lines = []
    for surf in (st["surface"], info["surface"]):
        one = _OneSurface(sd, "spencer", circle=surf, non_circ=None,
                          num_slices=NUM_SLICES, composite=False, seismic2=False)
        lines.append(one.failure_surface)
    fig = plt.figure(figsize=(10.0, 4.6))
    plot_inputs(sd, mode="shared", show_title=False, show_legend=False, fig=fig,
                frame="content")
    ax = fig.axes[0]
    labels = (f"Static critical circle (F = {st['FS']:.2f}; on it, k$_y$ = {st['k_y']:.3f})",
              f"Circle with the minimum k$_y$ ({info['k_y']:.3f})")
    handles = []
    for line, color, ls, lab in zip(lines, (BLUE, ORANGE), ("-", "--"), labels):
        x, y = zip(*line.coords)
        handles += ax.plot(x, y, color=color, lw=2.4, ls=ls, label=lab)
    # The sky over the toe, left of the slope face, is empty: the legend goes there.
    ax.legend(handles=handles, loc="upper left", fontsize=9, frameon=False)
    _save(fig, "seismic_ky_surfaces.png")
    print(f"  Fig. 10.7 soil 2 R envelope: static F {st['FS']:.3f}, k_y on it "
          f"{st['k_y']:.4f}; minimum k_y {info['k_y']:.4f} "
          f"({(info['k_y'] / st['k_y'] - 1) * 100:+.1f}%); circles {st['surface']} "
          f"and {info['surface']}")


# ---------------------------------------------------------------------------
# DWB Table 10.3 and 10.4
# ---------------------------------------------------------------------------

def envelopes():
    sd = load_slope_data(DWB_DRY)
    mats = {m["name"]: m for m in sd["materials"]}
    # The largest effective stress in the slope: its full height of soil.
    ys = [y for line in sd["profile_lines"] for _, y in line["coords"]]
    s_max = (max(ys) - min(ys)) * max(float(m.get("gamma") or 0) for m in mats.values())
    s_max = math.ceil(s_max / 1000.0) * 1000.0
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0))
    for n, ax in zip(range(1, 5), axes.flat):
        two, r = mats[f"Soil {n} two-stage"], mats[f"Soil {n} R envelope"]
        s = [0.0, s_max]
        d, psi = float(two["d"]), float(two["psi"])
        c, phi = float(two["c"]), float(two["phi"])
        cr, phir = float(r["c"]), float(r["phi"])
        ax.plot(s, [d + x * math.tan(math.radians(psi)) for x in s], color=BLUE, lw=2,
                label=f"K$_c$ = 1: d = {d:g} psf, ψ = {psi:g}°")
        ax.plot(s, [cr + x * math.tan(math.radians(phir)) for x in s], color=ORANGE, lw=2,
                ls="--", label=f"R: c$_R$ = {cr:g} psf, φ$_R$ = {phir:g}°")
        ax.plot(s, [c + x * math.tan(math.radians(phi)) for x in s], color=GREY, lw=1.2,
                ls=":", label=f"Drained: c′ = {c:g} psf, φ′ = {phi:g}°")
        ax.set_title(f"Soil {n}", fontsize=11)
        ax.set_xlim(0, s_max)
        ax.set_ylim(0, None)
        ax.grid(True, color="0.9")
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(fontsize=8.5, frameon=False, loc="upper left")
        ax.set_xlabel("σ′$_{fc}$ or σ (psf)")
        ax.set_ylabel("τ$_{ff}$ or τ (psf)")
    fig.tight_layout()
    _save(fig, "seismic_envelopes.png")


def _tagged_fs():
    """XSLOPE's Table 10.4 values: the expected_fs of the verification page's
    DWB-10.4 test tags, keyed (water, soil, "single" or "kc1")."""
    text = open(PUBLISHED, encoding="utf-8").read()
    out = {}
    for fs, water, soil, kind in re.findall(
            r"expected_fs=([\d.]+)[^>]*benchmark=DWB-10\.4-(dry|submerged)-(\d)-(single|kc1)\b",
            text):
        out[(water, int(soil), kind)] = float(fs)
    assert len(out) == 16, sorted(out)
    return out


def single_two_stage():
    xs = _tagged_fs()
    fig, ax = plt.subplots(figsize=(6.0, 5.4))
    top = 3.8
    # Equal aspect, so a line of slope frac is drawn at atan(frac) on the page.
    for frac, lab in ((1.0, "single = two-stage"), (0.9, "90 percent"), (0.8, "80 percent")):
        ax.plot([0, top], [0, frac * top], color="0.75", lw=1, ls="-" if frac == 1 else "--")
        x = 2.45
        ax.text(x, frac * x + 0.04, lab, fontsize=9, color="0.4", ha="center", va="bottom",
                rotation=math.degrees(math.atan(frac)), rotation_mode="anchor")
    for water, marker in (("dry", "o"), ("submerged", "s")):
        two = [xs[(water, n, "kc1")] for n in range(1, 5)]
        one = [xs[(water, n, "single")] for n in range(1, 5)]
        ax.plot(two, one, marker, color=BLUE, ms=8, ls="none",
                label=f"XSLOPE, {water}")
        book = DWB_T10_4[water]
        ax.plot([b[1] for b in book], [b[0] for b in book], marker, mfc="none",
                mec=ORANGE, mew=1.6, ms=13, ls="none",
                label=f"Duncan et al. (2014), {water}")
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_aspect("equal")
    ax.set_xlabel("Two-stage factor of safety")
    ax.set_ylabel("Single-stage factor of safety")
    ax.grid(True, color="0.93")
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    _save(fig, "seismic_single_two_stage.png")
    for water in ("dry", "submerged"):
        print("  " + water + ": " + ", ".join(
            f"soil {n} {xs[(water, n, 'single')] / xs[(water, n, 'kc1')]:.3f}"
            for n in range(1, 5)))


FIGURES = {"ky_curve": ky_curve, "newmark": newmark, "newmark_steps": newmark_steps,
           "ky_surfaces": ky_surfaces,
           "envelopes": envelopes, "single_two_stage": single_two_stage}

if __name__ == "__main__":
    for name in sys.argv[1:] or list(FIGURES):
        FIGURES[name]()
