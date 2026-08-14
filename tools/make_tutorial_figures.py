"""Regenerate the tutorial figures in ``docs/tutorials/images/``.

One producer for every committed tutorial image that a script can make:

* **Worksheet renders** — the filled ``main`` / ``mat`` / ``profile`` / ``circles``
  sheets of the tutorial's own model, through the same ``tools/render_sheet.py``
  renderer that draws the input-template page's sheet images. A tutorial's Excel
  path tells the reader which cells to fill; the render shows the sheet they should
  end up with, in this problem's numbers rather than the template's blanks.
* **Plot checkpoints** — ``plot_inputs`` / ``plot_solution`` /
  ``plot_circular_search_results`` at the states the tutorial asks the reader to
  compare their screen against, following the conventions of
  ``benchmarks/build_lem_figures.py`` (Agg backend, capture whatever the plot
  function would show, one PNG per state).
* **Placeholders** — a visible TODO card for each figure that only a hand capture
  can produce (the Studio main window, the assistant dock). The page references
  these like any other image, so it renders complete and the gap is visible on it
  rather than tracked somewhere else.

The method matters: a figure is drawn with the same method the tutorial's prose
quotes, so the factor of safety printed on the image is the number in the sentence
beside it. LEM-1's arc is the owner's: Spencer's search on the model as built
(1.276), Bishop's lower answer on a circle Spencer cannot solve (1.215, its crest
slices in tension), then the tension crack at 2c/γ = 8 ft — after which Spencer,
Bishop and M-P land on the same circle and the same 1.084.

Studio's dialog captures are a separate producer, because they need Qt:
``tools/capture_tutorial_screenshots.py``.

Run:  PYTHONPATH=. python3 tools/make_tutorial_figures.py            # everything
      PYTHONPATH=. python3 tools/make_tutorial_figures.py sheets     # one group
      PYTHONPATH=. python3 tools/make_tutorial_figures.py lem01      # by name
"""

from __future__ import annotations

import contextlib
import copy
import importlib.util
import io
import math
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402

from xslope.fileio import load_slope_data                            # noqa: E402
from xslope.generators import generate_starting_circles              # noqa: E402
from xslope.search import circular_search, file_search_window        # noqa: E402
from xslope.plot import (plot_inputs, plot_solution,                 # noqa: E402
                         plot_circular_search_results)

OUT_DIR = os.path.join(REPO_ROOT, "docs", "tutorials", "images")

#: The tutorial's model is the sample's model — one file, two pages (the tutorial
#: builds it, ``docs/lem/samples.md`` catalogues it). Nothing is copied.
LEM01 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_simple_embankment.xlsx")
LEM01_SLICES = 40


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def capture(name, fn, *args, **kwargs):
    """Run a plotting function and save whatever it would show to ``name``.

    The plot functions end in ``plt.show()``; under Agg that draws nothing, so
    ``show`` is swapped for a save. Same idiom as ``benchmarks/build_lem_figures.py``.
    """
    path = os.path.join(OUT_DIR, name)
    saved = []

    def _show(*a, **k):
        plt.gcf().savefig(path, dpi=200, bbox_inches="tight")
        saved.append(path)
        plt.close("all")

    orig = plt.show
    plt.show = _show
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fn(*args, **kwargs)
        if not saved:
            plt.gcf().savefig(path, dpi=200, bbox_inches="tight")
            plt.close("all")
    finally:
        plt.show = orig
    print("-> %s" % name)
    return path


def _render_sheet_module():
    """Import ``tools/render_sheet.py`` as a module (it is a script, not a package)."""
    path = os.path.join(REPO_ROOT, "tools", "render_sheet.py")
    spec = importlib.util.spec_from_file_location("render_sheet", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def render(name, src, sheet, rows=None, cols=None, tab_strip=False):
    # Tutorial captures drop the sheet-tab strip by default: with it, a narrow sheet
    # is padded with empty grid columns out to the strip's fixed width, and the
    # owner's review found the padding made the captures hard to read. The exception
    # is a figure whose subject IS the workbook's shape rather than one sheet's
    # cells — Tutorial 0's tour of the template — where the strip is what the reader
    # is being shown.
    mod = _render_sheet_module()
    out = os.path.join(OUT_DIR, name)
    mod.render_sheet(src, sheet, out, rows=rows, cols=cols, tab_strip=tab_strip)
    print("-> %s  (%s!%s)" % (name, os.path.relpath(src, REPO_ROOT), sheet))
    return out


def placeholder(name, title, lines, width=1200):
    """A visible TODO card standing in for a capture only a person can take.

    Drawn rather than left as a broken image reference: the page has to render
    complete for the layout around the missing figure to be reviewable, and a
    reader who reaches an un-captured state should be told so on the page.

    The height follows the text. A card is a gap marker, not a figure, and one
    sized to a screenshot's aspect outweighs the real figures around it on the
    page — so the box closes under the last line it carries.
    """
    from PIL import Image, ImageDraw
    from matplotlib import font_manager

    def font(px, weight="normal"):
        f = font_manager.findfont(font_manager.FontProperties(
            family="DejaVu Sans", weight=weight))
        from PIL import ImageFont
        return ImageFont.truetype(f, px)

    pad, bar, title_px, line_px, lead = 22, 50, 26, 19, 30
    y_title = bar + pad
    y_first = y_title + title_px + pad
    h = y_first + lead * len(lines) + pad - (lead - line_px)

    im = Image.new("RGB", (width, h), (247, 247, 249))
    d = ImageDraw.Draw(im)
    d.rectangle([6, 6, width - 7, h - 7], outline=(196, 148, 40), width=3)
    d.rectangle([6, 6, width - 7, bar], fill=(250, 226, 160),
                outline=(196, 148, 40), width=3)
    d.text((pad + 6, 16), "SCREENSHOT TO BE CAPTURED", font=font(19, "bold"),
           fill=(105, 76, 12))
    d.text((pad + 6, y_title), title, font=font(title_px, "bold"), fill=(40, 40, 46))
    y = y_first
    for line in lines:
        d.text((pad + 6, y), line, font=font(line_px), fill=(70, 70, 80))
        y += lead
    out = os.path.join(OUT_DIR, name)
    im.save(out)
    print("-> %s  (placeholder, %dx%d)" % (name, width, h))
    return out


# --------------------------------------------------------------------------- #
# LEM-1 — Simple Embankment
# --------------------------------------------------------------------------- #
def lem01_sheets():
    """The four worksheets the tutorial's Excel path fills, as the reader leaves them.

    Column windows, where given, are the real table edge: ``mat`` is 42 columns
    wide and only its strength band is filled here, and ``main``'s content stops at
    the Water-loads row — everything below it on that sheet is the hidden source
    data behind the dropdowns, which is not part of what a reader fills in.

    The ``mat`` window ends at Z because that is where the sheet's own row-9 band
    header (*Shear Strength/Stiffness*, merged C9:Z9) ends. A narrower window makes
    the renderer look up a merge anchor outside its own frame, which is a KeyError
    rather than a quiet crop — the same constraint the input-template manifest
    documents for its three ``mat`` views.
    """
    render("lem01_sheet_main.png", LEM01, "main", rows=(1, 24), cols="A:D")
    render("lem01_sheet_mat.png", LEM01, "mat", rows=(9, 13), cols="A:Z")
    # profile spans two neighbouring tables so the capture reads like the real
    # sheet — Profile Line #2 sits empty beside the filled #1, as the reader
    # sees it (the owner's review caught a one-table crop as "not the sheet").
    render("lem01_sheet_profile.png", LEM01, "profile", rows=(1, 12), cols="A:H")
    render("lem01_sheet_circles.png", LEM01, "circles", rows=(1, 6), cols="A:H")


def lem01_plots():
    """The states the tutorial asks the reader to compare their own screen against."""
    sd = load_slope_data(LEM01)

    # Where the three build paths rejoin: the finished model, before anything is
    # run. Geometry AND the starting circle, because that is the state all three
    # paths leave the reader in, and the page shows this figure at the point it
    # claims they now hold the same model.
    capture("lem01_inputs_geometry.png", plot_inputs, sd,
            title="Slope Geometry and Inputs")

    # The tutorial's run/exploration arc, per the owner's review: no single-circle
    # detour — Spencer's search on the model as built, the anomaly it exposes
    # (Bishop finds a lower circle Spencer cannot solve, its crest slices in
    # tension), then the tension crack added and Spencer re-run, at which point
    # every method lands on the same circle and the same number.
    def search(model, method):
        with contextlib.redirect_stdout(io.StringIO()):
            fs_cache, _, path, circles = circular_search(
                model, method, num_slices=LEM01_SLICES, diagnostic=False,
                **file_search_window(model))
        return fs_cache, path, circles

    fs_cache, path, circles = search(sd, "spencer")
    crit_sp = fs_cache[0]
    capture("lem01_search.png", plot_circular_search_results, sd, fs_cache, path,
            circle_cache=circles)
    capture("lem01_solution_search.png", plot_solution, sd, crit_sp["slices"],
            crit_sp["failure_surface"], crit_sp["solver_result"])

    # Bishop's lower answer on the same crackless model — the crest tension zone
    # (red bars) is drawn by the plot itself and is the anomaly the page reads.
    fs_cache_b, _, _ = search(sd, "bishop")
    crit_b = fs_cache_b[0]
    capture("lem01_solution_bishop.png", plot_solution, sd, crit_b["slices"],
            crit_b["failure_surface"], crit_b["solver_result"])

    # The cracked model: tension crack at the theoretical depth 2c/gamma = 8 ft,
    # dry. Spencer solves every trial and the methods stop disagreeing.
    sc = copy.deepcopy(sd)
    sc["tcrack_depth"] = 8.0
    sc["tcrack_water"] = 0.0
    fs_cache_c, path_c, circles_c = search(sc, "spencer")
    crit_c = fs_cache_c[0]
    capture("lem01_solution_cracked.png", plot_solution, sc, crit_c["slices"],
            crit_c["failure_surface"], crit_c["solver_result"])

    print("   spencer %.4f · bishop %.4f · cracked spencer %.4f (%d slices)"
          % (crit_sp["FS"], crit_b["FS"], crit_c["FS"], LEM01_SLICES))


def lem01_placeholders():
    """The two LEM-1 figures no script can take.

    A full main-window capture is the owner's — it carries the real window chrome
    at the project's fixed capture size — The assistant figure (lem01_assistant.png) is a
    live-session hand capture — the owner's, 2026-08-11, taken after the skill
    hardening round — since a provider conversation cannot be scripted.
    """
    # lem01_studio_canvas.png is now produced by capture_tutorial_screenshots.py
    # (offscreen main-window grab — the hand-capture convention is retired).


# --------------------------------------------------------------------------- #
# Tutorial 0 — Building Models Three Ways
# --------------------------------------------------------------------------- #
#: The blank master template the docs hand out — the file Tutorial 0's Excel path
#: starts from. Rendered blank on purpose: the reader has not filled anything in
#: yet, and the figure's subject is the workbook's shape, not a model.
TEMPLATE = os.path.join(REPO_ROOT, "docs/inputs/input_template.xlsx")


def t0_template():
    """The template's ``main`` worksheet as it is downloaded, out to column G.

    The window is wider than LEM-1's ``A:D`` because the sheet carries its own
    **Sheet / Description** table in ``F7:G21`` — the workbook naming and
    describing every worksheet in itself, which is what the tutorial's paragraph
    beside this figure is about. The tab strip would say the same thing less well:
    it names the sheets without describing them, and it pads a narrow capture out
    to its own fixed width with empty grid columns, which the owner's review of the
    LEM-1 captures rejected.

    Rows stop at 24 (the Surface-family row) for the same reason LEM-1's do —
    below it the sheet holds the hidden source data behind the dropdowns, which is
    not part of what a reader fills in.
    """
    render("t0_template_main.png", TEMPLATE, "main", rows=(1, 24), cols="A:G")


# --------------------------------------------------------------------------- #
# LEM-2 — Loads on the Crest
# --------------------------------------------------------------------------- #
#: LEM-2's model is LEM-1's plus one distributed load, built by
#: ``tools/build_lem02_model.py``. The tutorial's other states are edits the page
#: has the reader make to it, so they are built here in memory rather than as
#: files — exactly as LEM-1 builds its cracked model.
LEM02 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_crest_surcharge.xlsx")
LEM02_SLICES = 40

#: The 750 psf stockpile moved off the level crest and onto the 1:1 face, between
#: elevations 5 and 15. Direction is the whole subject there: perpendicular to a
#: 45° face carries a horizontal component equal to the load itself, while the
#: same stockpile's weight acts straight down. On the crest the two words describe
#: the same force, which is why the comparison cannot be made on the built model.
FACE_LOAD = [[{"X": 5.0, "Y": 5.0, "Normal": 750.0},
              {"X": 15.0, "Y": 15.0, "Normal": 750.0}]]

#: The same total force as the crest surcharge — 750 psf × 10 ft = 7500 lb/ft —
#: gathered onto the single point at the middle of the strip.
LINE_LOAD = [{"x": 30.0, "y": 20.0, "P": 7500.0, "angle": -90.0,
              "label": "footing"}]

LEM02_K = 0.15                     #: the seismic coefficient the page runs
LEM02_TARGET_FS = 1.5              #: the design target the page solves for


@contextlib.contextmanager
def _hold_show():
    """Let a plot function return without saving, so the caller can draw on it.

    ``capture`` saves at ``plt.show()``; a figure that gets annotations after the
    plot function has drawn it has to suppress that first show and let capture's
    fallback save the finished figure instead.
    """
    orig = plt.show
    plt.show = lambda *a, **k: None
    try:
        yield
    finally:
        plt.show = orig


def _lem02_problem_sketch(model):
    """The page's opening sketch: the section, the surcharge, and the numbers.

    ``frame="content"`` crops the panel to the model rather than padding the
    section out to the figure's aspect. The strength and the load are read off
    the model rather than written here, so the sketch cannot print a property
    the file does not carry.
    """
    with _hold_show():
        plot_inputs(model, title="Slope Geometry and Inputs", frame="content")
    ax = plt.gcf().axes[0]
    mat = model["materials"][0]
    ax.text(30.0, 9.0,
            "γ = %g pcf\nc = %g psf\nφ = %g" % (mat["gamma"], mat["c"],
                                                          mat["phi"]),
            ha="center", va="center", fontsize=11)
    band = model["dloads"][0]
    x0, x1 = band[0]["X"], band[-1]["X"]
    q = max(pt["Normal"] for pt in band)
    crest = max(pt["Y"] for pt in band)
    top = crest + q / model["gamma_water"]        # how tall plot_dloads draws it
    # The strip's width is dimensioned BELOW the loaded surface, where nothing
    # else is drawn: above it the arrows run to the top of the frame.
    ax.annotate("", xy=(x0, crest - 1.6), xytext=(x1, crest - 1.6),
                arrowprops=dict(arrowstyle="<->", color="0.25", linewidth=1.0))
    ax.text(0.5 * (x0 + x1), crest - 2.2, "%g ft" % (x1 - x0),
            ha="center", va="top", fontsize=10, color="0.25")
    ax.text(x1 + 1.5, 0.5 * (crest + top), "%g psf" % q,
            ha="left", va="center", fontsize=11)


def _lem02_search(model, method="spencer"):
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path, circles = circular_search(
            model, method, num_slices=LEM02_SLICES, diagnostic=False,
            **file_search_window(model))
    return fs_cache, path, circles


def lem02_sheets():
    """The one worksheet LEM-2's Excel path fills.

    Only ``dloads``: LEM-2 starts from the file LEM-1 finished, so ``main``,
    ``mat``, ``profile`` and ``circles`` are already right and already
    photographed on that page.

    The window runs to column H — the empty gap column and the whole of the
    still-blank Distributed Load #2 beside the filled #1, which is the sheet's
    four-columns-per-load shape a reader adding a second load has to see. It is
    also the narrowest window the renderer will take: #2's header is merged
    F4:H4, and a frame ending at F sends the renderer looking for a merge anchor
    outside it, which is a KeyError rather than a quiet crop (the same constraint
    LEM-1's ``mat`` window documents).

    Rows start at 4, the block header. Above it the sheet carries one merged
    formula cell whose text appears only when the workbook's water loads are
    automatic, and a freshly written file has no cached result for it — so rows
    1-3 render as three empty bands rather than as the note they hold in Excel.
    """
    render("lem02_sheet_dloads.png", LEM02, "dloads", rows=(4, 10), cols="A:H")


def lem02_plots():
    """The states the page compares against, each run with the method it quotes.

    The arc: the surcharge on the model LEM-1 left (Spencer, and every other
    method, on one circle), the same total force concentrated into a line load,
    the same load on the face read both ways round, the seismic coefficient, and
    the design sweep that says what cohesion would hold the target.
    """
    sd = load_slope_data(LEM02)

    # The problem, at the top of the page: the section and the load on it, with the
    # starting circle dropped. A trial circle is a step toward the answer rather
    # than part of the question, and the figure a reader decides the page by should
    # carry only what they are being asked about.
    problem = copy.deepcopy(sd)
    problem["circles"], problem["circular"] = [], False
    capture("lem02_problem.png", _lem02_problem_sketch, problem)

    # Where the three paths rejoin: the same model with the circle back.
    capture("lem02_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    fs_cache, _, _ = _lem02_search(sd)
    crit = fs_cache[0]
    capture("lem02_solution_load.png", plot_solution, sd, crit["slices"],
            crit["failure_surface"], crit["solver_result"])

    # The same resultant as a line load: one point instead of a 10 ft strip.
    ll = copy.deepcopy(sd)
    ll["dloads"], ll["dload_dirs"] = [], []
    ll["line_loads"] = copy.deepcopy(LINE_LOAD)
    fs_ll, _, _ = _lem02_search(ll)
    crit_ll = fs_ll[0]
    capture("lem02_solution_lload.png", plot_solution, ll, crit_ll["slices"],
            crit_ll["failure_surface"], crit_ll["solver_result"])

    # Direction, where it bites: the same load on the 1:1 face, read both ways.
    face = {}
    for word in ("normal", "vertical"):
        fd = copy.deepcopy(sd)
        fd["dloads"] = copy.deepcopy(FACE_LOAD)
        fd["dload_dirs"] = [word]
        fs_f, _, _ = _lem02_search(fd)
        face[word] = fs_f[0]
        capture("lem02_face_%s.png" % word, plot_solution, fd, face[word]["slices"],
                face[word]["failure_surface"], face[word]["solver_result"])

    # A second kind of demand on the same model.
    sk = copy.deepcopy(sd)
    sk["k_seismic"] = LEM02_K
    fs_k, _, _ = _lem02_search(sk)
    crit_k = fs_k[0]
    capture("lem02_solution_seismic.png", plot_solution, sk, crit_k["slices"],
            crit_k["failure_surface"], crit_k["solver_result"])

    # The design sweep: the cohesion that would hold the target under the load.
    # Drawn from the same run the page quotes, so the crossing on the figure is
    # the number in the sentence beside it — and drawn the way Studio's Design
    # view draws it (SweepCanvas.render_design): plot_sensitivity for the curve,
    # then annotate_design_crossing for the answer. A design sweep whose figure
    # showed only the curve would leave its own result off the picture.
    from xslope.sensitivity import design
    from xslope.plot import plot_sensitivity, annotate_design_crossing
    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = design(copy.deepcopy(sd), param="mat:soil:c", low=500.0,
                         high=1200.0, steps=8, target_fs=LEM02_TARGET_FS,
                         method="spencer", search=True, num_slices=LEM02_SLICES)
    if not ok:
        raise SystemExit("LEM-2 design sweep failed: %s" % (res,))

    def _design_figure(summary):
        fig = plot_sensitivity(summary["df"], target_fs=LEM02_TARGET_FS)
        annotate_design_crossing(fig.axes[0], LEM02_TARGET_FS, summary)

    capture("lem02_design.png", _design_figure, res)

    print("   surcharge %.4f · line load %.4f · face normal %.4f · face vertical "
          "%.4f · k=%.2f %.4f · c(FS %.1f) = %.1f"
          % (crit["FS"], crit_ll["FS"], face["normal"]["FS"], face["vertical"]["FS"],
             LEM02_K, crit_k["FS"], LEM02_TARGET_FS, res["crossing"]))


# --------------------------------------------------------------------------- #
# LEM-3 — A Layered Slope
# --------------------------------------------------------------------------- #
#: LEM-3's model is the layered sample — one file, two pages, as LEM-1's is.
LEM03 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_simple_mult_layers.xlsx")
LEM03_SLICES = 40

#: The foundation's cohesion in the page's what-if, against the 800 psf the file
#: carries. Below the fill's 400 psf, so the weaker layer is the deep one and the
#: mechanism the second starting circle describes is the one that controls.
LEM03_WEAK_C = 300.0


def _lem03_search(model, method="spencer"):
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path, circles = circular_search(
            model, method, num_slices=LEM03_SLICES, diagnostic=False,
            **file_search_window(model))
    return fs_cache, path, circles


def lem03_sheets():
    """The three worksheets LEM-3's Excel path fills.

    ``profile`` and ``circles`` take LEM-1's windows on the same sheets — through
    the empty Line #3 beside the two filled lines, and through the empty rows under
    the two circles — so a reader comparing the two pages sees the same frame
    twice.

    ``mat`` does not, and the reason is legibility at the page's width. LEM-1's
    window starts at row 9, which forces the frame out to column Z: row 9 carries
    the sheet's *Shear Strength/Stiffness* band header merged C9:Z9, and a frame
    that ends inside a merge is a KeyError in the renderer rather than a quiet
    crop. Starting at row 10 — the column-name row — drops that merge and lets the
    frame stop at O, the ``u`` column, which is the last one this problem fills.
    The result is 15 named columns across the figure's width instead of 26, and
    the two material rows are readable at the size the page renders them.
    """
    render("lem03_sheet_mat.png", LEM03, "mat", rows=(10, 13), cols="A:O")
    render("lem03_sheet_profile.png", LEM03, "profile", rows=(1, 12), cols="A:H")
    render("lem03_sheet_circles.png", LEM03, "circles", rows=(1, 6), cols="A:H")


def lem03_plots():
    """The states LEM-3 compares against, each from the run whose numbers it quotes.

    The arc is the layer contact: the model as delivered (Spencer's search settling
    on a circle tangent to the top of the foundation), and the same model with the
    foundation softened below the fill, where the minimum moves to the base of the
    foundation and the mechanism the second starting circle describes is the one
    that controls.
    """
    sd = load_slope_data(LEM03)

    capture("lem03_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    fs_cache, path, circles = _lem03_search(sd)
    crit = fs_cache[0]
    capture("lem03_search.png", plot_circular_search_results, sd, fs_cache, path,
            circle_cache=circles)
    capture("lem03_solution.png", plot_solution, sd, crit["slices"],
            crit["failure_surface"], crit["solver_result"])

    # The what-if: the foundation weaker than the fill it carries. Everything else
    # is the delivered model, so the only thing that moved is which layer is weak.
    weak = copy.deepcopy(sd)
    weak["materials"][1]["c"] = LEM03_WEAK_C
    fs_weak, _, _ = _lem03_search(weak)
    crit_w = fs_weak[0]
    capture("lem03_solution_weak.png", plot_solution, weak, crit_w["slices"],
            crit_w["failure_surface"], crit_w["solver_result"])

    print("   as delivered %.4f (tangent depth %.3f) · foundation at c = %g psf "
          "%.4f (tangent depth %.3f)"
          % (crit["FS"], crit["Depth"], LEM03_WEAK_C, crit_w["FS"], crit_w["Depth"]))


# --------------------------------------------------------------------------- #
# LEM-5 — A Weak Layer, Non-Circular
# --------------------------------------------------------------------------- #
#: LEM-5's model is the non-circular sample — one file, two pages, as LEM-1's and
#: LEM-3's are.
LEM05 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_noncircular.xlsx")
LEM05_SLICES = 40


def _lem05_solve(model, non_circ, method="spencer"):
    """One surface, one method, no search — the page's held-surface comparisons."""
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = generate_slices(model, non_circ=non_circ,
                                  num_slices=LEM05_SLICES)
        if not ok:
            raise SystemExit("LEM-5: slicing failed — %s" % (res,))
        slice_df, surface = res
        result = solve_selected(method, slice_df)
    if not isinstance(result, dict):
        raise SystemExit("LEM-5: %s failed — %s" % (method, result))
    return slice_df, surface, result


def lem05_sheets():
    """The four worksheets LEM-5's Excel path fills.

    ``non-circ`` is the one this page is about, and its window runs to column F
    rather than stopping at the three columns the reader types: E2:F5 is the
    sheet's own legend for the Movement column — *Free / Horiz / Fixed* against
    what each one does — which is the vocabulary the step beside the figure
    teaches.

    ``mat`` takes LEM-3's frame (row 10 down, through the ``u`` column) with one
    column more, because this model's fourth material fills ``ru``'s neighbour
    and a window that stopped at O would cut the pore-pressure pair in half.
    ``profile`` runs to L: four lines at three columns each, the fourth ending
    there.
    """
    render("lem05_sheet_mat.png", LEM05, "mat", rows=(10, 15), cols="A:P")
    render("lem05_sheet_profile.png", LEM05, "profile", rows=(1, 12), cols="A:L")
    render("lem05_sheet_piezo.png", LEM05, "piezo", rows=(1, 7), cols="A:E")
    render("lem05_sheet_noncirc.png", LEM05, "non-circ", rows=(1, 8), cols="A:F")


def lem05_plots():
    """LEM-5's arc: the search first, then the comparisons it frames.

    The page runs the automated search before anything else, so the first two
    figures are that search — the trial surfaces it walked through and the
    critical one it settled on, Spencer at 40 slices, which is the method and the
    slice count the prose quotes. Everything after them is a comparison against
    it: the circular search on the same section, and the weak-zone generator's own
    proposal solved as entered, which is the surface the page uses to say a
    generated shape is viable rather than critical.

    The generator's figure is deliberately a SINGLE-SURFACE solve. Its number is
    the one the page reads against the searched minimum, and a search from that
    seed would answer about a third surface again.
    """
    sd = load_slope_data(LEM05)

    capture("lem05_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    # The first analysis the page runs: Auto search, Spencer, from the four points
    # the model carries. The search plot is where the Movement column shows its
    # work — the entered surface among the trials, the critical one in red.
    from xslope.plot import plot_noncircular_search_results
    from xslope.search import noncircular_search
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path = noncircular_search(sd, "spencer",
                                               num_slices=LEM05_SLICES,
                                               diagnostic=False)
    if not fs_cache:
        raise SystemExit("LEM-5: the non-circular search found no valid surface")
    crit = fs_cache[0]
    capture("lem05_search.png", plot_noncircular_search_results, sd, fs_cache, path)
    capture("lem05_solution.png", plot_solution, sd, crit["slices"],
            crit["failure_surface"], crit["solver_result"])

    # The generator's proposal, solved as entered — the audit the Studio path asks
    # the reader to make, and the page's evidence that a generated surface is a
    # starting shape rather than an answer.
    from xslope.generators import generate_noncircular_surface
    gen = generate_noncircular_surface(sd, report=True)
    if not gen["surface"]:
        raise SystemExit("LEM-5: the weak-zone generator built nothing — %s"
                         % (gen["reason"],))
    gs, gsurf, gres = _lem05_solve(sd, gen["surface"])
    capture("lem05_solution_generated.png", plot_solution, sd, gs, gsurf, gres)

    # What a circle gets on the same section — search against search, which is the
    # only fair pairing. The file defines no circles (it is a non-circular model),
    # so the question is asked with the starting circles the geometry itself
    # proposes, refined by the ordinary circular search.
    from xslope.generators import generate_starting_circles
    circ = copy.deepcopy(sd)
    circ["circles"] = generate_starting_circles(circ)
    circ["circular"] = True
    with contextlib.redirect_stdout(io.StringIO()):
        fs_circ, _, _, _ = circular_search(circ, "spencer",
                                           num_slices=LEM05_SLICES,
                                           diagnostic=False)
    best = fs_circ[0]
    capture("lem05_solution_circle.png", plot_solution, circ, best["slices"],
            best["failure_surface"], best["solver_result"])

    print("   searched %.4f · generated as entered %.4f · best circle %.4f "
          "(depth %.3f)"
          % (crit["FS"], gres["FS"], best["FS"], best["Depth"]))


# --------------------------------------------------------------------------- #
# LEM-4 — Water in the Slope
# --------------------------------------------------------------------------- #
#: LEM-4's model is the piezometric-line sample — one file, two pages, as every
#: tutorial's is. The page never writes a variant file: its dry and pinned-circle
#: states are edits the reader makes and undoes, so each is rebuilt here in
#: memory.
LEM04 = os.path.join(REPO_ROOT,
                     "docs/lem/files/xslope_method_slices_problem.xlsx")
LEM04_SLICES = 40

#: The circle the search finds, to the two decimals the page prints and the
#: reader types back into the model. Every figure of the water comparison is
#: this surface under a different water assumption, so the page's dry and wet
#: readings differ by the water and by nothing else.
LEM04_CIRCLE = dict(Xo=182.37, Yo=88.32, Depth=26.90, R=88.32 - 26.90)


def _lem04_solve(model, method="spencer", circle=None):
    """One circle, one method, no search: what this page runs after the search.

    The circle defaults to the one the search finds, which is the surface the
    page holds while the water assumptions change around it.
    """
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = generate_slices(model, circle=circle or LEM04_CIRCLE,
                                  num_slices=LEM04_SLICES)
        if not ok:
            raise SystemExit("LEM-4: slicing failed — %s" % (res,))
        slice_df, surface = res
        result = solve_selected(method, slice_df)
    if not isinstance(result, dict):
        raise SystemExit("LEM-4: %s failed — %s" % (method, result))
    return slice_df, surface, result


def lem04_sheets():
    """The four worksheets LEM-4's Excel path fills.

    ``piezo`` is the one this page is about, and its window runs to G rather than
    stopping at the two columns the reader types: G4:G6 is the sheet's own legend
    for the Type cell — *piezo / phreatic* — which is the choice the step beside
    the figure explains. The rows run to 12 because the line takes eight points.

    ``mat`` takes LEM-5's frame (row 10 down, through column P) for the same
    reason LEM-5 chose it: the u / ru pair in O:P is half the subject, and this
    page also needs both unit-weight columns — C (**g**) and D (**gsat**), both
    filled here — because the weight split is one of the things the page reads.
    ``profile`` runs to I: three lines at three columns each, the third ending
    at H with its spacer.
    """
    render("lem04_sheet_mat.png", LEM04, "mat", rows=(10, 13), cols="A:P")
    render("lem04_sheet_profile.png", LEM04, "profile", rows=(1, 12), cols="A:I")
    render("lem04_sheet_piezo.png", LEM04, "piezo", rows=(1, 12), cols="A:G")
    render("lem04_sheet_circles.png", LEM04, "circles", rows=(1, 6), cols="A:H")


def lem04_plots():
    """The states LEM-4 compares against, in the order the page walks them.

    The search runs first, from the file's own seed, and what it finds is the
    surface the rest of the page is about: every remaining figure is that one
    circle solved by Spencer at 40 slices — wet (as the file ships) and dry
    (u = none in all three materials).
    """
    sd = load_slope_data(LEM04)

    capture("lem04_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    # What Auto search finds on this model: run from the file's own seed circle,
    # it settles on a deep circle through the foundation clay.
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path, circles = circular_search(sd, "spencer",
                                                     num_slices=LEM04_SLICES)
    crit = fs_cache[0]
    capture("lem04_search.png", plot_circular_search_results, sd, fs_cache, path,
            circle_cache=circles)

    # The found circle wet — the model exactly as the build leaves it, on the
    # surface the page has the reader enter after reading it off the search.
    w_slices, w_surface, w_result = _lem04_solve(sd)
    capture("lem04_solution_wet.png", plot_solution, sd, w_slices, w_surface,
            w_result)

    # The same circle dry: every material's pore-pressure option set to none —
    # the reader's three-cell edit, made in memory here. Nothing else moves, so
    # the two figures differ by the pore pressure and by nothing else.
    dry = copy.deepcopy(sd)
    for m in dry["materials"]:
        m["u"] = "none"
    d_slices, d_surface, d_result = _lem04_solve(dry)
    capture("lem04_solution_dry.png", plot_solution, dry, d_slices, d_surface,
            d_result)

    # What the γ_sat column is worth on this circle: the same surface with the
    # saturated weights withheld, which is the comparison the weight-split step
    # reports.
    blank = copy.deepcopy(sd)
    for m in blank["materials"]:
        m["gamma_sat"] = None
    b_slices, _, b_result = _lem04_solve(blank)
    print("   search finds %.4f (Xo %.2f Yo %.2f depth %.2f) · that circle wet "
          "%.4f (ΣW %.0f) · dry %.4f · γ_sat withheld %.4f (ΣW %.0f)"
          % (crit["FS"], crit["Xo"], crit["Yo"], crit["Depth"], w_result["FS"],
             w_slices["w"].sum(), d_result["FS"], b_result["FS"],
             b_slices["w"].sum()))


# --------------------------------------------------------------------------- #
# LEM-6 — Polygon Geometry
# --------------------------------------------------------------------------- #
#: LEM-6's model is the polygon sample — one file, two pages, as every tutorial's
#: is. Its two states beyond the delivered file (a circle pushed below the base,
#: and a softened foundation) are edits the page has the reader make, so both are
#: built in memory here rather than as files.
LEM06 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_sloping_bottom.xlsx")
LEM06_SLICES = 40

#: A circle 1.2 ft deeper than the deepest one that fits inside the domain — the
#: file's second starting circle is tangent to the dipping base at Depth
#: −10.7887, so this one crosses it. It is the page's demonstration of both
#: halves of the composite option: refused without it, truncated with it.
LEM06_DEEP = dict(Xo=20.0, Yo=40.0, Depth=-12.0, R=52.0)

#: The foundation's cohesion in the page's what-if, against the 800 psf the file
#: carries. Below the fill's 400 psf, which is what puts the critical mechanism
#: on the dipping base instead of on the contact above it.
LEM06_WEAK_C = 300.0


def _lem06_search(model, method="spencer", composite=False):
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path, circles = circular_search(
            model, method, num_slices=LEM06_SLICES, diagnostic=False,
            composite=composite, **file_search_window(model))
    return fs_cache, path, circles


def _lem06_solve(model, circle, composite=False, method="spencer"):
    """One circle, no search — the page's composite demonstration.

    Returns ``None`` when the surface cannot be built, because that refusal is
    itself one of the two states the page shows.
    """
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = generate_slices(model, circle=circle, num_slices=LEM06_SLICES,
                                  composite=composite)
        if not ok:
            return None, res
        slice_df, surface = res
        result = solve_selected(method, slice_df)
    if not isinstance(result, dict):
        raise SystemExit("LEM-6: %s failed — %s" % (method, result))
    return (slice_df, surface, result), None


def lem06_sheets():
    """The three worksheets LEM-6's Excel path fills.

    ``polygon`` is the one this page is about. Its window runs to H — three
    polygon blocks at three columns each — so the two filled zones are shown
    beside an empty third, the frame LEM-3 and LEM-4 use on the profile sheet.
    Rows run to 14 so the last vertex row sits above a blank one rather than on
    the frame's edge, which is what tells a reader the block is not full.

    The window starts at row 3 rather than row 1: A2 carries the sheet's
    instruction line, which is text long enough to force the frame out to column
    AB, and 28 columns of empty grid leave the eight the reader fills unreadable
    at the page's width. What A2 says is short enough to quote in the prose.

    ``mat`` and ``circles`` take LEM-3's frames; this section carries the same
    two soils, and its circles sheet fills the same four columns.
    """
    render("lem06_sheet_mat.png", LEM06, "mat", rows=(10, 13), cols="A:O")
    render("lem06_sheet_polygon.png", LEM06, "polygon", rows=(3, 14), cols="A:H")
    render("lem06_sheet_circles.png", LEM06, "circles", rows=(1, 6), cols="A:H")


def _lem06_problem_sketch(model):
    """The page's opening sketch: the two zones, their strengths, and the base.

    The starting circles are dropped from the model the sketch is drawn on — the
    sketch states the problem, and the circles are an input the reader has not
    entered yet. The strengths are read off the file rather than written here, so
    the sketch cannot print a property the model does not carry.
    """
    bare = copy.deepcopy(model)
    bare["circles"], bare["circular"] = [], False
    with _hold_show():
        plot_inputs(bare, title="Slope Geometry and Inputs", frame="content")
    ax = plt.gcf().axes[0]
    embankment, foundation = model["materials"][0], model["materials"][1]
    # Each label sits inside the zone it names, where that zone is thickest: the
    # embankment under its level crest, the foundation at the left edge where the
    # base has dipped away from it and nothing else is drawn.
    ax.text(80.0, 10.0, "%s\nγ = %g pcf\nc = %g psf\nφ = %g"
            % (embankment["name"], embankment["gamma"], embankment["c"],
               embankment["phi"]),
            ha="center", va="center", fontsize=10)
    ax.text(-25.0, -7.5, "%s\nγ = %g pcf\nc = %g psf\nφ = %g"
            % (foundation["name"], foundation["gamma"], foundation["c"],
               foundation["phi"]),
            ha="center", va="center", fontsize=10)
    # The two ends of the dipping base, called out below the hatching that draws it.
    ax.text(-49.0, -17.5, "El. −15", ha="left", va="top", fontsize=10, color="0.25")
    ax.text(119.0, -8.0, "El. −5", ha="right", va="top", fontsize=10, color="0.25")


def lem06_plots():
    """The states LEM-6 reads, in the order the page walks them.

    The search runs first, on the model as delivered, and settles well above the
    dipping base. The two figures after it are what makes the base do something:
    a circle pushed below it and truncated against it (the composite option), and
    the same section with the foundation softened below the fill, where the
    critical mechanism drops onto the base itself.
    """
    sd = load_slope_data(LEM06)

    capture("lem06_problem.png", _lem06_problem_sketch, sd)
    capture("lem06_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    fs_cache, path, circles = _lem06_search(sd)
    crit = fs_cache[0]
    capture("lem06_search.png", plot_circular_search_results, sd, fs_cache, path,
            circle_cache=circles)
    capture("lem06_solution.png", plot_solution, sd, crit["slices"],
            crit["failure_surface"], crit["solver_result"])

    # The circle that will not fit, truncated at the base. The same circle with
    # composite off is the page's refusal, and it produces no figure.
    refused, why = _lem06_solve(sd, LEM06_DEEP, composite=False)
    built, _ = _lem06_solve(sd, LEM06_DEEP, composite=True)
    c_slices, c_surface, c_result = built
    capture("lem06_solution_composite.png", plot_solution, sd, c_slices,
            c_surface, c_result)

    # The what-if: the foundation weaker than the fill it carries, which is what
    # sends the critical surface down to the base the rest of the page describes.
    weak = copy.deepcopy(sd)
    weak["materials"][1]["c"] = LEM06_WEAK_C
    fs_weak, _, _ = _lem06_search(weak)
    crit_w = fs_weak[0]
    capture("lem06_solution_weak.png", plot_solution, weak, crit_w["slices"],
            crit_w["failure_surface"], crit_w["solver_result"])

    print("   as delivered %.4f (Xo %.2f Yo %.2f depth %.3f) · circle at depth %g "
          "refused (%s) / truncated %.4f · foundation at c = %g psf %.4f "
          "(depth %.3f)"
          % (crit["FS"], crit["Xo"], crit["Yo"], crit["Depth"],
             LEM06_DEEP["Depth"], why, c_result["FS"], LEM06_WEAK_C,
             crit_w["FS"], crit_w["Depth"]))


# --------------------------------------------------------------------------- #
# LEM-8 — A Reinforced Slope
# --------------------------------------------------------------------------- #
#: LEM-8's model is the reinforced sample — one file, two pages. The page's two
#: comparisons are edits made in memory: the same section with the reinforcement
#: taken out, and the same lines with a longer pullout length.
LEM08 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_reinforce.xlsx")
LEM08_SLICES = 40


def _lem08_search(model, method="spencer"):
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path, circles = circular_search(
            model, method, num_slices=LEM08_SLICES, diagnostic=False,
            **file_search_window(model))
    return fs_cache, path, circles


def lem08_sheets():
    """The five worksheets LEM-8's Excel path fills.

    ``reinforce`` is the one this page is about. Its window runs from the
    column-name row (row 2) to column R, the last input column — the six lines
    plus two blank rows under them, so the block reads as unfilled below the
    last line. Columns beyond R hold the Type preset's lookup table, which is
    the sheet's own machinery rather than an input.

    ``mat`` and ``profile`` take LEM-3's frames — the same two-material shape,
    and the same window through the empty Line #3 beside the two filled ones —
    ``dloads`` takes LEM-2's four-columns-per-load frame, and ``circles`` the
    frame LEM-3 and LEM-6 use.
    """
    render("lem08_sheet_mat.png", LEM08, "mat", rows=(10, 13), cols="A:O")
    render("lem08_sheet_profile.png", LEM08, "profile", rows=(1, 13), cols="A:H")
    render("lem08_sheet_dloads.png", LEM08, "dloads", rows=(4, 10), cols="A:H")
    # LEM problem: the render stops before the FEM-only columns (Tres, E, Area).
    render("lem08_sheet_reinforce.png", LEM08, "reinforce", rows=(2, 10), cols="A:O")
    render("lem08_sheet_circles.png", LEM08, "circles", rows=(1, 5), cols="A:H")


def lem08_plots():
    """The states LEM-8 reads, in the order the page walks them.

    The search on the model as delivered comes first. The figure beside it is the
    same search run on a copy with the reinforcement lines removed — the page's
    measure of what the six layers are worth, searched rather than borrowed, so
    each answer sits on its own critical circle.
    """
    sd = load_slope_data(LEM08)

    capture("lem08_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    fs_cache, path, circles = _lem08_search(sd)
    crit = fs_cache[0]
    capture("lem08_search.png", plot_circular_search_results, sd, fs_cache, path,
            circle_cache=circles)
    capture("lem08_solution.png", plot_solution, sd, crit["slices"],
            crit["failure_surface"], crit["solver_result"])

    # The comparison: the same section, the same search, no reinforcement.
    bare = copy.deepcopy(sd)
    bare["reinforcement_lines"], bare["reinforce_lines"] = [], []
    fs_bare, _, _ = _lem08_search(bare)
    crit_b = fs_bare[0]
    capture("lem08_solution_bare.png", plot_solution, bare, crit_b["slices"],
            crit_b["failure_surface"], crit_b["solver_result"])

    print("   reinforced %.4f (Xo %.2f Yo %.2f depth %.3f, ΣP %.0f) · "
          "unreinforced %.4f (Xo %.2f Yo %.2f depth %.3f)"
          % (crit["FS"], crit["Xo"], crit["Yo"], crit["Depth"],
             crit["slices"]["p"].sum(), crit_b["FS"], crit_b["Xo"],
             crit_b["Yo"], crit_b["Depth"]))


#: The line lengths LEM-8's length study searches, in feet. The face end of every
#: geogrid is held and the back end moved, so the length IS the variable.
LEM08_LENGTHS = (10, 15, 20, 25, 30, 40)


def _lem08_at_length(model, length):
    """A copy of the model with every geogrid re-cut to ``length``, face end held.

    Every line in this problem is horizontal, so moving the back end to
    ``x1 + length`` changes the length and nothing else. The capacity envelope is
    rebuilt from the new endpoints rather than carried over — the pullout ramps
    are measured from the ends, so a line whose end moved has a different envelope
    even though Tmax, Lp1 and Lp2 are untouched.
    """
    from xslope.fileio import build_reinforce_lines

    out = copy.deepcopy(model)
    for r in out["reinforcement_lines"]:
        r["x2"] = r["x1"] + length
    out["reinforce_lines"] = build_reinforce_lines(out["reinforcement_lines"])
    return out


def lem08_lengths():
    """LEM-8's line-length study: a Spencer search per length, and the figure the
    plateau is read off.

    Each length is its own search — the critical circle is free to move, and below
    20 ft it does. The numbers printed here are the page's table: the factor of
    safety, the tension the found surface actually mobilizes, and how many of the
    six lines it crosses at all, which is what separates the two regimes.

    The figure is the SHORTEST length rather than the longest. At 40 ft the
    critical circle is the one ``lem08_solution.png`` already shows — identical
    center, identical radius — so a shot of it would repeat a figure the page
    carries. At 10 ft the surface is the one the reader cannot otherwise picture:
    it passes behind the back ends, crossing only two of the six lines and both of
    those within a couple of feet of a tip.
    """
    from shapely.geometry import LineString

    from xslope.fileio import reinforce_available_tension

    sd = load_slope_data(LEM08)
    for length in LEM08_LENGTHS:
        model = _lem08_at_length(sd, length)
        crit = _lem08_search(model)[0][0]
        crossed, tension = 0, []
        for r in model["reinforcement_lines"]:
            hit = crit["failure_surface"].intersection(
                LineString([(r["x1"], r["y1"]), (r["x2"], r["y2"])]))
            if hit.is_empty:
                continue
            crossed += 1
            pt = hit if hit.geom_type == "Point" else list(hit.geoms)[0]
            d1 = math.hypot(pt.x - r["x1"], pt.y - r["y1"])
            d2 = math.hypot(pt.x - r["x2"], pt.y - r["y2"])
            tension.append((pt.x, min(d1, d2),
                            reinforce_available_tension(d1, d2, r["t_max"],
                                                        r["lp1"], r["lp2"],
                                                        r.get("tend1", 0.0),
                                                        r.get("tend2", 0.0))))
        print("   %2d ft  FS %.4f  ΣT %5.0f  crossings %d  (Xo %.2f Yo %.2f "
              "depth %.3f)  nearest end %s"
              % (length, crit["FS"], crit["slices"]["p"].sum(), crossed,
                 crit["Xo"], crit["Yo"], crit["Depth"],
                 " ".join("%.1f" % d for _x, d, _t in tension)))
        if length == LEM08_LENGTHS[0]:
            capture("lem08_solution_short.png", plot_solution, model,
                    crit["slices"], crit["failure_surface"],
                    crit["solver_result"])
        if length == LEM08_LENGTHS[-1]:
            capture("lem08_solution_long.png", plot_solution, model,
                    crit["slices"], crit["failure_surface"],
                    crit["solver_result"])


# --------------------------------------------------------------------------- #
# LEM-9 — A Tieback Wall
# --------------------------------------------------------------------------- #
#: LEM-9's model is the tieback-wall verification problem — one file, two pages
#: (the tutorial builds it, ``docs/verification/rocscience.md`` VP49 catalogues
#: it against the SNAILZ manual and Slide). Nothing is copied, and nothing here
#: writes to it: the page's comparisons are edits made on in-memory copies.
LEM09 = os.path.join(REPO_ROOT,
                     "docs/verification/files/rocscience/vp049.xlsx")
#: The search runs at Studio's default slice count; the wedge the manual gives is
#: solved at the 60 the verification row locks, so the page's two numbers are each
#: produced the way their own sentence describes them.
LEM09_SLICES = 40
LEM09_WEDGE_SLICES = 60
LEM09_METHOD = "janbu"


def _lem09_search(model, method=LEM09_METHOD):
    from xslope.search import noncircular_search
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path = noncircular_search(model, method,
                                               num_slices=LEM09_SLICES,
                                               diagnostic=False)
    if not fs_cache:
        raise SystemExit("LEM-9: the non-circular search found no valid surface")
    return fs_cache, path


def _lem09_wedge(model, method=LEM09_METHOD):
    """The manual's wedge, solved as entered — one surface, no search."""
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = generate_slices(model, non_circ=model["non_circ"],
                                  num_slices=LEM09_WEDGE_SLICES)
        if not ok:
            raise SystemExit("LEM-9: slicing failed — %s" % (res,))
        slice_df, surface = res
        result = solve_selected(method, slice_df)
    if not isinstance(result, dict):
        raise SystemExit("LEM-9: %s failed on the wedge — %s" % (method, result))
    return slice_df, surface, result


def _lem09_problem_sketch(model):
    """The page's opening sketch: the wall, the two layers, the two tiebacks.

    The strengths are read off the model, so the sketch cannot print a property
    the file does not carry, and the wall height is measured from the profile
    line's own vertical run rather than written here.
    """
    with _hold_show():
        plot_inputs(model, title="Slope Geometry and Inputs", frame="content")
    ax = plt.gcf().axes[0]
    for mat, (x, y) in zip(model["materials"], ((135.0, 62.0), (95.0, 20.0))):
        ax.text(x, y, "%s\nγ = %g pcf\nc = %g psf\nφ = %g" % (
            mat["name"], mat["gamma"], mat["c"], mat["phi"]),
            ha="center", va="center", fontsize=10)
    face = [pt for pt in model["profile_lines"][1]["coords"] if pt[0] == 0.0]
    top = max(y for _, y in face)
    ax.annotate("", xy=(-6.0, 0.0), xytext=(-6.0, top),
                arrowprops=dict(arrowstyle="<->", color="0.25", linewidth=1.0))
    ax.text(-7.5, 0.5 * top, "%g ft wall" % top, ha="right", va="center",
            rotation=90, fontsize=10, color="0.25")
    for line in model["reinforcement_lines"]:
        ax.text(line["x2"] + 2.0, line["y2"], "%s lb/ft" % f"{line['t_max']:,.0f}",
                ha="left", va="center", fontsize=10)


def lem09_sheets():
    """The five worksheets LEM-9's Excel path fills.

    ``reinforce`` and ``piles`` are the two this page is about, and both windows
    run to the last input column so the reader sees the columns they are told to
    leave alone as well as the ones they fill. ``non-circ`` runs to F for the
    sheet's own Movement legend, as LEM-5's does. ``mat`` and ``profile`` take
    the two-material frames LEM-3 and LEM-8 use.
    """
    render("lem09_sheet_mat.png", LEM09, "mat", rows=(10, 13), cols="A:O")
    render("lem09_sheet_profile.png", LEM09, "profile", rows=(1, 14), cols="A:H")
    render("lem09_sheet_noncirc.png", LEM09, "non-circ", rows=(1, 6), cols="A:F")
    render("lem09_sheet_reinforce.png", LEM09, "reinforce", rows=(2, 6), cols="A:O")
    # LEM problem: stop before the FEM-only pile columns (E, I, Area, Fixity).
    render("lem09_sheet_piles.png", LEM09, "piles", rows=(2, 6), cols="A:M")


def lem09_plots():
    """LEM-9's arc: the search first, then the two things it is read against.

    The automated search on the model as built comes first — Janbu at 40 slices,
    the method the wall's published comparison uses — and its two figures are the
    surfaces the search walked and the critical one it kept. The wedge the
    reference manual gives is a SINGLE-SURFACE solve at the verification row's 60
    slices, announced as such on the page. The last figure is the same search on
    a copy with the tiebacks removed, so the page's measure of what they are
    worth sits on its own critical surface rather than on borrowed geometry.
    """
    from xslope.plot import plot_noncircular_search_results

    sd = load_slope_data(LEM09)

    # The sketch a reader decides the page by carries the question, not a step
    # toward the answer: the entered wedge is dropped from it, as LEM-2 drops its
    # starting circle.
    problem = copy.deepcopy(sd)
    problem["non_circ"] = []
    capture("lem09_problem.png", _lem09_problem_sketch, problem)

    capture("lem09_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    fs_cache, path = _lem09_search(sd)
    crit = fs_cache[0]
    capture("lem09_search.png", plot_noncircular_search_results, sd, fs_cache, path)
    capture("lem09_solution.png", plot_solution, sd, crit["slices"],
            crit["failure_surface"], crit["solver_result"])

    ws, wsurf, wres = _lem09_wedge(sd)
    capture("lem09_solution_wedge.png", plot_solution, sd, ws, wsurf, wres)

    # The comparison: the same wall, the same search, no tiebacks.
    bare = copy.deepcopy(sd)
    bare["reinforcement_lines"], bare["reinforce_lines"] = [], []
    fs_bare, _ = _lem09_search(bare)
    crit_b = fs_bare[0]
    capture("lem09_solution_bare.png", plot_solution, bare, crit_b["slices"],
            crit_b["failure_surface"], crit_b["solver_result"])

    # The circular cross-check: the wedge family is a modeling choice, and the
    # page tests it — one generated seed, the same Janbu, on circles.
    circ = copy.deepcopy(sd)
    with contextlib.redirect_stdout(io.StringIO()):
        circ["circles"] = generate_starting_circles(circ)
        fs_c, _, _, _ = circular_search(circ, LEM09_METHOD,
                                        num_slices=LEM09_SLICES,
                                        diagnostic=False)
    crit_c = fs_c[0]
    capture("lem09_solution_circle.png", plot_solution, circ, crit_c["slices"],
            crit_c["failure_surface"], crit_c["solver_result"])

    print("   searched %.4f (ΣW %.0f, ΣT_x %.0f, pile %.0f) · manual's wedge "
          "%.4f · no tiebacks %.4f"
          % (crit["FS"], crit["slices"]["w"].sum(),
             crit["slices"]["pa_cx"].sum(), crit["slices"]["h_pile"].sum(),
             wres["FS"], crit_b["FS"]))


# --------------------------------------------------------------------------- #
# LEM-10 — Finding the Global Minimum
# --------------------------------------------------------------------------- #
#: LEM-10's model is the multiple-local-minima sample — one file, two pages. The
#: page's comparison is a seeding change, not a model change: the same section
#: searched from the generated circle set and from the deep circle the file
#: carries, so each answer sits on the surface its own seed converged to.
LEM10 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_mult_min_KEY.xlsx")
LEM10_VP75 = os.path.join(REPO_ROOT,
                          "docs/verification/files/rocscience/vp075.xlsx")
LEM10_SLICES = 40


def _lem10_search(model, circles=None, **kwargs):
    """Spencer's search on the model, optionally re-seeded with ``circles``."""
    m = copy.deepcopy(model)
    if circles is not None:
        m["circles"] = circles
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path, circle_cache = circular_search(
            m, "spencer", num_slices=LEM10_SLICES, diagnostic=False,
            **dict(file_search_window(m), **kwargs))
    return fs_cache, path, circle_cache


def lem10_plots():
    """The two answers LEM-10 puts side by side, each on its own search.

    The first search is seeded with the generator's embankment circle — tangent
    to the top of the foundation, the shallower of the two mechanisms this
    section holds. It collapses onto the infinite-slope sliver, so its figure is
    the search plot: the story is the family of trial circles shrinking onto the
    face, not the surface that won. The generator's set as a whole reaches the
    same sliver, but the skim circle in it has a center 500 ft above the section
    and drawing that search puts the model in a corner of the frame.

    The second is seeded with the deep circle the file carries, tangent to the
    limiting depth, and finds the foundation mechanism. That one is drawn as a
    solution, because the mass it moves is the point.
    """
    sd = load_slope_data(LEM10)

    capture("lem10_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    generated = generate_starting_circles(sd, report=True)["circles"]
    shallow = [c for c in generated if abs(c["Depth"]) < 1e-9]
    fs_gen, path_gen, circles_gen = _lem10_search(sd, shallow)
    crit_g = fs_gen[0]
    capture("lem10_search_shallow.png", plot_circular_search_results, sd,
            fs_gen, path_gen, circle_cache=circles_gen)

    fs_deep, _, _ = _lem10_search(sd)
    crit_d = fs_deep[0]
    capture("lem10_solution_deep.png", plot_solution, sd, crit_d["slices"],
            crit_d["failure_surface"], crit_d["solver_result"])

    # The tools section's middle answer: the embankment seed with the 5 ft
    # surficial filter on — the search stays in the fill and returns the 33 ft
    # circle, deeper than the sliver and shallower than the foundation.
    fs_f5, _, _ = _lem10_search(sd, shallow, min_slip_depth=5.0)
    crit_f = fs_f5[0]
    capture("lem10_solution_filter5.png", plot_solution, sd, crit_f["slices"],
            crit_f["failure_surface"], crit_f["solver_result"])

    # The second slope: VP75, the James Bay dyke — the corpus's local-minimum
    # showcase. The generated per-layer seeds settle in the fill; the grid
    # finds the deep circle through all three clays. Both at the Run dialog's
    # default 40 slices, which is the flow the page walks.
    vp75 = load_slope_data(LEM10_VP75)
    capture("lem10_vp75_inputs.png", plot_inputs, vp75,
            title="Slope Geometry and Inputs")
    gen75 = generate_starting_circles(vp75)
    # The trap the page shows is the SINGLE plausible seed: the per-layer circle
    # tangent to the marine-clay top (the one an engineer might place, expecting
    # the failure in the fill), searched alone with the 2 m surficial filter.
    single = [c for c in gen75 if abs(c["Depth"] - 15.0) < 1e-6]
    with contextlib.redirect_stdout(io.StringIO()):
        fs_1, _, _, _ = circular_search(dict(vp75, circles=single), "spencer",
                                        num_slices=40, min_slip_depth=2.0,
                                        diagnostic=False)
        fs_all, _, _, _ = circular_search(dict(vp75, circles=gen75), "spencer",
                                          num_slices=40, min_slip_depth=2.0,
                                          diagnostic=False)
        fs_grid, _, _, _ = circular_search(copy.deepcopy(vp75), "spencer",
                                           num_slices=40, seed="grid",
                                           diagnostic=False)
    crit_1, crit_all, crit_grid = fs_1[0], fs_all[0], fs_grid[0]
    capture("lem10_vp75_single.png", plot_solution, vp75,
            crit_1["slices"], crit_1["failure_surface"],
            crit_1["solver_result"])
    capture("lem10_vp75_grid.png", plot_solution, vp75,
            crit_grid["slices"], crit_grid["failure_surface"],
            crit_grid["solver_result"])
    print("   vp75 single-seed %.4f · full set+filter %.4f · grid %.4f"
          % (crit_1["FS"], crit_all["FS"], crit_grid["FS"]))

    print("   shallow seed %.4f (Xo %.2f Yo %.2f depth %.3f, ΣW %.0f) · "
          "deep seed %.4f (Xo %.2f Yo %.2f depth %.3f, ΣW %.0f) · "
          "filter-5 %.4f (Xo %.2f Yo %.2f depth %.3f)"
          % (crit_g["FS"], crit_g["Xo"], crit_g["Yo"], crit_g["Depth"],
             crit_g["slices"]["w"].sum(), crit_d["FS"], crit_d["Xo"],
             crit_d["Yo"], crit_d["Depth"], crit_d["slices"]["w"].sum(),
             crit_f["FS"], crit_f["Xo"], crit_f["Yo"], crit_f["Depth"]))


GROUPS = {
    "t0_template": t0_template,
    "lem01_sheets": lem01_sheets,
    "lem01_plots": lem01_plots,
    "lem01_placeholders": lem01_placeholders,
    "lem02_sheets": lem02_sheets,
    "lem02_plots": lem02_plots,
    "lem03_sheets": lem03_sheets,
    "lem03_plots": lem03_plots,
    "lem05_sheets": lem05_sheets,
    "lem05_plots": lem05_plots,
    "lem04_sheets": lem04_sheets,
    "lem04_plots": lem04_plots,
    "lem06_sheets": lem06_sheets,
    "lem06_plots": lem06_plots,
    "lem08_sheets": lem08_sheets,
    "lem08_plots": lem08_plots,
    "lem08_lengths": lem08_lengths,
    "lem09_sheets": lem09_sheets,
    "lem09_plots": lem09_plots,
    "lem10_plots": lem10_plots,
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    os.makedirs(OUT_DIR, exist_ok=True)
    names = [n for n in GROUPS
             if not argv or any(a in n for a in argv)]
    if not names:
        print("no figure group matching %s; known groups: %s"
              % (argv, ", ".join(sorted(GROUPS))))
        return 1
    for name in names:
        print("== %s" % name)
        GROUPS[name]()
    print("\nwrote %d group(s) to docs/tutorials/images/" % len(names))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
