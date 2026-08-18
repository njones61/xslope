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


def render(name, src, sheet, rows=None, cols=None, tab_strip=False,
           identity_cols=None):
    # Tutorial captures drop the sheet-tab strip by default: with it, a narrow sheet
    # is padded with empty grid columns out to the strip's fixed width, and the
    # owner's review found the padding made the captures hard to read. The exception
    # is a figure whose subject IS the workbook's shape rather than one sheet's
    # cells — Tutorial 0's tour of the template — where the strip is what the reader
    # is being shown.
    #
    # identity_cols re-shows a leading band at the left of a distant window, the
    # split-view idiom the input-template manifest uses for the mat sheet: a capture
    # of the seepage columns forty across is unreadable without the ID and the name
    # of the row they belong to.
    mod = _render_sheet_module()
    out = os.path.join(OUT_DIR, name)
    mod.render_sheet(src, sheet, out, rows=rows, cols=cols, tab_strip=tab_strip,
                     identity_cols=identity_cols)
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
LEM10_JB = os.path.join(REPO_ROOT, "docs/lem/files/xslope_james_bay.xlsx")
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
    # Part B runs on the tutorial's own copy of VP75: the same dyke with a
    # single mid-depth seed on the circles sheet — the circle an engineer might
    # place expecting the failure under the crest — so the reader opens and
    # runs without editing anything.
    vp75 = load_slope_data(LEM10_JB)
    capture("lem10_vp75_inputs.png", plot_inputs, vp75,
            title="Slope Geometry and Inputs")
    gen75 = generate_starting_circles(vp75)
    with contextlib.redirect_stdout(io.StringIO()):
        fs_1, _, _, _ = circular_search(copy.deepcopy(vp75), "spencer",
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


# --------------------------------------------------------------------------- #
# LEM-7 — Strength Options Beyond Mohr-Coulomb
# --------------------------------------------------------------------------- #
#: Part A is Baker's compacted-clay slope (VP44) — the file carries the power
#: curve, and the page's one edit swaps it for the Mohr-Coulomb envelope fitted to
#: the same triaxial data. Part B is Low's layered undrained slope (VP23), whose
#: lower layer carries a strength that grows with depth; the page's edit replaces
#: that profile with a single constant. Both are tutorial copies of the corpus
#: models, so each figure is the verification answer the reader reproduces.
LEM07_A = os.path.join(REPO_ROOT, "docs/lem/files/xslope_baker_clay.xlsx")
LEM07_B = os.path.join(REPO_ROOT, "docs/lem/files/xslope_low_clay.xlsx")
#: The counterpoint pair — Baker's London clay (VP61), results only, no download.
LEM07_L_POW = os.path.join(REPO_ROOT, "docs/verification/files/rocscience/vp061a.xlsx")
LEM07_L_MC = os.path.join(REPO_ROOT, "docs/verification/files/rocscience/vp061b.xlsx")
#: Slice counts are the corpus tags': VP44 runs Spencer at 40, VP23 runs at 50.
LEM07_A_SLICES = 40
LEM07_B_SLICES = 50
#: Part B's edit: the constant an engineer reaches for when a linear su profile is
#: replaced by one number — the average of the layer's 15 kN/m² top and 30 kN/m²
#: bottom.
LEM07_B_CONST = 22.5


def _lem07_search(model, method, num_slices, **kwargs):
    """The page's search: the model's own circle, refined, nothing else set."""
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, _, _ = circular_search(
            copy.deepcopy(model), method, num_slices=num_slices, diagnostic=False,
            **dict(file_search_window(model), **kwargs))
    return fs_cache[0]


def _lem07_solution(name, model, crit):
    capture(name, plot_solution, model, crit["slices"], crit["failure_surface"],
            crit["solver_result"])


def _lem07_reading(crit, model=None, mat=None):
    """The numbers the page quotes off one search: FS, the circle, the surface.

    ``mat`` is a 1-based material id; when given, the length-weighted mean of the
    strength actually mobilized in that material — the honest single number for a
    layer whose strength varies along the surface — is reported with the mean base
    elevation it was mobilized at.
    """
    s = crit["slices"]
    out = ("FS %.6f · Xo %.4f Yo %.4f · tangent elev %.4f · L %.2f · ΣW %.1f"
           % (crit["FS"], crit["Xo"], crit["Yo"], crit["Depth"],
              s["dl"].sum(), s["w"].sum()))
    if mat is not None:
        m = s[s["mat"] == mat]
        if len(m):
            out += ("  |  mat %d: %.2f of the surface, mean strength %.4f at mean "
                    "base elev %.4f" % (mat, m["dl"].sum(),
                                        (m["c"] * m["dl"]).sum() / m["dl"].sum(),
                                        (m["y_cb"] * m["dl"]).sum() / m["dl"].sum()))
    return out


def lem07_plots():
    """Both parts' searches, each figure drawn with the method the prose quotes.

    Part A's two answers come off the same file: as downloaded (the power curve)
    and after the page's one edit (the fitted Mohr-Coulomb envelope). The London
    clay pair is the counterpoint — the same 43° slope with a data set that
    includes low-stress measurements, where the two envelopes nearly agree — and
    is shown as results only, so its models are read from the corpus rather than
    copied into the tutorial's file set.

    Part B's two answers likewise come off one file: the lower layer's strength
    profile as shipped, and that profile flattened to a constant. The interesting
    quantity there is not only the factor of safety but the tangent elevation, so
    the diagnostics print it for every run.
    """
    # ---- Part A: Baker example 1, compacted Israeli clays ------------------ #
    a = load_slope_data(LEM07_A)
    capture("lem07_baker_inputs.png", plot_inputs, a,
            title="Slope Geometry and Inputs")

    crit_pow = _lem07_search(a, "spencer", LEM07_A_SLICES)
    _lem07_solution("lem07_baker_pow.png", a, crit_pow)

    # The edit the page teaches: option pow -> mc, with Baker's fitted envelope.
    # The power-curve coefficients are left on the row exactly as a reader would
    # leave them — inert under mc, and the answer is identical either way.
    mc = dict(a["materials"][0])
    mc.update({"option": "mc", "c": 11.64, "phi": 24.7})
    a_mc = dict(a, materials=[mc])
    crit_mc = _lem07_search(a_mc, "spencer", LEM07_A_SLICES)
    _lem07_solution("lem07_baker_mc.png", a_mc, crit_mc)

    # The counterpoint: London clay, both envelopes, results only.
    london = {}
    for tag, path in (("pow", LEM07_L_POW), ("mc", LEM07_L_MC)):
        d = load_slope_data(path)
        crit = _lem07_search(d, "spencer", LEM07_A_SLICES)
        _lem07_solution("lem07_london_%s.png" % tag, d, crit)
        london[tag] = crit

    # ---- Part B: Low's layered undrained slope ---------------------------- #
    b = load_slope_data(LEM07_B)
    capture("lem07_low_inputs.png", plot_inputs, b,
            title="Slope Geometry and Inputs")

    crit_cp = _lem07_search(b, "bishop", LEM07_B_SLICES)
    _lem07_solution("lem07_low_cp.png", b, crit_cp)

    const = dict(b["materials"][2])
    const.update({"option": "mc", "c": LEM07_B_CONST})
    b_const = dict(b, materials=list(b["materials"][:2]) + [const])
    crit_const = _lem07_search(b_const, "bishop", LEM07_B_SLICES)
    _lem07_solution("lem07_low_const.png", b_const, crit_const)

    print("   A power  %s" % _lem07_reading(crit_pow))
    print("   A M-C    %s" % _lem07_reading(crit_mc))
    print("   London   power %.6f · M-C %.6f"
          % (london["pow"]["FS"], london["mc"]["FS"]))
    print("   B c/p    %s" % _lem07_reading(crit_cp, mat=3))
    print("   B su=%.1f %s" % (LEM07_B_CONST, _lem07_reading(crit_const, mat=3)))


# --------------------------------------------------------------------------- #
# LEM-11 — Reliability (open-and-run, one edit)
#
# One model, three analyses of it: the deterministic search, the Taylor series over
# the same search, and a Monte Carlo campaign on the surface that search found. The
# page's edit halves the standard deviation the variance Pareto measures as dominant,
# so the last two figures are the same Monte Carlo histogram before and after.
# --------------------------------------------------------------------------- #
LEM11 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_reliability.xlsx")
LEM11_METHOD = "spencer"
#: The page's edit: s(c) halved, 100 -> 50 psf.
LEM11_SIGMA_C = 50.0


def _lem11_sigma_c(model, value):
    """The model with the clay's s(c) set to ``value`` — the page's one edit."""
    mats = [dict(m) for m in model["materials"]]
    mats[0]["sigma_c"] = value
    return dict(model, materials=mats)


def _lem11_reliability(model, engine):
    """One reliability run, quiet. ``engine`` is 'taylor', 'mc' or 'rs'; all are
    called with the dialog's own defaults, so the figures carry what Studio's Run
    produces."""
    from xslope.reliability import (reliability_mc, reliability_rs,
                                    reliability_taylor)

    fn = {"taylor": reliability_taylor, "mc": reliability_mc,
          "rs": reliability_rs}[engine]
    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = fn(copy.deepcopy(model), LEM11_METHOD, debug_level=-1)
    if not ok:
        raise RuntimeError("lem11 %s reliability failed: %s" % (engine, res))
    return res


def lem11_plots():
    """The submerged clay slope's deterministic answer and the two estimators of
    what it is worth.

    The deterministic figure is a solution rather than a search plot: the surface
    the search settles on is the one the whole page is about, because both
    estimators evaluate the uncertainty on it. The Taylor figure is that surface
    with the F+ / F- perturbation surfaces over it — the 1 + 2N searches the method
    runs, drawn — and the Pareto beside it is where the dominant standard deviation
    is measured rather than guessed. The two histograms are the same Monte Carlo
    campaign before and after the page's edit, so the narrowing is read off one pair
    of pictures.
    """
    from xslope.plot import (plot_mc_rank_correlation, plot_reliability_histogram,
                             plot_reliability_results, plot_variance_pareto)
    from xslope.sensitivity import mc_rank_correlation, variance_contribution

    sd = load_slope_data(LEM11)
    capture("lem11_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, _, _ = circular_search(
            copy.deepcopy(sd), LEM11_METHOD, num_slices=40, diagnostic=False,
            **file_search_window(sd))
    crit = fs_cache[0]
    capture("lem11_solution.png", plot_solution, sd, crit["slices"],
            crit["failure_surface"], crit["solver_result"])

    taylor = _lem11_reliability(sd, "taylor")
    capture("lem11_taylor.png", plot_reliability_results, sd, taylor)

    with contextlib.redirect_stdout(io.StringIO()):
        ok, var = variance_contribution(copy.deepcopy(sd), method=LEM11_METHOD)
    if not ok:
        raise RuntimeError("lem11 variance contribution failed: %s" % var)
    capture("lem11_variance.png", plot_variance_pareto, var)
    with contextlib.redirect_stdout(io.StringIO()):
        ok_rank, rank = mc_rank_correlation(copy.deepcopy(sd), method=LEM11_METHOD, search=True)
    if not ok_rank:
        raise RuntimeError("lem11 rank correlation failed: %s" % rank)
    capture("lem11_rank.png", plot_mc_rank_correlation, rank)

    mc = _lem11_reliability(sd, "mc")
    capture("lem11_mc.png", plot_reliability_histogram, mc)

    rs = _lem11_reliability(sd, "rs")
    capture("lem11_rs.png", plot_reliability_histogram, rs)

    tight = _lem11_sigma_c(sd, LEM11_SIGMA_C)
    taylor_t = _lem11_reliability(tight, "taylor")
    mc_t = _lem11_reliability(tight, "mc")
    capture("lem11_mc_tightened.png", plot_reliability_histogram, mc_t)

    s = crit["slices"]
    print("   deterministic  FS %.6f · Xo %.4f Yo %.4f · tangent elev %.4f · "
          "L %.2f · ΣW %.1f"
          % (crit["FS"], crit["Xo"], crit["Yo"], crit["Depth"],
             s["dl"].sum(), s["w"].sum()))
    for tag, r in (("s(c) = 100", taylor), ("s(c) = %g" % LEM11_SIGMA_C, taylor_t)):
        print("   Taylor %-12s F_MLV %.6f · sigma_F %.6f · COV_F %.6f · beta_ln "
              "%.6f · R %.4f%% · Pf %.4f%%"
              % (tag, r["F_MLV"], r["sigma_F"], r["COV_F"], r["beta_ln"],
                 r["reliability"] * 100, r["prob_failure"] * 100))
        for p in r["param_info"]:
            print("      %-6s mlv %.4f sigma %.4f · F+ %.6f · F- %.6f · dF %.6f"
                  % (p["param"], p["mlv"], p["std"], p["F_plus"], p["F_minus"],
                     p["delta_F"]))
    for b in var["bars"]:
        print("   variance %-14s %.4f%% of Var(FS)  (cumulative %.4f%%)"
              % (b["label"], b["pct"], b["cumulative"]))
    for tag, r in (("s(c) = 100", mc), ("s(c) = %g" % LEM11_SIGMA_C, mc_t)):
        print("   MC     %-12s F_MLV %.6f · mean %.6f · sigma_F %.6f · beta_n "
              "%.6f · beta_ln %.6f · Pf %.4f%% (%d of %d below 1)"
              % (tag, r["F_MLV"], r["mean_FS"], r["sigma_F"], r["beta_normal"],
                 r["beta_ln"], r["pf_empirical"] * 100,
                 round(r["pf_empirical"] * r["n_valid"]), r["n_valid"]))
    # The second-order term the first-order Taylor series drops by construction:
    # the mean of a parameter's own F+ / F- pair against F_MLV. A parameter FS is
    # linear in contributes nothing; the curvature is what the MC mean picks up.
    for p in taylor["param_info"]:
        print("   curvature %-6s (F+ + F-)/2 - F_MLV = %+.6f"
              % (p["param"], (p["F_plus"] + p["F_minus"]) / 2 - taylor["F_MLV"]))
    print("   MC mean - F_MLV = %+.6f" % (mc["mean_FS"] - taylor["F_MLV"]))


# --------------------------------------------------------------------------- #
# LEM-12 — Piles (open-and-explore)
#
# The model ships with H blank on both piles, so every run below computes the pile
# force from the diameter and the spacing by Ito & Matsui, at the depth each trial
# surface reaches at each pile. The page's comparisons are therefore of two kinds,
# and the producer keeps them apart: SEARCHES, where the surface is free to move
# and the force moves with it, and HELD solves on the auto search's own critical
# circle, where the one variable named in the prose is the only thing that changes.
# Nothing here writes to the file — every edit is made on an in-memory copy.
# --------------------------------------------------------------------------- #
LEM12 = os.path.join(REPO_ROOT, "docs/lem/files/xslope_piles.xlsx")
LEM12_METHOD = "spencer"
LEM12_SLICES = 40
#: The page's specified-force edit: the two forces the auto run develops at its own
#: critical surface, to the digits the report's slice table prints them at.
LEM12_H = (2540.7, 1827.0)
#: The page's spacing edit: 6 ft out to 12 ft, S/D from 3 to 6.
LEM12_S = 12.0
#: The held spacing sweep, in feet. D = 2, so this is S/D from 1.5 to 8 — either
#: side of the 2-8 band Ito & Matsui is derived for.
LEM12_SPACINGS = (3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0)


def _lem12_piles(model, H=None, **over):
    """``model`` with every pile row overridden. ``H`` is a per-pile sequence."""
    rows = []
    for i, p in enumerate(model["pile_lines"]):
        row = dict(p, **over)
        if H is not None:
            row["H"] = H[i]
        rows.append(row)
    return dict(model, pile_lines=rows)


def _lem12_search(model, **kwargs):
    """The page's search: Spencer from the circle the file carries, 40 slices."""
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, _, _ = circular_search(
            copy.deepcopy(model), LEM12_METHOD, num_slices=LEM12_SLICES,
            diagnostic=False, **dict(file_search_window(model), **kwargs))
    return fs_cache[0]


def _lem12_held(model, circle):
    """One solve on a stated circle — no search. Returns ``(FS, slice_df)``."""
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected

    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = generate_slices(model, circle=circle, num_slices=LEM12_SLICES)
        if not ok:
            return None, res
        slice_df, _ = res
        result = solve_selected(LEM12_METHOD, slice_df)
    if isinstance(result, str):
        return None, result
    return result["FS"], slice_df


def _lem12_reading(crit):
    """FS, the circle, the surface and the pile force the page quotes off a search."""
    s = crit["slices"]
    nz = s[s["h_pile"] > 0]
    return ("FS %.6f · Xo %.4f Yo %.4f R %.4f · tangent elev %.4f · L %.2f · "
            "ΣW %.1f · ΣH %.4f (%s)"
            % (crit["FS"], crit["Xo"], crit["Yo"], crit["Yo"] - crit["Depth"],
               crit["Depth"], s["dl"].sum(), s["w"].sum(), nz["h_pile"].sum(),
               " + ".join("%.4f" % v for v in nz["h_pile"])))


def lem12_plots():
    """The pile-stabilized slope's five searches and the two held studies.

    The first search is the file exactly as it downloads — H blank on both piles,
    so the force comes from Ito & Matsui at whatever depth each trial surface
    reaches. Both the search plot and the solution are drawn: the trial family is
    where the force is being recomputed, and the solution is where the two pile
    crossings are marked.

    The other four searches are each one edit away from it, and each is drawn
    because its critical surface is somewhere else: the piles taken out, the
    spacing widened to 12 ft, the force stated instead of computed, and the same
    model searched from a grid of seeds rather than from the circle on the sheet.

    The two held studies — the spacing sweep and the structural caps — change one
    input on the auto search's own critical circle and are printed rather than
    drawn, because what they measure is a column of numbers.
    """
    from xslope.ito_matsui import (compute_ito_matsui_force_and_moment_arm,
                                   intersect_pile_with_materials,
                                   ito_matsui_coefficients)

    sd = load_slope_data(LEM12)
    capture("lem12_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path, cache = circular_search(
            copy.deepcopy(sd), LEM12_METHOD, num_slices=LEM12_SLICES,
            diagnostic=False, **file_search_window(sd))
    crit = fs_cache[0]
    capture("lem12_search.png", plot_circular_search_results, sd, fs_cache, path,
            circle_cache=cache)
    capture("lem12_solution.png", plot_solution, sd, crit["slices"],
            crit["failure_surface"], crit["solver_result"])
    deep = {"Xo": crit["Xo"], "Yo": crit["Yo"], "Depth": crit["Depth"],
            "R": crit["Yo"] - crit["Depth"]}

    bare = _lem12_search(dict(sd, pile_lines=[]))
    capture("lem12_solution_nopiles.png", plot_solution, dict(sd, pile_lines=[]),
            bare["slices"], bare["failure_surface"], bare["solver_result"])

    wide = _lem12_piles(sd, S=LEM12_S)
    crit_w = _lem12_search(wide)
    capture("lem12_solution_wide.png", plot_solution, wide, crit_w["slices"],
            crit_w["failure_surface"], crit_w["solver_result"])

    stated = _lem12_piles(sd, H=LEM12_H)
    crit_h = _lem12_search(stated)
    capture("lem12_solution_statedh.png", plot_solution, stated, crit_h["slices"],
            crit_h["failure_surface"], crit_h["solver_result"])

    grid = _lem12_search(sd, seed="grid")
    capture("lem12_solution_bypass.png", plot_solution, sd, grid["slices"],
            grid["failure_surface"], grid["solver_result"])

    print("   auto      %s" % _lem12_reading(crit))
    print("   no piles  %s" % _lem12_reading(bare))
    print("   S = %-5g %s" % (LEM12_S, _lem12_reading(crit_w)))
    print("   stated H  %s" % _lem12_reading(crit_h))
    print("   grid seed %s" % _lem12_reading(grid))

    # What the two piles are worth on ONE surface — the auto search's own critical
    # circle, so the only thing changing is which piles are present.
    for tag, rows in (("neither", []), ("pile 1", sd["pile_lines"][:1]),
                      ("pile 2", sd["pile_lines"][1:]), ("both", sd["pile_lines"])):
        fs, _ = _lem12_held(dict(sd, pile_lines=list(rows)), deep)
        print("   held, %-8s FS %.6f" % (tag, fs))

    # The stated force against the computed one, on two stated circles: the auto
    # search's own, where they agree, and the shallower circle the stated-force
    # search settles on, where the frozen number is a force for the wrong depth.
    shallow = {"Xo": crit_h["Xo"], "Yo": crit_h["Yo"], "Depth": crit_h["Depth"],
               "R": crit_h["Yo"] - crit_h["Depth"]}
    for label, circle in (("deep", deep), ("shallow", shallow)):
        for tag, model in (("computed", sd), ("stated", stated),
                           ("no piles", dict(sd, pile_lines=[]))):
            fs, df = _lem12_held(model, circle)
            nz = df[df["h_pile"] > 0] if fs is not None else None
            print("   held %-8s %-9s FS %.6f · ΣH %.4f"
                  % (label, tag, fs, nz["h_pile"].sum() if nz is not None else 0))

    # The spacing sweep, held on the deep circle: the arching coefficients, the
    # uncapped soil force per pile, and what survives the moment cap.
    gs = sd["ground_surface"].coords
    xs = [p[0] for p in gs]
    ys = [p[1] for p in gs]
    for S in LEM12_SPACINGS:
        model = _lem12_piles(sd, S=S)
        fs, df = _lem12_held(model, deep)
        nz = df[df["h_pile"] > 0]
        A1, A2 = ito_matsui_coefficients(S - model["pile_lines"][0]["D_pile"],
                                         model["pile_lines"][0]["D_pile"],
                                         sd["materials"][0]["phi"])
        parts = []
        for pile, (_, ps) in zip(model["pile_lines"], nz.iterrows()):
            y_g = _interp(ps["x_pile"], xs, ys)
            segs = intersect_pile_with_materials(ps["x_pile"], y_g, ps["y_pile"],
                                                 model["polygons"], model["materials"])
            _H, F, L_m = compute_ito_matsui_force_and_moment_arm(pile["D_pile"], S, segs)
            parts.append("z %.3f F_pile %.1f L_m %.3f M/L_m %.1f"
                         % (y_g - ps["y_pile"], F, L_m,
                            pile["M_cap"] / L_m if L_m else float("inf")))
        print("   S = %-5g S/D %.1f · A1 %.4f A2 %.4f · %s · ΣH %.4f · FS %.6f"
              % (S, S / model["pile_lines"][0]["D_pile"], A1, A2, " | ".join(parts),
                 nz["h_pile"].sum(), fs))

    # The structural caps, held on the same circle: which of the two governs, and
    # what the uncapped Ito & Matsui force alone would have read.
    for tag, over in (("as shipped", {}), ("V_cap only", {"M_cap": None}),
                      ("M_cap only", {"V_cap": None}),
                      ("neither cap", {"V_cap": None, "M_cap": None})):
        fs, df = _lem12_held(_lem12_piles(sd, **over), deep)
        if fs is None:
            print("   caps %-12s no solution: %s" % (tag, df))
            continue
        nz = df[df["h_pile"] > 0]
        print("   caps %-12s FS %.6f · ΣH %.4f (%s)"
              % (tag, fs, nz["h_pile"].sum(),
                 " + ".join("%.1f" % v for v in nz["h_pile"])))


def _interp(x, xs, ys):
    """Ground elevation at ``x`` along the ground surface's own vertices."""
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            if xs[i + 1] == xs[i]:
                return ys[i]
            t = (x - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + t * (ys[i + 1] - ys[i])
    return ys[-1] if x > xs[-1] else ys[0]


# --------------------------------------------------------------------------- #
# SEEP-1 — Seepage Under a Sheetpile (built from scratch, one model, many meshes)
#
# The model is confined: two specified-head boundaries, no exit face, one soil. That
# makes every number on the page a property of the mesh and of the boundary values
# alone, and the producer is organized around exactly that. Four states are drawn —
# the inputs, the mesh the tutorial's first run uses, that run's flow net, and the
# same flow net on a mesh auto-refined at the sheetpile tip — and three studies are
# printed: the discharge against the element size, the discharge against the element
# type, and the discharge against the conductivity.
#
# The tip is the reason the page exists. A sheetpile modeled as a slot in the ground
# surface leaves a re-entrant corner at its toe, where the hydraulic gradient is
# singular; the discharge therefore keeps falling as the mesh refines, and every
# quoted discharge names the mesh it was computed on.
# --------------------------------------------------------------------------- #
SEEP01 = os.path.join(REPO_ROOT, "docs/seep/files/xslope_clay_blanket.xlsx")
#: The tutorial's first mesh: tri3 (the seepage default), at the target size Studio's
#: own auto-size produces on this 50 m section — width / 100 divisions.
SEEP01_SIZE = 0.5
SEEP01_ELEMENT = "tri3"
#: The refinement the page turns on, in the Build mesh dialog's own terms: local
#: element size = target size / factor near model features, the sheetpile slot's toe
#: among them.
SEEP01_REFINE = 4
#: Element sizes for the mesh study, coarse to fine. 0.4167 is width / 120 — the size
#: the sample page's regression tag meshes at, and the reason its catalogued discharge
#: is neither of the neighbouring rows.
SEEP01_SIZES = (2.0, 1.0, 50.0 / 120.0, 0.5, 0.25, 0.125)
#: The size the element-type comparison is made at. Coarser than the page's own mesh,
#: so the difference between the element orders is not buried under a mesh fine enough
#: to hide it.
SEEP01_TYPE_SIZE = 1.0
#: Contour count for the flow net. 9 head drops is what makes the channel count come
#: out at a whole number here, so the drawn net is a flow net and not a rounding of one.
SEEP01_LEVELS = 10
#: The conductivity sweep, in the Parametric dialog's own From / To / Steps: from the
#: model's own 30 m/day up to ten times it. The first swept value is the model's own,
#: so both series start from the same measured point.
SEEP01_K_LOW, SEEP01_K_HIGH, SEEP01_K_STEPS = 30.0, 300.0, 10
#: The design target the sweep is asked to locate, a discharge in the model's units.
SEEP01_TARGET_Q = 100.0
#: The upstream heads the head-drop sweep uses, against a downstream head of 10 —
#: the other factor of q = k·Δh·Nf/Nd, checked the same way the conductivity is.
SEEP01_UPSTREAM_HEADS = (11.0, 12.0, 13.0, 14.0, 15.0)
#: The sheetpile toe, and the upstream edge of the clay blanket — the model's two
#: re-entrant corners, and the two places the gradient is measured.
SEEP01_TIP = (30.0, 7.0)
SEEP01_BLANKET_EDGE = (20.0, 10.0)


def seep01_sheets():
    """The four worksheets SEEP-1's Excel path fills.

    ``main`` and ``profile`` take LEM-1's own windows on those sheets, so a reader
    coming from that page sees the same frame twice — through the Water-loads row on
    one, and through the empty Profile Line #2 beside the filled #1 on the other.

    ``mat`` cannot: this problem fills the seepage band at columns AG–AP and nothing
    in the strength band at all, and a window wide enough to hold both would be forty
    columns across. It is rendered as the split view the input-template manifest uses
    for the same band — the seepage columns, with the mat ID and name re-shown at the
    left so the row is identifiable. Row 10 rather than 9 starts the window under the
    sheet's merged *Seepage* band header, which a frame ending inside would raise on.

    ``seep bc`` runs through the empty BC #3 beside the two filled blocks, for the
    reason LEM-1's profile window does: the empty neighbor is what the sheet looks
    like, and cropping it makes the capture read as a table rather than a worksheet.
    """
    render("seep01_sheet_main.png", SEEP01, "main", rows=(1, 24), cols="A:D")
    render("seep01_sheet_mat.png", SEEP01, "mat", rows=(10, 12), cols="AG:AP",
           identity_cols="A:B")
    render("seep01_sheet_profile.png", SEEP01, "profile", rows=(1, 15), cols="A:H")
    render("seep01_sheet_seep_bc.png", SEEP01, "seep bc", rows=(1, 8), cols="A:L")


def _seep01_mesh(model, size=SEEP01_SIZE, element_type=SEEP01_ELEMENT, **kwargs):
    """One mesh of the model's material polygons, quiet."""
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)

    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(
            get_material_polygons(model), size, element_type,
            size_regions=extract_size_regions(model), **kwargs)


def _seep01_solve(model, mesh):
    """One steady seepage solve on a built mesh. Returns ``(seep_data, solution)``."""
    from xslope.seep import build_seep_data, run_seepage_analysis

    with contextlib.redirect_stdout(io.StringIO()):
        seep_data = build_seep_data(mesh, model)
        solution = run_seepage_analysis(seep_data, tol=1e-4)
    return seep_data, solution


def _seep01_run(model, **kwargs):
    """Mesh and solve in one step, for the studies that only want numbers."""
    mesh = _seep01_mesh(model, **kwargs)
    seep_data, solution = _seep01_solve(model, mesh)
    return mesh, seep_data, solution


def _seep01_material(model, **fields):
    """``model`` with its one soil's seepage fields overridden, in memory only."""
    mats = [dict(m) for m in model["materials"]]
    mats[0].update(fields)
    return dict(model, materials=mats)


def _seep01_gradient_near(mesh, solution, point, radius=0.5):
    """The largest hydraulic-gradient magnitude within ``radius`` of ``point``, and
    how many nodes the mesh puts there.

    A re-entrant corner has no finite gradient, so this is a mesh measurement as much
    as a flow one: it says how hard the mesh is resolving the corner, which is what
    the page's refinement section compares.
    """
    import numpy as np

    nodes = np.asarray(mesh["nodes"])
    dist = np.linalg.norm(nodes - np.asarray(point), axis=1)
    near = dist < radius
    i_mag = np.asarray(solution["i_mag"])
    return float(i_mag[near].max()), int(near.sum()), float(np.sort(dist)[1])


def _seep01_tip_zoom(panels, half=1.6):
    """The meshes of ``panels`` — ``(label, mesh)`` pairs — side by side around the
    sheetpile toe.

    One figure rather than two, because the comparison IS the subject: the same
    window, the same axes and the same element edges, so the only difference on the
    page is the thing the refinement changed. The window is square and centered on the
    toe, and each panel carries the node count its mesh puts inside it.
    """
    import numpy as np

    x0, y0 = SEEP01_TIP
    fig, axes = plt.subplots(1, len(panels), figsize=(12, 5.2))
    for ax, (label, mesh) in zip(axes, panels):
        nodes = np.asarray(mesh["nodes"])
        tris = [e[:3] for e in np.asarray(mesh["elements"])]
        ax.triplot(nodes[:, 0], nodes[:, 1], tris, color="#5a6b7a", lw=0.5)
        inside = ((np.abs(nodes[:, 0] - x0) <= half)
                  & (np.abs(nodes[:, 1] - y0) <= half))
        ax.plot([29.9, x0, 30.1], [10.0, y0, 10.0], color="#b03030", lw=2.0,
                solid_capstyle="round")
        ax.plot([x0], [y0], "o", color="#b03030", ms=5)
        ax.set_xlim(x0 - half, x0 + half)
        ax.set_ylim(y0 - half, y0 + half)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x (m)")
        ax.set_title("%s — %d nodes in this window" % (label, int(inside.sum())))
    axes[0].set_ylabel("y (m)")
    fig.suptitle("The mesh at the sheetpile toe")
    fig.tight_layout()


def _seep01_q_vs_k(both, k1_only, target=None, crossing=None):
    """The discharge against conductivity, two sweeps on one pair of axes.

    Both series start at the model's own 30 m/day, so the figure's whole content is
    how they diverge above it. Scaling both principal conductivities together leaves
    the head field untouched and makes the discharge exactly proportional to k;
    sweeping k1 alone changes the anisotropy ratio, and with it the head field, so
    that series bends away from the straight one. Linear axes rather than log,
    because a straight line through the origin is what proportional LOOKS like and a
    log plot would draw both series straight.

    The proportional reference is drawn wide and the computed series narrow on top of
    it, so the series that is proportional reads as markers centered in a gray band
    rather than as a blue line hiding a gray one underneath.
    """
    ks = [k for k, _ in both]
    ratio = both[0][1] / both[0][0]

    fig, ax = plt.subplots(figsize=(9, 5.4))
    span = [0.0, max(ks)]
    ax.plot(span, [ratio * k for k in span], "-", color="#b8c1cb", lw=5.0,
            solid_capstyle="butt", label="q proportional to k")
    ax.plot(ks, [q for _, q in both], "o-", color="#1f6fb4", lw=1.2, ms=6,
            label="k₁ and k₂ scaled together")
    ax.plot([k for k, _ in k1_only], [q for _, q in k1_only], "s-", color="#c1663a",
            lw=1.4, ms=5, label="k₁ swept, k₂ held at 30 m/day")
    if target is not None:
        ax.axhline(target, color="#7a8592", lw=0.9, ls="--")
        ax.annotate("target q = %g" % target, (0.0, target), xytext=(4, 4),
                    textcoords="offset points", color="#5b646f", fontsize=9)
        if crossing is not None:
            ax.plot([crossing], [target], "*", color="#c1663a", ms=13)
            ax.annotate("k₁ = %.1f" % crossing, (crossing, target), xytext=(6, -14),
                        textcoords="offset points", color="#c1663a", fontsize=9)
    ax.set_xlim(0, max(ks) * 1.02)
    ax.set_ylim(0, None)
    ax.set_xlabel("hydraulic conductivity (m/day)")
    ax.set_ylabel("total discharge q (m³/day per m)")
    ax.set_title("Discharge against hydraulic conductivity")
    ax.grid(True, color="#e3e7eb", lw=0.6)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()


def seep01_plots():
    """The confined sheetpile problem: four figures and three studies.

    The figures are the four states the page asks the reader to compare their own
    screen against — the inputs as the boundary conditions leave them, the first
    mesh with its boundary nodes marked, the flow net that mesh produces, and the
    flow net after the tip refinement. The gradient field is drawn as well, because
    it is what says WHERE the mesh has work to do, and the two meshes are drawn
    together at the toe, because the refinement is otherwise invisible at section
    scale.

    The three studies are printed rather than tabulated in a figure: what the
    discharge does as the element size falls, what it does as the element type
    changes at one size, and what it does as the conductivity moves. Every number
    the page quotes is on one of these lines.
    """
    import numpy as np

    from xslope.plot_seep import plot_seep_data, plot_seep_solution

    sd = load_slope_data(SEEP01)
    # frame="content" because the section is five times as long as it is deep: the
    # default fill frame pads the data limits out to the figure's aspect and buries
    # the model in empty sky. show_mesh=False because this figure is the INPUTS —
    # the model as it stands before the reader has built a mesh at all, whatever
    # mesh the loader found beside the file.
    capture("seep01_inputs.png", plot_inputs, sd, mode="seep",
            title="Seepage Model Inputs", frame="content", show_mesh=False)

    mesh = _seep01_mesh(sd)
    seep_data, solution = _seep01_solve(sd, mesh)
    capture("seep01_mesh.png", plot_seep_data, seep_data, show_bc=True)
    # The solution as the view first shows it — the Display panel's defaults
    # (20 contour levels, Filled contours off) — and the same state with the fill
    # toggled on. The page walks the reader through both before setting the level
    # count to SEEP01_LEVELS for the flow-net reading.
    capture("seep01_solution_arrival.png", plot_seep_solution, seep_data, solution,
            levels=20, base_mat=1, fill_contours=False, mesh=False)
    capture("seep01_solution_filled.png", plot_seep_solution, seep_data, solution,
            levels=20, base_mat=1, fill_contours=True, mesh=False)
    capture("seep01_solution.png", plot_seep_solution, seep_data, solution,
            levels=SEEP01_LEVELS, base_mat=1, fill_contours=False, mesh=False)
    capture("seep01_gradient.png", plot_seep_solution, seep_data, solution,
            levels=SEEP01_LEVELS, variable="i_mag", mesh=False, flowlines=False)

    head = np.asarray(solution["head"])
    hdrop = float(head.max() - head.min())
    k = float(sd["materials"][0]["k1"])
    drops = SEEP01_LEVELS - 1
    channels = max(round(solution["flowrate"] * drops / (k * hdrop)) + 1, 2) - 1
    print("   base mesh   size %.4g %s · %d nodes · %d elements · q %.4f"
          % (SEEP01_SIZE, SEEP01_ELEMENT, len(mesh["nodes"]),
             len(mesh["elements"]), solution["flowrate"]))
    print("   head        %.4f to %.4f (drop %.4f) · u %.3f to %.3f · confined=%s "
          "converged=%s"
          % (head.min(), head.max(), hdrop, np.min(solution["u"]),
             np.max(solution["u"]), not solution.get("unconfined"),
             solution.get("converged")))
    print("   flow net    %d drops · %d channels · q = k·Δh·Nf/Nd = %.4f"
          % (drops, channels, k * hdrop * channels / drops))
    print("   60 m of wall q %.4f" % (60.0 * solution["flowrate"]))

    # Whether the section was drawn wide enough is a measurement, not a habit: if
    # the ends are far enough from the wall the head is flat there and nothing is
    # moving, which is what the no-flow ends assume.
    xs = np.asarray(mesh["nodes"])[:, 0]
    i_mag = np.asarray(solution["i_mag"])
    for label, near in (("upstream", xs <= 5.0), ("downstream", xs >= 45.0)):
        print("   %-10s end, outer 5 m: max |i| %.4f · head %.4f to %.4f"
              % (label, i_mag[near].max(), head[near].min(), head[near].max()))

    # The mesh study. The discharge falls monotonically, because the two re-entrant
    # corners have no finite gradient and a coarse mesh cannot see how much of the
    # flow is being forced through them.
    print("   -- element size, %s" % SEEP01_ELEMENT)
    for size in sorted(SEEP01_SIZES, reverse=True):
        m, _, s = _seep01_run(sd, size=size)
        print("   size %-7.4g %6d nodes %6d elements · q %.4f"
              % (size, len(m["nodes"]), len(m["elements"]), s["flowrate"]))

    # The element-type study, at one size, so the only variable is the order of the
    # element and the node count that comes with it.
    print("   -- element type, size %.4g" % SEEP01_TYPE_SIZE)
    for etype in ("tri3", "tri6", "quad4", "quad8", "quad9"):
        m, _, s = _seep01_run(sd, size=SEEP01_TYPE_SIZE, element_type=etype)
        print("   %-6s %6d nodes %6d elements · q %.4f"
              % (etype, len(m["nodes"]), len(m["elements"]), s["flowrate"]))

    # The tip refinement, against the uniform mesh it starts from and against the
    # uniform mesh that reaches the same answer the expensive way.
    mesh_r = _seep01_mesh(sd, refine_factor=SEEP01_REFINE)
    seep_r, solution_r = _seep01_solve(sd, mesh_r)
    capture("seep01_solution_refined.png", plot_seep_solution, seep_r, solution_r,
            levels=SEEP01_LEVELS, base_mat=1, fill_contours=False, mesh=False)
    capture("seep01_tip.png", _seep01_tip_zoom,
            [("Uniform, target size %.4g m" % SEEP01_SIZE, mesh),
             ("Refined near features, factor %g" % SEEP01_REFINE, mesh_r)])
    print("   -- refinement at the toe")
    for label, m, s in (("uniform", mesh, solution), ("refined", mesh_r, solution_r)):
        for name, point in (("toe", SEEP01_TIP), ("blanket edge", SEEP01_BLANKET_EDGE)):
            i_max, count, nearest = _seep01_gradient_near(m, s, point)
            print("   %-8s %-13s max |i| %.4f · %3d nodes within 0.5 m · nearest "
                  "node %.4f m" % (label, name, i_max, count, nearest))
        print("   %-8s %6d nodes %6d elements · q %.4f"
              % (label, len(m["nodes"]), len(m["elements"]), s["flowrate"]))

    # The conductivity sweep, twice over the same values: the whole conductivity
    # tensor scaled, and k1 alone. The second is run through xslope.sensitivity's own
    # design(), on the same mesh, with the bounds, step count and target the page
    # dictates into the Parametric dialog — so the crossing the page quotes is the
    # one that dialog reports rather than a number computed a private way.
    from xslope.sensitivity import design

    print("   -- discharge against conductivity, size %.4g %s"
          % (SEEP01_SIZE, SEEP01_ELEMENT))
    ks = list(np.linspace(SEEP01_K_LOW, SEEP01_K_HIGH, SEEP01_K_STEPS))
    both = []
    base_head = None
    for kk in ks:
        _, s = _seep01_solve(_seep01_material(sd, k1=kk, k2=kk), mesh)
        drift = (0.0 if base_head is None
                 else float(np.abs(np.asarray(s["head"]) - base_head).max()))
        base_head = np.asarray(s["head"]) if base_head is None else base_head
        both.append((kk, s["flowrate"]))
        print("   k1 = k2 = %-6g q %12.6f · q/k %.8f · max head change from the "
              "first row %.2e" % (kk, s["flowrate"], s["flowrate"] / kk, drift))

    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = design(dict(sd, mesh=mesh), param="seep:%s:k1"
                         % sd["materials"][0]["name"],
                         low=SEEP01_K_LOW, high=SEEP01_K_HIGH, steps=SEEP01_K_STEPS,
                         target_fs=SEEP01_TARGET_Q, mode="seep",
                         seep_opts={"bc": 1, "tol": 1e-4})
    if not ok:
        print("   design sweep refused: %s" % res)
        return
    # design() evaluates the model's own unmodified value as well as the swept ones,
    # and here that value IS the first swept point, so the row appears twice.
    swept = (res["df"][res["df"]["success"]]
             .drop_duplicates(subset="value").sort_values("value"))
    k1_only = [(float(r.value), float(r.fs)) for r in swept.itertuples()
               if SEEP01_K_LOW <= r.value <= SEEP01_K_HIGH]
    for kk, q in k1_only:
        print("   k1 = %-6g k2 = 30  q %12.6f · q/k1 %.6f" % (kk, q, q / kk))
    print("   design      %s · bracketed %s · %s"
          % (res["message"], res["bracketed"], res["direction"]))
    capture("seep01_q_vs_k.png", _seep01_q_vs_k, both, k1_only,
            target=SEEP01_TARGET_Q, crossing=res["crossing"])

    # The other half of q = k·Δh·Nf/Nd, measured the same way: the head drop, with
    # the conductivity left alone. The page states this one rather than drawing it,
    # so it is printed only.
    print("   -- discharge against the head drop, size %.4g %s"
          % (SEEP01_SIZE, SEEP01_ELEMENT))
    for upstream in SEEP01_UPSTREAM_HEADS:
        model = copy.deepcopy(sd)
        model["seepage_bc"]["specified_heads"][0]["head"] = upstream
        _, s = _seep01_solve(model, mesh)
        drop = upstream - model["seepage_bc"]["specified_heads"][1]["head"]
        print("   upstream head %-5g drop %-4g q %12.6f · q/Δh %.8f"
              % (upstream, drop, s["flowrate"], s["flowrate"] / drop))


# --------------------------------------------------------------------------- #
# SEEP-2 — Unconfined Seepage Through a Zoned Dam (open and explore)
#
# The Johnson Reservoir dam of docs/seep/samples.md #5, meshed and solved exactly
# as run_tests.py::run_seep_test solves it — tri3 at the ground-surface width over
# 120 — so the discharge this page quotes is the discharge that page's lock and the
# SEEP2D cross-check were recorded from.
#
# Everything here is a comparison, because the page is about choices rather than
# about one answer: the flow net drawn against each of the three material zones in
# turn, the three unsaturated conductivity models run on the same mesh, the linear
# front's floor swept until it reaches the van Genuchten answer, the convergence
# histories of all four side by side, and the core's conductivity raised until the
# seepage face on the downstream slope has something to do.
# --------------------------------------------------------------------------- #
SEEP02 = os.path.join(REPO_ROOT, "docs/seep/files/xslope_johnson_res.xlsx")
#: The mesh, in the Build mesh dialog's own controls: tri3, auto-sized at 120
#: divisions across the 750 ft section, which is the 6.25 ft target size the sample
#: page's regression tag and the SEEP2D comparison were both computed on.
SEEP02_DIVISIONS = 120
SEEP02_ELEMENT = "tri3"
#: Contour count for the flow net, the same 20 the sample figures are drawn at.
SEEP02_LEVELS = 20
#: The flow net's base material, 1-based into the mat sheet: 3 is the foundation,
#: which is what flownet_base_material picks for this model and what the sample
#: figure was drawn with. The page draws all three to show why.
SEEP02_BASE_MAT = 3
#: van Genuchten α and n per material, after Carsel & Parrish (1988) by texture —
#: sandy clay loam for the shell, clay for the core, clay loam for the foundation —
#: with α converted from the paper's 1/cm to this model's 1/ft by ×30.48.
SEEP02_VG = ((1.798, 1.48), (0.244, 1.09), (0.579, 1.31))
#: Gardner a and n per material, least-squares fits to each material's own van
#: Genuchten curve in log kr over 0.01 to 100 ft of suction. Gardner has no texture
#: table to read off — its parameters arrive with an imported model or from fitted
#: measurements — so fitting them to the van Genuchten curve is what makes the
#: three-model comparison a comparison of models rather than of soils. The producer
#: prints the misfit of these pinned pairs, so the calibration claim is measured.
SEEP02_GARD = ((115.5, 2.29), (128.2, 1.03), (52.8, 1.61))
#: The linear front's floor, swept down from its shipped 0.01 toward the 1e-4 the
#: other two models floor at. This is the test of the page's explanation for why the
#: linear front passes more water than the other two.
SEEP02_KR0_SWEEP = (0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001)
#: The linear-front pair that makes this dam hard: a floor four decades down reached
#: over 10 ft of suction, which is the shape of a van Genuchten curve drawn with
#: straight lines. It does not converge inside the default iteration ceiling.
SEEP02_HARD = {"kr0": 1e-4, "h0": -10.0}
SEEP02_HARD_MAX_ITER = 1000
#: Core conductivities for the seepage-face study, from the model's own 0.001 ft/day
#: up to the shell's 1.0 — a core, a poor core, a fill and no core at all.
SEEP02_CORE_KS = (0.001, 0.01, 0.1, 1.0)
#: Solve tolerances, to check the claim that the flow-closure condition makes the
#: converged discharge independent of the head tolerance it is asked for.
SEEP02_TOLS = (1e-3, 1e-4, 1e-5, 1e-6)
#: Vertical sections the flow is measured across. 370 is the dam centerline, through
#: the core and its cutoff key; 500 is in the downstream shell, where the phreatic
#: surface has dropped well below the slope and the unsaturated zone is thickest.
SEEP02_CENTERLINE = 370.0
SEEP02_DOWNSTREAM = 500.0
#: The bottom of the core's cutoff key, the elevation the underseepage passes below.
SEEP02_KEY_TOE = 60.0
#: Stations the phreatic surface is read at, upstream toe to downstream toe.
SEEP02_STATIONS = (220.0, 260.0, 300.0, 330.0, 355.0, 370.0, 390.0, 420.0,
                   460.0, 500.0, 540.0)


def _seep02_top(x):
    """Elevation of the top of the section at ``x`` — reservoir floor, upstream
    slope, crest, downstream slope, tailwater floor — from the model's own outer
    profile line."""
    if x <= 200.0:
        return 100.0
    if x <= 320.0:
        return 100.0 + (x - 200.0) * 60.0 / 120.0
    if x <= 360.0:
        return 160.0 + (x - 320.0) * 20.0 / 40.0
    if x <= 380.0:
        return 180.0
    if x <= 550.0:
        return 180.0 - (x - 380.0) * 80.0 / 170.0
    return 100.0


def _seep02_mesh(model, divisions=SEEP02_DIVISIONS, element_type=SEEP02_ELEMENT):
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)

    xs = [x for x, _ in model["ground_surface"].coords]
    size = (max(xs) - min(xs)) / divisions
    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(
            get_material_polygons(model), size, element_type,
            size_regions=extract_size_regions(model))


def _seep02_solve(model, mesh, tol=1e-4, max_iter=400):
    """One steady seepage solve, returning ``(seep_data, solution, log)``.

    The log is kept rather than discarded because on an unconfined problem it is
    half the result: the iteration count, the relaxation the solver fell back to,
    how many exit-face nodes ended up active, and — when the run does not
    converge — the sentence that says so.
    """
    from xslope.seep import build_seep_data, run_seepage_analysis

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        seep_data = build_seep_data(mesh, model)
        solution = run_seepage_analysis(seep_data, tol=tol, max_iter=max_iter)
    return seep_data, solution, buf.getvalue()


def _seep02_log_stats(log):
    """The unconfined iteration, read off the solver's own log."""
    import re

    done = re.search(r"Converged in (\d+) iterations \(residual = ([0-9.e+-]+), "
                     r"closure = ([0-9.e+-]+)", log)
    relax = [float(r) for r in re.findall(r"relax = ([0-9.]+)", log)]
    active = re.findall(r"(\d+)/(\d+) exit face active", log)
    tol = re.search(r"Convergence tolerance: ([0-9.e+-]+)", log)
    return {
        "converged": done is not None,
        "iterations": int(done.group(1)) if done else None,
        "residual": float(done.group(2)) if done else None,
        "closure": float(done.group(3)) if done else None,
        "relax_min": min(relax) if relax else None,
        "active": (int(active[-1][0]), int(active[-1][1])) if active else None,
        "tol_scaled": float(tol.group(1)) if tol else None,
    }


def _seep02_history(log):
    """``(iteration, residual, closure)`` for every sweep the log reported.

    The per-sweep line is printed only on the first three sweeps, every fifth
    sweep after that, and whenever the exit-face set changes, so a run that
    converges on a sweep matching none of those never prints its last one. The
    terminal ``Converged in N iterations`` line carries the same two numbers, so it
    is appended here: without it a converging series stops in mid-descent, above
    the thresholds it actually crossed.
    """
    import re

    rows = re.findall(r"Iteration (\d+): residual = ([0-9.e+-]+), "
                      r"closure = ([0-9.e+-]+)", log)
    out = [(int(i), float(r), float(c)) for i, r, c in rows]
    done = re.search(r"Converged in (\d+) iterations \(residual = ([0-9.e+-]+), "
                     r"closure = ([0-9.e+-]+)", log)
    if done is not None:
        last = (int(done.group(1)), float(done.group(2)), float(done.group(3)))
        if not out or out[-1][0] != last[0]:
            out.append(last)
    return out


def _seep02_active_history(log):
    """``(iteration, active nodes)`` for every sweep the log reported.

    The exit-face line prints on any sweep the set changed, so this carries every
    change even though it does not carry every sweep.
    """
    import re

    rows = re.findall(r"Iteration (\d+): .*? (\d+)/\d+ exit face active", log)
    return [(int(i), int(n)) for i, n in rows]


#: The solver's relaxation ladder (xslope/seep.py): the full step is taken through
#: sweep 20, and the factor drops at each of these sweep counts thereafter. Reading
#: the factor off the ladder rather than off the log is what makes it right on a run
#: whose last sweeps were not printed.
SEEP02_RELAX_LADDER = ((20, 0.5), (40, 0.2), (60, 0.1), (80, 0.05),
                       (100, 0.02), (120, 0.01))


def _seep02_relax_at(iteration):
    """The relaxation factor the ladder gives on sweep ``iteration``."""
    relax = 1.0
    for after, value in SEEP02_RELAX_LADDER:
        if iteration > after:
            relax = value
    return relax


def _seep02_materials(model, per_mat):
    """``model`` with each material's seepage fields overridden, in memory only.

    ``per_mat`` is one dict per material, in mat-sheet order; an empty dict leaves
    that material as the file carries it.
    """
    mats = [dict(m) for m in model["materials"]]
    for mat, fields in zip(mats, per_mat):
        mat.update(fields)
    return dict(model, materials=mats)


def _seep02_interp(mesh, values):
    """A linear interpolator over the tri3 mesh, for reading a nodal field along a
    line the mesh has no nodes on."""
    import matplotlib.tri as mtri
    import numpy as np

    nodes = np.asarray(mesh["nodes"])
    tris = [e[:3] for e in np.asarray(mesh["elements"])]
    triang = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tris)
    return mtri.LinearTriInterpolator(triang, np.asarray(values))


def _seep02_phreatic(mesh, solution, stations=SEEP02_STATIONS):
    """The elevation of the phreatic surface at each station.

    The phreatic surface is the pressure-head zero contour, so this walks a
    vertical line from the base up and returns the elevation where the pressure
    head crosses zero, interpolating between the two samples that straddle it.
    Reading it off an interpolator rather than off the nodes matters here: the
    three unsaturated models move this surface by less than a foot, and a
    nearest-node reading on a 6.25 ft mesh cannot see a difference that small.

    A column that is saturated all the way to the ground surface has no
    crossing inside the domain — on this dam that is the submerged upstream
    face, where the p = 0 level is the reservoir surface standing above the
    ground. Those stations return NaN so the drawn line starts where the
    reservoir level meets the face, the same convention as the solution plot's
    own phreatic line. (Clamping to the topmost sample instead traced the
    upstream face from the toe — the owner caught it on the published figure.)
    """
    import numpy as np

    ipsi = _seep02_interp(mesh, np.asarray(solution["u"]) / 62.4)
    out = []
    for x in stations:
        ys = np.linspace(0.02, _seep02_top(x) - 0.02, 4000)
        psi = np.ma.filled(ipsi(np.full(len(ys), x), ys), np.nan)
        ok = ~np.isnan(psi)
        ys, psi = ys[ok], psi[ok]
        wet = np.where(psi >= 0.0)[0]
        if not len(wet):
            out.append(float("nan"))
            continue
        i = wet[-1]
        if i + 1 < len(psi):
            out.append(float(ys[i] + (ys[i + 1] - ys[i]) * psi[i]
                             / (psi[i] - psi[i + 1])))
        else:
            out.append(float("nan"))
    return out


def _seep02_section_flow(mesh, solution, x, split=None):
    """The discharge crossing the vertical line at ``x``, and the share of it below
    ``split``.

    Darcy's own horizontal velocity integrated up the line, rather than a difference
    of stream-function values: the stream function is a companion solve whose
    additive constant is not pinned on an unconfined problem, while this integral
    reproduces the reported total discharge to within the discretization.
    """
    import numpy as np

    trapezoid = getattr(np, "trapezoid", np.trapz)
    ivx = _seep02_interp(mesh, np.asarray(solution["velocity"])[:, 0])
    ys = np.linspace(1e-3, _seep02_top(x) - 1e-3, 20001)
    vx = np.ma.filled(ivx(np.full(len(ys), x), ys), 0.0)
    total = float(trapezoid(vx, ys))
    if split is None:
        return total, None
    below = ys <= split
    return total, float(trapezoid(vx[below], ys[below]))


def _seep02_unsaturated_flow(mesh, solution, x=SEEP02_DOWNSTREAM):
    """The share of the flow crossing ``x`` that travels above the phreatic
    surface — the answer to what the unsaturated zone is carrying, which is the
    quantity the choice of ``kr`` model acts on."""
    import numpy as np

    trapezoid = getattr(np, "trapezoid", np.trapz)
    ivx = _seep02_interp(mesh, np.asarray(solution["velocity"])[:, 0])
    top = _seep02_top(x)
    y_phreatic = _seep02_phreatic(mesh, solution, (x,))[0]
    ys = np.linspace(1e-3, top - 1e-3, 20001)
    vx = np.ma.filled(ivx(np.full(len(ys), x), ys), 0.0)
    dry = ys > y_phreatic
    total = float(trapezoid(vx, ys))
    above = float(trapezoid(vx[dry], ys[dry]))
    return total, above, 100.0 * above / total, y_phreatic


def _seep02_outline():
    """The dam's own outline and its two internal zone boundaries, for the figures
    that draw a phreatic surface on the section rather than through the plotting
    module."""
    shell = [(200, 100), (320, 160), (360, 180), (380, 180), (550, 100)]
    core = [(320, 100), (360, 165), (380, 165), (420, 100)]
    key = [(320, 100), (360, 60), (380, 60), (420, 100)]
    ground = [(0, 100), (750, 100)]
    base = [(0, 0), (750, 0), (750, 100)]
    return shell, core, key, ground, base


def _seep02_section_axes(ax):
    """One frame for every figure that draws the section by hand: the dam, the core
    and its key, the foundation, at equal aspect and one scale."""
    shell, core, key, ground, base = _seep02_outline()
    for pts, color, lw in ((ground, "#8a939c", 0.9), (base, "#8a939c", 0.9),
                           (shell, "#4a5560", 1.6), (core, "#7a5a3a", 1.3),
                           (key, "#7a5a3a", 1.3)):
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, lw=lw)
    ax.set_xlim(0, 750)
    ax.set_ylim(0, 195)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (ft)")
    ax.set_ylabel("elevation (ft)")


def _seep02_figsize(mesh):
    """The sample figures' own sizing rule — fixed width, height from the meshed
    domain's aspect — so a SEEP-2 flow net and a samples-page flow net of the same
    dam come out at the same scale."""
    import numpy as np

    nodes = np.asarray(mesh["nodes"])
    width, height = np.ptp(nodes[:, 0]), np.ptp(nodes[:, 1])
    return (11.0, max(2.6, 11.0 * 0.80 * (height / width) + 2.05))


def _seep02_stack(name, panels, dpi=200):
    """Several already-rendered figures stacked into one image.

    The base-material comparison has to show three flow nets drawn by
    ``plot_seep_solution`` itself — the point being what that function does with the
    argument — and it draws one figure at a time. So each panel is rendered on its
    own and the three are laid up here, which keeps every panel the tool's own
    output rather than a redrawing of it.
    """
    import numpy as np

    images = [plt.imread(p) for p in panels]
    heights = [im.shape[0] / float(im.shape[1]) for im in images]
    fig, axes = plt.subplots(len(images), 1, figsize=(11.0, 11.0 * sum(heights)))
    for ax, im in zip(np.atleast_1d(axes), images):
        ax.imshow(im)
        ax.set_axis_off()
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0.02)
    out = os.path.join(OUT_DIR, name)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    # The panels were scratch: the page carries the stack, so leaving three
    # near-duplicates of it in the image directory would put figures nothing links
    # to under version control.
    for path in panels:
        os.remove(path)
    print("-> %s  (%d panels)" % (name, len(images)))
    return out


def _seep02_kr_curves(materials, vg, gard):
    """Each material's three relative-conductivity curves on its own axes.

    Log kr against suction on a log axis, because the models differ by decades and
    across decades: the linear front is a straight line on linear axes and a cliff
    on these, which is the shape that matters to a solver. One panel per material
    rather than one crowded panel, because the three soils' curves overlap.
    """
    import numpy as np

    from xslope.seep import kr_frontal_vec, kr_gardner_vec, kr_vg_vec

    suction = np.logspace(-2, 2, 500)
    psi = -suction
    fig, axes = plt.subplots(1, len(materials), figsize=(13, 4.2), sharey=True)
    for ax, mat, (a, n), (ga, gn) in zip(axes, materials, vg, gard):
        ax.loglog(suction, kr_frontal_vec(psi, mat["kr0"], mat["h0"]),
                  color="#1f6fb4", lw=2.0,
                  label="lf  kr₀ = %g, h₀ = %g ft" % (mat["kr0"], mat["h0"]))
        ax.loglog(suction, kr_vg_vec(psi, a, n, 1e-4), color="#c1663a", lw=1.6,
                  label="vg  a = %g, n = %g" % (a, n))
        ax.loglog(suction, kr_gardner_vec(psi, ga, gn, 1e-4), color="#3f8f5a",
                  lw=1.6, ls="--", label="gard  a = %g, n = %g" % (ga, gn))
        ax.set_xlabel("suction −ψ (ft)")
        ax.set_title(mat["name"])
        ax.grid(True, which="both", color="#e8ebee", lw=0.5)
        ax.legend(loc="lower left", frameon=False, fontsize=8)
    axes[0].set_ylabel("relative conductivity $k_r$")
    axes[0].set_ylim(5e-5, 2.0)
    fig.suptitle("The three unsaturated models on this dam's three soils")
    fig.tight_layout()


def _seep02_anchors(sdta, sol):
    """Where the phreatic surface touches the dam: the entry point on the
    upstream face and the exit point on the seepage face.

    The entry is a boundary fact — the highest fixed-head node standing at the
    reservoir head, which is where the water line meets the face — and the exit
    is the highest exit-face node still at positive pressure. Both are read
    from the solve so a run that moves them (the core sweep moves the exit)
    carries its own. The drawn curves are anchored on these so they touch the
    dam at both ends, as the solution plot's own phreatic line does.
    """
    import numpy as np

    bct = np.asarray(sdta["bc_type"])
    nds = np.asarray(sdta["nodes"])
    p = np.asarray(sol["u"]) / 62.4
    head = p + nds[:, 1]
    fixed = np.where(bct == 1)[0]
    res = fixed[head[fixed] >= head[fixed].max() - 1e-6]
    entry = (float(nds[res[np.argmax(nds[res, 1])], 0]),
             float(nds[res[np.argmax(nds[res, 1])], 1]))
    face = np.where(bct == 2)[0]
    wet = face[p[face] > -1e-6]
    exit_pt = (float(nds[wet[np.argmax(nds[wet, 1])], 0]),
               float(nds[wet[np.argmax(nds[wet, 1])], 1])) if len(wet) else None
    return entry, exit_pt


def _seep02_anchored(ys, entry, exit_pt, stations=SEEP02_STATIONS):
    """The drawable polyline: entry point, the finite stations, exit point."""
    import numpy as np

    pts = [entry] + [(x, y) for x, y in zip(stations, ys)
                     if np.isfinite(y)] + ([exit_pt] if exit_pt else [])
    return [x for x, _ in pts], [y for _, y in pts]


def _seep02_phreatic_figure(series, stations=SEEP02_STATIONS):
    """The phreatic surfaces of the three unsaturated models on the section. At
    section scale the three surfaces are one line, and that IS the result; the
    page quotes the per-station maxima rather than charting them."""
    fig, ax = plt.subplots(figsize=(11, 5.2))
    _seep02_section_axes(ax)
    colors = {"lf": "#1f6fb4", "vg": "#c1663a", "gard": "#3f8f5a"}
    styles = {"lf": "-", "vg": "--", "gard": ":"}
    for label, ys, entry, exit_pt in series:
        xs_a, ys_a = _seep02_anchored(ys, entry, exit_pt, stations)
        ax.plot(xs_a, ys_a, styles[label], color=colors[label], lw=2.0,
                label="%s" % label)
    ax.legend(loc="upper right", frameon=False)
    ax.set_title("The phreatic surface under each unsaturated model")
    fig.tight_layout()


def _seep02_kr0_figure(sweep, references):
    """The discharge against the linear front's floor, with the two other models'
    answers drawn across it.

    The floor is what the page's explanation names, so the figure is the test of
    it: if the explanation is right the swept series lands on the reference lines
    as the floor reaches theirs, and if it is wrong it does not.
    """
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.semilogx([k for k, _ in sweep], [q for _, q in sweep], "o-",
                color="#1f6fb4", lw=1.4, ms=6, label="linear front, h₀ = −1 ft")
    # The two reference lines land within 0.0012 of each other, so their labels are
    # stacked rather than placed on the lines they belong to.
    for i, ((label, q), color) in enumerate(zip(references, ("#c1663a", "#3f8f5a"))):
        ax.axhline(q, color=color, lw=1.2, ls="--")
        ax.annotate("%s: q = %.4f" % (label, q), (1.05e-4, q),
                    xytext=(0, 8 + 12 * i), textcoords="offset points",
                    color=color, fontsize=9)
    ax.set_xlabel("linear-front floor $kr_0$")
    ax.set_ylabel("total discharge q (ft³/day per ft)")
    ax.set_title("Discharge against the floor of the relative-conductivity curve")
    ax.grid(True, which="both", color="#e8ebee", lw=0.6)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()


def _seep02_convergence_figure(histories, tol_scaled, closure_tol=1e-3):
    """Head change and flow closure against iteration, for every model the page
    runs, with the two thresholds they are tested against drawn on.

    Both panels on log axes and both thresholds drawn, because the lesson is that
    the two conditions fail in different places: the hard run's head change is
    inside its tolerance for hundreds of sweeps while its flow closure is not.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.6))
    colors = {"lf": "#1f6fb4", "vg": "#c1663a", "gard": "#3f8f5a",
              "lf, kr₀ = 1e−4, h₀ = −10 ft": "#8e44ad"}
    for label, rows in histories:
        color = colors.get(label, "#5a6b7a")
        ax1.semilogy([i for i, _, _ in rows], [r for _, r, _ in rows],
                     color=color, lw=1.4, label=label)
        ax2.semilogy([i for i, _, _ in rows], [c for _, _, c in rows],
                     color=color, lw=1.4, label=label)
    ax1.axhline(tol_scaled, color="#8a939c", lw=1.0, ls="--")
    ax1.annotate("head tolerance %.3g" % tol_scaled, (1, tol_scaled),
                 xytext=(4, 4), textcoords="offset points", color="#5b646f",
                 fontsize=9)
    ax2.axhline(closure_tol, color="#8a939c", lw=1.0, ls="--")
    ax2.annotate("closure tolerance %g" % closure_tol, (1, closure_tol),
                 xytext=(4, 4), textcoords="offset points", color="#5b646f",
                 fontsize=9)
    ax1.set_title("Head change per sweep")
    # The solver's head test is on the RELATIVE change — max nodal change over the
    # largest head in the field — so the axis and the threshold are dimensionless.
    ax1.set_ylabel("relative head change ‖Δh‖∞ / max|h|")
    ax2.set_title("Flow closure per sweep")
    ax2.set_ylabel("closure, fraction of inflow")
    for ax in (ax1, ax2):
        ax.set_xlabel("iteration")
        ax.grid(True, which="both", color="#e8ebee", lw=0.5)
    ax1.legend(loc="upper right", frameon=False, fontsize=8)
    fig.tight_layout()


def _seep02_core_figure(series, stations=SEEP02_STATIONS):
    """The phreatic surface and the seepage-face exit point at four core
    conductivities, on one section.

    One figure rather than four, because what the reader needs to see is the
    surface climbing the downstream slope as the core stops being a cutoff, and
    the exit point climbing with it.
    """
    fig, ax = plt.subplots(figsize=(11, 4.6))
    _seep02_section_axes(ax)
    colors = ("#1f6fb4", "#3f8f5a", "#c9a227", "#c1663a")
    for (k, ys, entry, exit_pt, _q, _n), color in zip(series, colors):
        xs_a, ys_a = _seep02_anchored(ys, entry, exit_pt, stations)
        ax.plot(xs_a, ys_a, "-", color=color, lw=1.8,
                label="core k = %g ft/day" % k)
        if exit_pt is not None:
            ax.plot([exit_pt[0]], [exit_pt[1]], "o", color=color, ms=7,
                    markeredgecolor="white", markeredgewidth=0.8)
    ax.legend(loc="upper right", frameon=False)
    ax.set_title("The phreatic surface and its exit point as the core's "
                 "conductivity rises")
    fig.tight_layout()


def seep02_plots():
    """The zoned dam: the model, the flow net, and the four studies the page runs.

    Printed rather than tabulated in a figure: the base run and where its flow
    goes, the flow-net channel count each base material asks for, the three
    unsaturated models against each other, the floor sweep that explains the
    difference between them, the tolerance ladder, the run that does not converge,
    and the core sweep that makes the seepage face grow. Every number the page
    quotes is on one of these lines.
    """
    import numpy as np

    from xslope.plot_seep import (flownet_base_material, plot_seep_data,
                                  plot_seep_solution)

    sd = load_slope_data(SEEP02)
    capture("seep02_inputs.png", plot_inputs, sd, mode="seep",
            title="Seepage Model Inputs", frame="content", show_mesh=False)

    # The mesh that travels with the file, before the page builds its own. It is
    # quadratic, because the same workbook is also run for stability, and the page
    # says what it gives so a reader who presses Run before meshing knows why their
    # number is not the catalogued one.
    companion = sd.get("mesh")
    if companion is not None:
        _, sol_c, log_c = _seep02_solve(sd, companion)
        st_c = _seep02_log_stats(log_c)
        print("   companion   %d nodes · %d elements · element type %s · q %.4f · "
              "%s iterations"
              % (len(companion["nodes"]), len(companion["elements"]),
                 sorted(set(int(t) for t in companion["element_types"])),
                 sol_c["flowrate"], st_c["iterations"]))

    mesh = _seep02_mesh(sd)
    figsize = _seep02_figsize(mesh)
    seep_data, solution, log = _seep02_solve(sd, mesh)
    stats = _seep02_log_stats(log)
    capture("seep02_mesh.png", plot_seep_data, seep_data, figsize=figsize,
            show_bc=True)
    capture("seep02_solution.png", plot_seep_solution, seep_data, solution,
            figsize=figsize, levels=SEEP02_LEVELS, base_mat=SEEP02_BASE_MAT,
            fill_contours=True, mesh=False)
    capture("seep02_pressure.png", plot_seep_solution, seep_data, solution,
            figsize=figsize, levels=SEEP02_LEVELS, variable="u", mesh=False,
            flowlines=False)

    nodes = np.asarray(mesh["nodes"])
    head = np.asarray(solution["head"])
    psi = np.asarray(solution["u"]) / 62.4
    print("   mesh        %d nodes · %d elements · %s at width/%d = %.4g ft"
          % (len(nodes), len(mesh["elements"]), SEEP02_ELEMENT, SEEP02_DIVISIONS,
             750.0 / SEEP02_DIVISIONS))
    print("   base run    q %.4f · head %.3f to %.3f · u %.1f to %.1f psf"
          % (solution["flowrate"], head.min(), head.max(),
             np.min(solution["u"]), np.max(solution["u"])))
    print("   iteration   %s" % stats)
    print("   suction     min ψ %.2f ft at (%.1f, %.1f)"
          % (psi.min(), nodes[np.argmin(psi), 0], nodes[np.argmin(psi), 1]))

    # Where the flow goes. A cutoff core that reaches 40 ft into the foundation is
    # not a cutoff, and the centerline section says by how much.
    total, below = _seep02_section_flow(mesh, solution, SEEP02_CENTERLINE,
                                        SEEP02_KEY_TOE)
    print("   centerline  x %g: q %.4f · below elevation %g %.4f (%.1f%%)"
          % (SEEP02_CENTERLINE, total, SEEP02_KEY_TOE, below, 100.0 * below / total))
    for x in (400.0, 450.0, SEEP02_DOWNSTREAM):
        tot, above, share, y_ph = _seep02_unsaturated_flow(mesh, solution, x)
        # How far the section integral misses the reported total discharge is the
        # resolution of the share beside it: a share read off an integral that is
        # itself several percent out cannot be quoted to the tenth.
        print("   unsaturated x %g: q %.4f (%+.1f%% of the reported total) · above "
              "the phreatic surface at %.2f ft %.4f (%.1f%%)"
              % (x, tot, 100.0 * (tot - solution["flowrate"]) / solution["flowrate"],
                 y_ph, above, share))

    # The head the core drops, read either side of it at three elevations.
    ihead = _seep02_interp(mesh, head)
    for elev in (110.0, 130.0, 150.0):
        up, down = float(ihead(315.0, elev)), float(ihead(425.0, elev))
        print("   core drop   elevation %g: %.2f ft upstream, %.2f ft downstream, "
              "drop %.2f of 60" % (elev, up, down, up - down))

    # The exit face: how much of it the water actually uses.
    bc_type = np.asarray(seep_data["bc_type"])
    face = np.where(bc_type == 2)[0]
    active = face[psi[face] > -1e-6]
    top = face[np.argmax(nodes[face, 1])]
    print("   exit face   %d nodes from (%.1f, %.2f) to (%.1f, %.2f) · %d active · "
          "highest wet (%.1f, %.2f)"
          % (len(face), nodes[face, 0].max(), nodes[face, 1].min(),
             nodes[top, 0], nodes[top, 1], len(active),
             nodes[active[np.argmax(nodes[active, 1])], 0],
             nodes[active[np.argmax(nodes[active, 1])], 1]))

    # ---- the flow net's base material ------------------------------------- #
    # Three nets of the same solution, one per zone, because the argument's whole
    # effect is the channel count and the count is what the three panels differ by.
    print("   -- flow-net base material, %d contour levels" % SEEP02_LEVELS)
    drops = SEEP02_LEVELS - 1
    hdrop = float(head.max() - head.min())
    panels = []
    for i, mat in enumerate(sd["materials"], 1):
        k = math.sqrt(float(mat["k1"]) * float(mat["k2"]))
        channels = solution["flowrate"] * drops / (k * hdrop)
        print("   base_mat %d %-11s k %-8g Nf = q·Nd/(k·Δh) = %8.2f · %d φ contour "
              "levels requested"
              % (i, mat["name"], k, channels, max(round(channels) + 1, 2)))
        panels.append(capture("seep02_base_mat_%d.png" % i, plot_seep_solution,
                              seep_data, solution, figsize=figsize,
                              levels=SEEP02_LEVELS, base_mat=i,
                              fill_contours=False, mesh=False))
    print("   flownet_base_material picks %d"
          % flownet_base_material(seep_data, solution, levels=SEEP02_LEVELS))
    _seep02_stack("seep02_base_mat.png", panels)

    # ---- the three unsaturated models ------------------------------------- #
    from xslope.seep import kr_gardner_vec, kr_vg_vec

    suction = np.logspace(-2, 2, 600)
    print("   -- Gardner fitted to van Genuchten, per material")
    for mat, (a, n), (ga, gn) in zip(sd["materials"], SEEP02_VG, SEEP02_GARD):
        misfit = float(np.sqrt(np.mean(
            (np.log10(kr_gardner_vec(-suction, ga, gn, 1e-4))
             - np.log10(kr_vg_vec(-suction, a, n, 1e-4))) ** 2)))
        print("   %-11s vg a %-6g n %-5g · gard a %-6g n %-5g · rms log10 kr %.3f"
              % (mat["name"], a, n, ga, gn, misfit))
    capture("seep02_kr_models.png", _seep02_kr_curves, sd["materials"],
            SEEP02_VG, SEEP02_GARD)

    # The two places the curves part company, as numbers rather than as a picture:
    # near the surface, where the linear front is still on its straight line, and
    # deep, where it has stopped falling and the other two have not.
    from xslope.seep import kr_frontal_vec

    print("   -- kr at four suctions, shell parameters")
    probe = np.array([-0.5, -1.0, -5.0, -20.0])
    print("   %-6s %s" % ("ψ (ft)", " ".join("%9.1f" % p for p in probe)))
    for label, values in (
            ("lf", kr_frontal_vec(probe, sd["materials"][0]["kr0"],
                                  sd["materials"][0]["h0"])),
            ("vg", kr_vg_vec(probe, *SEEP02_VG[0], 1e-4)),
            ("gard", kr_gardner_vec(probe, *SEEP02_GARD[0], 1e-4))):
        print("   %-6s %s" % (label, " ".join("%9.5f" % v for v in values)))

    models = [
        ("lf", sd),
        ("vg", _seep02_materials(sd, [dict(unsat="vg", vg_a=a, vg_n=n)
                                      for a, n in SEEP02_VG])),
        ("gard", _seep02_materials(sd, [dict(unsat="gard", vg_a=a, vg_n=n)
                                        for a, n in SEEP02_GARD])),
    ]
    print("   -- the three unsaturated models on one mesh")
    phreatics, histories, refs = [], [], []
    for label, model in models:
        sdta, sol, mlog = _seep02_solve(model, mesh)
        st = _seep02_log_stats(mlog)
        ys = _seep02_phreatic(mesh, sol)
        _tot, _above, share, y_ph = _seep02_unsaturated_flow(mesh, sol)
        phreatics.append((label, ys, *_seep02_anchors(sdta, sol)))
        histories.append((label, _seep02_history(mlog)))
        if label != "lf":
            refs.append((label, sol["flowrate"]))
        # The relaxation factor comes off the LADDER at the sweep the run finished
        # on, not off the log: the log prints a sweep only every fifth one after the
        # sixth, so a run finishing on sweep 23 never prints sweeps 21-23 and the
        # scraped minimum reads 1.00 for a run that spent three sweeps at 0.5.
        print("   %-5s q %.4f · %3d iterations · relax at the last sweep %.2f "
              "(scraped from the log: %.2f) · %d sweeps past the full step · "
              "%d/%d exit face active · unsaturated share at x %g %.1f%% · "
              "phreatic there %.2f ft"
              % (label, sol["flowrate"], st["iterations"],
                 _seep02_relax_at(st["iterations"]), st["relax_min"],
                 max(0, st["iterations"] - SEEP02_RELAX_LADDER[0][0]),
                 st["active"][0], st["active"][1],
                 SEEP02_DOWNSTREAM, share, y_ph))
    for i, (a, ys_a, _e1, _x1) in enumerate(phreatics):
        for b, ys_b, _e2, _x2 in phreatics[i + 1:]:
            print("   %-4s vs %-4s phreatic surface differs by at most %.2f ft "
                  "over the %d stations"
                  % (a, b, np.nanmax(np.abs(np.asarray(ys_a)
                                            - np.asarray(ys_b))),
                     len(SEEP02_STATIONS)))
    print("   stations    %s" % " ".join("%g" % x for x in SEEP02_STATIONS))
    for label, ys, _entry, _exit in phreatics:
        print("   %-5s       %s" % (label, " ".join("%6.2f" % y for y in ys)))
    capture("seep02_phreatic_models.png", _seep02_phreatic_figure, phreatics)

    # ---- what the difference between them actually is --------------------- #
    print("   -- the linear front's floor, swept toward the other models'")
    sweep = []
    for kr0 in sorted(SEEP02_KR0_SWEEP, reverse=True):
        model = _seep02_materials(sd, [dict(unsat="lf", kr0=kr0, h0=-1.0)] * 3)
        _, sol, slog = _seep02_solve(model, mesh)
        st = _seep02_log_stats(slog)
        sweep.append((kr0, sol["flowrate"]))
        print("   kr0 %-8g q %.4f · %s iterations · relax at the last sweep %.2f "
              "(scraped from the log: %.2f) · %s"
              % (kr0, sol["flowrate"], st["iterations"],
                 _seep02_relax_at(st["iterations"]), st["relax_min"],
                 "converged" if st["converged"] else "DID NOT CONVERGE"))
    capture("seep02_kr0_sweep.png", _seep02_kr0_figure, sweep, refs)

    # The shape of the curve above the floor, at one floor — the other half of the
    # two-parameter linear front, so the page can say which of them matters.
    print("   -- the linear front's reference suction, at kr0 = 0.01")
    for h0 in (-0.5, -1.0, -2.0, -5.0, -10.0):
        model = _seep02_materials(sd, [dict(unsat="lf", kr0=0.01, h0=h0)] * 3)
        _, sol, slog = _seep02_solve(model, mesh)
        st = _seep02_log_stats(slog)
        print("   h0 %-6g q %.4f · %s iterations" % (h0, sol["flowrate"],
                                                     st["iterations"]))

    # ---- convergence ------------------------------------------------------ #
    print("   -- the solve tolerance, and what the converged discharge does with it")
    for tol in SEEP02_TOLS:
        _, sol, tlog = _seep02_solve(sd, mesh, tol=tol)
        st = _seep02_log_stats(tlog)
        print("   tol %-8g scaled to %.4g · q %.6f · %s iterations · relax at the "
              "last sweep %.2f" % (tol, st["tol_scaled"], sol["flowrate"],
                                   st["iterations"],
                                   _seep02_relax_at(st["iterations"])))

    print("   -- the run that does not converge inside the default ceiling")
    hard = _seep02_materials(sd, [dict(unsat="lf", **SEEP02_HARD)] * 3)
    _, sol_hard, hlog = _seep02_solve(hard, mesh)
    tail = [ln for ln in hlog.strip().splitlines()
            if ln.startswith("Iteration 400") or ln.startswith("Warning")
            or ln.startswith("WARNING")]
    print("   default cap q %.4f · converged %s"
          % (sol_hard["flowrate"], sol_hard.get("converged")))
    # WHICH of the three conditions failed is the lesson, so the sweep the head
    # test started passing on is measured rather than eyeballed off the figure.
    # The log prints every fifth sweep after the sixth, so this is the first
    # REPORTED sweep inside the head tolerance rather than the first one — which is
    # the right granularity for a page that quotes the log.
    inside = [i for i, r, _ in _seep02_history(hlog) if r < stats["tol_scaled"]]
    print("   head change inside the %.4g tolerance by reported sweep %s, and "
          "stayed there" % (stats["tol_scaled"], inside[0] if inside else "never"))
    # The band the flow closure orbits in, rather than an eyeballed pair of round
    # numbers off the figure, and the sweep the exit-face set stopped moving on.
    closures = [c for _, _, c in _seep02_history(hlog)]
    print("   closure     oscillates between %.4g and %.4g over the %d reported "
          "sweeps" % (min(closures), max(closures), len(closures)))
    for label, alog in (("base run", log), ("hard run", hlog)):
        act = _seep02_active_history(alog)
        settled = act[-1][1]
        first = next(i for i, (_it, n) in enumerate(act)
                     if all(m == settled for _jt, m in act[i:]))
        print("   %s   exit face %s ... settled to %d from sweep %d"
              % (label, " ".join("%d" % n for _it, n in act[:6]), settled,
                 act[first][0]))
    for line in tail:
        print("     %s" % line)
    histories.append(("lf, kr₀ = 1e−4, h₀ = −10 ft", _seep02_history(hlog)))
    _, sol_long, llog = _seep02_solve(hard, mesh, max_iter=SEEP02_HARD_MAX_ITER)
    st_long = _seep02_log_stats(llog)
    print("   cap %d     q %.4f · %s iterations · relax down to %.2f · converged %s"
          % (SEEP02_HARD_MAX_ITER, sol_long["flowrate"], st_long["iterations"],
             st_long["relax_min"], sol_long.get("converged")))
    print("   the truncated answer is %.2f%% from the converged one"
          % (100.0 * abs(sol_hard["flowrate"] - sol_long["flowrate"])
             / sol_long["flowrate"]))
    capture("seep02_convergence.png", _seep02_convergence_figure, histories,
            stats["tol_scaled"])

    # ---- the seepage face, given something to do -------------------------- #
    print("   -- the core's conductivity, and the seepage face it leaves")
    core_series = []
    for k in SEEP02_CORE_KS:
        model = _seep02_materials(sd, [{}, dict(k1=k, k2=k), {}])
        sdta, sol, clog = _seep02_solve(model, mesh)
        st = _seep02_log_stats(clog)
        bct = np.asarray(sdta["bc_type"])
        p = np.asarray(sol["u"]) / 62.4
        face = np.where(bct == 2)[0]
        wet = face[p[face] > -1e-6]
        entry, exit_pt = _seep02_anchors(sdta, sol)
        core_series.append((k, _seep02_phreatic(mesh, sol), entry, exit_pt,
                            sol["flowrate"], len(wet)))
        print("   core k %-8g q %8.4f · %2d/%d exit-face nodes wet · highest wet "
              "node %s · %s iterations"
              % (k, sol["flowrate"], len(wet), len(face), exit_pt, st["iterations"]))
    capture("seep02_core_sweep.png", _seep02_core_figure, core_series)


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
    "lem07_plots": lem07_plots,
    "lem08_sheets": lem08_sheets,
    "lem08_plots": lem08_plots,
    "lem08_lengths": lem08_lengths,
    "lem09_sheets": lem09_sheets,
    "lem09_plots": lem09_plots,
    "lem10_plots": lem10_plots,
    "lem11_plots": lem11_plots,
    "lem12_plots": lem12_plots,
    "seep01_sheets": seep01_sheets,
    "seep01_plots": seep01_plots,
    "seep02_plots": seep02_plots,
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
