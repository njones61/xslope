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
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402

from xslope.fileio import load_slope_data                            # noqa: E402
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

#: Where the page pulls the entry point to, to show what a steep toe wedge does.
#: 1 ft of horizontal run against 5 ft of drop is a 78.7 degree leading segment —
#: still sliceable, still solvable, and the answer it returns is nonsense.
LEM05_SLIVER_X = -1.0


def _lem05_solve(model, non_circ, method="spencer"):
    """One surface, one method, no search: what this page runs everywhere."""
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
    """The states LEM-5 compares against — every one a single surface, solved once.

    No search anywhere in this group. The page teaches the surface a reader
    enters, so each figure is that surface (or a deliberately damaged version of
    it) run through Spencer's method at 40 slices, which is the method and the
    slice count the prose quotes.

    The arc: the model as delivered, the generator's own proposal beside it, the
    best circle this section admits, and the same surface with its entry point
    pulled in until the leading segment is a sliver.
    """
    sd = load_slope_data(LEM05)

    capture("lem05_inputs.png", plot_inputs, sd, title="Slope Geometry and Inputs")

    slices, surface, result = _lem05_solve(sd, sd["non_circ"])
    capture("lem05_solution.png", plot_solution, sd, slices, surface, result)

    # The generator's proposal, solved the same way — the audit the Studio path
    # asks the reader to make against the surface they entered by hand.
    from xslope.generators import generate_noncircular_surface
    gen = generate_noncircular_surface(sd, report=True)
    if not gen["surface"]:
        raise SystemExit("LEM-5: the weak-zone generator built nothing — %s"
                         % (gen["reason"],))
    gs, gsurf, gres = _lem05_solve(sd, gen["surface"])
    capture("lem05_solution_generated.png", plot_solution, sd, gs, gsurf, gres)

    # What a circle gets on the same section. The file defines no circles — it is
    # a non-circular model — so the question is asked with the starting circles
    # the geometry itself proposes, refined by the ordinary circular search.
    from xslope.generators import generate_starting_circles
    circ = copy.deepcopy(sd)
    circ["circles"] = generate_starting_circles(circ)
    circ["circular"] = True
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, _, _ = circular_search(circ, "spencer",
                                            num_slices=LEM05_SLICES,
                                            diagnostic=False)
    best = fs_cache[0]
    capture("lem05_solution_circle.png", plot_solution, circ, best["slices"],
            best["failure_surface"], best["solver_result"])

    # The hazard: the same four points with the entry pulled in to x = -1, which
    # stands the leading segment up at 78.7 degrees. Drawn rather than described,
    # because the tell is the shape of the leading slices.
    sliver = copy.deepcopy(sd["non_circ"])
    sliver[0]["X"] = LEM05_SLIVER_X
    ss, ssurf, sres = _lem05_solve(sd, sliver)
    capture("lem05_sliver.png", plot_solution, sd, ss, ssurf, sres)

    print("   as entered %.4f · generated %.4f · best circle %.4f (depth %.3f) · "
          "entry at x = %g %.4f"
          % (result["FS"], gres["FS"], best["FS"], best["Depth"],
             LEM05_SLIVER_X, sres["FS"]))


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
