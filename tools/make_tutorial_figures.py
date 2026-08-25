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
from xslope.plot import (declared_unit_labels,                       # noqa: E402
                         plot_inputs, plot_solution,
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
    # The problem figure is the hand drawing crest_surcharge.fodg (private repo).

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



def lem06_plots():
    """The states LEM-6 reads, in the order the page walks them.

    The search runs first, on the model as delivered, and settles well above the
    dipping base. The two figures after it are what makes the base do something:
    a circle pushed below it and truncated against it (the composite option), and
    the same section with the foundation softened below the fill, where the
    critical mechanism drops onto the base itself.
    """
    sd = load_slope_data(LEM06)

    # The problem figure is the hand drawing sloping_bottom.fodg (private repo).
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
    column-name row (row 2) to column Q, the last column a limit equilibrium run
    reads — the six lines plus two blank rows under them, so the block reads as
    unfilled below the last line. Columns R:T (Tres, E, Area) are read only by
    the finite element engine, and the Type preset's lookup table beyond them is
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
    render("lem08_sheet_reinforce.png", LEM08, "reinforce", rows=(2, 10), cols="A:Q")
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
              "depth %.3f)  each crossing (x, ft to nearer end, T) %s"
              % (length, crit["FS"], crit["slices"]["p"].sum(), crossed,
                 crit["Xo"], crit["Yo"], crit["Depth"],
                 " ".join("(%.2f, %.1f, %.0f)" % (x, d, t)
                          for x, d, t in tension)))
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
    # LEM problem: stop before the FEM-only pile columns (E, I, Area, Head, Tip).
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
    # The problem figure is the hand drawing tieback_wall.fodg (private repo).

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


def _stack_panels(name, panels, dpi=200):
    """Several already-rendered figures stacked into one image.

    SEEP-2's base-material comparison has to show three flow nets drawn by
    ``plot_seep_solution`` itself — the point being what that function does with the
    argument — and it draws one figure at a time. So each panel is rendered on its
    own and the three are laid up here, which keeps every panel the tool's own
    output rather than a redrawing of it. SEEP-3's frame series is the same shape of
    problem, one instant per panel, and uses the same helper.
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
    _stack_panels("seep02_base_mat.png", panels)

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


# --------------------------------------------------------------------------- #
# SEEP-3 — Transient Seepage: Reservoir Drawdown Through a Cored Earth Dam
#
# The first page on which the answer is a sequence rather than a field. The model is
# the cored earth dam of docs/seep/samples.md #8, in the tutorial's own sidecar-free
# copy: a granular shell around a compacted clay core, a full pool at el 18 held
# against the upstream face, a tailwater at el 2, and a reservoir drawn down to the
# tailwater over 45 days.
#
# Every figure here is one of two runs. The steady solve at full pool is the initial
# condition — the state the dam is in before anything moves, and the state the
# transient march starts from. The march itself produces the frame series and the
# time histories, and it is run through the same ``run_transient_seepage`` call the
# Studio runner makes, on the mesh the page's Build Mesh step produces.
#
# The lesson the figures have to carry is the LAG: the pool falls in 45 days, the
# shell follows it within days, and the core does not. The frame series shows it as
# a shape and the history shows it as numbers, so both are drawn from the same
# march rather than from two.
# --------------------------------------------------------------------------- #
SEEP03 = os.path.join(REPO_ROOT,
                      "docs/tutorials/files/xslope_earth_dam_drawdown.xlsx")
#: The mesh, in the Build mesh dialog's own controls: tri3 (the seepage default),
#: auto-sized at 64 divisions across the 110 m section — the 1.71875 m target the
#: sample page's own transient figures are computed at.
SEEP03_DIVISIONS = 64
SEEP03_ELEMENT = "tri3"
#: Contour count. Twelve, the count the sample's frame panels use: a transient frame
#: is not a flow net (see below), so the level count is a readability choice rather
#: than the flow-channel arithmetic SEEP-1 and SEEP-2 do with it.
SEEP03_LEVELS = 12
#: The flow net's base material for the STEADY figure, 1-based into the mat sheet:
#: 2 is the core, which is what the steady sibling sample (Problem 3) is drawn with.
SEEP03_BASE_MAT = 2
#: Per-panel figure size for the stacked frame series, sized to this dam's aspect —
#: the same panel size the sample's own series figure uses.
SEEP03_PANEL_SIZE = (8.6, 2.55)
#: The instants the frame series is drawn at, and what each one is. Full pool is the
#: initial condition; 15 is mid-drawdown, with the pool already a few metres below the
#: interior surface; 47 is the end of the drawdown, where the lag is largest; 120 is
#: the recovery, where the interior mound is draining toward the new steady state.
SEEP03_FRAMES = ((0.0, "full pool — the initial condition"),
                 (15.0, "mid drawdown"),
                 (47.0, "end of drawdown — the largest lag"),
                 (120.0, "recovery toward the new steady state"))
#: Vertical stations the phreatic surface is read at, upstream toe to downstream toe.
#: 46 / 54.5 / 63 straddle the core (its base runs from x = 46 to x = 63), so the
#: three of them measure what the core is holding against what the shell has let go.
SEEP03_STATIONS = (10.0, 20.0, 30.0, 40.0, 46.0, 50.0, 54.5, 59.0, 63.0, 70.0,
                   80.0, 90.0, 100.0)
#: The three nodes the head history is read at, as target points the mesh is
#: searched near. One inside the clay core and one in the upstream shell at the SAME
#: elevation, so the two traces are a like-for-like comparison and their crossing is
#: the lag itself; the third is in the downstream shell, which the drawdown reaches
#: last. The producer prints the node ids and the coordinates it actually found.
SEEP03_NODES = (("clay core", (54.5, 9.0), "#b5460f"),
                ("upstream shell", (30.0, 9.0), "#1f6fb4"),
                ("downstream shell", (75.0, 6.0), "#3f8f5a"))


def _seep03_mesh(model, divisions=SEEP03_DIVISIONS, element_type=SEEP03_ELEMENT):
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)

    xs = [x for x, _ in model["ground_surface"].coords]
    size = (max(xs) - min(xs)) / divisions
    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(
            get_material_polygons(model), size, element_type,
            size_regions=extract_size_regions(model))


def _seep03_steady_model(model):
    """``model`` with the reservoir boundary pinned at its own t = 0 level, in memory
    only — the initial condition as a steady problem.

    The transient solver's initial condition IS a steady solve at the t = 0 boundary
    configuration, and at t = 0 the whole submerged face stands below full pool, so a
    submerged-only reservoir at el 18 and a fixed head of 18 are the same boundary.
    Pinning it here is what lets the page run the initial condition as an ordinary
    steady analysis and read a flowrate off it, which a series-bound value has none of.
    """
    bc = {k: (list(v) if isinstance(v, list) else v)
          for k, v in model["seepage_bc"].items()}
    heads = [dict(h) for h in bc["specified_heads"]]
    heads[0] = dict(heads[0], kind="head",
                    head=float(model["tseep"]["series"]["pool"][0]))
    bc["specified_heads"] = heads
    return dict(model, seepage_bc=bc, tseep=None)


def _seep03_solve(model, mesh, tol=1e-4, max_iter=400):
    """One steady seepage solve, returning ``(seep_data, solution, log)`` — the same
    shape SEEP-2's solver helper returns, and for the same reason: on an unconfined
    problem the iteration history is half the result."""
    from xslope.seep import build_seep_data, run_seepage_analysis

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        seep_data = build_seep_data(mesh, model)
        solution = run_seepage_analysis(seep_data, tol=tol, max_iter=max_iter)
    return seep_data, solution, buf.getvalue()


def _seep03_march(model, mesh):
    """The transient march, returning ``(seep_data, solution, frames, log)``.

    ``run_transient_seepage`` is the call Studio's ``SeepRunner`` makes, and the
    per-frame plottable dicts are rebuilt with the same ``_transient_frame_solution``
    the runner uses, so the frames this producer draws are the frames the reader's
    play bar shows.
    """
    from xslope.seep import (build_seep_data, build_tseep_data,
                             run_transient_seepage, _transient_frame_solution)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        seep_data = build_seep_data(mesh, model, seep_bc=1)
        solution = run_transient_seepage(seep_data, build_tseep_data(model),
                                         verbose=True)
    unconfined = bool(solution.get("unconfined", False))
    frames = [_transient_frame_solution(seep_data, fr["head"], fr["u"],
                                        fr.get("phi"), fr.get("inflow"),
                                        fr.get("outflow"), unconfined,
                                        time=fr["time"])
              for fr in solution["frames"]]
    return seep_data, solution, frames, buf.getvalue()


def _seep03_top(model, x):
    """Elevation of the top of the section at ``x``, off the model's own outer
    profile line — so a station's search column cannot run past the ground."""
    pts = list(model["ground_surface"].coords)
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return pts[-1][1]


def _seep03_phreatic(model, mesh, solution, stations=SEEP03_STATIONS):
    """The elevation of the phreatic surface at each station, read off a linear
    interpolator rather than off the nodes.

    Same construction as SEEP-2's: walk a vertical column and return the elevation
    where the pressure head crosses zero. A column that is saturated to the ground
    surface has no crossing inside the domain and returns NaN — on this dam that is
    the submerged upstream face at full pool, where the p = 0 level is the reservoir
    standing above the ground.
    """
    import matplotlib.tri as mtri
    import numpy as np

    nodes = np.asarray(mesh["nodes"])
    tris = [e[:3] for e in np.asarray(mesh["elements"])]
    triang = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tris)
    ipsi = mtri.LinearTriInterpolator(triang, np.asarray(solution["u"]) / 62.4)
    out = []
    for x in stations:
        ys = np.linspace(0.02, _seep03_top(model, x) - 0.02, 4000)
        psi = np.ma.filled(ipsi(np.full(len(ys), x), ys), np.nan)
        ok = ~np.isnan(psi)
        ys, psi = ys[ok], psi[ok]
        wet = np.where(psi >= 0.0)[0]
        if not len(wet) or wet[-1] + 1 >= len(psi):
            out.append(float("nan"))
            continue
        i = wet[-1]
        out.append(float(ys[i] + (ys[i + 1] - ys[i]) * psi[i]
                         / (psi[i] - psi[i + 1])))
    return out


def _seep03_pool(model, t):
    """The pool level at time ``t``, read off the model's own series with the tseep
    interpolation rule (linear between breakpoints, held beyond the ends)."""
    import numpy as np

    return float(np.interp(t, model["tseep"]["times"],
                           model["tseep"]["series"]["pool"]))


def _seep03_pick_nodes(mesh, seep_data):
    """The nearest mesh node to each target point, with the material it belongs to.

    The targets are written as points rather than as node numbers because the node
    numbering is the mesher's and changes with the element size; the page quotes the
    coordinates this returns, which are a property of the dam.
    """
    import numpy as np

    nodes = np.asarray(mesh["nodes"])
    elem_mat = np.asarray(seep_data["element_materials"])
    node_mat = {}
    for elem, mat in zip(np.asarray(seep_data["elements"]), elem_mat):
        for n in elem[:3]:
            node_mat.setdefault(int(n), set()).add(int(mat))
    picked = []
    for label, point, color in SEEP03_NODES:
        i = int(np.argmin(np.linalg.norm(nodes - np.asarray(point), axis=1)))
        picked.append((label, i, tuple(float(v) for v in nodes[i]), color,
                       sorted(node_mat.get(i, set()))))
    return picked


def _seep03_schedule(model):
    """The pool series the reservoir boundary follows, as the reader enters it —
    three breakpoints and a duration.

    A figure of its own rather than a panel under the section: the section is a
    sketch of the dam a reader meets before anything is built, while the schedule
    belongs with the transient dialog that carries these same three rows.
    """
    times = list(model["tseep"]["times"])
    pool = list(model["tseep"]["series"]["pool"])
    duration = float(model["tseep"]["duration"])

    fig, ax = plt.subplots(figsize=(8.6, 2.8))
    ax.axvspan(times[1], times[-1], color="#eef4f9", zorder=0)
    ax.plot(times + [duration], pool + [pool[-1]], color="#2b7bb0", lw=2.2,
            zorder=3)
    ax.plot(times, pool, "o", ms=6, mfc="white", mec="#2b7bb0", mew=1.8,
            zorder=4, label="schedule breakpoints")
    ax.annotate("held at full pool", xy=(times[1], pool[1]),
                xytext=(times[1] + 26.0, pool[1] + 1.6), fontsize=9,
                color="0.3", ha="left",
                arrowprops=dict(arrowstyle="-|>", color="0.45", lw=0.9))
    ax.annotate("drawn down over %g days" % (times[-1] - times[1]),
                xy=(0.5 * (times[1] + times[-1]), 0.5 * (pool[1] + pool[-1])),
                xytext=(times[-1] + 34.0, 0.66 * pool[0]), fontsize=9,
                color="#2b7bb0", ha="left",
                arrowprops=dict(arrowstyle="-|>", color="#2b7bb0", lw=0.9))
    ax.annotate("held for the rest of the %g-day run,\nwhile the dam relaxes"
                % duration, xy=(0.80 * duration, pool[-1]),
                xytext=(0.46 * duration, pool[-1] + 2.6), fontsize=9,
                color="0.3", ha="left",
                arrowprops=dict(arrowstyle="-|>", color="0.45", lw=0.9))
    ax.set_xlim(-8.0, duration + 8.0)
    ax.set_ylim(0.0, pool[0] + 5.0)
    ax.set_xlabel("time (day)")
    ax.set_ylabel("pool elevation (%s)" % declared_unit_labels(model)["length"])
    ax.grid(alpha=0.22)
    ax.legend(loc="upper right", frameon=False, fontsize=8.5)
    fig.tight_layout()


def _seep03_frame_panel(seep_data, frame, label, vmin, vmax, pool):
    """One saved instant, drawn through ``plot_seep_solution`` itself.

    Flow lines are off on every panel: a flow net requires divergence-free
    through-flow, and a transient state's flow comes out of released storage, so no
    stream function exists to draw equal-drop channels from. What each panel does
    show is the head field, the phreatic surface, and the instantaneous water levels
    — so the pool visibly falls through the series while the interior surface lags.
    The color scale is pinned across the series so the frames read as one story.
    """
    from xslope.plot_seep import plot_seep_solution

    with _hold_show():
        plot_seep_solution(seep_data, frame, figsize=SEEP03_PANEL_SIZE,
                           levels=SEEP03_LEVELS, phreatic=True, flowlines=False,
                           vectors=False, mesh=False, pad_frac=0.04,
                           cmap="Spectral_r", vmin=vmin, vmax=vmax,
                           show_title=False, show_legend=False,
                           show_bc_levels=True)
    _u = declared_unit_labels(seep_data)
    plt.gcf().axes[0].set_title("t = %g %s   ·   pool el %.1f %s   —   %s"
                                % (frame["time"], _u["time"], pool, _u["length"],
                                   label), fontsize=11, pad=5)


def _seep03_history_figure(model, times, series, picked):
    """Total head against time at the three nodes, with the pool schedule over them.

    One axes rather than three: the whole point is that the shell's trace and the
    core's trace CROSS — the shell starts above the core and ends below it — and a
    crossing cannot be read off panels side by side.
    """
    import numpy as np

    duration = float(model["tseep"]["duration"])
    t0, t1 = model["tseep"]["times"][1], model["tseep"]["times"][-1]
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.axvspan(t0, t1, color="#eef4f9", zorder=0)
    # The pool is the driver, not a fourth reading, so it is drawn in neutral gray:
    # the three colored traces are the dam's response to this one dashed line.
    ax.plot(times, [_seep03_pool(model, t) for t in times], color="#3f4a55",
            lw=2.0, ls=(0, (5, 3)), zorder=2, label="pool (the reservoir schedule)")
    for (label, _i, coords, color, _mats), values in zip(picked, series):
        ax.plot(times, values, "-o", color=color, lw=1.8, ms=4, zorder=3,
                label="%s  (%.1f, %.1f)" % (label, coords[0], coords[1]))

    # Where the core's trace passes ABOVE the shell's, interpolated between the two
    # saved frames that straddle it rather than read off the picture.
    core, shell = np.asarray(series[0]), np.asarray(series[1])
    gap = core - shell
    cross = np.where((gap[:-1] < 0) & (gap[1:] >= 0))[0]
    if len(cross):
        j = int(cross[0])
        tc = (times[j] + (times[j + 1] - times[j]) * (-gap[j])
              / (gap[j + 1] - gap[j]))
        ax.axvline(tc, color="0.6", lw=0.9, ls=(0, (2, 3)), zorder=1)
        ax.annotate("the core's head passes above the shell's\nat t ≈ %.0f day" % tc,
                    xy=(tc, float(np.interp(tc, times, core))),
                    xytext=(tc + 55.0, 0.80 * max(core.max(), shell.max())),
                    fontsize=9, color="0.25", ha="left", va="center",
                    arrowprops=dict(arrowstyle="-|>", color="0.45", lw=0.9))
    ax.annotate("drawdown", xy=(0.5 * (t0 + t1), 0.98), xycoords=("data",
                                                                  "axes fraction"),
                fontsize=8.5, color="#2b7bb0", ha="center", va="top")
    ax.set_xlim(-8.0, duration + 8.0)
    ax.set_xlabel("time (day)")
    ax.set_ylabel("total head (%s)" % declared_unit_labels(model)["length"])
    ax.grid(alpha=0.22)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_title("The core holds its head after the shell has let go", fontsize=11.5)
    fig.tight_layout()


def seep03_plots():
    """The transient dam: the problem, the model, the mesh, the initial condition,
    the frame series, and the histories that put numbers on the lag.

    Printed rather than drawn: the mesh, the steady discharge and where its phreatic
    surface sits, the saved-frame schedule, the mass-balance ledger the solver
    reports, the boundary-flow decay, and the head at the three history nodes at
    every saved instant. Every number the page quotes is on one of these lines.
    """
    import numpy as np

    from xslope.plot_seep import plot_seep_data, plot_seep_solution

    sd = load_slope_data(SEEP03)
    _u = declared_unit_labels(sd)
    tseep = sd["tseep"]
    print("   schedule    pool %s at t %s · duration %g · save_interval %g · "
          "save_times %s"
          % (tseep["series"]["pool"], tseep["times"], tseep["duration"],
             tseep["save_interval"], tseep["save_times"]))
    print("   staging     stage_1 %s · stage_2 %s · stability_time %s  (SEEP-3 is a "
          "seepage page; the staging belongs to the stability tutorial)"
          % (tseep.get("stage_1"), tseep.get("stage_2"),
             tseep.get("stability_time")))

    capture("seep03_schedule.png", _seep03_schedule, sd)
    capture("seep03_inputs.png", plot_inputs, sd, mode="seep",
            title="Seepage Model Inputs", frame="content", show_mesh=False)

    mesh = _seep03_mesh(sd)
    figsize = _seep02_figsize(mesh)
    nodes = np.asarray(mesh["nodes"])
    _width = max(nodes[:, 0]) - min(nodes[:, 0])
    print("   mesh        %d nodes · %d elements · %s at width/%d = %.5g %s"
          % (len(nodes), len(mesh["elements"]), SEEP03_ELEMENT, SEEP03_DIVISIONS,
             _width / SEEP03_DIVISIONS, _u["length"]))

    # ---- the initial condition, as an ordinary steady run ------------------ #
    steady_model = _seep03_steady_model(sd)
    seep_data, steady, log = _seep03_solve(steady_model, mesh)
    stats = _seep02_log_stats(log)
    capture("seep03_mesh.png", plot_seep_data, seep_data, figsize=figsize,
            show_bc=True)
    capture("seep03_steady.png", plot_seep_solution, seep_data, steady,
            figsize=figsize, levels=SEEP03_LEVELS, base_mat=SEEP03_BASE_MAT,
            fill_contours=True, mesh=False)
    head = np.asarray(steady["head"])
    psi = np.asarray(steady["u"]) / sd["gamma_water"]
    print("   steady      q %.6f %s · head %.3f to %.3f %s · u %.1f to %.1f %s"
          % (steady["flowrate"], _u["flowrate"], head.min(), head.max(),
             _u["length"], np.min(steady["u"]), np.max(steady["u"]),
             _u["stress"]))
    print("   iteration   %s" % stats)
    face = np.where(np.asarray(seep_data["bc_type"]) == 2)[0]
    wet = face[psi[face] > -1e-6]
    print("   exit face   %d nodes from (%.2f, %.2f) to (%.2f, %.2f) · %d wet%s"
          % (len(face), nodes[face][np.argmin(nodes[face, 1])][0],
             nodes[face, 1].min(),
             nodes[face][np.argmax(nodes[face, 1])][0], nodes[face, 1].max(),
             len(wet),
             "" if len(wet) else " — no seepage face at full pool; the phreatic "
             "surface reaches the downstream boundary at the tailwater"))
    ph0 = _seep03_phreatic(sd, mesh, steady)
    print("   stations    %s" % " ".join("%6g" % x for x in SEEP03_STATIONS))
    print("   steady ph   %s" % " ".join("%6.2f" % y for y in ph0))

    # ---- the march ---------------------------------------------------------- #
    tseep_data, solution, frames, tlog = _seep03_march(sd, mesh)
    times = [float(f["time"]) for f in frames]
    print("   march       %d saved frames at t = %s · converged %s · %d steps"
          % (len(frames), " ".join("%g" % t for t in times),
             solution["converged"], len(solution["dt_history"])))
    mb = solution["mass_balance"]
    print("   mass bal    cumulative inflow %.5g · final stored change %.5g · "
          "final closure %.3e"
          % (mb["cumulative_inflow"], mb["final_stored_change"],
             mb["final_closure"]))
    for row in mb["per_frame"]:
        print("   t %-6g    stored change %10.4f · cumulative inflow %10.4f · "
              "closure %.2e" % (row["time"], row["stored_change"],
                                row["cumulative_inflow"], row["closure"]))
    print("   solver log (the lines the page quotes):")
    for line in tlog.strip().splitlines():
        if "frame saved" in line or line.startswith("Transient"):
            print("     %s" % line.strip())

    outflow = [float(f.get("outflow") or 0.0) for f in frames]
    inflow = [float(f.get("inflow") or 0.0) for f in frames]
    peak = max(outflow)
    print("   %-8s %10s %10s %10s"
          % ("t (%s)" % _u["time"], "pool", "inflow", "outflow"))
    for t, qi, qo in zip(times, inflow, outflow):
        print("   %-8g %10.2f %10.5f %10.5f" % (t, _seep03_pool(sd, t), qi, qo))
    print("   outflow     peaks at %.5f on t = %g, and is %.5f at t = %g "
          "(%.2f%% of the peak)"
          % (peak, times[outflow.index(peak)], outflow[-1], times[-1],
             100.0 * outflow[-1] / peak))

    # ---- the frame series --------------------------------------------------- #
    all_head = np.concatenate([np.asarray(f["head"], float) for f in frames])
    vmin, vmax = float(all_head.min()), float(all_head.max())
    print("   color scale pinned across the series to head %.3f – %.3f %s"
          % (vmin, vmax, _u["length"]))
    panels = []
    for t, label in SEEP03_FRAMES:
        i = int(np.argmin(np.abs(np.asarray(times) - t)))
        panels.append(capture("seep03_frame_t%g.png" % times[i],
                              _seep03_frame_panel, tseep_data, frames[i], label,
                              vmin, vmax, _seep03_pool(sd, times[i])))
    _stack_panels("seep03_frames.png", panels)

    # ---- the phreatic surface, frame by frame -------------------------------- #
    print("   phreatic surface elevation (%s) by station and saved instant"
          % _u["length"])
    print("   %-8s %-6s %s" % ("t (%s)" % _u["time"], "pool",
                               " ".join("%6g" % x for x in SEEP03_STATIONS)))
    for f in frames:
        ys = _seep03_phreatic(sd, mesh, f)
        print("   %-8g %-6.2f %s"
              % (f["time"], _seep03_pool(sd, f["time"]),
                 " ".join("%6.2f" % y for y in ys)))

    # ---- the head histories -------------------------------------------------- #
    picked = _seep03_pick_nodes(mesh, tseep_data)
    for label, i, coords, _color, mats in picked:
        print("   node        %-17s -> %d at (%.3f, %.3f), material %s"
              % (label, i, coords[0], coords[1], mats))
    series = [[float(np.asarray(f["head"])[i]) for f in frames]
              for _l, i, _c, _col, _m in picked]
    pressures = [[float(np.asarray(f["u"])[i]) / sd["gamma_water"] for f in frames]
                 for _l, i, _c, _col, _m in picked]
    print("   %-8s %-8s %s" % ("t (%s)" % _u["time"], "pool",
                               " ".join("%18s" % l for l, *_ in picked)))
    for j, t in enumerate(times):
        print("   %-8g %-8.2f %s"
              % (t, _seep03_pool(sd, t),
                 " ".join("%18.4f" % s[j] for s in series)))
    print("   pressure head ψ = u/γw (%s) at the same nodes" % _u["length"])
    for j, t in enumerate(times):
        print("   %-8g %-8.2f %s"
              % (t, _seep03_pool(sd, t),
                 " ".join("%18.4f" % p[j] for p in pressures)))
    capture("seep03_history.png", _seep03_history_figure, sd, times, series,
            picked)


# --------------------------------------------------------------------------- #
# SEEP-4 — infiltration and flux boundaries
# --------------------------------------------------------------------------- #
# The Fredlund & Rahardjo (1993) earth dam under steady rainfall: 12 m high, 4 m
# crest, symmetric 2:1 faces, 52 m base, reservoir at 10 m, a 12 m horizontal
# drain at the downstream toe.  It is verification case 4 of GW6, and the
# tutorial's file pair is built from that corpus file by
# ``tools/build_dam_infiltration.py``.
#
# The lesson is a CONTRAST, so every figure and every number here is measured
# twice on one mesh: once with the rain and once without it.  The dry run is the
# completed model with its flux boundaries removed, which is bit-for-bit the same
# model as GW6 case 1 (gw006a) — identical geometry, identical soil — so the
# comparison is the manual's own pair of cases rather than a variation invented
# for the page.  Both runs are meshed and solved at the corpus discretization the
# verification page locks (tri3, 1.0 m, max_iter 2000), so a number printed here
# is a number that page can be checked against.
#
# Head figures are drawn on ONE pinned color scale, computed over both solutions,
# so the dry and wet panels are readable side by side: a contour that moves has
# moved, not been rescaled.
# --------------------------------------------------------------------------- #
SEEP04 = os.path.join(REPO_ROOT,
                      "docs/tutorials/files/xslope_dam_infiltration.xlsx")
#: The mesh, in the Build mesh dialog's own controls, and the solver's iteration
#: cap — the corpus discretization the GW6 verification tags run at.
SEEP04_SIZE = 1.0
SEEP04_ELEMENT = "tri3"
SEEP04_MAX_ITER = 2000
#: The manual's line 1-1: the crest centerline, where every GW6 case publishes its
#: pressure-head profile.  The dam is 24–28 m across the crest, so 26 is its axis.
SEEP04_LINE_X = 26.0
#: Elevations the line 1-1 profile is tabulated at — Fig 6.18's own marker spacing,
#: base to crest.  0.05 rather than 0 because the base node itself sits on the
#: no-flow boundary; it is the station the verification tag locks.
SEEP04_LINE_Y = (0.05, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
                 12.0)
#: Vertical stations the phreatic surface is read at, upstream toe to downstream
#: toe.  20 is the waterline, 24 and 28 are the crest shoulders, and 40 is the
#: upstream end of the toe drain, so the set brackets every feature that moves.
SEEP04_STATIONS = (4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 26.0, 28.0, 32.0, 36.0,
                   40.0, 44.0, 48.0)
#: Contour count for the head figures.
SEEP04_LEVELS = 12
#: The vendor's rain rate — a VERTICAL Darcy velocity, in m/s — and the 2:1 face
#: slope it falls on, both restated here so the projection the page teaches is
#: arithmetic the producer prints rather than a number copied into prose.
SEEP04_RAIN = 1e-8
SEEP04_FACE_RUN, SEEP04_FACE_RISE = 2.0, 1.0
#: The VENDOR's infiltration extents, for the comparison the page closes on. They
#: stop one vendor mesh edge short of the waterline and one short of the toe, so
#: they are properties of the vendor's mesh rather than of the dam; the tutorial
#: file draws the same rain on the geometry instead, waterline (20, 10) to toe
#: (52, 0). The corpus file gw006d keeps these, exactly as the vendor wrote them.
SEEP04_VENDOR_VERTICES = ((22.0, 11.0), (24.0, 12.0), (28.0, 12.0), (50.0, 1.0))
#: A flux block laid entirely inside the reservoir head boundary, used to prove the
#: collision rule: every node it loads is a specified-head node, so the answer must
#: not move. Wholly below the waterline at (20, 10), on the submerged upstream face.
SEEP04_SUBMERGED_BLOCK = ((12.0, 6.0), (20.0, 10.0))
#: Multiples of the file's rain rate the sweep solves at.  They span an order of
#: magnitude and a half around the vendor's rate, which is the middle of the set,
#: so the sweep brackets it rather than extending from it.
SEEP04_SWEEP_FACTORS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0)
#: Stations the swept free surfaces are drawn at: half-metre spacing from just
#: downstream of the waterline to the drain's upstream lip at x = 40, then tenth-
#: metre spacing out over the drain.  Denser than the thirteen the tables use,
#: because these are curves rather than readings, and denser again past the lip
#: because that is where the surfaces separate — the dry one turns down at the lip
#: itself while the wet ones stay up over the first metre or two of drain.
SEEP04_SWEEP_STATIONS = (tuple(20.5 + 0.5 * i for i in range(40))
                         + tuple(round(40.1 + 0.1 * i, 1) for i in range(25)))
#: The dam in profile — the ground surface the file carries, closed along the
#: impermeable base.  Drawn by hand because the sweep figure is a section with
#: curves on it, not a solution plot.
SEEP04_SECTION = ((0.0, 0.0), (24.0, 12.0), (28.0, 12.0), (52.0, 0.0))
#: One color per rain rate, cool to warm as the rate rises, so the six curves read
#: in order without the legend.  Same family the other seepage figures draw from.
SEEP04_SWEEP_COLORS = ("#2166ac", "#4393c3", "#3f8f5a", "#c9a227", "#c1663a",
                       "#96201f")


def _seep04_mesh(model, size=SEEP04_SIZE, element_type=SEEP04_ELEMENT):
    """One mesh of the model's material polygons, quiet — the corpus mesh."""
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)

    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(
            get_material_polygons(model), size, element_type,
            size_regions=extract_size_regions(model))


def _seep04_solve(model, mesh, max_iter=SEEP04_MAX_ITER):
    """One steady seepage solve on a built mesh, with the solver's log captured.

    ``tol=1e-5`` is the tolerance the page's own verification tags solve at
    (``run_seep_head_test``), so the heads printed here are the heads those tags
    compare against.
    """
    from xslope.seep import build_seep_data, run_seepage_analysis

    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        seep_data = build_seep_data(mesh, model)
        solution = run_seepage_analysis(seep_data, tol=1e-5, max_iter=max_iter)
    return seep_data, solution, log.getvalue()


def _seep04_dry(model):
    """The completed model with the rain switched off, in memory.

    Everything else — the reservoir head, the toe-drain exit face, the soil, the
    section — is untouched, so the difference between this run and the completed
    file's is the infiltration and nothing else.
    """
    bc = dict(model["seepage_bc"])
    bc["specified_fluxes"] = []
    return dict(model, seepage_bc=bc)


def _seep04_line(mesh, solution, x=SEEP04_LINE_X, elevations=SEEP04_LINE_Y):
    """Pressure head up the vertical line at ``x``, at the given elevations.

    Read off a linear interpolator over the mesh rather than off the nearest node,
    because the two runs are compared station by station and the difference at the
    base is smaller than the 1 m element.
    """
    import numpy as np

    ihead = _seep02_interp(mesh, np.asarray(solution["head"], dtype=float))
    xs = np.full(len(elevations), float(x))
    heads = np.ma.filled(ihead(xs, np.asarray(elevations, dtype=float)), np.nan)
    return np.asarray(heads), np.asarray(heads) - np.asarray(elevations)


def _seep04_top(x):
    """The dam's surface elevation at ``x`` — 2:1 faces up to a 4 m crest at 12 m."""
    if x <= 24.0:
        return x / 2.0
    if x <= 28.0:
        return 12.0
    return (52.0 - x) / 2.0


def _seep04_phreatic(mesh, solution, stations=SEEP04_STATIONS):
    """The phreatic-surface elevation at each station — the pressure-head zero
    crossing, walked up a vertical line and interpolated between the samples that
    straddle it, the same reading SEEP-2's producer takes.

    A column saturated all the way to the dam surface has no crossing inside the
    domain and returns NaN; on this dam that is the submerged upstream face.
    """
    import numpy as np

    ihead = _seep02_interp(mesh, np.asarray(solution["head"], dtype=float))
    out = []
    for x in stations:
        top = _seep04_top(x)
        ys = np.linspace(0.01, top - 0.01, 4000)
        psi = np.ma.filled(ihead(np.full(len(ys), float(x)), ys), np.nan) - ys
        ok = ~np.isnan(psi)
        ys, psi = ys[ok], psi[ok]
        wet = np.where(psi >= 0.0)[0]
        if not len(wet) or wet[-1] + 1 >= len(psi):
            out.append(float("nan"))
            continue
        i = wet[-1]
        out.append(float(ys[i] + (ys[i + 1] - ys[i]) * psi[i]
                         / (psi[i] - psi[i + 1])))
    return out


def _seep04_log_line(log, needle):
    """The last line of a solver log containing ``needle``, stripped."""
    hits = [ln.strip() for ln in log.splitlines() if needle in ln]
    return hits[-1] if hits else ""


def _seep04_exit_active(log):
    """The solver's own final report of how much of the exit face is draining.

    The unconfined iteration prints ``n/m exit face active`` on every logged
    sweep; the last one is the state the converged solution is in.
    """
    line = _seep04_log_line(log, "exit face active")
    if not line:
        return ""
    return line.split("relax")[-1].split(",")[-1].strip()


def _seep04_vendor_extent(model, vertices=SEEP04_VENDOR_VERTICES):
    """The completed model with the rain redrawn on the VENDOR's extents.

    The verification file stops one vendor mesh edge short of the waterline and
    one short of the toe — extents that are properties of the vendor's mesh, not
    of the dam.  Rebuilt here at the tutorial file's own rates so the only
    difference between the two runs is where the boundary starts and stops.
    """
    rates = {}
    for f in model["seepage_bc"]["specified_fluxes"]:
        rates["level" if f["coords"][0][1] == f["coords"][-1][1]
              else "sloping"] = f["flux"]
    blocks = [{"flux": rates["level"] if a[1] == b[1] else rates["sloping"],
               "coords": [a, b]}
              for a, b in zip(vertices, vertices[1:])]
    bc = dict(model["seepage_bc"], specified_fluxes=blocks)
    return dict(model, seepage_bc=bc)


def _seep04_submerged_block(model, coords=SEEP04_SUBMERGED_BLOCK):
    """The completed model with one extra flux block laid ENTIRELY inside the
    reservoir head boundary.

    Every node it loads is a specified-head node, so if the Dirichlet rows really
    do discard the flux they carry, this model must solve to exactly the same
    answer as the one without it.  That is the falsifiable form of the collision
    rule the page teaches, and the producer prints the comparison.
    """
    rate = max(f["flux"] for f in model["seepage_bc"]["specified_fluxes"])
    bc = dict(model["seepage_bc"])
    bc["specified_fluxes"] = list(bc["specified_fluxes"]) + [
        {"flux": rate, "coords": list(coords)}]
    return dict(model, seepage_bc=bc)


def _seep04_scaled_rain(model, factor):
    """The completed model with every flux block multiplied by ``factor``.

    All three blocks are scaled together, so the projection the file carries — the
    slope segments at cos(atan 1/2) of the crest rate — survives the scaling and
    the boundary still represents one vertical rain rate.  The workbook is never
    touched; ``factor = 0`` is the dry model.
    """
    bc = dict(model["seepage_bc"])
    bc["specified_fluxes"] = [dict(f, flux=f["flux"] * factor)
                              for f in model["seepage_bc"]["specified_fluxes"]]
    return dict(model, seepage_bc=bc)


def _seep04_free_surface(mesh, solution, stations=SEEP04_SWEEP_STATIONS):
    """The free surface at each station, and whether the column reached it.

    Returns ``(elevations, saturated)``.  Where the pressure head crosses zero
    inside the dam the elevation is that crossing — the phreatic surface.  Where
    the whole column is at or above zero the water has reached the dam surface and
    there is no crossing to find; the elevation returned is the dam surface itself
    and the flag is True, because the free surface has emerged rather than
    vanished.  A column that is dry throughout returns NaN.

    ``_seep04_phreatic`` returns NaN for both of those cases, which is the right
    answer for a table of phreatic elevations and the wrong one for a curve.
    """
    import numpy as np

    ihead = _seep02_interp(mesh, np.asarray(solution["head"], dtype=float))
    ys_out, sat_out = [], []
    for x in stations:
        top = _seep04_top(x)
        ys = np.linspace(0.01, top - 0.01, 4000)
        psi = np.ma.filled(ihead(np.full(len(ys), float(x)), ys), np.nan) - ys
        ok = ~np.isnan(psi)
        ys, psi = ys[ok], psi[ok]
        wet = np.where(psi >= 0.0)[0]
        if not len(wet):
            ys_out.append(float("nan"))
            sat_out.append(False)
        elif wet[-1] + 1 >= len(psi):
            ys_out.append(top)
            sat_out.append(True)
        else:
            i = wet[-1]
            ys_out.append(float(ys[i] + (ys[i + 1] - ys[i]) * psi[i]
                                / (psi[i] - psi[i + 1])))
            sat_out.append(False)
    return ys_out, sat_out


def _seep04_base_landing(mesh, solution, start=40.0, stop=46.0, step=0.002,
                         offset=0.01):
    """The x at which the free surface reaches the base, out along the drain.

    Read as the furthest x whose pressure head is still non-negative a hair above
    the base.  The hair is needed: the active drain nodes are held at head = y, so
    ψ is identically zero along the base itself and there is no sign change to find
    there.  ``offset`` is 1% of the 1 m element, and the answer moves by at most
    0.05 m as it is varied over 0.005–0.1 m — a tenth of one element, invisible at
    section scale.
    """
    import numpy as np

    ihead = _seep02_interp(mesh, np.asarray(solution["head"], dtype=float))
    xs = np.arange(start, stop, step)
    psi = np.ma.filled(ihead(xs, np.full(len(xs), offset)), np.nan) - offset
    ok = np.where(np.isfinite(psi) & (psi >= 0.0))[0]
    return float(xs[ok[-1]]) if len(ok) else float(start)


def _seep04_section_axes(ax):
    """One frame for the hand-drawn section: the dam, its impermeable base, the
    reservoir, and the toe drain, at equal aspect."""
    xs = [p[0] for p in SEEP04_SECTION]
    ys = [p[1] for p in SEEP04_SECTION]
    ax.plot(xs, ys, color="#4a5560", lw=1.6)
    ax.plot([SEEP04_SECTION[0][0], SEEP04_SECTION[-1][0]], [0.0, 0.0],
            color="#8a939c", lw=1.0)
    # The reservoir standing on the upstream face, with the water-table symbol at
    # its surface, and the toe drain that every free surface below drains into.
    head = 10.0
    ax.plot([0.0, 20.0], [head, head], color="#6ba8d6", lw=1.1)
    ax.plot([10.0], [head], marker="v", ms=7, color="#6ba8d6",
            markeredgecolor="#6ba8d6")
    ax.plot([40.0, 52.0], [0.0, 0.0], color="#3f8f5a", lw=3.0,
            solid_capstyle="butt")
    ax.annotate("toe drain", xy=(46.0, 0.0), xytext=(46.0, 1.2),
                color="#3f8f5a", fontsize=8, ha="center")
    ax.set_xlim(-1.0, 53.0)
    ax.set_ylim(-0.6, 13.9)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("elevation (m)")


def _seep04_sweep_figure(series, entry=(20.0, 10.0)):
    """The free surface on the dam section at each rain rate, all six on one axes.

    Each curve is anchored at both ends.  Upstream it starts at the reservoir
    waterline, which pins it whatever the rain does.  Downstream it ends ON the
    base, at the point where that run's own saturation runs out: the dry dam turns
    down at the drain's upstream lip, and the wetter the dam the further out over
    the drain the surface stays up before coming down.  The last stretch is steep
    because a free surface meeting a horizontal drain meets it at a right angle —
    the entry condition the Kozeny parabola is built on — and how steep is not
    resolved by a 1 m mesh, so it is drawn as what the stations report rather than
    smoothed into something the solution does not say.

    A rate whose curve lies on the dam surface has saturated the dam to its
    boundary; it is drawn dashed, because there the specified-flux boundary is
    offering more water than the soil can take and the run has stopped being a
    model of rain falling on a slope.
    """
    import numpy as np

    fig, ax = plt.subplots(figsize=(11.0, 4.4))
    _seep04_section_axes(ax)
    emerged = None
    for (label, ys, sat, land), color in zip(series, SEEP04_SWEEP_COLORS):
        pts = [entry] + [(x, y) for x, y in zip(SEEP04_SWEEP_STATIONS, ys)
                         if np.isfinite(y)]
        # The terminus is on the base, at whichever of the two readings of "the
        # surface is down" lies further out.  They can differ by up to an element:
        # over an active drain the base is held at atmospheric pressure and can
        # finish unsaturated while a saturated body is still draining onto it from
        # above, so the base scan stops before the column walk does.  Taking the
        # further of the two ends the curve where the saturation actually runs out
        # and keeps it from doubling back on itself.
        pts.append((max(pts[-1][0], land), 0.0))
        on_surface = any(sat)
        style = dict(lw=2.4, dashes=(4.5, 2.0)) if on_surface else dict(lw=1.8)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", color=color,
                label=label, **style)
        if on_surface:
            emerged = color
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=6,
              frameon=False, fontsize=9, columnspacing=1.6, handlelength=2.2)
    if emerged is not None:
        ax.annotate("dashed: the dam is saturated to its surface, so the\n"
                    "specified-flux boundary would pond in reality",
                    xy=(22.6, 11.3), xytext=(2.0, 12.7), color=emerged,
                    fontsize=8, ha="left", va="center",
                    arrowprops=dict(arrowstyle="->", color=emerged, lw=0.9,
                                    shrinkA=6, shrinkB=4))
    fig.tight_layout()
    plt.show()


def seep04_plots():
    """The infiltrated dam: the model, the mesh, the two solutions and the profile
    that compares them.

    Printed rather than drawn: the projection arithmetic behind the two flux rates,
    the mesh, each run's discharge and mass balance, the solver's iteration count
    and exit-face activity, the line 1-1 profile station by station with the rise
    the rain causes, and the phreatic surface at thirteen stations for each run.
    Every number the page quotes is on one of these lines.
    """
    import numpy as np

    from xslope.plot_seep import plot_seep_data, plot_seep_solution

    sd = load_slope_data(SEEP04)
    _u = declared_unit_labels(sd)
    mat = sd["materials"][0]
    print("   soil        %s · k %g/%g %s/s · %s α=%g n=%g · units %s"
          % (mat["name"], mat["k1"], mat["k2"], _u["length"], mat["unsat"],
             mat["vg_a"], mat["vg_n"], sd["unit_system"]))

    # ---- the boundary set, and the projection behind the two rates ---------- #
    bc = sd["seepage_bc"]
    for h in bc["specified_heads"]:
        print("   head        %g %s over %s" % (h["head"], _u["length"],
                                                h["coords"]))
    slope = math.degrees(math.atan2(SEEP04_FACE_RISE, SEEP04_FACE_RUN))
    projected = SEEP04_RAIN * SEEP04_FACE_RUN / math.hypot(SEEP04_FACE_RUN,
                                                           SEEP04_FACE_RISE)
    print("   projection  vertical rain %g %s/s on a %g:%g face (%.4f° from "
          "horizontal) has normal component %g × %g/√%g = %.8e"
          % (SEEP04_RAIN, _u["length"], SEEP04_FACE_RUN, SEEP04_FACE_RISE, slope,
             SEEP04_RAIN, SEEP04_FACE_RUN,
             SEEP04_FACE_RUN ** 2 + SEEP04_FACE_RISE ** 2, projected))
    footprint = 0.0
    for f in bc["specified_fluxes"]:
        (x0, y0), (x1, y1) = f["coords"][0], f["coords"][-1]
        length = math.hypot(x1 - x0, y1 - y0)
        footprint += abs(x1 - x0)
        print("   flux        %.8e over (%g, %g)→(%g, %g) · length %.4f %s · "
              "q·L %.6e · horizontal footprint %g %s"
              % (f["flux"], x0, y0, x1, y1, length, _u["length"],
                 f["flux"] * length, abs(x1 - x0), _u["length"]))
    total_in = sum(f["flux"] * math.hypot(f["coords"][-1][0] - f["coords"][0][0],
                                          f["coords"][-1][1] - f["coords"][0][1])
                   for f in bc["specified_fluxes"])
    print("   rain total  Σ q·L %.6e %s³/s per %s over a %g %s horizontal "
          "footprint = %.6e (rain × footprint %.6e)"
          % (total_in, _u["length"], _u["length"], footprint, _u["length"],
             total_in, SEEP04_RAIN * footprint))
    print("   exit face   %s" % (bc["exit_face"],))

    # ---- one mesh, two runs ------------------------------------------------- #
    mesh = _seep04_mesh(sd)
    nodes = np.asarray(mesh["nodes"])
    print("   mesh        %d nodes · %d elements · %s at %g %s"
          % (len(nodes), len(mesh["elements"]), SEEP04_ELEMENT, SEEP04_SIZE,
             _u["length"]))

    # The mesher pins a node at every seepage-BC vertex
    # (xslope/mesh.py::add_seep_bc_points_to_polygons), so a model whose boundary
    # polylines end anywhere but a corner of the section gets a mesh of its own.
    # The dry run's vertices are a subset of the completed model's, so it shares
    # this mesh; the vendor-extent variant ends at (22, 11) and (50, 1) and does
    # not, so it is meshed on its own terms — which is what the corpus file
    # gw006d does, and the comparison has to be model against model, not model
    # against a mesh built for something else. The collision test deliberately
    # STAYS on this mesh: it adds one load and nothing else, so the answer moving
    # would have only one possible cause.
    runs = {}
    ven_model = _seep04_vendor_extent(sd)
    for name, model, on in (("dry", _seep04_dry(sd), mesh),
                            ("wet", sd, mesh),
                            ("vendor", ven_model, _seep04_mesh(ven_model)),
                            ("collide", _seep04_submerged_block(sd), mesh)):
        seep_data, solution, log = _seep04_solve(model, on)
        runs[name] = (seep_data, solution, log, on)
        head = np.asarray(solution["head"], dtype=float)
        fn = np.asarray(seep_data["flux_nodal"], dtype=float)
        applied = float(np.sum(fn))
        # What the Dirichlet rows throw away: the load on specified-head nodes,
        # and on the exit-face nodes that finished the solve draining.
        held = np.asarray(seep_data["bc_type"]) == 1
        drained = ((np.asarray(seep_data["bc_type"]) == 2)
                   & (head <= np.asarray(seep_data["nodes"])[:, 1] + 1e-9))
        discarded = float(np.sum(fn[held | drained]))
        print("   %-7s run %d nodes · q %.6e %s³/s per %s · head %.4f–%.4f %s · "
              "converged %s · Σ assembled flux load %.6e (of which %.6e lands "
              "on head / draining nodes and is discarded)"
              % (name, len(np.asarray(seep_data["nodes"])), solution["flowrate"],
                 _u["length"], _u["length"], head.min(), head.max(),
                 _u["length"], solution.get("converged"), applied, discarded))
        print("        %s" % _seep04_log_line(log, "Converged in"))
        print("        %s" % _seep04_log_line(log, "Flow closure check"))
        print("        exit face at closure: %s" % _seep04_exit_active(log))
        # Where the water crosses the boundary, from the nodal reactions: the
        # reservoir head boundary and the toe drain, each as a signed total
        # (+ into the domain).  This is what the rain moves.
        q_react = np.asarray(solution["q"], dtype=float)
        print("        boundary flow: reservoir head %+.6e · toe drain %+.6e "
              "%s³/s per %s"
              % (float(np.sum(q_react[held])), float(np.sum(q_react[drained])),
                 _u["length"], _u["length"]))

    dry_data, dry_sol, _dry_log, _ = runs["dry"]
    wet_data, wet_sol, _wet_log, _ = runs["wet"]

    # The dry run is GW6 case 1 (gw006a) — same section, same soil, no rain. Solved
    # here rather than asserted, so the page's "this is the manual's case 1" holds.
    gw006a = os.path.join(REPO_ROOT,
                          "docs/verification/files/rocscience_gw/gw006a.xlsx")
    if os.path.exists(gw006a):
        ref = load_slope_data(gw006a)
        _rd, ref_sol, _rl = _seep04_solve(ref, _seep04_mesh(ref))
        print("   dry == gw006a (GW6 case 1): q %.6e vs %.6e, Δ %.3e"
              % (dry_sol["flowrate"], ref_sol["flowrate"],
                 dry_sol["flowrate"] - ref_sol["flowrate"]))

    # ---- the collision rule, measured -------------------------------------- #
    # A flux block laid wholly inside the reservoir head boundary loads only
    # specified-head nodes.  Those rows are overwritten with the prescribed head
    # (xslope/seep.py::_dirichlet_system: b is seeded with the flux vector, then
    # b[dirichlet] = head), so the load is discarded and the answer must not move.
    col_data, col_sol, _col_log, _ = runs["collide"]
    col_head = np.asarray(col_sol["head"], dtype=float)
    wet_head = np.asarray(wet_sol["head"], dtype=float)
    extra = (float(np.sum(np.asarray(col_data["flux_nodal"], dtype=float)))
             - float(np.sum(np.asarray(wet_data["flux_nodal"], dtype=float))))
    print("   collision   an extra flux block over %s (wholly inside the head "
          "boundary) assembles %.6e of load"
          % (list(SEEP04_SUBMERGED_BLOCK), extra))
    print("               q %.6e vs %.6e (Δ %.3e) · max |Δhead| %.3e %s "
          "→ the load on a specified-head node is discarded"
          % (col_sol["flowrate"], wet_sol["flowrate"],
             col_sol["flowrate"] - wet_sol["flowrate"],
             float(np.max(np.abs(col_head - wet_head))), _u["length"]))
    shared = int(np.argmin((nodes[:, 0] - 20.0) ** 2 + (nodes[:, 1] - 10.0) ** 2))
    print("               shared waterline node %d at (%.3f, %.3f): bc_type %d "
          "(1 = specified head) · assembled load %.6e · solved head %.6f %s "
          "(prescribed 10)"
          % (shared, nodes[shared, 0], nodes[shared, 1],
             int(np.asarray(wet_data["bc_type"])[shared]),
             float(np.asarray(wet_data["flux_nodal"], dtype=float)[shared]),
             wet_head[shared], _u["length"]))

    # ---- the extent comparison the page closes on --------------------------- #
    ven_data, ven_sol, _ven_log, ven_mesh = runs["vendor"]
    print("   extents     tutorial (geometry: waterline→toe) q %.6e vs vendor "
          "(mesh-edge-tied) q %.6e · Δ %+.3e (%+.2f%%)"
          % (wet_sol["flowrate"], ven_sol["flowrate"],
             wet_sol["flowrate"] - ven_sol["flowrate"],
             100.0 * (wet_sol["flowrate"] - ven_sol["flowrate"])
             / ven_sol["flowrate"]))
    _wl = sum(f["flux"] * math.hypot(f["coords"][-1][0] - f["coords"][0][0],
                                     f["coords"][-1][1] - f["coords"][0][1])
              for f in sd["seepage_bc"]["specified_fluxes"])
    _vl = sum(f["flux"] * math.hypot(f["coords"][-1][0] - f["coords"][0][0],
                                     f["coords"][-1][1] - f["coords"][0][1])
              for f in _seep04_vendor_extent(sd)["seepage_bc"]["specified_fluxes"])
    _wa = float(np.sum(np.asarray(wet_data["flux_nodal"], dtype=float)))
    _va = float(np.sum(np.asarray(ven_data["flux_nodal"], dtype=float)))
    print("               rain offered Σq·L %.6e over a 32 m footprint vs %.6e "
          "over 28 m · Δ %+.3e" % (_wl, _vl, _wl - _vl))
    print("               rain ASSEMBLED %.6e (%.2f%% of offered) vs %.6e "
          "(%.2f%%) — each on its own mesh, which pins a node at every BC "
          "vertex, so neither loses length at its ends"
          % (_wa, 100.0 * _wa / _wl, _va, 100.0 * _va / _vl))
    # What the vendor extents cost when they are NOT the extents the mesh was
    # built around: an edge carries load only when both its corners lie on the
    # polyline, so endpoints landing part-way along an edge drop it whole.
    _off_data, _off_sol, _off_log = _seep04_solve(ven_model, mesh)
    _off_a = float(np.sum(np.asarray(_off_data["flux_nodal"], dtype=float)))
    print("               the same vendor extents on the tutorial file's mesh "
          "(no node at (22, 11) or (50, 1)) assemble only %.6e (%.2f%% of "
          "offered) and give q %.6e — the loss is why a flux boundary is drawn "
          "on the geometry"
          % (_off_a, 100.0 * _off_a / _vl, _off_sol["flowrate"]))

    # ---- figures ------------------------------------------------------------ #
    # frame="content": the section is four times as long as it is deep, and the
    # default fill frame pads it out to the figure aspect and buries the dam.
    capture("seep04_inputs.png", plot_inputs, sd, mode="seep",
            title="Seepage Model Inputs", frame="content", show_mesh=False)
    figsize = _seep02_figsize(mesh)
    # The mesh figure shows the reader's OWN state at that point in the build:
    # head and exit face entered, no rain yet — so it is drawn from the dry
    # model's seep_data (same mesh; the completed model's flux vertices are all
    # already-pinned points), and no flux markers appear before the reader has
    # built a flux boundary.
    capture("seep04_mesh.png", plot_seep_data, dry_data, figsize=figsize,
            show_bc=True)

    # One color scale over BOTH solutions, so the two panels compare directly.
    both = np.concatenate([np.asarray(dry_sol["head"], dtype=float),
                           np.asarray(wet_sol["head"], dtype=float)])
    vmin, vmax = float(both.min()), float(both.max())
    print("   color scale pinned across both runs to head %.4f–%.4f %s"
          % (vmin, vmax, _u["length"]))
    for name, (data, sol, _log, _m) in (("dry", runs["dry"]), ("wet", runs["wet"])):
        capture("seep04_%s.png" % name, plot_seep_solution, data, sol,
                figsize=figsize, levels=SEEP04_LEVELS, base_mat=1,
                fill_contours=True, phreatic=True, flowlines=True, mesh=False,
                vmin=vmin, vmax=vmax)

    # ---- line 1-1, the manual's own comparison ------------------------------ #
    dry_h, dry_psi = _seep04_line(mesh, dry_sol)
    wet_h, wet_psi = _seep04_line(mesh, wet_sol)
    ven_h, ven_psi = _seep04_line(ven_mesh, ven_sol)
    print("   line 1-1 (x = %g %s): pressure head ψ = h − y"
          % (SEEP04_LINE_X, _u["length"]))
    print("   %-8s %10s %10s %10s %10s %10s %10s %10s"
          % ("y", "h dry", "h wet", "ψ dry", "ψ wet", "Δψ rain",
             "ψ vendor", "Δψ extent"))
    for y, hd, hw, pd_, pw, pv in zip(SEEP04_LINE_Y, dry_h, wet_h, dry_psi,
                                      wet_psi, ven_psi):
        print("   %-8g %10.4f %10.4f %10.4f %10.4f %+10.4f %10.4f %+10.4f"
              % (y, hd, hw, pd_, pw, pw - pd_, pv, pw - pv))

    # ---- what the rain does to the interior and to the free surface --------- #
    rise = (np.asarray(wet_sol["head"], dtype=float)
            - np.asarray(dry_sol["head"], dtype=float))
    i_max = int(np.argmax(rise))
    print("   head rise   mean %+.4f %s · max %+.4f at node %d (%.2f, %.2f)"
          % (rise.mean(), _u["length"], rise[i_max], i_max,
             nodes[i_max, 0], nodes[i_max, 1]))
    dry_ph = _seep04_phreatic(mesh, dry_sol)
    wet_ph = _seep04_phreatic(mesh, wet_sol)
    ven_ph = _seep04_phreatic(ven_mesh, ven_sol)
    print("   phreatic surface elevation (%s) by station" % _u["length"])
    print("   %-10s %s" % ("station", " ".join("%6g" % x
                                               for x in SEEP04_STATIONS)))
    for label, ys in (("dry", dry_ph), ("wet", wet_ph), ("vendor", ven_ph)):
        print("   %-10s %s" % (label, " ".join("%6.2f" % y for y in ys)))
    print("   %-10s %s" % ("Δ", " ".join("%+6.2f" % (w - d)
                                         for d, w in zip(dry_ph, wet_ph))))

    # ---- the rain sweep ----------------------------------------------------- #
    # The file's rates, scaled in memory, on the one mesh: what more rain does to
    # the discharge, to the water table under the crest, and to how much of the
    # dam stays unsaturated.  Every rate is reported as q/k, because that ratio —
    # not the rate itself — is what decides whether the soil can carry the rain
    # away, and the sweep runs up to the rate at which it cannot.
    k = float(sd["materials"][0]["k1"])
    print("   rain sweep at q/k = q ÷ k = q ÷ %g %s/s (all three flux blocks "
          "scaled together, in memory)" % (k, _u["length"]))
    print("   %-8s %-12s %-7s %-14s %-14s %-9s %-9s %-8s %s"
          % ("factor", "q (m/s)", "q/k", "Q", "reservoir", "WT x=26",
             "Δ from dry", "lands x", "exit face"))
    sweep, sweep_q, sweep_wt, sweep_res = [], [], [], []
    i26 = SEEP04_STATIONS.index(SEEP04_LINE_X)
    for factor in SEEP04_SWEEP_FACTORS:
        model = _seep04_scaled_rain(sd, factor)
        s_data, s_sol, s_log = _seep04_solve(model, mesh)
        ys, sat = _seep04_free_surface(mesh, s_sol)
        table_ph = _seep04_phreatic(mesh, s_sol)
        wt = table_ph[i26]
        q = float(s_sol["flowrate"])
        rain = SEEP04_RAIN * factor
        # How much of the drain's outflow the reservoir is supplying: the nodal
        # reactions on the specified-head boundary.  The rest is rain.
        res = float(np.sum(np.asarray(s_sol["q"], dtype=float)[
            np.asarray(s_data["bc_type"]) == 1]))
        # Where this run's free surface comes down onto the base, out along the
        # drain.  The drain's lip is at x = 40, so a landing further out is
        # saturated drain length: the dam is discharging over that much of it.
        land = _seep04_base_landing(mesh, s_sol)
        sweep.append(("q/k = %g" % round(rain / k, 6), ys, sat, land))
        sweep_q.append(q)
        sweep_wt.append(wt)
        sweep_res.append(res)
        print("   %-8g %-12.4e %-7g %-14.6e %-14.6e %-9s %-9s %-8.3f %s"
              % (factor, rain, round(rain / k, 6), q, res,
                 "surface" if np.isnan(wt) else "%.4f" % wt,
                 "—" if np.isnan(wt) or np.isnan(sweep_wt[0])
                 else "%+.4f" % (wt - sweep_wt[0]),
                 land, _seep04_exit_active(s_log)))
        n_sat = sum(1 for f in sat if f)
        if n_sat:
            xs_sat = [x for x, f in zip(SEEP04_SWEEP_STATIONS, sat) if f]
            print("        saturated to the dam surface over x = %g–%g %s "
                  "(%d of %d stations) — the flux boundary is offering more "
                  "than the soil can accept"
                  % (min(xs_sat), max(xs_sat), _u["length"], n_sat, len(sat)))
    # Linear in q or not, from the sweep's own numbers: a linear response has the
    # same rise per unit of rain everywhere, so the per-factor increments are the
    # test rather than the totals.
    print("   response per unit factor (a linear response repeats a constant):")
    for i in range(1, len(SEEP04_SWEEP_FACTORS)):
        f0, f1 = SEEP04_SWEEP_FACTORS[i - 1], SEEP04_SWEEP_FACTORS[i]
        dq = (sweep_q[i] - sweep_q[i - 1]) / (f1 - f0)
        dres = (sweep_res[i] - sweep_res[i - 1]) / (f1 - f0)
        dwt = (sweep_wt[i] - sweep_wt[i - 1]) / (f1 - f0)
        print("        factor %-5g → %-5g   ΔQ/Δfactor %.6e   "
              "Δreservoir/Δfactor %+.6e   ΔWT/Δfactor %s"
              % (f0, f1, dq, dres, "n/a" if np.isnan(dwt) else "%+.4f" % dwt))
    # The rain each run offers, against what the drain actually gains: the two
    # differ because the rising mound throttles the reservoir, so the drain never
    # gains a full unit of rain per unit of rain applied.
    offered = sum(f["flux"] * math.hypot(f["coords"][-1][0] - f["coords"][0][0],
                                         f["coords"][-1][1] - f["coords"][0][1])
                  for f in sd["seepage_bc"]["specified_fluxes"])
    i0 = SEEP04_SWEEP_FACTORS.index(0.0)
    i1 = SEEP04_SWEEP_FACTORS.index(1.0)
    print("        rain offered at factor 1 is %.6e; the drain gains %.6e of it "
          "because the reservoir gives up %.6e"
          % (offered, sweep_q[i1] - sweep_q[i0], sweep_res[i0] - sweep_res[i1]))
    capture("seep04_rain_sweep.png", _seep04_sweep_figure, sweep)


# --------------------------------------------------------------------------- #
# FEM-1 — Strength Reduction Basics
# --------------------------------------------------------------------------- #
#: FEM-1's file pair, written by ``tools/build_ssrm_embankment.py`` from the
#: Griffiths & Lane Example 1 corpus model.  The starter carries no elastic
#: constants, which is why the limit-equilibrium figures are drawn on IT — the
#: page's opening claim is that a factor of safety by slices needs neither E nor
#: nu, and drawing that state from the file that lacks them is the claim itself.
FEM01_START = os.path.join(REPO_ROOT,
                           "docs/tutorials/files/xslope_ssrm_embankment_start.xlsx")
FEM01_DONE = os.path.join(REPO_ROOT,
                          "docs/tutorials/files/xslope_ssrm_embankment.xlsx")
#: The method the page's limit-equilibrium number is quoted from.  Spencer
#: because it satisfies both equilibrium conditions, which is the closest
#: limit-equilibrium statement of the same problem the finite element run solves.
FEM01_METHOD = "spencer"
FEM01_SLICES = 40
#: The strength-reduction run, at the settings the page's final run is made at:
#: bracket [1.0, 2.0], tolerance 0.01, and 12,000 iterations a trial.
#: ``non_convergence`` rather than the API default ``hybrid`` because that is what
#: Studio's Run FEM dialog runs — its Failure criterion list offers no hybrid
#: entry, and its runner falls back to non-convergence.  On this model the two
#: criteria return the same factor of safety trial for trial.
FEM01_CRITERION = "non_convergence"
FEM01_F_MIN, FEM01_F_MAX = 1.0, 2.0
FEM01_TOLERANCE = 0.01
#: 12,000 — now ``solve_ssrm``'s own default, and the budget the page shows in
#: **Max iterations per trial**.  The figures have to show the state the page's own
#: instructions produce.  At the former 3,000 default the near-critical trials were
#: cut off mid-work and the answer read 1.3477; at 12,000 they finish and it
#: plateaus at 1.3633.  A budget below what the model needs is no longer a verdict
#: either — the engine extends it while the residual is still falling — so this
#: number now sets where the extension starts, not where the trial dies.
FEM01_MAX_ITERATIONS = 12000


def _fem01_mesh(model):
    """The tutorial's mesh, built at the element type and size the COMPLETED FILE
    declares rather than at numbers restated here — so the picture and the
    workbook's main sheet cannot disagree."""
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)

    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(get_material_polygons(model),
                                        model["target_size"],
                                        model["element_type"],
                                        size_regions=extract_size_regions(model))


def _fem01_search(model, method=FEM01_METHOD):
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, path, circles = circular_search(
            model, method, num_slices=FEM01_SLICES, diagnostic=False,
            **file_search_window(model))
    return fs_cache, path, circles


def fem01_plots():
    """The strength-reduction arc: the section, the limit-equilibrium answer, the
    mesh, and the finite element run's three result panels at failure.

    Printed rather than drawn: the limit-equilibrium search's critical circle, the
    mesh counts, the whole bisection walk trial by trial with its iteration
    counts, the factor of safety and the bracket it came from, and the two states
    the last pair of figures compares — the last trial that reached equilibrium
    and the captured collapse.  Every number the page quotes is on one of these
    lines.

    Every results panel is drawn autoscaled, exactly as Studio draws it — the
    converged state's figure matches the reader's screen, and the comparison
    between the states is carried by the color bar ranges, not a shared scale.
    """
    import time

    import numpy as np

    from xslope.fem import build_fem_data, solve_ssrm
    from xslope.plot_fem import (plot_fem_data, plot_fem_results,
                                 shared_panel_scales)

    # ---- the section, and the answer that needs no elastic constants -------- #
    start = load_slope_data(FEM01_START)
    _u = declared_unit_labels(start)
    mat = start["materials"][0]
    print("   starter     %s · γ %g %s · c %g %s · φ %g° · option %s · E %r · ν %r"
          % (mat["name"], mat["gamma"], _u["unit_weight"], mat["c"], _u["stress"],
             mat["phi"], mat["option"], mat["E"], mat["nu"]))
    print("   geometry    %s · one starting circle %s"
          % (start["ground_surface"], start["circles"]))

    capture("fem01_inputs.png", plot_inputs, start,
            title="Slope Geometry and Inputs")

    fs_cache, _path, _circles = _fem01_search(start)
    crit = fs_cache[0]
    capture("fem01_lem_solution.png", plot_solution, start, crit["slices"],
            crit["failure_surface"], crit["solver_result"])
    xs, ys = zip(*list(crit["failure_surface"].coords))
    print("   LEM         %s FS %.4f on Xo %.4f Yo %.4f R %.4f (tangent depth "
          "%.4f) · entry (%.3f, %.3f) exit (%.3f, %.3f) · %d candidates"
          % (FEM01_METHOD, crit["FS"], crit["Xo"], crit["Yo"],
             crit["Yo"] - crit["Depth"], crit["Depth"], xs[0], ys[0], xs[-1],
             ys[-1], len(fs_cache)))

    # ---- the mesh the strength reduction runs on ---------------------------- #
    done = load_slope_data(FEM01_DONE)
    dmat = done["materials"][0]
    print("   completed   same soil with E %g %s · ν %g · declares %s at target "
          "size %g %s" % (dmat["E"], _u["stress"], dmat["nu"],
                          done["element_type"], done["target_size"], _u["length"]))
    mesh = _fem01_mesh(done)
    nodes = np.asarray(mesh["nodes"])
    fem_data = build_fem_data(done, mesh)
    print("   mesh        %d nodes · %d elements · every element %s · side BC %s"
          % (len(nodes), len(mesh["elements"]), done["element_type"],
             done.get("side_bc") or "rollers (the shipped default; the file "
             "declares none)"))
    capture("fem01_mesh.png", plot_fem_data, fem_data)

    # ---- the strength reduction --------------------------------------------- #
    log = io.StringIO()
    t0 = time.time()
    with contextlib.redirect_stdout(log):
        result = solve_ssrm(fem_data, F_min=FEM01_F_MIN, F_max=FEM01_F_MAX,
                            tolerance=FEM01_TOLERANCE, debug_level=1,
                            failure_criterion=FEM01_CRITERION,
                            max_iterations=FEM01_MAX_ITERATIONS)
    wall = time.time() - t0
    print("   SSRM        FS %.4f from the bracket [%.4f, %.4f] (width %.4f) "
          "after %d bisection steps · %.1f s wall"
          % (result["FS"], result["final_interval"][0], result["final_interval"][1],
             result["interval_width"], result["iterations_ssrm"], wall))
    for tr in result["trials"]:
        print("        F %.4f  %-6s  %-13s  %s iterations"
              % (tr["F"], tr.get("role"), tr.get("verdict"), tr.get("iterations")))

    last, fail = result["last_solution"], result.get("failure_solution")
    for label, field in (("last converged", last), ("at failure", fail)):
        if field is None:
            print("   %-12s none" % label)
            continue
        u = np.asarray(field["displacements"]).reshape(-1, 2)
        umag = np.hypot(u[:, 0], u[:, 1])
        i = int(np.argmax(umag))
        ue = np.asarray(field["displacements_elastic"]).reshape(-1, 2)
        print("   %-11s F %.4f · %s in %d iterations (%s) · max|u| %.6f %s at "
              "(%.2f, %.2f), %.2f× the elastic %.6f"
              % (label, field.get("F", float("nan")),
                 "equilibrium" if field.get("converged") else "no equilibrium",
                 field.get("iterations"), field.get("exit_reason"), umag.max(),
                 _u["length"], nodes[i, 0], nodes[i, 1],
                 umag.max() / np.hypot(ue[:, 0], ue[:, 1]).max(),
                 np.hypot(ue[:, 0], ue[:, 1]).max()))

    # The three panels the FEM Results view offers, each drawn on its own so the
    # page can place them where its prose reaches them. Every one is at the
    # captured at-failure field, which is the state the mechanism is visible in.
    for name, kind in (("fem01_shear_strain.png", "shear_strain"),
                       ("fem01_deformed.png", "deformation"),
                       ("fem01_displacement_vectors.png", "displace_vector")):
        capture(name, plot_fem_results, fem_data, last, plot_type=kind,
                fs=result["FS"], failure_solution=fail, field_state="failure")

    # The comparison the page's failure definition rests on: the same section at
    # the last trial that reached equilibrium. Autoscaled to its own range, the
    # way Studio draws every field state — the reader's screen and this figure
    # agree, and the eightyfold gap between the two states' color bar maxima is
    # the comparison the page teaches.
    own = shared_panel_scales(fem_data, [last])
    both = shared_panel_scales(fem_data, [last, fail])
    print("   ranges      converged shear strain %s–%s · at failure %s–%s"
          % (own["vmin"], own["vmax"], both["vmin"], both["vmax"]))
    capture("fem01_shear_strain_converged.png", plot_fem_results, fem_data, last,
            plot_type="shear_strain", fs=result["FS"], failure_solution=fail,
            field_state="converged")


# --------------------------------------------------------------------------- #
# FEM-2 — Reinforcement: LEM against FEM
# --------------------------------------------------------------------------- #
#: FEM-2's file pair, written by ``tools/build_reinforced_slope_tutorial.py`` from
#: the LEM reinforced-slope model (docs/lem/samples.md problem 9, the model
#: Tutorial LEM-8 builds) and its finite element counterpart
#: (docs/fem/samples.md problem 1).  The starter carries the limit-equilibrium
#: model only — the six geogrid lines, but no elastic constants, no mesh and no
#: FEM reinforcement data; the completed file adds the soils' E and nu, Tres, the
#: bar modulus and the bar area, and declares the mesh.  Nothing here runs the FEM
#: on the starter: its three uses below are the inputs plot and the Spencer search.
FEM02_START = os.path.join(REPO_ROOT,
                           "docs/tutorials/files/xslope_reinforced_slope_start.xlsx")
FEM02_DONE = os.path.join(REPO_ROOT,
                          "docs/tutorials/files/xslope_reinforced_slope.xlsx")
#: The method the page's limit-equilibrium number is quoted from, and the slice
#: count the sample page's own locked search uses.
FEM02_METHOD = "spencer"
FEM02_SLICES = 40
#: The strength-reduction run, at the settings the page's final runs are made at.
#: ``non_convergence`` rather than the API default ``hybrid`` because that is what
#: Studio's Run FEM dialog runs.
FEM02_CRITERION = "non_convergence"
FEM02_F_MIN, FEM02_F_MAX = 1.0, 2.0
FEM02_TOLERANCE = 0.01
#: 12,000, which is also ``solve_ssrm``'s default budget.  On this model the
#: number no longer decides the answer: a trial that spends its budget while still
#: making progress is extended rather than failed, so both runs return the same
#: factors of safety (1.5586 elastic-perfectly-plastic, 1.5117 peak-residual) from
#: a 3,000 budget and from 12,000, with the same brackets and the same converged
#: field.  What the budget still decides is how far the captured failure state has
#: developed, which is why the producers ask for the larger one.
FEM02_MAX_ITERATIONS = 12000
#: The line the bar-force profiles are drawn for: the most heavily loaded of the
#: six in both runs (peak 800 lb/ft elastic-perfectly-plastic, 798 lb/ft
#: peak-residual, both at the last converged trial).
FEM02_PROFILE_LINE = 5
#: The residual capacities the sweep asks for, highest first.  ``None`` is the
#: blank cell — no post-peak drop at all.  The sweep runs down to zero because
#: that is where the answer moves: Tres = Tmax is the same run as a blank cell,
#: 600 and 400 give one lower factor of safety, and only Tres = 0 — the bar
#: tearing and carrying nothing — gives a third.
FEM02_TRES_SWEEP = (None, 800.0, 600.0, 400.0, 0.0)


def _fem02_mesh(model):
    """The tutorial's mesh, at the element type and target size the COMPLETED FILE
    declares, built the way Studio's Build Mesh dialog builds it — the
    reinforcement lines carried in as constraint lines so the bars land on mesh
    edges — with **Refine thin zones unticked**.

    The dialog's own default for that box is ticked, and on this section it is not
    inert: the shell is a 1.19 ft facing band, thinner than one element at the
    target size, so the refinement drives the local size there to 0.33 ft and the
    mesh from 2,101 elements to 5,096.  That mesh is a different model's worth of
    answer (peak-residual 1.4414 against 1.5117, elastic-perfectly-plastic 1.5117
    against 1.5586), so the page has to tell the reader which box to clear, and
    this producer draws the mesh the page's own numbers come from.
    """
    from xslope.mesh import (build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_size_regions, get_material_polygons)

    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(model)
    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(
            get_material_polygons(model, reinf_lines=lines),
            model["target_size"], model["element_type"], lines=lines or None,
            element_size_1d=model.get("element_size_1d"),
            size_regions=extract_size_regions(model))


def _fem02_solve(model, mesh, t_res, max_iterations=FEM02_MAX_ITERATIONS):
    """One strength-reduction run on this model with the six lines' residual
    capacity set to ``t_res`` (``None`` leaves the cells blank, which the loader
    reads as "no post-peak drop").  Returns ``(fem_data, result, log)``."""
    import copy as _copy

    from xslope.fem import build_fem_data, solve_ssrm

    sd = _copy.deepcopy(model)
    for line in sd["reinforcement_lines"]:
        line["t_res"] = float("nan") if t_res is None else float(t_res)
    fem_data = build_fem_data(sd, mesh)
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        result = solve_ssrm(fem_data, F_min=FEM02_F_MIN, F_max=FEM02_F_MAX,
                            tolerance=FEM02_TOLERANCE, debug_level=1,
                            failure_criterion=FEM02_CRITERION,
                            max_iterations=max_iterations)
    return fem_data, result, log.getvalue()


def _fem02_report(label, fem_data, result):
    """Print everything the page can quote from one run: the bracket walk, the two
    states, and what each of the six lines is doing at the last converged one."""
    import numpy as np

    from xslope import fem_details

    print("   %-14s FS %.4f from [%.6f, %.6f] (width %.6f) after %d bisection "
          "steps" % (label, result["FS"], result["final_interval"][0],
                     result["final_interval"][1], result["interval_width"],
                     result["iterations_ssrm"]))
    for tr in result["trials"]:
        print("        F %.4f  %-6s  %-13s  %s iterations"
              % (tr["F"], tr.get("role"), tr.get("verdict"), tr.get("iterations")))
    nodes = np.asarray(fem_data["nodes"])
    last, fail = result["last_solution"], result.get("failure_solution")
    for state, field in (("last converged", last), ("at failure", fail)):
        if field is None:
            print("        %-14s none" % state)
            continue
        u = np.asarray(field["displacements"]).reshape(-1, 2)
        umag = np.hypot(u[:, 0], u[:, 1])
        ue = np.asarray(field["displacements_elastic"]).reshape(-1, 2)
        i = int(np.argmax(umag))
        print("        %-14s F %.4f · %s in %s iterations (%s) · max|u| %.6f at "
              "(%.2f, %.2f), %.2f× elastic"
              % (state, field.get("F", float("nan")),
                 "equilibrium" if field.get("converged") else "no equilibrium",
                 field.get("iterations"), field.get("exit_reason"), umag.max(),
                 nodes[i, 0], nodes[i, 1],
                 umag.max() / np.hypot(ue[:, 0], ue[:, 1]).max()))
    mats_1d = np.asarray(fem_data.get("element_materials_1d", []), dtype=int)
    piles = np.asarray(fem_data.get("pile_elem_mask", np.zeros(len(mats_1d))),
                       dtype=bool)
    for line_id in sorted(set(mats_1d[~piles].tolist())):
        prof = fem_details.reinforcement_profile(
            fem_data, last, int(line_id), field_state="converged",
            failure_solution=fail)
        if not len(prof["s"]):
            continue
        # Two different "peaks", and the page needs both: the greatest force
        # anywhere on the bar, and the force at the point working hardest against
        # its own capacity — which on a bar with pullout ramps is usually an end
        # element carrying a fraction of Tmax against a fraction of the capacity.
        print("        line %d  max force %.1f · hardest-worked %.1f at s %.2f "
              "(utilization %.3f) · %s · %d element(s) at capacity, %d softened"
              % (line_id, float(prof["T"].max()), prof["peak_T"], prof["peak_s"],
                 prof["peak_utilization"], prof["status"],
                 int(prof["failed"].sum()), int(prof["softened"].sum())))


def fem02_plots():
    """The reinforcement arc: the section, the limit-equilibrium answer on the
    same six lines, the mesh the bars are meshed into, and the two strength
    reduction runs the page compares — elastic-perfectly-plastic (Tres blank) and
    peak-residual (Tres = 600 lb/ft).

    Printed rather than drawn: the critical circle, the mesh counts, both bracket
    walks trial by trial, the two states of each run, and every line's peak force
    and state at the last converged trial.
    """
    import time

    import numpy as np

    from xslope.plot_fem import plot_fem_data, plot_fem_results

    # ---- the section, and the answer the slices already gave ---------------- #
    start = load_slope_data(FEM02_START)
    _u = declared_unit_labels(start)
    for mat in start["materials"]:
        print("   starter     %-6s γ %g %s · c %g %s · φ %g° · E %r · ν %r"
              % (mat["name"], mat["gamma"], _u["unit_weight"], mat["c"],
                 _u["stress"], mat["phi"], mat["E"], mat["nu"]))
    r0 = start["reinforcement_lines"][0]
    print("   starter     %d reinforcement lines · Tmax %g · Lp %g/%g · "
          "Tres %r · E %r · Area %r"
          % (len(start["reinforcement_lines"]), r0["t_max"], r0["lp1"], r0["lp2"],
             r0["t_res"], r0["E"], r0["area"]))

    capture("fem02_inputs.png", plot_inputs, start,
            title="Slope Geometry and Inputs")

    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, _path, _circles = circular_search(
            start, FEM02_METHOD, num_slices=FEM02_SLICES, diagnostic=False,
            **file_search_window(start))
    crit = fs_cache[0]
    capture("fem02_lem_solution.png", plot_solution, start, crit["slices"],
            crit["failure_surface"], crit["solver_result"])
    xs, ys = zip(*list(crit["failure_surface"].coords))
    print("   LEM         %s FS %.4f on Xo %.4f Yo %.4f R %.4f · entry "
          "(%.3f, %.3f) exit (%.3f, %.3f) · %d candidates"
          % (FEM02_METHOD, crit["FS"], crit["Xo"], crit["Yo"],
             crit["Yo"] - crit["Depth"], xs[0], ys[0], xs[-1], ys[-1],
             len(fs_cache)))

    # ---- the mesh the bars are meshed into ---------------------------------- #
    done = load_slope_data(FEM02_DONE)
    d0 = done["reinforcement_lines"][0]
    print("   completed   Tres %g · bar E %g %s · Area %g %s²/%s · declares %s at "
          "target size %g %s"
          % (d0["t_res"], d0["E"], _u["stress"], d0["area"], _u["length"],
             _u["length"], done["element_type"], done["target_size"],
             _u["length"]))
    mesh = _fem02_mesh(done)
    print("   mesh        %d nodes · %d elements · %d bar elements · every "
          "element %s"
          % (len(mesh["nodes"]), len(mesh["elements"]),
             len(mesh.get("elements_1d", [])), done["element_type"]))
    from xslope.fem import build_fem_data
    capture("fem02_mesh.png", plot_fem_data, build_fem_data(done, mesh))

    # ---- the two runs -------------------------------------------------------- #
    runs = {}
    for tag, t_res in (("epp", None), ("pr", d0["t_res"])):
        t0 = time.time()
        fem_data, result, _log = _fem02_solve(done, mesh, t_res)
        print("   %s run · Tres %s · %.1f s wall (capture included)"
              % ("elastic-perfectly-plastic" if t_res is None else "peak-residual",
                 "blank" if t_res is None else "%g" % t_res, time.time() - t0))
        _fem02_report(tag, fem_data, result)
        runs[tag] = (fem_data, result)

    for tag, name in (("epp", "fem02_shear_strain_epp.png"),
                      ("pr", "fem02_shear_strain_pr.png")):
        fem_data, result = runs[tag]
        capture(name, plot_fem_results, fem_data, result["last_solution"],
                plot_type="shear_strain", fs=result["FS"],
                failure_solution=result.get("failure_solution"),
                field_state="failure")
    fem_data, result = runs["pr"]
    for name, kind in (("fem02_deformed_pr.png", "deformation"),
                       ("fem02_displacement_vectors_pr.png", "displace_vector")):
        capture(name, plot_fem_results, fem_data, result["last_solution"],
                plot_type=kind, fs=result["FS"],
                failure_solution=result.get("failure_solution"),
                field_state="failure")

    # ---- the bar the page reads along --------------------------------------- #
    # Studio's own 1D Details drawing, one figure per run: the two cannot share a
    # figure because ``plot_reinforcement_detail`` lays its panels out with
    # ``tight_layout``, which a matplotlib SubFigure does not carry.
    from xslope import fem_details
    from xslope.plot_fem_details import plot_reinforcement_detail

    for tag, name in (("epp", "fem02_bar_profile_epp.png"),
                      ("pr", "fem02_bar_profile_pr.png")):
        fem_data, result = runs[tag]
        prof = fem_details.reinforcement_profile(
            fem_data, result["last_solution"], FEM02_PROFILE_LINE,
            slope_data=done, field_state="converged",
            failure_solution=result.get("failure_solution"))
        print("   %-4s line %d · max force %.1f %s · hardest-worked %.1f at "
              "s %.2f %s (utilization %.3f) · %s · %d at capacity, %d softened"
              % (tag, FEM02_PROFILE_LINE, float(prof["T"].max()),
                 _u["force_per_len"], prof["peak_T"], prof["peak_s"],
                 _u["length"], prof["peak_utilization"], prof["status"],
                 int(prof["failed"].sum()), int(prof["softened"].sum())))
        capture(name, plot_reinforcement_detail, prof)

    # Where each engine puts the mechanism, measured rather than described: the
    # limit-equilibrium circle's entry and exit against the centroid of the
    # elements carrying the most viscoplastic shear strain.
    fem_data, result = runs["pr"]
    fail = result.get("failure_solution") or result["last_solution"]
    strain = np.asarray(fem_details._mechanism_field(
        fem_data, result["last_solution"], fail))
    cent = np.asarray(fem_details._element_centroids(fem_data))
    hot = strain >= 0.5 * np.nanmax(strain)
    w = strain[hot]
    top = int(np.argmax(cent[hot, 1]))              # where the band reaches the crest
    print("   mechanism   %d of %d elements above half the peak shear strain "
          "(%.4f); strain-weighted centroid (%.2f, %.2f); x %.2f–%.2f, "
          "y %.2f–%.2f; highest band element at (%.2f, %.2f)"
          % (int(hot.sum()), len(strain), float(np.nanmax(strain)),
             float((cent[hot, 0] * w).sum() / w.sum()),
             float((cent[hot, 1] * w).sum() / w.sum()),
             float(cent[hot, 0].min()), float(cent[hot, 0].max()),
             float(cent[hot, 1].min()), float(cent[hot, 1].max()),
             float(cent[hot][top, 0]), float(cent[hot][top, 1])))
    print("   LEM circle  entry (%.3f, %.3f) exit (%.3f, %.3f) · center "
          "(%.3f, %.3f) R %.3f — the surface the slices found, for comparison"
          % (xs[0], ys[0], xs[-1], ys[-1], crit["Xo"], crit["Yo"],
             crit["Yo"] - crit["Depth"]))


#: The overburden-dependent pullout law entered on FEM-2's alternative run:
#: adhesion zero, and the interface angle FHWA-NHI-10-024's default pullout
#: friction factor gives for a geogrid in sand — δ = arctan(F*·α) with
#: F* = (2/3)·tan φ′ and the geogrid scale-effect correction α = 0.8.  At
#: φ′ = 37° that is arctan(0.404) = 22.0°.
FEM02_LAW_ADHESION = 0.0
FEM02_LAW_DELTA = 22.0
#: The line the law's envelope is drawn along.  Line 2 rather than line 5: it is
#: the line whose hardest-worked point moves from the buried end to the face when
#: the law replaces the constant ramps, so its profile shows both the curve and
#: what the curve costs.
FEM02_LAW_PROFILE_LINE = 2


def fem02_pullout_law():
    """The same model with pullout resistance read from the overburden instead of
    from a stated development length: Adhesion 0, Delta 22°.

    One elastic-perfectly-plastic strength-reduction run and one Spencer search
    under the law, against the constant-law search on the same file, plus the
    bar-force profile that shows the curved capacity envelope.
    """
    import copy as _copy
    import time

    from xslope.fem import build_fem_data
    from xslope.fileio import ensure_reinforce_pullout
    from xslope import fem_details
    from xslope.plot_fem_details import plot_reinforcement_detail

    done = load_slope_data(FEM02_DONE)
    start = load_slope_data(FEM02_START)
    _u = declared_unit_labels(done)

    # ---- the limit-equilibrium half: the same search under each law ---------- #
    def _search(model):
        with contextlib.redirect_stdout(io.StringIO()):
            fs_cache, _, _p, _c = circular_search(
                model, FEM02_METHOD, num_slices=FEM02_SLICES, diagnostic=False,
                **file_search_window(model))
        return fs_cache[0]

    law_start = _copy.deepcopy(start)
    for line in law_start["reinforcement_lines"]:
        line["adhesion"] = FEM02_LAW_ADHESION
        line["delta"] = FEM02_LAW_DELTA
    ensure_reinforce_pullout(law_start)

    for tag, model in (("constant", start), ("law", law_start)):
        crit = _search(model)
        xs, ys = zip(*list(crit["failure_surface"].coords))
        print("   LEM %-9s %s FS %.4f on Xo %.4f Yo %.4f R %.4f · entry "
              "(%.3f, %.3f) exit (%.3f, %.3f) · ΣP %.0f %s"
              % (tag, FEM02_METHOD, crit["FS"], crit["Xo"], crit["Yo"],
                 crit["Yo"] - crit["Depth"], xs[0], ys[0], xs[-1], ys[-1],
                 crit["slices"]["p"].sum(), _u["force_per_len"]))

    # ---- the finite element half: one run, blank Tres ------------------------ #
    # The mesh is built from the model as the page meshes it, and the law is
    # applied to the copy that goes into ``build_fem_data``. Meshing the law's
    # own model instead would put a node at every one of its 41 stored tension
    # points per line — 240 bar elements against 60 — and the comparison would
    # be between two meshes rather than between two capacity laws.
    mesh = _fem02_mesh(done)
    law_done = _copy.deepcopy(done)
    for line in law_done["reinforcement_lines"]:
        line["adhesion"] = FEM02_LAW_ADHESION
        line["delta"] = FEM02_LAW_DELTA
    ensure_reinforce_pullout(law_done)
    print("   mesh        %d nodes · %d elements · %d bar elements"
          % (len(mesh["nodes"]), len(mesh["elements"]),
             len(mesh.get("elements_1d", []))))
    t0 = time.time()
    fem_data, result, _log = _fem02_solve(law_done, mesh, None)
    print("   law run · Tres blank · %.1f s wall" % (time.time() - t0))
    _fem02_report("law-epp", fem_data, result)

    prof = fem_details.reinforcement_profile(
        fem_data, result["last_solution"], FEM02_LAW_PROFILE_LINE,
        slope_data=law_done, field_state="converged",
        failure_solution=result.get("failure_solution"))
    print("   law  line %d · max force %.1f %s · hardest-worked %.1f at s %.2f "
          "%s (utilization %.3f) · %s"
          % (FEM02_LAW_PROFILE_LINE, float(prof["T"].max()), _u["force_per_len"],
             prof["peak_T"], prof["peak_s"], _u["length"],
             prof["peak_utilization"], prof["status"]))
    capture("fem02_bar_profile_law.png", plot_reinforcement_detail, prof)

    # The mechanism under the law, drawn and measured against the stated-length
    # run on the same mesh: where the band crosses each line, and where it
    # reaches the face. A capacity law that is weakest at the face can move the
    # band there even when the factor of safety barely moves.
    from xslope.plot_fem import plot_fem_results
    capture("fem02_shear_strain_law.png", plot_fem_results, fem_data, result,
            plot_type="shear_strain", fs=result["FS"],
            failure_solution=result.get("failure_solution"),
            field_state="failure")
    import numpy as np
    fem_con, res_con, _log = _fem02_solve(done, mesh, None)
    print("   constant-law run for comparison · FS %.4f" % res_con["FS"])
    for tag, fd, res, model in (("constant", fem_con, res_con, done),
                                ("law", fem_data, result, law_done)):
        fail = res.get("failure_solution") or res["last_solution"]
        strain = np.asarray(fem_details._mechanism_field(fd, res["last_solution"], fail))
        cent = np.asarray(fem_details._element_centroids(fd))
        hot = strain >= 0.5 * np.nanmax(strain)
        w = strain[hot]
        print("   %-8s mechanism · %d elements above half the peak (%.3f) · "
              "centroid (%.2f, %.2f) · x %.1f–%.1f, y %.1f–%.1f"
              % (tag, int(hot.sum()), float(np.nanmax(strain)),
                 float((cent[hot, 0] * w).sum() / w.sum()),
                 float((cent[hot, 1] * w).sum() / w.sum()),
                 cent[hot, 0].min(), cent[hot, 0].max(),
                 cent[hot, 1].min(), cent[hot, 1].max()))
        for ln in range(1, len(model["reinforcement_lines"]) + 1):
            pr = fem_details.reinforcement_profile(
                fd, res["last_solution"], ln, slope_data=model,
                field_state="converged", failure_solution=res.get("failure_solution"))
            print("        line %d · band %s–%s · hardest-worked at s %.2f (%.3f) · %s"
                  % (ln, pr.get("band_lo"), pr.get("band_hi"), pr["peak_s"],
                     pr["peak_utilization"], pr["status"]))

    # What the two laws allow along that line, so the page's claim about WHERE the
    # bond is critical is a measurement rather than a reading of the figure.
    from xslope.fileio import reinforce_available_tension
    r_law = law_done["reinforcement_lines"][FEM02_LAW_PROFILE_LINE - 1]
    r_con = done["reinforcement_lines"][FEM02_LAW_PROFILE_LINE - 1]
    L = math.hypot(r_con["x2"] - r_con["x1"], r_con["y2"] - r_con["y1"])
    print("   envelope along line %d (s from the face end, %s):"
          % (FEM02_LAW_PROFILE_LINE, _u["length"]))
    for s in (1.0, 2.0, 3.0, 4.0, 5.0, 10.0, L - 3.0, L - 2.0, L - 1.0):
        con = reinforce_available_tension(s, L - s, r_con["t_max"], r_con["lp1"],
                                          r_con["lp2"], r_con.get("tend1", 0.0),
                                          r_con.get("tend2", 0.0))
        law = reinforce_available_tension(s, L - s, r_law["t_max"], r_law["lp1"],
                                          r_law["lp2"], r_law.get("tend1", 0.0),
                                          r_law.get("tend2", 0.0),
                                          pullout=r_law.get("_pullout_profile"))
        print("        s %5.1f   constant %6.1f   law %6.1f" % (s, con, law))


def fem02_tres_sweep():
    """FS against the residual capacity typed into Tres, with the
    limit-equilibrium answer drawn across it.

    Its own group because it is five strength-reduction runs — about ten
    minutes — and nothing else on the page needs them.
    """
    import numpy as np

    done = load_slope_data(FEM02_DONE)
    mesh = _fem02_mesh(done)
    start = load_slope_data(FEM02_START)
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _, _p, _c = circular_search(
            start, FEM02_METHOD, num_slices=FEM02_SLICES, diagnostic=False,
            **file_search_window(start))
    lem_fs = fs_cache[0]["FS"]

    values, factors = [], []
    for t_res in FEM02_TRES_SWEEP:
        fem_data, result, _log = _fem02_solve(done, mesh, t_res)
        values.append(t_res)
        factors.append(result["FS"])
        print("   Tres %-6s FS %.4f  interval [%.6f, %.6f]"
              % ("blank" if t_res is None else "%g" % t_res, result["FS"],
                 result["final_interval"][0], result["final_interval"][1]))

    def _draw():
        _u = declared_unit_labels(done)
        numeric = sorted((v, f) for v, f in zip(values, factors) if v is not None)
        blank = [f for v, f in zip(values, factors) if v is None]
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        # The measured answer lands on three tiers, not on a curve: sloped
        # connectors between the points would draw a gradient the runs do not
        # have.  Steps hold each measured factor of safety flat between its
        # neighbors and put the change where the answer actually changes.
        ax.plot([v for v, _ in numeric], [f for _, f in numeric],
                drawstyle="steps-mid", marker="o", color="#1f4e79",
                linewidth=1.8, label="Tres entered")
        if blank:
            ax.axhline(blank[0], color="#7a5195", linestyle="--", linewidth=1.5,
                       label="Tres blank (elastic-perfectly-plastic)")
        ax.axhline(lem_fs, color="#c0392b", linestyle=":", linewidth=1.5,
                   label="Spencer (limit equilibrium)")
        ax.set_xlabel("Residual capacity Tres (%s)" % _u["force_per_len"])
        ax.set_ylabel("Factor of safety")
        ax.grid(True, alpha=0.3)
        # The measured curve steps twice and leaves a clear corner, so the legend
        # goes inside where matplotlib finds room for it rather than eating a
        # third of the figure width in a reserved column beside the axes.
        ax.legend(loc="best", fontsize=9, framealpha=0.9)
        fig.tight_layout()

    capture("fem02_tres_sweep.png", _draw)


# --------------------------------------------------------------------------- #
# FEM-3 — Piles: LEM against FEM
# --------------------------------------------------------------------------- #
#: FEM-3's two models, both on the same slope.  The discrete row is the pile
#: sample problem — one model, two sample pages: docs/lem/samples.md problem 10
#: locks Spencer's search on it and docs/fem/samples.md problem 2 locks a
#: strength-reduction run on it.  The sample's locked run is made over its own
#: bracket on the mesh committed beside its finite element copy; the tutorial
#: runs below are made at the settings the page has the reader enter, so they
#: build their own mesh from the limit-equilibrium copy the page links.
FEM03_PILES = os.path.join(REPO_ROOT, "docs/lem/files/xslope_piles.xlsx")
#: The continuous wall: the pair written by ``tools/build_pile_wall_tutorial.py``
#: from that same slope — the starter with no structural line in it, and the same
#: file with one PZ-27 sheet pile row.  Neither carries a mesh; the page builds
#: one, and so does every run here.
FEM03_WALL_START = os.path.join(REPO_ROOT,
                                "docs/tutorials/files/xslope_pile_wall_start.xlsx")
FEM03_WALL_DONE = os.path.join(REPO_ROOT,
                               "docs/tutorials/files/xslope_pile_wall.xlsx")
#: The limit-equilibrium method and slice count the piles model's own locked
#: search runs at (docs/lem/samples.md's circular_search tag).
FEM03_METHOD = "spencer"
FEM03_SLICES = 40
#: The mesh every figure on the page is drawn on: quadratic triangles at 2 ft,
#: with the thin-zone refinement off, which is what the page has the reader set
#: in Build Mesh.
FEM03_ELEMENT_TYPE = "tri6"
FEM03_TARGET_SIZE = 2.0
#: The optional refinement step of the wall half: a 0.5 ft 1D element size, which
#: also refines the soil the beam is embedded in.
FEM03_REFINED_1D = 0.5
#: The strength-reduction settings the page runs every trial at, which are
#: Studio's Run FEM dialog as it opens except for the bracket: ``non_convergence``
#: rather than the API default ``hybrid``, a 12,000-iteration budget, and a
#: 1.0–2.0 bracket, which has to reach 2.0 because the socketed runs stand at 1.6
#: and above.
FEM03_CRITERION = "non_convergence"
FEM03_F_MIN, FEM03_F_MAX = 1.0, 2.0
FEM03_TOLERANCE = 0.01
FEM03_MAX_ITERATIONS = 12000
#: The spacings the sweep runs, spanning S/D = 1.5 to 6 on the model's 2 ft
#: shafts.  6 ft is the file's own value and the locked case; 3 and 12 ft halve
#: and double it, which quarters and quadruples what the finite element model
#: sees, since spacing reaches it only as the divisor on EA and EI.
FEM03_SPACINGS = (3.0, 6.0, 12.0)


def _fem03_search(model):
    """One Spencer search on the piles model at the sample's own slice count."""
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _conv, _path, _cache = circular_search(
            model, FEM03_METHOD, num_slices=FEM03_SLICES, diagnostic=False,
            **file_search_window(model))
    return fs_cache[0]


def _fem03_reading(crit):
    """The one-line reading of a search: factor of safety, circle, and the pile
    force that reached the slices."""
    piles = (crit["slices"].attrs.get("pile_report") or []
             if crit.get("slices") is not None else [])
    total = sum(r["H_width"] for r in piles)
    return ("FS %.4f on Xo %.3f Yo %.3f R %.3f (deepest y %.3f) · ΣH %.1f "
            "per unit width from %d row(s)"
            % (crit["FS"], crit["Xo"], crit["Yo"], crit["Yo"] - crit["Depth"],
               crit["Depth"], total, len(piles)))


def _fem03_piles_model(spacing=None, head=None, tip=None, pile_E=None,
                       caps=True, path=FEM03_PILES):
    """The piles model with one thing changed on both rows.

    Every argument left ``None`` leaves the file's own cell alone, so a call with
    no arguments is the model as it is shipped: 2 ft shafts at 6 ft centers, both
    ends free to rotate, 46,000 lb of shear capacity and 60,000 lb·ft of moment
    capacity each.
    """
    sd = load_slope_data(path)
    for pile in sd["pile_lines"]:
        if spacing is not None:
            pile["S"] = float(spacing)
        if head is not None:
            pile["head_fixity"] = head
        if tip is not None:
            pile["tip_fixity"] = tip
        if pile_E is not None:
            pile["E"] = float(pile_E)
        if not caps:
            pile["V_cap"] = None
            pile["M_cap"] = None
    return sd


def _fem03_mesh(model, element_size_1d=None):
    """The mesh Studio's Build Mesh dialog builds for this model at the page's
    settings: quadratic triangles at 2 ft with the pile or wall lines carried in
    as constraint lines, so the beam elements lie on element edges and share
    their nodes with the soil on both sides.
    """
    from xslope.mesh import (build_mesh_from_polygons,
                             extract_constraint_line_geometry,
                             extract_point_constraints,
                             extract_size_regions, get_material_polygons)

    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(model)
    with contextlib.redirect_stdout(io.StringIO()):
        return build_mesh_from_polygons(
            get_material_polygons(model, reinf_lines=lines),
            FEM03_TARGET_SIZE, FEM03_ELEMENT_TYPE, lines=lines or None,
            element_size_1d=element_size_1d,
            point_constraints=extract_point_constraints(model),
            size_regions=extract_size_regions(model))


def _fem03_solve(model, mesh):
    """One strength-reduction run at the page's settings.  Returns
    ``(fem_data, result, solution, seconds)``, where ``solution`` is the last
    converged field with the captured mechanism attached — the pair the results
    plots and the detail profiles are both read from."""
    import time

    from xslope.fem import build_fem_data, solve_ssrm

    fem_data = build_fem_data(model, mesh)
    t0 = time.time()
    with contextlib.redirect_stdout(io.StringIO()):
        result = solve_ssrm(fem_data, F_min=FEM03_F_MIN, F_max=FEM03_F_MAX,
                            tolerance=FEM03_TOLERANCE, debug_level=0,
                            capture_failure_state=True,
                            failure_criterion=FEM03_CRITERION,
                            max_iterations=FEM03_MAX_ITERATIONS)
    seconds = time.time() - t0
    solution = dict(result["last_solution"])
    solution["failure_solution"] = result.get("failure_solution")
    return fem_data, result, solution, seconds


def _fem03_mechanism(fem_data, solution):
    """Where the shear strain concentrates at the captured mechanism: the peak,
    how many elements stand above half of it, and where they are."""
    import numpy as np

    from xslope import fem_details

    fail = solution.get("failure_solution") or solution
    strain = np.asarray(fem_details._mechanism_field(fem_data, solution, fail))
    cent = np.asarray(fem_details._element_centroids(fem_data))
    peak = float(np.nanmax(strain))
    hot = strain >= 0.5 * peak
    w = strain[hot]
    j = int(np.nanargmax(strain))
    return ("peak shear strain %.3f at element centroid (%.1f, %.1f) · %d of %d "
            "elements above half the peak · strain-weighted centroid (%.2f, %.2f) "
            "· x %.2f–%.2f, y %.2f–%.2f"
            % (peak, cent[j, 0], cent[j, 1], int(hot.sum()), len(strain),
               float((cent[hot, 0] * w).sum() / w.sum()),
               float((cent[hot, 1] * w).sum() / w.sum()),
               float(cent[hot, 0].min()), float(cent[hot, 0].max()),
               float(cent[hot, 1].min()), float(cent[hot, 1].max())))


def _fem03_beam(field, spacing, m_cap, v_cap):
    """What the beam elements carry in one field: the peaks per unit width of
    section, what those are per member, and how many elements have yielded."""
    import numpy as np

    mom = np.asarray(field.get("forces_pile_moment", []), dtype=float)
    shear = np.asarray(field.get("forces_pile_lateral", []), dtype=float)
    y_m = np.asarray(field.get("yielded_pile_M", []), dtype=bool)
    y_v = np.asarray(field.get("yielded_pile_V", []), dtype=bool)
    m_peak = float(np.abs(mom).max()) if mom.size else float("nan")
    v_peak = float(np.abs(shear).max()) if shear.size else float("nan")
    return ("peak moment %.0f per unit width (%.0f per member, %s) · "
            "peak shear %.0f per unit width (%.0f per member, %s) · "
            "%d element(s) yielded in bending, %d in shear"
            % (m_peak, m_peak * spacing,
               ("%.0f%% of Mcap" % (100.0 * m_peak * spacing / m_cap)) if m_cap
               else "Mcap blank",
               v_peak, v_peak * spacing,
               ("%.0f%% of Vcap" % (100.0 * v_peak * spacing / v_cap)) if v_cap
               else "Vcap blank",
               int(y_m.sum()) if y_m.size else 0,
               int(y_v.sum()) if y_v.size else 0))


def _fem03_report(label, model, mesh, fem_data, result, solution, seconds,
                  beam=True):
    """Everything the page can quote from one run."""
    row = model["pile_lines"][0] if model["pile_lines"] else None
    spacing = float(row["S"]) if row else 1.0
    print("   %-22s FS %.4f from [%.6f, %.6f] in %d trials · %d nodes, %d "
          "elements, %d beam elements · %.0f s"
          % (label, result["FS"], result["final_interval"][0],
             result["final_interval"][1], len(result["trials"]),
             len(mesh["nodes"]), len(mesh["elements"]),
             len(mesh.get("elements_1d", [])), seconds))
    if not (beam and row):
        return
    fail = solution.get("failure_solution") or solution
    print("      converged  %s"
          % _fem03_beam(solution, spacing, row["M_cap"], row["V_cap"]))
    print("      at failure %s"
          % _fem03_beam(fail, spacing, row["M_cap"], row["V_cap"]))
    print("      mechanism  %s" % _fem03_mechanism(fem_data, solution))


def _fem03_profiles(model, fem_data, solution, state="converged"):
    """The per-member profile Studio's 1D Details panel draws, read as numbers."""
    from xslope import fem_details

    out = []
    fail = solution.get("failure_solution")
    for i, row in enumerate(model["pile_lines"]):
        prof = fem_details.pile_profile(fem_data, solution, i, slope_data=model,
                                        field_state=state, failure_solution=fail)
        head_y = row["y1"]
        print("      %-16s %d beam elements · peak moment %.0f per unit width at "
              "depth %.2f (el. %.2f) · moment %.0f at the head and %.0f at the "
              "toe · peak shear %.0f at depth %.2f · head deflection %.4f · %s"
              % (row["label"], prof["n_elements"], prof["max_moment"],
                 prof["max_moment_depth"], head_y - prof["max_moment_depth"],
                 float(prof["moment"][0]), float(prof["moment"][-1]),
                 prof["max_shear"], prof["max_shear_depth"],
                 float(prof["u_lateral"][0]), prof["status"]))
        out.append(prof)
    return out


def fem03_piles():
    """The discrete row: the section both engines are given, the limit
    equilibrium answer with the two rows crossing the critical circle, the mesh
    the page builds, and the strength-reduction mechanism the finite element
    model produces on the file as it is shipped — both rows free to rotate at
    head and tip.
    """
    from xslope.plot_fem import plot_fem_data, plot_fem_results

    sd = _fem03_piles_model()
    _u = declared_unit_labels(sd)
    mat = sd["materials"][0]
    print("   model       %s: γ %g %s · c %g %s · φ %g° · E %g %s · ν %g · u %s"
          % (mat["name"], mat["gamma"], _u["unit_weight"], mat["c"], _u["stress"],
             mat["phi"], mat["E"], _u["stress"], mat["nu"], mat["u"]))
    for pile in sd["pile_lines"]:
        print("   pile row    %-9s (%g, %g) to (%g, %g) · D %g %s · S %g %s · "
              "E %g %s · Vcap %g · Mcap %g · head %s · tip %s · %s"
              % (pile["label"], pile["x1"], pile["y1"], pile["x2"], pile["y2"],
                 pile["D_pile"], _u["length"], pile["S"], _u["length"],
                 pile["E"], _u["stress"], pile["V_cap"], pile["M_cap"],
                 pile["head_fixity"], pile["tip_fixity"], pile["appl"]))

    capture("fem03_inputs_piles.png", plot_inputs, sd,
            title="Slope Geometry and Inputs")

    crit = _fem03_search(sd)
    capture("fem03_lem_solution_piles.png", plot_solution, sd, crit["slices"],
            crit["failure_surface"], crit["solver_result"])
    print("   LEM         %s" % _fem03_reading(crit))
    for rec in (crit["slices"].attrs.get("pile_report") or []):
        print("      row %-9s crosses at y %.3f · Ito & Matsui %.0f per pile · "
              "%s → %.0f per pile = %.1f %s applied"
              % (rec["label"], rec["y"], rec["F_soil"], rec["governed"],
                 rec["F_used"], rec["H_width"], _u["force_per_len"]))
    xs, ys = zip(*list(crit["failure_surface"].coords))
    print("   LEM circle  entry (%.3f, %.3f) exit (%.3f, %.3f)"
          % (xs[0], ys[0], xs[-1], ys[-1]))

    mesh = _fem03_mesh(sd)
    fem_data, result, solution, seconds = _fem03_solve(sd, mesh)
    capture("fem03_mesh_piles.png", plot_fem_data, fem_data)
    capture("fem03_fem_shear_piles.png", plot_fem_results, fem_data, solution,
            plot_type="shear_strain", fs=result["FS"],
            failure_solution=solution.get("failure_solution"),
            field_state="failure")
    _fem03_report("FEM as shipped", sd, mesh, fem_data, result, solution, seconds)
    _fem03_profiles(sd, fem_data, solution)


def fem03_tip():
    """What the flat spacing line is actually measuring: the shafts' toe.

    Five runs on the pile model, each changing one thing from the file as it is
    shipped, and all of them at the page's own settings on the page's own mesh.
    The first four are the diagnostic — stiffness up, stiffness down, the
    structural capacities cleared, the heads restrained — and the fifth is the
    answer they point at: the shafts end on the rigid base with a free toe, so
    they swing about it, and fixing that toe is what moves the factor of safety.
    The tip-fixed run also draws the mechanism figure, which is a different
    failure from the one the shipped file produces.
    """
    from xslope.plot_fem import plot_fem_results

    E0 = load_slope_data(FEM03_PILES)["pile_lines"][0]["E"]
    cases = [
        ("pile E ×100", dict(pile_E=E0 * 100.0), None),
        ("pile E ÷100", dict(pile_E=E0 / 100.0), None),
        ("Vcap and Mcap cleared", dict(caps=False), None),
        ("heads fixed, tips free", dict(head="fixed"), None),
        ("tips fixed", dict(tip="fixed"), "fem03_fem_shear_piles_fixed.png"),
    ]
    for label, kwargs, figure in cases:
        sd = _fem03_piles_model(**kwargs)
        mesh = _fem03_mesh(sd)
        fem_data, result, solution, seconds = _fem03_solve(sd, mesh)
        if figure:
            capture(figure, plot_fem_results, fem_data, solution,
                    plot_type="shear_strain", fs=result["FS"],
                    failure_solution=solution.get("failure_solution"),
                    field_state="failure")
        _fem03_report(label, sd, mesh, fem_data, result, solution, seconds)
        _fem03_profiles(sd, fem_data, solution)


def fem03_spacing():
    """The page's falsifiable test: the same two pile rows at 3, 6 and 12 ft
    spacing, put through both engines, and through the finite element engine
    twice — once with the shafts' toes free, as the file ships, and once with
    them fixed.

    The limit equilibrium runs are full searches, because spacing changes which
    surface governs and a held circle would hide half the effect.  The
    strength-reduction runs change nothing but the S cell, which reaches the
    model only as the divisor on EA and EI.
    """
    sd = load_slope_data(FEM03_PILES)
    _u = declared_unit_labels(sd)
    row = sd["pile_lines"][0]
    D, E = row["D_pile"], row["E"]
    area, inertia = math.pi * D ** 2 / 4.0, math.pi * D ** 4 / 64.0
    print("   section     D %g %s → A %.4f %s² · I %.4f %s⁴ · E %g %s"
          % (D, _u["length"], area, _u["length"], inertia, _u["length"], E,
             _u["stress"]))

    lem, free, fixed = {}, {}, {}
    for spacing in FEM03_SPACINGS:
        model = _fem03_piles_model(spacing=spacing)
        crit = _fem03_search(model)
        lem[spacing] = crit["FS"]
        print("   LEM  S %-5g %s" % (spacing, _fem03_reading(crit)))

    for tip, store in (("free", free), ("fixed", fixed)):
        for spacing in FEM03_SPACINGS:
            model = _fem03_piles_model(spacing=spacing, tip=tip)
            mesh = _fem03_mesh(model)
            fem_data, result, solution, seconds = _fem03_solve(model, mesh)
            store[spacing] = result["FS"]
            print("   FEM  tips %s, S %-5g EA/S %.4g · EI/S %.4g"
                  % (tip, spacing, E * area / spacing, E * inertia / spacing))
            _fem03_report("     ", model, mesh, fem_data, result, solution,
                          seconds)

    _fem03_spacing_figure(lem, fixed, free, _u)


def _fem03_spacing_figure(lem, fixed, free, _u):
    """The sweep on one axis: the limit equilibrium curve and the two
    strength-reduction lines, each point labeled with its own answer.

    The labels sit above the limit equilibrium curve and below both flat lines,
    which is what keeps the 6 ft column readable: the limit equilibrium point and
    the tip-fixed point are 0.05 apart there, and a label above each would print
    one over the other.
    """
    def _draw():
        fig, ax = plt.subplots(figsize=(7.0, 4.6))
        xs = list(FEM03_SPACINGS)
        series = [
            (lem, "LEM (Spencer, Ito & Matsui)", "#1f77b4", "o-", 9),
            (fixed, "FEM (SSRM, pile tips fixed)", "#9467bd", "^-", -16),
            (free, "FEM (SSRM, pile tips free)", "#d62728", "s-", -16),
        ]
        for data, label, color, style, dy in series:
            ax.plot(xs, [data[s] for s in xs], style, color=color,
                    linewidth=2.0, markersize=7, label=label)
            for s in xs:
                ax.annotate("%.3f" % data[s], (s, data[s]),
                            textcoords="offset points", xytext=(0, dy),
                            ha="center", fontsize=9, color=color)
        ax.set_xlabel("Pile spacing S (%s)" % _u["length"])
        ax.set_ylabel("Factor of safety")
        ax.set_title("Two rows of 2 ft shafts: what each engine does with spacing")
        ax.set_xticks(xs)
        ax.set_xlim(min(xs) - 1.0, max(xs) + 1.0)
        # Room for the annotation rows, measured off the data rather than padded
        # by a guess: the lower labels sit under the two flat lines and would
        # otherwise print on the axis.
        lo = min(min(d.values()) for d, *_ in series)
        hi = max(max(d.values()) for d, *_ in series)
        ax.set_ylim(lo - 0.12 * (hi - lo), hi + 0.10 * (hi - lo))
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", frameon=False, fontsize=9)
        fig.tight_layout()

    capture("fem03_spacing_sweep.png", _draw)


def fem03_wall():
    """The continuous wall: the slope the reader starts from, the mesh the wall
    line forces, the two mechanisms its toe condition produces, and the internal
    actions only the finite element engine can give.

    Four runs, all at the page's settings: the starter with no wall in it, the
    wall with its tip free, the wall with its tip fixed, and the tip-fixed wall
    again on a 0.5 ft 1D element size.  A fifth run clears Mcap on the tip-free
    wall, which is the page's check that the cap never binds there.
    """
    from xslope.plot_fem import plot_fem_data, plot_fem_results
    from xslope.plot_fem_details import plot_pile_detail

    start = load_slope_data(FEM03_WALL_START)
    _u = declared_unit_labels(start)
    mat = start["materials"][0]
    print("   model       %s: γ %g %s · c %g %s · φ %g° · E %g %s · ν %g"
          % (mat["name"], mat["gamma"], _u["unit_weight"], mat["c"], _u["stress"],
             mat["phi"], mat["E"], _u["stress"], mat["nu"]))
    print("   starter     %d pile row(s) · %s at target size %g %s · SSRM "
          "bracket [%g, %g]"
          % (len(start["pile_lines"]), start["element_type"],
             start["target_size"], _u["length"], start["ssrm_f_min"],
             start["ssrm_f_max"]))
    wall = load_slope_data(FEM03_WALL_DONE)["pile_lines"][0]
    print("   wall row    %s: (%g, %g) to (%g, %g) · E %g %s · Area %g %s²/%s · "
          "I %g %s⁴/%s · S %g · D %r · Vcap %r · Mcap %g · head %s · tip %s"
          % (wall["label"], wall["x1"], wall["y1"], wall["x2"], wall["y2"],
             wall["E"], _u["stress"], wall["area"], _u["length"], _u["length"],
             wall["I"], _u["length"], _u["length"], wall["S"], wall["D_pile"],
             wall["V_cap"], wall["M_cap"], wall["head_fixity"],
             wall["tip_fixity"]))
    print("   wall EA %.4g · EI %.4g" % (wall["E"] * wall["area"],
                                         wall["E"] * wall["I"]))

    capture("fem03_inputs_wall.png", plot_inputs, start,
            title="Slope Geometry and Inputs")

    # The bare slope, both ways — the baseline both halves of the page measure
    # against, drawn so the reader sees the two mechanisms before any member.
    crit = _fem03_search(start)
    capture("fem03_lem_solution_bare.png", plot_solution, start, crit["slices"],
            crit["failure_surface"], crit["solver_result"])
    print("   LEM (bare)  %s" % _fem03_reading(crit))
    mesh = _fem03_mesh(start)
    fem_data, result, solution, seconds = _fem03_solve(start, mesh)
    _fem03_report("no wall", start, mesh, fem_data, result, solution, seconds)
    capture("fem03_fem_shear_bare.png", plot_fem_results, fem_data, solution,
            plot_type="shear_strain", fs=result["FS"],
            failure_solution=result.get("failure_solution"),
            field_state="failure")

    #: (label, tip, 1D element size, keep Mcap, mesh figure, shear-strain figure,
    #: profile figure)
    runs = [
        ("wall, tip free", "free", None, True, "fem03_mesh_wall.png",
         "fem03_wall_shear.png", "fem03_wall_profiles.png"),
        ("wall, tip fixed", "fixed", None, True, None,
         "fem03_wall_shear_fixed.png", "fem03_wall_profiles_fixed.png"),
        ("wall, tip fixed, 1D 0.5", "fixed", FEM03_REFINED_1D, True, None, None,
         "fem03_wall_profiles_refined.png"),
        ("wall, tip free, Mcap blank", "free", None, False, None, None, None),
    ]
    for label, tip, size_1d, keep_cap, mesh_fig, shear_fig, prof_fig in runs:
        sd = load_slope_data(FEM03_WALL_DONE)
        sd["pile_lines"][0]["tip_fixity"] = tip
        if not keep_cap:
            sd["pile_lines"][0]["M_cap"] = None
        mesh = _fem03_mesh(sd, element_size_1d=size_1d)
        fem_data, result, solution, seconds = _fem03_solve(sd, mesh)
        if mesh_fig:
            capture(mesh_fig, plot_fem_data, fem_data)
        if shear_fig:
            capture(shear_fig, plot_fem_results, fem_data, solution,
                    plot_type="shear_strain", fs=result["FS"],
                    failure_solution=solution.get("failure_solution"),
                    field_state="failure")
        _fem03_report(label, sd, mesh, fem_data, result, solution, seconds)
        profs = _fem03_profiles(sd, fem_data, solution)
        if prof_fig:
            capture(prof_fig, plot_pile_detail, profs[0])


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
    "seep03_plots": seep03_plots,
    "seep04_plots": seep04_plots,
    "fem01_plots": fem01_plots,
    "fem02_plots": fem02_plots,
    "fem02_pullout_law": fem02_pullout_law,
    "fem02_tres_sweep": fem02_tres_sweep,
    "fem03_piles": fem03_piles,
    "fem03_tip": fem03_tip,
    "fem03_spacing": fem03_spacing,
    "fem03_wall": fem03_wall,
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
