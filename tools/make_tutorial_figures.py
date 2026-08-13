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
    """The template's ``main`` worksheet as it is downloaded, with the tab strip.

    The strip is the point of this one: it names every worksheet in the workbook,
    which is what the tutorial's paragraph beside it is describing. Same row/column
    window as LEM-1's filled ``main`` capture — content stops at the Water-loads
    row, below which the sheet holds the hidden source data behind the dropdowns.
    """
    render("t0_template_main.png", TEMPLATE, "main", rows=(1, 24), cols="A:D",
           tab_strip=True)


GROUPS = {
    "t0_template": t0_template,
    "lem01_sheets": lem01_sheets,
    "lem01_plots": lem01_plots,
    "lem01_placeholders": lem01_placeholders,
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
