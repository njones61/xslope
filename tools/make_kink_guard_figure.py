"""Render docs/lem/images/noncircular_kink_guard.png for the search page.

VP17, with the non-circular search started from eight free points on its
critical circle (the start test/noncircular_kink_check.py uses), searched twice
with Spencer's method: on the left with the bend-back limit turned off
(``max_kink_back=None``), on the right with the default limit. The left panel is
the false minimum the limit exists to refuse.

Run from the repo root:
    python tools/make_kink_guard_figure.py
"""
import contextlib
import io
import os
import sys
import warnings

warnings.filterwarnings("ignore")
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

from xslope.fileio import load_slope_data  # noqa: E402
from xslope.plot import plot_solution  # noqa: E402
from xslope.search import MAX_KINK_BACK, noncircular_search  # noqa: E402

sys.path.insert(0, os.path.join(ROOT, "test"))
from noncircular_kink_check import NUM_SLICES, START, VP17  # noqa: E402

OUT = os.path.join(ROOT, "docs", "lem", "images", "noncircular_kink_guard.png")


def panel(max_kink_back, path, panel_size=(8.0, 5.0), dpi=150):
    sd = load_slope_data(VP17)
    sd["non_circ"] = [{"X": x, "Y": y, "Movement": "Free"} for x, y in START]
    with contextlib.redirect_stdout(io.StringIO()):
        cache, _conv, _path = noncircular_search(sd, "spencer", num_slices=NUM_SLICES,
                                                 max_kink_back=max_kink_back)
    best = cache[0]
    fig = plt.figure(figsize=panel_size)
    plot_solution(sd, best["slices"], best["failure_surface"], best["solver_result"],
                  fig=fig, show_title=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return best["FS"]


if __name__ == "__main__":
    parts = [OUT.replace(".png", "_off.png"), OUT.replace(".png", "_on.png")]
    f_off = panel(None, parts[0])
    f_on = panel(MAX_KINK_BACK, parts[1])
    imgs = [Image.open(p) for p in parts]
    h = min(im.height for im in imgs)
    imgs = [im.resize((round(im.width * h / im.height), h), Image.LANCZOS) for im in imgs]
    combo = Image.new("RGB", (sum(im.width for im in imgs) + 20, h), "white")
    x = 0
    for im in imgs:
        combo.paste(im, (x, 0))
        x += im.width + 20
    combo.save(OUT)
    for p in parts:
        os.remove(p)
    print(f"limit off F = {f_off:.3f}, limit on F = {f_on:.3f} -> {OUT}")
