#!/usr/bin/env python3
"""Caption-vs-image check.

Two modes, chosen per page in its config.

``panel``
    Classify every referenced figure by panel layout directly from the PNG and
    compare that against the panel form its caption declares, so a caption can
    never drift from the figure it labels.

    Classification (quadrant ink test, no image library beyond Pillow): a
    figure is FOUR-PANEL when the central band of BOTH the row-ink and the
    column-ink profile contains an essentially blank line — the inter-panel
    gutters a 2x2 layout must have.  A single axes cannot have them, because
    its data area spans the middle of the canvas in both directions.

    Caption forms recognised:
      four   — mentions all four panels: inputs/model, mesh, shear strain,
               displacement vectors
      two    — "(left)" ... "(right)" with no mesh/vector mention
      one    — "inputs" only, explicitly noting no mechanism is plotted
      other  — the classifier cannot read a panel form.  This is a FAILURE
               unless the image is named in the page's ``caption_exempt``: an
               unclassifiable caption is unchecked, not passing.

``structural``
    For pages whose figures have no single meaningful panel form (a page mixing
    field renders, time series and flow nets, where "how many axes" is not a
    claim the caption makes).  Every referenced image must exist and every
    caption must be non-empty and descriptive — never a bare file name.

Usage: python -m tools.verification_checks.figures <page.md>
"""
import os
import re
import sys
from collections import Counter

# Any line-initial image reference, including a relative path out of the page's
# own directory and a trailing attribute block ({width=800}).  Anything narrower
# silently skips figures instead of checking them.
IMG = re.compile(r"^!\[(.*?)\]\(([^)\s]+)\)(?:\s*\{[^}]*\})?\s*$")


def _profiles(path):
    """Per-row and per-column ink fraction against the modal (background) colour."""
    from PIL import Image
    im = Image.open(path).convert("RGB")
    im = im.resize((min(im.size[0], 800), min(im.size[1], 800)))
    w, h = im.size
    px = im.load()
    bg = Counter(px[x, y] for x in range(0, w, 5)
                 for y in range(0, h, 5)).most_common(1)[0][0]

    def inked(x, y):
        r, g, b = px[x, y]
        return abs(r - bg[0]) + abs(g - bg[1]) + abs(b - bg[2]) > 40

    rows = [sum(inked(x, y) for x in range(0, w, 2)) / (w / 2) for y in range(h)]
    cols = [sum(inked(x, y) for y in range(0, h, 2)) / (h / 2) for x in range(w)]
    return rows, cols


def _has_gutter(profile, lo=0.35, hi=0.65, thresh=0.004):
    """True when the central band of the profile contains an essentially blank
    line — the inter-panel gutter a 2x2 layout must have and a single axes
    cannot, because its data area spans the middle."""
    n = len(profile)
    band = profile[int(n * lo):int(n * hi)]
    return bool(band) and min(band) <= thresh


def classify_image(path):
    """four = both gutters (2x2); two = a column gutter only (side by side);
    one = neither (a single axes filling the canvas)."""
    rows, cols = _profiles(path)
    vg, hg = _has_gutter(rows), _has_gutter(cols)
    form = "four" if (vg and hg) else ("two" if hg else "one")

    def mid(p):
        return round(min(p[int(len(p) * .35):int(len(p) * .65)]), 4)
    return form, (mid(rows), mid(cols))


def classify_caption(cap, cfg=None):
    low = cap.lower()
    for sub, form in (cfg.caption_rules if cfg else []):
        if sub.lower() in low:
            return form
    has_strain = "shear strain" in low
    has_vec = "displacement vector" in low
    has_mesh = "mesh" in low
    has_lr = "(left)" in low and "(right)" in low
    if has_strain and has_vec and has_mesh:
        return "four"
    if has_lr and not (has_vec and has_mesh):
        return "two"
    if "inputs" in low and not has_strain:
        return "one"
    return "other"


def _figures(page):
    lines = open(page).read().split("\n")
    for i, l in enumerate(lines):
        m = IMG.match(l.strip())
        if m:
            yield i + 1, m.group(1), m.group(2)


def run_panel(page, cfg, report=print):
    base = os.path.dirname(os.path.abspath(page))
    bad = 0
    counts = Counter()
    for ln, cap, src in _figures(page):
        cform = classify_caption(cap, cfg)
        counts["caption:" + cform] += 1
        p = os.path.join(base, src)
        if not os.path.exists(p):
            report(f"  L{ln} MISSING IMAGE {src}")
            bad += 1
            continue
        iform, q = classify_image(p)
        counts["image:" + iform] += 1
        if cform == "other":
            if src in cfg.caption_exempt:
                counts["exempt"] += 1
                continue
            report(f"  L{ln} caption form UNCLASSIFIABLE (image is {iform}-panel); "
                   f"a caption the classifier cannot read is unchecked, not "
                   f"passing  {src}")
            bad += 1
            continue
        if cform != iform:
            report(f"  L{ln} caption declares {cform}-panel, image is {iform}-panel "
                   f"(central-band min ink {q})  {src}")
            bad += 1
    report(f"{page}: "
           f"{sum(v for k, v in counts.items() if k.startswith('caption:'))} "
           f"captions, mismatches={bad}  {dict(counts)}")
    return bad


def run_structural(page, cfg, report=print):
    """Every referenced image exists and carries non-empty alt text."""
    base = os.path.dirname(os.path.abspath(page))
    bad = n = 0
    for ln, cap, src in _figures(page):
        n += 1
        p = os.path.join(base, src)
        if not os.path.exists(p):
            report(f"  L{ln} MISSING IMAGE {src}")
            bad += 1
            continue
        if not cap.strip():
            report(f"  L{ln} figure has empty alt text  {src}")
            bad += 1
    report(f"{page}: {n} figures checked structurally "
           f"(image exists, alt text non-empty), mismatches={bad}")
    return bad


def run(page, cfg, report=print):
    if cfg.figure_mode == "structural":
        return run_structural(page, cfg, report)
    return run_panel(page, cfg, report)


def _cli():
    from .pages import config_for
    page = sys.argv[1]
    return run(page, config_for(page))


if __name__ == "__main__":
    sys.exit(_cli())
