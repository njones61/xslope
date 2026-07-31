"""Standing checks on what a double-click on the Studio Inputs canvas resolves to
when the click lands inside a material zone.

A polygon-based file has no profile sheet: its geometry IS the polygon sheet, so a
zone's interior is the most direct way to reach the row that defines it. The click
therefore opens the polygon editor on that zone. Getting this wrong is silent — a
wrong row opens a real editor on real geometry, and the user edits the wrong zone —
so the checks assert the RESOLVED ROW, not merely that something opened.

What is guarded here:

  A. PRECEDENCE — a line, point, or vertex within the pick tolerance still wins.
     A piezometric line drawn across a zone's interior, a failure circle whose arc
     crosses it, and a polygon boundary itself all beat the interior fallback; the
     interior is consulted ONLY when nothing is within tolerance. A click near a
     boundary shared by two zones resolves through the edge test, which is why it
     still answers from OUTSIDE the domain where no interior contains the point.
  B. THE INTERIOR — a click inside a zone, with nothing near, resolves to
     ('polygons', row) where row is that zone's position in slope_data['polygons'],
     which is the row the editor builds it at.
  C. OVERLAP — where zones nest or overlap, the SMALLEST containing zone wins. A
     zone wholly inside a larger one is otherwise unreachable by an interior click,
     and the innermost zone is what a click most specifically identifies.
  D. PROFILE FILES ARE UNCHANGED — a profile-based file's zones are derived from
     its profile lines and own no polygon-sheet row, so an interior click there
     still resolves to that zone's material, as it did before. A click inside no
     zone at all resolves to nothing, on either kind of file.
  E. THE WIRE — the canvas `picked` signal reaches MainWindow._on_canvas_pick, which
     hands the category and row to edit_category. edit_category is replaced by a
     recorder, so no modal dialog opens and the assertion is about the arguments.
  F. THE ROW LANDS — PolygonEditor.build(..., select=row) opens with that row
     current, and the row's vertices are the picked zone's. This is the check that
     fails if the polygon editor's merged list (material zones, then SSR overlays,
     then refine regions) ever stops putting the material zones first, which is the
     offset the resolved row assumes.

Legs A-D are file-less and Qt-free; E and F open real Studio widgets offscreen and
skip cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import io
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: A polygon-based section (no profile sheet) for the Studio legs.
POLY_FILE = os.path.join(_REPO, "docs/seep/files/xslope_levee_poly.xlsx")

#: The pick tolerance MplCanvas derives from a ~12 px offset, in data units. The
#: synthetic sections below are ~100 units wide, so 1.0 is a realistic radius.
TOL = 1.0


def _quiet(fn, *args, **kwargs):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        return fn(*args, **kwargs)


def _zones():
    """Two stacked material zones sharing the boundary y = 50, as a polygon sheet
    loads them: lower 0 <= y <= 50, upper 50 <= y <= 80, both 0 <= x <= 100."""
    from shapely.geometry import Polygon
    lower = Polygon([(0, 0), (100, 0), (100, 50), (0, 50)])
    upper = Polygon([(0, 50), (100, 50), (100, 80), (0, 80)])
    return [{"polygon": lower, "mat_id": 0}, {"polygon": upper, "mat_id": 1}]


def _poly_data(**extra):
    """A polygon-based slope_data: material zones and no profile sheet."""
    d = {"polygons": _zones(), "materials": [{"name": "lower"}, {"name": "upper"}]}
    d.update(extra)
    return d


# --------------------------------------------------------------- A. precedence
def test_precedence():
    """Anything within the tolerance still wins; the interior is the last resort."""
    fails = []
    from studio.picking import pick_category

    # A piezometric line drawn ACROSS the lower zone's interior. A click on it is a
    # click inside the zone too — the line has to win, or the piezo line becomes
    # unpickable everywhere it crosses a zone, which is everywhere it exists.
    d = _poly_data(piezo_line=[(0, 25), (100, 25)])
    hit = pick_category(d, 50.0, 25.3, TOL)
    if hit != ("piezo", None):
        fails.append(f"click on a piezo line inside a zone resolved to {hit!r}, "
                     "want ('piezo', None) — the interior fallback pre-empted a feature")
    # Same file, same zone, away from the line: now the interior answers.
    hit = pick_category(d, 50.0, 40.0, TOL)
    if hit != ("polygons", 0):
        fails.append(f"click inside the lower zone away from the piezo line resolved "
                     f"to {hit!r}, want ('polygons', 0)")

    # A failure circle whose arc crosses the same interior.
    d = _poly_data(circles=[{"Xo": 50.0, "Yo": 80.0, "R": 40.0}])
    hit = pick_category(d, 50.0, 40.2, TOL)          # arc passes through (50, 40)
    if hit != ("circles", 0):
        fails.append(f"click on a circle arc inside a zone resolved to {hit!r}, "
                     "want ('circles', 0)")

    # Reinforcement inside the zone.
    d = _poly_data(reinforcement_lines=[{"x1": 10.0, "y1": 30.0,
                                         "x2": 90.0, "y2": 30.0}])
    hit = pick_category(d, 50.0, 30.4, TOL)
    if hit != ("reinforce", 0):
        fails.append(f"click on reinforcement inside a zone resolved to {hit!r}, "
                     "want ('reinforce', 0)")

    # The boundary the two zones SHARE is a polygon edge, and the edge test — not
    # the interior fallback — is what answers there. Probed from just outside the
    # domain, where no zone contains the point at all: a non-None answer can only
    # have come from the edge branch.
    d = _poly_data()
    hit = pick_category(d, 100.4, 50.0, TOL)
    if hit is None or hit[0] != "polygons":
        fails.append(f"click within tolerance of a shared polygon edge, outside the "
                     f"domain, resolved to {hit!r}, want a ('polygons', row) edge hit")
    return fails


# ----------------------------------------------------------------- B. interior
def test_interior():
    """A click inside a zone, nothing near, opens the polygon editor on that zone."""
    fails = []
    from studio.picking import pick_category

    d = _poly_data()
    for (x, y), want in (((50.0, 25.0), 0), ((50.0, 65.0), 1),
                         ((5.0, 10.0), 0), ((95.0, 75.0), 1)):
        hit = pick_category(d, x, y, TOL)
        if hit != ("polygons", want):
            fails.append(f"click at ({x}, {y}) resolved to {hit!r}, "
                         f"want ('polygons', {want})")

    # The resolved row must genuinely be the zone under the click, not an accident
    # of list order: check containment for every zone from its own interior.
    for i, zone in enumerate(d["polygons"]):
        c = zone["polygon"].representative_point()
        hit = pick_category(d, c.x, c.y, 1e-9)
        if hit != ("polygons", i):
            fails.append(f"a point inside zone {i} resolved to {hit!r}, "
                         f"want ('polygons', {i})")
    return fails


# ------------------------------------------------------------------ C. overlap
def test_overlap():
    """Nested / overlapping zones: the smallest container wins."""
    fails = []
    from shapely.geometry import Polygon
    from studio.picking import pick_category

    big = Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
    small = Polygon([(40, 40), (60, 40), (60, 60), (40, 60)])
    # Listed big-first and small-first: the answer must be the small zone either
    # way, so it cannot be coming from list order.
    for order, want in (((big, small), 1), ((small, big), 0)):
        d = {"polygons": [{"polygon": p, "mat_id": i} for i, p in enumerate(order)],
             "materials": [{"name": "a"}, {"name": "b"}]}
        hit = pick_category(d, 50.0, 50.0, 1e-9)      # deep inside both
        if hit != ("polygons", want):
            fails.append(f"click inside nested zones (small at row {want}) resolved "
                         f"to {hit!r}, want ('polygons', {want})")

    # Partial overlap, not nesting: in the shared lens the smaller zone wins; in
    # the part of the big zone outside the lens the big zone answers.
    other = Polygon([(30, 30), (70, 30), (70, 70), (30, 70)])
    d = {"polygons": [{"polygon": big, "mat_id": 0}, {"polygon": other, "mat_id": 1}],
         "materials": [{"name": "a"}, {"name": "b"}]}
    hit = pick_category(d, 50.0, 50.0, 1e-9)
    if hit != ("polygons", 1):
        fails.append(f"click in the overlap of two zones resolved to {hit!r}, "
                     "want ('polygons', 1) — the smaller zone")
    hit = pick_category(d, 10.0, 10.0, 1e-9)
    if hit != ("polygons", 0):
        fails.append(f"click inside only the larger zone resolved to {hit!r}, "
                     "want ('polygons', 0)")
    return fails


# --------------------------------------------- D. profile files / nothing under
def test_profile_and_empty():
    """Profile-based files and clicks over empty space keep their old answers."""
    fails = []
    from studio.picking import pick_category

    # A profile-based file also carries `polygons` — the loader derives them from
    # the profile lines — but they own no polygon-sheet row, so an interior click
    # still edits the zone's material, exactly as before.
    zones = _zones()
    d = {"polygons": zones, "materials": [{"name": "lower"}, {"name": "upper"}],
         "profile_lines": [{"coords": [(0, 80), (100, 80)]}], "max_depth": 0.0}
    hit = pick_category(d, 50.0, 25.0, TOL)
    if hit != ("materials", 0):
        fails.append(f"interior click on a profile-based file resolved to {hit!r}, "
                     "want ('materials', 0) — unchanged behaviour")
    hit = pick_category(d, 50.0, 65.0, TOL)
    if hit != ("materials", 1):
        fails.append(f"interior click on the upper zone of a profile-based file "
                     f"resolved to {hit!r}, want ('materials', 1)")

    # Nothing under the click: still nothing, on either kind of file.
    for label, data in (("polygon-based", _poly_data()), ("profile-based", d)):
        hit = pick_category(data, 500.0, 500.0, TOL)
        if hit is not None:
            fails.append(f"click over empty space on a {label} file resolved to "
                         f"{hit!r}, want None")
    # An empty document resolves to nothing rather than raising.
    if pick_category({}, 0.0, 0.0, TOL) is not None:
        fails.append("click on an empty document did not resolve to None")
    return fails


# ---------------------------------------------------------------------- E. wire
def test_wire():
    """The canvas pick signal reaches edit_category with the category and row."""
    fails = []
    from studio.main_window import MainWindow

    mw = MainWindow()
    try:
        _quiet(mw.open_path, POLY_FILE)
        polys = mw.doc.slope_data.get("polygons") or []
        if not polys:
            return [f"{os.path.basename(POLY_FILE)} loaded with no polygons — "
                    "the wire leg has nothing to click inside"]
        if mw.doc.slope_data.get("profile_lines"):
            return [f"{os.path.basename(POLY_FILE)} is profile-based — "
                    "the wire leg needs a polygon-based file"]

        seen = []
        mw.edit_category = lambda category, select=None: seen.append((category, select))

        # Click deep inside each zone, through the canvas signal the double-click
        # emits, with a tolerance small enough that no boundary is in reach.
        for i, zone in enumerate(polys):
            geom = zone["polygon"]
            c = geom.representative_point()
            seen.clear()
            mw.canvas.picked.emit(c.x, c.y, 1e-9)
            if seen != [("polygons", i)]:
                fails.append(f"a canvas pick inside zone {i} of "
                             f"{os.path.basename(POLY_FILE)} reached edit_category as "
                             f"{seen!r}, want [('polygons', {i})]")

        # Far outside the section: nothing opens.
        xs = [x for z in polys for (x, _) in z["polygon"].exterior.coords]
        ys = [y for z in polys for (_, y) in z["polygon"].exterior.coords]
        seen.clear()
        mw.canvas.picked.emit(max(xs) + 1000.0, max(ys) + 1000.0, 1e-9)
        if seen:
            fails.append(f"a canvas pick far outside the section opened {seen!r}, "
                         "want nothing")
    finally:
        mw.close()
    return fails


# ------------------------------------------------------------------ F. the row
def test_row_lands():
    """The resolved row is the row the polygon editor opens on, with those vertices."""
    fails = []
    from studio.editors import CATEGORY_EDITORS
    from studio.main_window import MainWindow
    from studio.picking import pick_category

    editor = CATEGORY_EDITORS.get("polygons")
    if editor is None:
        return ["no 'polygons' category editor registered"]

    mw = MainWindow()
    try:
        _quiet(mw.open_path, POLY_FILE)
        data = mw.doc.slope_data
        polys = data.get("polygons") or []
        if not polys:
            return [f"{os.path.basename(POLY_FILE)} loaded with no polygons"]

        for i, zone in enumerate(polys):
            c = zone["polygon"].representative_point()
            hit = pick_category(data, c.x, c.y, 1e-9)
            if hit != ("polygons", i):
                fails.append(f"zone {i} of {os.path.basename(POLY_FILE)} resolved to "
                             f"{hit!r}, want ('polygons', {i})")
                continue
            dlg = editor.build(data, mw, select=hit[1])
            try:
                row = dlg.list.currentRow()
                if row != i:
                    fails.append(f"the polygon editor opened on row {row} for a click "
                                 f"inside zone {i}")
                    continue
                # ...and that row's vertices are this zone's (the merged list puts
                # material zones first; an overlay offset would show up here).
                want = list(zone["polygon"].exterior.coords)
                if len(want) >= 2 and want[0] == want[-1]:
                    want = want[:-1]
                got = dlg.result_lines()[row]["coords"]
                if [tuple(round(v, 6) for v in c) for c in got] != \
                        [tuple(round(v, 6) for v in c) for c in want]:
                    fails.append(f"the polygon editor row for zone {i} carries "
                                 f"{len(got)} vertices that are not that zone's")
            finally:
                dlg.deleteLater()
    finally:
        mw.close()
    return fails


def run():
    legs = [("interior click -> polygons editor", test_interior),
            ("features win over the interior", test_precedence),
            ("overlapping zones: smallest wins", test_overlap),
            ("profile files and empty space unchanged", test_profile_and_empty)]
    try:
        import PySide6                                   # noqa: F401
        legs += [("canvas pick -> edit_category", test_wire),
                 ("resolved row is the editor row", test_row_lands)]
    except Exception:
        print("  (PySide6 not installed — Studio legs skipped)")

    failures = []
    for name, fn in legs:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("polygon interior picking (Studio):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll polygon picking checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
