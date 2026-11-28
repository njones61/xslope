"""The seepage plot's contour lines can be turned off.

plot_seep_solution(contour_lines=False) draws no contour lines (the artist with gid
HEAD_CONTOURS) and no "... contour" legend entry, while flow lines and the phreatic
surface stay; the default still draws both. Studio's seepage Display panel carries a
"Contour lines" box, ticked by default, that feeds the option (steady and transient).

Run from the repo root:  PYTHONPATH=. python test/seep_contour_lines_check.py
"""

import os
import sys

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")

MODEL = os.path.join(_REPO, "docs", "seep", "files", "xslope_earth_dam2.xlsx")


def _solve():
    import contextlib
    import io
    from xslope.fileio import load_slope_data
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_size_regions)
    from xslope.seep import build_seep_data, run_seepage_analysis
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(MODEL)
        mesh = sd.get("mesh")
        if mesh is None:
            xs = [x for x, _ in sd["ground_surface"].coords]
            mesh = build_mesh_from_polygons(get_material_polygons(sd),
                                            (max(xs) - min(xs)) / 40, "tri6",
                                            size_regions=extract_size_regions(sd))
        seep = build_seep_data(mesh, sd)
        sol = run_seepage_analysis(seep, tol=1e-4)
    return seep, sol


def _drawn(fig):
    gids = set()
    labels = []
    for ax in fig.axes:
        for a in ax.get_children():
            if a.get_gid():
                gids.add(a.get_gid())
        leg = ax.get_legend()
        if leg is not None:
            labels += [t.get_text() for t in leg.get_texts()]
    for leg in fig.legends:
        labels += [t.get_text() for t in leg.get_texts()]
    return gids, labels


def main():
    import contextlib
    import io
    import matplotlib.pyplot as plt
    from xslope.plot_seep import plot_seep_solution
    failures = []
    seep, sol = _solve()

    for on in (True, False):
        fig = plt.figure(figsize=(10, 6))
        with contextlib.redirect_stdout(io.StringIO()):
            plot_seep_solution(seep, sol, fig=fig, fill_contours=False,
                               contour_lines=on)
        gids, labels = _drawn(fig)
        has_lines = "HEAD_CONTOURS" in gids
        has_key = any(lab.endswith("contour") for lab in labels)
        print(f"  contour_lines={on}: lines drawn {has_lines}, legend entry {has_key}, "
              f"flow line key {'Flow line' in labels}")
        if has_lines != on or has_key != on:
            failures.append(f"contour_lines={on}: lines {has_lines}, key {has_key}")
        if "Flow line" not in labels:
            failures.append(f"contour_lines={on}: flow lines lost")
        plt.close(fig)

    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)   # noqa: F841
    from studio.display_panels import SeepDisplayPanel
    for transient in (False, True):
        panel = SeepDisplayPanel([{"name": "a"}], transient=transient)
        box = getattr(panel, "contour_lines", None)
        if box is None:
            failures.append(f"transient={transient}: no Contour lines box")
            continue
        default = panel.options().get("contour_lines")
        box.setChecked(False)
        off = panel.options().get("contour_lines")
        print(f"  panel transient={transient}: default {default}, unticked {off}")
        if default is not True or off is not False:
            failures.append(f"transient={transient}: options {default} -> {off}")

    if failures:
        print("FAIL:\n  " + "\n  ".join(failures))
        sys.exit(1)
    print("Seepage contour-line checks passed.")


if __name__ == "__main__":
    main()
