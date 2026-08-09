"""Regenerate the transient history figures for docs/seep/samples.md.

Both transient samples are drawn the same way and by the same function: the run
is marched once and :func:`xslope.plot_seep.plot_transient_history` draws it —
the reservoir level, the water table at an interior station and the top of the
seepage face above, the boundary inflow and outflow below, with the drawdown
shaded across both.

What this script adds to that figure is what belongs to the documentation rather
than to the plot: the station each sample is read at, the elevation window, and
the two sentence titles that say what the sample shows. The station is a choice —
the function's own default is the column that lags the reservoir most — so each
sample names the one its text is written about.

    PYTHONPATH=. python3 docs/seep/images/tseep_history.py            # both
    PYTHONPATH=. python3 docs/seep/images/tseep_history.py johnson    # one
"""
import contextlib
import io
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", ".."))

from xslope.fileio import load_slope_data                                # noqa: E402
from xslope.mesh import build_mesh_from_polygons, get_material_polygons  # noqa: E402
from xslope.plot_seep import plot_transient_history                      # noqa: E402
from xslope.seep import (build_seep_data, build_tseep_data,              # noqa: E402
                         run_transient_seepage)

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = os.path.join(HERE, "..", "files")

#: Per sample: the model, the mesh the sample is solved on, the interior station
#: the water table is followed at, and the elevation window the panel is drawn in.
#: ``sat_tol`` is how far below zero a node's pressure head may sit and still
#: count as saturated — a tenth of a foot, well inside one element on either mesh.
SAMPLES = {
    "earth_dam": dict(
        xlsx="xslope_earth_dam_tseep.xlsx",
        stem="earth_dam_tseep",
        target=110.0 / 64.0,
        station=30.0,
        station_halfwidth=2.5,
        sat_tol=0.1,
        elev_lim=(0, 20),
        titles=("Phreatic surface lags the reservoir; the exit point migrates "
                "down the face",
                "Inflow → 0 as the face drains; outflow spikes on released "
                "storage, then decays"),
    ),
    "johnson": dict(
        xlsx="xslope_johnson_res_tseep.xlsx",
        stem="johnson_res_tseep",
        target=10.0,
        station=280.0,
        station_halfwidth=15.0,
        sat_tol=0.1,
        elev_lim=(95, 170),
        titles=("Phreatic surface lags the reservoir; the exit point migrates "
                "down the face",
                "Inflow → 0 as the face drains; outflow spikes on released "
                "storage, then decays"),
    ),
}


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def solve(cfg):
    """``(seep_data, transient_solution)`` for one sample, marched once."""
    slope_data = load_slope_data(os.path.join(FILES, cfg["xlsx"]))
    mesh = _quiet(build_mesh_from_polygons, get_material_polygons(slope_data),
                  cfg["target"], "tri3")
    seep_data = _quiet(build_seep_data, mesh, slope_data)
    solution = _quiet(run_transient_seepage, seep_data,
                      build_tseep_data(slope_data), verbose=False)
    return seep_data, solution


def draw(cfg, seep_data, solution):
    fig = plt.figure(figsize=(7.4, 6.4))
    history = plot_transient_history(
        seep_data, solution, fig=fig, station=cfg["station"],
        station_halfwidth=cfg["station_halfwidth"], sat_tol=cfg["sat_tol"],
        elev_lim=cfg["elev_lim"], show_title=False)
    ax_elev, ax_flow = fig.axes
    for ax, title in zip((ax_elev, ax_flow), cfg["titles"]):
        ax.set_title(title, fontsize=10.5)

    # The interior station can settle BELOW the drawn-down pool: at quasi-
    # equilibrium the water table slopes from the low pool toward the downstream
    # exit, so a column inside the section drains under the upstream level.
    # Annotated where it happens, so the trace is not read as an error.
    phreatic, level = history["phreatic"], history["level"]
    if phreatic[-1] < level[-1] - history["sat_tol"]:
        ax_elev.annotate("interior surface drains\ntoward the downstream exit",
                         xy=(history["times"][-1], phreatic[-1]),
                         xytext=(0.60, 0.26), textcoords="axes fraction",
                         fontsize=8.0, color="#1f4e79", ha="left", va="center",
                         arrowprops=dict(arrowstyle="-|>", color="#1f4e79",
                                         lw=0.9))
    fig.tight_layout()
    path = os.path.join(HERE, f"{cfg['stem']}_history.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path, history


def main(argv):
    for name in (argv[1:] or list(SAMPLES)):
        cfg = SAMPLES[name]
        seep_data, solution = solve(cfg)
        path, history = draw(cfg, seep_data, solution)
        print(f"wrote {path}")
        print(f"  {len(seep_data['nodes'])} nodes, "
              f"{len(solution['frames'])} saved frames, "
              f"converged={solution['converged']}")
        for key in ("level", "phreatic", "exit_point", "inflow", "outflow"):
            values = history[key]
            print(f"  {key:11s} " + " ".join(f"{v:.3f}" for v in values))


if __name__ == "__main__":
    main(sys.argv)
