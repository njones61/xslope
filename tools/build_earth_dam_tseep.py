"""Build the transient (reservoir-drawdown) seepage SAMPLE workbooks for
``docs/seep/samples.md``.

Two worked transient samples share this one deterministic builder — each is the
transient variant of an existing *steady* sample, built from that committed base
file so the geometry never diverges. The recipe is identical for both: load the
base, retype the upstream head boundary as a submerged-only ``reservoir`` bound to
a falling ``pool`` series, attach per-zone storage properties and a v18 ``tseep``
sheet, declare units + a day time base, and round-trip through the packaged v18
template with ``save_slope_data_to_xlsx``.

  ``earth_dam``  -> docs/seep/files/xslope_earth_dam_tseep.xlsx
      Transient variant of the homogeneous cored earth dam of
      [Problem 3](../seep/samples.md) (``xslope_earth_dam1.xlsx``): the SAME
      cross-section, zones and downstream boundaries, driven by a falling
      upstream reservoir.

  ``johnson``    -> docs/seep/files/xslope_johnson_res_tseep.xlsx
      Transient variant of the zoned Johnson Reservoir dam of
      [Problem 5](../seep/samples.md#johnson-reservoir) (``xslope_johnson_res.xlsx``):
      shell over a low-permeability clay core carried down into the foundation.
      The zones are the story here — the core drains far slower than the shells,
      so it holds an elevated interior head long after the shells have emptied.

Physics choices (all stated in the docs entries):

  * Time base = DAY, so the conductivities must be expressed in units that share
    the schedule's time unit (the transient march balances storage against
    ``div(k grad h)``). The Johnson base file is already in **ft/day** — its steady
    discharge is the 1.958 ft^3/day/ft SEEP2D benchmark — so its k is reused
    verbatim. The homogeneous earth-dam base is in ft/yr, so representative
    granular-shell / clay-core conductivities are set in ft/day instead.
  * Storage (mat ``Ss``/``Sy``), from the transient.md storage tables, assigned by
    material type: granular shell (sand) ``Ss = 1e-4`` /ft, ``Sy = 0.22``; clay core
    ``Ss = 1e-3`` /ft, ``Sy = 0.03``; the Johnson silty-sand foundation
    ``Ss = 2e-4`` /ft, ``Sy = 0.15``. With the linear-front law and ``|h0| = 1`` ft the
    drainable-band storage is ~``Sy``, so the phreatic surface drains at the
    unconfined rate ``k/Sy``.
  * Schedule: the pool is held at full pool briefly, drawn down to the tailwater
    datum over ~45 days, then held while the field relaxes toward a new steady
    state (the boundary outflow decays to ~1% of its drawdown peak). ``stage_1``
    (full pool) and ``stage_2`` (end of drawdown, the maximum interior lag) are the
    critical rapid-drawdown pair.

Run:  PYTHONPATH=. python3 tools/build_earth_dam_tseep.py            # build both
      PYTHONPATH=. python3 tools/build_earth_dam_tseep.py johnson    # build one
"""

from __future__ import annotations

import os
import sys

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           default_template_path)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(REPO_ROOT, "docs", "seep", "files")

# --- representative storage by material type (1/ft, dimensionless) ----------- #
# From the transient.md storage tables; |h0| = 1 so the drainable band storage ~ Sy.
SAND_STORAGE = dict(Ss=1.0e-4, Sy=0.22)   # granular shell (fine sand)
CLAY_STORAGE = dict(Ss=1.0e-3, Sy=0.03)   # compacted clay core
SILT_STORAGE = dict(Ss=2.0e-4, Sy=0.15)   # silty-sand foundation


# --------------------------------------------------------------------------- #
#  Sample specifications
# --------------------------------------------------------------------------- #
# Each spec: base file, output file, per-material property updates (indexed to the
# base file's material order), the index of the upstream head boundary to retype as
# a reservoir, and the tseep control block (times/series/run controls/stage pair).
SAMPLES = {
    # Homogeneous cored earth dam (Problem 3).  Base k is in ft/yr, so k is reset to
    # representative ft/day values (shell ~ fine sand, core ~ silty clay).
    "earth_dam": dict(
        base="xslope_earth_dam1.xlsx",
        out="xslope_earth_dam_tseep.xlsx",
        mat_updates=[
            dict(k1=0.75, k2=0.25, **SAND_STORAGE),   # 0: shell
            dict(k1=0.012, k2=0.005, **CLAY_STORAGE),  # 1: core
        ],
        reservoir_idx=0,
        tseep=dict(
            times=[0.0, 2.0, 47.0],
            series={"pool": [18.0, 18.0, 2.0]},   # full pool -> tailwater over 45 d
            duration=360.0,
            save_interval=60.0,
            save_times=[2.0, 15.0, 30.0, 47.0, 80.0],
            stage_1=0.0,
            stage_2=47.0,
        ),
    ),
    # Zoned Johnson Reservoir dam (Problem 5).  Base k is already in ft/day (its
    # steady discharge IS the 1.958 SEEP2D benchmark), so k is reused verbatim and
    # only storage is added.  The reservoir is drawn from full pool (el 160) to the
    # tailwater datum (el 100) over 45 days; the low-k core paces the relaxation, so
    # the run is long (1000 d) before the outflow decays to ~1% of its peak.
    "johnson": dict(
        base="xslope_johnson_res.xlsx",
        out="xslope_johnson_res_tseep.xlsx",
        mat_updates=[
            dict(**SAND_STORAGE),   # 0: shell   (k = 1.0 ft/day, base)
            dict(**CLAY_STORAGE),   # 1: core    (k = 0.001 ft/day, base)
            dict(**SILT_STORAGE),   # 2: foundation (k = 0.1 ft/day, base)
        ],
        reservoir_idx=0,
        tseep=dict(
            times=[0.0, 5.0, 50.0],
            series={"pool": [160.0, 160.0, 100.0]},  # full pool -> tailwater over 45 d
            duration=1000.0,
            save_interval=200.0,
            save_times=[5.0, 35.0, 50.0, 80.0, 150.0, 300.0],
            stage_1=0.0,
            stage_2=50.0,
        ),
    ),
}


def build(name):
    spec = SAMPLES[name]
    sd = load_slope_data(os.path.join(FILES, spec["base"]))

    # per-zone storage (+ ft/day conductivities where the base file is not in /day)
    for i, upd in enumerate(spec["mat_updates"]):
        sd["materials"][i].update(upd)

    # declared units + a day time base (required for a transient run)
    sd["unit_system"] = "imperial"
    sd["time_unit"] = "day"

    # retype the upstream head boundary as a submerged-only reservoir driven by the
    # 'pool' series; drawn up the full submerged face, so a node the falling level
    # leaves above the waterline converts to an exit face.  The downstream tailwater
    # head and the exit face are left unchanged.
    series_name = next(iter(spec["tseep"]["series"]))
    head = sd["seepage_bc"]["specified_heads"][spec["reservoir_idx"]]
    head["kind"] = "reservoir"
    head["head"] = series_name

    sd["tseep"] = dict(spec["tseep"])

    out = os.path.join(FILES, spec["out"])
    save_slope_data_to_xlsx(sd, out, template=default_template_path())
    print(f"wrote {os.path.relpath(out, REPO_ROOT)}")
    return out


def main(argv):
    names = argv[1:] or list(SAMPLES)
    for name in names:
        if name not in SAMPLES:
            raise SystemExit(f"unknown sample {name!r}; choices: {list(SAMPLES)}")
        build(name)


if __name__ == "__main__":
    main(sys.argv)
