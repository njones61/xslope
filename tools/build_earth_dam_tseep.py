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
      Transient variant of the cored earth dam of
      [Problem 3](../seep/samples.md) (``xslope_earth_dam1.xlsx``): the SAME
      cross-section, zones and downstream boundaries, driven by a falling
      upstream reservoir. Metric, like its base.

  ``johnson``    -> docs/seep/files/xslope_johnson_res_tseep.xlsx
      Transient variant of the zoned Johnson Reservoir dam of
      [Problem 5](../seep/samples.md#johnson-reservoir) (``xslope_johnson_res.xlsx``):
      shell over a low-permeability clay core carried down into the foundation.
      The zones are the story here — the core drains far slower than the shells,
      so it holds an elevated interior head long after the shells have emptied.

Physics choices (all stated in the docs entries):

  * Each sample keeps its base file's UNIT SYSTEM — the earth dam is metric, the
    Johnson dam is imperial — so a conversion of one never reaches the other.
  * Time base = DAY, so the conductivities must be expressed in units that share
    the schedule's time unit (the transient march balances storage against
    ``div(k grad h)``). The Johnson base file is already in **ft/day** — its steady
    discharge is the 1.958 ft^3/day/ft SEEP2D benchmark — so its k is reused
    verbatim. The earth-dam base is in m/yr, so representative granular-shell /
    clay-core conductivities are set in **m/day** instead. They are deliberately
    not the steady sketch's 46/18 m/yr: at that rate the shell drains at
    ``k2/Sy`` = 0.22 m/yr and could not follow a 45-day drawdown at all, so the
    model would have no transient story to tell.
  * Storage (mat ``Ss``/``Sy``), from the transient.md storage tables, assigned by
    material type in the model's OWN length unit — ``Ss`` is a 1/length, so the
    metric values are the imperial ones times ~3.28: granular shell (sand)
    ``Ss = 3e-4`` /m or ``1e-4`` /ft, ``Sy = 0.22``; clay core ``Ss = 3e-3`` /m or
    ``1e-3`` /ft, ``Sy = 0.03``; the Johnson silty-sand foundation
    ``Ss = 2e-4`` /ft, ``Sy = 0.15``. With the linear-front law and a small
    ``|h0|`` the drainable-band storage is ~``Sy``, so the phreatic surface drains
    at the unconfined rate ``k/Sy``.
  * Schedule: the pool is held at full pool briefly, drawn down to the tailwater
    datum over ~45 days, then held while the field relaxes toward a new steady
    state (the boundary outflow decays to ~1% of its drawdown peak). ``stage_1``
    (full pool) and ``stage_2`` (end of drawdown, the maximum interior lag) are the
    critical rapid-drawdown pair.

The same recipe also builds Tutorial SEEP-3's file PAIR into docs/tutorials/files/,
from the same base file and the same storage/conductivity choices, so the tutorial
and the sample can never describe different soils:

  ``seep03_start``  -> docs/tutorials/files/xslope_earth_dam_drawdown_start.xlsx
      The STARTER. Geometry, zones, units and materials (storage included) of the
      transient earth dam, with the ``seep bc`` sheet and the ``tseep`` sheet left
      blank: the reader builds every boundary condition and the whole schedule by
      hand, and the file loads STEADY until they do.

  ``seep03``        -> docs/tutorials/files/xslope_earth_dam_drawdown.xlsx
      The COMPLETED file — the skip-ahead download. (A basename distinct from the
      sample's on purpose: the tutorial file differs from the sample by design,
      and a shared name invites someone to "fix" the difference.) Identical to the
      starter plus
      the boundary set (reservoir bound to the ``pool`` series, tailwater head, exit
      face) and the transient schedule. It carries NO ``stage_1``/``stage_2`` and no
      ``stability_time``: SEEP-3 is a seepage page, and the rapid-drawdown staging
      those flags exist for belongs to the stability tutorial that follows it. Every
      other physics choice is the sample's.

Run:  PYTHONPATH=. python3 tools/build_earth_dam_tseep.py            # build all
      PYTHONPATH=. python3 tools/build_earth_dam_tseep.py johnson    # build one
      PYTHONPATH=. python3 tools/build_earth_dam_tseep.py seep03     # the pair
"""

from __future__ import annotations

import os
import sys

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           default_template_path)
from xslope.units import GAMMA_W

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(REPO_ROOT, "docs", "seep", "files")
TUTORIAL_FILES = os.path.join(REPO_ROOT, "docs", "tutorials", "files")

# --- representative storage by material type --------------------------------- #
# From the transient.md storage tables. Ss is a 1/length, so it is quoted in the
# model's own length unit and the metric values are the imperial ones times ~3.28;
# Sy is a drainable water content and is dimensionless, so it is shared. With the
# linear-front law the drainable band storage is ~Sy.
SAND_STORAGE_SI = dict(Ss=3.0e-4, Sy=0.22)   # granular shell (fine sand), 1/m
CLAY_STORAGE_SI = dict(Ss=3.0e-3, Sy=0.03)   # compacted clay core, 1/m
SAND_STORAGE = dict(Ss=1.0e-4, Sy=0.22)      # granular shell (fine sand), 1/ft
CLAY_STORAGE = dict(Ss=1.0e-3, Sy=0.03)      # compacted clay core, 1/ft
SILT_STORAGE = dict(Ss=2.0e-4, Sy=0.15)      # silty-sand foundation, 1/ft


# --------------------------------------------------------------------------- #
#  Sample specifications
# --------------------------------------------------------------------------- #
# Each spec: base file, output file, per-material property updates (indexed to the
# base file's material order), the index of the upstream head boundary to retype as
# a reservoir, and the tseep control block (times/series/run controls/stage pair).
SAMPLES = {
    # Cored earth dam (Problem 3).  Metric; base k is in m/yr, so k is reset to
    # representative m/day values (shell ~ fine sand, core ~ compacted clay).
    "earth_dam": dict(
        base="xslope_earth_dam1.xlsx",
        out="xslope_earth_dam_tseep.xlsx",
        unit_system="si",
        mat_updates=[
            dict(k1=1.5, k2=0.5, **SAND_STORAGE_SI),     # 0: shell
            dict(k1=0.012, k2=0.005, **CLAY_STORAGE_SI),  # 1: core
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
        unit_system="imperial",
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


def _soil(spec):
    """The base model with this sample's soil made transient: per-zone storage,
    per-day conductivities where the base file is not in /day, the base file's own
    declared unit system and a day time base (required for a transient run).

    The unit system is named per sample rather than assumed, so the metric earth dam
    and the imperial Johnson dam can share this one call without either inheriting
    the other's declaration.

    Everything downstream of this — the sample, the tutorial starter and the
    tutorial's completed file — is the same soil, because it is the same call.
    """
    sd = load_slope_data(os.path.join(FILES, spec["base"]))
    for i, upd in enumerate(spec["mat_updates"]):
        sd["materials"][i].update(upd)
    sd["unit_system"] = spec["unit_system"]
    sd["gamma_water"] = GAMMA_W[spec["unit_system"]]
    sd["time_unit"] = "day"
    return sd


def _bind_reservoir(sd, spec):
    """Retype the upstream head boundary as a submerged-only reservoir driven by the
    series the schedule defines.

    Drawn up the full submerged face, so a node the falling level leaves above the
    waterline converts to an exit face. The downstream tailwater head and the exit
    face are left unchanged.
    """
    series_name = next(iter(spec["tseep"]["series"]))
    head = sd["seepage_bc"]["specified_heads"][spec["reservoir_idx"]]
    head["kind"] = "reservoir"
    head["head"] = series_name
    return sd


def _write(sd, directory, filename):
    out = os.path.join(directory, filename)
    save_slope_data_to_xlsx(sd, out, template=default_template_path())
    print(f"wrote {os.path.relpath(out, REPO_ROOT)}")
    return out


def build(name):
    spec = SAMPLES[name]
    sd = _bind_reservoir(_soil(spec), spec)
    sd["tseep"] = dict(spec["tseep"])
    return _write(sd, FILES, spec["out"])


# --------------------------------------------------------------------------- #
#  Tutorial SEEP-3's file pair (docs/tutorials/files/)
# --------------------------------------------------------------------------- #
#: Both tutorial files are the ``earth_dam`` sample's soil. The pair differs from
#: the sample in exactly two ways, both deliberate: the files carry no companion
#: sidecars (the reader meshes and solves), and the completed file carries no
#: rapid-drawdown staging — SEEP-3 stops at the pore-pressure field, and the
#: ``stage_1`` / ``stage_2`` / ``stability_time`` flags belong to the stability
#: tutorial that consumes it.
SEEP03 = "earth_dam"
SEEP03_START_OUT = "xslope_earth_dam_drawdown_start.xlsx"
SEEP03_OUT = "xslope_earth_dam_drawdown.xlsx"
#: The tseep controls the completed file keeps, in the schedule's own order. Naming
#: them here rather than deleting the stage keys afterwards means a stage flag added
#: to the SAMPLE never leaks into the tutorial file by omission.
SEEP03_TSEEP_KEYS = ("times", "series", "duration", "save_interval", "save_times")


def build_seep03_start():
    """The starter: the transient soil, no boundary conditions, no schedule.

    Both sheets are left blank rather than partly filled. A ``seep bc`` sheet with
    the tailwater already in it would make the page's boundary section a correction
    exercise, and a ``tseep`` sheet with the controls already in it would load the
    file transient before the reader has entered a series — so the file loads
    STEADY, with nothing to solve, which is the state the page starts from.
    """
    sd = _soil(SAMPLES[SEEP03])
    sd["seepage_bc"] = {"specified_heads": [], "specified_fluxes": [],
                        "exit_face": []}
    sd["seepage_bc2"] = {"specified_heads": [], "specified_fluxes": [],
                         "exit_face": []}
    sd["has_seepage_bc2"] = False
    sd["tseep"] = None
    return _write(sd, TUTORIAL_FILES, SEEP03_START_OUT)


def build_seep03():
    """The completed file: the sample's boundary set and schedule, no staging."""
    spec = SAMPLES[SEEP03]
    sd = _bind_reservoir(_soil(spec), spec)
    sd["tseep"] = {k: spec["tseep"][k] for k in SEEP03_TSEEP_KEYS}
    return _write(sd, TUTORIAL_FILES, SEEP03_OUT)


TUTORIALS = {
    "seep03_start": build_seep03_start,
    "seep03_done": build_seep03,
}


def main(argv):
    names = argv[1:] or (list(SAMPLES) + list(TUTORIALS))
    if names == ["seep03"]:
        names = list(TUTORIALS)
    for name in names:
        if name in SAMPLES:
            build(name)
        elif name in TUTORIALS:
            TUTORIALS[name]()
        else:
            raise SystemExit(f"unknown target {name!r}; choices: "
                             f"{list(SAMPLES) + list(TUTORIALS)}")


if __name__ == "__main__":
    main(sys.argv)
