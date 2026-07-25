"""Build the transient (drawdown) earth-dam seepage SAMPLE workbook
``docs/seep/files/xslope_earth_dam_tseep.xlsx`` -> docs/seep/samples.md,
"Earth Dam — Reservoir Drawdown (Transient)".

This is the transient variant of the homogeneous cored earth dam of
[Problem 3](../seep/samples.md) (``xslope_earth_dam1.xlsx``): the SAME cross-section,
zones, and downstream boundaries, driven by a falling upstream reservoir. It is
built deterministically from that committed base file so the geometry never
diverges: load the base, retype the upstream head boundary as a submerged-only
``reservoir`` bound to a drawdown series, attach storage properties and a v18
``tseep`` sheet, declare units + a day time base, and round-trip through the
packaged v18 template with ``save_slope_data_to_xlsx``.

Physics choices (all stated in the docs entry):

  * Time base = DAY, so the conductivities are expressed in **ft/day** — the
    transient march balances storage against ``div(k grad h)`` and k must share the
    schedule's time unit (the same rule the SEEP/W corpus builders follow: k in
    m/s with time in seconds). The base file's ft/yr numbers are not reused
    verbatim; representative granular-shell / clay-core conductivities are set in
    ft/day (shell K ~ 0.75 ft/day ~ 3e-6 m/s, a fine sand; core K ~ 0.012 ft/day
    ~ 4e-8 m/s, a compacted silty clay), preserving the shell>>core contrast.
  * Storage (mat Ss/Sy), from the transient.md storage tables: shell (sand)
    Ss = 1e-4 /ft, Sy = 0.22; core (clay) Ss = 1e-3 /ft, Sy = 0.03. With the
    linear-front law and |h0| = 1 ft the drainable band storage is ~Sy, so the
    phreatic surface drains at the unconfined rate k/Sy.
  * Schedule: full pool (el 18) held briefly, drawn down to the tailwater level
    (el 2) over 45 days (t = 2 -> 47), then held. Duration 360 days lets the field
    reach quasi-equilibrium (the boundary outflow decays to ~1.5% of its drawdown
    peak; see the figure script). stage_1 = 0 (full-pool steady state) and
    stage_2 = 47 (end of drawdown, the maximum interior lag) are the critical
    rapid-drawdown pair.

Run:  PYTHONPATH=. python3 tools/build_earth_dam_tseep.py
"""

from __future__ import annotations

import os

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           default_template_path)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(REPO_ROOT, "docs", "seep", "files", "xslope_earth_dam1.xlsx")
OUT = os.path.join(REPO_ROOT, "docs", "seep", "files", "xslope_earth_dam_tseep.xlsx")

# --- representative properties (ft/day, 1/ft, dimensionless) ----------------- #
SHELL = dict(k1=0.75, k2=0.25, Ss=1.0e-4, Sy=0.22)   # granular shell (fine sand)
CORE = dict(k1=0.012, k2=0.005, Ss=1.0e-3, Sy=0.03)  # clay core (silty clay)

# --- reservoir drawdown series + run controls (time = day) ------------------- #
FULL_POOL = 18.0     # crest-region pool = the base file's upstream head
LOW_POOL = 2.0       # drawn down to the downstream tailwater datum
HOLD = 2.0           # brief full-pool hold before drawdown starts (day)
DRAWDOWN_END = 47.0  # HOLD + 45-day drawdown window
DURATION = 360.0     # long enough to reach quasi-equilibrium


def build():
    sd = load_slope_data(BASE)

    # storage + ft/day conductivities on the two zones (order: shell, core)
    sd["materials"][0].update(SHELL)
    sd["materials"][1].update(CORE)

    # declared units + a day time base (required for a transient run)
    sd["unit_system"] = "imperial"
    sd["time_unit"] = "day"

    # retype the upstream head boundary as a submerged-only reservoir driven by the
    # 'pool' series; drawn up the full submerged face (toe (0,0) -> el 18), so a
    # node the falling level leaves above the waterline converts to an exit face.
    heads = sd["seepage_bc"]["specified_heads"]
    heads[0]["kind"] = "reservoir"
    heads[0]["head"] = "pool"      # names the tseep series
    # heads[1] (downstream tailwater, el 2) and the exit face are left unchanged.

    # v18 tseep sheet: one shared time axis, the 'pool' reservoir series, run
    # controls, and the rapid-drawdown stage pair.
    sd["tseep"] = {
        "times": [0.0, HOLD, DRAWDOWN_END],
        "series": {"pool": [FULL_POOL, FULL_POOL, LOW_POOL]},
        "duration": DURATION,
        "save_interval": 60.0,
        "save_times": [2.0, 15.0, 30.0, 47.0, 80.0],
        "stage_1": 0.0,
        "stage_2": DRAWDOWN_END,
    }

    save_slope_data_to_xlsx(sd, OUT, template=default_template_path())
    print(f"wrote {os.path.relpath(OUT, REPO_ROOT)}")
    return OUT


if __name__ == "__main__":
    build()
