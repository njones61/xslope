"""Build the rapid-drawdown LEM worked-example workbook for
``docs/lem/rapid.md`` (the "Rapid Drawdown from a Transient Solution" section).

This is the *stability* companion to the Johnson Reservoir transient **seepage**
sample (``docs/seep/files/xslope_johnson_res_tseep.xlsx``, samples §9). It is built
deterministically from that committed seepage file so the cross-section, zones,
storage, reservoir schedule and ``stage_1``/``stage_2`` times never diverge from the
seepage sample — the two files are the same dam seen from the seepage side and the
stability side. The seepage sample is left untouched (kept clean); the rapid-drawdown
LEM inputs live here in a ``_rapid`` variant.

Two additions turn the seepage sample into a rapid-drawdown model:

  * **Undrained core.** The compacted-clay core (material 2) is given the Duncan
    d/psi pair of its ``Kc = 1`` (isotropically consolidated, undrained) envelope,
    ``d = 250`` psf, ``psi = 14`` deg — below the drained ``c' = 400`` psf,
    ``phi' = 18`` deg. The free-draining sand shell and the more permeable silty-sand
    foundation are analysed drained (no d/psi), the textbook shell-versus-core
    dichotomy the rapid.md page sets out: the low-permeability core is the zone that
    cannot drain over a 45-day drawdown, so only it is carried undrained into Stage 2.
  * **The critical upstream circle.** A single circle on the upstream (reservoir) face
    — the slope the drawdown destabilises — located as the critical surface by a
    rapid-drawdown circular search over the transient-staged fields. It toes near the
    upstream base (x ~ 192), daylights just past the crest (x ~ 415) and cuts twelve
    slices of core, so the undrained zone governs Stage 2.

Everything else is inherited from the seepage sample and needs no change:

  * The **Stage 1 water load** is the full-pool distributed load already on the
    seepage file (the reservoir standing against the upstream face, 3744 psf at the
    toe tapering to 0 at the el-160 pool line). The **Stage 2 water load is empty**:
    the schedule draws the pool all the way down to the el-100 tailwater datum, so at
    ``stage_2`` no water remains on the upstream slope — the removal of that load is
    the destabilising essence of rapid drawdown.
  * ``u = seep`` on every material: the two stage pore-pressure fields are supplied at
    run time by ``stage_transient_for_drawdown`` from the transient frames at
    ``stage_1`` (t = 0, full pool) and ``stage_2`` (t = 50, end of drawdown). No
    ``_seep.csv`` / ``_seep2.csv`` companions are shipped — the transient staging path
    builds both fields in memory.
  * Declared imperial units and a ``day`` time base, and the tseep sheet
    (schedule + storage + ``stage_1 = 0`` / ``stage_2 = 50``), all carried over verbatim.

Run:  PYTHONPATH=. python3 tools/build_johnson_rapid.py
"""

from __future__ import annotations

import os

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,
                           default_template_path)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEEP_FILES = os.path.join(REPO_ROOT, "docs", "seep", "files")
LEM_FILES = os.path.join(REPO_ROOT, "docs", "lem", "files")

BASE = os.path.join(SEEP_FILES, "xslope_johnson_res_tseep.xlsx")
OUT = os.path.join(LEM_FILES, "xslope_johnson_res_rapid.xlsx")

# Undrained core Kc=1 (CU) envelope, below the drained c'=400/phi'=18.
CORE_D, CORE_PSI = 250.0, 14.0

# Critical upstream circle (located by a rapid-drawdown circular search over the
# transient-staged fields). Yo - Depth = R keeps the tangent elevation consistent.
CIRCLE = dict(Xo=275.0, Yo=235.0, R=160.0, Depth=75.0)


def build():
    sd = load_slope_data(BASE)

    # undrained compacted-clay core (material index 1); shell + foundation stay drained
    sd["materials"][1]["d"] = CORE_D
    sd["materials"][1]["psi"] = CORE_PSI

    # the critical upstream (reservoir-face) circle. max_depth is left at the base
    # file's value (the full foundation depth) — the fixed circle carries its own
    # tangent Depth, and overriding max_depth would clip the domain through the core.
    sd["circular"] = True
    sd["circles"] = [dict(CIRCLE)]

    out = OUT
    save_slope_data_to_xlsx(sd, out, template=default_template_path())
    print(f"wrote {os.path.relpath(out, REPO_ROOT)}")
    return out


if __name__ == "__main__":
    build()
