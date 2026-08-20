"""Restore the steady "Earth Dam with Core" sample workbooks to SI.

The model is the owner's SEEP2D tutorial dam, and it has always been metric: the
source sketch (``docs/seep/images/earth_dam1.png``) dimensions it 22 m high and
110 m long, with an 18 m pool over a 2 m tailwater, and labels the shell
``kx = 46 / ky = 18 m/yr`` over a ``4.5 / 1.8 m/yr`` core. The committed
workbooks declared Imperial units against those same numbers and carried a
shell ``k1 = 56`` that matches no label on the sketch. This script restores the
declaration and the conductivity to the model.

What it changes, and nothing else:

  * ``main`` Units selector  -> SI, Time selector -> ``yr``. The geometry numbers
    are already metres, so no coordinate is touched. The time base is declared
    because ``xslope.units.labels()`` builds the ``k`` and flow-rate unit strings
    from it -- without it a metric model plots a flow rate with no unit at all.
  * ``main`` unit weight of water -> 9.81 kN/m3. The cell held 62.4 pcf, which an
    SI declaration rejects as a 536% divergence from canonical, and which would
    have reported the pore pressures in the wrong unit by a factor of 6.4. The SI
    siblings (``xslope_clay_blanket.xlsx``) carry 9.81 in the same cell.
  * shell ``k1`` 56 -> 46 (the sketch's label), ``k2`` 18 unchanged;
    core ``k1`` 4.5 / ``k2`` 1.8 unchanged. All in m/yr.
  * ``h0`` -1 -> -0.3 for both materials: the linear-front reference suction is a
    length, and -1 ft is -0.3 m. ``kr0`` stays 0.01.
  * van Genuchten variant only: the ``a`` column is Carsel & Parrish alpha
    converted to the model's length unit, so the 1/ft values become 1/m --
    sandy loam 0.075/cm -> 7.5 /m (shell), loam 0.036/cm -> 3.6 /m (core).
    ``n`` is dimensionless and unchanged.

Every write goes through ``xslope.fileio``: load, mutate the dict, save back
through the packaged template. Run once; re-running is idempotent.

    PYTHONPATH=. python3 tools/convert_earth_dam_to_si.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from xslope.fileio import (load_slope_data, save_slope_data_to_xlsx,   # noqa: E402
                           default_template_path)
from xslope.units import GAMMA_W                                       # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(REPO_ROOT, "docs", "seep", "files")

#: Carsel & Parrish alpha in 1/cm x 100 -> 1/m (see docs/seep/overview.md):
#: sandy loam 0.075 -> 7.5, loam 0.036 -> 3.6. Written as the exact decimals the
#: table implies rather than a float product, so the cell reads 3.6 and not
#: 3.5999999999999996.
VG_ALPHA_PER_M = {"sandy loam": 7.5, "loam": 3.6}

#: (filename, per-material updates indexed to the file's material order).
#: Shell = granular fill, sandy loam unsaturated character; core = compacted clay,
#: loam unsaturated character. k in m/yr.
TARGETS = [
    ("xslope_earth_dam1.xlsx", [
        dict(k1=46.0, k2=18.0, kr0=0.01, h0=-0.3),               # 0: shell
        dict(k1=4.5, k2=1.8, kr0=0.01, h0=-0.3),                 # 1: core
    ]),
    ("xslope_earth_dam1_vg.xlsx", [
        dict(k1=46.0, k2=18.0, kr0=0.01, h0=-0.3,
             vg_a=VG_ALPHA_PER_M["sandy loam"], vg_n=1.89),      # 0: shell
        dict(k1=4.5, k2=1.8, kr0=0.01, h0=-0.3,
             vg_a=VG_ALPHA_PER_M["loam"], vg_n=1.56),            # 1: core
    ]),
]


def convert(filename, mat_updates):
    path = os.path.join(FILES, filename)
    sd = load_slope_data(path)
    for i, upd in enumerate(mat_updates):
        sd["materials"][i].update(upd)
    sd["unit_system"] = "si"
    sd["time_unit"] = "yr"
    sd["gamma_water"] = GAMMA_W["si"]
    save_slope_data_to_xlsx(sd, path, template=default_template_path())
    print(f"wrote {os.path.relpath(path, REPO_ROOT)}  "
          f"[{sd['unit_system']} / {sd['time_unit']}, "
          f"gamma_w={sd['gamma_water']}]")
    for m in sd["materials"]:
        print(f"    {m['name']:<6} k1={m['k1']} k2={m['k2']} "
              f"kr0={m['kr0']} h0={m['h0']} "
              f"a={m.get('vg_a')} n={m.get('vg_n')} unsat={m.get('unsat')}")


def main():
    for filename, mat_updates in TARGETS:
        convert(filename, mat_updates)


if __name__ == "__main__":
    main()
