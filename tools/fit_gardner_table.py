"""Fit Gardner (power form) parameters to the Carsel & Parrish van Genuchten curves.

Produces the fitted-parameter table in docs/seep/overview.md's Gardner section.
Gardner has no published texture table of its own, so each USDA texture's
(a, n) pair here is a least-squares fit — in log10 kr over 0.01 to 100 ft
(0.003 to 30.5 m) of suction, with both curves clipped at kr = 1e-4 — to that
texture's van Genuchten curve from the Carsel & Parrish table on the same page.
The same procedure calibrates the Gardner runs in Tutorial SEEP-2.

``FIT_CLIP`` is a property of the fit, not of the solver: it bounds the band the
two laws are matched over, so the dry tail — where both curves are decades below
anything that carries water and the power form cannot follow the Mualem
integral — does not dominate a least-squares fit in log kr. It is independent of
``seep.kr_min``, which is a numerical guard on the conduction matrix, and the
shipped pairs below are the ones fitted over this band.

`n` is unit-invariant; `a` carries units of (1/length)^n and converts exactly
by a_ft = a_m * 0.3048**n. The fit is done in meters and converted.

    python tools/fit_gardner_table.py    # prints the markdown rows
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scipy.optimize import curve_fit

from xslope.seep import kr_gardner_vec, kr_vg_vec

#: Carsel & Parrish (1988) van Genuchten alpha (1/cm) and n, as tabulated in
#: docs/seep/overview.md.
CARSEL_PARRISH = [
    ("Sand", 0.145, 2.68), ("Loamy sand", 0.124, 2.28),
    ("Sandy loam", 0.075, 1.89), ("Loam", 0.036, 1.56),
    ("Silt", 0.016, 1.37), ("Silt loam", 0.020, 1.41),
    ("Sandy clay loam", 0.059, 1.48), ("Clay loam", 0.019, 1.31),
    ("Silty clay loam", 0.010, 1.23), ("Sandy clay", 0.027, 1.23),
    ("Silty clay", 0.005, 1.09), ("Clay", 0.008, 1.09),
]
FIT_CLIP = 1e-4
#: Same physical suction band as Tutorial SEEP-2's fits: 0.01–100 ft.
SUCTION_M = np.logspace(np.log10(0.01 * 0.3048), np.log10(100 * 0.3048), 200)


def fit_texture(alpha_cm, n_vg):
    psi = -SUCTION_M
    y = np.log10(kr_vg_vec(psi, alpha_cm * 100.0, n_vg, FIT_CLIP))
    logs = np.log10(SUCTION_M)

    def model(ls, a, n):
        return np.log10(kr_gardner_vec(-(10 ** ls), a, n, FIT_CLIP))

    (a_m, n), _ = curve_fit(model, logs, y, p0=((alpha_cm * 100.0) ** n_vg, n_vg),
                            maxfev=20000)
    rms = float(np.sqrt(np.mean((model(logs, a_m, n) - y) ** 2)))
    return a_m, n, a_m * 0.3048 ** n, rms


def main():
    print("| Soil texture | `a` (ψ in m) | `a` (ψ in ft) | `n` | RMS log₁₀k_r |")
    print("| --- | --- | --- | --- | --- |")
    for name, alpha_cm, n_vg in CARSEL_PARRISH:
        a_m, n, a_ft, rms = fit_texture(alpha_cm, n_vg)
        print(f"| {name} | {a_m:.3g} | {a_ft:.3g} | {n:.2f} | {rms:.2f} |")


if __name__ == "__main__":
    main()
