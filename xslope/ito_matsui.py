# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Ito & Matsui (1975) method for computing lateral force on passive piles in slopes.

The method models soil between adjacent piles as being in plastic equilibrium
and provides closed-form equations for the distributed lateral pressure on piles
as a function of depth. The total force per pile is obtained by integrating
the pressure from the ground surface down to the failure surface.

Reference:
    Ito, T. and Matsui, T. (1975). "Methods to Estimate Lateral Force Acting
    on Stabilizing Piles." Soils and Foundations, 15(4), 43-59.
"""

import numpy as np
from shapely.geometry import LineString


def ito_matsui_coefficients(D1, D, phi_deg):
    """
    Compute dimensionless Ito & Matsui coefficients A1 and A2.

    These coefficients define the distributed lateral force p(z) on a pile:
        p(z) = c * A1 + gamma * z * A2

    Parameters
    ----------
    D1 : float
        Clear spacing between piles (S - D).
    D : float
        Pile diameter (or width).
    phi_deg : float
        Soil friction angle in degrees.

    Returns
    -------
    A1 : float
        Cohesion coefficient.
    A2 : float
        Overburden coefficient.
    """
    # --- Phi = 0 baseline: Ito & Matsui (1975) Eq. (23) ---
    # The phi=0 case is derived independently from the c-phi case using a
    # different plastic equilibrium formulation. It serves as the physical
    # lower bound for all friction angles, since adding friction to a cohesive
    # soil can only increase arching (and thus pile force), never decrease it.
    S = D1 + D
    A1_phi0 = S * (3.0 * np.log(S / D1) + (D / D1) * np.tan(np.pi / 8) - 2.0) + 2.0 * D1
    A2_phi0 = D

    if phi_deg < 0.01:
        return A1_phi0, A2_phi0

    # --- General c-phi soil: Ito & Matsui (1975) Eqs. (12)-(14) ---
    phi = np.radians(phi_deg)
    Nphi = np.tan(np.pi / 4 + phi / 2) ** 2  # passive earth pressure coefficient
    sN = np.sqrt(Nphi)
    tp = np.tan(phi)

    # Geometric amplification factor (Eq. 5/11)
    R = (S / D1) ** (sN * tp + Nphi - 1)

    # Exponential plastic flow term (Eq. 10/11)
    E = np.exp(D / D1 * Nphi * tp * np.tan(np.pi / 8 + phi / 4))

    # A2: overburden coefficient, from Eq. (14) (c=0 case of Eq. 13)
    A2 = (S * R * E - D1) / Nphi

    # A1: cohesion coefficient, from Eq. (13)
    F = (2 * tp + 2 * sN + 1 / sN) / (sN * tp + Nphi - 1)
    A1 = (S * R * (E - 2 * sN * tp - 1) / (Nphi * tp)
          + S * F * (R - 1)
          + 2 * D1 / sN)

    # Apply phi=0 floor as safety net. The correct Eq. (13) converges
    # smoothly to Eq. (23) as phi->0, so this should rarely activate.
    A1 = max(A1, A1_phi0)
    A2 = max(A2, A2_phi0)

    return A1, A2


def intersect_pile_with_materials(pile_x, y_ground, y_failure, polygons, materials):
    """
    Find soil layers along the pile between the ground surface and failure surface.

    Uses the unified material-zone polygon representation (see
    plans/plan_polygons.md): each polygon's vertical extent at the pile location
    gives a layer band, the same way slice weights are computed.

    Parameters
    ----------
    pile_x : float
        X-coordinate of the pile (assumed vertical for layer intersection).
    y_ground : float
        Y-coordinate of the ground surface at the pile location.
    y_failure : float
        Y-coordinate of the failure surface at the pile location.
    polygons : list
        slope_data['polygons'] — list of dicts with 'polygon' (shapely Polygon)
        and 'mat_id' (0-based) keys.
    materials : list
        List of material property dicts with 'c', 'phi', 'gamma' keys.

    Returns
    -------
    segments : list of dict
        Each dict contains:
        - 'z_top': depth from ground surface to top of segment
        - 'z_bot': depth from ground surface to bottom of segment
        - 'c': cohesion
        - 'phi': friction angle (degrees)
        - 'gamma': unit weight
    """
    # Build the per-polygon top/bottom boundary tables (shared with slice.py).
    from .slice import _build_polygon_edges
    poly_edges = _build_polygon_edges(polygons)

    # For each polygon, its vertical extent at pile_x defines the layer band; clip
    # that to [failure surface, ground surface] (polygons are ordered top-to-bottom).
    segments = []
    for idx, pe in enumerate(poly_edges):
        if not (pe['xmin'] - 1e-9 <= pile_x <= pe['xmax'] + 1e-9):
            continue
        poly_top = np.interp(pile_x, pe['txs'], pe['tys'])
        poly_bot = np.interp(pile_x, pe['bxs'], pe['bys'])
        if np.isnan(poly_top) or np.isnan(poly_bot):
            continue

        overlap_top = min(y_ground, poly_top)
        overlap_bot = max(y_failure, poly_bot)
        h = overlap_top - overlap_bot
        if h <= 0:
            continue

        # Convert y-coordinates to depths from ground surface (z increases downward)
        z_top = y_ground - overlap_top
        z_bot = y_ground - overlap_bot

        mat_id = pe['mat_id']
        if mat_id is not None and 0 <= mat_id < len(materials):
            mat_index = mat_id
        else:
            mat_index = idx
        mat = materials[mat_index]

        segments.append({
            'z_top': z_top,
            'z_bot': z_bot,
            'c': mat.get('c', 0.0),
            'phi': mat.get('phi', 0.0),
            'gamma': mat.get('gamma', 0.0),
        })

    return segments


def compute_ito_matsui_force(D_pile, S, segments):
    """
    Compute total pile force per unit width using Ito & Matsui theory.

    Integrates p(z) = c*A1 + gamma*z*A2 piecewise over each soil layer
    from ground surface to failure surface, then divides by spacing S.

    Parameters
    ----------
    D_pile : float
        Pile diameter.
    S : float
        Center-to-center spacing.
    segments : list of dict
        Soil layer segments from intersect_pile_with_materials().

    Returns
    -------
    H : float
        Force per unit width of slope (F_pile / S).
    F_pile : float
        Total force per pile.
    """
    D1 = S - D_pile  # clear spacing between piles

    if D1 <= 0:
        return 0.0, 0.0

    F_pile = 0.0
    for seg in segments:
        z_top = seg['z_top']
        z_bot = seg['z_bot']
        c = seg['c']
        phi = seg['phi']
        gamma = seg['gamma']

        A1, A2 = ito_matsui_coefficients(D1, D_pile, phi)

        # Piecewise integration: F_j = c*A1*(z_bot - z_top) + gamma*A2*(z_bot^2 - z_top^2)/2
        dz = z_bot - z_top
        F_j = c * A1 * dz + gamma * A2 * (z_bot ** 2 - z_top ** 2) / 2.0
        F_pile += F_j

    H = F_pile / S if S > 0 else 0.0
    return H, F_pile


def compute_ito_matsui_force_and_moment_arm(D_pile, S, segments):
    """
    Compute pile force and moment arm from Ito & Matsui pressure distribution.

    The moment arm L_m is the distance from the failure surface (bottom of
    the pressure diagram) to the centroid of the pressure resultant. Used
    for the M_cap structural capacity check: F_pile_max = M_cap / L_m.

    Parameters
    ----------
    D_pile : float
        Pile diameter.
    S : float
        Center-to-center spacing.
    segments : list of dict
        Soil layer segments from intersect_pile_with_materials().

    Returns
    -------
    H : float
        Force per unit width of slope (F_pile / S).
    F_pile : float
        Total force per pile.
    L_m : float
        Moment arm from failure surface to pressure centroid.
        Returns 0.0 if F_pile is zero.
    """
    D1 = S - D_pile

    if D1 <= 0 or not segments:
        return 0.0, 0.0, 0.0

    F_pile = 0.0
    moment_about_base = 0.0

    # Total depth = max z_bot across all segments (failure surface depth)
    z_max = max(seg['z_bot'] for seg in segments)

    for seg in segments:
        z_top = seg['z_top']
        z_bot = seg['z_bot']
        c = seg['c']
        phi = seg['phi']
        gamma = seg['gamma']

        A1, A2 = ito_matsui_coefficients(D1, D_pile, phi)

        # Force for this segment
        dz = z_bot - z_top
        F_j = c * A1 * dz + gamma * A2 * (z_bot ** 2 - z_top ** 2) / 2.0
        F_pile += F_j

        # Moment of p(z) about the failure surface (z_max).
        # Lever arm at depth z is (z_max - z).
        # Integral of p(z)*(z_max - z) dz from z_top to z_bot:
        #   constant term: c*A1 * [z_max*dz - (z_bot^2 - z_top^2)/2]
        #   linear term:   gamma*A2 * [z_max*(z_bot^2 - z_top^2)/2 - (z_bot^3 - z_top^3)/3]
        M_j = (c * A1 * (z_max * dz - (z_bot ** 2 - z_top ** 2) / 2.0)
               + gamma * A2 * (z_max * (z_bot ** 2 - z_top ** 2) / 2.0
                               - (z_bot ** 3 - z_top ** 3) / 3.0))
        moment_about_base += M_j

    H = F_pile / S if S > 0 else 0.0
    L_m = moment_about_base / F_pile if F_pile > 0 else 0.0

    return H, F_pile, L_m
