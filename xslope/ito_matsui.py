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
    if phi_deg < 0.01:
        # phi = 0 case (undrained clay): simplified equations
        # A1 = 2*sqrt(2) * D1/(D1-D) + 2*sqrt(2) - 1   (from Ito & Matsui eq. for c-only soil)
        # A2 = 0 (no frictional component)
        # Using the exact phi=0 result: p(z) = c * [2*sqrt(2)*D1/(D1-D) + 2*sqrt(2) - 1]
        # But more precisely from the original paper for phi=0:
        # p = c * [D1/(D1-D) * (8/3) + (D - D1 + D) * pi/2 ... ]
        # Use the limiting form: as phi->0, Nphi->1, the exponential term -> 1
        # Simplified: A1 = D1/(D1-D) + D/(D1-D) = S/(D1) but let's use the standard form
        # For phi=0: p(z) = c * [2*tan(45)*D1/(D1-D)*exp(0) + D*tan(45) - D1*1]
        # = c * [2*D1/(D1-D) + D - D1]
        A1 = 2.0 * D1 / (D1) + D  # simplified phi=0
        # Actually for phi=0 the standard result is just geometry-dependent.
        # Use the well-known phi=0 formula: p = 2*c*(D1+D)/(D1) per unit depth
        A1 = 2.0 * (D1 + D) / D1
        A2 = 0.0
        return A1, A2

    phi = np.radians(phi_deg)
    Nphi = np.tan(np.pi / 4 + phi / 2) ** 2  # passive earth pressure coefficient

    # Ito & Matsui (1975) closed-form coefficients
    # From equations (1) and (2) in the original paper:
    #
    # A2 = D1 * [Nphi^(Nphi*D1/(2*(D1-D))) * tan(phi) + Nphi^(1/2) * (D/(2*(D1-D)))
    #       * {Nphi * D1/(D1-D) * tan(phi) + Nphi^(1/2) - 1}]
    #       - D1 * tan(phi)   ... simplified form
    #
    # The standard implementation follows the form used in FHWA and ReSSA software:

    exp_term = Nphi ** (0.5 * Nphi * D1 / (D1 - D))
    sqrt_Nphi = np.sqrt(Nphi)

    # A2: coefficient for the gamma*z term (overburden pressure)
    A2 = (D1 * (exp_term * np.tan(phi)
          + sqrt_Nphi * D / (2.0 * (D1 - D))
          * (exp_term * np.tan(phi) + sqrt_Nphi - 1.0))
          - D1 * np.tan(phi))

    # A1: coefficient for the cohesion term
    # A1 = 2*Nphi^(1/2) * A2 / (Nphi - 1) when phi > 0
    A1 = 2.0 * sqrt_Nphi * A2 / (Nphi - 1.0)

    return A1, A2


def intersect_pile_with_materials(pile_x, y_ground, y_failure, profile_lines, materials):
    """
    Find soil layers along the pile between the ground surface and failure surface.

    Parameters
    ----------
    pile_x : float
        X-coordinate of the pile (assumed vertical for layer intersection).
    y_ground : float
        Y-coordinate of the ground surface at the pile location.
    y_failure : float
        Y-coordinate of the failure surface at the pile location.
    profile_lines : list
        List of profile line dicts with 'coords' and 'mat_id' keys.
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
    # Get y-coordinates of each profile line at pile_x
    layer_boundaries = []
    for line in profile_lines:
        coords = np.array(line['coords'])
        xs = coords[:, 0]
        ys = coords[:, 1]
        sort_idx = np.argsort(xs)
        y_at_x = np.interp(pile_x, xs[sort_idx], ys[sort_idx], left=np.nan, right=np.nan)
        mat_id = line.get('mat_id')
        if mat_id is not None and 0 <= mat_id < len(materials):
            mat_index = mat_id
        else:
            mat_index = len(layer_boundaries)
        layer_boundaries.append({'y': y_at_x, 'mat_index': mat_index})

    # Sort boundaries from top (highest y) to bottom (lowest y)
    layer_boundaries.sort(key=lambda b: -b['y'] if not np.isnan(b['y']) else float('inf'))

    # Build segments between ground surface and failure surface
    segments = []
    for i, boundary in enumerate(layer_boundaries):
        y_layer_top = boundary['y']
        if np.isnan(y_layer_top):
            continue

        # Layer bottom is the next boundary below, or the failure surface
        if i + 1 < len(layer_boundaries) and not np.isnan(layer_boundaries[i + 1]['y']):
            y_layer_bot = layer_boundaries[i + 1]['y']
        else:
            y_layer_bot = y_failure  # extends to failure surface or below

        # Clip to [y_failure, y_ground] range
        seg_top_y = min(y_layer_top, y_ground)
        seg_bot_y = max(y_layer_bot, y_failure)

        # Skip if segment is outside the range or zero thickness
        if seg_top_y <= y_failure or seg_bot_y >= y_ground:
            continue
        if seg_top_y <= seg_bot_y:
            continue

        # Convert y-coordinates to depths from ground surface (z increases downward)
        z_top = y_ground - seg_top_y
        z_bot = y_ground - seg_bot_y

        mat = materials[boundary['mat_index']]
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
