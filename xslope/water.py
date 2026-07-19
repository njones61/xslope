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

"""Water-load synthesis shared by the vendor importers.

Ponded water is the one thing every commercial slope-stability program stores
implicitly — SLOPE/W, Slide2 and RS2 all derive the weight of water standing on
the slope from the water surface, never as an object of its own. xslope needs it
as an explicit distributed load, or the reservoir is simply lost from the
statics. This module holds the one geometric operation that recovers it, so
geostudio.py, slide2.py and rs2.py all convert water the same way:

  - material_above_ground_dload(ground, upper_line, unit_weight): the weight of
    whatever fills the gap between a line and the ground, as dload blocks. Both
    ponded water (upper_line = the water surface, unit_weight = gamma_w) and a
    GeoStudio-style surcharge (upper_line = the fill's top, unit_weight = the
    fill's unit weight) are this same object.
  - ponded_water_dload(ground, piezo_line, gamma_water): the ponded-water name
    for the above, for readers where the intent is only ever water.
  - _y_on(line, x): elevation of a left-to-right polyline at x.

The dload convention is xslope's own: each block is a list of ``{'X','Y','Normal'}``
points along the loaded surface, ``Normal`` the pressure (perpendicular
intensity) at that point, tapering to zero at a water/ground crossing so the load
ends where the water does. It was validated exact against the dry-buoyant oracle
(a fully submerged still-water slope reads identically to the same slope run dry
with gamma' = gamma_sat - gamma_w).
"""

from shapely.geometry import LineString


def _y_on(line, x):
    """Elevation of a left-to-right polyline at x, or None if x is off its ends."""
    if line is None or line.is_empty:
        return None
    coords = list(line.coords)
    if not coords or x < coords[0][0] or x > coords[-1][0]:
        return None
    for (x1, y1), (x2, y2) in zip(coords, coords[1:]):
        if x1 <= x <= x2:
            if x2 == x1:
                return max(y1, y2)
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    return None


def material_above_ground_dload(ground_surface, upper_line, unit_weight):
    """The weight of whatever fills the gap between a line and the ground, as dloads.

    Two GeoStudio things are this same object, and both would otherwise be lost:

    * **Ponded water.** GeoStudio has no ponded-water object -- where the piezometric
      line rises above the ground, SLOPE/W simply carries the water's weight. Pass the
      piezo line and gamma_w.
    * **A surcharge.** GeoStudio's ``<Surcharge>`` is a body of fill between a drawn
      line and the ground, and its ``<Pressure>`` is the fill's UNIT WEIGHT, not a
      pressure -- so the load varies with how deep the fill is. (Verified against
      SLOPE/W's own per-slice surcharge forces: a 20 kN/m3 fill under a line running
      1 m to 2 m above flat ground gives 300 kN/m, not 20 x 10 = 200.)

    Returns a list of dload blocks -- one per stretch where the line is above ground --
    each a list of ``{'X','Y','Normal'}``, with ``Normal`` the pressure
    ``unit_weight * depth``, zero at the crossing so the load tapers out correctly.
    """
    if not upper_line or ground_surface is None or ground_surface.is_empty:
        return []
    piezo = LineString(upper_line)
    gamma_water = unit_weight

    # Sample where either line has a vertex, so no break in either is missed.
    xs = sorted({x for x, _ in ground_surface.coords} | {x for x, _ in piezo.coords})
    samples = []
    for x in xs:
        yg, yp = _y_on(ground_surface, x), _y_on(piezo, x)
        if yg is not None and yp is not None:
            samples.append((x, yg, yp - yg))          # depth of water over the ground
    if not samples:
        return []

    blocks, run = [], []
    for i, (x, yg, d) in enumerate(samples):
        if d > 1e-9:
            if not run and i > 0:                     # entering the water: add the waterline
                px, pyg, pd = samples[i - 1]
                t = -pd / (d - pd) if (d - pd) else 0.0
                xw = px + t * (x - px)
                run.append({"X": xw, "Y": _y_on(ground_surface, xw) or pyg, "Normal": 0.0})
            run.append({"X": x, "Y": yg, "Normal": gamma_water * d})
        elif run:                                     # leaving the water: close at the waterline
            px, pyg, pd = samples[i - 1]
            t = pd / (pd - d) if (pd - d) else 0.0
            xw = px + t * (x - px)
            run.append({"X": xw, "Y": _y_on(ground_surface, xw) or yg, "Normal": 0.0})
            if len(run) >= 2:
                blocks.append(run)
            run = []
    if len(run) >= 2:
        blocks.append(run)
    return blocks


def ponded_water_dload(ground_surface, piezo_line, gamma_water):
    """Water standing above the ground surface, as a distributed load.

    A thin name for :func:`material_above_ground_dload` -- ponded water is just the
    material above the ground being water.
    """
    return material_above_ground_dload(ground_surface, piezo_line, gamma_water)
