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

"""The force at each crossing of a reinforcement line by a failure surface.

:func:`xslope.slice.generate_slices` records every point where the surface
crosses a reinforcement line on ``slice_df.attrs['reinf_crossings']``: the slice
whose base it falls on, the point, the tension the line's capacity envelope
makes available there, and which term of the envelope sets it. This module
turns those records into the one table the report and Studio both print, so the
two show the same rows in the same words.

A passive line is mobilized with the soil, and the force the solvers apply for
it is the available tension divided by the factor of safety. The table states
both numbers, so a reader can see the division.
"""

#: Column headings, before the unit suffix the caller adds to the quantities.
HEADERS = ("Line", "Slice", "x", "y", "Appl", "Available T", "Applied force",
           "Limited by")

#: Which quantity each heading carries, for the unit suffix.
QUANTITIES = (None, None, "length", "length", None, "force_per_len",
              "force_per_len", None)

#: What each heading means, for a table legend or a tooltip.
DEFINITIONS = (
    "The reinforcement line, by its label.",
    "The slice whose base the crossing is on, numbered from the left edge of "
    "the failure surface as in the slice table.",
    "Horizontal coordinate of the crossing.",
    "Elevation of the crossing.",
    "Active lines act at their full available tension. Passive lines are "
    "mobilized with the soil, so their force is divided by the factor of safety.",
    "The tension the line's capacity envelope allows at the crossing, per unit "
    "width.",
    "The force applied to the sliding mass at the crossing, per unit width: the "
    "available tension for an active line, and the available tension divided "
    "by the factor of safety for a passive one.",
    "Which term of the capacity envelope sets the available tension: Tmax, the "
    "line's own capacity, or pullout from an end, where the pullout resistance "
    "between that end and the crossing, plus any end anchorage Tend, develops "
    "less than Tmax.",
)


def limit_text(limit, line=None):
    """The words for one envelope term (:func:`xslope.fileio.reinforce_limiting_term`).

    ``line`` is the model's own reinforcement line, read for its end anchorage:
    an end that develops its force partly through an anchorage says so.
    """
    line = line or {}

    def end(n):
        tend = line.get(f"tend{n}") or 0.0
        try:
            anchored = float(tend) > 0.0
        except (TypeError, ValueError):
            anchored = False
        return (f"Tend{n} + pullout from end {n}" if anchored
                else f"Pullout from end {n}")

    if limit == "tmax":
        return "Tmax"
    if limit == "end1":
        return end(1)
    if limit == "end2":
        return end(2)
    if limit == "ends":
        return "Pullout from both ends"
    return ""


def line_names(slope_data):
    """The display name of every reinforcement line, in the order defined: the
    names the report's reinforcement table and the FEM details print
    (:func:`xslope.fem_details.display_labels`)."""
    from .fem_details import display_labels
    # The rule generate_slices reads the lines by: an explicit list is the
    # model's statement, the empty list included, and the derived point lists
    # stand in only for a model that carries no such key.
    lines = slope_data.get("reinforcement_lines")
    if lines is not None:
        return display_labels([r.get("label") for r in lines], "Line")
    legacy = slope_data.get("reinforce_lines") or []
    return [f"Line {i + 1}" for i in range(len(legacy))]


def crossings(slice_df, slope_data, fs=None):
    """The crossings of one solved surface, or ``None`` for a model with no
    reinforcement lines.

    Returns ``{"rows": [...], "uncrossed": [names], "passive": bool}``. Each row
    is a dict with the line's name and index, the slice number, the point
    (``x``, ``y``), ``appl``, the available tension ``T``, the ``applied`` force
    (``T``, or ``T / fs`` for a passive line; ``None`` for a passive line when
    ``fs`` is not known), and ``limit``, the words for the term that sets ``T``.
    Rows are ordered by line, then along the surface.
    """
    names = line_names(slope_data)
    attrs = getattr(slice_df, "attrs", None) or {}
    records = attrs.get("reinf_crossings")
    if not names or records is None:
        return None
    lines = slope_data.get("reinforcement_lines") or []
    try:
        fs = float(fs) if fs is not None else None
    except (TypeError, ValueError):
        fs = None
    if fs is not None and not fs > 0.0:
        fs = None
    rows = []
    for rec in sorted(records, key=lambda r: (r["line"], r["slice"])):
        idx = rec["line"]
        line = lines[idx] if 0 <= idx < len(lines) else {}
        passive = rec.get("appl") == "passive"
        T = float(rec.get("T") or 0.0)
        applied = (T / fs if fs is not None else None) if passive else T
        rows.append({
            "name": names[idx] if 0 <= idx < len(names) else f"Line {idx + 1}",
            "index": idx, "slice": int(rec["slice"]),
            "x": float(rec["x"]), "y": float(rec["y"]),
            "appl": "passive" if passive else "active",
            "T": T, "applied": applied,
            "limit": limit_text(rec.get("limit"), line),
        })
    crossed = {r["index"] for r in rows}
    uncrossed = [n for i, n in enumerate(names) if i not in crossed]
    return {"rows": rows, "uncrossed": uncrossed,
            "passive": any(r["appl"] == "passive" for r in rows)}


def crossing_table(slice_df, slope_data, fs=None, units=None):
    """``(headers, rows, found)`` for printing, or ``None`` without reinforcement.

    ``units`` is :func:`xslope.units.labels` for the model; the quantities'
    headings carry its units. ``rows`` are lists of strings; ``found`` is the
    dict :func:`crossings` returned, for the sentences around the table.
    """
    found = crossings(slice_df, slope_data, fs)
    if found is None:
        return None
    units = units or {}
    headers = []
    for head, qty in zip(HEADERS, QUANTITIES):
        unit = units.get(qty) if qty else ""
        headers.append(f"{head} ({unit})" if unit else head)
    rows = []
    for r in found["rows"]:
        rows.append([
            r["name"], str(r["slice"]), f"{r['x']:.2f}", f"{r['y']:.2f}",
            r["appl"].capitalize(), f"{r['T']:,.1f}",
            "" if r["applied"] is None else f"{r['applied']:,.1f}",
            r["limit"],
        ])
    return headers, rows, found
