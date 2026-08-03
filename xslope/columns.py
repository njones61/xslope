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

"""What every slice-table column means, in one declaration.

``slice_df`` carries about eighty columns. Most of them are the solvers' own
working storage; a handful are the numbers an engineer checks a hand calculation
against. Nothing said which was which, so the report had no basis for choosing
what to print and the documentation had no public definition of ``n_eff``.

This module is that declaration. Each :class:`Column` states the label a reader
sees, the sentence that defines it, the physical quantity it carries (so its
header can be labelled in the model's own units), how it is formatted, and
whether it belongs in a report at all. The report's slice table and its legend
are both generated from it, and a documentation page can generate a column
reference from the same rows.

The quantity vocabulary is :func:`xslope.units.labels`' own key set, plus:

* ``"deg"`` — an angle, always degrees (the package's internal convention);
* ``"moment"`` — a moment per unit width, formed from the declared units;
* ``""`` — dimensionless, or an index/count.

Nothing here computes anything. It reads a column and says what it is.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Column:
    """One ``slice_df`` column, declared.

    Attributes
    ----------
    key : str
        The column name in ``slice_df``.
    label : str
        The heading a reader sees. Kept close to the symbol used in the method
        derivations so a table can be checked against the equations.
    description : str
        One sentence, for the table legend and the documentation.
    quantity : str
        Which physical quantity it carries — a :func:`xslope.units.labels` key,
        ``"deg"``, ``"moment"``, or ``""`` for dimensionless.
    fmt : str
        A ``str.format`` spec for the value.
    report : bool
        True when the column belongs in a report's slice table. Everything else
        is solver working storage, reachable from Python but not printed.
    """

    key: str
    label: str
    description: str
    quantity: str = ""
    fmt: str = "{:.2f}"
    report: bool = False


#: Every declared column, in the order a slice table prints them. Columns absent
#: from this registry are internal by definition — a new solver column shows up
#: nowhere until someone declares what it means.
SLICE_COLUMNS = (
    Column("slice #", "Slice", "Slice number, counted from the left edge of the "
           "failure surface.", "", "{:.0f}", True),
    Column("x_c", "x_c", "Horizontal coordinate of the slice mid-point.",
           "length", "{:.2f}", True),
    Column("y_cb", "y_b", "Elevation of the failure surface at the slice "
           "mid-point.", "length", "{:.2f}", True),
    Column("y_ct", "y_t", "Elevation of the ground surface at the slice "
           "mid-point.", "length", "{:.2f}", True),
    Column("dx", "Δx", "Slice width.", "length", "{:.2f}", True),
    Column("dl", "Δl", "Length of the slice base along the failure surface.",
           "length", "{:.2f}", True),
    Column("alpha", "α", "Inclination of the slice base, positive when the base "
           "rises toward the toe.", "deg", "{:.2f}", True),
    Column("w", "W", "Weight of the slice, per unit thickness.",
           "force_per_len", "{:.1f}", True),
    Column("dload", "Q", "Resultant of the distributed load acting on the top of "
           "the slice, per unit thickness.", "force_per_len", "{:.1f}", True),
    Column("beta", "β", "Inclination of the distributed-load resultant from the "
           "horizontal.", "deg", "{:.2f}", True),
    Column("kw", "kW", "Horizontal seismic force on the slice, per unit "
           "thickness.", "force_per_len", "{:.1f}", True),
    Column("t", "T_c", "Water force in the tension crack, per unit thickness.",
           "force_per_len", "{:.1f}", True),
    Column("p", "P", "Reinforcement force acting on the slice base, resolved "
           "tangent to the failure surface, per unit thickness.",
           "force_per_len", "{:.1f}", True),
    Column("h_pile", "H_p", "Pile resistance mobilised at the slice base, per "
           "unit thickness.", "force_per_len", "{:.1f}", True),
    Column("hw", "h_w", "Height of water above the slice base at its mid-point.",
           "length", "{:.2f}", True),
    Column("u", "u", "Pore water pressure on the slice base at its mid-point.",
           "stress", "{:.1f}", True),
    Column("mat", "Mat", "Index of the material the slice base runs through, "
           "numbered as on the materials sheet.", "", "{:.0f}", True),
    Column("c", "c", "Cohesion of the base material, in the strength "
           "formulation the analysis is run under.", "stress", "{:.1f}", True),
    Column("phi", "φ", "Friction angle of the base material.", "deg",
           "{:.2f}", True),
    Column("n_eff", "N'", "Effective normal force on the slice base, per unit "
           "thickness, as solved.", "force_per_len", "{:.1f}", True),
    Column("z", "Z", "Interslice force on the right-hand side of the slice, per "
           "unit thickness.", "force_per_len", "{:.1f}", True),

    # Declared, but not printed: geometry the plots already carry, per-layer
    # heights whose count varies with the model, and the solvers' scratch space.
    Column("x_l", "x_L", "Horizontal coordinate of the slice's left edge.", "length"),
    Column("x_r", "x_R", "Horizontal coordinate of the slice's right edge.", "length"),
    Column("y_lb", "y_Lb", "Elevation of the failure surface at the left edge.", "length"),
    Column("y_lt", "y_Lt", "Elevation of the ground surface at the left edge.", "length"),
    Column("y_rb", "y_Rb", "Elevation of the failure surface at the right edge.", "length"),
    Column("y_rt", "y_Rt", "Elevation of the ground surface at the right edge.", "length"),
    Column("y_cg", "y_cg", "Elevation of the slice's centre of gravity.", "length"),
    Column("qL", "q_L", "Distributed-load intensity at the slice's left edge.", "stress"),
    Column("qR", "q_R", "Distributed-load intensity at the slice's right edge.", "stress"),
    Column("d_x", "d_x", "Horizontal coordinate of the distributed-load resultant.", "length"),
    Column("d_y", "d_y", "Elevation of the distributed-load resultant.", "length"),
    Column("y_t", "y_Tc", "Elevation of the tension-crack water force's line of action.", "length"),
    Column("piezo_y", "y_piezo", "Elevation of the piezometric surface at the slice mid-point.", "length"),
    Column("theta", "θ", "Interslice force inclination on the right-hand side of the slice.", "deg"),
    Column("c_suction", "c_suc", "Apparent cohesion contributed by matric suction.", "stress"),
    Column("r", "R", "Radius of the trial circle (circular surfaces only).", "length"),
    Column("xo", "X_o", "Horizontal coordinate of the trial circle's centre.", "length"),
    Column("yo", "Y_o", "Elevation of the trial circle's centre.", "length"),
)

#: Keyed by column name, for a direct lookup.
BY_KEY = {c.key: c for c in SLICE_COLUMNS}


def report_columns():
    """The columns a report's slice table prints, in table order."""
    return tuple(c for c in SLICE_COLUMNS if c.report)


def unit_label(column, unit_labels):
    """The unit string for one column, or ``""`` when it has none to state.

    ``unit_labels`` is a :func:`xslope.units.labels` dict, or None for a model
    that declares no unit system — in which case every header stays bare, the
    same way the plots do.
    """
    q = column.quantity
    if not q:
        return ""
    if q == "deg":
        return "deg"
    if not unit_labels:
        return ""
    if q == "moment":
        length, force = unit_labels.get("length", ""), unit_labels.get("force_per_len", "")
        return f"{force}·{length}" if length and force else ""
    return unit_labels.get(q, "")


def header(column, unit_labels=None):
    """``"W (lb/ft)"`` — the column's label with its unit, when it has one."""
    unit = unit_label(column, unit_labels)
    return f"{column.label} ({unit})" if unit else column.label


def format_value(column, value):
    """One cell, formatted by the column's own spec. Blank for a missing value."""
    if value is None:
        return ""
    try:
        if value != value:                       # NaN
            return ""
    except TypeError:
        return str(value)
    try:
        return column.fmt.format(float(value))
    except (TypeError, ValueError):
        return str(value)


def slice_table(slice_df, unit_labels=None, drop_empty=True):
    """The report's slice table: ``(headers, rows, legend)``.

    ``legend`` is one ``(label, description)`` pair per printed column, which is
    what the table's legend paragraph is written from — the definitions travel
    with the numbers instead of living in a manual.

    ``drop_empty`` removes a report-worthy column that is entirely zero or blank
    for this surface, so a model with no seismic load, no distributed load and no
    reinforcement does not print three columns of zeros. Geometry and strength
    columns are always kept: a zero there is a reading, not an absence.
    """
    import numpy as np

    always = {"slice #", "x_c", "y_cb", "y_ct", "dx", "dl", "alpha", "w",
              "mat", "c", "phi", "u", "n_eff"}
    cols = []
    for c in report_columns():
        if c.key not in slice_df.columns:
            continue
        if drop_empty and c.key not in always:
            values = np.asarray(slice_df[c.key].values, dtype="object")
            nums = [float(v) for v in values
                    if isinstance(v, (int, float, np.floating, np.integer))
                    and v == v]
            if not nums or all(abs(v) < 1e-12 for v in nums):
                continue
        cols.append(c)

    headers = [header(c, unit_labels) for c in cols]
    rows = [[format_value(c, row[c.key]) for c in cols]
            for _i, row in slice_df.iterrows()]
    legend = [(header(c, unit_labels), c.description) for c in cols]
    return headers, rows, legend
