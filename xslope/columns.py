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
    computed : bool
        True when the column is written by the REPORT rather than by the solver —
        the per-slice terms of the factor of safety equation, which exist only in
        the table the calculations section builds. A solved ``slice_df`` is not
        expected to carry one.
    gated_by : str
        The column whose value decides whether this one is printed, for a column
        that describes a force rather than carrying it. ``drop_empty`` asks
        whether a column is all zeros, which is the right question for a force
        and the wrong one for the angle it acts at: the angle is geometry, stored
        for every slice whether or not a load is applied there, so it survives
        the test on its own and prints an inclination for a force that is not on
        the model. A gated column is printed only when the column it is gated by
        carries a value, so the pair appears and disappears together.
    """

    key: str
    label: str
    description: str
    quantity: str = ""
    fmt: str = "{:.2f}"
    report: bool = False
    computed: bool = False
    gated_by: str = ""


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
    # D, not Q: every derivation page and every equation the report prints calls
    # the distributed-load resultant D, and Q is Spencer's interslice resultant.
    # A table that labelled this Q left the reader of "D cos β" with no column to
    # find it in and a Q that meant something else two columns over.
    Column("dload", "D", "Resultant of the distributed load acting on the top of "
           "the slice, per unit thickness.", "force_per_len", "{:.1f}", True),
    # Gated on D. This angle is the inclination of the slice's top edge, which
    # the slicer computes for every slice whether a distributed load is applied
    # there or not, so it is not zero on a model that carries none — a model
    # loaded by a line load alone printed β = 84.55° beside a D column that had
    # been dropped for being zero, an inclination belonging to no force.
    Column("beta", "β", "Inclination of the distributed-load resultant from the "
           "horizontal.", "deg", "{:.2f}", True, gated_by="dload"),
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

    # The per-slice terms of the factor of safety equation, written by the report's
    # Calculations section (:func:`xslope.report.calculations`) for the method it
    # documents. They are in the table so that the summed values printed in that
    # section can be checked term by term against the slices they came from — which
    # is the whole point of showing the calculation. A method contributes only its
    # own terms, so no slice table carries them all.
    Column("m_res", "M_R", "Resisting moment this slice contributes about the "
           "center of rotation, (c·Δl + N'·tan φ)·a_S.", "moment", "{:.1f}",
           True, True),
    Column("m_drv", "M_D", "Net driving moment this slice contributes about the "
           "center of rotation.", "moment", "{:.1f}", True, True),
    Column("f_res", "F_R", "Resisting force this slice contributes to the "
           "horizontal balance, (c·Δl + N'·tan φ)·cos α.", "force_per_len",
           "{:.1f}", True, True),
    Column("f_drv", "F_D", "Net driving force this slice contributes to the "
           "horizontal balance.", "force_per_len", "{:.1f}", True, True),
    # Spencer's two force sums, ahead of the Q they are the inputs to, so that a
    # reader checking a row reads it left to right in the order the equation
    # uses it.
    Column("F_h", "F_h", "Resultant of all the forces on the slice except the "
           "base normal, the base shear and the interslice forces — horizontal "
           "component, per unit thickness.", "force_per_len", "{:.1f}",
           True, True),
    Column("F_v", "F_v", "Resultant of all the forces on the slice except the "
           "base normal, the base shear and the interslice forces — vertical "
           "component, per unit thickness.", "force_per_len", "{:.1f}",
           True, True),
    Column("q_s", "Q_s", "Resultant of the interslice forces on the slice "
           "(Spencer's Q), per unit thickness.", "force_per_len", "{:.1f}",
           True, True),
    Column("y_q", "y_Q", "Elevation of the line of action of Q_s.", "length",
           "{:.2f}", True, True),

    # Declared, but not printed: geometry the plots already carry, per-layer
    # heights whose count varies with the model, and the solvers' scratch space.
    Column("x_l", "x_L", "Horizontal coordinate of the slice's left edge.", "length"),
    Column("x_r", "x_R", "Horizontal coordinate of the slice's right edge.", "length"),
    Column("y_lb", "y_Lb", "Elevation of the failure surface at the left edge.", "length"),
    Column("y_lt", "y_Lt", "Elevation of the ground surface at the left edge.", "length"),
    Column("y_rb", "y_Rb", "Elevation of the failure surface at the right edge.", "length"),
    Column("y_rt", "y_Rt", "Elevation of the ground surface at the right edge.", "length"),
    Column("y_cg", "y_cg", "Elevation of the slice's center of gravity.", "length"),
    Column("qL", "q_L", "Distributed-load intensity at the slice's left edge.", "stress"),
    Column("qR", "q_R", "Distributed-load intensity at the slice's right edge.", "stress"),
    Column("d_x", "d_x", "Horizontal coordinate of the distributed-load resultant.", "length"),
    Column("d_y", "d_y", "Elevation of the distributed-load resultant.", "length"),
    Column("y_t", "y_Tc", "Elevation of the tension-crack water force's line of action.", "length"),
    Column("piezo_y", "y_piezo", "Elevation of the piezometric surface at the slice mid-point.", "length"),
    Column("theta", "θ", "Interslice force inclination on the right-hand side of the slice.", "deg"),
    Column("c_suction", "c_suc", "Apparent cohesion contributed by matric suction.", "stress"),
    Column("r", "R", "Radius of the trial circle (circular surfaces only).", "length"),
    Column("xo", "X_o", "Horizontal coordinate of the trial circle's center.", "length"),
    Column("yo", "Y_o", "Elevation of the trial circle's center.", "length"),
)

#: Keyed by column name, for a direct lookup.
BY_KEY = {c.key: c for c in SLICE_COLUMNS}


def report_columns():
    """The columns a report's slice table prints, in table order."""
    return tuple(c for c in SLICE_COLUMNS if c.report)


def solver_columns():
    """The printed columns a solved ``slice_df`` is expected to carry."""
    return tuple(c for c in report_columns() if not c.computed)


def computed_columns():
    """The printed columns the report itself writes — the per-slice terms of the
    factor of safety equation."""
    return tuple(c for c in report_columns() if c.computed)


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
        # A moment per unit width of section, spelled the one way the report
        # spells it — force x length per length, as the pile capacities are
        # (:func:`xslope.fem_details.unit_labels`). "lb/ft·ft" is the same
        # quantity written so that it reads as force per area.
        length = unit_labels.get("length", "")
        force = (unit_labels.get("force_per_len", "") or "").split("/")[0]
        return f"{force}·{length} per {length}" if length and force else ""
    return unit_labels.get(q, "")


def header(column, unit_labels=None):
    """``"W (lb/ft)"`` — the column's label with its unit, when it has one."""
    unit = unit_label(column, unit_labels)
    return f"{column.label} ({unit})" if unit else column.label


# ---------------------------------------------------------------------------
# How a table cell reads
#
# A report table's justification is not written table by table. It is read off
# the cells, by one policy applied to every table the report builds: a column of
# numbers is centered, and everything else begins at the left, where a reader's
# eye starts a line. The one other centered case is a column of symbols, which is
# what a nomenclature's first column is and what nothing else in the report is.
#
# The classification lives here because both ends of the pipeline need it — the
# report to justify a column, the renderer to keep a number off two lines — and
# because a policy stated twice is a policy that will disagree with itself.
# ---------------------------------------------------------------------------

def is_number(text):
    """Is this cell a number and nothing else?

    Read off the printed string, not the value behind it: by the time a cell
    exists the value has been through :func:`format_value` and the string is all
    there is. A blank is not a number — a column of forces with an empty cell in
    it is still a column of forces.

    A number written with thousands separators is a number. "1,550.0" is what
    the report prints a force as, and reading it as prose cost it everything a
    number is given: its cell was free to break mid-number across two lines, and
    its column was measured as though it could.
    """
    s = str(text).strip()
    if not s:
        return False
    try:
        float(s.replace(",", "") if "," in s else s)
    except (TypeError, ValueError):
        return False
    return True


#: The longest a cell can be and still be read as a symbol rather than as a word.
#: "N' (lb/ft)" is the longest the report writes.
SYMBOL_MAX_CHARS = 16


def is_symbol(text):
    """Is this cell a symbol — one term of an equation, with its unit if it has
    one?

    A symbol is a single character, or carries something no ordinary word does: a
    Greek letter, a prime, a subscript. That is what keeps the test off a column
    of one-word English — "Warning", "Sand" — which is a column of words however
    short its entries are.
    """
    import re

    s = str(text).strip()
    if not s or len(s) > SYMBOL_MAX_CHARS:
        return False
    core = re.sub(r"\s*\([^()]*\)$", "", s)             # a trailing unit
    if not core or any(c.isspace() for c in core):
        return False
    return len(core) == 1 or any(not (c.isascii() and c.isalnum()) for c in core)


def column_alignment(cells):
    """The justification one column's cells ask for — ``"c"`` or ``"l"``.

    Empty cells are ignored: they say nothing about what the column holds.
    """
    values = [str(c).strip() for c in cells]
    values = [v for v in values if v]
    if not values:
        return "l"
    if all(is_number(v) for v in values) or all(is_symbol(v) for v in values):
        return "c"
    return "l"


def infer_alignment(headers, rows):
    """One alignment letter per column, read off the body cells.

    Headers are not consulted. "Factor of safety" is a phrase set over a column
    of numbers, and it is the numbers that decide; the header then takes its own
    column's alignment, so the two can never disagree.
    """
    n = len(headers)
    return [column_alignment([row[j] for row in rows if j < len(row)])
            for j in range(n)]


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


#: How a factor of safety is written everywhere in a report.
FS_FMT = "{:.3f}"

#: Significant digits carried by a summed quantity in the Calculations section.
#:
#: Derived from :data:`FS_FMT`, not chosen. A factor of safety printed to three
#: decimals must be reproducible from the printed sums to within one unit in that
#: last digit — an absolute tolerance of 1e-3, so a relative tolerance of 1e-3/F.
#: A quotient of two values each rounded to ``s`` significant digits carries a
#: relative error of at most 10^(1-s), so the reproduction needs
#: 10^(1-s) <= 1e-3/F. Six digits satisfies that with an order of magnitude in
#: hand for every factor of safety up to F = 100, which is past anything a slope
#: reports. (The reproduction itself is checked, on the printed strings, by
#: ``test/report_check.py``; this constant is why the check can pass.)
SUM_DIGITS = 6


def format_sum(value, digits=SUM_DIGITS):
    """A summed quantity, at the precision the factor of safety can be rebuilt
    from. Fixed-point for values a reader can hold, exponent notation past that."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return ""
    if v != v:
        return ""
    if v == 0:
        return "0"
    from math import floor, log10
    exponent = floor(log10(abs(v)))
    if -4 <= exponent < digits + 3:
        decimals = max(0, digits - 1 - int(exponent))
        text = f"{v:.{decimals}f}"
        # A trailing run of zeros past the decimal point is precision the value
        # does not have; the digits before the point are all significant.
        if "." in text:
            text = text.rstrip("0").rstrip(".")
        return text
    return f"{v:.{digits - 1}e}"


def format_fs(value):
    """A factor of safety, in the one format the whole report uses."""
    try:
        return FS_FMT.format(float(value))
    except (TypeError, ValueError):
        return ""


def format_residual(value):
    """An equilibrium residual: scientific notation, because the number being
    small is the whole statement."""
    try:
        return f"{float(value):.3e}"
    except (TypeError, ValueError):
        return ""


#: The columns whose SUM over the slices is a number, and which a slice table
#: therefore totals at its foot. Extensive quantities only: forces, moments, and
#: the two lengths that partition the surface. A column left out of this set is
#: left out because adding it up gives nothing — the mean elevation of a slice
#: base, the sum of fifteen friction angles — and a total under such a column
#: would be a number a reader has to work out how to disbelieve.
#:
#: The factor of safety terms (``m_res``/``m_drv``, ``f_res``/``f_drv``) are the
#: reason this exists: the calculations section asks the reader to accept a
#: quotient of two sums, and those two sums belong at the foot of the table the
#: terms came from, where they can be checked. Spencer's ``q_s`` is here for the
#: same reason in reverse — its total is the equilibrium residual, and seeing it
#: come out at essentially zero is the method closing its own books.
#:
#: Spencer's ``F_h`` and ``F_v`` are deliberately NOT here. They are per-slice
#: inputs to Q, and no step of the derivation forms either sum: the equilibrium
#: statement is ΣQ = 0, and a total under F_h would be a number with no equation
#: to check it against.
TOTALLED = frozenset({
    "dx", "dl", "w", "dload", "kw", "t", "p", "h_pile", "n_eff",
    "m_res", "m_drv", "f_res", "f_drv", "q_s",
})

#: Columns kept even when they are entirely zero: a zero there is a reading
#: about the slice, not the absence of a feature.
#:
#: ``F_h`` is here for a third reason: on a slope with no horizontal load it is
#: zero on every slice, and it is still an operand of the Q printed beside it. A
#: reader checking a row needs the zero to be on the page. It reaches only the
#: tables that carry it — a slice table with no F_h column is unaffected.
ALWAYS_PRINTED = frozenset({
    "slice #", "x_c", "y_cb", "y_ct", "dx", "dl", "alpha", "w",
    "mat", "c", "phi", "u", "n_eff", "F_h", "F_v",
})


def selected_columns(slice_df, drop_empty=True):
    """The registry columns a slice table prints for this surface, in order.

    Split out so the table and its totals row cannot disagree about which
    columns are in the table.

    A column with a ``gated_by`` is tested on the column that gates it rather
    than on itself: β is the inclination of the distributed load, and the angle
    is stored on every slice whether or not a load acts there. Tested on its own
    values it survives on a model that applies no distributed load at all.
    """
    import numpy as np

    def empty(key):
        if key not in slice_df.columns:
            return True
        values = np.asarray(slice_df[key].values, dtype="object")
        nums = [float(v) for v in values
                if isinstance(v, (int, float, np.floating, np.integer))
                and v == v]
        return not nums or all(abs(v) < 1e-12 for v in nums)

    cols = []
    for c in report_columns():
        if c.key not in slice_df.columns:
            continue
        if drop_empty and c.key not in ALWAYS_PRINTED and empty(c.gated_by or c.key):
            continue
        cols.append(c)
    return cols


def slice_table(slice_df, unit_labels=None, drop_empty=True):
    """The report's slice table: ``(headers, rows, legend)``.

    ``legend`` is one ``(label, description)`` pair per printed column, which is
    what the table's legend paragraph is written from — the definitions travel
    with the numbers instead of living in a manual.

    ``drop_empty`` removes a report-worthy column that is entirely zero or blank
    for this surface, so a model with no seismic load, no distributed load and no
    reinforcement does not print three columns of zeros. Geometry and strength
    columns are always kept: a zero there is a reading, not an absence.

    The totals row is :func:`slice_totals`, kept separate so this function's
    three-value return stays what every caller already unpacks.
    """
    cols = selected_columns(slice_df, drop_empty)
    headers = [header(c, unit_labels) for c in cols]
    rows = [[format_value(c, row[c.key]) for c in cols]
            for _i, row in slice_df.iterrows()]
    legend = [(header(c, unit_labels), c.description) for c in cols]
    return headers, rows, legend


def slice_totals(slice_df, drop_empty=True, label="Total"):
    """The slice table's totals row, aligned with :func:`slice_table`'s columns.

    One cell per printed column: the column's sum where it has one
    (:data:`TOTALLED`), ``label`` in the first cell, and blank everywhere else.
    Formatted with the column's own format, so a total reads in the same units
    and to the same precision as the values above it.
    """
    import numpy as np

    cols = selected_columns(slice_df, drop_empty)
    if not cols:
        return []
    out = []
    for i, c in enumerate(cols):
        if i == 0:
            out.append(label)
        elif c.key in TOTALLED:
            values = np.asarray(slice_df[c.key].values, dtype=float)
            total = float(np.nansum(values))
            out.append(format_value(c, total))
        else:
            out.append("")
    return out
