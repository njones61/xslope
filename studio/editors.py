"""Input editors — modal dialogs launched from the Inputs tree (Phase 2).

Each editable ``slope_data`` category has a small ``CategoryEditor`` that builds a
dialog from the current data and writes the edited values back. Most categories
reuse the generic ``TableEditorDialog`` (a column spec + add/remove rows);
scalars use ``FormEditorDialog``. Editors only expose the commonly-edited fields
of each record and **preserve any unshown keys** (e.g. material reliability
sigmas), so a round-trip through an editor never drops data.
"""

from __future__ import annotations

import copy
import math

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor, QIcon, QKeySequence, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QButtonGroup, QCheckBox, QComboBox, QDialog,
    QDialogButtonBox, QFormLayout, QFrame, QGridLayout, QGroupBox, QHBoxLayout, QHeaderView,
    QLabel, QLineEdit, QListWidget, QListWidgetItem, QMessageBox, QPushButton,
    QScrollArea, QSplitter, QStackedWidget, QTableWidget, QTableWidgetItem,
    QTabWidget, QVBoxLayout, QWidget,
)

# Point-to-polyline distance shared with the Inputs-canvas hit-tester — reused by the
# preview-pane pick resolvers so a click resolves against the same geometry.
from .picking import _line_dist
# The polygon-overlay vocabulary lives with the loader; the Studio imports it rather
# than restating it, so the Type words, the sentinel codes and their display strings
# can never drift. Same for the reinforcement support-type presets: the table the
# Type column fills Dir/Appl from is the loader's, which is the sheet's.
from xslope.fileio import (POLYGON_TYPE_WORDS, REINFORCE_TYPE_PRESETS,
                           SEARCH_WINDOW_KEYS, SSR_ZONE_LABELS, SSR_ZONE_SENTINELS)

# Column "usage" tags: which analysis a field applies to. Header text is colored
# to mirror the input template's header coloring (red = LEM-specific inputs,
# blue = FEM-specific inputs); seepage/reliability extend the same idea.
USAGE_COLOR = {"lem": "#c00000", "fem": "#0432ff", "seep": "#2e7d32", "rel": "#9a7d0a"}
USAGE_NAME = {"lem": "LEM only", "fem": "FEM only",
              "seep": "Seepage only", "rel": "Reliability"}
# Short labels for the show/hide toggles at the top of an editor.
USAGE_TOGGLE_LABEL = {"lem": "LEM", "fem": "FEM", "seep": "Seepage", "rel": "Reliability"}


# --------------------------------------------------------------------------- #
# Numeric display
#
# Significant digits an editor shows a stored number to. Inputs are typed, computed
# from typed geometry (a circle's R from its intercept point) or generated, and the
# computed ones arrive as full-precision floats: 41.23105625617661 is the radius from
# (10, 40) to (0, 0), and printing it that way is how a coordinate table turns into a
# wall of noise. Six digits is more than any of these inputs carries physically, and
# the stored value is never touched — Field.read_text hands back the number the cell
# was FILLED from whenever the text is still the rendering of it, so a rounding can
# only ever be seen, never saved.
#
# (The transient dialog's own _tseep_fmt renders at 10 digits instead. It edits its
# values in place with no stored original to fall back on, so its display has to be
# lossless; here the original is kept, which is what buys the shorter one.)
#
# The digits are capped; the NOTATION is chosen by magnitude. A plain %g switches to
# an exponent as soon as a value needs more than the significant digits it is allowed
# -- which turns a Young's modulus of 2088500 into 2.0885e+06, an engineering input
# rewritten as physics notation. 101 of the 312 corpus models carry at least one such
# value. Inside the window below, the value is written out in full; outside it -- a
# hydraulic conductivity of 7e-05, a modulus in the billions -- an exponent is what
# anyone would write by hand anyway.
# --------------------------------------------------------------------------- #
_DISPLAY_SIG_DIGITS = 6
#: Magnitudes written in full rather than with an exponent. The lower bound is where
#: leading zeros start to outnumber digits; the upper is where a written-out integer
#: stops being readable at a glance.
_DISPLAY_FIXED_RANGE = (1e-4, 1e9)


def _display_number(v):
    """A stored value rendered for a numeric editor widget.

    ``None`` and NaN (the templates' two spellings of "unset") render blank; ints
    render whole; floats render at :data:`_DISPLAY_SIG_DIGITS` significant digits with
    no trailing zeros (``40.0`` -> ``"40"``), written out in full over
    :data:`_DISPLAY_FIXED_RANGE` (``2088500.0`` -> ``"2088500"``, never
    ``"2.0885e+06"``) and with an exponent outside it. A non-number is left to
    ``str``."""
    if v is None:
        return ""
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return str(v)
    if isinstance(v, float) and v != v:            # NaN = unset, never the text "nan"
        return ""
    if isinstance(v, int):
        return str(v)
    if v == 0:
        return "0"
    low, high = _DISPLAY_FIXED_RANGE
    if not low <= abs(v) < high:
        return f"{v:.{_DISPLAY_SIG_DIGITS}g}"
    # Significant digits, spelled out: the decimals left over once the digits before
    # the point have taken their share (none left over for a value in the millions,
    # which is exactly the case a %g would have put an exponent on).
    decimals = max(0, _DISPLAY_SIG_DIGITS - 1 - math.floor(math.log10(abs(v))))
    text = f"{v:.{decimals}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text or "0"


def _unedited(text, stored):
    """True when ``text`` is still exactly what :func:`_display_number` rendered for
    ``stored`` — the widget is showing a display rounding nobody has typed over, so
    the stored value is what must be kept.

    :meth:`Field.read_text` is this rule for a field in a table or a form; this is
    the same rule for the handful of scalar edits that have no Field behind them (the
    max-depth elevation, a seepage BC's head and flux)."""
    return (isinstance(stored, (int, float)) and not isinstance(stored, bool)
            and text == _display_number(stored))


# --------------------------------------------------------------------------- #
# Field spec + dialogs
# --------------------------------------------------------------------------- #
class Field:
    """One editable column/field. kind ∈ {'float','optfloat','int','str','choice','bool'}.

    'optfloat' is an optional float: a blank cell reads back as None (not 0.0), so
    optional fields like a pile's H or capacities stay unset instead of becoming
    zero (which the loader would reject).

    'bool' is a YES/blank flag, edited through the same combo widget as 'choice' but
    stored as a real Python bool — matching what the loader produces for a YES/blank
    worksheet column. Writing the literal strings back would make the loader reject
    the saved file, so the conversion happens at BOTH ends here."""

    _BLANK = {"float": 0.0, "optfloat": None, "int": 0, "str": "", "choice": "",
              "bool": False}
    # Kinds edited through a QComboBox rather than a line edit.
    COMBO_KINDS = ("choice", "bool")
    # Kinds whose value is a number, and so is display-rounded on the way out and
    # preserved on the way back in (see to_text / read_text).
    NUMERIC_KINDS = ("float", "optfloat", "int")
    _BOOL_CHOICES = ["", "YES"]

    def __init__(self, key, header, kind="float", choices=None, default=None,
                 usage=None, applies=None, tooltip=None, unit=None):
        self.key = key
        self.header = header
        self.kind = kind
        # Optional units-plan (phase 4) tag: the name of a ``xslope.units.labels()``
        # dict key (e.g. "stress", "unit_weight", "k"). When the project declares a
        # unit system the form builder appends the corresponding label to the header
        # (" (kPa)"); None leaves the field unlabeled. Purely display — the stored key
        # and the parsed value are untouched.
        self.unit = unit
        # Optional hover help — shown on the table column header and on the list-view
        # cell label/edit. None (the default) leaves the field un-annotated.
        self.tooltip = tooltip
        if kind == "bool" and not choices:
            choices = self._BOOL_CHOICES
        self.choices = [str(c) for c in (choices or [])]
        # `applies`: the set of analyses this field is used in — drives the
        # show/hide column toggles; None = universal (geometry/identity, always
        # shown). `usage`: the single analysis a field is *specific* to, used for
        # the mirrored header color — only set when the field belongs to exactly
        # one analysis (shared fields like c/φ stay uncolored but still hide when
        # neither of their analyses is toggled on).
        if applies is None and usage is not None:
            applies = {usage}
        self.applies = frozenset(applies) if applies is not None else None
        self.usage = usage if usage is not None else (
            next(iter(self.applies)) if (self.applies and len(self.applies) == 1) else None)
        if default is None:
            default = (self.choices[0] if kind == "choice" and self.choices
                       else self._BLANK[kind])
        self.default = default

    def to_text(self, val):
        """``val`` rendered for this field's editor widget (the inverse of
        :meth:`from_text` for the combo kinds).

        Numeric kinds are rendered at :data:`_DISPLAY_SIG_DIGITS` significant digits
        rather than through ``str``, which prints a float's full repr: a radius the
        editor itself computed from an intercept point is stored as
        41.23105625617661 and has no business being *shown* that way. The stored
        value is untouched — :meth:`read_text` is the other half of that promise."""
        if self.kind == "bool":
            return "YES" if val else ""
        if self.kind in self.NUMERIC_KINDS:
            return _display_number(val)
        return "" if val is None else str(val)

    def read_text(self, text, stored):
        """``text`` parsed back to a value — unless it is still exactly what
        :meth:`to_text` rendered for ``stored``, in which case ``stored`` comes back
        untouched.

        This is what makes a display rounding a *display* rounding. A cell showing
        41.2311 for a stored 41.23105625617661 must save the stored number, not the
        rounding; only text the user actually changed may replace the value. Every
        editor reads its widgets back through here, so no editor can round a value
        into the model by being opened and OK'd."""
        if self.kind in self.NUMERIC_KINDS and _unedited(text, stored):
            return stored
        return self.from_text(text)

    def from_text(self, text):
        text = (text or "").strip()
        if self.kind == "bool":
            return text.lower() in ("yes", "true", "1")
        if self.kind == "float":
            try:
                return float(text)
            except ValueError:
                return 0.0
        if self.kind == "optfloat":
            if text == "" or text.lower() == "none":
                return None
            try:
                return float(text)
            except ValueError:
                return None
        if self.kind == "int":
            try:
                return int(float(text))
            except ValueError:
                return 0
        return text

    def display_header(self, unit_labels=None):
        """The column/field header, with a unit suffix when this field declares a
        ``unit`` and ``unit_labels`` supplies a non-empty label for it.

        ``unit_labels`` is a ``xslope.units.labels()`` dict, or ``None`` for a project
        that declares no unit system. When it is ``None``, or the field has no ``unit``,
        or the label string is empty (e.g. a time-bearing unit with no declared time
        base), the bare header is returned unchanged -- so an undeclared project renders
        byte-identically to today."""
        return _with_unit(self.header, self, unit_labels)


def _with_unit(base, field, unit_labels):
    """``"<base> (<label>)"`` when ``field`` declares a ``unit`` and ``unit_labels`` has
    a non-empty label for it; ``base`` unchanged otherwise (the undeclared/legacy path,
    kept byte-identical to today)."""
    if field is not None and getattr(field, "unit", None) and unit_labels:
        suffix = unit_labels.get(field.unit)
        if suffix:
            return f"{base} ({suffix})"
    return base


def _unit_labels_for(slope_data):
    """The ``xslope.units.labels()`` dict for a project's declared unit system, or
    ``None`` when it declares none.

    Returning ``None`` (rather than the all-empty labels dict) lets every dialog's
    ``unit_labels=None`` default mean "no suffixes" with a single truthiness check, and
    keeps an undeclared project's forms byte-identical to today. When a system is
    declared but no time unit is, the time-bearing labels (k, flowrate) come back empty
    and simply produce no suffix -- the time unit is never guessed."""
    from xslope.units import labels, normalize_unit_system
    system = normalize_unit_system((slope_data or {}).get("unit_system"))
    if system is None:
        return None
    return labels(system, (slope_data or {}).get("time_unit"))


# --------------------------------------------------------------------------- #
# Per-element vs per-unit-width labeling for the reinforcement / pile editors.
#
# The loader divides every reinforcement/pile capacity term (and Area) by the row's
# out-of-plane Spacing (S), so everything internal and every reported output is per
# unit width of slope. The AMBIGUITY is input-side only: a discrete support (a nail
# or anchor, Spacing set) is entered per element and divided by Spacing; a
# geosynthetic (Spacing blank or 1) is entered already per unit width. The list-view
# form makes that explicit per row — the force/capacity/area labels re-word from
# "per element" to "per unit width" as the row's Spacing/S goes set -> blank, and a
# unit string joins the wording only when the project declares a unit system.
# --------------------------------------------------------------------------- #
def _spacing_is_per_element(text):
    """True when a Spacing/S entry means the row's capacities are entered PER ELEMENT.

    Blank (or an explicit ``1``) is the geosynthetic case — the entries are already
    per unit width, so no division happens and the label reads "per unit width".
    Any other positive value is a discrete-support spacing the loader divides by, so
    the entries are per element."""
    t = (text or "").strip()
    if t == "" or t.lower() == "none":
        return False
    try:
        v = float(t)
    except ValueError:
        return False
    return v > 0 and v != 1.0


def _reinf_capacity_unit(unit_labels, quantity, per_element):
    """Unit string for a spacing-scaled reinforcement/pile quantity, or ``""`` when the
    project declares no unit system.

    ``quantity`` is ``"force"`` ([F] / [F/L]), ``"moment"`` ([F·L] / [F·L/L]), or
    ``"area"`` ([L²] / [L²/L]); ``per_element`` picks the per-element vs per-unit-width
    form. Derived from ``xslope.units.labels()``'s ``force_per_len`` (e.g. "kN/m") and
    ``length`` (e.g. "m"), since ``labels()`` carries no bare force/area/moment key --
    a per-unit-width force IS ``force_per_len``, and its numerator is the per-element
    force. The per-unit-width forms keep the ``/L`` explicit (``"m²/m"``, ``"kN·m/m"``)
    rather than dimensionally simplifying, so the label announces the convention."""
    if not unit_labels:
        return ""
    fpl = unit_labels.get("force_per_len", "")      # "kN/m" / "lb/ft"
    length = unit_labels.get("length", "")           # "m" / "ft"
    force = fpl.split("/")[0] if fpl else ""          # "kN" / "lb"
    if quantity == "force":
        return force if per_element else fpl
    if quantity == "moment":
        if not (force and length):
            return ""
        moment = f"{force}·{length}"             # "kN·m"
        return moment if per_element else f"{moment}/{length}"
    if quantity == "area":
        if not length:
            return ""
        return f"{length}²" if per_element else f"{length}²/{length}"
    return ""


def _dynamic_capacity_label(base_header, quantity, unit_labels, per_element):
    """``"<header> (per element[, <unit>])"`` / ``"<header> (per unit width[, <unit>])"``.

    The convention wording always shows; the unit string joins it only when a unit
    system is declared (``unit_labels`` truthy). This is the live label the list-view
    form shows for a spacing-scaled field, rebuilt whenever the row's Spacing changes."""
    convention = "per element" if per_element else "per unit width"
    unit = _reinf_capacity_unit(unit_labels, quantity, per_element)
    inner = f"{convention}, {unit}" if unit else convention
    return f"{base_header} ({inner})"


class FormEditorDialog(QDialog):
    """Simple key/value form for scalar parameters."""

    def __init__(self, title, fields, values, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._fields = fields
        self._edits = {}
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self._values = dict(values)
        for f in fields:
            edit = QLineEdit(f.to_text(values.get(f.key, f.default)))
            if f.kind in Field.NUMERIC_KINDS:
                edit.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self._edits[f.key] = edit
            form.addRow(f.header, edit)
        layout.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    def result_values(self):
        return {f.key: f.read_text(self._edits[f.key].text(),
                                   self._values.get(f.key, f.default))
                for f in self._fields}


# --------------------------------------------------------------------------- #
# Table clipboard — a block of cells in, a block of cells out
#
# The inputs of a slope model are tables of numbers, and they nearly always exist
# somewhere else first: a worksheet, a page of coordinates, a set of vertices out
# of a report. Retyping such a block is where a model gains a digit it was never
# given, so every editor table takes a block from the clipboard and fills itself
# from it — the same block a spreadsheet copies, tabs between columns and newlines
# between rows.
#
# Tabs and newlines ONLY: that is what Excel puts on the clipboard and what a table
# copied out of a browser carries, and runs of spaces are what a column of typed
# numbers is padded with, not what separates them.
# --------------------------------------------------------------------------- #
def _parse_clipboard_block(text):
    """Clipboard ``text`` as a list of rows of cell strings, or ``[]`` for nothing
    pastable.

    Tabs separate columns, newlines separate rows (all three line endings), and
    trailing newlines — which every spreadsheet appends to a copied block — are
    dropped rather than read as an empty last row. Rows are left RAGGED: a row with
    fewer cells than its neighbours fills fewer columns rather than blanking the
    rest, so a block whose last column is empty on some rows leaves those cells
    alone."""
    if not text:
        return []
    body = text.replace("\r\n", "\n").replace("\r", "\n").rstrip("\n")
    if body == "":
        return []
    return [line.split("\t") for line in body.split("\n")]


#: Why a cell in a pasted block was not written, in the order the status line lists
#: them. Each is a thing the block asked for that the table will not do.
_PASTE_SKIP_REASONS = ("past the last column", "not a listed choice", "read-only")


def _paste_skip_phrase(counts):
    """"2 past the last column, 1 read-only" — the skipped cells named by reason,
    or "" when nothing was skipped."""
    parts = [f"{counts[r]} {r}" for r in _PASTE_SKIP_REASONS if counts.get(r)]
    return ", ".join(parts)


def _plural(n, noun):
    return f"{n} {noun}" if n == 1 else f"{n} {noun}s"


def _spreadsheet_manners(table):
    """Cell-level selection and type-to-edit on a QTableWidget: typing into the
    current cell starts the edit, a second click opens it, Tab moves along the
    row. Single click never opens an editor, so block paste stays on the grid."""
    table.setSelectionBehavior(QAbstractItemView.SelectItems)
    table.setEditTriggers(
        QAbstractItemView.AnyKeyPressed | QAbstractItemView.EditKeyPressed
        | QAbstractItemView.DoubleClicked | QAbstractItemView.SelectedClicked)


class _ClipboardTable(QTableWidget):
    """The editors' table widget, with Copy and Paste over a block of cells.

    Subclassed rather than wired to a QShortcut on purpose: a key event reaches the
    focused widget first, so while a cell is open for editing its line edit keeps
    its own Ctrl/Cmd+V (pasting text into the cell being typed in) and this handler
    runs only when the TABLE itself has the keystroke. A shortcut would take the
    keystroke off the line edit and paste a block over the cell instead.

    ``QKeySequence.Paste`` rather than a literal Ctrl+V so the platform's own
    binding is what works — Cmd+V on macOS, Ctrl+V elsewhere."""

    def __init__(self, rows, cols, owner):
        super().__init__(rows, cols)
        self._owner = owner

    def keyPressEvent(self, event):
        if event.matches(QKeySequence.Copy):
            self._owner.copy_selection()
            event.accept()
            return
        if event.matches(QKeySequence.Paste):
            self._owner.paste_clipboard()
            event.accept()
            return
        super().keyPressEvent(event)


class _EditableTable(QWidget):
    """A table over a list of dict records with Add/Remove rows. Unshown keys are
    preserved. Reused standalone (TableEditorDialog) and per-tab (TabbedTableEditorDialog)."""

    def __init__(self, fields, rows, new_row, parent=None, swatch_state=None,
                 on_change=None, on_select=None, dim_rule=None, unit_labels=None,
                 preset_spec=None, dim_on_edit=False):
        super().__init__(parent)
        self._fields = fields
        self._new_row = new_row
        # Optional preset rule: one choice column FILLS others (reinforcement's Type
        # -> Dir/Appl). See _apply_preset_row for what picking one does.
        self._preset = preset_spec or None
        # Units-plan (phase 4) label dict, or None when the project declares no unit
        # system. Only affects the displayed column headers (see below); None keeps
        # them byte-identical to today.
        self._unit_labels = unit_labels
        # Optional live-edit hooks (used by the editor previews): on_change fires on
        # any data edit (cell text, combo, add/remove); on_select on a row-selection
        # change. Suppressed while the table populates so construction is silent.
        self._on_change = on_change
        self._on_select = on_select
        # Optional per-row disable rule: row_dict -> set of field keys to gray out
        # (kept read-only, value retained). The Materials editor uses it to mirror
        # the mat-sheet conditional formatting for an option=elastic row.
        self._dim_rule = dim_rule
        # Whether a TYPED cell can change which cells apply. The materials rule is
        # driven by combos, which re-derive on their own; the reinforcement rule is
        # driven by whether Adhesion and Delta carry values, so it has to re-derive
        # when one is typed or cleared.
        self._dim_on_edit = bool(dim_on_edit)
        self._suppress_notify = True
        self._bases = [dict(r) for r in rows]  # keep originals to preserve extra keys
        # Optional leading display-color swatch column (Materials editor). It is a
        # *display* column, not a data field: kept as the appended LOGICAL column
        # (len(fields)) so every field/column index elsewhere — result_rows,
        # apply_usage_filter, callers indexing by field position — is unchanged, then
        # moved to VISUAL position 0 so it reads as a leading swatch. Slot-keyed via
        # the shared _MaterialColorState; committed to the style delta on OK.
        self._swatch = swatch_state
        ncols = len(fields) + (1 if swatch_state is not None else 0)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.table = _ClipboardTable(len(rows), ncols, self)
        self.table.setHorizontalHeaderLabels(
            [f.display_header(self._unit_labels) for f in fields])
        # Spreadsheet manners: click selects the CELL (not the row — a row-header
        # click still selects whole rows for Remove), typing into the current cell
        # starts the edit, a second click on the current cell opens it, and Tab
        # commits and moves along the row. Block paste stays on the table itself:
        # a single click never opens an editor, so Ctrl+V lands on the grid.
        self.table.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.table.setEditTriggers(
            QAbstractItemView.AnyKeyPressed | QAbstractItemView.EditKeyPressed
            | QAbstractItemView.DoubleClicked | QAbstractItemView.SelectedClicked)
        for j, f in enumerate(fields):  # color usage-tagged headers (red=LEM, blue=FEM, …)
            hi = self.table.horizontalHeaderItem(j)
            if hi is None:
                continue
            if f.usage:
                hi.setForeground(QColor(USAGE_COLOR[f.usage]))
                fnt = hi.font()
                fnt.setBold(True)
                hi.setFont(fnt)
            if f.tooltip:
                hi.setToolTip(f.tooltip)
        if swatch_state is not None:
            col = len(fields)
            hdr = QTableWidgetItem("\U0001F3A8")   # palette glyph — no data header text
            hdr.setToolTip("Display color on the Inputs plot — not a 'mat' worksheet "
                           "column. Click a swatch to change it; pick “Default” "
                           "to reset to the palette color.")
            self.table.setHorizontalHeaderItem(col, hdr)
            self.table.horizontalHeader().moveSection(col, 0)   # -> leading (visual)
        layout.addWidget(self.table)
        for i, row in enumerate(rows):
            self._set_row(i, row)
        if swatch_state is not None:
            self._rebuild_swatches()
        self._apply_dim_all()      # gray inapplicable cells (e.g. an elastic row)
        # Notifications are wired only AFTER the initial population, and item edits
        # go through a suppress flag, so building the table fires nothing.
        self.table.itemChanged.connect(self._on_item_changed)
        self.table.itemSelectionChanged.connect(self._emit_select)
        self._suppress_notify = False

        bar = QHBoxLayout()
        add = QPushButton("Add row")
        add.clicked.connect(self._add_row)
        rem = QPushButton("Remove selected")
        rem.clicked.connect(self._remove_rows)
        bar.addWidget(add)
        bar.addWidget(rem)
        bar.addStretch(1)
        layout.addLayout(bar)
        # What the last paste did, written into the pane under the buttons rather
        # than into a box to dismiss: a paste that filled fewer cells than the block
        # held is something to READ against the rows it landed in. Hidden until
        # there is one, so a dialog that has not been pasted into measures and lays
        # out exactly as it did before.
        self.paste_summary = QLabel()
        self.paste_summary.setObjectName("paste_summary")
        self.paste_summary.setWordWrap(True)
        self.paste_summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.paste_summary.setVisible(False)
        layout.addWidget(self.paste_summary)
        # Every table is fitted to its own content the moment it is populated, so no
        # editor opens on Qt's default section width — which is the same number for
        # an x-coordinate and a support-type combo, and therefore both too wide for
        # one and too narrow for the other.
        self.fit_columns()

    # ---------------------------------------------------------------- #
    # Column fit
    # ---------------------------------------------------------------- #
    def fit_columns(self):
        """Size every column to its own widest content; returns the width the
        VISIBLE columns need in the table's viewport.

        The field kinds are what tell a numeric column (floored at the widest
        ordinary number, so an empty table is as wide as a full one) from a text or
        choice column (floored at its header). A display column that is not a field
        — the Materials swatch — is not numeric either, and floors on its own
        header glyph."""
        numeric = {j for j, f in enumerate(self._fields)
                   if f.kind in Field.NUMERIC_KINDS}
        return _fit_columns(self.table, numeric)

    def columns_width(self):
        """The width this pane needs to show every visible column: the fitted
        columns plus everything they sit inside — the row-number gutter, the frame,
        and the scroll bar the table reserves.

        Measured from the widgets' metrics rather than read off the current
        geometry, so it is the same answer before the pane is ever shown as after.
        sizeHint(), not width(), for the gutter: it is only as wide as the row count
        needs, and before layout its width() has not caught up with the rows in it —
        a 40-row table's gutter is wider than a 1-row table's, and measuring the
        stale one puts a column back over the edge."""
        from PySide6.QtWidgets import QStyle

        return (self.fit_columns()
                + self.table.verticalHeader().sizeHint().width()
                + 2 * self.table.frameWidth()
                + self.style().pixelMetric(QStyle.PM_ScrollBarExtent))

    # ---------------------------------------------------------------- #
    # Clipboard
    # ---------------------------------------------------------------- #
    def _clipboard_columns(self):
        """The data columns a copy or paste walks, left to right AS THE TABLE SHOWS
        THEM.

        Hidden columns are not walked: a usage toggle is the user saying which
        parameters they are working with, and a block pasted over a table with the
        FEM columns folded away belongs to the columns in front of them. The
        Materials swatch is skipped for a different reason — it is a display color,
        not a field, so it is not a column a block can carry. Visual order, because
        that swatch is moved to the front of the table and "rightward" means what
        the eye means by it."""
        header = self.table.horizontalHeader()
        return sorted((j for j in range(len(self._fields))
                       if not self.table.isColumnHidden(j)),
                      key=header.visualIndex)

    def _cell_text(self, i, j):
        """The text column ``j`` of row ``i`` currently shows."""
        w = self.table.cellWidget(i, j)
        if isinstance(w, QComboBox):
            return w.currentText()
        it = self.table.item(i, j)
        return it.text() if it is not None else ""

    def copy_selection(self):
        """The selected block onto the clipboard as tab-separated text — the shape
        this table's own paste reads back, and the shape a spreadsheet takes.

        These tables select whole rows (the selection drives the preview's
        emphasis), so a selection is a run of rows across every column on show; with
        no selection at all, the current cell alone is copied."""
        cols = self._clipboard_columns()
        if not cols:
            return
        chosen = {(ix.row(), ix.column()) for ix in self.table.selectedIndexes()}
        if chosen:
            rows = sorted({r for r, _ in chosen})
            cols = [c for c in cols if any((r, c) in chosen for r in rows)]
        else:
            r, c = self.table.currentRow(), self.table.currentColumn()
            if r < 0 or c not in cols:
                return
            rows, cols = [r], [c]
        block = "\n".join("\t".join(self._cell_text(r, c) for c in cols) for r in rows)
        QApplication.clipboard().setText(block)

    def paste_clipboard(self):
        """Fill the table from a tab-separated block on the clipboard, starting at
        the current cell and running right and down.

        ROWS GROW to fit the block; COLUMNS DO NOT. A block is a set of records, and
        a table that would not lengthen for one could only ever be pasted into after
        the rows had been added by hand — which is most of the typing the paste is
        there to avoid, and the whole of it when the table starts empty. Columns are
        the opposite case: the columns are the fields this input HAS, and a block
        wider than they are is a block that came from somewhere else. Its extra
        columns are dropped and counted rather than shifted onto fields they do not
        name.

        Every value goes in the way a typed one does — the cell's text for a number
        (so the field's own parser reads it, and a pasted rendering of the stored
        value is kept the way an untouched cell is), the matching entry for a choice
        column, matched without regard to case. Text that names no choice, and a
        column held read-only by the row's own state, leave their cells as they
        were and are counted. One change notification for the whole block, so a
        preview redraws for the finished table rather than for every cell of a
        half-filled one."""
        block = _parse_clipboard_block(QApplication.clipboard().text())
        cols = self._clipboard_columns()
        if not block or not cols:
            return
        anchor_col = self.table.currentColumn()
        start = cols.index(anchor_col) if anchor_col in cols else 0
        start_row = max(0, self.table.currentRow())
        skipped = {r: 0 for r in _PASTE_SKIP_REASONS}
        width = 0

        prev = self._suppress_notify
        self._suppress_notify = True
        try:
            grow = start_row + len(block) - self.table.rowCount()
            for _ in range(max(0, grow)):
                i = self.table.rowCount()
                self.table.insertRow(i)
                base = self._new_row()
                self._bases.append(base)
                self._set_row(i, base)
            for dr, cells in enumerate(block):
                for dc, text in enumerate(cells):
                    if start + dc >= len(cols):
                        skipped["past the last column"] += 1
                        continue
                    width = max(width, dc + 1)
                    reason = self._set_cell_text(start_row + dr, cols[start + dc], text)
                    if reason:
                        skipped[reason] += 1
            self._rebuild_swatches()
            self._apply_dim_all()
        finally:
            self._suppress_notify = prev
        self._show_paste_summary(len(block), width, skipped)
        self._emit_change()

    def _set_cell_text(self, i, j, text):
        """Write one pasted cell. Returns "" when it landed, else the reason it did
        not — the wording the status line counts under."""
        f = self._fields[j]
        text = (text or "").strip()
        w = self.table.cellWidget(i, j)
        if isinstance(w, QComboBox):
            if not w.isEnabled():
                return "read-only"
            for k in range(w.count()):
                if w.itemText(k).strip().lower() == text.lower():
                    w.setCurrentIndex(k)
                    # A pasted Type fills Dir/Appl exactly as a picked one does —
                    # called rather than left to the index-changed signal, which does
                    # not fire when the block names the type the cell already holds.
                    # (Idempotent when it did fire.) A block that also carries Dir or
                    # Appl still wins: those columns come after this one, so they are
                    # written over the fill.
                    if self._preset is not None and f.key == self._preset["driver"]:
                        self._apply_preset_row(w)
                    return ""
            return "not a listed choice"
        it = self.table.item(i, j)
        if it is None:
            it = QTableWidgetItem("")
            if f.kind in Field.NUMERIC_KINDS:
                it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.table.setItem(i, j, it)
        if not (it.flags() & Qt.ItemIsEditable):
            return "read-only"
        it.setText(text)
        return ""

    def _show_paste_summary(self, rows, cols, skipped):
        """"Pasted 6 rows × 4 columns" — plus what the block asked for that the
        table did not do, and why."""
        text = f"Pasted {_plural(rows, 'row')} × {_plural(cols, 'column')}"
        phrase = _paste_skip_phrase(skipped)
        if phrase:
            text += f"; {_plural(sum(skipped.values()), 'cell')} skipped ({phrase})"
        self.paste_summary.setText(text + ".")
        self.paste_summary.setVisible(True)

    def _emit_change(self):
        if not self._suppress_notify and self._on_change is not None:
            self._on_change()

    def _emit_select(self):
        if not self._suppress_notify and self._on_select is not None:
            self._on_select()

    def selected_row(self):
        """Index of the currently selected row, or -1 if none."""
        return self.table.currentRow()

    def select_row(self, i):
        """Programmatically select row ``i`` (used by the preview-pane pick). Fires
        the normal selection-changed path, so the preview re-renders with the new
        emphasis — the visual acknowledgement of the pick."""
        if 0 <= i < self.table.rowCount():
            self.table.selectRow(i)

    def column_at(self, widget):
        """Logical column index for a focused ``widget``: the current column when
        the table itself holds focus, else the column of the cell widget that owns
        ``widget`` (a combo, or the swatch button), else the current column. Lets
        the help strip name the field under focus in the table view."""
        tbl = self.table
        if widget is tbl:
            return tbl.currentColumn()
        w = widget
        while w is not None and w is not tbl:
            for r in range(tbl.rowCount()):
                for c in range(tbl.columnCount()):
                    if tbl.cellWidget(r, c) is w:
                        return c
            w = w.parentWidget()
        return tbl.currentColumn()

    def apply_usage_filter(self, enabled):
        """Show only columns whose ``applies`` set intersects ``enabled`` (a set of
        analysis tags). Universal columns (applies=None) are always shown. Hidden
        columns keep their data, so toggling never drops values on save."""
        for j, f in enumerate(self._fields):
            visible = f.applies is None or bool(f.applies & enabled)
            self.table.setColumnHidden(j, not visible)

    def _set_row(self, i, row):
        for j, f in enumerate(self._fields):
            val = row.get(f.key, f.default)
            if f.kind in Field.COMBO_KINDS:
                combo = QComboBox()
                combo.addItems(f.choices)
                _txt = f.to_text(val)
                if _txt in f.choices:
                    combo.setCurrentText(_txt)
                # A preset driver fills its dependent columns BEFORE the general
                # notify, so one edit is one redraw of the finished row. Connected
                # here (after the initial setCurrentText above) so populating a table
                # from a file never re-fills a Dir/Appl the file overrode; `activated`
                # as well as `currentIndexChanged` because re-picking the SAME Type is
                # the user re-asserting the preset over an override, and that fires no
                # index change.
                if self._preset is not None and f.key == self._preset["driver"]:
                    combo.activated.connect(
                        lambda *_, w=combo: self._apply_preset_row(w))
                    combo.currentIndexChanged.connect(
                        lambda *_, w=combo: self._apply_preset_row(w))
                # A combo edit (e.g. the strength 'option') can change which cells are
                # applicable, so re-evaluate the disable rule for every row, then
                # notify. Re-dimming all rows keeps this correct across add/remove
                # without capturing a stale row index.
                combo.currentIndexChanged.connect(lambda *_: self._on_combo_changed())
                self.table.setCellWidget(i, j, combo)
            else:
                it = QTableWidgetItem(f.to_text(val))
                if f.kind in Field.NUMERIC_KINDS:
                    # Numbers right-align, text left — the spreadsheet
                    # convention, so magnitudes line up down a column.
                    it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.table.setItem(i, j, it)

    def replace_rows(self, rows):
        """Swap the whole table for ``rows``, as one undoable edit of the dialog.

        Used by a generator button, which does not add a row to what is there --
        it *replaces* the surface with one derived from the model. Notifications
        are suppressed while the table repopulates and fired once at the end, so
        the preview redraws for the finished set rather than for every row of a
        half-built one."""
        self._suppress_notify = True
        try:
            self.table.setRowCount(0)
            self._bases = [dict(r) for r in rows]
            self.table.setRowCount(len(rows))
            for i, row in enumerate(rows):
                self._set_row(i, row)
            self._rebuild_swatches()
            self._apply_dim_all()
        finally:
            self._suppress_notify = False
        self._emit_change()

    def _rebuild_swatches(self):
        """(Re)build the leading display-color swatch button for every row, bound to
        the row's CURRENT index. Rebuilt after any add/remove so each button targets
        its live slot — colors are slot-keyed (they follow the row position, matching
        the canvas), so after a removal the remaining rows show their slot's color."""
        if self._swatch is None:
            return
        from .styles_dialog import ColorButton, MATERIAL_PALETTE
        col = len(self._fields)
        for i in range(self.table.rowCount()):
            btn = ColorButton(self._swatch.resolved_hex(i),
                              default_hex=self._swatch.default_hex(i),
                              palette=MATERIAL_PALETTE)
            btn.colorChanged.connect(lambda h, idx=i: self._swatch.set(idx, h))
            self.table.setCellWidget(i, col, btn)

    def _apply_preset_row(self, combo):
        """Fill this row's dependent columns from the preset the driver combo now
        names — the reinforcement sheet's ``Dir``/``Appl`` formulas, in the table.

        Picking a Type SETS Dir and Appl; typing over either afterwards keeps what
        was typed, because nothing rewrites them until the Type is picked again. A
        driver value with no preset (a blank Type — a generic tensile line) fills
        nothing: the sheet's formula leaves those cells empty for it, and the values
        already in them are the loader's defaults for exactly that line.

        Row is resolved from the widget rather than captured, so a row inserted or
        removed above this one cannot make an old index fill the wrong line."""
        if self._preset is None:
            return
        driver_col = self._key_column(self._preset["driver"])
        i = next((r for r in range(self.table.rowCount())
                  if self.table.cellWidget(r, driver_col) is combo), -1)
        if i < 0:
            return
        preset = self._preset["presets"].get(combo.currentText().strip().lower())
        if preset is None:
            return
        for key, value in zip(self._preset["fills"], preset):
            col = self._key_column(key)
            if col is not None:
                self._set_cell_text(i, col, value)

    def _key_column(self, key):
        return next((j for j, f in enumerate(self._fields) if f.key == key), None)

    def _on_combo_changed(self):
        self._apply_dim_all()
        self._emit_change()

    def _row_values(self, i):
        """Current widget values for row ``i`` as a dict — enough for the dim rule to
        read the row's ``option`` (and any other driver) at signal time."""
        vals = {}
        for j, f in enumerate(self._fields):
            w = self.table.cellWidget(i, j)
            if isinstance(w, QComboBox):
                vals[f.key] = w.currentText()
            else:
                it = self.table.item(i, j)
                vals[f.key] = it.text() if it is not None else ""
        return vals

    def _on_item_changed(self, item):
        """A typed cell may change which cells apply — the reinforcement pullout
        law is chosen by whether Adhesion and Delta carry values, not by a combo
        — so the row's graying is re-derived from the edit before the change is
        announced. (The materials rule is combo-driven and re-derives there.)"""
        if (self._dim_rule is not None and self._dim_on_edit
                and not self._suppress_notify):
            prev = self._suppress_notify
            self._suppress_notify = True
            try:
                self._apply_dim_row(item.row())
            finally:
                self._suppress_notify = prev
        self._emit_change()

    def _apply_dim_all(self):
        if self._dim_rule is None:
            return
        prev = self._suppress_notify
        self._suppress_notify = True     # flag edits from setFlags fire nothing
        try:
            for i in range(self.table.rowCount()):
                self._apply_dim_row(i)
        finally:
            self._suppress_notify = prev

    def _apply_dim_row(self, i):
        """Gray (disable, keep value) the fields the dim rule marks inapplicable for
        this row; restore the rest to normal editable/enabled state."""
        dim = self._dim_rule(self._row_values(i))
        for j, f in enumerate(self._fields):
            dimmed = f.key in dim
            w = self.table.cellWidget(i, j)
            if isinstance(w, QComboBox):
                w.setEnabled(not dimmed)
            else:
                it = self.table.item(i, j)
                if it is None:
                    continue
                if dimmed:
                    it.setFlags(Qt.ItemIsSelectable)
                else:
                    it.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled
                                | Qt.ItemIsEditable)

    def _add_row(self):
        i = self.table.rowCount()
        self.table.insertRow(i)
        base = self._new_row()
        self._bases.append(base)
        self._set_row(i, base)
        self._rebuild_swatches()
        self._apply_dim_all()
        self._emit_change()

    def _remove_rows(self):
        for r in sorted({idx.row() for idx in self.table.selectedIndexes()}, reverse=True):
            self.table.removeRow(r)
            if r < len(self._bases):
                self._bases.pop(r)
        self._rebuild_swatches()
        self._apply_dim_all()
        self._emit_change()

    def result_rows(self):
        out = []
        for i in range(self.table.rowCount()):
            base = dict(self._bases[i]) if i < len(self._bases) else self._new_row()
            for j, f in enumerate(self._fields):
                widget = self.table.cellWidget(i, j)
                if isinstance(widget, QComboBox):
                    base[f.key] = f.from_text(widget.currentText())
                else:
                    item = self.table.item(i, j)
                    # read_text, not from_text: an untouched cell hands back the value
                    # the row was built from, so the display rounding never lands in
                    # the record (see Field.read_text).
                    base[f.key] = f.read_text(item.text() if item else "",
                                              base.get(f.key))
            out.append(base)
        return out


def _help_label(text):
    lbl = QLabel(text)
    lbl.setWordWrap(True)
    return lbl


# --------------------------------------------------------------------------- #
# Context-sensitive help strip — a persistent, focus-driven one/two-line pane at
# the bottom of an editor dialog that shows the currently-focused field's help
# text. Built as a reusable widget (_HelpStrip) + a wiring helper (attach_help)
# so any editor can adopt it by supplying a field->text mapping and a resolver
# that names the field for a focused widget; the Materials editor is the pilot.
# --------------------------------------------------------------------------- #
class _HelpStrip(QFrame):
    """A subtle, fixed-height (~2 line) strip of dim help text with a top border.

    Empty by default; ``set_help`` swaps in a field's help string (or clears it
    when focus is nowhere useful). Word-wrapped so a longer string flows onto the
    second line rather than being cut; kept concise so two lines always suffice."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFrameShape(QFrame.NoFrame)
        # Separation from the editing area: a single top hairline. Strong blue
        # text (Norm: make it pop) — distinct from the red/blue header-legend
        # colors, readable on light and dark palettes.
        self.setStyleSheet(
            "_HelpStrip { border-top: 1px solid palette(mid); }"
            "_HelpStrip QLabel { color: #1a56b0; }")
        lay = QHBoxLayout(self)
        lay.setContentsMargins(8, 4, 8, 4)
        lay.setSpacing(0)
        self._label = QLabel("")
        self._label.setWordWrap(True)
        self._label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        f = self._label.font()
        if f.pointSizeF() > 0:
            f.setPointSizeF(f.pointSizeF() * 0.9)   # a touch smaller than body text
        self._label.setFont(f)
        lay.addWidget(self._label, 1)
        # Fix the height at two text lines (+ the layout margins) so the strip is
        # always present and never reflows the dialog as text changes length.
        fm = self._label.fontMetrics()
        self.setFixedHeight(fm.lineSpacing() * 2 + 8 + 2)

    def set_help(self, text):
        self._label.setText(text or "")


def attach_help(dialog, mapping, resolver):
    """Wire a context-sensitive :class:`_HelpStrip` into ``dialog``.

    ``mapping`` is the field-key -> help-text dict (the single source of truth for
    the dialog's help). ``resolver`` is ``focused_widget -> key | None``: it names
    the field a focused widget belongs to (or returns None when focus is on
    something without help). The strip is inserted at the bottom of the dialog's
    layout — just above the button box if the last item is one — and updated on
    every ``QApplication.focusChanged`` whose new focus lands inside the dialog.

    Returns the created strip; also stored on ``dialog._help_strip`` so callers can
    push updates from signals a focus change doesn't cover (e.g. a table's
    ``currentCellChanged`` when arrow keys move the current cell without moving
    keyboard focus off the table). Reusable by any editor without refactoring:
    supply a mapping and a resolver for that editor's widgets."""
    strip = _HelpStrip(dialog)
    lay = dialog.layout()
    idx = lay.count()
    # Sit above a trailing button box (the conventional last row) if present.
    last = lay.itemAt(idx - 1) if idx else None
    if last is not None and isinstance(last.widget(), QDialogButtonBox):
        idx -= 1
    lay.insertWidget(idx, strip)
    dialog._help_strip = strip

    def _update(widget):
        if widget is None:
            return
        # Only react to focus inside this dialog; leave the strip as-is otherwise
        # (e.g. focus moving to another window shouldn't blank the help).
        w = widget
        while w is not None:
            if w is dialog:
                break
            w = w.parentWidget()
        if w is None:
            return
        key = resolver(widget)
        strip.set_help(mapping.get(key, "") if key else "")

    app = QApplication.instance()
    if app is not None:
        slot = lambda _old, now: _update(now)
        app.focusChanged.connect(slot)
        # Drop the global connection when the dialog dies so the slot can't fire
        # against a deleted dialog.
        dialog.destroyed.connect(lambda *_: app.focusChanged.disconnect(slot))
    return strip


def _table_help_resolver(get_table, get_fields):
    """Build an ``attach_help`` resolver for a table-based dialog: derive the
    logical column from the focused widget (via ``_EditableTable.column_at``) and
    map it to that column's field key. ``get_table``/``get_fields`` are callables
    (not the values themselves) so the resolver keeps working when the active
    table is rebuilt or swapped — e.g. a tabbed dialog whose "current" table
    changes with the active tab."""
    def resolver(widget):
        table = get_table()
        if table is None:
            return None
        fields = get_fields()
        col = table.column_at(widget)
        if col is None or not (0 <= col < len(fields)):
            return None
        return fields[col].key
    return resolver


def _wire_cell_help(strip, table, fields, mapping):
    """Push ``table``'s current-cell help into ``strip`` on every
    ``currentCellChanged``. Arrow-key navigation moves the table's current cell
    without moving keyboard focus off the table, so the ``focusChanged``-driven
    strip (``attach_help``) alone would miss it — this covers that gap, mirroring
    the Materials table's ``currentCellChanged`` wiring. No-op if either
    ``strip`` or ``table`` is None (e.g. help wasn't requested for this dialog)."""
    if strip is None or table is None:
        return

    def _push(cr, cc, pr, pc):
        key = fields[cc].key if 0 <= cc < len(fields) else None
        strip.set_help(mapping.get(key, "") if key else "")

    table.table.currentCellChanged.connect(_push)


def _usage_legend(fields):
    """A colored legend for any usage-tagged columns, or None if there are none."""
    present = [u for u in ("lem", "fem", "seep", "rel")
               if any(getattr(f, "usage", None) == u for f in fields)]
    if not present:
        return None
    spans = [f'<b><span style="color:{USAGE_COLOR[u]}">&#9632; {USAGE_NAME[u]}</span></b>'
             for u in present]
    lbl = QLabel("Header colors:&nbsp;&nbsp; " + "&nbsp;&nbsp;&nbsp; ".join(spans))
    lbl.setTextFormat(Qt.RichText)
    return lbl


def _ok_cancel(dialog, layout):
    bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    bb.accepted.connect(dialog.accept)
    bb.rejected.connect(dialog.reject)
    layout.addWidget(bb)


# --------------------------------------------------------------------------- #
# Editor preview drawing — the full cross-section with the object being edited
# HIGHLIGHTED, shared by the profile-line / polygon / starting-circle editors
# (Materials list-view plots generalized). One emphasis style + one dim style,
# consistent across editors: emphasis is the single deliberate color choice; "dim"
# is just reduced alpha on each feature's own (material/style-derived) color, so the
# dimmed context still reads as itself. Each drawer renders from the editor's PENDING
# rows; the debounced PreviewPane keeps the last good frame if a half-typed row makes
# a drawer raise.
# --------------------------------------------------------------------------- #
_PREVIEW_EMPH = "#ff6d00"        # emphasis color for the selected object (saturated
                                 # orange — stands out against the tab10 material
                                 # palette and the red trial surfaces alike)
_PREVIEW_EMPH_LW = 3.0           # emphasized line width
_PREVIEW_DIM_ALPHA = 0.35        # alpha applied to non-selected context geometry


def _doc_style(parent):
    """The project's resolved style delta reached via the parent window's document,
    so a preview renders with the same material colors / feature styles as the Inputs
    plot. None (defaults) in headless/round-trip use where there is no document."""
    doc = getattr(parent, "doc", None)
    return getattr(doc, "style", None) if doc is not None else None


def _finish_preview_axes(ax):
    """Match plot_inputs' final touches so the preview frames the section the same
    way: equal aspect (data-lim adjustable, so circles read round), no grid, a touch
    of headroom. Titles/legends are omitted — the preview is a focused thumbnail."""
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(False)
    y0, y1 = ax.get_ylim()
    if y1 > y0:
        pad = 0.05 * (y1 - y0)
        ax.set_ylim(y0, y1 + pad)


# v21 polygon Type vocabulary for the polygon editor's Type combo: the template's
# own words (the polygon sheet's row-5 validation list) paired with the loader's
# internal kinds. Both halves are imported rather than restated, so the words the
# Studio offers and the kinds it stores can never drift from the file's.
POLYGON_TYPE_ITEMS = list(POLYGON_TYPE_WORDS.items())   # [(word, kind), ...]
# Display wording for a non-material kind, used in the item list and the preview.
# 'refine' has no SSR label — it is a meshing overlay, not an analysis one.
POLYGON_KIND_LABELS = dict(SSR_ZONE_LABELS, refine="Refine")
# kind -> the plot feature whose color the preview borrows, so a zone previews in
# the same hue plot_inputs draws it in.
_ZONE_PREVIEW_FEATURE = {"reduce": "ssr_zone_reduce", "hold": "ssr_zone_hold",
                         "hold_elastic": "ssr_zone_elastic"}
# A refine region has no material and no analysis meaning, so it borrows no feature
# color: a muted neutral keeps it legible without reading as a soil zone.
_REFINE_PREVIEW_COLOR = "#8d6e63"


def _material_color(style, mat_id, fallback_idx, kind="material"):
    """Preview color for one polygon: an overlay's feature color when ``kind`` names
    one, otherwise the material palette entry for ``mat_id``."""
    from xslope.style import feature_style, material_style
    if kind == "refine":
        return _REFINE_PREVIEW_COLOR
    if kind in _ZONE_PREVIEW_FEATURE:
        return feature_style(style, _ZONE_PREVIEW_FEATURE[kind]).get("color", "black")
    idx = mat_id if mat_id is not None else fallback_idx
    return material_style(style, idx)["color"]


def _draw_profile_preview(ax, lines, selected, max_depth, slope_data, style):
    """Preview for the profile-line editor: ONLY the PENDING profile lines plus the
    max-depth base line — nothing else. The selected line is bold (emphasis color)
    with vertex markers; the others keep their material color, thin and dimmed. The
    max-depth base stays because this same dialog edits max_depth. The other input
    overlays (piezo / loads / reinforcement / piles) are deliberately NOT drawn here:
    Norm asked this preview stay a clean profile-lines-only view — they cluttered the
    picture and aren't what's being edited. Zone fills are NOT re-derived either
    (rebuilding them from a half-edited line is fragile), so this is a line preview
    (see the caption)."""
    from xslope.plot import plot_max_depth
    from xslope.style import resolve_style
    rstyle = resolve_style(style)
    if max_depth is not None:
        drawable = [{"coords": ln["coords"]} for ln in lines if ln.get("coords")]
        if drawable:
            try:
                plot_max_depth(ax, drawable, max_depth, style=rstyle)
            except Exception:
                pass
    for i, ln in enumerate(lines):
        coords = ln.get("coords") or []
        if not coords:
            continue
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        color = _material_color(rstyle, ln.get("mat_id"), i)
        if i == selected:
            ax.plot(xs, ys, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW,
                    marker="o", markersize=5, markerfacecolor=_PREVIEW_EMPH,
                    markeredgecolor="white", zorder=20)
        else:
            ax.plot(xs, ys, color=color, linewidth=1.2, alpha=_PREVIEW_DIM_ALPHA,
                    zorder=6)
    _finish_preview_axes(ax)


def _draw_polygon_preview(ax, polys, selected, slope_data, style):
    """Preview for the polygon editor: the PENDING material zones. The selected zone
    is filled (its material color, higher opacity + hatch) with a bold emphasis edge
    and vertex markers; the others are dimmed fills. The domain base gives context
    (uncluttered, per the preview rule). Rings are closed for display."""
    from xslope.plot import plot_domain_base
    from xslope.style import resolve_style
    rstyle = resolve_style(style)
    try:
        plot_domain_base(ax, slope_data.get("domain_polygon"), style=rstyle)
    except Exception:
        pass
    for i, pg in enumerate(polys):
        coords = pg.get("coords") or []
        if len(coords) < 2:
            continue
        xs = [c[0] for c in coords] + [coords[0][0]]   # close the ring for display
        ys = [c[1] for c in coords] + [coords[0][1]]
        color = _material_color(rstyle, pg.get("mat_id"), i,
                                pg.get("kind") or "material")
        if i == selected:
            ax.fill(xs, ys, facecolor=color, alpha=0.55, hatch="//",
                    edgecolor=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW, zorder=20)
            ax.plot(xs, ys, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW,
                    marker="o", markersize=4, markerfacecolor=_PREVIEW_EMPH,
                    markeredgecolor="white", zorder=21)
        else:
            ax.fill(xs, ys, facecolor=color, alpha=0.12, zorder=4)
            ax.plot(xs, ys, color=color, linewidth=1.0, alpha=_PREVIEW_DIM_ALPHA,
                    zorder=5)
    _finish_preview_axes(ax)


def _circle_radius_depth(c):
    """(Xo, Yo, R, Depth) from a pending circle row, resolving R/Depth from the row's
    Option EXACTLY as CirclesEditor.apply does — so the preview tracks the option live
    (a Radius/Intercept row previews from R, a Depth row from Depth)."""
    def f(key):
        try:
            return float(c.get(key, 0) or 0)
        except (TypeError, ValueError):
            return 0.0
    Xo, Yo = f("Xo"), f("Yo")
    opt = str(c.get("Option", "Depth"))
    if opt == "Radius":
        R = f("R"); depth = Yo - R
    elif opt == "Intercept":
        R = ((f("Xi") - Xo) ** 2 + (f("Yi") - Yo) ** 2) ** 0.5
        depth = Yo - R
    else:
        depth = f("Depth"); R = Yo - depth
    return Xo, Yo, R, depth


def _circle_arc(slope_data, Xo, Yo, R, depth):
    """The clipped failure-surface arc for a circle, as a list of (x, y), or None —
    exactly what the engine draws for it. Shared by the circles drawer and its pick
    resolver so both hit-test the same geometry. Guarded: returns None on any error."""
    ground_surface = slope_data.get("ground_surface")
    if ground_surface is None or getattr(ground_surface, "is_empty", False):
        return None
    try:
        from shapely.geometry import LineString
        from xslope.slice import generate_failure_surface
        ok, res = generate_failure_surface(
            ground_surface, circular=True,
            circle={"Xo": Xo, "Yo": Yo, "Depth": depth, "R": R},
            tcrack_depth=slope_data.get("tcrack_depth", 0))
        if not ok:
            return None
        clipped = res[4]
        if not isinstance(clipped, LineString):
            clipped = LineString(clipped)
        return list(clipped.coords)
    except Exception:
        return None


#: Color for the search-window overlay in the circles preview — a limit on WHERE a
#: search may run is neither a trial surface (red) nor the selected circle (orange),
#: so it gets its own hue rather than borrowing one of theirs.
_WINDOW_COLOR = "#00838f"


def _ground_band(slope_data, x_lo, x_hi, n=48):
    """The ground surface between ``x_lo`` and ``x_hi`` as (xs, ys), clipped to the
    section's own x range, or None when the range falls outside it entirely."""
    import numpy as np
    gs = slope_data.get("ground_surface")
    if gs is None or getattr(gs, "is_empty", False):
        return None
    coords = list(gs.coords)
    gx = [p[0] for p in coords]
    gy = [p[1] for p in coords]
    if gx[0] > gx[-1]:
        gx, gy = gx[::-1], gy[::-1]
    lo = max(min(x_lo, x_hi), min(gx))
    hi = min(max(x_lo, x_hi), max(gx))
    if not hi > lo:
        return None
    xs = np.linspace(lo, hi, n)
    return xs, np.interp(xs, gx, gy)


def _draw_search_window(ax, slope_data, window, annot_line):
    """The search window over the section: the entry and exit ranges as bars lying on
    the ground surface, the center box as a dashed rectangle.

    Only a limit the search will actually apply is drawn — a range with one end filled
    is not a window and the engine ignores it, so drawing it would show a constraint
    that is not there. The center box is an annotation artist (like the circle centers
    it confines) because it usually sits well above the section and must not inflate
    the framed view."""
    for lo, hi, label in (("entry_x_min", "entry_x_max", "entry"),
                          ("exit_x_min", "exit_x_max", "exit")):
        if window.get(lo) is None or window.get(hi) is None:
            continue
        band = _ground_band(slope_data, window[lo], window[hi])
        if band is None:
            continue
        xs, ys = band
        ax.plot(xs, ys, color=_WINDOW_COLOR, linewidth=4.0, alpha=0.75,
                solid_capstyle="butt", zorder=8)
        text = ax.annotate(label, (xs[len(xs) // 2], ys[len(ys) // 2]),
                           textcoords="offset points", xytext=(0, 6), ha="center",
                           fontsize=7, color=_WINDOW_COLOR, zorder=9)
        text.set_in_layout(False)
    box = ("center_box_x_min", "center_box_x_max",
           "center_box_y_min", "center_box_y_max")
    if all(window.get(k) is not None for k in box):
        x0, x1, y0, y1 = (window[k] for k in box)
        annot_line([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0],
                   color=_WINDOW_COLOR, linewidth=1.2, linestyle="--")


def _draw_circles_preview(ax, circles, selected, slope_data, style, window=None):
    """Preview for the starting-circles editor: the full cross-section (base geometry
    + overlays) with the PENDING circles over it. Each circle's clipped failure arc is
    drawn as the engine does; the selected one is bold (emphasis color) with a center
    marker, a radius line and a depth line, the others faint. The center/radius are
    annotation-layer artists (in_layout=False, clipped) so a center far above the
    section can't inflate the framed view — matching plot_circles.

    ``window`` is the pending search window (the editor's group, live as it is typed),
    drawn under the circles: the entry/exit ranges as bars on the ground surface, the
    center box dashed."""
    from matplotlib.lines import Line2D
    from xslope.plot import plot_base_geometry
    from xslope.style import resolve_style
    rstyle = resolve_style(style)
    plot_base_geometry(ax, slope_data, labels=False, style=rstyle)

    def _annot_line(xs, ys, **kw):
        ln = Line2D(xs, ys, **kw)
        ln.set_in_layout(False)          # keep an off-section center out of the
        ln.set_clip_box(ax.bbox)         # tight-bbox layout + view autoscale
        ln.set_clip_on(True)
        ax.add_artist(ln)

    if window:
        _draw_search_window(ax, slope_data, window, _annot_line)

    for i, c in enumerate(circles):
        Xo, Yo, R, depth = _circle_radius_depth(c)
        if R <= 0:
            continue
        emph = (i == selected)
        arc = _circle_arc(slope_data, Xo, Yo, R, depth)
        if arc:
            ax_ = [p[0] for p in arc]; ay_ = [p[1] for p in arc]
            if emph:
                ax.plot(ax_, ay_, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW,
                        zorder=20)
            else:
                ax.plot(ax_, ay_, color="red", linestyle="--", linewidth=1.0,
                        alpha=_PREVIEW_DIM_ALPHA, zorder=6)
        if emph:
            # Center +, radius line (down to the circle bottom) and depth line, all
            # as annotation artists so they never blow up the framed section.
            _annot_line([Xo], [Yo], marker="+", color=_PREVIEW_EMPH, linestyle="None",
                        markersize=12, markeredgewidth=2)
            _annot_line([Xo, Xo], [Yo, depth], color=_PREVIEW_EMPH, linewidth=1.2)
            _annot_line([Xo - R, Xo + R], [depth, depth], color=_PREVIEW_EMPH,
                        linewidth=1.0, linestyle=":")
    _finish_preview_axes(ax)


# --------------------------------------------------------------------------- #
# Overlay-feature previews (reinforcement / piles / line loads / non-circular /
# distributed loads / piezometric lines). These edit features that layer OVER the
# base geometry, so each preview draws the FULL cross-section as backdrop — the base
# geometry plus every other overlay — and then draws the feature being edited itself,
# the selected object emphasized and the rest dimmed. Same emphasis/dim vocabulary as
# the geometry previews above.
# --------------------------------------------------------------------------- #
def _draw_section_context(ax, slope_data, style, include=()):
    """Backdrop for a feature preview — Norm's rule: keep it uncluttered, just
    enough to confirm placement. Draws the base geometry (profile/polygon layers
    + max-depth/domain base) and ONLY the overlay features named in ``include``.
    The one standing inclusion: the dloads editor passes ('piezo',), since a
    distributed load often represents a water load and reads against the water
    line. Guarded so a half-edited model can never blank the backdrop."""
    from xslope import plot as _p
    try:
        _p.plot_base_geometry(ax, slope_data, labels=False, style=style)
    except Exception:
        pass
    extras = {
        "piezo": lambda: _p.plot_piezo_line(ax, slope_data, style=style),
        "dloads": lambda: _p.plot_dloads(ax, slope_data, style=style),
    }
    for name in include:
        fn = extras.get(name)
        if fn is not None:
            try:
                fn()
            except Exception:
                pass


def _xy(row, kx, ky):
    """(x, y) floats from a pending row, or None if either is blank/half-typed."""
    try:
        return float(row.get(kx, 0) or 0), float(row.get(ky, 0) or 0)
    except (TypeError, ValueError):
        return None


def _draw_reinforcement_preview(ax, rows, selected, slope_data, style):
    """Preview for the reinforcement editor: the base geometry ONLY (profile/polygon
    layers + domain base — no loads, surfaces, or other overlays; Norm: they get in
    the way here) with each pending line from (x1,y1)→(x2,y2). The selected line is
    bold (emphasis color) with endpoint markers; the others keep the reinforcement
    color, thin and dimmed. Pullout-length positions — lp1 measured from end 1 and
    lp2 from end 2, where the available tension reaches t_max — are dotted on every
    line, mirroring the standard reinforcement plot's tension breakpoints: they are
    part of the geometry being edited."""
    import math as _math
    from xslope.plot import plot_base_geometry
    from xslope.style import resolve_style, feature_style
    rstyle = resolve_style(style)
    try:
        plot_base_geometry(ax, slope_data, labels=False, style=rstyle)
    except Exception:
        pass
    base = feature_style(rstyle, "reinforcement").get("color", "darkgray")
    for i, r in enumerate(rows):
        p1, p2 = _xy(r, "x1", "y1"), _xy(r, "x2", "y2")
        if p1 is None or p2 is None:
            continue
        xs, ys = [p1[0], p2[0]], [p1[1], p2[1]]
        emph = i == selected
        if emph:
            ax.plot(xs, ys, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW,
                    marker="o", markersize=6, markerfacecolor=_PREVIEW_EMPH,
                    markeredgecolor="white", zorder=20)
        else:
            ax.plot(xs, ys, color=base, linewidth=2.0, alpha=_PREVIEW_DIM_ALPHA,
                    zorder=6)
        # Pullout-length breakpoints along the line (only when 0 < lp < length).
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        length = _math.hypot(dx, dy)
        if length <= 0:
            continue
        ux, uy = dx / length, dy / length
        for lp_key, (ox, oy), sgn in (("lp1", p1, +1), ("lp2", p2, -1)):
            try:
                lp = float(r.get(lp_key, 0) or 0)
            except (TypeError, ValueError):
                continue
            if 0 < lp < length:
                px, py = ox + sgn * ux * lp, oy + sgn * uy * lp
                if emph:
                    ax.plot([px], [py], marker="o", markersize=7,
                            markerfacecolor="white", markeredgecolor=_PREVIEW_EMPH,
                            markeredgewidth=1.8, zorder=21)
                else:
                    ax.plot([px], [py], marker="o", markersize=4, color=base,
                            alpha=_PREVIEW_DIM_ALPHA, zorder=7)
    _finish_preview_axes(ax)


def _draw_piles_preview(ax, rows, selected, slope_data, style):
    """Preview for the piles editor: each pending pile from (x1,y1)→(x2,y2) over the
    full section. The selected pile is bold (emphasis color) with a square cap marker
    at its upper end and a downward-triangle tip marker at its lower end; the others
    keep the pile color, thinner and dimmed."""
    from xslope.style import resolve_style, feature_style
    rstyle = resolve_style(style)
    _draw_section_context(ax, slope_data, rstyle)
    base = feature_style(rstyle, "piles").get("color", "green")
    for i, p in enumerate(rows):
        p1, p2 = _xy(p, "x1", "y1"), _xy(p, "x2", "y2")
        if p1 is None or p2 is None:
            continue
        xs, ys = [p1[0], p2[0]], [p1[1], p2[1]]
        if i == selected:
            ax.plot(xs, ys, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW + 1.0,
                    solid_capstyle="butt", zorder=20)
            cap, tip = (p1, p2) if p1[1] >= p2[1] else (p2, p1)  # cap = upper end
            ax.plot([cap[0]], [cap[1]], marker="s", markersize=9, color=_PREVIEW_EMPH,
                    markeredgecolor="white", markeredgewidth=1.2, zorder=21)
            ax.plot([tip[0]], [tip[1]], marker="v", markersize=10, color=_PREVIEW_EMPH,
                    markeredgecolor="white", markeredgewidth=1.2, zorder=21)
        else:
            ax.plot(xs, ys, color=base, linewidth=3.0, alpha=_PREVIEW_DIM_ALPHA,
                    solid_capstyle="butt", zorder=6)
    _finish_preview_axes(ax)


def _draw_line_loads_preview(ax, rows, selected, slope_data, style):
    """Preview for the line-loads editor: each pending load as an arrow pointing IN
    the force direction, its head on the point of application (mirrors plot_line_loads
    exactly, but with the selected load emphasized and the others dimmed). Arrow-tail
    length is 6% of the model span, so the glyph and its label stay in view."""
    import numpy as np
    from xslope.style import resolve_style
    rstyle = resolve_style(style)
    # Trial surfaces don't inform load placement (and a far circle center can inflate
    # the frame), so the section + other overlays are context enough here.
    _draw_section_context(ax, slope_data, rstyle)
    gs = slope_data.get("ground_surface")
    if gs is not None and not getattr(gs, "is_empty", False):
        xs = [p[0] for p in gs.coords]
        span = max(xs) - min(xs)
    else:
        span = 100.0
    alen = 0.06 * span
    for i, ll in enumerate(rows):
        pt = _xy(ll, "x", "y")
        if pt is None:
            continue
        try:
            ang = np.radians(float(ll.get("angle", -90.0) or 0))
            P = float(ll.get("P", 0) or 0)
        except (TypeError, ValueError):
            continue
        x, y = pt
        tx, ty = x - np.cos(ang) * alen, y - np.sin(ang) * alen
        emph = (i == selected)
        color = _PREVIEW_EMPH if emph else "purple"
        lw = _PREVIEW_EMPH_LW if emph else 2.0
        alpha = 1.0 if emph else _PREVIEW_DIM_ALPHA
        ax.annotate("", xy=(x, y), xytext=(tx, ty), annotation_clip=False,
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, alpha=alpha))
        ax.plot([tx], [ty], linestyle="None")   # keep the arrow (+ label) framed
        label = str(ll.get("label") or "L")
        ax.annotate(f"{label}={P:.0f}", (tx, ty), textcoords="offset points",
                    xytext=(4, 4), fontsize=8, color=color, alpha=alpha,
                    fontweight="bold" if emph else "normal")
    _finish_preview_axes(ax)


def _draw_seep_bc_preview(ax, entries, selected, slope_data, style, set_no=1):
    """Preview for the seep BC editor: the base geometry with every boundary in the
    set drawn as its polyline — the selected entry bold (emphasis color) with point
    markers, the others dimmed in the same hue families the main canvas uses (set-1
    head navy / reservoir azure / flux green / exit face red; set 2 rose / seagreen /
    orangered). ``entries`` mirrors the widget's list order (heads … fluxes … exit
    face) as dicts with 'kind' and 'coords'."""
    from xslope.plot import _SEEP_BC_RESERVOIR, _SEEP_BC_SET2_HEAD
    from xslope.style import resolve_style, feature_style
    rstyle = resolve_style(style)
    _draw_section_context(ax, slope_data, rstyle)
    if int(set_no) == 1:
        colors = {"head": feature_style(rstyle, "seep_bc").get("color", "darkblue"),
                  "reservoir": _SEEP_BC_RESERVOIR[0],
                  "flux": feature_style(rstyle, "seep_flux").get("color", "darkgreen"),
                  "exit": feature_style(rstyle, "seep_exit_face").get("color", "red")}
    else:
        colors = {"head": _SEEP_BC_SET2_HEAD[0], "reservoir": _SEEP_BC_SET2_HEAD[0],
                  "flux": "seagreen", "exit": "orangered"}
    styles = {"head": "--", "reservoir": "--", "flux": "-.", "exit": "--"}
    for i, e in enumerate(entries):
        coords = [c for c in (e.get("coords") or []) if c is not None]
        if not coords:
            continue
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        kind = str(e.get("kind") or "head").strip().lower()
        if i == selected:
            ax.plot(xs, ys, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW,
                    marker="o", markersize=5, markerfacecolor=_PREVIEW_EMPH,
                    markeredgecolor="white", zorder=20)
        else:
            ax.plot(xs, ys, color=colors.get(kind, "darkblue"),
                    linestyle=styles.get(kind, "--"), linewidth=2.0,
                    alpha=_PREVIEW_DIM_ALPHA, zorder=6)
    _finish_preview_axes(ax)


_NCPT_MARKER = {"Fixed": "s", "Horiz": "D", "Free": "o"}   # per movement type


def _draw_noncirc_preview(ax, rows, selected, slope_data, style):
    """Preview for the non-circular editor: the base geometry ONLY (profile/polygon
    layers + max-depth base — no piezo or other overlays; Norm's declutter) with the
    pending polyline (all rows) drawn bold, each vertex marked with a per-movement
    glyph — □ Fixed, ◇ Horiz-only, ○ Free — and the selected vertex enlarged and
    filled. Horiz vertices also carry a small ↔ direction glyph. Points are ordered
    left→right as entered."""
    import numpy as np
    from xslope.plot import plot_base_geometry
    from xslope.style import resolve_style
    rstyle = resolve_style(style)
    try:
        plot_base_geometry(ax, slope_data, labels=False, style=rstyle)
    except Exception:
        pass
    pts = []
    for r in rows:
        pt = _xy(r, "X", "Y")
        if pt is not None:
            pts.append((pt[0], pt[1], str(r.get("Movement", "Free"))))
    if len(pts) >= 2:
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=_PREVIEW_EMPH,
                linewidth=_PREVIEW_EMPH_LW, zorder=18)
    xs = [p[0] for p in pts]
    glyph = ((max(xs) - min(xs)) * 0.035) if len(xs) >= 2 and max(xs) > min(xs) else 1.0
    for i, (x, y, mv) in enumerate(pts):
        emph = (i == selected)
        # Horiz vertices carry a ↔ direction glyph UNDER the marker, so an open
        # (non-selected) marker shows it through and the caption legend explains it.
        if mv == "Horiz":
            ax.annotate("", xy=(x + glyph, y), xytext=(x - glyph, y), zorder=19,
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle="<->", color=_PREVIEW_EMPH, lw=1.5))
        ax.plot([x], [y], marker=_NCPT_MARKER.get(mv, "o"),
                markersize=11 if emph else 7, color=_PREVIEW_EMPH,
                markerfacecolor=_PREVIEW_EMPH if emph else "white",
                markeredgecolor=_PREVIEW_EMPH, markeredgewidth=1.6, zorder=20)
    _finish_preview_axes(ax)


def _draw_dload_block(ax, block, color, gamma_w, emph, vertical=False):
    """Draw one distributed-load block: its surface polyline, the load-profile line
    (surface offset by Normal/γ_w, i.e. the equivalent water depth) and pressure
    arrows sampled along it. Emphasized blocks use the emphasis color, full opacity
    and a heavier weight; the rest dim to context. Raises nothing the caller doesn't
    guard — a half-typed point skips the whole block.

    ``vertical`` is the v21 Direction: a vertical block is DRAWN vertical — arrows
    straight up off the loaded surface rather than leaning with it — matching
    plot_dloads, so the preview shows the direction the user just picked."""
    import numpy as np
    pts = []
    for p in block:
        try:
            pts.append((float(p["X"]), float(p["Y"]), float(p["Normal"])))
        except (TypeError, ValueError, KeyError):
            return
    if len(pts) < 2 or gamma_w <= 0:
        return
    ec = _PREVIEW_EMPH if emph else color
    lw = _PREVIEW_EMPH_LW if emph else 1.6
    alpha = 0.95 if emph else _PREVIEW_DIM_ALPHA
    z = 20 if emph else 6
    xs = [p[0] for p in pts]
    ax.plot(xs, [p[1] for p in pts], color=ec, linewidth=lw, alpha=alpha, zorder=z)
    total_dx = max(xs) - min(xs)
    step = (total_dx / 14) if total_dx > 0 else 1.0   # ~14 arrows across the block
    top_xs, top_ys = [], []
    for i in range(len(pts) - 1):
        x1, y1, n1 = pts[i]
        x2, y2, n2 = pts[i + 1]
        dx, dy = x2 - x1, y2 - y1
        seg = np.hypot(dx, dy)
        if seg == 0:
            continue
        if vertical:
            perp_x, perp_y = 0.0, 1.0          # gravity surcharge: straight up
        else:
            perp_x, perp_y = -dy / seg, dx / seg  # segment dir +90° (away from soil)
        n = max(1, int(round(abs(dx) / step))) if step else 1
        for t in np.linspace(0, 1, n + 1):
            x, y = x1 + t * dx, y1 + t * dy
            h = (n1 + t * (n2 - n1)) / gamma_w
            ax_, ay_ = x + perp_x * h, y + perp_y * h
            top_xs.append(ax_)
            top_ys.append(ay_)
            if h > 1e-6:
                ax.annotate("", xy=(x, y), xytext=(ax_, ay_), annotation_clip=False,
                            arrowprops=dict(arrowstyle="-|>", color=ec, alpha=alpha,
                                            lw=lw * 0.6, shrinkA=0, shrinkB=0))
    if top_xs:
        ax.plot(top_xs, top_ys, color=ec, linewidth=lw, alpha=alpha, zorder=z)


def _draw_dloads_preview(ax, set1_blocks, set2_blocks, active_set, selected_block,
                         slope_data, style, set_dirs=None):
    """Preview for the distributed-loads editor: both sets on the section (set 1 and
    set 2 in their distinct feature colors, so a rapid-drawdown pair reads apart), the
    selected block of the ACTIVE tab emphasized and everything else dimmed. Follows the
    active set tab and the selected load within it."""
    from xslope.style import resolve_style, feature_style
    from xslope.units import require_gamma_water
    rstyle = resolve_style(style)
    _draw_section_context(ax, slope_data, rstyle, include=("piezo",))
    # A Studio project always carries gamma_water (seeded on new/loaded projects), so
    # read it loudly here rather than masking a missing value with a hardcoded default
    # that could preview arrows at the wrong scale for the model's unit system.
    gamma_w = require_gamma_water(slope_data, "distributed-load preview")
    c1 = feature_style(rstyle, "dloads").get("color", "purple")
    c2 = feature_style(rstyle, "dloads2").get("color", "orange")
    dirs = set_dirs or ([], [])
    for si, (blocks, color) in enumerate(((set1_blocks, c1), (set2_blocks, c2))):
        sd = dirs[si] if si < len(dirs) else []
        for bi, block in enumerate(blocks or []):
            emph = (si == active_set and bi == selected_block)
            vert = str(sd[bi] if bi < len(sd) else "normal").lower() == "vertical"
            try:
                _draw_dload_block(ax, block, color, gamma_w, emph, vertical=vert)
            except Exception:
                pass
    _finish_preview_axes(ax)


def _draw_piezo_preview(ax, rows_per_tab, active_tab, slope_data, style):
    """Preview for the piezometric-lines editor: BOTH piezo lines on the section, the
    one whose tab is active bold (emphasis color) with vertex markers, the other in
    its own color, thin and dimmed. Line 2 is the rapid-drawdown / second table."""
    from xslope.style import resolve_style, feature_style
    rstyle = resolve_style(style)
    _draw_section_context(ax, slope_data, rstyle)
    colors = [feature_style(rstyle, "piezo_line").get("color", "b"),
              feature_style(rstyle, "piezo_line2").get("color", "skyblue")]
    for ti, rows in enumerate(rows_per_tab):
        pts = [p for p in (_xy(r, "x", "y") for r in rows) if p is not None]
        if not pts:
            continue
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        if ti == active_tab:
            ax.plot(xs, ys, color=_PREVIEW_EMPH, linewidth=_PREVIEW_EMPH_LW,
                    marker="o", markersize=5, markerfacecolor=_PREVIEW_EMPH,
                    markeredgecolor="white", zorder=20)
        else:
            ax.plot(xs, ys, color=colors[ti] if ti < len(colors) else "gray",
                    linewidth=1.6, marker="o", markersize=3,
                    alpha=_PREVIEW_DIM_ALPHA, zorder=6)
    _finish_preview_axes(ax)


# --------------------------------------------------------------------------- #
# Preview-pane PICK resolvers — the inverse of the drawers above. Each maps a
# clicked (x, y) + a data-unit tolerance (both from the canvas, tol already scaled
# from ~8 screen px through the current zoom) to a selection target, hit-testing
# the SAME pending geometry its drawer renders. They return None for a click beyond
# tolerance (empty space), so the wiring leaves the selection untouched. Point-to-
# segment/polyline distance is the shared `_line_dist`; a two-point list is a
# single segment, so vertex/segment/polyline picks all go through it.
# --------------------------------------------------------------------------- #
from shapely.geometry import Point as _Point, Polygon as _Polygon


def _pick_circles(rows, x, y, tol, slope_data):
    """Nearest starting-circle row to the click — by its clipped failure arc (or the
    full circle if the arc can't be built) or its center marker. Row index or None."""
    pt = _Point(x, y)
    best_i, best_d = None, float("inf")
    for i, c in enumerate(rows):
        Xo, Yo, R, depth = _circle_radius_depth(c)
        if R <= 0:
            continue
        d = pt.distance(_Point(Xo, Yo))                         # center marker
        arc = _circle_arc(slope_data, Xo, Yo, R, depth)
        if arc:
            d = min(d, _line_dist(pt, arc))                     # the drawn arc
        else:
            d = min(d, abs(((x - Xo) ** 2 + (y - Yo) ** 2) ** 0.5 - R))
        if d < best_d:
            best_i, best_d = i, d
    return best_i if best_d <= tol else None


def _pick_line_rows(rows, x, y, tol):
    """Nearest 2-endpoint line row (reinforcement / piles): the segment or either
    endpoint (all select the same row). Row index within tol, else None."""
    pt = _Point(x, y)
    best_i, best_d = None, float("inf")
    for i, r in enumerate(rows):
        p1, p2 = _xy(r, "x1", "y1"), _xy(r, "x2", "y2")
        if p1 is None or p2 is None:
            continue
        d = _line_dist(pt, [p1, p2])
        if d < best_d:
            best_i, best_d = i, d
    return best_i if best_d <= tol else None


def _pick_line_loads(rows, x, y, tol, slope_data):
    """Nearest line-load row — its arrow (tail→application point). Mirrors the arrow
    geometry the drawer builds (tail length = 6% of the model span). Row or None."""
    import numpy as np
    pt = _Point(x, y)
    gs = slope_data.get("ground_surface")
    if gs is not None and not getattr(gs, "is_empty", False):
        xs = [p[0] for p in gs.coords]
        span = max(xs) - min(xs)
    else:
        span = 100.0
    alen = 0.06 * span
    best_i, best_d = None, float("inf")
    for i, ll in enumerate(rows):
        p = _xy(ll, "x", "y")
        if p is None:
            continue
        try:
            ang = np.radians(float(ll.get("angle", -90.0) or 0))
        except (TypeError, ValueError):
            continue
        px, py = p
        tx, ty = px - np.cos(ang) * alen, py - np.sin(ang) * alen
        d = _line_dist(pt, [(tx, ty), (px, py)])
        if d < best_d:
            best_i, best_d = i, d
    return best_i if best_d <= tol else None


def _pick_noncirc(rows, x, y, tol):
    """Non-circular surface: nearest vertex row, or (failing that) the nearest
    segment's start-vertex row. Row index within tol, else None."""
    pt = _Point(x, y)
    pts, idx = [], []
    for i, r in enumerate(rows):
        p = _xy(r, "X", "Y")
        if p is not None:
            pts.append(p); idx.append(i)
    vbest_i, vbest_d = None, float("inf")
    for k, p in enumerate(pts):
        d = pt.distance(_Point(p))
        if d < vbest_d:
            vbest_i, vbest_d = idx[k], d
    if vbest_d <= tol:
        return vbest_i
    sbest_i, sbest_d = None, float("inf")
    for k in range(len(pts) - 1):
        d = _line_dist(pt, [pts[k], pts[k + 1]])
        if d < sbest_d:
            sbest_i, sbest_d = idx[k], d
    return sbest_i if sbest_d <= tol else None


def _pick_matgeom_lines(lines, x, y, tol, closed):
    """Profile lines / polygons. Returns ``(feature, row|None)``: a vertex within tol
    → (that line, that vertex row); else a segment/edge within tol → (that line,
    None); else, for polygons (``closed``), a click inside a zone → (that zone,
    None). None when nothing is within tolerance. Vertices win over edges/interior,
    mirroring the emphasis the drawer gives them."""
    pt = _Point(x, y)
    vfeat, vrow, vd = None, None, float("inf")
    efeat, ed = None, float("inf")
    for fi, ln in enumerate(lines):
        coords = ln.get("coords") or []
        if not coords:
            continue
        for ri, c in enumerate(coords):
            d = pt.distance(_Point(c[0], c[1]))
            if d < vd:
                vfeat, vrow, vd = fi, ri, d
        ring = list(coords) + ([coords[0]] if closed and len(coords) >= 2 else [])
        if len(ring) >= 2:
            d = _line_dist(pt, ring)
            if d < ed:
                efeat, ed = fi, d
    if vfeat is not None and vd <= tol:
        return (vfeat, vrow)
    if efeat is not None and ed <= tol:
        return (efeat, None)
    if closed:                              # interior: the smallest containing zone
        cfeat, carea = None, float("inf")
        for fi, ln in enumerate(lines):
            coords = ln.get("coords") or []
            if len(coords) < 3:
                continue
            try:
                poly = _Polygon(coords)
                if poly.contains(pt) and poly.area < carea:
                    cfeat, carea = fi, poly.area
            except Exception:
                pass
        if cfeat is not None:
            return (cfeat, None)
    return None


def _pick_piezo(rows_per_tab, x, y, tol):
    """Piezometric lines. Returns ``(tab, row|None)``: the tab of the nearest line and
    the nearest vertex row on it (or None for a bare segment hit). None beyond tol.
    A hit on the OTHER tab's line switches the active tab (its line 'becomes active')."""
    pt = _Point(x, y)
    vtab, vrow, vd = None, None, float("inf")
    stab, srow, sd = None, None, float("inf")
    for ti, rows in enumerate(rows_per_tab):
        pts, ridx = [], []
        for i, r in enumerate(rows):
            p = _xy(r, "x", "y")
            if p is not None:
                pts.append(p); ridx.append(i)
        for k, p in enumerate(pts):
            d = pt.distance(_Point(p))
            if d < vd:
                vtab, vrow, vd = ti, ridx[k], d
        for k in range(len(pts) - 1):
            d = _line_dist(pt, [pts[k], pts[k + 1]])
            if d < sd:
                stab, srow, sd = ti, ridx[k], d
    if vtab is not None and vd <= tol:
        return (vtab, vrow)
    if stab is not None and sd <= tol:
        return (stab, srow)
    return None


def _pick_dloads(set_blocks, x, y, tol):
    """Distributed loads. ``set_blocks = (blocks_set1, blocks_set2)``. Returns
    ``(set, block)`` for the nearest block outline (its surface polyline) within tol,
    else None. A hit in the other set switches the active tab to it."""
    pt = _Point(x, y)
    best = (None, None, float("inf"))
    for si, blocks in enumerate(set_blocks):
        for bi, block in enumerate(blocks or []):
            coords = []
            ok = True
            for p in block:
                try:
                    coords.append((float(p["X"]), float(p["Y"])))
                except (TypeError, ValueError, KeyError):
                    ok = False
                    break
            if not ok or len(coords) < 2:
                continue
            d = _line_dist(pt, coords)
            if d < best[2]:
                best = (si, bi, d)
    return (best[0], best[1]) if best[2] <= tol else None


# --------------------------------------------------------------------------- #
# Fitting a table to its own content
#
# A table dialog that opens with a column past its right edge is a table dialog with
# a hidden input, and one that gives every column the width of its widest is a dialog
# most of whose width is blank: an x-coordinate does not need the room a
# "geosynthetic" combo does, and seventeen columns sized for that combo is a
# seventeen-hundred pixel dialog. So each column is measured on ITS OWN widest thing
# — header text, cell text, and the cell WIDGETS a combo column is made of, whose
# width no item hint knows about — and given that, plus a cushion.
#
# A numeric column is floored at the widest ordinary number instead of at its header,
# so a table opened EMPTY is as wide as the same table opened full and typing the
# first row never moves a column; a text column floors on its header, which is the
# widest thing it is guaranteed to hold.
#
# Everything is measured in the widget's own font, so the fit follows the font, the
# display and the platform instead of a pixel guess made on one of them.
# --------------------------------------------------------------------------- #

#: Rows of empty table to leave visible below the last one, so a short table still
#: looks like something you can add to.
_TABLE_SPARE_ROWS = 3
#: Preview-pane depth, in lines of the dialog's own text. Deep enough to read a
#: cross-section in; it scales with the font rather than fixing a pixel height.
_PREVIEW_MIN_LINES = 16
#: A negative number written to the display precision — the widest ordinary thing a
#: numeric cell holds. It floors a NUMERIC column's width, so a table opened EMPTY is
#: as wide as the same table opened full, and typing the first row does not need a
#: resize.
_NUMBER_SAMPLE = "-" + "0" * _DISPLAY_SIG_DIGITS + "."
#: A free-text column's floor — room for a name a little longer than "Line 10".
_TEXT_SAMPLE = "0" * 12


def _column_widths(table, numeric):
    """The width each column needs to show its own widest content, header included.

    ``numeric`` is the set of column indexes holding numbers; those are the ones
    floored at :data:`_NUMBER_SAMPLE`."""
    from PySide6.QtWidgets import QComboBox, QStyle

    header = table.horizontalHeader()
    fm = table.fontMetrics()
    cushion = fm.horizontalAdvance("00")                    # one em-ish, both sides
    number = fm.horizontalAdvance(_NUMBER_SAMPLE) + cushion
    # A free-text column (a label, a name) floors on a sample wider than its
    # header so a short entry does not leave it cramped.
    text_floor = fm.horizontalAdvance(_TEXT_SAMPLE) + cushion
    widths = []
    for c in range(table.columnCount()):
        w = max(table.sizeHintForColumn(c), header.sectionSizeHint(c))
        free_text = combo = False
        for r in range(table.rowCount()):
            cell = table.cellWidget(r, c)
            if isinstance(cell, QComboBox):
                # Just wide enough for the longest choice it offers plus the
                # arrow: a combo's own sizeHint pads well beyond that.
                longest = max((fm.horizontalAdvance(cell.itemText(i))
                               for i in range(cell.count())), default=0)
                arrow = cell.style().pixelMetric(QStyle.PM_ScrollBarExtent)
                w = max(w, longest + arrow + cushion)
                combo = True
            elif cell is not None:
                w = max(w, cell.sizeHint().width())
            else:
                # Measured row by row rather than left to sizeHintForColumn, which
                # only looks at the rows currently in the viewport: the longest
                # label in a forty-row table is usually not one of the first
                # twenty, and a fit that had not seen it would cut it off.
                item = table.item(r, c)
                if item is not None and item.text():
                    w = max(w, fm.horizontalAdvance(item.text()))
                    if c not in numeric:
                        free_text = True
        if not combo:                # a combo column's cushion is already in
            w += cushion
        if c in numeric:
            w = max(w, number)
        elif free_text:
            w = max(w, text_floor)
        widths.append(max(w, header.minimumSectionSize()))
    return widths


def _fit_columns(table, numeric):
    """Size each column to its own content (:func:`_column_widths`). Returns the
    viewport width the VISIBLE columns need.

    The sections stay Interactive at the fitted width rather than stretching to fill
    the viewport: a fitted column is already as wide as the widest thing in it, so
    sharing out spare width only pads it, and a dialog dragged narrower would take
    that width back out of every column at once. Left alone, the fit survives a
    resize in both directions — the spare width is blank space on the right, and a
    dialog dragged narrower than its columns scrolls, which is the one direction a
    table has a scroll bar for. The user keeps the column drag Interactive means."""
    header = table.horizontalHeader()
    # Forget the previous fit before measuring: a minimum section size left standing
    # is reported back as the sections' own hint, so re-fitting would ratchet the
    # columns wider every time it ran.
    header.setMinimumSectionSize(-1)
    header.setSectionResizeMode(QHeaderView.Interactive)
    header.setStretchLastSection(False)
    widths = _column_widths(table, numeric)
    for c, w in enumerate(widths):
        header.resizeSection(c, w)
    return sum(w for c, w in enumerate(widths) if not table.isColumnHidden(c))


def _dialog_width_for(dialog, editable):
    """The width ``dialog`` needs to show every visible column of the table pane
    ``editable``: the pane's fitted width plus the dialog's own margins, and never
    narrower than the rest of the dialog already needs — the view toggle and usage
    checkboxes above the table, the buttons below it. A four-column table does not
    get to fold the toggle bar."""
    margins = dialog.layout().contentsMargins()
    return max(editable.columns_width() + margins.left() + margins.right(),
               dialog.layout().minimumSize().width())


def _set_dialog_width(dialog, want, cap=None, grow_only=False):
    """Take ``dialog`` to ``want`` pixels wide, never past ``cap`` or the screen.

    ``grow_only`` for a dialog already on screen: showing a hidden column (a usage
    toggle ticked back on) or switching to a wider view must not leave a column past
    the right edge, but neither may it take back width the user gave the dialog
    themselves. It is applied after the caps, so a dialog someone has already made
    wider than the screen (a headless capture does exactly that) is left alone."""
    if cap is not None:
        want = min(want, cap)
    screen = dialog.screen() or QApplication.primaryScreen()
    if screen is not None:
        want = min(want, screen.availableGeometry().width())
    if grow_only:
        want = max(want, dialog.width())
    if want != dialog.width():
        dialog.resize(want, dialog.height())


def _fit_dialog_to_columns(dialog, editable, cap=None, grow_only=False):
    """Set ``dialog``'s width to what the fitted columns of ``editable`` need."""
    _set_dialog_width(dialog, _dialog_width_for(dialog, editable),
                      cap=cap, grow_only=grow_only)


def _grow_dialog_to(dialog, want):
    """Widen ``dialog`` to ``want`` if it is narrower than that; never shrink it."""
    _set_dialog_width(dialog, want, grow_only=True)


class TableEditorDialog(QDialog):
    """Editable table over a list of dict records.

    ``usage_toggles`` (a list of analysis tags, e.g. ``["lem", "fem"]``) adds a row
    of checkboxes that show/hide the columns specific to each analysis, so the user
    sees only the inputs relevant to what they're doing. The toggle state persists
    per dialog. When omitted, a static color legend is shown instead.

    ``preview_draw`` (a hook ``draw(ax, rows, selected_index)``) attaches a live
    preview pane on the right behind a splitter — the highlight follows the selected
    table row. This keeps a pure-table editor a table (no list-view rewrite); only the
    editors that pass a hook (currently starting circles) grow a preview.

    ``preview_below`` stacks that preview UNDER the table instead of beside it, and
    sizes the dialog from the table's measured content (see :func:`_fit_columns`). A
    wide table and a preview cannot both have the width they need side by side, and
    the table is the half that loses columns off its right edge when they compete —
    so a table with more columns than a pane's width can hold takes the full width
    and the section, which is wide and short anyway, takes the space below it.

    ``extra_widget`` is a widget of settings that belong to the table as a whole
    rather than to any row of it (the circles editor's search window), placed under
    the table and above the buttons. It is reached back as ``dlg.extra``. Two optional
    hooks on it are honoured: ``set_on_change(cb)`` wires it to the live preview, and
    ``validate() -> message`` refuses OK with that message rather than letting the
    dialog save something the loader would reject.

    ``generate`` (a hook ``propose() -> (rows, message, reason)``) adds a button that
    derives the whole table from the model. It follows the same contract as a
    preflight remedy: a button that cannot run is DIMMED with the reason in its
    tooltip rather than live and failing when pressed, and what it will do is stated
    before it does it. The proposed rows land in the table, not in the document, so
    the preview shows them and Cancel still discards them."""

    def __init__(self, title, fields, rows, new_row, parent=None, help_text=None,
                 usage_toggles=None, preview_draw=None, preview_caption=None,
                 pick_resolve=None, field_help=None, unit_labels=None,
                 generate=None, preview_below=False, extra_widget=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._title = title
        self._fields = fields
        self._unit_labels = unit_labels
        self._preview_draw = preview_draw
        self._preview = None
        # ``pick_resolve(x, y, tol, rows) -> row_index | None`` maps a preview click
        # to the table row to select (None = beyond tolerance, leave selection).
        self._pick_resolve = pick_resolve
        # Optional field-key -> help-text mapping for the context-sensitive help
        # strip (see attach_help). None (the default) leaves the dialog without one.
        self._field_help = field_help
        self.resize(min(1200, 160 + 110 * len(fields)), 460)
        layout = QVBoxLayout(self)
        if help_text:
            layout.addWidget(_help_label(help_text))
        on_change = self._schedule_preview if preview_draw is not None else None
        on_select = self._schedule_preview if preview_draw is not None else None
        self._editable = _EditableTable(fields, rows, new_row,
                                        on_change=on_change, on_select=on_select,
                                        unit_labels=self._unit_labels)
        if usage_toggles:
            layout.addLayout(self._build_toggle_bar(usage_toggles))
        else:
            legend = _usage_legend(fields)
            if legend:
                layout.addWidget(legend)
        if preview_draw is not None:
            from .canvas import PreviewPane
            self._preview = PreviewPane(
                lambda ax: self._preview_draw(ax, self._editable.result_rows(),
                                              self._editable.selected_row()),
                caption=preview_caption)
            if self._pick_resolve is not None:
                self._preview.clicked.connect(self._on_preview_click)
            split = QSplitter(Qt.Vertical if preview_below else Qt.Horizontal)
            split.addWidget(self._editable)
            split.addWidget(self._preview)
            self._split = split
            layout.addWidget(split, 1)
            if preview_below:
                # The table keeps the height its rows ask for and the preview absorbs
                # everything else, so dragging the dialog taller grows the picture.
                split.setStretchFactor(0, 0)
                split.setStretchFactor(1, 1)
            else:
                split.setStretchFactor(0, 1)
                split.setStretchFactor(1, 1)
                split.setSizes([560, 500])
                self.resize(min(1280, 620 + 110 * len(fields)), 520)
        else:
            layout.addWidget(self._editable)
        if generate is not None:
            layout.addLayout(self._build_generate_bar(generate))
        self.extra = extra_widget
        if extra_widget is not None:
            layout.addWidget(extra_widget)
            if self._preview is not None and hasattr(extra_widget, "set_on_change"):
                extra_widget.set_on_change(self._schedule_preview)
        _ok_cancel(self, layout)
        if usage_toggles:
            self._apply_toggles()      # set initial column visibility
        if self._preview is not None:
            self._preview.refresh_now()
        if self._field_help is not None:
            attach_help(self, self._field_help,
                       _table_help_resolver(lambda: self._editable, lambda: self._fields))
            _wire_cell_help(self._help_strip, self._editable, self._fields, self._field_help)
        # Sized from its own table where the table IS the dialog's width — the
        # preview stacked below it. A preview BESIDE the table splits the width with
        # it, so there the table's fit is not the dialog's width and the opening
        # size stays the split one.
        self._content_sized = bool(preview_below)
        if preview_below:
            # Last, so the measurement sees every strip the dialog ended up with.
            self._size_to_content()

    def accept(self):
        """OK, unless the extra widget says what it holds cannot be saved — a window
        the loader would refuse to read back is stopped here, with the pair named,
        rather than at the next open of the file."""
        extra = getattr(self, "extra", None)
        problem = extra.validate() if hasattr(extra, "validate") else ""
        if problem:
            QMessageBox.warning(self, self._title, problem)
            return
        super().accept()

    def _table_pane_height(self, rows_shown=None):
        """The height the table half needs for ``rows_shown`` rows (its own row count
        by default), plus a few spare rows, its header and the Add/Remove bar."""
        table = self._editable.table
        vh = table.verticalHeader()
        if rows_shown is None:
            rows = sum(table.rowHeight(r) for r in range(table.rowCount()))
        else:
            rows = rows_shown * vh.defaultSectionSize()
        spare = _TABLE_SPARE_ROWS * vh.defaultSectionSize()
        chrome = 2 * table.frameWidth()
        # The Add/Remove bar is whatever the pane's own hint holds beyond the table's.
        bar = max(0, self._editable.sizeHint().height() - table.sizeHint().height())
        return table.horizontalHeader().height() + rows + spare + chrome + bar

    def _content_width(self):
        """Dialog width that fits every column: the table pane's own fitted width
        (:meth:`_EditableTable.columns_width`) plus the dialog's margins, and never
        narrower than the rows above and below the table — the toggle bar, the
        buttons — already need."""
        return _dialog_width_for(self, self._editable)

    def _size_to_content(self):
        """Open at the size the content asks for, capped by the screen.

        Width comes from the table (every column fitted, plus the scroll bar and
        frame it lives inside); height from the table's own rows plus a preview deep
        enough to read a section in and whatever strips the dialog carries. Nothing
        here is a remembered pixel count -- change the font, the columns or the row
        count and the dialog opens to match.

        Each pane also gets a MINIMUM, which is what keeps it honest afterwards: the
        dialog cannot be dragged small enough to squeeze the preview into a strip, and
        anything that appears later (the generator's summary line) grows the layout's
        own minimum, so the dialog grows with it rather than taking the room out of
        the picture. Where the two cannot both be satisfied -- a table of forty
        circles is taller than any screen -- the preview keeps its minimum and the
        table scrolls, which is the pane that has a scroll bar for the purpose."""
        width = self._content_width()
        pane = self._table_pane_height()
        # The caption wraps at the width it is GIVEN, which is the dialog's width less
        # the layout margins -- reserving at the full width buys one line too few and
        # cuts the last line off along the bottom edge.
        margins = self.layout().contentsMargins()
        caption_width = width - margins.left() - margins.right()
        preview = (_PREVIEW_MIN_LINES * self.fontMetrics().height()
                   + self._preview.reserve_caption(caption_width))
        self._preview.setMinimumHeight(preview)
        self._editable.setMinimumHeight(self._table_pane_height(rows_shown=0))
        # A QSplitter does not carry its children's minimums into its own, so a strip
        # that appears later (the generator's summary) is free to squeeze the split
        # rather than grow the dialog -- and what gets cut is the bottom of the
        # preview's caption. Stating the minimum here puts the panes into the layout's
        # own minimum, so the dialog grows for the new strip instead.
        self._split.setMinimumHeight(self._editable.minimumHeight() + preview
                                     + self._split.handleWidth())
        # Everything that is not the splitter: help text, legend, generate bar,
        # buttons, help strip. sizeHint() knows them all; the splitter's own hint is
        # replaced by the two pane heights measured above.
        chrome = self.sizeHint().height() - self._split.sizeHint().height()
        height = chrome + pane + preview

        screen = self.screen() or QApplication.primaryScreen()
        if screen is not None:
            avail = screen.availableGeometry()
            width = min(width, avail.width())
            height = min(height, avail.height())
        self.resize(width, height)
        split_height = max(0, height - chrome)
        pane = max(0, min(pane, split_height - preview))
        self._split.setSizes([pane, split_height - pane])

    def _refit_columns(self):
        """Re-fit the columns to content that arrived after the dialog opened (a
        generated set of rows), widening the dialog if the new content needs it.

        Only ever wider: a dialog that shrank itself while the user was working in it
        would be taking room away from them to save room they did not ask to save.

        Only for a dialog that was sized from its content in the first place. The
        fit resets the columns to the new content's widths, which discards any width
        the user dragged a column to; in a dialog whose whole width was measured to
        fit its columns that is the fit doing its job, while in one sized some other
        way it would only be overwriting the user's own arrangement."""
        if not getattr(self, "_content_sized", False):
            return
        self._grow_to_layout_minimum()
        needed = self._content_width()
        if needed > self.width():
            self._grow(needed)
        self._grow_to_fit_columns()

    def _grow_to_layout_minimum(self):
        """Grow tall enough for everything the layout now holds, measured at the
        width the dialog actually has.

        A strip that appears after the dialog was sized — the generator's summary
        line — needs room that nothing gave it, and what pays for it is the bottom of
        the preview's caption. ``minimumSizeHint`` is no help: it measures every
        wrapped label at the dialog's own MINIMUM width, where the same text takes
        more lines than it does at the width on screen, so it answers a question
        nobody asked. Summing the rows at the real width is the measurement that
        matches what is drawn."""
        layout = self.layout()
        margins = layout.contentsMargins()
        width = self.width() - margins.left() - margins.right()
        want = margins.top() + margins.bottom()
        want += layout.spacing() * max(0, layout.count() - 1)
        for i in range(layout.count()):
            item = layout.itemAt(i)
            if item.hasHeightForWidth():
                want += item.heightForWidth(width)
                continue
            widget = item.widget()
            want += (max(widget.sizeHint().height(), widget.minimumHeight())
                     if widget is not None else item.sizeHint().height())
        if want > self.height():
            self.resize(self.width(), want)

    def _grow(self, width):
        """Widen to ``width``, never past the screen."""
        screen = self.screen() or QApplication.primaryScreen()
        avail = screen.availableGeometry().width() if screen is not None else width
        self.resize(min(width, avail), self.height())

    def showEvent(self, event):
        super().showEvent(event)
        if getattr(self, "_content_sized", False):
            self._grow_to_layout_minimum()
            self._grow_to_fit_columns()

    def _grow_to_fit_columns(self):
        """Close any gap between what the columns need and what the viewport gives
        them, now that the widgets have been laid out.

        The estimate made before the first layout is close but not exact: the
        row-number gutter widens with the row count, and a style adds its own header
        margins only when it lays the header out. Both are small and both are in the
        direction that hides a column, so the last word on the fit belongs to the
        widgets once they have real geometry."""
        table = self._editable.table
        deficit = table.horizontalHeader().length() - table.viewport().width()
        if deficit > 0:
            self._grow(self.width() + deficit)

    def _build_generate_bar(self, spec):
        """The generator button, dimmed with its own reason when it cannot run.

        Same contract as a preflight remedy: a button that cannot succeed is dimmed
        with the reason in its tooltip rather than live and failing when pressed.
        ``spec`` carries ``label``, ``available``, ``reason``, ``tooltip``,
        ``propose(parent) -> (rows, message, reason)``, the noun for what it makes
        (``unit``, e.g. "circle"), and ``append`` — whether adding the generated rows
        to the ones already there is a sensible thing to offer. It is for a table of
        independent records (trial circles) and not for one that IS a single object
        (the points of one non-circular surface).

        The bar carries the generator's summary line under the button, because the
        summary is the audit: it names which face each circle came from, how it was
        sized, and what it threw away — the things a reader checks the generated rows
        against. It has to stay readable while they do that, so it is written into
        the dialog rather than into a box they have to dismiss first."""
        self._generate = spec
        bar = QVBoxLayout()
        row = QHBoxLayout()
        btn = QPushButton(spec.get("label") or "Generate…")
        btn.setObjectName("generate_button")
        available = bool(spec.get("available", True))
        btn.setEnabled(available)
        btn.setToolTip(spec.get("reason") if not available
                       else (spec.get("tooltip") or ""))
        btn.clicked.connect(self._run_generate)
        self.generate_button = btn
        row.addWidget(btn)
        row.addStretch(1)
        bar.addLayout(row)
        self.generate_summary = QLabel()
        self.generate_summary.setObjectName("generate_summary")
        self.generate_summary.setWordWrap(True)
        self.generate_summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.generate_summary.setVisible(False)
        bar.addWidget(self.generate_summary)
        return bar

    def _generate_unit(self, n):
        """"3 circles" / "1 point" — the generator says what it makes."""
        noun = self._generate.get("unit") or "point"
        return f"{n} {noun}" if n == 1 else f"{n} {noun}s"

    def _ask_generate_disposition(self, existing, proposed, message):
        """What to do with rows already in the table: ``"replace"``, ``"append"`` or
        ``None`` for cancel.

        Asked only when there ARE rows, and never answered for the user: work that is
        already in the table is the user's, and a generator that quietly threw it away
        would be a generator nobody could risk pressing."""
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Question)
        box.setWindowTitle("Replace what is already here?")
        box.setText(f"This table already holds {self._generate_unit(existing)}. "
                    f"The generator built {self._generate_unit(proposed)}.")
        box.setInformativeText(message)
        buttons = QMessageBox.Apply | QMessageBox.Cancel
        if self._generate.get("append"):
            buttons |= QMessageBox.Yes
        box.setStandardButtons(buttons)
        box.button(QMessageBox.Apply).setText("Replace")
        if self._generate.get("append"):
            box.button(QMessageBox.Yes).setText("Append")
        box.setDefaultButton(QMessageBox.Apply)
        answer = box.exec()
        if answer == QMessageBox.Apply:
            return "replace"
        if answer == QMessageBox.Yes and self._generate.get("append"):
            return "append"
        return None

    def _run_generate(self):
        """Ask the generator for rows, ask the user what to do with the ones already
        in the table, then do it and show what was done.

        The generated rows land in the TABLE, not in the document: the preview
        redraws on them, the user can still edit them, and Cancel discards them the
        way it discards any other edit made in this dialog."""
        rows, message, reason = self._generate["propose"](self)
        if not rows:
            if reason:
                QMessageBox.information(self, "Nothing was generated", reason)
            return
        existing = self._editable.result_rows()
        how = "replace"
        if existing:                       # nothing to lose = nothing to ask about
            how = self._ask_generate_disposition(len(existing), len(rows), message)
            if how is None:
                return
        self._editable.replace_rows(existing + rows if how == "append" else rows)
        self._refit_columns()          # generated numbers may be wider than typed ones
        self._show_generate_summary(message, len(rows), how if existing else "fill")
        self._schedule_preview()

    def _show_generate_summary(self, message, n, how):
        """Write the generator's own summary into the dialog, under the button.

        The summary is a strip that was not there a moment ago, and a dialog whose
        height does not change has to find its room somewhere: the somewhere is the
        preview, which is the only thing in the dialog that stretches. So the dialog
        grows by exactly what the layout says it now needs more of -- measured as the
        difference in the layout's own total hint across the label appearing, which
        counts the label's wrapped height and the spacing around it and needs no
        pixel figure of its own. The preview keeps the area it had."""
        lead = {"append": "Added", "replace": "Replaced the table with"}.get(how,
                                                                            "Generated")
        label = self.generate_summary
        layout = self.layout()
        before = layout.totalSizeHint().height()
        label.setText(f"{lead} {self._generate_unit(n)}: {message}".rstrip(". ") + ".")
        label.setVisible(True)
        # A wrapped label's height is only knowable once its width is: reserve it, or
        # the layout budgets one line for a summary that takes two and the difference
        # is taken out of the preview below it.
        margins = layout.contentsMargins()
        wrapped = label.heightForWidth(self.width() - margins.left() - margins.right())
        if wrapped > 0:
            label.setMinimumHeight(wrapped)
        layout.invalidate()
        grown = layout.totalSizeHint().height() - before
        if getattr(self, "_content_sized", False):
            self._size_to_content()        # re-fit around the strip that just appeared
        elif grown > 0:
            screen = self.screen() or QApplication.primaryScreen()
            avail = (screen.availableGeometry().height() if screen is not None
                     else self.height() + grown)
            self.resize(self.width(), min(self.height() + grown, avail))

    def _schedule_preview(self, *_):
        if self._preview is not None:
            self._preview.schedule()

    def _on_preview_click(self, x, y, tol):
        """A click in the preview: resolve it to a table row and select that row
        (the selection-changed path then re-renders the preview with the new
        emphasis). Beyond tolerance the resolver returns None → no-op."""
        row = self._pick_resolve(x, y, tol, self._editable.result_rows())
        if row is not None:
            self._editable.select_row(row)

    def _build_toggle_bar(self, tags):
        from PySide6.QtWidgets import QCheckBox
        from PySide6.QtCore import QSettings
        s = QSettings("XSlope", "XSlope Studio")
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Show parameters for:"))
        self._toggles = {}
        for t in tags:
            cb = QCheckBox(USAGE_TOGGLE_LABEL[t])
            default = (t != "rel")     # reliability σ columns hidden by default (niche)
            cb.setChecked(bool(s.value(f"editor_toggles/{self._title}/{t}",
                                       default, type=bool)))
            cb.setStyleSheet(f"color:{USAGE_COLOR[t]}; font-weight:bold;")
            cb.toggled.connect(self._apply_toggles)
            self._toggles[t] = cb
            bar.addWidget(cb)
        bar.addStretch(1)
        return bar

    def _apply_toggles(self):
        from PySide6.QtCore import QSettings
        enabled = {t for t, cb in self._toggles.items() if cb.isChecked()}
        self._editable.apply_usage_filter(enabled)
        s = QSettings("XSlope", "XSlope Studio")
        for t, cb in self._toggles.items():
            s.setValue(f"editor_toggles/{self._title}/{t}", cb.isChecked())

    def result_rows(self):
        return self._editable.result_rows()


class TabbedTableEditorDialog(QDialog):
    """A tabbed dialog, one editable table per tab (e.g. the two piezo lines).

    ``preview_draw`` (a hook ``draw(ax, rows_per_tab, active_tab)`` where
    ``rows_per_tab`` is each tab's live rows) attaches a live preview pane on the
    right behind a splitter; the highlight follows the ACTIVE tab. Each tab's table
    drives the preview via its ``on_change`` hook, and a tab switch reschedules it."""

    def __init__(self, title, tabs, parent=None, help_text=None,
                 preview_draw=None, preview_caption=None, pick_resolve=None,
                 field_help=None, unit_labels=None):
        # tabs: list of (tab_title, fields, rows, new_row)
        super().__init__(parent)
        self.setWindowTitle(title)
        self._unit_labels = unit_labels
        self._preview_draw = preview_draw
        self._preview = None
        # ``pick_resolve(x, y, tol, rows_per_tab) -> (tab, row|None) | None`` maps a
        # preview click to the tab + points-table row to select.
        self._pick_resolve = pick_resolve
        # Optional field-key -> help-text mapping for the context-sensitive help
        # strip; each tab keeps its OWN fields list (usually identical across tabs,
        # e.g. piezo's x/y), so the resolver looks up the ACTIVE tab's fields.
        self._field_help = field_help
        self._tab_fields = [fields for _, fields, _, _ in tabs]
        max_cols = max(len(fields) for _, fields, _, _ in tabs)
        self.resize(min(1200, 200 + 110 * max_cols), 480)
        layout = QVBoxLayout(self)
        if help_text:
            layout.addWidget(_help_label(help_text))
        self._tabs = QTabWidget()
        self._editables = []
        on_change = self._schedule_preview if preview_draw is not None else None
        for tab_title, fields, rows, new_row in tabs:
            et = _EditableTable(fields, rows, new_row, on_change=on_change,
                                unit_labels=self._unit_labels)
            self._tabs.addTab(et, tab_title)
            self._editables.append(et)
        if preview_draw is not None:
            from .canvas import PreviewPane
            self._tabs.currentChanged.connect(self._schedule_preview)
            self._preview = PreviewPane(
                lambda ax: self._preview_draw(
                    ax, [et.result_rows() for et in self._editables],
                    self._tabs.currentIndex()),
                caption=preview_caption)
            if self._pick_resolve is not None:
                self._preview.clicked.connect(self._on_preview_click)
            split = QSplitter(Qt.Horizontal)
            split.addWidget(self._tabs)
            split.addWidget(self._preview)
            split.setStretchFactor(0, 1)
            split.setStretchFactor(1, 1)
            split.setSizes([560, 500])
            layout.addWidget(split, 1)
            self.resize(min(1280, 620 + 110 * max_cols), 520)
        else:
            layout.addWidget(self._tabs)
        _ok_cancel(self, layout)
        if self._preview is not None:
            self._preview.refresh_now()
        if self._field_help is not None:
            attach_help(self, self._field_help,
                       _table_help_resolver(self._active_table, self._active_fields))
            for et, flds in zip(self._editables, self._tab_fields):
                _wire_cell_help(self._help_strip, et, flds, self._field_help)
            self._tabs.currentChanged.connect(self._push_active_col_help)

    def _active_table(self):
        i = self._tabs.currentIndex()
        return self._editables[i] if 0 <= i < len(self._editables) else None

    def _active_fields(self):
        i = self._tabs.currentIndex()
        return self._tab_fields[i] if 0 <= i < len(self._tab_fields) else []

    def _push_active_col_help(self, *_):
        """On a tab switch, refresh the strip for the newly-active tab's current
        cell — a tab change doesn't itself move keyboard focus, so focusChanged
        wouldn't otherwise fire."""
        strip = getattr(self, "_help_strip", None)
        if strip is None:
            return
        table = self._active_table()
        fields = self._active_fields()
        col = table.table.currentColumn() if table is not None else -1
        key = fields[col].key if 0 <= col < len(fields) else None
        strip.set_help(self._field_help.get(key, "") if key else "")

    def _schedule_preview(self, *_):
        if self._preview is not None:
            self._preview.schedule()

    def _on_preview_click(self, x, y, tol):
        """A click in the preview: resolve to (tab, row), switch to that tab (its
        line 'becomes active') and select the row. The tab-change / selection paths
        re-render the preview. Beyond tolerance the resolver returns None → no-op."""
        hit = self._pick_resolve(x, y, tol,
                                 [et.result_rows() for et in self._editables])
        if hit is None:
            return
        tab, row = hit
        if 0 <= tab < self._tabs.count():
            self._tabs.setCurrentIndex(tab)
            if row is not None:
                self._editables[tab].select_row(row)

    def result_rows(self, index):
        return self._editables[index].result_rows()


class _BlockListWidget(QWidget):
    """Master/detail over a list of blocks, where each block is a list of dict
    rows: a list of blocks (left) + the selected block's row table (right).
    Used for distributed loads (and reusable for other block-structured inputs).

    A block may also carry ONE per-block property (``prop_spec``) edited by a combo
    above the row table — the distributed loads' v21 Direction. The property is held
    ON the block record, not in a parallel list keyed by position: a parallel list is
    what silently shifts every direction onto the wrong load the moment a block is
    deleted or reordered, and the record makes that structurally impossible."""

    def __init__(self, blocks, fields, new_row, block_label="Load", parent=None,
                 on_change=None, unit_labels=None, prop_spec=None, values=None):
        # prop_spec: {"label", "items": [(display, value), ...], "help_key"} for the
        #   per-block combo, or None for no per-block property.
        # values: the per-block property values, positional, as they arrive from the
        #   model (which stores them as a parallel list). They are zipped onto the
        #   block records HERE, once, and travel with them from then on.
        super().__init__(parent)
        self._fields = fields
        self._new_row = new_row
        self._unit_labels = unit_labels
        self._block_label = block_label
        # Optional live-edit hook (used by the dloads preview): fires on any point
        # edit (via the inner table), block add/remove, or block selection change.
        self._on_change = on_change
        self._prop_spec = prop_spec
        self._prop_default = (prop_spec["items"][0][1] if prop_spec else None)
        _vals = list(values or [])
        self._blocks = [{"rows": [dict(r) for r in blk],
                         "prop": (_vals[i] if i < len(_vals) else self._prop_default)}
                        for i, blk in enumerate(blocks or [])]
        self._cur = -1
        self.table = None
        # Context-sensitive help: set by the owning dialog AFTER construction (it
        # needs the dialog's layout/button box built first — see attach_help), then
        # applied retroactively via _wire_help() to whichever table exists then and
        # to every one built afterward (_load rebuilds the table per block).
        self._help_strip = None
        self._field_help = None

        body = QHBoxLayout(self)
        body.setContentsMargins(0, 0, 0, 0)
        left = QVBoxLayout()
        body.addLayout(left)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._on_select)
        left.addWidget(self.list)
        lb = QHBoxLayout()
        b_add = QPushButton(f"Add {block_label.lower()}")
        b_add.clicked.connect(self._add_block)
        b_rem = QPushButton("Remove")
        b_rem.clicked.connect(self._remove_block)
        lb.addWidget(b_add)
        lb.addWidget(b_rem)
        left.addLayout(lb)
        right = QVBoxLayout()
        body.addLayout(right, 1)

        self.prop_combo = None
        if prop_spec is not None:
            prow = QHBoxLayout()
            prow.addWidget(QLabel(prop_spec["label"]))
            self.prop_combo = QComboBox()
            for display, _val in prop_spec["items"]:
                self.prop_combo.addItem(display)
            self.prop_combo.currentIndexChanged.connect(self._on_prop_changed)
            prow.addWidget(self.prop_combo)
            prow.addStretch(1)
            right.addLayout(prow)

        self._holder = QVBoxLayout()
        right.addLayout(self._holder, 1)

        self._refresh_list()
        if self._blocks:
            self.list.setCurrentRow(0)

    # --- per-block property -------------------------------------------------
    def _prop_value_at(self, index):
        items = self._prop_spec["items"] if self._prop_spec else []
        return items[index][1] if 0 <= index < len(items) else self._prop_default

    def _prop_index(self, value):
        items = self._prop_spec["items"] if self._prop_spec else []
        for i, (_d, v) in enumerate(items):
            if v == value:
                return i
        return 0

    def _label(self, i):
        base = f"{self._block_label} {i + 1}"
        if self._prop_spec is None:
            return base
        # The property is shown in the LIST, not just in the detail pane, so a load's
        # direction is readable while deleting or re-picking blocks — the operation
        # that used to shift directions onto the wrong load. The list shows the
        # file's own word; the combo carries the longer explanatory wording.
        return f"{base}  ({self._blocks[i]['prop']})"

    def _on_prop_changed(self, idx):
        if 0 <= self._cur < len(self._blocks):
            self._blocks[self._cur]["prop"] = self._prop_value_at(idx)
            item = self.list.item(self._cur)
            if item:
                item.setText(self._label(self._cur))
        self._notify()

    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i in range(len(self._blocks)):
            self.list.addItem(self._label(i))
        self.list.blockSignals(False)

    def _commit_current(self):
        if 0 <= self._cur < len(self._blocks) and self.table is not None:
            self._blocks[self._cur]["rows"] = self.table.result_rows()
            if self.prop_combo is not None:
                self._blocks[self._cur]["prop"] = self._prop_value_at(
                    self.prop_combo.currentIndex())

    def _load(self, idx):
        if self.table is not None:
            self.table.setParent(None)
            self.table = None
        if not (0 <= idx < len(self._blocks)):
            return
        self.table = _EditableTable(self._fields, self._blocks[idx]["rows"],
                                    self._new_row, on_change=self._notify,
                                    unit_labels=self._unit_labels)
        self._holder.addWidget(self.table)
        self._wire_help()
        if self.prop_combo is not None:
            self.prop_combo.blockSignals(True)
            self.prop_combo.setCurrentIndex(self._prop_index(self._blocks[idx]["prop"]))
            self.prop_combo.blockSignals(False)

    def help_key_for_widget(self, widget):
        """The field key ``widget`` edits — the per-block property combo's key, or
        the current block table's column key, or None."""
        if self.prop_combo is not None and widget is self.prop_combo:
            return self._prop_spec.get("help_key")
        if self.table is None:
            return None
        col = self.table.column_at(widget)
        return self._fields[col].key if 0 <= col < len(self._fields) else None

    def _wire_help(self):
        _wire_cell_help(self._help_strip, self.table, self._fields, self._field_help or {})

    def _notify(self):
        if self._on_change is not None:
            self._on_change()

    def _on_select(self, idx):
        self._commit_current()
        self._cur = idx
        self._load(idx)
        self._notify()

    def _add_block(self):
        self._commit_current()
        self._blocks.append({"rows": [], "prop": self._prop_default})
        self._refresh_list()
        self.list.setCurrentRow(len(self._blocks) - 1)
        self._notify()

    def _remove_block(self):
        idx = self.list.currentRow()
        if idx < 0:
            return
        # One pop removes the block AND its property together — the whole reason the
        # property lives on the record.
        self._blocks.pop(idx)
        self._cur = -1
        self._refresh_list()
        if self._blocks:
            self.list.setCurrentRow(min(idx, len(self._blocks) - 1))
        else:
            self._load(-1)
        self._notify()

    def selected_block(self):
        """Index of the currently selected block, or -1 if none."""
        return self.list.currentRow()

    def pending_blocks(self):
        """A snapshot of every block's ROWS with the in-progress table's LIVE rows
        folded in, without mutating internal state — the preview's view of the
        current edit."""
        out = [list(b["rows"]) for b in self._blocks]
        if 0 <= self._cur < len(out) and self.table is not None:
            out[self._cur] = self.table.result_rows()
        return out

    def pending_values(self):
        """The per-block property values with the in-progress combo folded in — the
        preview's counterpart to :meth:`pending_blocks`, same length, same order, so
        the two can be zipped by the drawing code."""
        if self._prop_spec is None:
            return []
        out = [b["prop"] for b in self._blocks]
        if 0 <= self._cur < len(out) and self.prop_combo is not None:
            out[self._cur] = self._prop_value_at(self.prop_combo.currentIndex())
        return out

    def result_blocks(self):
        self._commit_current()
        return [list(b["rows"]) for b in self._blocks]

    def result_values(self):
        """The per-block property values, positional — the parallel list the model
        stores. Built from the same records ``result_blocks`` returns and in the same
        order, so the two can never disagree about which value belongs to which
        block. Empty when the widget has no per-block property."""
        self._commit_current()
        if self._prop_spec is None:
            return []
        return [b["prop"] for b in self._blocks]


# --------------------------------------------------------------------------- #
# Per-category editors
# --------------------------------------------------------------------------- #
class CategoryEditor:
    label = ""

    def build(self, slope_data, parent):
        raise NotImplementedError

    def apply(self, slope_data, dlg):
        raise NotImplementedError


# main-sheet Units selector (v18): combo label <-> canonical slope_data value.
# Blank label = undeclared (None) -> gamma_w is inferred from magnitude on load and
# time-bearing quantities stay unlabeled. Item order matches the combo, so a combo
# index maps straight into this list.
_UNIT_SYSTEM_ITEMS = [("", None), ("SI", "si"), ("Imperial", "imperial")]

# Time selector legend tokens, verbatim from the v18 template list $D$54:$D$57
# ("sec", not "s"). Blank = time unit undeclared.
_TIME_UNIT_ITEMS = ["", "sec", "min", "hr", "day", "yr"]
# Tension SRF (v19, main!D17) is a TRI-STATE: blank means "unspecified", which is
# not the same as NO. Blank leaves the engine default (reduce the cap) in force and
# keeps the saved file silent on the question; NO pins the static cap.
_TENSION_SRF_ITEMS = [("", None), ("YES", True), ("NO", False)]

WATER_LOADS_HELP = (
    "Who supplies the weight of water standing on the slope. auto derives the "
    "ponded-water load at solve time from the model's own water definition (the "
    "seepage boundary conditions, otherwise the piezometric line); manual leaves "
    "it to the dloads sheets. Switching a model that already draws its reservoir "
    "as load blocks to auto counts the water twice -- the Model checks warn when "
    "that happens.")

K0_HELP = (
    "At-rest lateral earth pressure coefficient K0 for the FEM initial stress "
    "state. Blank (the default) starts from zero stress and switches gravity on in "
    "one step, so the initial lateral stress is whatever elasticity gives — "
    "sigma_h = nu/(1-nu)·sigma_v, about 0.43 at nu = 0.3, set by the STIFFNESS "
    "rather than by the soil. Enter a value and the initial stress is built from "
    "the overburden instead: sigma_h = K0·sigma_v in-plane and out-of-plane, then "
    "iterated to equilibrium. Compacted fills and overconsolidated clays run at "
    "K0 = 1 and above; under-confining a thin structural column (a reinforced-soil "
    "block) lowers its factor of safety.")

TENSION_SRF_HELP = (
    "Reduce the tensile-strength cap along with c and tan(phi) during an SSRM. "
    "YES (the shipped setting) makes the factor of safety the factor on the WHOLE "
    "strength envelope, shear and tensile — RS2's tensilestrength_SRF = 1 and what "
    "Plaxis does. NO holds each cap at its authored value through the bisection. "
    "Blank leaves it to the solver, which reduces. Either way this is a no-op on a "
    "model that sets no t_cut and no global cutoff: there is no cap to reduce.")


class GlobalParamsDialog(QDialog):
    """Global-parameters form with the v18 Units and Time selectors.

    The two dropdowns sit above the numeric fields; picking a Units system autofills
    the unit-weight-of-water field with that system's canonical value (9.81 SI /
    62.4 Imperial) — typeable-over, so a seawater/brine override still stands. A
    blank Units selection leaves the entered value untouched. Unit-suffix field
    labels are a later phase and are deliberately not added here.
    """

    def __init__(self, numeric_fields, values, parent=None):
        from xslope.units import normalize_unit_system
        super().__init__(parent)
        self.setWindowTitle("Global parameters")
        self._numeric_fields = numeric_fields
        self._values = dict(values)
        self._edits = {}
        layout = QVBoxLayout(self)
        form = QFormLayout()

        # Units selector — canonical si/imperial/None <-> SI/Imperial/blank label.
        self._units = QComboBox()
        for lbl, _val in _UNIT_SYSTEM_ITEMS:
            self._units.addItem(lbl)
        cur_sys = normalize_unit_system(values.get("unit_system"))
        self._units.setCurrentIndex(
            next(i for i, (_l, v) in enumerate(_UNIT_SYSTEM_ITEMS) if v == cur_sys))
        form.addRow("Units", self._units)

        # Time selector — legend tokens verbatim; an unexpected stored token is
        # preserved as an extra entry rather than silently dropped.
        self._time = QComboBox()
        self._time.addItems(_TIME_UNIT_ITEMS)
        cur_time = values.get("time_unit")
        cur_time = str(cur_time) if cur_time else ""
        if cur_time and cur_time not in _TIME_UNIT_ITEMS:
            self._time.addItem(cur_time)
        self._time.setCurrentText(cur_time)
        form.addRow("Time", self._time)

        # Tension SRF (v19) — a tri-state selector, not a checkbox, because blank
        # (unspecified, engine default) must stay distinguishable from NO.
        self._tension_srf = QComboBox()
        for lbl, _val in _TENSION_SRF_ITEMS:
            self._tension_srf.addItem(lbl)
        _cur_srf = values.get("tension_srf")
        self._tension_srf.setCurrentIndex(
            next(i for i, (_l, v) in enumerate(_TENSION_SRF_ITEMS) if v is _cur_srf))
        self._tension_srf.setToolTip(TENSION_SRF_HELP)
        form.addRow("Tension SRF (FEM)", self._tension_srf)

        # Water loads (main!D23) -- auto derives the ponded-water load from
        # the model's own water definition; manual leaves it on the dloads
        # sheets. The loader always resolves the mode to one of the two, so
        # the combo carries no blank entry.
        self._water_loads = QComboBox()
        self._water_loads.addItems(["auto", "manual"])
        _cur_wl = str(values.get("water_loads") or "auto").strip().lower()
        if _cur_wl not in ("auto", "manual"):
            _cur_wl = "auto"
        self._water_loads.setCurrentText(_cur_wl)
        self._water_loads.setToolTip(WATER_LOADS_HELP)
        form.addRow("Water loads", self._water_loads)

        # Numeric fields (gamma_water first — the Units autofill targets it). The
        # header carries a unit suffix ("Unit weight of water (pcf)") when the project
        # declares a system; unlabeled otherwise. Reflects the STORED declaration at
        # open time; it does not live-update when the selector changes (a later pass).
        unit_labels = _unit_labels_for(values)
        for f in numeric_fields:
            _v = values.get(f.key, f.default)
            edit = QLineEdit(f.to_text(_v))
            edit.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            if f.tooltip:
                edit.setToolTip(f.tooltip)
            self._edits[f.key] = edit
            form.addRow(f.display_header(unit_labels), edit)

        layout.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

        # Wire the autofill only AFTER the initial combo state is set, so opening the
        # dialog on an existing project never clobbers its stored gamma_w override.
        self._units.currentIndexChanged.connect(self._autofill_gamma)

    def _autofill_gamma(self, index):
        from xslope.units import GAMMA_W
        system = (_UNIT_SYSTEM_ITEMS[index][1]
                  if 0 <= index < len(_UNIT_SYSTEM_ITEMS) else None)
        edit = self._edits.get("gamma_water")
        if system is not None and edit is not None:
            edit.setText(_display_number(GAMMA_W[system]))

    def result_values(self):
        out = {f.key: f.read_text(self._edits[f.key].text(),
                                  self._values.get(f.key, f.default))
               for f in self._numeric_fields}
        out["unit_system"] = _UNIT_SYSTEM_ITEMS[self._units.currentIndex()][1]
        t = self._time.currentText().strip()
        out["time_unit"] = t or None
        out["tension_srf"] = _TENSION_SRF_ITEMS[self._tension_srf.currentIndex()][1]
        out["water_loads"] = self._water_loads.currentText()
        return out


class GlobalEditor(CategoryEditor):
    label = "Global parameters"
    FIELDS = [
        Field("gamma_water", "Unit weight of water", unit="unit_weight"),
        Field("tcrack_depth", "Tension crack depth"),
        Field("tcrack_water", "Water in crack"),
        Field("k_seismic", "Seismic coefficient k"),
        # v19: K0 initial stress. optfloat so a blank stays None — "unspecified",
        # which is the gravity turn-on initialization, NOT K0 = 0.
        Field("k0", "K0 initial stress (FEM)", "optfloat", usage="fem",
              tooltip=K0_HELP),
    ]

    def build(self, slope_data, parent):
        return GlobalParamsDialog(self.FIELDS, slope_data, parent)

    def apply(self, slope_data, dlg):
        slope_data.update(dlg.result_values())
        # Rebuild the derived tension-crack line (loader does the same from depth).
        from shapely.geometry import LineString
        td = slope_data.get("tcrack_depth", 0)
        gs = slope_data.get("ground_surface")
        if td and td > 0 and gs is not None and not gs.is_empty:
            slope_data["tcrack_surface"] = LineString([(x, y - td) for x, y in gs.coords])
        else:
            slope_data["tcrack_surface"] = None


def _new_material():
    # Same key set load_slope_data produces (gsat blank -> None; the rest zeroed).
    return {"name": "", "gamma": 0.0, "gamma_sat": None, "option": "mc",
            "c": 0.0, "phi": 0.0, "cp": 0.0, "r_elev": 0.0, "d": 0.0, "psi": 0.0,
            "pow_a": 0.0, "pow_b": 0.0, "pow_c": 0.0, "pow_d": 0.0,
            "hb_sci": 0.0, "hb_gsi": 0.0, "hb_mi": 0.0, "hb_d": 0.0,
            "u": "none", "ru": 0.0,
            "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0, "sigma_cp": 0.0,
            "sigma_d": 0.0, "sigma_psi": 0.0, "k1": 0.0, "k2": 0.0, "alpha": 0.0,
            "unsat": "lf", "kr0": 0.0, "h0": 0.0, "vg_a": 0.0, "vg_n": 0.0,
            "Ss": None, "Sy": None,
            "t_cut": None, "phi_b": None, "s_cap": None,
            "E": 0.0, "nu": 0.0}


# --------------------------------------------------------------------------- #
# Materials editor — Table view (the Excel-mirroring table) + List view (a
# per-material form with strength/kr confirmation plots). Both bind to the SAME
# underlying rows, so switching mid-edit is lossless.
# --------------------------------------------------------------------------- #

# Remembered within the session (module-global); first open defaults to Table.
_LAST_MATERIALS_VIEW = "table"

# Base fields that carry a paired sigma_* (shown beside the value under Reliability).
_MAT_SIGMA = {"gamma": "sigma_gamma", "c": "sigma_c", "phi": "sigma_phi",
              "cp": "sigma_cp", "d": "sigma_d", "psi": "sigma_psi"}
# Strength option -> the option-specific fields shown in the list-view form (only
# the selected option's are visible; blank shows none, per bc0999d).
_MAT_OPTION_FIELDS = {"mc": ["c", "phi"], "cp": ["c", "cp", "r_elev"],
                      "pow": ["pow_a", "pow_b", "pow_c", "pow_d"],
                      "hb": ["hb_sci", "hb_gsi", "hb_mi", "hb_d"],
                      "elastic": [], "": []}
_MAT_ALL_OPTION_FIELDS = ["c", "phi", "cp", "r_elev", "pow_a", "pow_b", "pow_c",
                          "pow_d", "hb_sci", "hb_gsi", "hb_mi", "hb_d"]

# v17: the matric-suction pair phi_b/s_cap (red "LEM & FEM" block, cols Q:R right of
# ru) is read only for an effective-stress strength option (mc/pow/hb) paired with a
# signed-u pore-pressure model (piezo or seep). It is inert — and grayed, mirroring
# the mat-sheet CF on cols Q:R — for cp (total-stress Su already embodies field
# suction), for elastic (no strength at all), and whenever u delivers no signed
# pressure (none/ru/blank: ru is a positive ratio by construction, none/blank carry
# no suction). Both solvers now read it (in the FEM/SSRM the suction term is reduced
# by F alongside tan φ'); the graying rule is unchanged by that.
_MAT_SUCTION_DIM = frozenset(["phi_b", "s_cap"])

# v16: an option=elastic material cannot fail — the mat sheet grays every strength
# column, t_cut, the dilation pair, the matric-suction pair, u/ru and the strength
# standard deviations for such a row (CF ranges F..N, Q..R, AB..AF on $E="elastic").
# g/gsat, E/ν, s(γ) and the seepage block stay live. The editor mirrors that: these
# keys read-only/gray on both views when the row's option is elastic.
_MAT_ELASTIC_DIM = frozenset(
    _MAT_ALL_OPTION_FIELDS + ["d", "psi", "phi_b", "s_cap", "t_cut", "u", "ru",
                              "sigma_c", "sigma_phi", "sigma_cp",
                              "sigma_d", "sigma_psi"])


def _mat_dim_keys(row):
    """Field keys to gray for a material row.

    An elastic row grays its whole inert set. Otherwise the only option-coupled
    graying is the matric-suction pair phi_b/s_cap, which is inert for a total-stress
    strength option (cp) or a pore-pressure model that delivers no signed u — i.e. u
    not in {piezo, seep} — mirroring the template's L:M conditional formatting."""
    opt = str(row.get("option", "") or "").strip().lower()
    if opt == "elastic":
        return _MAT_ELASTIC_DIM
    u = str(row.get("u", "") or "").strip().lower()
    if opt == "cp" or u not in ("piezo", "seep"):
        return _MAT_SUCTION_DIM
    return frozenset()


# Unsaturated model -> its curve params. gard reuses vg's a/n columns (fileio.py).
_MAT_UNSAT_FIELDS = {"lf": ["kr0", "h0"], "vg": ["vg_a", "vg_n"],
                     "gard": ["vg_a", "vg_n"]}
_MAT_ALL_UNSAT_FIELDS = ["kr0", "h0", "vg_a", "vg_n"]


class _MaterialColorState:
    """Working per-material display-color overrides for the Materials editor, shared
    by both views and committed to the document's style delta on OK.

    Seeded from the incoming sparse style delta's material *color* entries; stays
    sparse — an index only appears when its color differs from the tab10 palette
    default (so "reset to default" simply drops the key). Keyed by material INDEX
    (slot), exactly as ``xslope.style`` / the canvas resolve material color
    (``str(mat_id)`` == index), so an override stays with the *slot*, not the
    material, when materials are added/removed/reordered — matching the Styles
    dialog and the Inputs plot. Hatch/alpha (owned by the Styles dialog) are not
    touched here; the OK merge carries them through untouched."""

    def __init__(self, style):
        from matplotlib.colors import to_hex
        self._over = {}          # idx -> hex string
        for k, ov in ((style or {}).get("materials") or {}).items():
            if isinstance(ov, dict) and ov.get("color"):
                try:
                    self._over[int(k)] = to_hex(ov["color"])
                except (ValueError, TypeError):
                    pass

    def default_hex(self, idx):
        from matplotlib.colors import to_hex
        from xslope.plot import get_material_color
        return to_hex(get_material_color(idx))

    def resolved_hex(self, idx):
        return self._over.get(idx) or self.default_hex(idx)

    def has_override(self, idx):
        return idx in self._over

    def set(self, idx, hex_color):
        """Record a color for ``idx``, staying sparse: a value equal to the palette
        default drops the override (that's how "Default" in the picker resets)."""
        if not hex_color or hex_color.lower() == self.default_hex(idx).lower():
            self._over.pop(idx, None)
        else:
            self._over[idx] = hex_color

    def reset(self, idx):
        self._over.pop(idx, None)


def _material_swatch(idx, color=None):
    """A color swatch QIcon for material ``idx``. When ``color`` (a hex / mpl color)
    is given it is drawn directly — the editor passes the pending override or the
    resolved default so a list entry tracks in-dialog color edits. Otherwise it falls
    back to the SAME palette the canvas zones use (``xslope.style.material_style`` ->
    tab10 fallback), so a list entry reads as the same color as its zone on the
    Inputs plot."""
    from matplotlib.colors import to_rgba
    if color is None:
        from xslope.style import material_style, resolve_style
        color = material_style(resolve_style(None), idx)["color"]
    r, g, b, a = to_rgba(color)
    pm = QPixmap(16, 16)
    pm.fill(Qt.transparent)
    from PySide6.QtGui import QPainter, QPen
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing, True)
    p.setBrush(QColor.fromRgbF(r, g, b, a))
    p.setPen(QPen(QColor("#666666")))
    p.drawRoundedRect(1, 1, 13, 13, 3, 3)
    p.end()
    return QIcon(pm)


class _MaterialListView(QWidget):
    """List view for the materials editor: three splitter panes —
    [materials list | property form | confirmation plots]. Binds to a list of
    material dicts; edits write through to those dicts immediately, so
    ``result_rows`` (and a view switch) never lose an in-progress edit. All 36
    loader keys have a widget (hidden ones keep their value), so a round-trip
    preserves every field exactly as the table view does.

    Two independent rules decide what the form shows, combined in
    :meth:`_apply_visibility`: the MODEL rules (only the selected strength option's
    fields, ru only for u='ru', only the selected unsat model's curve parameters)
    and the USAGE toggles (a field whose ``applies`` set names no ticked analysis is
    not part of the problem being edited). The usage half mirrors the table's column
    filter field for field — same ``Field.applies`` tags, so the two views hide the
    same parameters — and a group whose every field hides, hides whole."""

    def __init__(self, fields, rows, new_row, reliability_on, parent=None,
                 color_state=None, unit_labels=None, enabled_usage=None):
        super().__init__(parent)
        self._field_by_key = {f.key: f for f in fields}
        self._unit_labels = unit_labels
        self._new_row = new_row
        self._rows = [dict(r) for r in rows]
        self._reliability_on = bool(reliability_on)
        # Ticked analyses (the toggle bar's set). None = built without a filter:
        # everything the analyses cover is shown, σ per the reliability flag.
        self._usage = (set(enabled_usage) if enabled_usage is not None
                       else {"lem", "fem", "seep"} | ({"rel"} if reliability_on else set()))
        self._color = color_state if color_state is not None else _MaterialColorState({})
        self._loading = False     # suppress color write-through while populating a row
        self._cur = -1
        self._edits = {}          # key -> QLineEdit / QComboBox
        self._edit_keys = {}      # focusable widget -> key (help-strip resolver)
        self._cell_widgets = {}   # key -> the labeled cell QWidget (for show/hide)
        self._cell_keys = {}      # the cell QWidget -> key (pair bookkeeping)
        self._model_hidden = set()   # keys the model rules hide (option/u/unsat)
        self._pairs = []          # (side-by-side wrapper, [keys]) — hides when empty
        self._group_boxes = []    # (QGroupBox, [keys]) — hides when every field hides
        self._sigma_widgets = {}  # base key -> (sigma label, sigma edit)

        self._plot_timer = QTimer(self)
        self._plot_timer.setSingleShot(True)
        self._plot_timer.setInterval(160)     # debounce so plots don't lag typing
        self._plot_timer.timeout.connect(self._refresh_plots)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self._build_list_pane())
        splitter.addWidget(self._build_form_pane())
        splitter.addWidget(self._build_plots_pane())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 1)
        splitter.setSizes([200, 470, 530])

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(splitter)

        self._refresh_list()
        self._update_sigma_visibility()
        self._apply_visibility()
        self._apply_plot_visibility()
        if self._rows:
            self.list.setCurrentRow(0)
        else:
            self._load(-1)

    # --- panes -----------------------------------------------------------
    def _build_list_pane(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._on_select)
        v.addWidget(self.list, 1)
        bar = QHBoxLayout()
        add = QPushButton("Add")
        add.clicked.connect(self._add)
        rem = QPushButton("Remove")
        rem.clicked.connect(self._remove)
        bar.addWidget(add)
        bar.addWidget(rem)
        v.addLayout(bar)
        return w

    def _make_edit(self, key):
        f = self._field_by_key[key]
        if f.kind in Field.COMBO_KINDS:
            w = QComboBox()
            w.addItems(f.choices)
        else:
            w = QLineEdit()
            if f.kind in Field.NUMERIC_KINDS:
                w.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            w.editingFinished.connect(self._on_edit)
        self._edits[key] = w
        self._edit_keys[w] = key      # reverse lookup for the help strip
        return w

    def help_key_for_widget(self, widget):
        """The material field key the focused ``widget`` edits, or None. Climbs a
        few parent levels so a combo/color-button focus proxy still resolves."""
        w = widget
        for _ in range(4):
            if w is None:
                break
            key = self._edit_keys.get(w)
            if key is not None:
                return key
            w = w.parentWidget()
        return None

    def _cell(self, label, key, sigma=False, label_w=80):
        """A [label | value | (± σ | sigma value)] cell as a QWidget, registered by
        key for show/hide + load/commit. Widths are kept tight so even the widest
        row (a d|ψ pair, both carrying σ) fits the form pane without a horizontal
        scrollbar."""
        cell = QWidget()
        h = QHBoxLayout(cell)
        h.setContentsMargins(0, 2, 0, 2)
        h.setSpacing(4)
        lab = QLabel(label)
        lab.setMinimumWidth(label_w)
        h.addWidget(lab)
        edit = self._make_edit(key)
        edit.setMinimumWidth(44)
        h.addWidget(edit, 1)
        f = self._field_by_key.get(key)
        tip = getattr(f, "tooltip", None) if f is not None else None
        if tip:
            lab.setToolTip(tip)
            edit.setToolTip(tip)
        if sigma and key in _MAT_SIGMA:
            skey = _MAT_SIGMA[key]
            slab = QLabel("± σ")
            slab.setStyleSheet(f"color:{USAGE_COLOR['rel']};")
            sedit = self._make_edit(skey)
            sedit.setFixedWidth(48)
            h.addWidget(slab)
            h.addWidget(sedit)
            self._sigma_widgets[key] = (slab, sedit)
        self._register_cell(key, cell)
        return cell

    def _register_cell(self, key, cell):
        """Record a form cell under ``key`` for the show/hide rules (both
        directions: the pair/group bookkeeping asks a widget for its key)."""
        self._cell_widgets[key] = cell
        self._cell_keys[cell] = key
        return cell

    def _pair(self, cell_a, cell_b):
        """Two cells side by side, tracked so the row itself hides when both of its
        cells do — an empty half-height wrapper is still a gap in the form."""
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(cell_a, 1)
        h.addWidget(cell_b, 1)
        self._pairs.append((w, [self._cell_keys.get(cell_a),
                                self._cell_keys.get(cell_b)]))
        return w

    def _build_color_row(self):
        """The Identity display-color control: a swatch button (same picker as the
        Styles dialog, seeded with the resolved override-or-default) beside a small
        Reset button that drops the override so the palette default shows through."""
        from .styles_dialog import ColorButton, MATERIAL_PALETTE
        cell = QWidget()
        h = QHBoxLayout(cell)
        h.setContentsMargins(0, 2, 0, 2)
        h.setSpacing(4)
        lab = QLabel("Display color")
        lab.setMinimumWidth(80)
        lab.setToolTip("Color of this material's zone on the Inputs plot. Stored as a "
                       "style override, not a 'mat' property.")
        h.addWidget(lab)
        self._color_btn = ColorButton("#000000", palette=MATERIAL_PALETTE)
        self._edit_keys[self._color_btn] = "_swatch"    # help-strip resolver
        self._color_btn.colorChanged.connect(self._on_color_changed)
        h.addWidget(self._color_btn)
        reset = QPushButton("Reset")
        reset.setToolTip("Remove the override — show the default palette color.")
        reset.clicked.connect(self._on_color_reset)
        h.addWidget(reset)
        h.addStretch(1)
        # Registered under the help strip's swatch key: it carries no field, so it
        # has no usage tags and is always shown (like the material's name).
        return self._register_cell("_swatch", cell)

    def _on_color_changed(self, hexc):
        if self._loading or not (0 <= self._cur < len(self._rows)):
            return
        self._color.set(self._cur, hexc)
        self._refresh_list_item(self._cur)   # track the change in the list swatch

    def _on_color_reset(self):
        if not (0 <= self._cur < len(self._rows)):
            return
        self._color.reset(self._cur)
        self._loading = True
        self._color_btn.blockSignals(True)
        self._color_btn.set_hex(self._color.default_hex(self._cur))
        self._color_btn.blockSignals(False)
        self._loading = False
        self._refresh_list_item(self._cur)

    def _build_form_pane(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        body = QWidget()
        v = QVBoxLayout(body)

        # Identity — name plus the display color (swatch button + reset), the same
        # override the Styles dialog edits and the Inputs plot resolves.
        g = QGroupBox("Identity")
        gv = QVBoxLayout(g)
        gv.addWidget(self._cell("Name", "name"))
        gv.addWidget(self._build_color_row())
        self._group_boxes.append((g, ["name", "_swatch"]))
        v.addWidget(g)

        # Unit weights (γ | γsat side by side; γ carries its σ)
        g = QGroupBox("Unit weights")
        gv = QVBoxLayout(g)
        gv.addWidget(self._pair(self._cell("γ", "gamma", sigma=True, label_w=28),
                                self._cell("γ_sat", "gamma_sat", label_w=44)))
        self._group_boxes.append((g, ["gamma", "gamma_sat"]))
        v.addWidget(g)

        # Strength: option combo, then only the selected option's fields, then
        # dilation d/ψ and elastic E/ν.
        g = QGroupBox("Strength")
        gv = QVBoxLayout(g)
        opt_cell = QWidget()
        oh = QHBoxLayout(opt_cell)
        oh.setContentsMargins(0, 2, 0, 2)
        olab = QLabel("Model (option)")
        olab.setMinimumWidth(96)
        oh.addWidget(olab)
        self._opt_combo = self._make_edit("option")
        self._opt_combo.currentIndexChanged.connect(self._on_option_changed)
        oh.addWidget(self._opt_combo, 1)
        gv.addWidget(self._register_cell("option", opt_cell))
        for key in _MAT_ALL_OPTION_FIELDS:
            gv.addWidget(self._cell(self._label_for(key), key,
                                    sigma=key in _MAT_SIGMA))
        gv.addSpacing(4)
        gv.addWidget(self._pair(self._cell("d", "d", sigma=True, label_w=52),
                                self._cell("ψ", "psi", sigma=True, label_w=18)))
        # v16: tensile cutoff (file order t_cut, E, nu), then the elastic E/ν pair.
        gv.addWidget(self._cell("t_cut", "t_cut", label_w=52))
        gv.addWidget(self._pair(self._cell("E", "E", label_w=18),
                                self._cell("ν", "nu", label_w=18)))
        # Mixed-band group: the option fields and d/ψ are LEM(+FEM), t_cut/E/ν FEM,
        # so per-field hiding decides the group rather than one tag.
        self._group_boxes.append((g, ["option"] + _MAT_ALL_OPTION_FIELDS
                                  + ["d", "psi", "t_cut", "E", "nu"]))
        v.addWidget(g)

        # Pore pressure (ru only when u = 'ru')
        g = QGroupBox("Pore pressure")
        gv = QVBoxLayout(g)
        u_cell = QWidget()
        uh = QHBoxLayout(u_cell)
        uh.setContentsMargins(0, 2, 0, 2)
        ulab = QLabel("Model (u)")
        ulab.setMinimumWidth(96)
        uh.addWidget(ulab)
        self._u_combo = self._make_edit("u")
        self._u_combo.currentIndexChanged.connect(self._on_u_changed)
        uh.addWidget(self._u_combo, 1)
        self._u_cell = u_cell            # grayed whole for an elastic material
        gv.addWidget(self._register_cell("u", u_cell))
        gv.addWidget(self._cell("ru", "ru"))
        # v17: matric-suction pair (file order phi_b, s_cap; right of ru — LEM & FEM).
        # Placed here to mirror the template, where the pair moved right of ru and is
        # coupled to the pore-pressure model (matric suction = negative pore pressure).
        # Always shown, grayed when inert (cp/elastic, or u not signed) by
        # _update_suction_disable.
        gv.addWidget(self._pair(self._cell("φ_b", "phi_b", label_w=52),
                                self._cell("s_cap", "s_cap", label_w=44)))
        self._group_boxes.append((g, ["u", "ru", "phi_b", "s_cap"]))
        v.addWidget(g)

        # Conductivity: k1/k2/alpha, then unsat model + its curve params.
        g = QGroupBox("Conductivity")
        gv = QVBoxLayout(g)
        gv.addWidget(self._pair(self._cell("k1", "k1", label_w=22),
                                self._cell("k2", "k2", label_w=22)))
        gv.addWidget(self._cell("alpha", "alpha"))
        us_cell = QWidget()
        ush = QHBoxLayout(us_cell)
        ush.setContentsMargins(0, 2, 0, 2)
        uslab = QLabel("Unsat model")
        uslab.setMinimumWidth(96)
        ush.addWidget(uslab)
        self._unsat_combo = self._make_edit("unsat")
        self._unsat_combo.currentIndexChanged.connect(self._on_unsat_changed)
        ush.addWidget(self._unsat_combo, 1)
        gv.addWidget(self._register_cell("unsat", us_cell))
        for key in _MAT_ALL_UNSAT_FIELDS:
            gv.addWidget(self._cell(self._label_for(key), key))
        self._group_boxes.append((g, ["k1", "k2", "alpha", "unsat"]
                                  + _MAT_ALL_UNSAT_FIELDS))
        v.addWidget(g)

        v.addStretch(1)
        scroll.setWidget(body)
        return scroll

    def _build_plots_pane(self):
        from .canvas import MplCanvas
        pane = QSplitter(Qt.Vertical)
        self._strength_canvas = MplCanvas()
        self._kr_canvas = MplCanvas()
        pane.addWidget(self._strength_canvas)
        pane.addWidget(self._kr_canvas)
        pane.setSizes([320, 320])
        self._plots_pane = pane
        return pane

    # List view shows friendly symbols where the table mirrors the terse 'mat'
    # sheet headers (Norm: "use φ in the list view"). Only display labels — the
    # underlying keys/headers are untouched.
    _FRIENDLY = {"f": "φ", "psi": "ψ", "s(f)": "σ(φ)", "s(g)": "σ(γ)",
                 "s(c)": "σ(c)", "s(c/p)": "σ(c/p)", "s(d)": "σ(d)",
                 "s(psi)": "σ(ψ)",
                 # The Hoek-Brown four, in the notation the criterion is published
                 # in: σci, GSI, mi, D. The hb_ prefix disambiguates the sheet
                 # columns from the power-curve ones; on a labeled field it only
                 # gets in the way of recognizing the parameter.
                 "hb_sci": "σci", "hb_gsi": "GSI", "hb_mi": "mi", "hb_d": "D"}

    def _label_for(self, key):
        f = self._field_by_key.get(key)
        header = f.header if f is not None else key
        # Friendly symbol first (φ/ψ/σ), then the unit suffix. No field is both
        # friendly-mapped AND unit-bearing, so order is immaterial; a declared unit
        # yields e.g. "g (pcf)", undeclared leaves the bare symbol unchanged.
        base = self._FRIENDLY.get(header, header)
        return _with_unit(base, f, self._unit_labels)

    # --- list ------------------------------------------------------------
    def _item_text(self, i):
        name = self._rows[i].get("name") or ""
        return f"{i + 1}.  {name}" if name else f"{i + 1}."

    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i in range(len(self._rows)):
            item = QListWidgetItem(_material_swatch(i, self._color.resolved_hex(i)),
                                   self._item_text(i))
            self.list.addItem(item)
        self.list.blockSignals(False)

    def _refresh_list_item(self, i):
        if 0 <= i < self.list.count():
            it = self.list.item(i)
            it.setText(self._item_text(i))
            it.setIcon(_material_swatch(i, self._color.resolved_hex(i)))

    # --- load / commit ---------------------------------------------------
    def _load(self, idx):
        self._cur = -1                       # suppress write-through while populating
        for key, w in self._edits.items():
            val = self._rows[idx].get(key) if (0 <= idx < len(self._rows)) else None
            if isinstance(w, QComboBox):
                w.blockSignals(True)
                txt = self._field_by_key[key].to_text(val)
                j = w.findText(txt)
                w.setCurrentIndex(j if j >= 0 else 0)
                w.blockSignals(False)
            else:
                w.blockSignals(True)
                # NaN is "unset" throughout the templates (e.g. t_res, E, area on
                # lines that don't carry them) — render blank, never the string "nan"
                # (to_text does that, and rounds the rest for display).
                w.setText(self._field_by_key[key].to_text(val))
                w.blockSignals(False)
        ok = 0 <= idx < len(self._rows)
        self.setEnabled(True)
        self._loading = True
        self._color_btn.blockSignals(True)
        if ok:
            self._color_btn.set_default(self._color.default_hex(idx))
            self._color_btn.set_hex(self._color.resolved_hex(idx))
        self._color_btn.setEnabled(ok)
        self._color_btn.blockSignals(False)
        self._loading = False
        if ok:
            self._cur = idx
            self._update_option_visibility()
            self._update_u_visibility()
            self._update_unsat_visibility()
            self._update_sigma_visibility()
            self._update_elastic_disable()
            self._update_suction_disable()
            self._refresh_plots()
        else:
            self._clear_plots()

    def _commit(self):
        if not (0 <= self._cur < len(self._rows)):
            return
        row = self._rows[self._cur]
        for key, w in self._edits.items():
            f = self._field_by_key[key]
            if isinstance(w, QComboBox):
                row[key] = f.from_text(w.currentText())
            else:
                # read_text: an untouched field re-commits the value it was loaded
                # with, so the display rounding is never written back (Field.read_text).
                row[key] = f.read_text(w.text(), row.get(key))

    # --- visibility ------------------------------------------------------
    def _usage_visible(self, key):
        """Is ``key`` part of a ticked analysis? Universal cells (the name, the color
        swatch — no ``applies`` tags) are always part of the problem; every other
        field answers with the same tags the table filters its columns by."""
        f = self._field_by_key.get(key)
        if f is None or f.applies is None:
            return True
        return bool(f.applies & self._usage)

    def _set_model_hidden(self, key, hidden):
        """Record the MODEL rules' verdict on ``key`` (the selected strength option /
        pore-pressure model / unsaturated model). Kept apart from the usage verdict so
        neither rule can un-hide what the other hid."""
        if hidden:
            self._model_hidden.add(key)
        else:
            self._model_hidden.discard(key)

    def _apply_visibility(self):
        """Show a cell only where both rules agree, then fold away the side-by-side
        rows and the groups left with nothing in them.

        isHidden, not isVisible: this runs on a form that has not been shown yet
        (the list view is built before the stack switches to it), where isVisible()
        is False for every cell and would collapse the whole form."""
        for key, cell in self._cell_widgets.items():
            cell.setVisible(key not in self._model_hidden
                            and self._usage_visible(key))
        for w, keys in self._pairs:
            w.setVisible(any(not self._cell_widgets[k].isHidden()
                             for k in keys if k in self._cell_widgets))
        for g, keys in self._group_boxes:
            g.setVisible(any(not self._cell_widgets[k].isHidden()
                             for k in keys if k in self._cell_widgets))

    def _apply_plot_visibility(self):
        """Gate the two confirmation plots on the bands they confirm: the strength
        envelope on LEM/FEM, the kr curve on Seepage. With neither band ticked the
        whole pane goes — a splitter drops a hidden pane's width rather than leaving
        a blank half-dialog."""
        strength = bool(self._usage & {"lem", "fem"})
        kr = "seep" in self._usage
        self._strength_canvas.setVisible(strength)
        self._kr_canvas.setVisible(kr)
        self._plots_pane.setVisible(strength or kr)

    def _update_option_visibility(self):
        opt = self._opt_combo.currentText().strip().lower()
        shown = set(_MAT_OPTION_FIELDS.get(opt, []))
        for key in _MAT_ALL_OPTION_FIELDS:
            self._set_model_hidden(key, key not in shown)
        self._apply_visibility()

    def _update_elastic_disable(self):
        """Gray the fields inert for an option=elastic material (mirrors the mat-sheet
        conditional formatting): t_cut, the dilation pair d/ψ and the pore-pressure
        model u/ru go read-only; g/gsat, E/ν and the seepage block stay live. The
        strength-option cells are already hidden for 'elastic' (empty option set), so
        only the always-shown cells need graying. Disabling a whole d/ψ cell also
        grays its s(d)/s(ψ) — matching the template, which grays those too."""
        elastic = self._opt_combo.currentText().strip().lower() == "elastic"
        for key in ("t_cut", "d", "psi", "ru"):
            w = self._cell_widgets.get(key)
            if w is not None:
                w.setEnabled(not elastic)
        self._u_cell.setEnabled(not elastic)

    def _update_suction_disable(self):
        """Gray the matric-suction pair phi_b/s_cap when it is inert (mirrors the
        mat-sheet CF on cols Q:R): grayed for cp (total-stress Su already embodies
        field suction) or elastic (no strength), and whenever u delivers no signed
        pressure (u not in {piezo, seep}). Reuses the table view's dim rule so the two
        views stay a single source of truth. Depends on both the option and u combos,
        so it is re-run whenever either changes as well as on load."""
        dim = _mat_dim_keys({"option": self._opt_combo.currentText(),
                             "u": self._u_combo.currentText()})
        for key in ("phi_b", "s_cap"):
            w = self._cell_widgets.get(key)
            if w is not None:
                w.setEnabled(key not in dim)

    def _update_u_visibility(self):
        self._set_model_hidden(
            "ru", self._u_combo.currentText().strip().lower() != "ru")
        self._apply_visibility()

    def _update_unsat_visibility(self):
        um = self._unsat_combo.currentText().strip().lower()
        shown = set(_MAT_UNSAT_FIELDS.get(um, []))
        for key in _MAT_ALL_UNSAT_FIELDS:
            self._set_model_hidden(key, key not in shown)
        self._apply_visibility()

    def _update_sigma_visibility(self):
        for _key, (slab, sedit) in self._sigma_widgets.items():
            slab.setVisible(self._reliability_on)
            sedit.setVisible(self._reliability_on)

    def apply_usage_filter(self, enabled):
        """Show only the fields, groups and plots the ticked analyses use — the
        list-view half of the toggle bar's contract, mirroring the table's column
        filter. Hidden cells keep their widgets and their values, so toggling never
        drops data on save. Reliability rides along as it always has: its σ readouts
        follow the 'rel' tag."""
        self._usage = set(enabled)
        self._reliability_on = "rel" in self._usage
        self._update_sigma_visibility()
        self._apply_visibility()
        self._apply_plot_visibility()
        self._refresh_plots()

    # --- events ----------------------------------------------------------
    def _on_select(self, idx):
        self._commit()
        self._load(idx)

    def _on_edit(self):
        self._commit()
        self._refresh_list_item(self._cur)   # name edits update the list label
        self._plot_timer.start()

    def _on_option_changed(self):
        self._commit()
        self._update_option_visibility()
        self._update_elastic_disable()
        self._update_suction_disable()       # option gates the suction pair too
        self._plot_timer.start()

    def _on_u_changed(self):
        self._commit()
        self._update_u_visibility()          # u affects neither plot; no refresh
        self._update_suction_disable()       # u gates the matric-suction pair

    def _on_unsat_changed(self):
        self._commit()
        self._update_unsat_visibility()
        self._plot_timer.start()             # kr plot depends on the unsat model

    def _add(self):
        self._commit()
        self._rows.append(self._new_row())
        self._cur = -1
        self._refresh_list()
        self.list.setCurrentRow(len(self._rows) - 1)

    def _remove(self):
        idx = self.list.currentRow()
        if not (0 <= idx < len(self._rows)):
            return
        self._rows.pop(idx)
        self._cur = -1
        self._refresh_list()
        if self._rows:
            self.list.setCurrentRow(min(idx, len(self._rows) - 1))
        else:
            self._load(-1)

    # --- plots -----------------------------------------------------------
    def _refresh_plots(self):
        from xslope.plot import plot_material_strength, plot_material_kr
        if not (0 <= self._cur < len(self._rows)):
            self._clear_plots()
            return
        mat = dict(self._rows[self._cur])    # snapshot for the deferred draw
        self._strength_canvas.render_axes(lambda ax: plot_material_strength(ax, mat))
        self._kr_canvas.render_axes(lambda ax: plot_material_kr(ax, mat))

    def _clear_plots(self):
        self._strength_canvas.clear()
        self._kr_canvas.clear()

    # --- external API (used by MaterialsDialog on view switch) -----------
    def set_rows(self, rows, reliability_on):
        self._reliability_on = bool(reliability_on)
        self._rows = [dict(r) for r in rows]
        self._cur = -1
        self._refresh_list()
        if self._rows:
            self.list.setCurrentRow(0)
        else:
            self._load(-1)

    def result_rows(self):
        self._commit()
        return [dict(r) for r in self._rows]


# Single source of truth for the Materials editor's field help — one short,
# semantics-first line per field, wording aligned with the v16 'mat' worksheet
# legend (docs/usage/input_template.md) and the existing t_cut tooltip. Shown in
# the context-sensitive help strip (both views) and reused as the t_cut hover
# tooltip. Keyed by the material field keys plus '_swatch' for the display-color
# control. Kept concise so the two-line strip never clips.
MATERIALS_HELP = {
    "name": "Material name — identity only; also labels this material's zone on the Inputs plot.",
    "gamma": "Moist/total unit weight (above the piezometric line).",
    "gamma_sat": "Saturated unit weight (below the piezometric line; blank = use g).",
    "option": "Strength model — determines which strength columns apply.",
    "c": "Cohesion intercept (mc: Mohr-Coulomb c; cp: undrained strength at r-elev).",
    "phi": "Friction angle φ, degrees (mc option).",
    "cp": "Rate of undrained-strength increase per unit elevation below r-elev (cp option).",
    "r_elev": "Reference elevation; at or above it the strength equals c (cp option).",
    "d": "Rapid-drawdown R-envelope cohesion intercept (LEM rapid-drawdown only).",
    "psi": "Rapid-drawdown R-envelope friction angle (LEM rapid-drawdown only).",
    "phi_b": ("Fredlund unsaturated friction angle φ_b, degrees — credits matric-"
              "suction strength above the water line (mc/pow/hb; LEM & FEM — in the "
              "FEM/SSRM the suction term is reduced by F alongside tan φ'). Blank = no "
              "suction strength (the default). Caution: with u=piezo hydrostatic "
              "suction is unbounded — set s_cap; with u=seep the FE field self-bounds."),
    "s_cap": ("Cap on credited matric suction, stress units; blank = uncapped. "
              "Essential with u=piezo (suction grows unbounded above the line); a "
              "backstop with u=seep, where the FE field self-bounds. LEM & FEM. "
              "Paired with phi_b."),
    # 393 chars — MEASURED to wrap in exactly two lines at the dialog's natural width
    # (the strip is fixed at two lines and clips beyond). Keep any edit at or under
    # this length.
    "t_cut": ("Tensile-strength cap (Rankine, major principal stress). FEM only. Blank "
              "is NOT 'no tension': Mohr-Coulomb grants c/tanφ (28 kPa at c=20, φ=35°; "
              "unbounded at φ=0), which strength reduction never touches — enough to "
              "hold a steep crest cut shut and inflate SSRM FS. A cap you set IS reduced "
              "with F, like c and tanφ (RS2/Plaxis do the same). 0 = no tension (can "
              "block reinforced-fill equilibrium)."),
    "E": "FEM elastic (Young's) modulus — with ν, the only mechanical property read for an elastic material.",
    "nu": "FEM Poisson's ratio — with E, the only mechanical property read for an elastic material.",
    "u": "Pore-pressure model: none, piezo (piezometric line), seep (seepage solution), or ru.",
    "ru": "Pore-pressure ratio; u = ru·σv on the soil column only (u = ru option).",
    "pow_a": "Power-curve envelope coefficient a in τ = a·(σ'n + d)^b + c (pow option).",
    "pow_b": "Power-curve exponent b; b = 1 collapses to Mohr-Coulomb (pow option).",
    "pow_c": "Power-curve additive strength term c (pow option).",
    "pow_d": "Power-curve normal-stress offset d (pow option).",
    "hb_sci": "σci — uniaxial compressive strength of the intact rock, in the model's "
              "stress units: 30 MPa is 30,000 kPa (hb option).",
    "hb_gsi": "GSI — Geological Strength Index, (0, 100] (hb option).",
    "hb_mi": "mi — intact Hoek-Brown constant, a rock-type property (hb option).",
    "hb_d": "D — disturbance factor, 0 (undisturbed) to 1 (blast-damaged) (hb option).",
    "sigma_gamma": "Reliability: standard deviation of g.",
    "sigma_c": "Reliability: standard deviation of c.",
    "sigma_phi": "Reliability: standard deviation of φ.",
    "sigma_cp": "Reliability: standard deviation of cp.",
    "sigma_d": "Reliability: standard deviation of d.",
    "sigma_psi": "Reliability: standard deviation of ψ.",
    "k1": "Seepage: major hydraulic conductivity (k1 = kx when alpha = 0).",
    "k2": "Seepage: minor hydraulic conductivity (k2 = ky when alpha = 0).",
    "alpha": "Seepage: orientation angle of the permeability tensor (degrees).",
    "unsat": "Seepage: unsaturated relative-permeability model — lf, vg, or gard.",
    "kr0": "Relative conductivity kr0 at suction head h0 (linear-front, unsat = lf).",
    "h0": "Suction head at which k = kr0 (linear-front, unsat = lf).",
    "vg_a": "Curve parameter a (vg: α in 1/length; gard: power-form a).",
    "vg_n": "Curve parameter n (van Genuchten / Gardner, unsat = vg or gard).",
    "Ss": "Specific storage — water released per unit volume of saturated soil per "
          "unit head drop (1/length). Read by a transient seepage run only; "
          "required on every material then, blank otherwise.",
    "Sy": "Specific yield — the drainable porosity: the volume fraction released "
          "as the water table falls (dimensionless). Read by a transient seepage "
          "run only; required on every material then, blank otherwise.",
    "_swatch": "Display color on the Inputs plot — a style override, not a 'mat' property.",
}


class MaterialsDialog(QDialog):
    """The materials editor dialog: a segmented Table/List view toggle over a
    QStackedWidget. Both views bind to ``self._rows`` (the single source of truth):
    switching harvests the active view's edits into it and rebuilds the target, so
    a switch mid-edit is lossless. The Reliability toggle drives σ columns (table)
    and σ fields (list) in both views."""

    def __init__(self, title, fields, rows, new_row, parent=None, help_text=None,
                 usage_toggles=None, style=None, doc=None, unit_labels=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._title = title
        self._fields = fields
        self._new_row = new_row
        self._unit_labels = unit_labels
        self._rows = [dict(r) for r in rows]
        self._mode = None
        self._table = None
        self._list_view = None
        # Display-color overrides live in the project style delta (NOT the material
        # dicts). Seed a working state from the incoming delta, shared by both views;
        # result_style() folds it back into the delta on OK. doc is the document to
        # commit to (None in headless/round-trip use — colors just aren't committed).
        self._doc = doc
        self._orig_style = copy.deepcopy(style or {})
        self._color = _MaterialColorState(self._orig_style)
        # The width the LIST view was laid out for, and the height both views open
        # at. The table view sizes itself from its own fitted columns instead (see
        # _fit_table_width), which on a filtered material sheet is narrower.
        self._designed_width = 1180
        self.resize(self._designed_width, 640)

        layout = QVBoxLayout(self)
        if help_text:
            layout.addWidget(_help_label(help_text))

        top = QHBoxLayout()
        seg_group = QButtonGroup(self)
        seg_group.setExclusive(True)
        self._seg = {}
        for mode, lbl in (("table", "Table view"), ("list", "List view")):
            b = QPushButton(lbl)
            b.setCheckable(True)
            b.clicked.connect(lambda _checked=False, m=mode: self._set_mode(m))
            seg_group.addButton(b)
            self._seg[mode] = b
            top.addWidget(b)
        top.addSpacing(24)
        top.addWidget(QLabel("Show parameters for:"))
        from PySide6.QtCore import QSettings
        s = QSettings("XSlope", "XSlope Studio")
        self._toggles = {}
        for t in (usage_toggles or []):
            cb = QCheckBox(USAGE_TOGGLE_LABEL[t])
            default = (t != "rel")
            cb.setChecked(bool(s.value(f"editor_toggles/{self._title}/{t}",
                                       default, type=bool)))
            cb.setStyleSheet(f"color:{USAGE_COLOR[t]}; font-weight:bold;")
            if t == "rel":
                cb.setToolTip("Reliability rides on the LEM/FEM parameters — "
                              "tick LEM or FEM to show each value's σ.")
            cb.toggled.connect(self._on_toggle)
            self._toggles[t] = cb
            top.addWidget(cb)
        top.addStretch(1)
        layout.addLayout(top)
        self._sync_rel_enabled()

        self._stack = QStackedWidget()
        self._table_holder = QWidget()
        self._table_lay = QVBoxLayout(self._table_holder)
        self._table_lay.setContentsMargins(0, 0, 0, 0)
        self._list_holder = QWidget()
        self._list_lay = QVBoxLayout(self._list_holder)
        self._list_lay.setContentsMargins(0, 0, 0, 0)
        self._stack.addWidget(self._table_holder)     # index 0
        self._stack.addWidget(self._list_holder)      # index 1
        layout.addWidget(self._stack, 1)

        _ok_cancel(self, layout)

        # Context-sensitive help strip: one source-of-truth mapping, one resolver
        # that names the field under focus in whichever view is active. Attached
        # before the first _set_mode so _build_table can wire its cell-tracking.
        attach_help(self, MATERIALS_HELP, self._help_key_for)

        initial = _LAST_MATERIALS_VIEW if _LAST_MATERIALS_VIEW in ("table", "list") else "table"
        self._set_mode(initial)

    # --- context-sensitive help ------------------------------------------
    def _help_key_for(self, widget):
        """The material field key the focused ``widget`` belongs to, or None.

        Table view: derive the column from the focused cell/cell-widget (falling
        back to the current column) and map it to its field. List view: defer to
        the list view's own widget->key lookup."""
        if self._mode == "table" and self._table is not None:
            col = self._table.column_at(widget)
            return self._help_key_for_col(col)
        if self._mode == "list" and self._list_view is not None:
            return self._list_view.help_key_for_widget(widget)
        return None

    def _help_key_for_col(self, col):
        """Field key for a logical table column, or '_swatch' for the trailing
        display-color column, or None for no/invalid column."""
        if col is None or col < 0:
            return None
        if 0 <= col < len(self._fields):
            return self._fields[col].key
        return "_swatch"        # the appended (logical-last) swatch column

    def _help_from_col(self, col):
        """Push the help for a table column into the strip (for currentCellChanged,
        which fires as arrow keys move the current cell without changing focus)."""
        strip = getattr(self, "_help_strip", None)
        if strip is not None:
            key = self._help_key_for_col(col)
            strip.set_help(MATERIALS_HELP.get(key, "") if key else "")

    # --- toggles ---------------------------------------------------------
    def _rel_applies(self):
        """Does Reliability have anything to qualify? A σ is the scatter of an
        LEM/FEM parameter (γ, c, φ, c/p, d, ψ), so with neither of those analyses
        ticked there is nothing for it to apply to."""
        return any(self._toggles[t].isChecked()
                   for t in ("lem", "fem") if t in self._toggles)

    def _sync_rel_enabled(self):
        """Gray the Reliability toggle while it has nothing to qualify. Its checked
        state is left alone — untouched in the widget and in QSettings — so ticking
        LEM or FEM back on restores the σ readouts exactly as they were."""
        cb = self._toggles.get("rel")
        if cb is not None:
            cb.setEnabled(self._rel_applies())

    def _reliability_on(self):
        cb = self._toggles.get("rel")
        return bool(cb is not None and cb.isChecked() and self._rel_applies())

    def _enabled_usage(self):
        on = {t for t, cb in self._toggles.items() if cb.isChecked()}
        if not self._rel_applies():
            on.discard("rel")
        return on

    def _on_toggle(self):
        from PySide6.QtCore import QSettings
        s = QSettings("XSlope", "XSlope Studio")
        for t, cb in self._toggles.items():
            s.setValue(f"editor_toggles/{self._title}/{t}", cb.isChecked())
        self._sync_rel_enabled()
        if self._mode == "table" and self._table is not None:
            self._table.apply_usage_filter(self._enabled_usage())
            # Columns that were just shown belong on the dialog, not past its right
            # edge — so the width follows the filter (wider only).
            self._fit_table_width()
        if self._mode == "list" and self._list_view is not None:
            self._list_view.apply_usage_filter(self._enabled_usage())

    # --- view switching --------------------------------------------------
    def _harvest(self):
        if self._mode == "table" and self._table is not None:
            self._rows = self._table.result_rows()
        elif self._mode == "list" and self._list_view is not None:
            self._rows = self._list_view.result_rows()

    def _build_table(self):
        if self._table is not None:
            self._table.setParent(None)
            self._table.deleteLater()
        self._table = _EditableTable(self._fields, self._rows, self._new_row,
                                     swatch_state=self._color, dim_rule=_mat_dim_keys,
                                     unit_labels=self._unit_labels)
        self._table_lay.addWidget(self._table)
        self._table.apply_usage_filter(self._enabled_usage())
        # Track the current cell so the help strip follows keyboard navigation
        # across columns (focus stays on the table, so focusChanged won't fire).
        if getattr(self, "_help_strip", None) is not None:
            self._table.table.currentCellChanged.connect(
                lambda cr, cc, pr, pc: self._help_from_col(cc))

    def _ensure_list(self):
        if self._list_view is None:
            self._list_view = _MaterialListView(self._fields, self._rows,
                                                self._new_row, self._reliability_on(),
                                                color_state=self._color,
                                                unit_labels=self._unit_labels,
                                                enabled_usage=self._enabled_usage())
            self._list_lay.addWidget(self._list_view)
        else:
            # Filter first: set_rows loads a material, and the load re-runs the
            # form's show/hide against whatever the toggles say now.
            self._list_view.apply_usage_filter(self._enabled_usage())
            self._list_view.set_rows(self._rows, self._reliability_on())

    def _set_mode(self, mode):
        if mode not in ("table", "list"):
            return
        if mode == self._mode:
            self._seg[mode].setChecked(True)
            return
        self._harvest()                      # pull current edits into self._rows
        strip = getattr(self, "_help_strip", None)
        if strip is not None:
            strip.set_help("")               # drop the other view's stale help text
        if mode == "table":
            self._build_table()
            self._stack.setCurrentIndex(0)
            self._fit_table_width()
        else:
            self._ensure_list()
            self._stack.setCurrentIndex(1)
            if self.isVisible():
                _grow_dialog_to(self, self._designed_width)
        self._mode = mode
        self._seg[mode].setChecked(True)
        global _LAST_MATERIALS_VIEW
        _LAST_MATERIALS_VIEW = mode

    def _fit_table_width(self):
        """Take the dialog to the width the table view's fitted columns need, up to
        the width the list view was laid out for — the two views share one dialog,
        and the material sheet has more columns than any screen has room for, so the
        table view's fit is a floor to open at rather than a width to insist on.
        Only ever wider once the dialog is on screen (see _fit_dialog_to_columns)."""
        if self._table is None:
            return
        _fit_dialog_to_columns(self, self._table, cap=self._designed_width,
                               grow_only=self.isVisible())

    def set_view_mode(self, mode):
        """Programmatic view switch (used by the round-trip guard)."""
        self._set_mode(mode)

    def result_rows(self):
        self._harvest()
        return self._rows

    def result_style(self):
        """The project style delta with material *display colors* reconciled from the
        editor's working overrides. Sparse: an override is written only where the
        color differs from the palette default; hatch/alpha (owned by the Styles
        dialog) are carried through untouched. Slot-keyed by material index, so
        entries for slots past the current material count are dropped. Applying with
        no color change returns a delta equal to the one passed in."""
        from matplotlib.colors import to_hex
        n = len(self.result_rows())
        out = copy.deepcopy(self._orig_style)
        mats = dict(out.get("materials") or {})
        new_mats = {}
        for idx in range(n):
            key = str(idx)
            entry = dict(mats.get(key) or {})
            if self._color.has_override(idx):
                new_hex = self._color.resolved_hex(idx)
                old = entry.get("color")
                # Keep the original color string when it's the same color (avoids
                # churning e.g. "red" -> "#ff0000" on an untouched override).
                if old is None or to_hex(old).lower() != new_hex.lower():
                    entry["color"] = new_hex
            else:
                entry.pop("color", None)
            if entry:
                new_mats[key] = entry
        if new_mats:
            out["materials"] = new_mats
        else:
            out.pop("materials", None)
        return out


class MaterialsEditor(CategoryEditor):
    label = "Materials"
    # Columns mirror the v17 'mat' worksheet IN FILE ORDER: name, g, gsat, option,
    # c, f, c/p, r-elev, d, psi, t_cut, E, nu, u, ru, phi_b, s_cap, pow_a..pow_d,
    # hb_sci/hb_gsi/hb_mi/hb_d, s(g), s(c), s(f), s(c/p), s(d), s(psi), k1, k2,
    # alpha, unsat, kr0, h0, a(vg_a), n(vg_n).  v16 inserted t_cut and moved E/nu up
    # beside the strength values, with u/ru following; v17's matric-suction pair
    # phi_b/s_cap now sits at cols Q/R, right of ru (relocated there when the feature
    # landed in the FEM too — the red "LEM & FEM" block, same class as c/φ).
    # `applies` tags mirror the template's analysis usage (input_template.md):
    # the strength block (g, gsat, option, c, f, c/p, r-elev, the pow_* power-curve
    # and hb_* Hoek-Brown envelope parameters, u, ru, phi_b/s_cap) is shared by LEM+FEM;
    # d/psi are rapid-drawdown (LEM); s(...) reliability; k1..vg_n seepage; t_cut/E/nu FEM.
    # gsat is optional (blank -> fall back to g), so it reads back as None when
    # left empty rather than 0.0. The alternate-envelope columns (pow_*/hb_*) are
    # always shown; option-driven show/hide grouping is a later UX pass.
    #
    # Every property whose ZERO is physically meaningful is an 'optfloat' -- gamma,
    # c, phi, d, psi, E, nu, ru, k1 -- so a blank cell comes back as None and the
    # model checks can tell "the user entered 0" from "the user left it blank".
    # A plain 'float' stamps 0.0 onto every blank cell of every row on OK (the table
    # rewrites all fields), which made a cohesionless material and an unfilled one
    # indistinguishable, a blank Poisson's ratio silently 0.0 (-33% on a measured
    # SSRM factor of safety), and a blank r_u a model with no pore pressure at all.
    # The writer blanks the cell for a None rather than stamping a 0.0 into the
    # sheet, so an unfilled property is not invented on the way out either. (The
    # LOADER still reads a blank numeric cell back as 0.0, unchanged here; what this
    # buys is that a model edited in Studio carries the distinction while it is being
    # edited, which is when the checks run and when the user can act on them.)
    LF = {"lem", "fem"}
    FIELDS = [
        Field("name", "name", "str"),
        Field("gamma", "γ", "optfloat", applies=LF, unit="unit_weight"),
        Field("gamma_sat", "γsat", "optfloat", applies=LF, unit="unit_weight"),
        # A BLANK option is valid for seep-only material rows (the loader keeps ''
        # via _choice; document._blank_material produces it for DXF imports). Offer
        # it as an empty combo entry so the editor round-trips it instead of
        # normalizing blank -> 'mc'. 'elastic' (v16): infinite strength — the row's
        # strength/t_cut/u cells gray out. Kept last so the default stays 'mc'.
        Field("option", "option", "choice", choices=["mc", "cp", "pow", "hb", "elastic", ""], applies=LF),
        Field("c", "c", "optfloat", applies=LF, unit="stress"),
        Field("phi", "φ", "optfloat", applies=LF),
        Field("cp", "c/p", applies=LF), Field("r_elev", "r-elev", applies=LF),
        Field("d", "d", "optfloat", usage="lem"),
        Field("psi", "psi", "optfloat", usage="lem"),
        # v16: tensile-strength cutoff (FEM only). optfloat so a blank cell stays
        # None (no cutoff), never 0.0 (which would mean "no tension").
        Field("t_cut", "t_cut", "optfloat", usage="fem", tooltip=MATERIALS_HELP["t_cut"]),
        Field("E", "E", "optfloat", usage="fem", unit="stress"),
        Field("nu", "n", "optfloat", usage="fem"),
        Field("u", "u", "choice", choices=["none", "piezo", "seep", "ru"], applies=LF),
        Field("ru", "ru", "optfloat", applies=LF),
        # v17: matric-suction pair (file order phi_b, s_cap; right of ru — the red
        # "LEM & FEM" block, mirroring c/φ, since both solvers now read it). optfloat
        # so a blank cell stays None — phi_b None = no suction strength (the default,
        # exactly the pre-v17 behavior); s_cap None = uncapped suction. In the
        # FEM/SSRM the suction term is reduced by F alongside tan φ'.
        Field("phi_b", "phi_b", "optfloat", applies=LF, tooltip=MATERIALS_HELP["phi_b"]),
        Field("s_cap", "s_cap", "optfloat", applies=LF, tooltip=MATERIALS_HELP["s_cap"],
              unit="stress"),
        Field("pow_a", "pow_a", applies=LF), Field("pow_b", "pow_b", applies=LF),
        Field("pow_c", "pow_c", applies=LF), Field("pow_d", "pow_d", applies=LF),
        Field("hb_sci", "hb_sci", applies=LF), Field("hb_gsi", "hb_gsi", applies=LF),
        Field("hb_mi", "hb_mi", applies=LF), Field("hb_d", "hb_d", applies=LF),
        Field("sigma_gamma", "s(g)", usage="rel"), Field("sigma_c", "s(c)", usage="rel"),
        Field("sigma_phi", "s(f)", usage="rel"), Field("sigma_cp", "s(c/p)", usage="rel"),
        Field("sigma_d", "s(d)", usage="rel"), Field("sigma_psi", "s(psi)", usage="rel"),
        Field("k1", "k1", "optfloat", usage="seep", unit="k"),
        Field("k2", "k2", usage="seep", unit="k"),
        Field("alpha", "alpha", usage="seep"),
        # Unsaturated model: lf (linear front -> kr0/h0), vg (van Genuchten) or gard
        # (Gardner); vg/gard share the vg_a/vg_n curve pair.
        Field("unsat", "unsat", "choice", choices=["lf", "vg", "gard"], usage="seep"),
        Field("kr0", "kr0", usage="seep"), Field("h0", "h0", usage="seep", unit="length"),
        Field("vg_a", "vg_a", usage="seep"), Field("vg_n", "vg_n", usage="seep"),
        # v18: transient storage (mat sheet seepage block, after the unsat curve
        # pair). optfloat so a blank stays None — the loader never defaults these
        # to 0 (a silent zero would drop the storage term), and preflight demands
        # them only when a tseep sheet is in use.
        Field("Ss", "Ss", "optfloat", usage="seep", unit="inv_length",
              tooltip=MATERIALS_HELP["Ss"]),
        Field("Sy", "Sy", "optfloat", usage="seep", tooltip=MATERIALS_HELP["Sy"]),
    ]

    def build(self, slope_data, parent):
        # The display-color swatch edits the project's style delta, which lives on the
        # document (reached via the parent window); headless/round-trip callers pass a
        # bare parent, so both are optional and colors simply aren't committed then.
        doc = getattr(parent, "doc", None)
        style = getattr(doc, "style", None) if doc is not None else None
        return MaterialsDialog(
            "Materials", self.FIELDS, slope_data.get("materials", []), _new_material, parent,
            unit_labels=_unit_labels_for(slope_data),
            help_text="Table view mirrors the 'mat' worksheet (row order = Mat ID "
                      "order). List view edits one material at a time as a form with "
                      "strength- and conductivity-model plots that confirm the "
                      "selected options. Both views edit the same rows, so switching "
                      "is lossless. The toggles hide the parameters an analysis does "
                      "not use — table columns, and list-view fields and plots; "
                      "'Reliability' adds each value's σ. The color "
                      "swatch sets the material's display color on the Inputs plot.",
            usage_toggles=["lem", "seep", "fem", "rel"], style=style, doc=doc)

    def apply(self, slope_data, dlg):
        slope_data["materials"] = dlg.result_rows()
        # Display colors are style deltas, not material data — commit them via the
        # document's style store (same path as the Styles dialog), only when changed,
        # so the sparse delta stays sparse and an unchanged apply is a no-op.
        doc = getattr(dlg, "_doc", None)
        if doc is not None:
            new_style = dlg.result_style()
            if new_style != doc.style:
                doc.set_style(new_style)


# --------------------------------------------------------------------------- #
# Line editors — Table view (the bulk-entry table) + List view (a per-line
# grouped form beside the live section preview). Shared by the reinforcement and
# pile editors: both edit a list of 2-endpoint line records. Both views bind to
# the SAME rows, so switching mid-edit is lossless (mirrors MaterialsDialog).
# --------------------------------------------------------------------------- #

# Last-used view per line editor, remembered within the session. Unlike materials
# (which mirrors the wide 'mat' sheet and opens on the table), these open on the
# LIST view: piles are almost always 1-3 rows, and a tiered wall's reinforcement is
# edited one line at a time — the table stays the bulk-entry path for the 15-20
# lines of a big wall, but the list is the primary editing surface.
_LAST_LINE_VIEW = {"reinforce": "list", "piles": "list"}


def _last_line_view(key):
    return _LAST_LINE_VIEW.get(key, "list")


def _set_last_line_view(key, mode):
    _LAST_LINE_VIEW[key] = mode


class _LineListView(QWidget):
    """List view for the 2-endpoint line editors (reinforcement / piles): three
    splitter panes — [items list | grouped property form | section preview]. Binds
    to a list of row dicts; edits write through immediately, so ``result_rows`` (and
    a view switch) never lose an in-progress edit. Every table field has a widget
    (grouped, all shown), and unshown keys (e.g. a pile's derived θ) ride along
    untouched, so a round-trip preserves every field exactly as the table view does.
    The preview is the SAME debounced ``PreviewPane`` the table view uses; a click in
    it selects the corresponding list item (the emphasis follows)."""

    def __init__(self, fields, rows, new_row, groups, item_label,
                 preview_draw, pick_resolve, preview_caption=None, parent=None,
                 unit_labels=None, dynamic_spec=None, preset_spec=None,
                 switch_spec=None):
        super().__init__(parent)
        self._field_by_key = {f.key: f for f in fields}
        # Optional per-line either/or switch inside one group (reinforcement's
        # Pullout: development length OR overburden). It is a VIEW of the row —
        # which pair carries values decides where it opens — not a stored field,
        # so nothing new reaches the file. See _build_switch.
        self._switch = switch_spec or None
        self._switch_combo = None
        self._switch_stash = {}   # row index -> {key: text} parked by a switch
        # Same preset rule the table view carries (reinforcement's Type -> Dir/Appl),
        # so the two views fill identically — see _apply_preset.
        self._preset = preset_spec or None
        self._unit_labels = unit_labels
        self._new_row = new_row
        self._rows = [dict(r) for r in rows]
        self._groups = groups
        self._item_label = item_label
        self._preview_draw = preview_draw
        self._pick_resolve = pick_resolve
        self._cur = -1
        self._edits = {}          # key -> QLineEdit / QComboBox
        self._edit_keys = {}      # focusable widget -> key (help-strip resolver)
        # Per-element/per-unit-width labeling (reinforcement/pile): the field that
        # drives the wording (Spacing / S), the map of scaled-field key -> quantity
        # ('force'/'moment'/'area'), and the QLabels to re-word live. Empty for the
        # editors that don't scale by spacing.
        spec = dynamic_spec or {}
        self._dyn_driver = spec.get("driver")
        self._dyn_fields = dict(spec.get("fields", {}))
        self._dyn_labels = {}     # key -> (QLabel, quantity)

        from .canvas import PreviewPane
        list_pane = self._build_list_pane()
        form_pane = self._build_form_pane()
        self._preview = PreviewPane(
            lambda ax: self._preview_draw(ax, self._rows, self._cur),
            caption=preview_caption)
        self._preview.clicked.connect(self._on_preview_click)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(list_pane)
        splitter.addWidget(form_pane)
        splitter.addWidget(self._preview)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 1)
        splitter.setSizes([200, 380, 470])

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(splitter)

        self._refresh_list()
        if self._rows:
            self.list.setCurrentRow(0)
        else:
            self._load(-1)
        self._preview.refresh_now()

    # --- panes -----------------------------------------------------------
    def _build_list_pane(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._on_select)
        v.addWidget(self.list, 1)
        bar = QHBoxLayout()
        add = QPushButton("Add")
        add.clicked.connect(self._add)
        rem = QPushButton("Remove")
        rem.clicked.connect(self._remove)
        bar.addWidget(add)
        bar.addWidget(rem)
        v.addLayout(bar)
        return w

    def _make_edit(self, key):
        f = self._field_by_key[key]
        if f.kind in Field.COMBO_KINDS:
            w = QComboBox()
            w.addItems(f.choices)               # blank-tolerant: a "" choice is a real
            if self._preset is not None and key == self._preset["driver"]:
                # Fill the dependent combos FIRST, so the commit below writes the row
                # with them already set. `activated` too: re-picking the same Type is
                # the user re-asserting the preset, and fires no index change.
                w.activated.connect(self._apply_preset)
                w.currentIndexChanged.connect(self._apply_preset)
            w.currentIndexChanged.connect(self._on_edit)   # (empty) entry, as in the table
        else:
            w = QLineEdit()
            if f.kind in Field.NUMERIC_KINDS:
                w.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            w.editingFinished.connect(self._on_edit)
        self._edits[key] = w
        self._edit_keys[w] = key      # reverse lookup for the help strip
        return w

    def help_key_for_widget(self, widget):
        """The field key ``widget`` edits, or None. Climbs a few parent levels so a
        combo's internal focus proxy still resolves (mirrors _MaterialListView)."""
        w = widget
        for _ in range(4):
            if w is None:
                break
            key = self._edit_keys.get(w)
            if key is not None:
                return key
            w = w.parentWidget()
        return None

    def _cell(self, key, label_w=58):
        f = self._field_by_key[key]
        cell = QWidget()
        h = QHBoxLayout(cell)
        h.setContentsMargins(0, 2, 0, 2)
        h.setSpacing(4)
        quantity = self._dyn_fields.get(key)
        if quantity is not None:
            # Spacing-scaled field: a live per-element/per-unit-width label, seeded to
            # the per-unit-width wording and re-worded by _relabel_dynamic once a row
            # (with its Spacing) is loaded.
            lab = QLabel(_dynamic_capacity_label(f.header, quantity,
                                                 self._unit_labels, False))
            self._dyn_labels[key] = (lab, quantity)
        else:
            lab = QLabel(f.display_header(self._unit_labels))
        lab.setMinimumWidth(label_w)
        if f.usage:                             # mirror the table's header coloring
            lab.setStyleSheet(f"color:{USAGE_COLOR[f.usage]};")
        tip = getattr(f, "tooltip", None)
        if tip:
            lab.setToolTip(tip)
        h.addWidget(lab)
        edit = self._make_edit(key)
        edit.setMinimumWidth(56)
        if tip:
            edit.setToolTip(tip)
        h.addWidget(edit, 1)
        self._cells[key] = cell
        return cell

    @staticmethod
    def _pair(cell_a, cell_b):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(cell_a, 1)
        h.addWidget(cell_b, 1)
        return w

    def _build_form_pane(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        body = QWidget()
        v = QVBoxLayout(body)
        self._cells = {}
        self._group_boxes = []
        for title, group_rows in self._groups:
            g = QGroupBox(title)
            gv = QVBoxLayout(g)
            if self._switch is not None and title == self._switch.get("group"):
                gv.addWidget(self._build_switch())
            for keys in group_rows:
                # A spacing-scaled field carries a long "per element / per unit width"
                # label, so a pair containing one is broken onto full-width single rows
                # (its label would otherwise be clipped by the half-width cell).
                if len(keys) == 2 and not any(k in self._dyn_fields for k in keys):
                    gv.addWidget(self._pair(self._cell(keys[0]), self._cell(keys[1])))
                else:
                    for k in keys:
                        gv.addWidget(self._cell(k))
            self._group_boxes.append((g, [k for row in group_rows for k in row]))
            v.addWidget(g)
        v.addStretch(1)
        scroll.setWidget(body)
        self._form_scroll = scroll
        # Re-word the spacing-scaled labels live as the Spacing/S field is typed (its
        # editingFinished already commits the value; textChanged is the immediate one).
        driver = self._edits.get(self._dyn_driver) if self._dyn_driver else None
        if driver is not None and hasattr(driver, "textChanged"):
            driver.textChanged.connect(lambda *_: self._relabel_dynamic())
        return scroll

    # --- either/or switch -------------------------------------------------
    def _build_switch(self):
        """The group's law selector: one combo naming each option and the pair of
        fields it governs.

        Nothing about it is stored. The row itself says which law is in force —
        for reinforcement, a filled Adhesion/Delta pair IS the overburden law —
        so the combo opens where ``resolve`` puts it and the fields it does not
        govern are grayed, values intact.
        """
        spec = self._switch
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 2, 0, 2)
        h.setSpacing(4)
        lab = QLabel(spec.get("label", "Mode"))
        lab.setMinimumWidth(58)
        combo = QComboBox()
        for opt in spec["options"]:
            combo.addItem(opt[1])
        tip = spec.get("tooltip")
        if tip:
            lab.setToolTip(tip)
            combo.setToolTip(tip)
        combo.currentIndexChanged.connect(self._on_switch_changed)
        h.addWidget(lab)
        h.addWidget(combo, 1)
        self._switch_combo = combo
        return w

    def _switch_index(self, row):
        want = self._switch["resolve"](row)
        for j, opt in enumerate(self._switch["options"]):
            if opt[0] == want:
                return j
        return 0

    def _apply_switch_state(self):
        """Enable the governing pair, gray the rest — labels and edits alike, so a
        dimmed field reads as inactive rather than merely unfocused."""
        if self._switch_combo is None:
            return
        active = self._switch_combo.currentIndex()
        for j, opt in enumerate(self._switch["options"]):
            for key in opt[2]:
                cell = self._cells.get(key)
                if cell is not None:
                    cell.setEnabled(j == active)

    def _on_switch_changed(self, *_):
        """Move the row onto the law the combo now names.

        Only a SELF-ASSERTING pair is cleared, and it is parked rather than
        dropped: values left in Adhesion and Delta would keep the overburden law
        in force whatever the combo said, so leaving that option stashes them and
        entering it again brings them back. The development lengths assert
        nothing while the other law runs, so they merely dim — their values stay
        in the cells, which is what the switch promises.
        """
        if self._switch_combo is None or self._cur < 0:
            return
        active = self._switch_combo.currentIndex()
        stash = self._switch_stash.setdefault(self._cur, {})
        for j, opt in enumerate(self._switch["options"]):
            keys, asserting = opt[2], opt[3]
            for key in keys:
                w = self._edits.get(key)
                if w is None or isinstance(w, QComboBox):
                    continue
                if j == active:
                    if not w.text().strip() and stash.get(key):
                        w.setText(stash.pop(key))
                elif asserting:
                    if w.text().strip():
                        stash[key] = w.text()
                    w.clear()
        self._apply_switch_state()
        self._commit()
        self._preview.schedule()

    def apply_usage_filter(self, enabled):
        """Hide the form cells of fields whose usage tag is not enabled — the
        list-view half of the toggle bar's contract, mirroring the table's
        column filter. A group whose every field hides, hides whole."""
        for key, cell in getattr(self, "_cells", {}).items():
            f = self._field_by_key[key]
            cell.setVisible((not f.usage) or (f.usage in enabled))
        for g, keys in getattr(self, "_group_boxes", []):
            # isHidden, not isVisible: the filter runs on a lazily built list view
            # before it is shown, where isVisible() is False for every cell and
            # would hide every group.
            g.setVisible(any(not self._cells[k].isHidden() for k in keys
                             if k in self._cells))

    def _relabel_dynamic(self):
        """Re-word every spacing-scaled label from the current Spacing/S entry: set ->
        'per element', blank/1 -> 'per unit width' (with the declared unit string
        joined when a system is set). A no-op for editors without scaled fields."""
        if not self._dyn_labels:
            return
        driver = self._edits.get(self._dyn_driver)
        text = driver.text() if driver is not None else ""
        per_element = _spacing_is_per_element(text)
        for key, (lab, quantity) in self._dyn_labels.items():
            f = self._field_by_key[key]
            lab.setText(_dynamic_capacity_label(f.header, quantity,
                                                self._unit_labels, per_element))

    # --- list ------------------------------------------------------------
    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i in range(len(self._rows)):
            self.list.addItem(QListWidgetItem(self._item_label(i, self._rows[i])))
        self.list.blockSignals(False)

    def _refresh_list_item(self, i):
        if 0 <= i < self.list.count():
            self.list.item(i).setText(self._item_label(i, self._rows[i]))

    # --- load / commit ---------------------------------------------------
    def _load(self, idx):
        self._cur = -1                       # suppress write-through while populating
        ok = 0 <= idx < len(self._rows)
        for key, w in self._edits.items():
            val = self._rows[idx].get(key) if ok else None
            w.blockSignals(True)
            if isinstance(w, QComboBox):
                txt = self._field_by_key[key].to_text(val)
                j = w.findText(txt)
                w.setCurrentIndex(j if j >= 0 else 0)
            else:
                # NaN is "unset" throughout the templates (e.g. t_res, E, area on
                # lines that don't carry them) — render blank, never the string "nan"
                # (to_text does that, and rounds the rest for display).
                w.setText(self._field_by_key[key].to_text(val))
            w.blockSignals(False)
        # Loading blocks the driver's textChanged, so re-word the scaled labels from
        # the freshly-loaded Spacing/S here.
        self._relabel_dynamic()
        # The either/or switch reads the freshly-loaded row, silently: seeding it
        # is not the user choosing a law, so it must not move any value.
        if self._switch_combo is not None:
            self._switch_combo.blockSignals(True)
            self._switch_combo.setCurrentIndex(
                self._switch_index(self._rows[idx]) if ok else 0)
            self._switch_combo.blockSignals(False)
            self._apply_switch_state()
        self._form_scroll.setEnabled(ok)
        if ok:
            self._cur = idx
        self._preview.schedule()

    def _commit(self):
        if not (0 <= self._cur < len(self._rows)):
            return
        row = self._rows[self._cur]
        for key, w in self._edits.items():
            f = self._field_by_key[key]
            if isinstance(w, QComboBox):
                row[key] = f.from_text(w.currentText())
            else:
                # read_text: an untouched field re-commits the value it was loaded
                # with, so the display rounding is never written back (Field.read_text).
                row[key] = f.read_text(w.text(), row.get(key))

    # --- events ----------------------------------------------------------
    def _on_select(self, idx):
        self._commit()
        self._load(idx)

    def _apply_preset(self, *_):
        """Fill the dependent combos from the preset the driver combo now names —
        the table view's ``_apply_preset_row``, on the one row this view shows.

        Picking a Type SETS Dir and Appl; editing either afterwards keeps what was
        entered, since nothing rewrites them until a Type is picked again. A driver
        value with no preset (blank Type) fills nothing. Silent while a row loads:
        ``_load`` blocks every widget's signals, so switching lines never re-fills a
        Dir/Appl the file overrode."""
        if self._preset is None or self._cur < 0:
            return
        driver = self._edits.get(self._preset["driver"])
        if driver is None:
            return
        preset = self._preset["presets"].get(driver.currentText().strip().lower())
        if preset is None:
            return
        for key, value in zip(self._preset["fills"], preset):
            w = self._edits.get(key)
            if isinstance(w, QComboBox):
                j = w.findText(value)
                if j >= 0:
                    w.setCurrentIndex(j)     # its own edit signal commits the row

    def _on_edit(self, *_):
        self._commit()
        self._refresh_list_item(self._cur)   # label/type edits update the list line
        self._preview.schedule()

    def _add(self):
        self._commit()
        self._rows.append(self._new_row())
        self._cur = -1
        self._refresh_list()
        self.list.setCurrentRow(len(self._rows) - 1)

    def _remove(self):
        idx = self.list.currentRow()
        if not (0 <= idx < len(self._rows)):
            return
        self._rows.pop(idx)
        self._cur = -1
        self._refresh_list()
        if self._rows:
            self.list.setCurrentRow(min(idx, len(self._rows) - 1))
        else:
            self._load(-1)

    def _on_preview_click(self, x, y, tol):
        """A click in the preview: resolve it to a row and select that list item (the
        selection path then re-renders the preview with the new emphasis)."""
        row = self._pick_resolve(x, y, tol, self._rows)
        if row is not None and 0 <= row < self.list.count():
            self.list.setCurrentRow(row)

    # --- external API (used by _LineEditorDialog on view switch) ---------
    def set_rows(self, rows):
        self._rows = [dict(r) for r in rows]
        self._cur = -1
        self._refresh_list()
        if self._rows:
            self.list.setCurrentRow(0)
        else:
            self._load(-1)
        self._preview.refresh_now()

    def result_rows(self):
        self._commit()
        return [dict(r) for r in self._rows]


class _LineEditorDialog(QDialog):
    """Table/List view toggle for the 2-endpoint line editors (reinforcement, piles),
    mirroring ``MaterialsDialog``. Both views bind to ``self._rows`` (the single source
    of truth): switching harvests the active view's edits into it and rebuilds the
    target, so a switch mid-edit is lossless. Each view carries the SAME kind of
    cross-section ``PreviewPane`` (right of the table; third pane of the list view)
    with click-to-select. The usage toggles show/hide table columns; the list view
    always shows every field, so they don't affect it (as in materials).

    These default to the LIST view (materials keeps its table default), remembering the
    last-used view for the session — see ``_LAST_LINE_VIEW``."""

    def __init__(self, title, fields, rows, new_row, groups, item_label,
                 preview_draw, pick_resolve, view_state, parent=None,
                 help_text=None, usage_toggles=None, preview_caption=None,
                 field_help=None, unit_labels=None, dynamic_spec=None,
                 preset_spec=None, dim_rule=None, switch_spec=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._title = title
        self._fields = fields
        self._unit_labels = unit_labels
        # Per-element/per-unit-width labeling spec for the list view (see
        # _LineListView); None for editors whose fields aren't scaled by spacing.
        self._dynamic_spec = dynamic_spec
        # Preset rule handed to BOTH views, so a Type picked in either fills the same
        # two columns from the same table (reinforcement; None elsewhere).
        self._preset_spec = preset_spec
        # Per-row graying rule (row dict -> keys to gray) and the list view's
        # either/or switch. Both express the SAME fact — which of two mutually
        # exclusive parameter sets a row is using — so the table's grayed pair and
        # the list's dimmed pair always agree.
        self._dim_rule = dim_rule
        self._switch_spec = switch_spec
        self._new_row = new_row
        self._groups = groups
        self._item_label = item_label
        self._preview_draw = preview_draw
        self._pick_resolve = pick_resolve
        self._preview_caption = preview_caption
        self._view_state = view_state
        # Field-key -> help-text mapping for the context-sensitive help strip
        # (see attach_help); shared by both the table and list views.
        self._field_help = field_help
        self._rows = [dict(r) for r in rows]
        self._mode = None
        self._table = None
        self._table_split = None
        self._table_preview = None
        self._list_view = None
        # The width the LIST view was laid out for, and the height both views open
        # at. The table view sizes itself from its own fitted columns instead (see
        # _set_mode) — narrower here, since a line's fields are a dozen numbers —
        # but never wider than the view it shares the dialog with.
        self._designed_width = 1200
        self.resize(self._designed_width, 620)

        layout = QVBoxLayout(self)
        if help_text:
            layout.addWidget(_help_label(help_text))

        top = QHBoxLayout()
        seg_group = QButtonGroup(self)
        seg_group.setExclusive(True)
        self._seg = {}
        for mode, lbl in (("table", "Table view"), ("list", "List view")):
            b = QPushButton(lbl)
            b.setCheckable(True)
            b.clicked.connect(lambda _checked=False, m=mode: self._set_mode(m))
            seg_group.addButton(b)
            self._seg[mode] = b
            top.addWidget(b)
        top.addSpacing(24)
        from PySide6.QtCore import QSettings
        s = QSettings("XSlope", "XSlope Studio")
        self._toggles = {}
        if usage_toggles:
            top.addWidget(QLabel("Show parameters for:"))
        for t in (usage_toggles or []):
            cb = QCheckBox(USAGE_TOGGLE_LABEL[t])
            cb.setChecked(bool(s.value(f"editor_toggles/{self._title}/{t}",
                                       True, type=bool)))
            cb.setStyleSheet(f"color:{USAGE_COLOR[t]}; font-weight:bold;")
            cb.toggled.connect(self._on_toggle)
            self._toggles[t] = cb
            top.addWidget(cb)
        top.addStretch(1)
        layout.addLayout(top)

        self._stack = QStackedWidget()
        self._table_holder = QWidget()
        self._table_lay = QVBoxLayout(self._table_holder)
        self._table_lay.setContentsMargins(0, 0, 0, 0)
        self._list_holder = QWidget()
        self._list_lay = QVBoxLayout(self._list_holder)
        self._list_lay.setContentsMargins(0, 0, 0, 0)
        self._stack.addWidget(self._table_holder)     # index 0
        self._stack.addWidget(self._list_holder)      # index 1
        layout.addWidget(self._stack, 1)

        _ok_cancel(self, layout)

        # Attached before the first _set_mode so _build_table can wire its
        # cell-tracking (mirrors MaterialsDialog).
        if self._field_help is not None:
            attach_help(self, self._field_help, self._help_key_for)

        last = _last_line_view(self._view_state)
        self._set_mode(last if last in ("table", "list") else "list")

    # --- context-sensitive help ------------------------------------------
    def _help_key_for(self, widget):
        """The field key the focused ``widget`` belongs to, or None.

        Table view: derive the column from the focused cell (falling back to the
        current column) and map it to its field. List view: defer to the list
        view's own widget->key lookup (the widget-map pattern)."""
        if self._mode == "table" and self._table is not None:
            col = self._table.column_at(widget)
            return self._help_key_for_col(col)
        if self._mode == "list" and self._list_view is not None:
            return self._list_view.help_key_for_widget(widget)
        return None

    def _help_key_for_col(self, col):
        if col is None or col < 0 or col >= len(self._fields):
            return None
        return self._fields[col].key

    def _help_from_col(self, col):
        """Push the help for a table column into the strip (for currentCellChanged,
        which fires as arrow keys move the current cell without changing focus)."""
        strip = getattr(self, "_help_strip", None)
        if strip is not None:
            key = self._help_key_for_col(col)
            strip.set_help((self._field_help or {}).get(key, "") if key else "")

    # --- toggles ---------------------------------------------------------
    def _enabled_usage(self):
        return {t for t, cb in self._toggles.items() if cb.isChecked()}

    def _on_toggle(self):
        from PySide6.QtCore import QSettings
        s = QSettings("XSlope", "XSlope Studio")
        for t, cb in self._toggles.items():
            s.setValue(f"editor_toggles/{self._title}/{t}", cb.isChecked())
        if self._table is not None:
            self._table.apply_usage_filter(self._enabled_usage())
            if self._mode == "table":
                # Columns that were just shown belong on the dialog, not past its
                # right edge — so the width follows the filter (wider only).
                self._fit_table_width()
        if self._list_view is not None:
            self._list_view.apply_usage_filter(self._enabled_usage())

    # --- view switching --------------------------------------------------
    def _harvest(self):
        if self._mode == "table" and self._table is not None:
            self._rows = self._table.result_rows()
        elif self._mode == "list" and self._list_view is not None:
            self._rows = self._list_view.result_rows()

    def _schedule_table_preview(self, *_):
        if self._table_preview is not None:
            self._table_preview.schedule()

    def _on_table_preview_click(self, x, y, tol):
        row = self._pick_resolve(x, y, tol, self._table.result_rows())
        if row is not None:
            self._table.select_row(row)

    def _build_table(self):
        from .canvas import PreviewPane
        if self._table_split is not None:      # rebuild from the shared rows
            self._table_split.setParent(None)
            self._table_split.deleteLater()
        self._table = _EditableTable(self._fields, self._rows, self._new_row,
                                     on_change=self._schedule_table_preview,
                                     on_select=self._schedule_table_preview,
                                     unit_labels=self._unit_labels,
                                     preset_spec=self._preset_spec,
                                     dim_rule=self._dim_rule,
                                     dim_on_edit=self._dim_rule is not None)
        self._table_preview = PreviewPane(
            lambda ax: self._preview_draw(ax, self._table.result_rows(),
                                          self._table.selected_row()),
            caption=self._preview_caption)
        self._table_preview.clicked.connect(self._on_table_preview_click)
        # The preview sits BELOW the table: these editors carry a dozen or more
        # columns, and a side-by-side preview forced the dialog wider than a
        # screen (owner ruling 2026-08-14). The table keeps the height its rows
        # ask for; the preview absorbs the rest, so a taller dialog grows the
        # picture — the circles editor's stacked pattern.
        split = QSplitter(Qt.Vertical)
        split.addWidget(self._table)
        split.addWidget(self._table_preview)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        self._table_preview.setMinimumHeight(220)
        # Open with every row visible: the table's opening share is its own
        # measured hint, not a constant, so six rows show six rows.
        split.setSizes([self._table.sizeHint().height(), 400])
        self._table_split = split
        self._table_lay.addWidget(split)
        self._table.apply_usage_filter(self._enabled_usage())
        self._table_preview.refresh_now()
        # Track the current cell so the help strip follows keyboard navigation
        # across columns (focus stays on the table, so focusChanged won't fire).
        if getattr(self, "_help_strip", None) is not None:
            self._table.table.currentCellChanged.connect(
                lambda cr, cc, pr, pc: self._help_from_col(cc))

    def _ensure_list(self):
        if self._list_view is None:
            self._list_view = _LineListView(
                self._fields, self._rows, self._new_row, self._groups,
                self._item_label, self._preview_draw, self._pick_resolve,
                preview_caption=self._preview_caption, unit_labels=self._unit_labels,
                dynamic_spec=self._dynamic_spec, preset_spec=self._preset_spec,
                switch_spec=self._switch_spec)
            self._list_lay.addWidget(self._list_view)
            # A lazily built list view starts under whatever the toggle bar
            # already says — the same filter the table is showing.
            if getattr(self, "_toggles", None):
                self._on_toggle()
        else:
            self._list_view.set_rows(self._rows)

    def _set_mode(self, mode):
        if mode not in ("table", "list"):
            return
        if mode == self._mode:
            self._seg[mode].setChecked(True)
            return
        self._harvest()                      # pull current edits into self._rows
        strip = getattr(self, "_help_strip", None)
        if strip is not None:
            strip.set_help("")               # drop the other view's stale help text
        if mode == "table":
            self._build_table()
            self._stack.setCurrentIndex(0)
            self._fit_table_width()
        else:
            self._ensure_list()
            self._stack.setCurrentIndex(1)
            if self.isVisible():
                _grow_dialog_to(self, self._designed_width)
        self._mode = mode
        self._seg[mode].setChecked(True)
        _set_last_line_view(self._view_state, mode)

    def _fit_table_width(self):
        """Take the dialog to the width the table view's fitted columns need.

        Bounded above by the width the list view was laid out for: the two views
        share one dialog, and a table that wanted more than that would be sizing the
        other view's form off the screen. Bounded below by nothing but the dialog's
        own controls — a table of a dozen number columns is a narrow dialog, and it
        should open as one.

        Only ever wider once the dialog is on screen: a switch back to the table is
        not a reason to take back width the user gave it."""
        if self._table is None:
            return
        _fit_dialog_to_columns(self, self._table, cap=self._designed_width,
                               grow_only=self.isVisible())

    def set_view_mode(self, mode):
        """Programmatic view switch (used by the round-trip guard)."""
        self._set_mode(mode)

    def result_rows(self):
        self._harvest()
        return self._rows


def _new_circle():
    return {"Xo": 0.0, "Yo": 0.0, "Option": "Depth", "Depth": 0.0,
            "Xi": 0.0, "Yi": 0.0, "R": 0.0}


CIRCLES_HELP = {
    "Xo": "X-coordinate of the circle center.",
    "Yo": "Y-coordinate of the circle center.",
    "Option": "How the circle's size is defined — Depth, Radius, or Intercept; "
              "only the matching field below is used.",
    "Depth": "Depth of the circle bottom below the center (R = Yo − Depth). Option = Depth.",
    "Xi": "X-coordinate of a point the circle passes through. Option = Intercept.",
    "Yi": "Y-coordinate of a point the circle passes through. Option = Intercept.",
    "R": "Circle radius, specified directly. Option = Radius.",
}


SEARCH_WINDOW_HELP = {
    "entry_x_min": "X range the failure surface's crest-side (higher-ground) endpoint "
                   "must fall in. A trial surface breaking outside it is rejected.",
    "exit_x_min": "X range the toe-side (lower-ground) endpoint must fall in. A trial "
                  "surface daylighting outside it is rejected.",
    "center_box_x_min": "Rectangle the circle centers are confined to. The refining "
                        "grid stays inside it, so the search cannot walk out. Applies "
                        "only when all four cells are filled.",
    "center_box_y_min": "Rectangle the circle centers are confined to. The refining "
                        "grid stays inside it, so the search cannot walk out. Applies "
                        "only when all four cells are filled.",
    "max_tangent_depth": "The lowest ELEVATION the circle's bottom (its tangent point) "
                         "may reach. A deeper trial surface is rejected.",
    "min_slip_depth": "Minimum depth below the ground surface a surface must reach. "
                      "Rejects shallow surficial “skin” mechanisms, whose factor "
                      "of safety is depth-independent on a cohesionless face and would "
                      "otherwise win.",
}
# Each grid row of the group: (label, key or None, key or None). The four ranges pair
# their two ends on one row; the two single limits fill the min column alone.
_SEARCH_WINDOW_ROWS = (
    ("Entry x", "entry_x_min", "entry_x_max"),
    ("Exit x", "exit_x_min", "exit_x_max"),
    ("Center box x", "center_box_x_min", "center_box_x_max"),
    ("Center box y", "center_box_y_min", "center_box_y_max"),
    ("Max tangent depth", "max_tangent_depth", None),
    ("Min slip depth", "min_slip_depth", None),
)
#: The four ranges, as (low key, high key) — what has to be increasing, and what has
#: to be filled at both ends to apply at all.
_SEARCH_WINDOW_PAIRS = tuple((lo, hi) for _lbl, lo, hi in _SEARCH_WINDOW_ROWS
                             if hi is not None)


class _SearchWindowGroup(QGroupBox):
    """The circles sheet's optional search window (J8:K17), edited as a group.

    Ten independent limits confining an automated circular search, every one of them
    optional: a blank field is a limit that is not applied, and an all-blank group is
    the unconstrained search. That is exactly what the loader produces for a blank
    block — no ``search_window`` key at all — so :meth:`result_values` returns only
    the filled keys and an empty dict for a group nobody filled in, and the editor
    drops the key rather than storing ten Nones.

    Each field is an ``optfloat`` :class:`Field`, so a value shown rounded for reading
    is written back at full precision unless it was actually typed over."""

    CAPTION = ("Optional limits confining an automated circular search. A blank field "
               "is a limit that is not applied; a range applies only when BOTH of its "
               "ends are filled.")

    def __init__(self, window, unit_labels=None, parent=None):
        super().__init__("Search window", parent)
        self._stored = dict(window or {})
        self._on_change = None
        self._fields, self._edits = {}, {}
        outer = QVBoxLayout(self)
        outer.addWidget(_help_label(self.CAPTION))
        grid = QGridLayout()
        grid.addWidget(QLabel("min"), 0, 1)
        grid.addWidget(QLabel("max"), 0, 2)
        for r, (label, lo, hi) in enumerate(_SEARCH_WINDOW_ROWS, start=1):
            tip = SEARCH_WINDOW_HELP[lo]
            name = QLabel(_with_unit(label, Field(lo, label, unit="length"), unit_labels))
            name.setToolTip(tip)
            grid.addWidget(name, r, 0)
            for col, key in ((1, lo), (2, hi)):
                if key is None:
                    continue
                field = Field(key, label, "optfloat", unit="length", tooltip=tip)
                edit = QLineEdit(field.to_text(self._stored.get(key)))
                edit.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
                edit.setToolTip(tip)
                edit.textChanged.connect(self._changed)
                self._fields[key] = field
                self._edits[key] = edit
                grid.addWidget(edit, r, col)
        grid.setColumnStretch(3, 1)
        outer.addLayout(grid)

    def set_on_change(self, callback):
        """Call ``callback`` whenever a field is typed in — how the dialog's live
        preview follows the window being edited."""
        self._on_change = callback

    def _changed(self, _text=None):
        if self._on_change is not None:
            self._on_change()

    def result_values(self):
        """The window as the loader would produce it: only the filled keys, in sheet
        order, each a float."""
        out = {}
        for key in SEARCH_WINDOW_KEYS:
            field = self._fields[key]
            value = field.read_text(self._edits[key].text(), self._stored.get(key))
            if value is not None:
                out[key] = float(value)
        return out

    def validate(self):
        """The message to refuse OK with, or "" when the window is savable.

        The one thing the loader rejects outright is a range that runs backwards, and
        a file saved with one would not open again — so it is caught here, where the
        user can still see which pair they typed."""
        window = self.result_values()
        for lo, hi in _SEARCH_WINDOW_PAIRS:
            if lo in window and hi in window and window[lo] > window[hi]:
                return (f"The search window has {lo} = {window[lo]:g} greater than "
                        f"{hi} = {window[hi]:g}. Every range must be increasing "
                        f"(min ≤ max).")
        return ""


def _circles_generate_spec(slope_data):
    """The Generate button's state and behaviour for the starting-circles editor.

    The generator is ``xslope.search.generate_starting_circles`` — the same function
    the Run LEM preflight offers as its ``generate_starting_circles`` remedy, so the
    button and the remedy cannot propose different circles. Whether it CAN run is
    settled before the dialog is drawn, from the model's own geometry: a section with
    no room for a circle to daylight leaves a dimmed button carrying that reason
    rather than a live one that fails when pressed."""
    from xslope.search import generate_starting_circles

    spec = {"label": "Generate starting circles…",
            "unit": "circle",
            # Each circle is its own trial surface, so adding to the ones already in
            # the table is a real thing to want: a generated set alongside a circle
            # the user typed is a wider search, not a contradiction.
            "append": True,
            "tooltip": "Derive a starting set from the slope's own geometry: for each "
                       "significant face, a circle through the toe and one at the base "
                       "of each layer, centered above the middle of the face at twice "
                       "the slope height."}
    probe = generate_starting_circles(slope_data, report=True)
    if not probe["circles"]:
        spec["available"] = False
        spec["reason"] = (f"No starting circle can be derived from this model: "
                          f"{probe['reason']}.")
        return spec
    spec["available"] = True

    def propose(parent):
        # Re-run against the model as it stands rather than reusing the probe above.
        result = generate_starting_circles(slope_data, report=True)
        if not result["circles"]:
            return [], "", result["reason"]
        rows = []
        for circle in result["circles"]:
            row = _new_circle()
            row.update(circle)
            # The generator sizes circles by the elevation of their lowest point,
            # which is what Option = Depth means on this sheet.
            row["Option"] = "Depth"
            rows.append(row)
        return rows, result["summary"], ""

    spec["propose"] = propose
    return spec


class CirclesEditor(CategoryEditor):
    label = "Circles"
    # Mirror the 'circles' worksheet: Xo, Yo, Option, Depth, Xi, Yi, R.
    FIELDS = [
        Field("Xo", "Xo"), Field("Yo", "Yo"),
        Field("Option", "Option", "choice", choices=["Depth", "Radius", "Intercept"]),
        Field("Depth", "Depth"), Field("Xi", "Xi"), Field("Yi", "Yi"), Field("R", "R"),
    ]

    def build(self, slope_data, parent):
        style = _doc_style(parent)
        window = _SearchWindowGroup(slope_data.get("search_window"),
                                    _unit_labels_for(slope_data))

        def preview(ax, rows, selected):
            _draw_circles_preview(ax, rows, selected, slope_data, style,
                                  window=window.result_values())

        return TableEditorDialog(
            "Circles", self.FIELDS, slope_data.get("circles", []), _new_circle, parent,
            help_text="Option sets how each circle's size is defined (only the matching "
                      "field is used):\n"
                      "  • Depth — elevation of the circle bottom at the center (R = Yo − Depth)\n"
                      "  • Radius — the circle radius R directly (Depth = Yo − R)\n"
                      "  • Intercept — a point (Xi, Yi) the circle passes through",
            preview_draw=preview,
            preview_caption="Preview shows the starting circles on the section "
                            "(selected circle bold with center, radius and depth "
                            "lines; others faint). Click a circle to select it. Any "
                            "search window is drawn with it: entry and exit ranges as "
                            "bars on the ground surface, the center box dashed.",
            pick_resolve=lambda x, y, tol, rows: _pick_circles(rows, x, y, tol, slope_data),
            field_help=CIRCLES_HELP,
            generate=_circles_generate_spec(slope_data),
            # Seven columns and a section drawing do not both fit across a dialog:
            # side by side, the table loses R off its right edge. Stacked, the table
            # gets the full width for its columns and the section — wide and short —
            # gets the room below.
            preview_below=True,
            # The search window limits where a SEARCH may run, so it belongs to the
            # circle table as a whole rather than to any one starting circle.
            extra_widget=window)

    def apply(self, slope_data, dlg):
        rows = dlg.result_rows()
        for c in rows:
            opt = str(c.get("Option", "Depth"))
            if opt == "Radius":
                c["Depth"] = c["Yo"] - c["R"]
            elif opt == "Intercept":
                c["R"] = ((c["Xi"] - c["Xo"]) ** 2 + (c["Yi"] - c["Yo"]) ** 2) ** 0.5
                c["Depth"] = c["Yo"] - c["R"]
            else:  # Depth
                c["R"] = c["Yo"] - c["Depth"]
        slope_data["circles"] = rows
        slope_data["circular"] = len(rows) > 0
        # Only the filled limits are stored, and a window nobody filled in leaves no
        # key behind at all — the state the loader produces for a blank J8:K17 block,
        # so an untouched model round-trips through this editor unchanged.
        window = dlg.extra.result_values()
        if window:
            slope_data["search_window"] = window
        else:
            slope_data.pop("search_window", None)


def _new_ncpt():
    return {"X": 0.0, "Y": 0.0, "Movement": "Free"}


NONCIRC_HELP = {
    "X": "X-coordinate of a point on the failure surface, listed left→right.",
    "Y": "Y-coordinate of a point on the failure surface.",
    "Movement": "Movement constraint used by the automated search — Free "
               "(unrestricted), Horiz (horizontal only), or Fixed (does not move).",
}


class WeakZoneDialog(QDialog):
    """Choose the material zone a generated non-circular surface should track.

    One dialog for two situations, which is the point of it. When no zone is
    clearly the weakest the generator has nothing to pre-select and this is how the
    question gets asked; when one *is*, this is how the user overrides it. The list
    is identical either way — every zone with its material colour, its computed
    mobilisable strength and its extent — so there is one thing to learn and one
    thing to build. Only the pre-selection and the wording above the list differ.
    """

    def __init__(self, zones, parent=None, headline="", preselect=0):
        super().__init__(parent)
        self.setWindowTitle("Choose the weak zone")
        self.resize(680, 340)
        layout = QVBoxLayout(self)
        layout.addWidget(_help_label(
            headline or "Choose the material zone the starting surface should track. "
                        "Strength is the shear strength each zone can mobilise at the "
                        "stress it actually carries, which is what makes zones of "
                        "different material models comparable."))

        self.table = QTableWidget(len(zones), 4, self)
        self.table.setHorizontalHeaderLabels(
            ["Zone", "Strength model", "Mobilisable strength", "Extent"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        head = self.table.horizontalHeader()
        head.setSectionResizeMode(0, QHeaderView.Stretch)
        for col in (1, 2, 3):
            head.setSectionResizeMode(col, QHeaderView.ResizeToContents)

        self._indices = []
        for row, zone in enumerate(zones):
            self._indices.append(zone.index)
            name = QTableWidgetItem(zone.name)
            name.setIcon(_material_swatch(zone.mat_id))
            self.table.setItem(row, 0, name)
            self.table.setItem(row, 1, QTableWidgetItem(str(zone.option or "mc")))
            self.table.setItem(row, 2, QTableWidgetItem(f"{zone.tau:.4g}"))
            self.table.setItem(row, 3, QTableWidgetItem(
                f"x {zone.x_range[0]:g} to {zone.x_range[1]:g}, "
                f"y {zone.y_range[0]:g} to {zone.y_range[1]:g}"))
        if zones:
            self.table.selectRow(max(0, min(preselect, len(zones) - 1)))
        layout.addWidget(self.table, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.table.doubleClicked.connect(self.accept)

    def chosen_zone(self):
        """The polygon index of the chosen zone."""
        row = self.table.currentRow()
        return self._indices[row if row >= 0 else 0]


def _noncirc_generate_spec(slope_data):
    """The Generate button's state and behaviour for the non-circular editor.

    Computed before anything is drawn, so a model the generator cannot serve
    produces a dimmed button carrying the reason rather than a live one that
    fails when pressed."""
    from xslope.generators import generate_noncircular_surface, rank_weak_zones

    zones = rank_weak_zones(slope_data)
    spec = {"label": "Generate from the weak zone…",
            "tooltip": "Build a starting surface that tracks the base of the "
                       "weakest material zone and ramps up to the ground at each "
                       "end — the shape a circular search cannot find."}
    if not zones:
        spec["available"] = False
        spec["reason"] = ("This model defines no material zones, so there is no "
                          "weak layer for a surface to track. Add the geometry on "
                          "the profile or polygon sheet first.")
        return spec
    if len(zones) == 1:
        spec["available"] = False
        spec["reason"] = (
            f"This model has one material zone ('{zones[0].name}'), so it has no "
            f"weak layer for a surface to track. A non-circular surface is for a "
            f"slope whose mechanism follows a weak zone; with one material the "
            f"mechanism follows the slope's own geometry, which is what a circular "
            f"search finds.")
        return spec
    spec["available"] = True

    def propose(parent):
        result = generate_noncircular_surface(slope_data, report=True)
        if result["surface"]:
            return result["surface"], result["summary"], ""
        if not result["candidates"] or len(result["candidates"]) < 2:
            return [], "", result["reason"]
        dlg = WeakZoneDialog(result["candidates"], parent,
                             headline=result["reason"] + ".")
        if not dlg.exec():
            return [], "", ""
        chosen = generate_noncircular_surface(slope_data, zone=dlg.chosen_zone(),
                                              report=True)
        if not chosen["surface"]:
            return [], "", chosen["reason"]
        return chosen["surface"], chosen["summary"], ""

    spec["propose"] = propose
    return spec


class NonCircEditor(CategoryEditor):
    label = "Non-circular surface"
    FIELDS = [Field("X", "X"), Field("Y", "Y"),
              Field("Movement", "Movement", "choice", choices=["Free", "Horiz", "Fixed"])]

    def build(self, slope_data, parent):
        style = _doc_style(parent)

        def preview(ax, rows, selected):
            _draw_noncirc_preview(ax, rows, selected, slope_data, style)

        return TableEditorDialog(
            "Non-circular surface", self.FIELDS, slope_data.get("non_circ", []), _new_ncpt, parent,
            help_text="Points ordered left→right. Entry/exit points use Movement='Free'.",
            preview_draw=preview,
            preview_caption="Preview shows the non-circular surface on the section "
                            "(selected vertex enlarged; ○ Free, ◇ Horiz-only, □ Fixed). "
                            "Click a vertex or the surface to select it.",
            pick_resolve=lambda x, y, tol, rows: _pick_noncirc(rows, x, y, tol),
            field_help=NONCIRC_HELP,
            generate=_noncirc_generate_spec(slope_data))

    def apply(self, slope_data, dlg):
        slope_data["non_circ"] = dlg.result_rows()


def _new_pt():
    return {"x": 0.0, "y": 0.0}


PIEZO_HELP = {
    "x": "X-coordinate of a point on this piezometric line, listed left→right.",
    "y": "Y-coordinate (elevation) of a point on this piezometric line.",
}


class PiezoEditor(CategoryEditor):
    label = "Piezometric lines"
    FIELDS = [Field("x", "x"), Field("y", "y")]

    def _rows(self, slope_data, key):
        return [{"x": x, "y": y} for (x, y) in (slope_data.get(key) or [])]

    def build(self, slope_data, parent):
        style = _doc_style(parent)

        def preview(ax, rows_per_tab, active_tab):
            _draw_piezo_preview(ax, rows_per_tab, active_tab, slope_data, style)

        return TabbedTableEditorDialog(
            "Piezometric lines",
            [("Line 1", self.FIELDS, self._rows(slope_data, "piezo_line"), _new_pt),
             ("Line 2 (rapid drawdown)", self.FIELDS, self._rows(slope_data, "piezo_line2"), _new_pt)],
            parent,
            help_text="Points ordered left→right. Line 2 is only used for rapid-drawdown "
                      "analysis (the drawdown / second water table).",
            preview_draw=preview,
            preview_caption="Preview shows both piezometric lines on the section (the "
                            "active tab's line bold with its points; the other dimmed). "
                            "Click a line to switch to its tab; click a vertex to select it.",
            pick_resolve=lambda x, y, tol, rpt: _pick_piezo(rpt, x, y, tol),
            field_help=PIEZO_HELP)

    def apply(self, slope_data, dlg):
        slope_data["piezo_line"] = [(r["x"], r["y"]) for r in dlg.result_rows(0)]
        slope_data["piezo_line2"] = [(r["x"], r["y"]) for r in dlg.result_rows(1)]


# --- distributed loads (two sets; each a list of point blocks) -------------- #
DLOAD_FIELDS = [Field("X", "X"), Field("Y", "Y"),
                Field("Normal", "Normal", unit="stress")]


def _new_dload_pt():
    return {"X": 0.0, "Y": 0.0, "Normal": 0.0}


def _apply_set_selection(widgets, tabs, select):
    """Activate the tab and list row encoded in ``select = (set_idx, row)`` — shared
    by the dload and seep-BC editors (both two-set tabbed master/detail dialogs), so
    a canvas double-click jumps straight to the picked item. No-op for other
    `select` shapes (e.g. None, or a plain index)."""
    if not (isinstance(select, tuple) and len(select) == 2):
        return
    s, row = select
    if s not in (0, 1):
        return
    tabs.setCurrentIndex(s)
    w = widgets[s]
    if row is not None and 0 <= row < w.list.count():
        w.list.setCurrentRow(row)


DLOADS_HELP = {
    "X": "X-coordinate of a point on the load's distribution line, listed left→right.",
    "Y": "Y-coordinate of a point on the load's distribution line.",
    "Normal": "Stress (force per unit area) at this point, applied in the load's "
              "Direction.",
    "dir": ("Which way this load pushes. Normal (blank in the file, and every load "
            "before template v21) applies it perpendicular to its own line — a "
            "pressure. Vertical applies it straight down whatever the line's slope, "
            "which is what a dead-weight surcharge does; on an inclined crest the two "
            "differ by a horizontal thrust the surcharge never has."),
}

# v21 distributed-load Direction: the template's two words, the wording the run
# dialogs use, and 'normal' first so it is the default a new load takes.
DLOAD_DIR_SPEC = {
    "label": "Direction:",
    "items": [("Normal (perpendicular to the line)", "normal"),
              ("Vertical (gravity surcharge)", "vertical")],
    "help_key": "dir",
}


class DloadsEditor(CategoryEditor):
    label = "Distributed loads"

    def build(self, slope_data, parent, select=None):
        dlg = QDialog(parent)
        dlg.setWindowTitle("Distributed loads")
        dlg.resize(640, 520)
        dlg._preview = None            # so the on_change closure is safe pre-build
        style = _doc_style(parent)
        layout = QVBoxLayout(dlg)
        layout.addWidget(_help_label(
            "Each load is a left→right series of points (≥2) plus a Direction — normal "
            "(a pressure perpendicular to its own line) or vertical (a gravity "
            "surcharge). Select a load to edit it. Set 2 is the second rapid-drawdown "
            "stage."))

        def schedule():
            if dlg._preview is not None:
                dlg._preview.schedule()

        tabs = QTabWidget()
        _ul = _unit_labels_for(slope_data)
        w1 = _BlockListWidget(slope_data.get("dloads"), DLOAD_FIELDS, _new_dload_pt,
                              "Load", on_change=schedule, unit_labels=_ul,
                              prop_spec=DLOAD_DIR_SPEC,
                              values=slope_data.get("dload_dirs"))
        w2 = _BlockListWidget(slope_data.get("dloads2"), DLOAD_FIELDS, _new_dload_pt,
                              "Load", on_change=schedule, unit_labels=_ul,
                              prop_spec=DLOAD_DIR_SPEC,
                              values=slope_data.get("dload2_dirs"))
        tabs.addTab(w1, "Set 1")
        tabs.addTab(w2, "Set 2 (rapid drawdown)")
        tabs.currentChanged.connect(lambda *_: schedule())

        def preview(ax):
            active = tabs.currentIndex()
            w = (w1, w2)[active] if active in (0, 1) else w1
            _draw_dloads_preview(ax, w1.pending_blocks(), w2.pending_blocks(),
                                 active, w.selected_block(), slope_data, style,
                                 set_dirs=(w1.pending_values(), w2.pending_values()))

        from .canvas import PreviewPane
        dlg._preview = PreviewPane(
            preview,
            caption="Preview shows both load sets on the section (set 1 / set 2 in "
                    "their own colors; the active tab's selected load emphasized, the "
                    "rest dimmed; a vertical load's arrows stand straight up). Click a "
                    "load block to select it (switching set if it belongs to the "
                    "other tab).")

        def on_pick(x, y, tol):
            hit = _pick_dloads((w1.pending_blocks(), w2.pending_blocks()), x, y, tol)
            if hit is None:
                return
            si, bi = hit
            tabs.setCurrentIndex(si)                 # -> the block's set becomes active
            w = (w1, w2)[si]
            if 0 <= bi < w.list.count():
                w.list.setCurrentRow(bi)

        dlg._preview.clicked.connect(on_pick)
        split = QSplitter(Qt.Horizontal)
        split.addWidget(tabs)
        split.addWidget(dlg._preview)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 1)
        split.setSizes([560, 500])
        layout.addWidget(split, 1)
        _ok_cancel(dlg, layout)

        def _help_key_for(widget):
            active = tabs.currentIndex()
            w = (w1, w2)[active] if active in (0, 1) else w1
            return w.help_key_for_widget(widget)

        strip = attach_help(dlg, DLOADS_HELP, _help_key_for)
        for w in (w1, w2):
            w._help_strip = strip
            w._field_help = DLOADS_HELP
            w._wire_help()          # the block table built during __init__ needs it too

        dlg._sets = (w1, w2)
        _apply_set_selection((w1, w2), tabs, select)
        dlg.resize(1160, 560)
        dlg._preview.refresh_now()
        return dlg

    def apply(self, slope_data, dlg):
        # Blocks and directions are written together, from the same records and in
        # the same order. Writing the blocks alone (which is what this did before the
        # Direction column existed) left dload_dirs at its pre-edit length and
        # positions, so deleting or reordering a load silently moved every later
        # direction onto the wrong load.
        w1, w2 = dlg._sets
        slope_data["dloads"] = w1.result_blocks()
        slope_data["dload_dirs"] = w1.result_values()
        slope_data["dloads2"] = w2.result_blocks()
        slope_data["dload2_dirs"] = w2.result_values()


# --- seep BC (two sets; each: a list of specified-head BCs + an exit face) --- #
XY_FIELDS = [Field("x", "x"), Field("y", "y")]


class _SeepBcSetWidget(QWidget):
    """One seepage BC set as a single master/detail: the left list holds each
    specified-head boundary, each specified-flux boundary, AND the exit face; the
    middle column shows one full-height point table plus a value field — "Head
    value:" for heads, "Flux value:" for fluxes, both hidden for the exit face.

    List order is heads … fluxes … exit face; picking.py mirrors this so a canvas
    double-click jumps to the right row.

    With ``slope_data`` the right side is a live section preview (the shared
    ``PreviewPane``) with the selected boundary emphasized; a click on a boundary in
    it selects the corresponding list row. The value/points column keeps a fixed
    narrow width — the points table is only x/y — and the preview absorbs the rest."""

    def __init__(self, bc, parent=None, unit_labels=None, constant_only=False,
                 slope_data=None, style=None, set_no=1):
        super().__init__(parent)
        bc = bc or {}
        self._unit_labels = unit_labels
        # Set 2 is the constant-steady rapid-drawdown boundary set: it can carry
        # neither a reservoir (submerged-only, time-varying) type nor a time-series
        # value (fileio rejects both at load time — there is one transient timeline and
        # it belongs to set 1). When True the Type selector is hidden (every head is a
        # plain Dirichlet) and a non-numeric value entered here is coerced to 0.
        self._constant_only = bool(constant_only)
        self._heads = [{"head": h.get("head", 0.0),
                        "kind": str(h.get("kind", "head")).strip().lower() or "head",
                        "coords": [tuple(c) for c in h.get("coords", [])]}
                       for h in (bc.get("specified_heads") or [])]
        self._fluxes = [{"flux": f.get("flux", 0.0),
                         "coords": [tuple(c) for c in f.get("coords", [])]}
                        for f in (bc.get("specified_fluxes") or [])]
        self._exit = [tuple(c) for c in (bc.get("exit_face") or [])]
        self._cur = -1
        self.table = None
        # Context-sensitive help: set by the owning dialog AFTER construction (it
        # needs the dialog's layout/button box built first — see attach_help), then
        # applied retroactively via _wire_help() to whichever table exists then and
        # to every one built afterward (_load rebuilds the table per selection).
        self._help_strip = None
        self._field_help = None

        left_pane = QWidget()
        left = QVBoxLayout(left_pane)
        left.setContentsMargins(0, 0, 0, 0)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._on_select)
        left.addWidget(self.list)
        lb = QHBoxLayout()
        b_add = QPushButton("Add head")
        b_add.clicked.connect(self._add_head)
        b_rem = QPushButton("Remove head")
        b_rem.clicked.connect(self._remove_head)
        lb.addWidget(b_add)
        lb.addWidget(b_rem)
        left.addLayout(lb)
        lb2 = QHBoxLayout()
        b_addf = QPushButton("Add flux")
        b_addf.clicked.connect(self._add_flux)
        b_remf = QPushButton("Remove flux")
        b_remf.clicked.connect(self._remove_flux)
        lb2.addWidget(b_addf)
        lb2.addWidget(b_remf)
        left.addLayout(lb2)

        # Head is a length; specified flux is a Darcy velocity (length/time), i.e. the
        # same dimension as k. Append the declared unit to each value label; undeclared
        # leaves the exact "Head value:" / "Flux value:" strings unchanged.
        head_text = "Head value:"
        flux_text = "Flux value:"
        if unit_labels:
            if unit_labels.get("length"):
                head_text = f"Head value ({unit_labels['length']}):"
            if unit_labels.get("k"):
                flux_text = f"Flux value ({unit_labels['k']}):"

        mid_pane = QWidget()
        right = QVBoxLayout(mid_pane)
        right.setContentsMargins(0, 0, 0, 0)
        # Head TYPE selector (plain Dirichlet head vs submerged-only reservoir).
        # Shown only for head rows; fluxes and the exit face hide it.
        trow = QHBoxLayout()
        self.type_label = QLabel("Type:")
        trow.addWidget(self.type_label)
        self.type_combo = QComboBox()
        self.type_combo.addItem("head")
        self.type_combo.addItem("reservoir")
        self.type_combo.setItemData(
            0, "held at the value/series at all nodes, all times", Qt.ToolTipRole)
        self.type_combo.setItemData(
            1, "draw the full face; nodes above the water level become seepage "
               "faces as it falls", Qt.ToolTipRole)
        trow.addWidget(self.type_combo, 1)
        right.addLayout(trow)
        # Set 2 has only plain heads — hide the whole Type row so reservoir can't be
        # picked (mirrors the fileio rule; keeps the two-solve rapid-drawdown set clean).
        if self._constant_only:
            self.type_label.setVisible(False)
            self.type_combo.setVisible(False)
        hrow = QHBoxLayout()
        self.head_label = QLabel(head_text)
        hrow.addWidget(self.head_label)
        self.head_edit = QLineEdit()
        self.head_edit.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        if self._constant_only:
            self.head_edit.setToolTip(
                "Set 2 is the constant, steady rapid-drawdown boundary set — enter a "
                "number. Time-varying (tseep series) values and reservoir boundaries "
                "belong on Set 1.")
        hrow.addWidget(self.head_edit, 1)
        right.addLayout(hrow)
        frow = QHBoxLayout()
        self.flux_label = QLabel(flux_text)
        frow.addWidget(self.flux_label)
        self.flux_edit = QLineEdit()
        self.flux_edit.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        frow.addWidget(self.flux_edit, 1)
        right.addLayout(frow)
        self._holder = QVBoxLayout()
        right.addLayout(self._holder, 1)

        body = QHBoxLayout(self)
        body.setContentsMargins(0, 0, 0, 0)
        self._preview = None
        if slope_data is not None:
            from .canvas import PreviewPane
            self._preview = PreviewPane(
                lambda ax: _draw_seep_bc_preview(
                    ax, self._pending_entries(), self.list.currentRow(),
                    slope_data, style, set_no),
                caption="Preview shows this set's boundaries on the section "
                        "(selected boundary bold; others dimmed). Click a boundary "
                        "to select it.")
            self._preview.clicked.connect(self._on_preview_click)
            # A type change recolors the selected boundary's family in the preview.
            self.type_combo.currentIndexChanged.connect(self._schedule_preview)
            splitter = QSplitter(Qt.Horizontal)
            splitter.addWidget(left_pane)
            splitter.addWidget(mid_pane)
            splitter.addWidget(self._preview)
            # The value/points column is just a value field over an x/y table, so it
            # keeps a fixed narrow share; dragging the dialog wider grows the preview.
            splitter.setStretchFactor(0, 0)
            splitter.setStretchFactor(1, 0)
            splitter.setStretchFactor(2, 1)
            splitter.setSizes([190, 260, 560])
            body.addWidget(splitter)
        else:
            body.addWidget(left_pane)
            body.addWidget(mid_pane, 1)

        self._refresh()
        self.list.setCurrentRow(0)
        if self._preview is not None:
            self._preview.refresh_now()

    # Index scheme: heads occupy [0, n_heads); fluxes [n_heads, n_heads+n_flux);
    # the exit face is the single trailing row.
    def _is_head(self, idx):
        return 0 <= idx < len(self._heads)

    def _is_flux(self, idx):
        return len(self._heads) <= idx < len(self._heads) + len(self._fluxes)

    def _flux_idx(self, idx):
        return idx - len(self._heads)

    def _is_exit(self, idx):
        return idx == len(self._heads) + len(self._fluxes)

    def _refresh(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i, h in enumerate(self._heads):
            self.list.addItem(f"Head {i + 1}  (h = {h['head']})")
        for i, f in enumerate(self._fluxes):
            self.list.addItem(f"Flux {i + 1}  (q = {f['flux']})")
        self.list.addItem("Exit face")
        self.list.blockSignals(False)

    def _commit(self):
        if self._cur < 0:
            return
        coords = [(r["x"], r["y"]) for r in self.table.result_rows()] if self.table else []
        if self._is_exit(self._cur):
            self._exit = coords
        elif self._is_flux(self._cur):
            f = self._fluxes[self._flux_idx(self._cur)]
            txt = self.flux_edit.text()
            if not _unedited(txt, f.get("flux")):   # untouched keeps the stored value
                try:
                    f["flux"] = float(txt or 0)
                except ValueError:
                    f["flux"] = 0.0
            f["coords"] = coords
        elif self._is_head(self._cur):
            txt = self.head_edit.text().strip()
            try:
                if not _unedited(txt, self._heads[self._cur].get("head")):
                    self._heads[self._cur]["head"] = float(txt or 0)
            except ValueError:
                # a non-numeric head value is a tseep series name (time-varying BC) —
                # allowed on set 1 only; set 2 is constant, so coerce it to 0.
                self._heads[self._cur]["head"] = 0.0 if self._constant_only else txt
            # Set 2 is plain-head only; its Type selector is hidden, so force "head".
            self._heads[self._cur]["kind"] = (
                "head" if self._constant_only else self.type_combo.currentText())
            self._heads[self._cur]["coords"] = coords

    def _load(self, idx):
        if self.table is not None:
            self.table.setParent(None)
            self.table = None
        if idx < 0:
            return
        is_head = self._is_head(idx)
        is_flux = self._is_flux(idx)
        # Set 2 never shows the Type selector (constant plain heads only).
        self.type_label.setVisible(is_head and not self._constant_only)
        self.type_combo.setVisible(is_head and not self._constant_only)
        self.head_label.setVisible(is_head)
        self.head_edit.setVisible(is_head)
        self.flux_label.setVisible(is_flux)
        self.flux_edit.setVisible(is_flux)
        if is_head:
            kind = self._heads[idx].get("kind", "head")
            self.type_combo.setCurrentText("reservoir" if kind == "reservoir" else "head")
            self.head_edit.setText(_display_number(self._heads[idx]["head"]))
            rows = [{"x": x, "y": y} for (x, y) in self._heads[idx]["coords"]]
        elif is_flux:
            f = self._fluxes[self._flux_idx(idx)]
            self.flux_edit.setText(_display_number(f["flux"]))
            rows = [{"x": x, "y": y} for (x, y) in f["coords"]]
        else:  # exit face
            rows = [{"x": x, "y": y} for (x, y) in self._exit]
        self.table = _EditableTable(XY_FIELDS, rows, _new_pt,
                                    on_change=self._schedule_preview)
        self._holder.addWidget(self.table)
        self._wire_help()

    def help_key_for_widget(self, widget):
        """The field key ``widget`` edits: 'head'/'flux' for the value edits, else
        the point table's column key (x/y), or None."""
        if widget is self.head_edit:
            return "head"
        if widget is self.flux_edit:
            return "flux"
        if self.table is not None:
            col = self.table.column_at(widget)
            if 0 <= col < len(XY_FIELDS):
                return XY_FIELDS[col].key
        return None

    def _wire_help(self):
        _wire_cell_help(self._help_strip, self.table, XY_FIELDS, self._field_help or {})

    def _on_select(self, idx):
        self._commit()
        self._cur = idx
        self._load(idx)
        self._schedule_preview()

    # --- live preview -----------------------------------------------------
    def _schedule_preview(self, *_):
        if self._preview is not None:
            self._preview.schedule()

    def _pending_entries(self):
        """The set's boundaries in list order (heads … fluxes … exit face) with the
        current selection's pending edits committed — what the preview draws."""
        self._commit()
        return ([{"kind": h.get("kind", "head"), "coords": h["coords"]}
                 for h in self._heads]
                + [{"kind": "flux", "coords": f["coords"]} for f in self._fluxes]
                + [{"kind": "exit", "coords": self._exit}])

    def _on_preview_click(self, x, y, tol):
        """Select the boundary nearest the preview click (within tolerance)."""
        from shapely.geometry import LineString, Point
        pt = Point(x, y)
        best, best_d = None, float("inf")
        for i, e in enumerate(self._pending_entries()):
            coords = [c for c in (e["coords"] or []) if c is not None]
            try:
                d = (LineString(coords).distance(pt) if len(coords) >= 2
                     else Point(coords[0]).distance(pt) if coords
                     else float("inf"))
            except Exception:
                continue
            if d < best_d:
                best, best_d = i, d
        if best is not None and best_d <= tol:
            self.list.setCurrentRow(best)

    def _add_head(self):
        self._commit()
        self._heads.append({"head": 0.0, "kind": "head", "coords": []})
        self._cur = -1                     # rows below shifted; avoid a stale commit
        self._refresh()
        self.list.setCurrentRow(len(self._heads) - 1)

    def _remove_head(self):
        idx = self.list.currentRow()
        if not self._is_head(idx):         # only a head row is removable here
            return
        self._heads.pop(idx)
        self._cur = -1
        self._refresh()
        self.list.setCurrentRow(min(idx, self.list.count() - 1))

    def _add_flux(self):
        self._commit()
        self._fluxes.append({"flux": 0.0, "coords": []})
        self._cur = -1
        self._refresh()
        self.list.setCurrentRow(len(self._heads) + len(self._fluxes) - 1)

    def _remove_flux(self):
        idx = self.list.currentRow()
        if not self._is_flux(idx):
            return
        self._fluxes.pop(self._flux_idx(idx))
        self._cur = -1
        self._refresh()
        self.list.setCurrentRow(min(idx, self.list.count() - 1))

    def result(self):
        self._commit()
        return {"specified_heads": [{"head": h["head"],
                                     "kind": h.get("kind", "head"),
                                     "coords": list(h["coords"])}
                                    for h in self._heads],
                "specified_fluxes": [{"flux": f["flux"], "coords": list(f["coords"])}
                                     for f in self._fluxes],
                "exit_face": list(self._exit)}


SEEPBC_HELP = {
    "head": "Total head at this boundary — the height of water above the datum "
           "(length units).",
    "flux": "Specified flux — the normal Darcy velocity across the boundary "
           "(length/time); positive = inflow (infiltration/recharge).",
    "x": "X-coordinate of a point on this boundary (or the exit face).",
    "y": "Y-coordinate of a point on this boundary (or the exit face).",
}


class SeepBcEditor(CategoryEditor):
    label = "Seep BC"

    def build(self, slope_data, parent, select=None):
        dlg = QDialog(parent)
        dlg.setWindowTitle("Seep BC")
        dlg.resize(1080, 560)
        layout = QVBoxLayout(dlg)
        layout.addWidget(_help_label(
            "Seepage boundary conditions: specified-head boundaries (each a head value + "
            "its points), specified-flux boundaries (each a flux value + its points), and "
            "an exit face (where water leaves the slope). Set 2 is the constant, steady "
            "rapid-drawdown set (the second seepage solution): it takes plain heads only "
            "— reservoir and time-varying (tseep series) boundaries belong on Set 1."))
        tabs = QTabWidget()
        _ul = _unit_labels_for(slope_data)
        _style = _doc_style(parent)
        w1 = _SeepBcSetWidget(slope_data.get("seepage_bc"), unit_labels=_ul,
                              slope_data=slope_data, style=_style, set_no=1)
        w2 = _SeepBcSetWidget(slope_data.get("seepage_bc2"), unit_labels=_ul,
                              constant_only=True,
                              slope_data=slope_data, style=_style, set_no=2)
        tabs.addTab(w1, "Set 1")
        tabs.addTab(w2, "Set 2 (rapid drawdown)")
        # A tab shown for the first time renders its preview at the real viewport
        # size (the hidden tab's first paint happened at zero size).
        tabs.currentChanged.connect(
            lambda i: (w1, w2)[i]._schedule_preview() if i in (0, 1) else None)
        layout.addWidget(tabs)
        _ok_cancel(dlg, layout)

        def _help_key_for(widget):
            active = tabs.currentIndex()
            w = (w1, w2)[active] if active in (0, 1) else w1
            return w.help_key_for_widget(widget)

        strip = attach_help(dlg, SEEPBC_HELP, _help_key_for)
        for w in (w1, w2):
            w._help_strip = strip
            w._field_help = SEEPBC_HELP
            w._wire_help()          # the point table built during __init__ needs it too

        dlg._sets = (w1, w2)
        _apply_set_selection((w1, w2), tabs, select)
        return dlg

    def apply(self, slope_data, dlg):
        w1, w2 = dlg._sets
        slope_data["seepage_bc"] = w1.result()
        slope_data["seepage_bc2"] = w2.result()
        bc2 = slope_data["seepage_bc2"]
        slope_data["has_seepage_bc2"] = bool(bc2.get("specified_heads") or bc2.get("exit_face"))


# --- piles ------------------------------------------------------------------ #
def _new_pile():
    return {"label": "Pile", "x1": 0.0, "y1": 0.0, "x2": 0.0, "y2": 0.0, "H": None,
            "theta_p": 0.0, "D_pile": None, "S": None, "E": None, "I": None,
            "area": None, "V_cap": None, "M_cap": None, "head_fixity": "free",
            "tip_fixity": "free", "appl": "active"}


# List-view form layout for a pile: every PilesEditor.FIELDS key grouped (Identity /
# Geometry / Capacity / Behavior). theta_p is NOT a field — it's derived from the axis
# on apply — so it has no widget and rides along as an unshown key.
_PILE_FORM_GROUPS = [
    ("Identity", [["label"]]),
    ("Geometry", [["x1", "y1"], ["x2", "y2"]]),
    ("Capacity / design", [["H"], ["D_pile", "S"], ["V_cap", "M_cap"],
                           ["E", "I"], ["area"]]),
    ("Behavior", [["appl"], ["head_fixity", "tip_fixity"]]),
]


def _pile_item_label(i, row):
    name = str(row.get("label") or "Pile")
    try:
        return f"{i + 1}. {name} @ x={float(row.get('x1', 0) or 0):g}"
    except (TypeError, ValueError):
        return f"{i + 1}. {name}"


# Wording mirrors the 'piles' worksheet section of input_template.md. Where a field
# is genuinely single-analysis (matching its column's usage tag / header color) the
# text says so explicitly; D and S feed both engines — the LEM force auto-computation
# AND the FEM section and its per-unit-width smear — so neither is tagged "LEM only"
# and each tooltip names both readers.
PILES_HELP = {
    "label": "Name used in error messages, summaries, and plots (optional).",
    "x1": "Pile top X-coordinate.",
    "y1": "Pile top Y-coordinate.",
    "x2": "Pile tip (bottom) X-coordinate.",
    "y2": "Pile tip (bottom) Y-coordinate.",
    "H": "Pile force per unit width of slope (force/length). Blank = auto-computed "
        "via Ito & Matsui from D and S (vertical piles only). LEM only; the FEM "
        "never reads a stated pile force.",
    "D_pile": "Pile diameter. Read by both engines: LEM needs it for the Ito & "
             "Matsui auto-computation of H, and FEM derives the section from it "
             "when I and Area are blank — I = pi·D^4/64 and Area = pi·D^2/4 for a "
             "solid circular pile.",
    "S": "Center-to-center pile spacing. In the LEM, spacing is physics: it sets "
        "the arching between piles (Ito & Matsui) and makes Vcap/Mcap per-pile. "
        "In the FEM, S converts per-pile properties to per-unit-width (EA/S, "
        "EI/S) — the correct plane-strain stiffness of the smeared row; what it "
        "cannot model is the soil arching between piles.",
    "E": "Young's modulus of the pile material; with I and Area gives the flexural "
        "(EI) and axial (EA) stiffness, each divided by S for the per-unit-width "
        "2D beam. FEM only.",
    "I": "Moment of inertia of a single pile (÷ S for per unit width). FEM only; "
        "auto-computed from D for a solid circular section when left blank.",
    "area": "Cross-sectional area of a single pile (÷ S for per unit width). FEM "
           "only; auto-computed from D for a solid circular section when left blank.",
    "V_cap": "Shear capacity of a single pile (force units); requires S. Read by "
            "both engines: the LEM checks the per-pile shear against it; the FEM "
            "clips the beam shear at Vcap ÷ S per unit width.",
    "M_cap": "Moment capacity of a single pile (force×length); requires S. Read by "
            "both engines: the LEM checks the per-pile moment against it; the FEM "
            "releases a plastic hinge where the beam moment reaches Mcap ÷ S per "
            "unit width.",
    "appl": "Force application — Active: H is an allowable force, not divided by "
           "FS (default). Passive: H is an ultimate capacity divided by FS. LEM only.",
    "head_fixity": "Restraint at the top of the pile: free (default; no connection "
                   "at the head), pinned (translation held, rotation free — tie-rods "
                   "or anchors), unrotated (rotation held, translation free — a cap "
                   "beam tying the heads), or fixed (both held — cap beam and "
                   "anchors). FEM only.",
    "tip_fixity": "Restraint at the bottom of the pile: free (default; the tip "
                  "moves with the soil around it, or sits on the model boundary), "
                  "pinned (translation held, rotation free — bearing on a hard "
                  "stratum inside the mesh), unrotated (rotation held, translation "
                  "free), or fixed (both held — socketed into rock). FEM only.",
}


class PilesEditor(CategoryEditor):
    label = "Piles"
    # tooltip=PILES_HELP[key] gives the table-header hover and the list-view
    # label/edit hover; the same dict feeds the context-sensitive help strip.
    # Field order mirrors the piles sheet's columns (Label, the endpoints, H, Appl,
    # D, S, Vcap, Mcap, then the FEM tail E, I, Area, Head, Tip) so a block copied from
    # the sheet or the docs' tables pastes straight in. The sheet's qp (θ) sits
    # between H and Appl and has no column here: θ is derived from the pile axis on
    # save, so a block spanning it goes in as two — the endpoints through H, then D
    # onward, which is how the tutorials print it.
    LF = {"lem", "fem"}
    FIELDS = [
        Field("label", "Label", "str", tooltip=PILES_HELP["label"]),
        Field("x1", "x1", tooltip=PILES_HELP["x1"]), Field("y1", "y1", tooltip=PILES_HELP["y1"]),
        Field("x2", "x2", tooltip=PILES_HELP["x2"]), Field("y2", "y2", tooltip=PILES_HELP["y2"]),
        Field("H", "H", "optfloat", usage="lem", tooltip=PILES_HELP["H"]),
        # Force application (v12, LEM only): active = allowable force applied as-is;
        # passive = ultimate capacity divided by FS (loader default 'active').
        Field("appl", "Appl", "choice", choices=["active", "passive"], usage="lem",
              tooltip=PILES_HELP["appl"]),
        # D and S are applies=LF, like the reinforcement editor's Spacing: the FEM
        # beam assembly derives I and Area from D when they are blank and divides
        # EA and EI by S, so neither usage toggle may hide them and neither header
        # is colored for a single engine.
        Field("D_pile", "D", "optfloat", applies=LF, tooltip=PILES_HELP["D_pile"]),
        Field("S", "S", "optfloat", applies=LF, tooltip=PILES_HELP["S"]),
        Field("V_cap", "Vcap", "optfloat", applies=LF, tooltip=PILES_HELP["V_cap"]),
        Field("M_cap", "Mcap", "optfloat", applies=LF, tooltip=PILES_HELP["M_cap"]),
        Field("E", "E", "optfloat", usage="fem", tooltip=PILES_HELP["E"]),
        Field("I", "I", "optfloat", usage="fem", tooltip=PILES_HELP["I"]),
        Field("area", "Area", "optfloat", usage="fem", tooltip=PILES_HELP["area"]),
        Field("head_fixity", "Head", "choice", choices=["free", "pinned", "unrotated", "fixed"], usage="fem",
              tooltip=PILES_HELP["head_fixity"]),
        Field("tip_fixity", "Tip", "choice", choices=["free", "pinned", "unrotated", "fixed"], usage="fem",
              tooltip=PILES_HELP["tip_fixity"]),
    ]

    def build(self, slope_data, parent):
        style = _doc_style(parent)

        def preview(ax, rows, selected):
            _draw_piles_preview(ax, rows, selected, slope_data, style)

        return _LineEditorDialog(
            "Piles", self.FIELDS, slope_data.get("pile_lines", []), _new_pile,
            _PILE_FORM_GROUPS, _pile_item_label, preview,
            lambda x, y, tol, rows: _pick_line_rows(rows, x, y, tol),
            view_state="piles", parent=parent, unit_labels=_unit_labels_for(slope_data),
            dynamic_spec={"driver": "S",
                          "fields": {"V_cap": "force", "M_cap": "moment",
                                     "area": "area"}},
            help_text="List view edits one pile at a time as a grouped form beside a "
                      "live section preview; the table view is available for bulk entry "
                      "of many piles. Both views edit the same rows, so switching is "
                      "lossless. Leave H blank for auto Ito & Matsui force. I / Area "
                      "auto-compute from D when blank. θ is auto-derived from the pile "
                      "axis. Vcap/Mcap require S (spacing). Appl: active = allowable "
                      "force; passive = ultimate capacity ÷ FS.",
            usage_toggles=["lem", "fem"],
            preview_caption="Preview shows the pile rows on the section (selected "
                            "row bold with □ cap and ▽ tip markers; others dimmed). "
                            "Click a pile row to select it.",
            field_help=PILES_HELP)

    def apply(self, slope_data, dlg):
        rows = dlg.result_rows()
        for p in rows:  # θ is auto-derived from the axis (loader does the same)
            dx, dy = p["x2"] - p["x1"], p["y2"] - p["y1"]
            p["theta_p"] = math.degrees(math.atan2(dx, -dy))
        slope_data["pile_lines"] = rows


# --- geometry: profile lines & polygons (master/detail) --------------------- #
def _set_derived_geometry(slope_data, polys):
    """Store material-zone polygons and rebuild the geometry derived from them
    (ground surface, domain polygon, tension-crack line) — the resync the loader
    performs after parsing. ``polys`` is the [{'polygon': Polygon, 'mat_id': int}]
    list both the profile and polygon editors converge on."""
    from shapely.geometry import LineString
    from xslope.fileio import build_ground_surface_from_polygons

    slope_data["polygons"] = polys
    gs, dom = (build_ground_surface_from_polygons(polys) if polys else (None, None))
    slope_data["ground_surface"] = gs
    slope_data["domain_polygon"] = dom
    td = slope_data.get("tcrack_depth", 0)
    if td and td > 0 and gs is not None and not gs.is_empty:
        slope_data["tcrack_surface"] = LineString([(x, y - td) for x, y in gs.coords])
    else:
        slope_data["tcrack_surface"] = None


def _parse_opt_size(text):
    """Parse an optional local-mesh-Size field. Returns ``(ok, value)``: a blank
    field is ``(True, None)`` — the global target size — a positive number is
    ``(True, float)``, and anything else is ``(False, None)``. Mirrors the loader's
    ``_opt_size_cell``, which rejects a non-numeric or non-positive Size rather than
    ignoring it, so the Studio and the file agree on what a valid Size is."""
    s = (text or "").strip()
    if not s:
        return True, None
    try:
        v = float(s)
    except (TypeError, ValueError):
        return False, None
    return (True, v) if v > 0 else (False, None)


def _normalize_polygon(p):
    """Coerce a polygon entry to the canonical {'polygon': Polygon, 'mat_id': int}.
    Accepts entries that carry a shapely 'polygon', a 'coords' list, or a bare
    shapely Polygon — so a model that built polygons with only 'coords' still
    resyncs."""
    from shapely.geometry import Polygon
    if isinstance(p, dict):
        poly = p.get("polygon")
        if poly is None and p.get("coords"):
            poly = Polygon(p["coords"])
        return {"polygon": poly, "mat_id": p.get("mat_id"), "size": p.get("size")}
    return {"polygon": p, "mat_id": None, "size": None}


def _resync_geometry(slope_data):
    """Rebuild polygons (and the geometry derived from them: ground surface,
    domain polygon, t-crack line) from the edited source — the resync the loader/
    design driver also do. Two source paths: profile_lines (build zone polygons
    from them), or polygons set directly (e.g. a zoned dam built without
    profile_lines), which are normalized and used as-is."""
    from shapely.geometry import Polygon
    from xslope.mesh import build_polygons

    profile_lines = slope_data.get("profile_lines") or []
    if profile_lines:
        polys = [{"polygon": Polygon(p["coords"]), "mat_id": p["mat_id"],
                  "size": p.get("size")}
                 for p in build_polygons(slope_data={"profile_lines": profile_lines,
                                                      "max_depth": slope_data.get("max_depth")})]
        _set_derived_geometry(slope_data, polys)
        return
    raw = slope_data.get("polygons") or []
    if raw:
        _set_derived_geometry(slope_data, [_normalize_polygon(p) for p in raw])


class MatGeometryDialog(QDialog):
    """Master/detail editor over a list of material-tagged coordinate sequences:
    the items list (left) + the selected item's material and vertex table (right).
    Shared by the profile-line editor (open polylines) and the polygon editor
    (closed rings); the two differ only in labels, help text, and how coords are
    extracted/rebuilt by the owning CategoryEditor."""

    XY = [Field("x", "x"), Field("y", "y")]

    def __init__(self, title, help_text, item_label, items, materials, parent=None,
                 select=None, max_depth=None, preview_draw=None, preview_caption=None,
                 slope_data=None, style=None, pick_resolve=None, field_help=None,
                 poly_types=False):
        # items: list of {"mat_id": int|None, "coords": [(x, y), ...],
        #                 "kind": str, "size": float|None}
        # poly_types: show the v21 Type combo (polygon sheet only — a profile line is
        #   always a material boundary, so the profile sheet has no Type cell).
        # select: row to pre-highlight (e.g. the double-clicked line); else first.
        # max_depth: when not None, show a "Max depth" field (profile sheet only —
        #   it has no meaning for polygon input); the bottom boundary elevation
        #   used when building zone polygons from the profile lines.
        # preview_draw: draw(ax, lines, selected_index, max_depth) — a per-editor
        #   hook that renders the full section with the selected item highlighted.
        #   When given, a live preview pane is added on the right behind a splitter.
        super().__init__(parent)
        self.setWindowTitle(title)
        self._item_label = item_label
        self._materials = materials
        # One record per item: material, v21 Type kind, v21 local mesh Size, vertices.
        # Every field the record carries is edited here, so nothing has to survive as
        # a pass-through — but the copy is still explicit, because a rebuild-from-
        # fields apply() is exactly what silently deletes a key nobody listed.
        self._lines = [{"mat_id": it.get("mat_id"), "size": it.get("size"),
                        "kind": it.get("kind") or "material",
                        "coords": [tuple(c) for c in it.get("coords", [])]}
                       for it in (items or [])]
        self._cur = -1
        self.table = None
        self._preview_draw = preview_draw
        self._preview = None
        # ``pick_resolve(x, y, tol, lines) -> (feature, row|None) | None`` maps a
        # preview click to the item (list row) and, for a vertex hit, its vertex row.
        self._pick_resolve = pick_resolve
        # Optional field-key -> help-text mapping ("x"/"y"/"mat_id"/"size", + "type"
        # on the polygon sheet and "max_depth" on the profile sheet) for the
        # context-sensitive help strip.
        self._field_help = field_help
        # Combo index -> mat_id, identity over the materials list. The mapping stays
        # explicit (rather than "index == mat_id") so the combo can never be read as
        # a material index it isn't.
        self._mat_ids = list(range(len(self._materials)))
        self._poly_types = bool(poly_types)

        main = QVBoxLayout(self)
        main.addWidget(_help_label(help_text))

        self._max_depth_edit = None
        if max_depth is not None:
            mdrow = QHBoxLayout()
            mdrow.addWidget(QLabel("Max depth (bottom boundary elevation):"))
            self._max_depth = max_depth
            self._max_depth_edit = QLineEdit(_display_number(max_depth))
            self._max_depth_edit.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self._max_depth_edit.setToolTip(
                "Elevation of the model's bottom boundary, used to build the zone "
                "polygons from the profile lines.")
            self._max_depth_edit.setMaximumWidth(90)
            self._max_depth_edit.textChanged.connect(self._schedule_preview)
            mdrow.addWidget(self._max_depth_edit)
            mdrow.addStretch(1)
            main.addLayout(mdrow)

        # Master/detail (list | material+vertex table), plus an optional preview pane
        # on the right behind a splitter — the materials list-view layout pattern.
        left_w = QWidget()
        left = QVBoxLayout(left_w)
        left.setContentsMargins(0, 0, 0, 0)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._on_select)
        left.addWidget(self.list)
        lbtns = QHBoxLayout()
        b_add = QPushButton(f"Add {item_label.lower()}")
        b_add.clicked.connect(self._add_line)
        b_rem = QPushButton(f"Remove {item_label.lower()}")
        b_rem.clicked.connect(self._remove_line)
        lbtns.addWidget(b_add)
        lbtns.addWidget(b_rem)
        left.addLayout(lbtns)

        right_w = QWidget()
        right = QVBoxLayout(right_w)
        right.setContentsMargins(0, 0, 0, 0)
        matrow = QHBoxLayout()

        # v21 Type (polygon sheet only). The words are the template's, and picking a
        # non-material one greys the Mat ID + material name exactly as the sheet's
        # row-7 echo formula does — an overlay has no material.
        self.type_combo = None
        if self._poly_types:
            matrow.addWidget(QLabel("Type:"))
            self.type_combo = QComboBox()
            for word, _kind in POLYGON_TYPE_ITEMS:
                self.type_combo.addItem(word)
            _type_help = (self._field_help or {}).get("type", "")
            if _type_help:
                self.type_combo.setToolTip(_type_help)
            self.type_combo.currentIndexChanged.connect(self._on_type_changed)
            matrow.addWidget(self.type_combo)

        self._mat_label = QLabel("Material:")
        matrow.addWidget(self._mat_label)
        self.mat_combo = QComboBox()
        for mid in self._mat_ids:
            m = materials[mid] if mid < len(materials) else {}
            self.mat_combo.addItem(f"{mid + 1}: {m.get('name', '')}")
        self.mat_combo.currentIndexChanged.connect(self._on_mat_changed)
        matrow.addWidget(self.mat_combo, 1)

        # v21 local mesh Size. Optional on every item; blank = the global target size.
        _size_help = (self._field_help or {}).get("size", "")
        matrow.addWidget(QLabel("Size:"))
        self.size_edit = QLineEdit()
        self.size_edit.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.size_edit.setPlaceholderText("global")
        # Fixed, not just capped: the material combo takes the row's stretch, and a
        # merely-capped Size edit collapses under it and scrolls its own value out of
        # sight — a declared refinement has to stay readable.
        self.size_edit.setFixedWidth(76)
        if _size_help:
            self.size_edit.setToolTip(_size_help)
        self.size_edit.textChanged.connect(self._on_size_changed)
        matrow.addWidget(self.size_edit)
        right.addLayout(matrow)
        self._holder = QVBoxLayout()
        right.addLayout(self._holder, 1)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_w)
        splitter.addWidget(right_w)
        if preview_draw is not None:
            from .canvas import PreviewPane
            self._preview = PreviewPane(
                lambda ax: self._preview_draw(ax, self._pending_lines(), self._cur,
                                              self.result_max_depth()),
                caption=preview_caption)
            if self._pick_resolve is not None:
                self._preview.clicked.connect(self._on_preview_click)
            splitter.addWidget(self._preview)
            splitter.setStretchFactor(0, 0)
            splitter.setStretchFactor(1, 0)
            splitter.setStretchFactor(2, 1)
            splitter.setSizes([190, 320, 470])
            self.resize(1060, 560)
        else:
            splitter.setStretchFactor(0, 0)
            splitter.setStretchFactor(1, 1)
            self.resize(680, 540)
        main.addWidget(splitter, 1)

        _ok_cancel(self, main)

        if self._field_help is not None:
            # Attached before the first _load (below), so the initial vertex table
            # picks up cell-tracking help immediately, exactly as the later ones do.
            attach_help(self, self._field_help, self._help_key_for)

        self._refresh_list()
        if self._lines:
            row = select if (select is not None and 0 <= select < len(self._lines)) else 0
            self.list.setCurrentRow(row)
        else:
            self._sync_type_enabled()      # empty list: still show a coherent form
        if self._preview is not None:
            self._preview.refresh_now()

    def _help_key_for(self, widget):
        """The field key ``widget`` edits: 'mat_id' for the material combo, 'type'
        for the Type combo (polygon only), 'size' for the Size edit, 'max_depth' for
        that edit (profile only), else the vertex table's column key (x/y), or None."""
        if widget is self.mat_combo:
            return "mat_id"
        if self.type_combo is not None and widget is self.type_combo:
            return "type"
        if widget is self.size_edit:
            return "size"
        if self._max_depth_edit is not None and widget is self._max_depth_edit:
            return "max_depth"
        if self.table is not None:
            col = self.table.column_at(widget)
            if 0 <= col < len(self.XY):
                return self.XY[col].key
        return None

    def _wire_help(self):
        _wire_cell_help(getattr(self, "_help_strip", None), self.table, self.XY,
                        self._field_help or {})

    def _schedule_preview(self, *_):
        if self._preview is not None:
            self._preview.schedule()

    def _on_preview_click(self, x, y, tol):
        """A click in the preview: resolve to (feature, row). Select that item in the
        list (switching feature loads its vertex table via _on_select); if the hit was
        on a vertex, also select that vertex row. Selection re-renders the preview."""
        hit = self._pick_resolve(x, y, tol, self._pending_lines())
        if hit is None:
            return
        feature, row = hit
        if not (0 <= feature < len(self._lines)):
            return
        if feature != self._cur:
            self.list.setCurrentRow(feature)     # -> _on_select rebuilds self.table
        if row is not None and self.table is not None:
            self.table.select_row(row)

    def _pending_lines(self):
        """The current items with the SELECTED one's coords/material taken live from
        the vertex table + material combo (they aren't committed to self._lines until
        a selection change / OK), so the preview reflects the in-progress edit."""
        lines = [{"mat_id": ln["mat_id"], "kind": ln.get("kind") or "material",
                  "coords": list(ln["coords"])}
                 for ln in self._lines]
        if 0 <= self._cur < len(lines) and self.table is not None:
            coords = [(r["x"], r["y"]) for r in self.table.result_rows()]
            mid = self._mat_id_at(self.mat_combo.currentIndex())
            lines[self._cur] = {"mat_id": mid if mid is not None else lines[self._cur]["mat_id"],
                                "kind": self._kind_at(),
                                "coords": coords}
        return lines

    def _kind_at(self):
        """The Type combo's current kind ('material' when there is no Type combo —
        a profile line is always a material boundary)."""
        if self.type_combo is None:
            return "material"
        i = self.type_combo.currentIndex()
        return POLYGON_TYPE_ITEMS[i][1] if 0 <= i < len(POLYGON_TYPE_ITEMS) else "material"

    def _type_index(self, kind):
        """Combo row for a kind, defaulting to 'material' for anything unrecognized."""
        for i, (_word, k) in enumerate(POLYGON_TYPE_ITEMS):
            if k == kind:
                return i
        return 0

    def _combo_index(self, mat_id):
        """Combo row for a mat_id, or 0 when it is unset / not offered."""
        try:
            return self._mat_ids.index(mat_id)
        except ValueError:
            return 0

    def _mat_id_at(self, index):
        """mat_id for a combo row, or None when the row is out of range."""
        if 0 <= index < len(self._mat_ids):
            return self._mat_ids[index]
        return None

    def _label(self, i):
        ln = self._lines[i]
        kind = ln.get("kind") or "material"
        size = ln.get("size")
        # A declared Size is shown in the list so a refinement is visible without
        # clicking through every item (and a refine region, which is NOTHING but its
        # size, always reads as something).
        tail = f", size {size:g}" if isinstance(size, (int, float)) else ""
        if kind != "material":
            return (f"{self._item_label} {i + 1}  "
                    f"({POLYGON_KIND_LABELS.get(kind, kind)}{tail})")
        mid = ln["mat_id"]
        if mid is not None and 0 <= mid < len(self._materials):
            return (f"{self._item_label} {i + 1}  (mat {mid + 1}: "
                    f"{self._materials[mid].get('name', '')}{tail})")
        return f"{self._item_label} {i + 1}  (mat ?{tail})"

    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i in range(len(self._lines)):
            self.list.addItem(self._label(i))
        self.list.blockSignals(False)

    def _commit_current(self):
        if 0 <= self._cur < len(self._lines) and self.table is not None:
            self._lines[self._cur]["coords"] = [(r["x"], r["y"]) for r in self.table.result_rows()]
            mid = self._mat_id_at(self.mat_combo.currentIndex())
            if mid is not None:
                self._lines[self._cur]["mat_id"] = mid
            self._lines[self._cur]["kind"] = self._kind_at()
            ok, size = _parse_opt_size(self.size_edit.text())
            if ok:
                self._lines[self._cur]["size"] = size

    def _load(self, idx):
        if self.table is not None:
            self.table.setParent(None)
            self.table = None
        if not (0 <= idx < len(self._lines)):
            return
        ln = self._lines[idx]
        rows = [{"x": x, "y": y} for (x, y) in ln["coords"]]
        self.table = _EditableTable(self.XY, rows, _new_pt,
                                    on_change=self._schedule_preview)
        self._holder.addWidget(self.table)
        self._wire_help()
        self.mat_combo.blockSignals(True)
        self.mat_combo.setCurrentIndex(self._combo_index(ln["mat_id"]))
        self.mat_combo.blockSignals(False)
        if self.type_combo is not None:
            self.type_combo.blockSignals(True)
            self.type_combo.setCurrentIndex(self._type_index(ln.get("kind")))
            self.type_combo.blockSignals(False)
        self.size_edit.blockSignals(True)
        _sz = ln.get("size")
        self.size_edit.setText("" if _sz is None else f"{float(_sz):g}")
        self.size_edit.blockSignals(False)
        self._sync_type_enabled()

    def _sync_type_enabled(self):
        """Grey the Mat ID + material name for any Type other than 'material',
        mirroring the polygon sheet (its row-7 name echo blanks on the same test):
        an overlay is not a soil zone, so it has no material."""
        is_mat = self._kind_at() == "material"
        self.mat_combo.setEnabled(is_mat)
        self._mat_label.setEnabled(is_mat)

    def _on_select(self, new_idx):
        self._commit_current()
        self._cur = new_idx
        self._load(new_idx)
        self._schedule_preview()

    def _on_mat_changed(self, idx):
        mid = self._mat_id_at(idx)
        if 0 <= self._cur < len(self._lines) and mid is not None:
            self._lines[self._cur]["mat_id"] = mid
            item = self.list.item(self._cur)
            if item:
                item.setText(self._label(self._cur))
        self._schedule_preview()

    def _on_type_changed(self, _idx):
        self._sync_type_enabled()
        if 0 <= self._cur < len(self._lines):
            self._lines[self._cur]["kind"] = self._kind_at()
            item = self.list.item(self._cur)
            if item:
                item.setText(self._label(self._cur))
        self._schedule_preview()

    def _on_size_changed(self, _text):
        ok, size = _parse_opt_size(self.size_edit.text())
        if ok and 0 <= self._cur < len(self._lines):
            self._lines[self._cur]["size"] = size
            item = self.list.item(self._cur)
            if item:
                item.setText(self._label(self._cur))

    def _add_line(self):
        self._commit_current()
        self._lines.append({"mat_id": 0, "coords": [], "kind": "material",
                            "size": None})
        self._refresh_list()
        self.list.setCurrentRow(len(self._lines) - 1)
        self._schedule_preview()

    def _remove_line(self):
        idx = self.list.currentRow()
        if idx < 0:
            return
        self._lines.pop(idx)
        self._cur = -1  # invalidate so the next select doesn't write to a stale index
        self._refresh_list()
        if self._lines:
            self.list.setCurrentRow(min(idx, len(self._lines) - 1))
        else:
            self._load(-1)
        self._schedule_preview()

    def result_lines(self):
        self._commit_current()
        # 'size' is emitted only when it is declared, so a blank Size stays absent
        # from the record rather than becoming an explicit None the writer would have
        # to special-case. 'kind' travels only on the polygon sheet, where the owning
        # editor splits the list back apart by it; a profile line has no Type, so
        # adding the key there would invent a field the record never had.
        out = []
        for ln in self._lines:
            item = {"coords": list(ln["coords"]), "mat_id": ln["mat_id"]}
            if self._poly_types:
                item["kind"] = ln.get("kind") or "material"
            if ln.get("size") is not None:
                item["size"] = ln["size"]
            out.append(item)
        return out

    def result_max_depth(self):
        """The edited max-depth value (float), or None if the field isn't shown."""
        if self._max_depth_edit is None:
            return None
        text = self._max_depth_edit.text()
        if _unedited(text, self._max_depth):     # untouched: keep the stored value
            return self._max_depth
        try:
            return float(text)
        except (TypeError, ValueError):
            return None

    def accept(self):
        """Validate the two Size rules before closing: a Size must be a positive
        number, and a 'refine' polygon must have one. Both are load-time errors in
        the engine — catching them here means the user fixes the item they are
        already looking at instead of meeting a traceback on the next Open."""
        from PySide6.QtWidgets import QMessageBox
        self._commit_current()
        title = self.windowTitle()
        for i, ln in enumerate(self._lines):
            name = f"{self._item_label} {i + 1}"
            if i == self._cur:
                ok, _ = _parse_opt_size(self.size_edit.text())
                if not ok:
                    QMessageBox.warning(
                        self, title,
                        f"{name}: Size must be a positive number, or blank for the "
                        "global target element size.")
                    return
            if (ln.get("kind") == "refine") and ln.get("size") is None:
                QMessageBox.warning(
                    self, title,
                    f"{name} has Type 'refine' but no Size. A refine region carries "
                    "no material and no analysis meaning — its only effect is the "
                    "local element size, so the Size is required.")
                self.list.setCurrentRow(i)
                return
        super().accept()


PROFILE_HELP = {
    "x": "X-coordinate of a point on this profile line, listed left→right.",
    "y": "Y-coordinate (elevation) of a point on this profile line.",
    "mat_id": "Material assigned to the layer below this line and above the next "
             "profile line down.",
    "max_depth": "Elevation of the model's bottom boundary (bedrock). Profile "
                "lines and the failure surface cannot go below it.",
    "size": "Optional target finite-element size along this line, used only when a "
           "mesh is generated. Blank = the global target size. A Size only ever "
           "refines: a value at or above the global size cannot coarsen the mesh.",
}


class ProfileEditor(CategoryEditor):
    label = "Profile lines"

    def build(self, slope_data, parent, select=None):
        style = _doc_style(parent)

        def preview(ax, lines, selected, max_depth):
            _draw_profile_preview(ax, lines, selected, max_depth, slope_data, style)

        return MatGeometryDialog(
            "Profile lines",
            "Each profile line is the top of a material layer, drawn left→right and "
            "ordered shallowest first. Select a line to edit its material and vertices.",
            "Line",
            slope_data.get("profile_lines") or [],
            slope_data.get("materials") or [], parent,
            select=select, max_depth=slope_data.get("max_depth") or 0.0,
            preview_draw=preview,
            preview_caption="Preview shows the pending profile lines over the current "
                            "model (selected line highlighted). Material zones aren't "
                            "re-filled until you save. Click a line or vertex to select it.",
            slope_data=slope_data, style=style,
            pick_resolve=lambda x, y, tol, lines: _pick_matgeom_lines(
                lines, x, y, tol, closed=False),
            field_help=PROFILE_HELP)

    def apply(self, slope_data, dlg):
        slope_data["profile_lines"] = dlg.result_lines()
        md = dlg.result_max_depth()
        if md is not None:
            slope_data["max_depth"] = md
        _resync_geometry(slope_data)  # rebuild polygons / ground surface / t-crack


POLYGON_HELP = {
    "x": "X-coordinate of a polygon vertex (CW or CCW order; the ring closes "
        "automatically — don't repeat the start point).",
    "y": "Y-coordinate of a polygon vertex.",
    # Each entry is <= the t_cut budget (393 chars — MEASURED to wrap in exactly two
    # lines at the dialog's natural width; the strip is fixed at two lines and clips
    # beyond).
    "mat_id": ("Material assigned to this closed zone. Greyed out for any Type other than 'material': an overlay is not a soil zone, so it has no material — the same rule the polygon sheet applies when it blanks the material-name echo."),
    "type": ("What kind of region this polygon is. 'material' (the default) is a soil zone. The three SSR types are FEM analysis OVERLAYS — never meshed, never sliced: 'ssr reduce' reduces only inside, 'ssr hold' holds full strength inside, 'ssr elastic' cannot yield inside. 'refine' is neither — it is a pure meshing region and needs a Size."),
    "size": ("Optional target finite-element size inside this polygon, used only when a mesh is generated. Blank = the global target size. Independent of Type: a material zone or an SSR overlay may carry one, and a 'refine' polygon is nothing but one. A Size only ever refines — a value at or above the global size cannot coarsen the mesh."),
}


class PolygonEditor(CategoryEditor):
    """Geometry editor for polygon-based files (no profile sheet). Each polygon is
    a closed material zone; the loader closes the ring implicitly, so the editor
    shows/stores only the distinct exterior vertices.

    The same list also carries the polygon sheet's OVERLAY rows — the SSR zones and
    the v21 refine regions — distinguished by the v21 Type. They live on
    ``slope_data['ssr_zones']`` / ``['refine_zones']``, not in ``polygons``: an SSR
    zone is an analysis overlay and a refine region is a meshing overlay, and neither
    may reach the mesher as geometry. build() merges the three lists for editing (in
    that order, which picking.py mirrors) and apply() splits them apart again."""

    label = "Polygons"

    def build(self, slope_data, parent, select=None):
        items = []
        for p in (slope_data.get("polygons") or []):
            coords = list(p["polygon"].exterior.coords)
            if len(coords) >= 2 and coords[0] == coords[-1]:
                coords = coords[:-1]                       # drop the closing duplicate
            # A negative Mat ID is the v20 sentinel encoding of an SSR overlay. The
            # loader splits those out for itself, so this only fires for a record
            # built by some other path; reading it as a kind keeps such a record from
            # being edited as "material -1".
            mid = p.get("mat_id")
            kind = SSR_ZONE_SENTINELS.get(mid, "material") if mid is not None else "material"
            items.append({"mat_id": None if kind != "material" else mid,
                          "kind": kind, "coords": coords, "size": p.get("size")})
        for z in (slope_data.get("ssr_zones") or []):
            kind = str(z.get("kind", "")).strip()
            if kind not in SSR_ZONE_LABELS:
                continue
            items.append({"mat_id": None, "kind": kind, "size": z.get("size"),
                          "coords": [tuple(c) for c in (z.get("polygon") or [])]})
        for r in (slope_data.get("refine_zones") or []):
            items.append({"mat_id": None, "kind": "refine", "size": r.get("size"),
                          "coords": [tuple(c) for c in (r.get("polygon") or [])]})
        style = _doc_style(parent)

        def preview(ax, polys, selected, _max_depth):
            _draw_polygon_preview(ax, polys, selected, slope_data, style)

        return MatGeometryDialog(
            "Polygons",
            "Each polygon is a closed region (the ring closes automatically, so list "
            "each vertex once) — a material zone, an SSR analysis overlay, or a mesh "
            "refinement region, set by its Type. Select a polygon to edit it.",
            "Polygon", items, slope_data.get("materials") or [], parent, select=select,
            preview_draw=preview,
            preview_caption="Preview shows the pending zones (selected zone filled and "
                            "outlined, others dimmed). Click a zone or vertex to select it.",
            slope_data=slope_data, style=style,
            pick_resolve=lambda x, y, tol, lines: _pick_matgeom_lines(
                lines, x, y, tol, closed=True),
            field_help=POLYGON_HELP, poly_types=True)

    def apply(self, slope_data, dlg):
        from shapely.geometry import Polygon
        polys, zones, refines = [], [], []
        for it in dlg.result_lines():
            coords = it["coords"]
            if len(coords) < 3:
                continue                                   # not a valid ring; skip
            kind = it.get("kind") or "material"
            if kind == "refine":
                # A refine region with no Size is refused by accept(), so one cannot
                # reach here from the dialog; drop it defensively rather than write a
                # record the loader would reject on the next Open.
                if it.get("size") is None:
                    continue
                refines.append({"polygon": [tuple(c) for c in coords],
                                "size": it["size"]})
            elif kind in SSR_ZONE_LABELS:
                zones.append({"kind": kind, "polygon": [tuple(c) for c in coords],
                              "label": SSR_ZONE_LABELS[kind],
                              "size": it.get("size")})
            else:
                polys.append({"polygon": Polygon(coords), "mat_id": it["mat_id"],
                              "size": it.get("size")})
        slope_data["ssr_zones"] = zones
        slope_data["refine_zones"] = refines
        _set_derived_geometry(slope_data, polys)


# --- reinforcement ---------------------------------------------------------- #
def _new_reinf():
    # A generic tensile line (blank type -> tangent/active via the loader presets,
    # spacing 1 -> no per-width division); mirrors the pre-v12 default row.
    return {"label": "", "x1": 0.0, "y1": 0.0, "x2": 0.0, "y2": 0.0,
            "t_max": 0.0, "t_res": 0.0,
            "lp1": 0.0, "lp2": 0.0, "E": 0.0, "area": 0.0,
            "type": "", "dir": "tangent", "appl": "active",
            "tend1": 0.0, "tend2": 0.0, "spacing": 1.0,
            # Blank, not zero: a new line uses the development-length law, and a
            # zero Adhesion with a zero Delta would be a real (and useless)
            # overburden law rather than the absence of one.
            "adhesion": float("nan"), "delta": float("nan")}


# List-view form layout for a reinforcement line: every ReinforcementEditor.FIELDS
# key grouped (Identity / Geometry / Capacity / Anchorage / Type). The Type combos
# are blank-tolerant exactly as the table enums are (a "" type is a real, selectable
# entry).
_REINF_FORM_GROUPS = [
    ("Identity", [["label"]]),
    ("Geometry", [["x1", "y1"], ["x2", "y2"]]),
    ("Capacity", [["t_max", "t_res"], ["E", "area"]]),
    ("Anchorage", [["lp1", "lp2"], ["adhesion", "delta"],
                   ["tend1", "tend2"], ["spacing"]]),
    ("Type", [["type"], ["dir", "appl"]]),
]


def _reinf_item_label(i, row):
    typ = str(row.get("type") or "line")
    try:
        return (f"{i + 1}. {typ} (x={float(row.get('x1', 0) or 0):g}"
                f"→{float(row.get('x2', 0) or 0):g})")
    except (TypeError, ValueError):
        return f"{i + 1}. {typ}"


# Wording mirrors the 'reinforce' worksheet section of input_template.md. Tmax,
# Lp1/Lp2, Adhesion/Delta, Tend1/Tend2 and Spacing form the capacity envelope used
# by BOTH LEM and FEM (fem.py caps the truss yield force at the same envelope) even
# though their header color is "LEM only" red — so they aren't tagged that way here.
# Lp1/Lp2 and Adhesion/Delta are two ways to state ONE thing, the pullout law: a
# filled Adhesion/Delta pair takes over and Lp1/Lp2 stop being read, which is what
# the list view's Pullout switch and the grayed pair in either view say. Type/Dir/
# Appl are truly LEM only (FEM ignores them); Tres/E/Area are truly FEM only.
REINFORCE_HELP = {
    "label": "Name used in error messages, summaries, and plots (optional).",
    "x1": "Start point X-coordinate.",
    "y1": "Start point Y-coordinate.",
    "x2": "End point X-coordinate.",
    "y2": "End point Y-coordinate.",
    "type": "Support-type preset — fills Dir and Appl automatically (Geosynthetic, "
           "Nail, Tieback, Anchor); blank = generic line. LEM only.",
    "dir": "Force direction at the slip surface — Tangent (flexible, e.g. "
          "geosynthetics; the default) or Axial (rigid, e.g. nails/tiebacks — the "
          "UTEXAS/UTEXASED convention). LEM only.",
    "appl": "Force application — Active: allowable force on the driving side, not "
           "divided by FS (default). Passive: ultimate capacity on the resisting "
           "side, divided by FS. LEM only.",
    "t_max": "Maximum tensile force the line can mobilize, per unit width (discrete "
            "supports: enter the per-element capacity with Spacing). Caps both the "
            "LEM force and the FEM yield force.",
    "t_res": "Residual tensile force the line retains after it ruptures (post-peak), "
            "per unit width (÷ Spacing for discrete supports). Capped by the pullout "
            "envelope: bond slip is perfectly plastic, so an element keeps carrying "
            "whatever its embedment develops. FEM only; blank = elastic-perfectly-"
            "plastic (holds capacity), 0 = brittle rupture (carries nothing).",
    "tend1": "End anchorage/connection capacity at end 1, per unit width (0 = "
            "friction only; ÷ Spacing for discrete supports).",
    "tend2": "End anchorage/connection capacity at end 2, per unit width (0 = "
            "friction only; ÷ Spacing for discrete supports).",
    "lp1": "Pullout bond length at end 1 — tapers the mobilized force toward that end.",
    "lp2": "Pullout bond length at end 2 — tapers the mobilized force toward that end.",
    "adhesion": "Adhesion = soil-reinforcement interface adhesion (blank = use Lp).",
    "delta": "Delta = soil-reinforcement interface friction angle (blank = use Lp).",
    "spacing": "Out-of-plane spacing for discrete supports (nails, tiebacks); leave "
              "blank or 1 for geosynthetics (already per unit width). All capacity "
              "terms and Area are divided by it, once, for both engines.",
    "E": "Elastic modulus of the reinforcement; with Area gives the axial stiffness "
        "EA (Area carries the ÷ Spacing). FEM only (models the line as a 1D truss).",
    "area": "Cross-sectional area of the reinforcement, per unit width (÷ Spacing). "
           "FEM only.",
}


class ReinforcementEditor(CategoryEditor):
    label = "Reinforcement"
    LF = {"lem", "fem"}
    # Support-type columns (v12). `type` choices are REINFORCE_TYPE_PRESETS' keys
    # (the loader's, which are the sheet's) plus a blank entry — a blank type is a
    # generic line whose Dir/Appl default via the presets; offering '' as an empty
    # combo entry lets a blank type round-trip unchanged (same treatment as the
    # materials blank option). Dir/Appl mirror the loader's accepted values.
    # tooltip=REINFORCE_HELP[key] gives the table-header hover and the list-view
    # label/edit hover; the same dict feeds the context-sensitive help strip.
    # Column order is the reinforce sheet's, Label first, so a block copied from the
    # sheet or the docs' tables pastes straight in.
    FIELDS = [
        Field("label", "Label", "str", tooltip=REINFORCE_HELP["label"]),
        Field("x1", "x1", tooltip=REINFORCE_HELP["x1"]),
        Field("y1", "y1", tooltip=REINFORCE_HELP["y1"]),
        Field("x2", "x2", tooltip=REINFORCE_HELP["x2"]),
        Field("y2", "y2", tooltip=REINFORCE_HELP["y2"]),
        Field("type", "Type", "choice",
              choices=[""] + list(REINFORCE_TYPE_PRESETS), applies=LF,
              tooltip=REINFORCE_HELP["type"]),
        Field("dir", "Dir", "choice", choices=["tangent", "axial"], applies=LF,
              tooltip=REINFORCE_HELP["dir"]),
        Field("appl", "Appl", "choice", choices=["active", "passive"], applies=LF,
              tooltip=REINFORCE_HELP["appl"]),
        # Field order past Type/Dir/Appl mirrors the reinforce sheet's columns
        # (Tmax, Lp1, Lp2, Tend1, Tend2, Spacing, then the FEM-only tail) so a
        # block copied from the sheet or the docs' tables pastes straight in.
        Field("t_max", "Tmax", usage="lem", tooltip=REINFORCE_HELP["t_max"]),
        Field("lp1", "Lp1", usage="lem", tooltip=REINFORCE_HELP["lp1"]),
        Field("lp2", "Lp2", usage="lem", tooltip=REINFORCE_HELP["lp2"]),
        # The overburden pullout law. applies=LF, like Spacing: both engines read
        # it, so neither usage toggle may hide it. The sheet colors the whole
        # envelope block one red; the editor leaves the both-engines members of
        # that block uncolored, as it already does for Spacing.
        # kind="optfloat", not "float": a cleared cell must come back BLANK, and a
        # blank pair is what selects the development-length law. A plain float
        # would read a cleared cell as 0.0 — a real, and refused, overburden law.
        Field("adhesion", "Adhesion", "optfloat", applies=LF, unit="stress",
              tooltip=REINFORCE_HELP["adhesion"]),
        Field("delta", "Delta", "optfloat", applies=LF,
              tooltip=REINFORCE_HELP["delta"]),
        Field("tend1", "Tend1", usage="lem", tooltip=REINFORCE_HELP["tend1"]),
        Field("tend2", "Tend2", usage="lem", tooltip=REINFORCE_HELP["tend2"]),
        Field("spacing", "Spacing", applies=LF, tooltip=REINFORCE_HELP["spacing"]),
        Field("t_res", "Tres", usage="fem", tooltip=REINFORCE_HELP["t_res"]),
        # E is a Young's modulus (stress). Area (a cross-section, length²) has no
        # xslope.units.labels() key -- the labels() contract is length/stress/
        # unit_weight/force_per_len/k/flowrate/time -- so it is deliberately left
        # unlabeled by the static header; the list view labels it dynamically as
        # per-element (L²) / per-unit-width (L²/L) instead.
        Field("E", "E", usage="fem", unit="stress", tooltip=REINFORCE_HELP["E"]),
        Field("area", "Area", usage="fem", tooltip=REINFORCE_HELP["area"]),
    ]
    # Which pullout law a line uses is not a stored field — the row says it. A
    # filled Adhesion/Delta pair IS the overburden law, and the development lengths
    # stop being read; anything else is the development-length law, with Adhesion
    # and Delta blank and free to be typed. One function answers that question for
    # the table's graying, the list view's switch, and the preflight note, so the
    # three cannot disagree.
    @staticmethod
    def pullout_mode(row):
        """'overburden' when this row carries both Adhesion and Delta, else 'lp'."""
        def _filled(v):
            if v is None or (isinstance(v, str) and not v.strip()):
                return False
            try:
                return float(v) == float(v)          # NaN is blank
            except (TypeError, ValueError):
                return False
        return ("overburden" if _filled(row.get("adhesion")) and _filled(row.get("delta"))
                else "lp")

    @classmethod
    def dim_keys(cls, row):
        """Field keys to gray for a reinforcement row.

        Only the pair the row is NOT using is grayed, and only when the other pair
        is actually in force: a line on the development-length law leaves Adhesion
        and Delta live, because graying an empty cell is graying out the only way
        to fill it.
        """
        if cls.pullout_mode(row) == "overburden":
            return frozenset(("lp1", "lp2"))
        return frozenset()

    SWITCH_SPEC = {
        "group": "Anchorage",
        "label": "Pullout",
        # (id, label, fields, self-asserting). Adhesion/Delta are self-asserting:
        # values left in them WOULD keep the overburden law in force, so leaving
        # that option parks them. Lp1/Lp2 are inert while the other law runs, so
        # they simply dim -- their values stay on screen and in the file.
        "options": [("lp", "Development length (Lp1, Lp2)", ("lp1", "lp2"), False),
                    ("overburden", "Overburden (Adhesion, Delta)",
                     ("adhesion", "delta"), True)],
        "resolve": lambda row: ReinforcementEditor.pullout_mode(row),
        "tooltip": "Which pullout law this line uses. Development length develops "
                   "the full capacity over Lp1/Lp2. Overburden develops it at "
                   "2·(Adhesion + σ′v·tan Delta) per unit "
                   "length, from the effective overburden at each point of the "
                   "line.",
    }

    # The sheet's Dir/Appl formulas, in the editor: picking a Type fills both from
    # the loader's own table; typing over either keeps it until a Type is picked
    # again. Handed to BOTH views (and used by a pasted Type), so the three ways a
    # reader can set a support type mean one thing.
    PRESET_SPEC = {"driver": "type", "fills": ("dir", "appl"),
                   "presets": REINFORCE_TYPE_PRESETS}

    def build(self, slope_data, parent):
        style = _doc_style(parent)

        def preview(ax, rows, selected):
            _draw_reinforcement_preview(ax, rows, selected, slope_data, style)

        return _LineEditorDialog(
            "Reinforcement", self.FIELDS, slope_data.get("reinforcement_lines", []),
            _new_reinf, _REINF_FORM_GROUPS, _reinf_item_label, preview,
            lambda x, y, tol, rows: _pick_line_rows(rows, x, y, tol),
            view_state="reinforce", parent=parent, unit_labels=_unit_labels_for(slope_data),
            dynamic_spec={"driver": "spacing",
                          "fields": {"t_max": "force", "t_res": "force",
                                     "tend1": "force", "tend2": "force",
                                     "area": "area"}},
            preset_spec=self.PRESET_SPEC,
            dim_rule=self.dim_keys,
            switch_spec=self.SWITCH_SPEC,
            help_text="List view edits one line at a time as a grouped form beside a "
                      "live section preview; the table view is available for bulk entry "
                      "of the many lines of a tiered wall. Both views edit the same rows, "
                      "so switching is lossless. A line develops its pullout capacity one "
                      "of two ways, chosen per line by the Pullout switch: over the "
                      "development lengths Lp1/Lp2 (0 = fully anchored), or from the "
                      "effective overburden through Adhesion and Delta. The pair not in "
                      "use is grayed, values intact. The LEM tension distribution on the "
                      "preview is derived from whichever is in force. Picking a Type "
                      "fills Dir and Appl "
                      "— change either afterwards to override it — and a blank Type is a "
                      "generic tensile line. Tend1/Tend2 are the "
                      "end-anchorage capacities; capacities and E/Area are per-unit-width "
                      "(Spacing divides discrete supports).",
            usage_toggles=["lem", "fem"],
            preview_caption="Preview shows the reinforcement lines on the section "
                            "(selected line bold with its endpoints and pullout-length "
                            "markers; others dimmed). Click a line to select it.",
            field_help=REINFORCE_HELP)

    def apply(self, slope_data, dlg):
        from xslope.fileio import build_reinforce_lines, attach_reinforce_pullout
        rows = dlg.result_rows()
        slope_data["reinforcement_lines"] = rows
        # Rebuild the LEM display/analysis format so the canvas reflects the edit.
        slope_data["reinforce_lines"] = build_reinforce_lines(rows)
        # An edited Adhesion/Delta changes the envelope's shape, not just its
        # numbers, so the pullout profiles are rebuilt here rather than waiting
        # for a solve — the canvas draws the curve the analysis will use.
        attach_reinforce_pullout(slope_data)


# --- line loads ------------------------------------------------------------- #
def _new_lload():
    # Same key set load_slope_data produces (fileio.py:1325-1330). Angle defaults
    # to -90° (straight down), mirroring the loader's blank-angle fallback.
    return {"x": 0.0, "y": 0.0, "P": 0.0, "angle": -90.0, "label": "Load"}


LLOADS_HELP = {
    "label": "Name used in error messages, summaries, and plots (optional).",
    "x": "X-coordinate of the point of application — must lie on the ground "
        "surface (snapped on save).",
    "y": "Y-coordinate of the point of application — must lie on the ground "
        "surface (snapped on save).",
    "P": "Force magnitude, per unit width of slope (force/length).",
    "angle": "Direction of the force from horizontal, degrees (−90 = straight "
            "down, the default).",
}


class LineLoadsEditor(CategoryEditor):
    label = "Line loads"
    # Loader keys (fileio.py:1325-1330): x, y, P (magnitude), angle, label. A line
    # load is a concentrated force per unit width acting at a point on the ground
    # surface; used by both LEM (slice.py) and FEM (fem.py), so no usage tags.
    FIELDS = [
        Field("label", "Label", "str"),
        Field("x", "x"), Field("y", "y"),
        Field("P", "P"),
        Field("angle", "Angle", default=-90.0),
    ]

    def build(self, slope_data, parent):
        style = _doc_style(parent)

        def preview(ax, rows, selected):
            _draw_line_loads_preview(ax, rows, selected, slope_data, style)

        return TableEditorDialog(
            "Line loads", self.FIELDS, slope_data.get("line_loads", []),
            _new_lload, parent,
            help_text="A concentrated line load (force per unit width) acting at a "
                      "point on the ground surface — e.g. the weight of a facing "
                      "plate. P is the magnitude (positive); Angle is the direction "
                      "in degrees (−90 = straight down). Each load is snapped onto "
                      "the ground surface on save, since the loader requires line "
                      "loads to act on the surface.",
            preview_draw=preview,
            preview_caption="Preview shows the line loads on the section (selected "
                            "load's arrow emphasized; others dimmed). The arrow points "
                            "in the force direction, head on the point of application. "
                            "Click a load's arrow to select it.",
            pick_resolve=lambda x, y, tol, rows: _pick_line_loads(rows, x, y, tol, slope_data),
            field_help=LLOADS_HELP)

    def apply(self, slope_data, dlg):
        rows = dlg.result_rows()
        # Mirror the loader (fileio.py:1315-1324), which requires a line load to act
        # ON the ground surface: snap each point onto the nearest surface location so
        # Studio can never author a load the loader would refuse. (The loader also
        # rejects points beyond a tolerance; snapping is the forgiving inverse — the
        # saved point always lands exactly on the surface.)
        gs = slope_data.get("ground_surface")
        if gs is not None and not gs.is_empty:
            from shapely.geometry import Point
            for ll in rows:
                try:
                    snapped = gs.interpolate(gs.project(Point(ll["x"], ll["y"])))
                    ll["x"], ll["y"] = float(snapped.x), float(snapped.y)
                except Exception:
                    pass
        slope_data["line_loads"] = rows


# --- transient seepage (the v18 'tseep' sheet) ----------------------------- #
# Five series columns (C..G), matching the v18 tseep sheet geometry and the parser/
# writer in fileio (_parse_tseep_sheet / save_slope_data_to_xlsx).
_TSEEP_N_SERIES = 5


def _tseep_optf(text):
    """Blank / 'none' -> None; otherwise the float. Mirrors the optfloat field so an
    unset control (duration, a series breakpoint) round-trips as None, not 0.0."""
    text = (text or "").strip()
    if text == "" or text.lower() == "none":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _tseep_fmt(v):
    """Display a stored float with enough precision to round-trip exactly, but clean
    for round numbers (100.0 -> '100'). Blank for None."""
    return "" if v is None else f"{v:.10g}"


TSEEP_HELP = {
    "duration": "Total simulated time of the run, in the model's declared Time unit. "
                "The march ends here — scheduled times beyond the duration are never "
                "reached.",
    "save_interval": "Spacing of the regularly-saved frames (one every interval, up to "
                     "the duration). The saved set is the UNION of this grid, the extra "
                     "save times, the stage times, and every series' own breakpoints.",
    "stage_1": "Rapid-drawdown STAGE 1 time: the saved frame whose pore pressures seed "
               "the first drawdown stage of an LEM/FEM run (u = seep). Leave BOTH "
               "stages blank for a plain transient run with no drawdown coupling.",
    "stage_2": "Rapid-drawdown STAGE 2 time — must be later than Stage 1. Set BOTH "
               "stage times, or neither.",
    "stability_time": "Which instant an LEM or FEM run with u = seep reads its pore "
                      "pressures from. It selects a frame out of the march; it does "
                      "not change the march. Leave it blank and a run reads the LAST "
                      "saved frame. The Run LEM and Run FEM dialogs set it too — this "
                      "is the same value.",
    "time": "The shared time axis for every series (ascending, in the Time unit). A "
            "REPEATED time is an instantaneous step: the series jumps to the new value "
            "at that instant (right-continuous), drawn as a vertical segment on the plot.",
    "series": "This series' value at the row's time. Blank = no breakpoint here; the "
              "series is linearly interpolated between its own defined points. Each "
              "named series drives any seep BC head/flux VALUE cell that contains its "
              "name (a time-varying boundary condition).",
    "series_name": "Name of this time series. A seep BC value cell holding this exact "
                   "name is driven by the series (a time-varying BC). A blank name "
                   "column is unused.",
    "save_times": "Extra explicit times to save a frame at, beyond the save-interval "
                  "grid — e.g. a specific instant of interest.",
}


class TransientDialog(QDialog):
    """Editor for the transient-seepage (``tseep``) inputs: run controls (duration,
    save interval, rapid-drawdown stage times, the stability time), the
    extra-save-times list, and the
    time-series table (a shared time axis plus up to five named series whose values
    drive the seep BC value cells that name them). A live plot beside the table draws
    every defined series versus time — markers at the breakpoints, linear between,
    stage times as reference lines — and updates as the tables are edited.

    There is no enable toggle: blank fields already mean "no transient data" — the
    writer emits a blank tseep sheet and the file loads steady, round-tripping
    untouched — and the steady-vs-transient RUN choice lives on the Run Seepage
    dialog. :meth:`result_tseep` therefore returns ``None`` for an all-blank editor,
    matching the fileio parser exactly (series aligned to the time axis with ``None``
    gaps; an all-blank sheet collapses to ``None``).

    Layout: the LEFT column is one time-data block — series names on top, then the
    series table and the extra-save-times list side by side. The RIGHT column stacks
    the five run-control fields (duration, save interval, stage times, stability
    time) at the top with
    the series plot below them, so the plot's canvas toolbar is never clipped."""

    def __init__(self, tseep, slope_data, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Transient seepage")
        self.resize(1080, 620)
        self._unit_labels = _unit_labels_for(slope_data)
        self._preview = None            # so an early schedule() is safe
        self._populating = True
        tseep = tseep or {}

        layout = QVBoxLayout(self)
        layout.addWidget(_help_label(
            "Transient (time-dependent) seepage. Define the run duration and save "
            "schedule, the optional rapid-drawdown stage times, and one or more named "
            "time series (a shared time axis with a value per series). A seep BC value "
            "cell that names a series is driven by it — a time-varying boundary "
            "condition. With no times or values entered the model stays steady. The "
            "plot previews the series as you edit."))

        # --- run controls (built here; placed at the TOP of the right column) ---
        controls_widget = QWidget()
        controls = QFormLayout(controls_widget)
        controls.setContentsMargins(0, 0, 0, 0)
        self._duration = QLineEdit(_tseep_fmt(tseep.get("duration")))
        self._duration.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self._save_interval = QLineEdit(_tseep_fmt(tseep.get("save_interval")))
        self._save_interval.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self._stage_1 = QLineEdit(_tseep_fmt(tseep.get("stage_1")))
        self._stage_1.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self._stage_2 = QLineEdit(_tseep_fmt(tseep.get("stage_2")))
        self._stage_2.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self._stability_time = QLineEdit(_tseep_fmt(tseep.get("stability_time")))
        tlab = f" ({self._unit_labels['time']})" if (self._unit_labels
                                                     and self._unit_labels.get("time")) else ""
        for lbl, w, key in (("Duration" + tlab, self._duration, "duration"),
                            ("Save interval" + tlab, self._save_interval, "save_interval"),
                            ("Stage 1 time" + tlab, self._stage_1, "stage_1"),
                            ("Stage 2 time" + tlab, self._stage_2, "stage_2"),
                            ("Stability time" + tlab, self._stability_time,
                             "stability_time")):
            w.setToolTip(TSEEP_HELP[key])
            controls.addRow(lbl, w)

        # --- main split: time-data block (left) | controls + plot (right) ---
        split = QSplitter(Qt.Horizontal)

        # LEFT — one time-data block: series names on top, then the series table and
        # the extra-save-times list side by side (save times to the right of the last
        # series column), so the whole left column reads as one time-data unit.
        series_group = QGroupBox("Time series")
        sgl = QVBoxLayout(series_group)
        names_row = QHBoxLayout()
        names_row.addWidget(QLabel("Series names:"))
        self._name_edits = []
        existing = list((tseep.get("series") or {}).keys())
        for i in range(_TSEEP_N_SERIES):
            if i < len(existing):
                name = existing[i]
            else:
                # Empty slots carry the input template's default header names
                # (tseep sheet row 2: t1..t5); a default that collides with a
                # loaded series name stays blank instead.
                default = f"t{i + 1}"
                name = default if default not in existing else ""
            e = QLineEdit(name)
            e.setToolTip(TSEEP_HELP["series_name"])
            e.setPlaceholderText(f"t{i + 1}")
            self._name_edits.append(e)
            names_row.addWidget(e)
        sgl.addLayout(names_row)

        tables_row = QHBoxLayout()

        series_col = QVBoxLayout()
        times = list(tseep.get("times") or [])
        series_vals = tseep.get("series") or {}
        n_rows = max(len(times) + 2, 8)
        self._series_table = QTableWidget(n_rows, 1 + _TSEEP_N_SERIES)
        _spreadsheet_manners(self._series_table)
        self._series_table.setMinimumHeight(220)   # show several breakpoints at once
        self._sync_series_headers()
        for r, t in enumerate(times):
            self._series_table.setItem(r, 0, QTableWidgetItem(_tseep_fmt(t)))
            for c, name in enumerate(existing):
                vals = series_vals.get(name) or []
                v = vals[r] if r < len(vals) else None
                if v is not None:
                    self._series_table.setItem(r, c + 1, QTableWidgetItem(_tseep_fmt(v)))
        series_col.addWidget(self._series_table, 1)
        srow = QHBoxLayout()
        add_s = QPushButton("Add row")
        add_s.clicked.connect(lambda: self._add_table_row(self._series_table))
        rem_s = QPushButton("Remove selected")
        rem_s.clicked.connect(lambda: self._remove_table_rows(self._series_table))
        srow.addWidget(add_s)
        srow.addWidget(rem_s)
        srow.addStretch(1)
        series_col.addLayout(srow)
        tables_row.addLayout(series_col, 1)

        # Extra save times — a TALL vertical list immediately to the right of the
        # series table (no longer a short bottom strip), so more rows are visible.
        save_group = QGroupBox("Extra save times")
        save_group.setMaximumWidth(170)
        svl = QVBoxLayout(save_group)
        save_times = list(tseep.get("save_times") or [])
        self._save_table = QTableWidget(max(len(save_times) + 1, 8), 1)
        _spreadsheet_manners(self._save_table)
        self._save_table.setHorizontalHeaderLabels(["save time"])
        self._save_table.horizontalHeaderItem(0).setToolTip(TSEEP_HELP["save_times"])
        for r, v in enumerate(save_times):
            self._save_table.setItem(r, 0, QTableWidgetItem(_tseep_fmt(v)))
        svl.addWidget(self._save_table, 1)
        vrow = QHBoxLayout()
        add_v = QPushButton("Add")
        add_v.clicked.connect(lambda: self._add_table_row(self._save_table))
        rem_v = QPushButton("Remove")
        rem_v.clicked.connect(lambda: self._remove_table_rows(self._save_table))
        vrow.addWidget(add_v)
        vrow.addWidget(rem_v)
        svl.addLayout(vrow)
        tables_row.addWidget(save_group)

        sgl.addLayout(tables_row, 1)
        split.addWidget(series_group)

        # RIGHT — run controls stacked at the top, series plot beneath them (so the
        # plot's Fit/+/- canvas toolbar stays unobstructed).
        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.addWidget(controls_widget)
        from .canvas import PreviewPane
        self._preview = PreviewPane(self._draw_plot, dxf=False,
                                    caption="Series value vs time. Markers are the "
                                            "defined breakpoints; stage times show as "
                                            "dashed reference lines.")
        rv.addWidget(self._preview, 1)
        split.addWidget(right)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 1)
        split.setSizes([680, 520])
        layout.addWidget(split, 1)

        _ok_cancel(self, layout)
        attach_help(self, TSEEP_HELP, self._help_resolver)
        self._series_table.currentCellChanged.connect(self._on_series_cell)

        # Wire live updates AFTER the initial population so building fires nothing.
        self._populating = False
        for w in (self._duration, self._save_interval, self._stage_1, self._stage_2):
            w.textChanged.connect(self._schedule)
        for i, e in enumerate(self._name_edits):
            e.textChanged.connect(self._sync_series_headers)
            e.textChanged.connect(self._schedule)
        self._series_table.itemChanged.connect(self._schedule)
        self._save_table.itemChanged.connect(self._schedule)

        self._preview.refresh_now()

    # --- helpers ------------------------------------------------------
    @staticmethod
    def _cell_text(table, r, c):
        it = table.item(r, c)
        return it.text() if it is not None else ""

    def _sync_series_headers(self, *_):
        labels = ["time"]
        for i, e in enumerate(self._name_edits):
            nm = e.text().strip()
            labels.append(nm if nm else f"(t{i + 1})")
        self._series_table.setHorizontalHeaderLabels(labels)
        h0 = self._series_table.horizontalHeaderItem(0)
        if h0 is not None:
            h0.setToolTip(TSEEP_HELP["time"])
        for c in range(1, 1 + _TSEEP_N_SERIES):
            hi = self._series_table.horizontalHeaderItem(c)
            if hi is not None:
                hi.setToolTip(TSEEP_HELP["series_name"])

    def _add_table_row(self, table):
        table.insertRow(table.rowCount())
        self._schedule()

    def _remove_table_rows(self, table):
        rows = sorted({ix.row() for ix in table.selectedIndexes()}, reverse=True)
        for r in rows:
            table.removeRow(r)
        self._schedule()

    def _schedule(self, *_):
        if not self._populating and self._preview is not None:
            self._preview.schedule()

    # --- help strip ----------------------------------------------------
    def _help_resolver(self, widget):
        direct = {self._duration: "duration",
                  self._save_interval: "save_interval", self._stage_1: "stage_1",
                  self._stage_2: "stage_2"}
        if widget in direct:
            return direct[widget]
        if widget in self._name_edits:
            return "series_name"
        w = widget
        while w is not None:
            if w is self._series_table:
                return "time" if self._series_table.currentColumn() == 0 else "series"
            if w is self._save_table:
                return "save_times"
            w = w.parentWidget()
        return None

    def _on_series_cell(self, row, col, *_):
        # Arrow-key nav moves the current cell without a focus change; push the help.
        if getattr(self, "_help_strip", None) is not None:
            key = "time" if col == 0 else "series"
            self._help_strip.set_help(TSEEP_HELP.get(key, ""))

    # --- live plot -----------------------------------------------------
    def _pending_rows(self):
        rows = []
        for r in range(self._series_table.rowCount()):
            t = _tseep_optf(self._cell_text(self._series_table, r, 0))
            if t is None:
                continue
            vals = [_tseep_optf(self._cell_text(self._series_table, r, c + 1))
                    for c in range(_TSEEP_N_SERIES)]
            rows.append((t, vals))
        return rows

    def _draw_plot(self, ax):
        ul = self._unit_labels
        t_unit = f" ({ul['time']})" if (ul and ul.get("time")) else ""
        len_unit = f" ({ul['length']})" if (ul and ul.get("length")) else ""
        ax.set_xlabel(f"Time{t_unit}")
        ax.set_ylabel(f"Series value{len_unit}")
        ax.grid(True, alpha=0.3)
        rows = self._pending_rows()
        names = [e.text().strip() for e in self._name_edits]
        any_series = False
        for c in range(_TSEEP_N_SERIES):
            nm = names[c]
            if not nm:
                continue
            pts = [(t, vals[c]) for t, vals in rows if vals[c] is not None]
            if not pts:
                continue
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            ax.plot(xs, ys, marker="o", ms=4, lw=1.5, label=nm)
            any_series = True
        # Rapid-drawdown stage times as dashed reference lines.
        for key, lbl in (("_stage_1", "Stage 1"), ("_stage_2", "Stage 2")):
            v = _tseep_optf(getattr(self, key).text())
            if v is not None:
                ax.axvline(v, color="crimson", ls="--", lw=1.0)
                ax.annotate(lbl, xy=(v, 1.0), xycoords=("data", "axes fraction"),
                            xytext=(3, -3), textcoords="offset points",
                            ha="left", va="top", fontsize=8, color="crimson")
        if any_series:
            ax.legend(fontsize=8, loc="best")
        else:
            ax.text(0.5, 0.5, "Enter time/value rows under a named series to plot it",
                    transform=ax.transAxes, ha="center", va="center", color="gray")

    # --- result --------------------------------------------------------
    def result_tseep(self):
        """Reconstruct the tseep dict (or ``None`` for an all-blank editor), matching
        the fileio parser: series aligned to the time axis with ``None`` gaps; a series
        registered only when it has a non-blank name AND at least one value."""
        rows = self._pending_rows()
        times = [t for t, _ in rows]
        names = [e.text().strip() for e in self._name_edits]
        series = {}
        for c in range(_TSEEP_N_SERIES):
            nm = names[c]
            col_vals = [vals[c] for _, vals in rows]
            if nm and any(v is not None for v in col_vals):
                series[nm] = col_vals
        save_times = []
        for r in range(self._save_table.rowCount()):
            v = _tseep_optf(self._cell_text(self._save_table, r, 0))
            if v is not None:
                save_times.append(v)
        duration = _tseep_optf(self._duration.text())
        save_interval = _tseep_optf(self._save_interval.text())
        stage_1 = _tseep_optf(self._stage_1.text())
        stage_2 = _tseep_optf(self._stage_2.text())
        stability_time = _tseep_optf(self._stability_time.text())
        # Mirror the parser's "present but all-blank -> None (steady)".
        if not (times or series or save_times
                or any(v is not None for v in (duration, save_interval, stage_1,
                                               stage_2, stability_time))):
            return None
        return {"times": times, "series": series, "duration": duration,
                "save_interval": save_interval, "save_times": save_times,
                "stage_1": stage_1, "stage_2": stage_2,
                "stability_time": stability_time}

    def accept(self):
        from PySide6.QtWidgets import QMessageBox
        ts = self.result_tseep()
        if ts is not None:
            times = ts["times"]
            if any(b < a for a, b in zip(times, times[1:])):
                QMessageBox.warning(self, "Transient seepage",
                                    "Times must be in ascending order (equal "
                                    "consecutive times are allowed — they define an "
                                    "instantaneous step).")
                return
            names = [e.text().strip() for e in self._name_edits if e.text().strip()]
            if len(names) != len(set(names)):
                QMessageBox.warning(self, "Transient seepage",
                                    "Series names must be unique.")
                return
            s1, s2 = ts["stage_1"], ts["stage_2"]
            if (s1 is None) != (s2 is None):
                QMessageBox.warning(self, "Transient seepage",
                                    "Set BOTH rapid-drawdown stage times, or neither.")
                return
            if s1 is not None and s2 is not None and s1 >= s2:
                QMessageBox.warning(self, "Transient seepage",
                                    "Stage 1 time must be less than Stage 2 time.")
                return
            st, dur = ts["stability_time"], ts["duration"]
            if st is not None and dur is not None and not (0 < st <= dur):
                QMessageBox.warning(
                    self, "Transient seepage",
                    "The stability time must be greater than 0 and no later than the "
                    "run duration — the march never reaches an instant outside it.")
                return
        super().accept()


class TransientEditor(CategoryEditor):
    label = "Transient"

    def build(self, slope_data, parent):
        return TransientDialog(slope_data.get("tseep"), slope_data, parent)

    def apply(self, slope_data, dlg):
        ts = dlg.result_tseep()
        if ts is None:
            slope_data.pop("tseep", None)
        else:
            slope_data["tseep"] = ts


CATEGORY_EDITORS = {
    "global": GlobalEditor(),
    "materials": MaterialsEditor(),
    "circles": CirclesEditor(),
    "non_circ": NonCircEditor(),
    "piezo": PiezoEditor(),
    "dloads": DloadsEditor(),
    "seep_bc": SeepBcEditor(),
    "piles": PilesEditor(),
    "reinforce": ReinforcementEditor(),
    "line_loads": LineLoadsEditor(),
    "profile": ProfileEditor(),
    "polygons": PolygonEditor(),
    "transient": TransientEditor(),
}
