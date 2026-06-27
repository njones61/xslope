"""Input editors — modal dialogs launched from the Inputs tree (Phase 2).

Each editable ``slope_data`` category has a small ``CategoryEditor`` that builds a
dialog from the current data and writes the edited values back. Most categories
reuse the generic ``TableEditorDialog`` (a column spec + add/remove rows);
scalars use ``FormEditorDialog``. Editors only expose the commonly-edited fields
of each record and **preserve any unshown keys** (e.g. material reliability
sigmas), so a round-trip through an editor never drops data.
"""

from __future__ import annotations

import math

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDialog, QDialogButtonBox, QFormLayout,
    QHBoxLayout, QLabel, QLineEdit, QListWidget, QPushButton, QTableWidget,
    QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

# Column "usage" tags: which analysis a field applies to. Header text is colored
# to mirror the input template's header coloring (red = LEM-specific inputs,
# blue = FEM-specific inputs); seepage/reliability extend the same idea.
USAGE_COLOR = {"lem": "#c00000", "fem": "#0432ff", "seep": "#2e7d32", "rel": "#9a7d0a"}
USAGE_NAME = {"lem": "LEM only", "fem": "FEM only",
              "seep": "Seepage only", "rel": "Reliability"}
# Short labels for the show/hide toggles at the top of an editor.
USAGE_TOGGLE_LABEL = {"lem": "LEM", "fem": "FEM", "seep": "Seepage", "rel": "Reliability"}


# --------------------------------------------------------------------------- #
# Field spec + dialogs
# --------------------------------------------------------------------------- #
class Field:
    """One editable column/field. kind ∈ {'float','optfloat','int','str','choice'}.

    'optfloat' is an optional float: a blank cell reads back as None (not 0.0), so
    optional fields like a pile's H or capacities stay unset instead of becoming
    zero (which the loader would reject)."""

    _BLANK = {"float": 0.0, "optfloat": None, "int": 0, "str": "", "choice": ""}

    def __init__(self, key, header, kind="float", choices=None, default=None,
                 usage=None, applies=None):
        self.key = key
        self.header = header
        self.kind = kind
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
            default = self.choices[0] if kind == "choice" and self.choices else self._BLANK[kind]
        self.default = default

    def from_text(self, text):
        text = (text or "").strip()
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


class FormEditorDialog(QDialog):
    """Simple key/value form for scalar parameters."""

    def __init__(self, title, fields, values, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._fields = fields
        self._edits = {}
        layout = QVBoxLayout(self)
        form = QFormLayout()
        for f in fields:
            edit = QLineEdit(str(values.get(f.key, f.default)))
            self._edits[f.key] = edit
            form.addRow(f.header, edit)
        layout.addLayout(form)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    def result_values(self):
        return {f.key: f.from_text(self._edits[f.key].text()) for f in self._fields}


class _EditableTable(QWidget):
    """A table over a list of dict records with Add/Remove rows. Unshown keys are
    preserved. Reused standalone (TableEditorDialog) and per-tab (TabbedTableEditorDialog)."""

    def __init__(self, fields, rows, new_row, parent=None):
        super().__init__(parent)
        self._fields = fields
        self._new_row = new_row
        self._bases = [dict(r) for r in rows]  # keep originals to preserve extra keys

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.table = QTableWidget(len(rows), len(fields))
        self.table.setHorizontalHeaderLabels([f.header for f in fields])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        for j, f in enumerate(fields):  # color usage-tagged headers (red=LEM, blue=FEM, …)
            if f.usage:
                hi = self.table.horizontalHeaderItem(j)
                if hi is not None:
                    hi.setForeground(QColor(USAGE_COLOR[f.usage]))
                    fnt = hi.font()
                    fnt.setBold(True)
                    hi.setFont(fnt)
        layout.addWidget(self.table)
        for i, row in enumerate(rows):
            self._set_row(i, row)

        bar = QHBoxLayout()
        add = QPushButton("Add row")
        add.clicked.connect(self._add_row)
        rem = QPushButton("Remove selected")
        rem.clicked.connect(self._remove_rows)
        bar.addWidget(add)
        bar.addWidget(rem)
        bar.addStretch(1)
        layout.addLayout(bar)

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
            if f.kind == "choice":
                combo = QComboBox()
                combo.addItems(f.choices)
                if str(val) in f.choices:
                    combo.setCurrentText(str(val))
                self.table.setCellWidget(i, j, combo)
            else:
                self.table.setItem(i, j, QTableWidgetItem("" if val is None else str(val)))

    def _add_row(self):
        i = self.table.rowCount()
        self.table.insertRow(i)
        base = self._new_row()
        self._bases.append(base)
        self._set_row(i, base)

    def _remove_rows(self):
        for r in sorted({idx.row() for idx in self.table.selectedIndexes()}, reverse=True):
            self.table.removeRow(r)
            if r < len(self._bases):
                self._bases.pop(r)

    def result_rows(self):
        out = []
        for i in range(self.table.rowCount()):
            base = dict(self._bases[i]) if i < len(self._bases) else self._new_row()
            for j, f in enumerate(self._fields):
                widget = self.table.cellWidget(i, j)
                if isinstance(widget, QComboBox):
                    base[f.key] = widget.currentText()
                else:
                    item = self.table.item(i, j)
                    base[f.key] = f.from_text(item.text() if item else "")
            out.append(base)
        return out


def _help_label(text):
    lbl = QLabel(text)
    lbl.setWordWrap(True)
    return lbl


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


class TableEditorDialog(QDialog):
    """Editable table over a list of dict records.

    ``usage_toggles`` (a list of analysis tags, e.g. ``["lem", "fem"]``) adds a row
    of checkboxes that show/hide the columns specific to each analysis, so the user
    sees only the inputs relevant to what they're doing. The toggle state persists
    per dialog. When omitted, a static color legend is shown instead."""

    def __init__(self, title, fields, rows, new_row, parent=None, help_text=None,
                 usage_toggles=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._title = title
        self.resize(min(1200, 160 + 110 * len(fields)), 460)
        layout = QVBoxLayout(self)
        if help_text:
            layout.addWidget(_help_label(help_text))
        self._editable = _EditableTable(fields, rows, new_row)
        if usage_toggles:
            layout.addLayout(self._build_toggle_bar(usage_toggles))
        else:
            legend = _usage_legend(fields)
            if legend:
                layout.addWidget(legend)
        layout.addWidget(self._editable)
        _ok_cancel(self, layout)
        if usage_toggles:
            self._apply_toggles()      # set initial column visibility

    def _build_toggle_bar(self, tags):
        from PySide6.QtWidgets import QCheckBox
        from PySide6.QtCore import QSettings
        s = QSettings("XSlope", "XSlope Studio")
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Show columns for:"))
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
    """A tabbed dialog, one editable table per tab (e.g. the two piezo lines)."""

    def __init__(self, title, tabs, parent=None, help_text=None):
        # tabs: list of (tab_title, fields, rows, new_row)
        super().__init__(parent)
        self.setWindowTitle(title)
        max_cols = max(len(fields) for _, fields, _, _ in tabs)
        self.resize(min(1200, 200 + 110 * max_cols), 480)
        layout = QVBoxLayout(self)
        if help_text:
            layout.addWidget(_help_label(help_text))
        self._tabs = QTabWidget()
        self._editables = []
        for tab_title, fields, rows, new_row in tabs:
            et = _EditableTable(fields, rows, new_row)
            self._tabs.addTab(et, tab_title)
            self._editables.append(et)
        layout.addWidget(self._tabs)
        _ok_cancel(self, layout)

    def result_rows(self, index):
        return self._editables[index].result_rows()


class _BlockListWidget(QWidget):
    """Master/detail over a list of blocks, where each block is a list of dict
    rows: a list of blocks (left) + the selected block's row table (right).
    Used for distributed loads (and reusable for other block-structured inputs)."""

    def __init__(self, blocks, fields, new_row, block_label="Load", parent=None):
        super().__init__(parent)
        self._fields = fields
        self._new_row = new_row
        self._block_label = block_label
        self._blocks = [[dict(r) for r in blk] for blk in (blocks or [])]
        self._cur = -1
        self.table = None

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
        self._holder = QVBoxLayout()
        body.addLayout(self._holder, 1)

        self._refresh_list()
        if self._blocks:
            self.list.setCurrentRow(0)

    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i in range(len(self._blocks)):
            self.list.addItem(f"{self._block_label} {i + 1}")
        self.list.blockSignals(False)

    def _commit_current(self):
        if 0 <= self._cur < len(self._blocks) and self.table is not None:
            self._blocks[self._cur] = self.table.result_rows()

    def _load(self, idx):
        if self.table is not None:
            self.table.setParent(None)
            self.table = None
        if not (0 <= idx < len(self._blocks)):
            return
        self.table = _EditableTable(self._fields, self._blocks[idx], self._new_row)
        self._holder.addWidget(self.table)

    def _on_select(self, idx):
        self._commit_current()
        self._cur = idx
        self._load(idx)

    def _add_block(self):
        self._commit_current()
        self._blocks.append([])
        self._refresh_list()
        self.list.setCurrentRow(len(self._blocks) - 1)

    def _remove_block(self):
        idx = self.list.currentRow()
        if idx < 0:
            return
        self._blocks.pop(idx)
        self._cur = -1
        self._refresh_list()
        if self._blocks:
            self.list.setCurrentRow(min(idx, len(self._blocks) - 1))
        else:
            self._load(-1)

    def result_blocks(self):
        self._commit_current()
        return [list(b) for b in self._blocks]


# --------------------------------------------------------------------------- #
# Per-category editors
# --------------------------------------------------------------------------- #
class CategoryEditor:
    label = ""

    def build(self, slope_data, parent):
        raise NotImplementedError

    def apply(self, slope_data, dlg):
        raise NotImplementedError


class GlobalEditor(CategoryEditor):
    label = "Global parameters"
    FIELDS = [
        Field("gamma_water", "Unit weight of water"),
        Field("tcrack_depth", "Tension crack depth"),
        Field("tcrack_water", "Water in crack"),
        Field("k_seismic", "Seismic coefficient k"),
    ]

    def build(self, slope_data, parent):
        return FormEditorDialog("Global parameters", self.FIELDS, slope_data, parent)

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
    return {"name": "", "gamma": 0.0, "option": "mc", "c": 0.0, "phi": 0.0,
            "cp": 0.0, "r_elev": 0.0, "d": 0.0, "psi": 0.0, "u": "none",
            "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0, "sigma_cp": 0.0,
            "sigma_d": 0.0, "sigma_psi": 0.0, "k1": 0.0, "k2": 0.0, "alpha": 0.0,
            "kr0": 0.0, "h0": 0.0, "E": 0.0, "nu": 0.0}


class MaterialsEditor(CategoryEditor):
    label = "Materials"
    # Columns mirror the 'mat' worksheet in order: name, g, option, c, f, c/p,
    # r-elev, d, psi, u, s(g), s(c), s(f), s(c/p), s(d), s(psi), k1, k2, alpha,
    # kr0, h0, E, n.
    # `applies` tags mirror the template's analysis usage (input_template.md):
    # strength (g, option, c, f, c/p, r-elev, u) is shared by LEM+FEM; d/psi are
    # rapid-drawdown (LEM); s(...) are reliability; k1..h0 seepage; E/n FEM.
    LF = {"lem", "fem"}
    FIELDS = [
        Field("name", "name", "str"),
        Field("gamma", "g", applies=LF),
        Field("option", "option", "choice", choices=["mc", "cp"], applies=LF),
        Field("c", "c", applies=LF), Field("phi", "f", applies=LF),
        Field("cp", "c/p", applies=LF), Field("r_elev", "r-elev", applies=LF),
        Field("d", "d", usage="lem"), Field("psi", "psi", usage="lem"),
        Field("u", "u", "choice", choices=["none", "piezo", "seep"], applies=LF),
        Field("sigma_gamma", "s(g)", usage="rel"), Field("sigma_c", "s(c)", usage="rel"),
        Field("sigma_phi", "s(f)", usage="rel"), Field("sigma_cp", "s(c/p)", usage="rel"),
        Field("sigma_d", "s(d)", usage="rel"), Field("sigma_psi", "s(psi)", usage="rel"),
        Field("k1", "k1", usage="seep"), Field("k2", "k2", usage="seep"),
        Field("alpha", "alpha", usage="seep"), Field("kr0", "kr0", usage="seep"),
        Field("h0", "h0", usage="seep"),
        Field("E", "E", usage="fem"), Field("nu", "n", usage="fem"),
    ]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Materials", self.FIELDS, slope_data.get("materials", []), _new_material, parent,
            help_text="Columns mirror the 'mat' worksheet. Row order = Mat ID order "
                      "(row 1 → Mat ID 1). Use the toggles to show only the columns "
                      "for your analysis; reliability σ columns are hidden unless "
                      "'Reliability' is on.",
            usage_toggles=["lem", "seep", "fem", "rel"])

    def apply(self, slope_data, dlg):
        slope_data["materials"] = dlg.result_rows()


def _new_circle():
    return {"Xo": 0.0, "Yo": 0.0, "Option": "Depth", "Depth": 0.0,
            "Xi": 0.0, "Yi": 0.0, "R": 0.0}


class CirclesEditor(CategoryEditor):
    label = "Circles"
    # Mirror the 'circles' worksheet: Xo, Yo, Option, Depth, Xi, Yi, R.
    FIELDS = [
        Field("Xo", "Xo"), Field("Yo", "Yo"),
        Field("Option", "Option", "choice", choices=["Depth", "Radius", "Intercept"]),
        Field("Depth", "Depth"), Field("Xi", "Xi"), Field("Yi", "Yi"), Field("R", "R"),
    ]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Circles", self.FIELDS, slope_data.get("circles", []), _new_circle, parent,
            help_text="Option sets how each circle's size is defined (only the matching "
                      "field is used):\n"
                      "  • Depth — elevation of the circle bottom at the center (R = Yo − Depth)\n"
                      "  • Radius — the circle radius R directly (Depth = Yo − R)\n"
                      "  • Intercept — a point (Xi, Yi) the circle passes through")

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


def _new_ncpt():
    return {"X": 0.0, "Y": 0.0, "Movement": "Free"}


class NonCircEditor(CategoryEditor):
    label = "Non-circular surface"
    FIELDS = [Field("X", "X"), Field("Y", "Y"),
              Field("Movement", "Movement", "choice", choices=["Free", "Horiz", "Fixed"])]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Non-circular surface", self.FIELDS, slope_data.get("non_circ", []), _new_ncpt, parent,
            help_text="Points ordered left→right. Entry/exit points use Movement='Free'.")

    def apply(self, slope_data, dlg):
        slope_data["non_circ"] = dlg.result_rows()


def _new_pt():
    return {"x": 0.0, "y": 0.0}


class PiezoEditor(CategoryEditor):
    label = "Piezometric lines"
    FIELDS = [Field("x", "x"), Field("y", "y")]

    def _rows(self, slope_data, key):
        return [{"x": x, "y": y} for (x, y) in (slope_data.get(key) or [])]

    def build(self, slope_data, parent):
        return TabbedTableEditorDialog(
            "Piezometric lines",
            [("Line 1", self.FIELDS, self._rows(slope_data, "piezo_line"), _new_pt),
             ("Line 2 (rapid drawdown)", self.FIELDS, self._rows(slope_data, "piezo_line2"), _new_pt)],
            parent,
            help_text="Points ordered left→right. Line 2 is only used for rapid-drawdown "
                      "analysis (the drawdown / second water table).")

    def apply(self, slope_data, dlg):
        slope_data["piezo_line"] = [(r["x"], r["y"]) for r in dlg.result_rows(0)]
        slope_data["piezo_line2"] = [(r["x"], r["y"]) for r in dlg.result_rows(1)]


# --- distributed loads (two sets; each a list of point blocks) -------------- #
DLOAD_FIELDS = [Field("X", "X"), Field("Y", "Y"), Field("Normal", "Normal")]


def _new_dload_pt():
    return {"X": 0.0, "Y": 0.0, "Normal": 0.0}


class DloadsEditor(CategoryEditor):
    label = "Distributed loads"

    def build(self, slope_data, parent):
        dlg = QDialog(parent)
        dlg.setWindowTitle("Distributed loads")
        dlg.resize(640, 520)
        layout = QVBoxLayout(dlg)
        layout.addWidget(_help_label(
            "Each load is a left→right series of points (≥2). Select a load to edit its "
            "points. Set 2 is the second rapid-drawdown stage."))
        tabs = QTabWidget()
        w1 = _BlockListWidget(slope_data.get("dloads"), DLOAD_FIELDS, _new_dload_pt, "Load")
        w2 = _BlockListWidget(slope_data.get("dloads2"), DLOAD_FIELDS, _new_dload_pt, "Load")
        tabs.addTab(w1, "Set 1")
        tabs.addTab(w2, "Set 2 (rapid drawdown)")
        layout.addWidget(tabs)
        _ok_cancel(dlg, layout)
        dlg._sets = (w1, w2)
        return dlg

    def apply(self, slope_data, dlg):
        w1, w2 = dlg._sets
        slope_data["dloads"] = w1.result_blocks()
        slope_data["dloads2"] = w2.result_blocks()


# --- seep BC (two sets; each: a list of specified-head BCs + an exit face) --- #
XY_FIELDS = [Field("x", "x"), Field("y", "y")]


class _SeepBcSetWidget(QWidget):
    """One seepage BC set as a single master/detail: the left list holds each
    specified-head boundary AND the exit face; the right side shows one
    full-height point table (plus a head-value field, hidden for the exit face)."""

    def __init__(self, bc, parent=None):
        super().__init__(parent)
        bc = bc or {}
        self._heads = [{"head": h.get("head", 0.0),
                        "coords": [tuple(c) for c in h.get("coords", [])]}
                       for h in (bc.get("specified_heads") or [])]
        self._exit = [tuple(c) for c in (bc.get("exit_face") or [])]
        self._cur = -1
        self.table = None

        body = QHBoxLayout(self)
        body.setContentsMargins(0, 0, 0, 0)
        left = QVBoxLayout()
        body.addLayout(left)
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

        right = QVBoxLayout()
        body.addLayout(right, 1)
        hrow = QHBoxLayout()
        self.head_label = QLabel("Head value:")
        hrow.addWidget(self.head_label)
        self.head_edit = QLineEdit()
        hrow.addWidget(self.head_edit, 1)
        right.addLayout(hrow)
        self._holder = QVBoxLayout()
        right.addLayout(self._holder, 1)

        self._refresh()
        self.list.setCurrentRow(0)

    def _is_exit(self, idx):
        return idx == len(self._heads)

    def _refresh(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i, h in enumerate(self._heads):
            self.list.addItem(f"Head {i + 1}  (h = {h['head']})")
        self.list.addItem("Exit face")
        self.list.blockSignals(False)

    def _commit(self):
        if self._cur < 0:
            return
        coords = [(r["x"], r["y"]) for r in self.table.result_rows()] if self.table else []
        if self._is_exit(self._cur):
            self._exit = coords
        elif self._cur < len(self._heads):
            try:
                self._heads[self._cur]["head"] = float(self.head_edit.text() or 0)
            except ValueError:
                self._heads[self._cur]["head"] = 0.0
            self._heads[self._cur]["coords"] = coords

    def _load(self, idx):
        if self.table is not None:
            self.table.setParent(None)
            self.table = None
        if idx < 0:
            return
        is_exit = self._is_exit(idx)
        self.head_label.setVisible(not is_exit)
        self.head_edit.setVisible(not is_exit)
        if is_exit:
            rows = [{"x": x, "y": y} for (x, y) in self._exit]
        else:
            self.head_edit.setText(str(self._heads[idx]["head"]))
            rows = [{"x": x, "y": y} for (x, y) in self._heads[idx]["coords"]]
        self.table = _EditableTable(XY_FIELDS, rows, _new_pt)
        self._holder.addWidget(self.table)

    def _on_select(self, idx):
        self._commit()
        self._cur = idx
        self._load(idx)

    def _add_head(self):
        self._commit()
        self._heads.append({"head": 0.0, "coords": []})
        self._refresh()
        self.list.setCurrentRow(len(self._heads) - 1)

    def _remove_head(self):
        idx = self.list.currentRow()
        if idx < 0 or self._is_exit(idx):  # the exit face is permanent
            return
        self._heads.pop(idx)
        self._cur = -1
        self._refresh()
        self.list.setCurrentRow(min(idx, len(self._heads)))

    def result(self):
        self._commit()
        return {"specified_heads": [{"head": h["head"], "coords": list(h["coords"])}
                                    for h in self._heads],
                "exit_face": list(self._exit)}


class SeepBcEditor(CategoryEditor):
    label = "Seep BC"

    def build(self, slope_data, parent):
        dlg = QDialog(parent)
        dlg.setWindowTitle("Seep BC")
        dlg.resize(640, 520)
        layout = QVBoxLayout(dlg)
        layout.addWidget(_help_label(
            "Seepage boundary conditions: a list of specified-head boundaries (each a head "
            "value + its points) and an exit face (where water leaves the slope). Set 2 is "
            "used for rapid-drawdown (the second seepage solution)."))
        tabs = QTabWidget()
        w1 = _SeepBcSetWidget(slope_data.get("seepage_bc"))
        w2 = _SeepBcSetWidget(slope_data.get("seepage_bc2"))
        tabs.addTab(w1, "Set 1")
        tabs.addTab(w2, "Set 2 (rapid drawdown)")
        layout.addWidget(tabs)
        _ok_cancel(dlg, layout)
        dlg._sets = (w1, w2)
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
            "area": None, "V_cap": None, "M_cap": None, "fixity": "free"}


class PilesEditor(CategoryEditor):
    label = "Piles"
    FIELDS = [
        Field("label", "Label", "str"),
        Field("x1", "x1"), Field("y1", "y1"), Field("x2", "x2"), Field("y2", "y2"),
        Field("H", "H", "optfloat"),
        Field("D_pile", "D", "optfloat", usage="lem"), Field("S", "S", "optfloat", usage="lem"),
        Field("E", "E", "optfloat", usage="fem"), Field("I", "I", "optfloat", usage="fem"),
        Field("area", "Area", "optfloat", usage="fem"),
        Field("V_cap", "Vcap", "optfloat", usage="lem"), Field("M_cap", "Mcap", "optfloat", usage="lem"),
        Field("fixity", "Fixity", "choice", choices=["free", "fixed"], usage="fem"),
    ]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Piles", self.FIELDS, slope_data.get("pile_lines", []), _new_pile, parent,
            help_text="Leave H blank for auto Ito & Matsui force. I / Area auto-compute "
                      "from D when blank. θ is auto-derived from the pile axis. Vcap/Mcap "
                      "require S (spacing).",
            usage_toggles=["lem", "fem"])

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
        return {"polygon": poly, "mat_id": p.get("mat_id")}
    return {"polygon": p, "mat_id": None}


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
        polys = [{"polygon": Polygon(p["coords"]), "mat_id": p["mat_id"]}
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
                 select=None, max_depth=None):
        # items: list of {"mat_id": int|None, "coords": [(x, y), ...]}
        # select: row to pre-highlight (e.g. the double-clicked line); else first.
        # max_depth: when not None, show a "Max depth" field (profile sheet only —
        #   it has no meaning for polygon input); the bottom boundary elevation
        #   used when building zone polygons from the profile lines.
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(680, 540)
        self._item_label = item_label
        self._materials = materials
        self._lines = [{"mat_id": it.get("mat_id"),
                        "coords": [tuple(c) for c in it.get("coords", [])]}
                       for it in (items or [])]
        self._cur = -1
        self.table = None

        main = QVBoxLayout(self)
        main.addWidget(_help_label(help_text))

        self._max_depth_edit = None
        if max_depth is not None:
            mdrow = QHBoxLayout()
            mdrow.addWidget(QLabel("Max depth (bottom boundary elevation):"))
            self._max_depth_edit = QLineEdit(str(max_depth))
            self._max_depth_edit.setToolTip(
                "Elevation of the model's bottom boundary, used to build the zone "
                "polygons from the profile lines.")
            mdrow.addWidget(self._max_depth_edit, 1)
            main.addLayout(mdrow)

        body = QHBoxLayout()
        main.addLayout(body, 1)

        left = QVBoxLayout()
        body.addLayout(left)
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

        right = QVBoxLayout()
        body.addLayout(right, 1)
        matrow = QHBoxLayout()
        matrow.addWidget(QLabel("Material:"))
        self.mat_combo = QComboBox()
        for i, m in enumerate(materials):
            self.mat_combo.addItem(f"{i + 1}: {m.get('name', '')}")
        self.mat_combo.currentIndexChanged.connect(self._on_mat_changed)
        matrow.addWidget(self.mat_combo, 1)
        right.addLayout(matrow)
        self._holder = QVBoxLayout()
        right.addLayout(self._holder, 1)

        _ok_cancel(self, main)

        self._refresh_list()
        if self._lines:
            row = select if (select is not None and 0 <= select < len(self._lines)) else 0
            self.list.setCurrentRow(row)

    def _label(self, i):
        mid = self._lines[i]["mat_id"]
        if mid is not None and 0 <= mid < len(self._materials):
            return f"{self._item_label} {i + 1}  (mat {mid + 1}: {self._materials[mid].get('name', '')})"
        return f"{self._item_label} {i + 1}  (mat ?)"

    def _refresh_list(self):
        self.list.blockSignals(True)
        self.list.clear()
        for i in range(len(self._lines)):
            self.list.addItem(self._label(i))
        self.list.blockSignals(False)

    def _commit_current(self):
        if 0 <= self._cur < len(self._lines) and self.table is not None:
            self._lines[self._cur]["coords"] = [(r["x"], r["y"]) for r in self.table.result_rows()]
            if self.mat_combo.currentIndex() >= 0:
                self._lines[self._cur]["mat_id"] = self.mat_combo.currentIndex()

    def _load(self, idx):
        if self.table is not None:
            self.table.setParent(None)
            self.table = None
        if not (0 <= idx < len(self._lines)):
            return
        ln = self._lines[idx]
        rows = [{"x": x, "y": y} for (x, y) in ln["coords"]]
        self.table = _EditableTable(self.XY, rows, _new_pt)
        self._holder.addWidget(self.table)
        self.mat_combo.blockSignals(True)
        mid = ln["mat_id"]
        self.mat_combo.setCurrentIndex(mid if (mid is not None and 0 <= mid < self.mat_combo.count()) else 0)
        self.mat_combo.blockSignals(False)

    def _on_select(self, new_idx):
        self._commit_current()
        self._cur = new_idx
        self._load(new_idx)

    def _on_mat_changed(self, idx):
        if 0 <= self._cur < len(self._lines) and idx >= 0:
            self._lines[self._cur]["mat_id"] = idx
            item = self.list.item(self._cur)
            if item:
                item.setText(self._label(self._cur))

    def _add_line(self):
        self._commit_current()
        self._lines.append({"mat_id": 0, "coords": []})
        self._refresh_list()
        self.list.setCurrentRow(len(self._lines) - 1)

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

    def result_lines(self):
        self._commit_current()
        return [{"coords": list(ln["coords"]), "mat_id": ln["mat_id"]} for ln in self._lines]

    def result_max_depth(self):
        """The edited max-depth value (float), or None if the field isn't shown."""
        if self._max_depth_edit is None:
            return None
        try:
            return float(self._max_depth_edit.text())
        except (TypeError, ValueError):
            return None


class ProfileEditor(CategoryEditor):
    label = "Profile lines"

    def build(self, slope_data, parent, select=None):
        return MatGeometryDialog(
            "Profile lines",
            "Each profile line is the top of a material layer, drawn left→right and "
            "ordered shallowest first. Select a line to edit its material and vertices.",
            "Line",
            slope_data.get("profile_lines") or [],
            slope_data.get("materials") or [], parent,
            select=select, max_depth=slope_data.get("max_depth") or 0.0)

    def apply(self, slope_data, dlg):
        slope_data["profile_lines"] = dlg.result_lines()
        md = dlg.result_max_depth()
        if md is not None:
            slope_data["max_depth"] = md
        _resync_geometry(slope_data)  # rebuild polygons / ground surface / t-crack


class PolygonEditor(CategoryEditor):
    """Geometry editor for polygon-based files (no profile sheet). Each polygon is
    a closed material zone; the loader closes the ring implicitly, so the editor
    shows/stores only the distinct exterior vertices."""

    label = "Polygons"

    def build(self, slope_data, parent, select=None):
        items = []
        for p in (slope_data.get("polygons") or []):
            coords = list(p["polygon"].exterior.coords)
            if len(coords) >= 2 and coords[0] == coords[-1]:
                coords = coords[:-1]                       # drop the closing duplicate
            items.append({"mat_id": p.get("mat_id"), "coords": coords})
        return MatGeometryDialog(
            "Polygons",
            "Each polygon is a closed material zone (the ring is closed automatically, "
            "so list each vertex once). Select a polygon to edit its material and vertices.",
            "Polygon", items, slope_data.get("materials") or [], parent, select=select)

    def apply(self, slope_data, dlg):
        from shapely.geometry import Polygon
        polys = []
        for it in dlg.result_lines():
            coords = it["coords"]
            if len(coords) < 3:
                continue                                   # not a valid ring; skip
            polys.append({"polygon": Polygon(coords), "mat_id": it["mat_id"]})
        _set_derived_geometry(slope_data, polys)


# --- reinforcement ---------------------------------------------------------- #
def _new_reinf():
    return {"x1": 0.0, "y1": 0.0, "x2": 0.0, "y2": 0.0, "t_max": 0.0, "t_res": 0.0,
            "lp1": 0.0, "lp2": 0.0, "E": 0.0, "area": 0.0}


class ReinforcementEditor(CategoryEditor):
    label = "Reinforcement"
    FIELDS = [
        Field("x1", "x1"), Field("y1", "y1"), Field("x2", "x2"), Field("y2", "y2"),
        Field("t_max", "Tmax", usage="lem"), Field("t_res", "Tres", usage="fem"),
        Field("lp1", "Lp1", usage="lem"), Field("lp2", "Lp2", usage="lem"),
        Field("E", "E", usage="fem"), Field("area", "Area", usage="fem"),
    ]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Reinforcement", self.FIELDS, slope_data.get("reinforcement_lines", []),
            _new_reinf, parent,
            help_text="Lp1/Lp2 are the pullout lengths at each end (0 = fully anchored). "
                      "The LEM tension distribution shown on the plot is derived from these.",
            usage_toggles=["lem", "fem"])

    def apply(self, slope_data, dlg):
        from xslope.fileio import build_reinforce_lines
        rows = dlg.result_rows()
        slope_data["reinforcement_lines"] = rows
        # Rebuild the LEM display/analysis format so the canvas reflects the edit.
        slope_data["reinforce_lines"] = build_reinforce_lines(rows)


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
    "profile": ProfileEditor(),
    "polygons": PolygonEditor(),
}
