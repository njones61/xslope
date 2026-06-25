"""Input editors — modal dialogs launched from the Inputs tree (Phase 2).

Each editable ``slope_data`` category has a small ``CategoryEditor`` that builds a
dialog from the current data and writes the edited values back. Most categories
reuse the generic ``TableEditorDialog`` (a column spec + add/remove rows);
scalars use ``FormEditorDialog``. Editors only expose the commonly-edited fields
of each record and **preserve any unshown keys** (e.g. material reliability
sigmas), so a round-trip through an editor never drops data.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDialog, QDialogButtonBox, QFormLayout,
    QHBoxLayout, QLabel, QLineEdit, QPushButton, QTableWidget, QTableWidgetItem,
    QVBoxLayout,
)


# --------------------------------------------------------------------------- #
# Field spec + dialogs
# --------------------------------------------------------------------------- #
class Field:
    """One editable column/field. kind ∈ {'float','int','str','choice'}."""

    _BLANK = {"float": 0.0, "int": 0, "str": "", "choice": ""}

    def __init__(self, key, header, kind="float", choices=None, default=None):
        self.key = key
        self.header = header
        self.kind = kind
        self.choices = [str(c) for c in (choices or [])]
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


class TableEditorDialog(QDialog):
    """Editable table over a list of dict records. Unshown keys are preserved."""

    def __init__(self, title, fields, rows, new_row, parent=None, help_text=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(min(1200, 160 + 110 * len(fields)), 460)
        self._fields = fields
        self._new_row = new_row
        self._bases = [dict(r) for r in rows]  # keep originals to preserve extra keys

        layout = QVBoxLayout(self)
        if help_text:
            lbl = QLabel(help_text)
            lbl.setWordWrap(True)
            layout.addWidget(lbl)

        self.table = QTableWidget(len(rows), len(fields))
        self.table.setHorizontalHeaderLabels([f.header for f in fields])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
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

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

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
                self.table.setItem(i, j, QTableWidgetItem(str(val)))

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
    FIELDS = [
        Field("name", "Name", "str"),
        Field("gamma", "γ"),
        Field("option", "Option", "choice", choices=["mc", "cp"]),
        Field("c", "c"), Field("phi", "φ"), Field("cp", "c/p"),
        Field("r_elev", "r-elev"), Field("d", "d"), Field("psi", "ψ"),
        Field("u", "u", "choice", choices=["none", "piezo", "seep"]),
        Field("E", "E"), Field("nu", "ν"),
        Field("k1", "k1"), Field("k2", "k2"), Field("alpha", "α"),
        Field("kr0", "kr0"), Field("h0", "h0"),
    ]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Materials", self.FIELDS, slope_data.get("materials", []), _new_material, parent,
            help_text="Row order = Mat ID order (row 1 → Mat ID 1). Reliability sigmas "
                      "are preserved but not shown here.")

    def apply(self, slope_data, dlg):
        slope_data["materials"] = dlg.result_rows()


def _new_circle():
    return {"Xo": 0.0, "Yo": 0.0, "Depth": 0.0, "R": 0.0}


class CirclesEditor(CategoryEditor):
    label = "Circles"
    FIELDS = [Field("Xo", "Xo"), Field("Yo", "Yo"), Field("Depth", "Depth (elev)")]

    def build(self, slope_data, parent):
        return TableEditorDialog(
            "Circles", self.FIELDS, slope_data.get("circles", []), _new_circle, parent,
            help_text="Depth is the elevation of the circle bottom; R = Yo − Depth.")

    def apply(self, slope_data, dlg):
        rows = dlg.result_rows()
        for c in rows:
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
    label = "Piezometric line"
    FIELDS = [Field("x", "x"), Field("y", "y")]

    def build(self, slope_data, parent):
        rows = [{"x": x, "y": y} for (x, y) in (slope_data.get("piezo_line") or [])]
        return TableEditorDialog("Piezometric line 1", self.FIELDS, rows, _new_pt, parent)

    def apply(self, slope_data, dlg):
        slope_data["piezo_line"] = [(r["x"], r["y"]) for r in dlg.result_rows()]


CATEGORY_EDITORS = {
    "global": GlobalEditor(),
    "materials": MaterialsEditor(),
    "circles": CirclesEditor(),
    "non_circ": NonCircEditor(),
    "piezo": PiezoEditor(),
}
