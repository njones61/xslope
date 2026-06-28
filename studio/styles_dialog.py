"""Styles dialog — edit the project's per-feature visual styles (plan §8a).

v1 covers **material** styles (fill color / hatch / alpha), which drive the
material zones and profile/polygon colors. Other feature styles (piezo, circles,
…) get added as they are threaded through the plot functions. Edits preview live
on the canvas; OK keeps them (and marks the project dirty → written to
``{stem}_styles.json`` on Save), Cancel restores the prior style.

The result is the sparse style-delta dict the document stores: only materials the
user actually changed from the defaults appear.
"""

from __future__ import annotations

import copy

from matplotlib.colors import to_hex
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QColorDialog, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QHeaderView, QLabel, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout,
)

from xslope.plot import get_material_color

# Curated soil-ish hatches over Matplotlib's built-in set (an approximation of
# USCS symbols, not publication-grade). Label → hatch string (None = solid fill).
HATCHES = [
    ("None", None),
    ("Dots", "...."),
    ("Dense dots", "......."),
    ("Circles", "oooo"),
    ("Horizontal", "----"),
    ("Vertical", "||||"),
    ("Diagonal /", "////"),
    ("Diagonal \\", "\\\\\\\\"),
    ("Cross", "++++"),
    ("X", "xxxx"),
]
DEFAULT_ALPHA = 0.6


class StylesDialog(QDialog):
    def __init__(self, materials, style, on_preview, parent=None):
        # materials: slope_data['materials']; style: current sparse delta dict;
        # on_preview(style_deltas): called live as edits are made (and on Cancel
        # with the original) so the canvas re-renders.
        super().__init__(parent)
        self.setWindowTitle("Styles")
        self.resize(520, 420)
        self._materials = materials or []
        self._orig = copy.deepcopy(style or {})
        self._on_preview = on_preview
        mat_overrides = (style or {}).get("materials") or {}

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "Per-material fill style. Color seeds the zone and profile/polygon lines; "
            "hatch and opacity apply to filled zones. Changes preview live."))

        self.table = QTableWidget(len(self._materials), 4, self)
        self.table.setHorizontalHeaderLabels(["Material", "Color", "Hatch", "Opacity"])
        self.table.verticalHeader().setVisible(False)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)
        for c in (1, 2, 3):
            hh.setSectionResizeMode(c, QHeaderView.ResizeToContents)

        self._rows = []   # (idx, color_btn, hatch_combo, alpha_spin)
        for idx, mat in enumerate(self._materials):
            ov = mat_overrides.get(str(idx)) or {}
            name = mat.get("name") if isinstance(mat, dict) else str(mat)
            item = QTableWidgetItem(name or f"Material {idx + 1}")
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(idx, 0, item)

            color = ov.get("color") or to_hex(get_material_color(idx))
            btn = QPushButton()
            btn.setToolTip("Pick fill color")
            self._set_btn_color(btn, color)
            btn.clicked.connect(lambda _c, b=btn: self._pick_color(b))
            self.table.setCellWidget(idx, 1, btn)

            combo = QComboBox()
            for label, h in HATCHES:
                combo.addItem(label, h)
            hi = max(0, combo.findData(ov.get("hatch")))
            combo.setCurrentIndex(hi)
            combo.currentIndexChanged.connect(self._emit)
            self.table.setCellWidget(idx, 2, combo)

            spin = QDoubleSpinBox()
            spin.setRange(0.0, 1.0)
            spin.setSingleStep(0.05)
            spin.setDecimals(2)
            spin.setValue(float(ov.get("alpha", DEFAULT_ALPHA)))
            spin.valueChanged.connect(self._emit)
            self.table.setCellWidget(idx, 3, spin)

            self._rows.append((idx, btn, combo, spin))
        layout.addWidget(self.table, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    # --- color button helpers -------------------------------------------
    def _set_btn_color(self, btn, hex_color):
        btn.setText(hex_color)
        btn.setStyleSheet(f"background-color: {hex_color};")
        btn._hex = hex_color

    def _pick_color(self, btn):
        col = QColorDialog.getColor(QColor(btn._hex), self, "Fill color")
        if col.isValid():
            self._set_btn_color(btn, col.name())
            self._emit()

    # --- result / preview ------------------------------------------------
    def result(self):
        """The sparse style-delta dict: only materials that differ from defaults
        (palette color, no hatch, 0.6 opacity). Non-material keys are preserved."""
        out = copy.deepcopy(self._orig)
        mats = {}
        for idx, btn, combo, spin in self._rows:
            entry = {}
            if btn._hex.lower() != to_hex(get_material_color(idx)).lower():
                entry["color"] = btn._hex
            if combo.currentData() is not None:
                entry["hatch"] = combo.currentData()
            if abs(spin.value() - DEFAULT_ALPHA) > 1e-9:
                entry["alpha"] = round(spin.value(), 3)
            if entry:
                mats[str(idx)] = entry
        if mats:
            out["materials"] = mats
        else:
            out.pop("materials", None)
        return out

    def _emit(self, *_):
        self._on_preview(self.result())

    def reject(self):
        self._on_preview(self._orig)     # restore the prior look
        super().reject()
