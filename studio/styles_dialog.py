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
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QColorDialog, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QGridLayout, QHeaderView, QLabel, QMenu, QTableWidget, QTableWidgetItem,
    QToolButton, QVBoxLayout, QWidget, QWidgetAction,
)

from xslope.plot import get_material_color

# Curated preset palette for the color popup — earth/soil tones, greens, blues,
# accents, neutrals (6 columns). "Custom…" falls back to the OS color dialog.
COLOR_PALETTE = [
    ("saddlebrown", "#8b4513"), ("sienna", "#a0522d"), ("peru", "#cd853f"),
    ("tan", "#d2b48c"), ("burlywood", "#deb887"), ("wheat", "#f5deb3"),
    ("goldenrod", "#daa520"), ("darkkhaki", "#bdb76b"), ("olive", "#808000"),
    ("darkolivegreen", "#556b2f"), ("khaki", "#f0e68c"), ("gold", "#ffd700"),
    ("forestgreen", "#228b22"), ("seagreen", "#2e8b57"), ("yellowgreen", "#9acd32"),
    ("darkseagreen", "#8fbc8f"), ("mediumseagreen", "#3cb371"), ("lightgreen", "#90ee90"),
    ("steelblue", "#4682b4"), ("cornflowerblue", "#6495ed"), ("cadetblue", "#5f9ea0"),
    ("lightskyblue", "#87cefa"), ("lightblue", "#add8e6"), ("navy", "#000080"),
    ("firebrick", "#b22222"), ("indianred", "#cd5c5c"), ("salmon", "#fa8072"),
    ("orange", "#ffa500"), ("plum", "#dda0dd"), ("teal", "#008080"),
    ("black", "#000000"), ("dimgray", "#696969"), ("gray", "#808080"),
    ("darkgray", "#a9a9a9"), ("lightgray", "#d3d3d3"), ("white", "#ffffff"),
]


class ColorButton(QToolButton):
    """A swatch button showing the current color; clicking pops up a preset palette
    grid plus a "Custom…" entry that opens the OS color dialog. Emits ``colorChanged``
    (hex) on any change."""

    colorChanged = Signal(str)

    def __init__(self, hex_color, parent=None):
        super().__init__(parent)
        self._hex = hex_color
        self.setFixedSize(60, 22)
        self.setCursor(Qt.PointingHandCursor)
        self._apply()
        self.clicked.connect(self._popup)

    @property
    def hex(self):
        return self._hex

    def set_hex(self, h):
        self._hex = h
        self._apply()
        self.colorChanged.emit(h)

    def _apply(self):
        c = QColor(self._hex)
        fg = "#000" if c.lightnessF() > 0.5 else "#fff"
        self.setText(self._hex)
        self.setStyleSheet(f"QToolButton{{background:{self._hex}; color:{fg}; "
                           f"border:1px solid #888; font-size:9px;}}")

    def _popup(self):
        menu = QMenu(self)
        grid = QWidget()
        g = QGridLayout(grid)
        g.setSpacing(2)
        g.setContentsMargins(6, 6, 6, 6)
        for i, (name, hexc) in enumerate(COLOR_PALETTE):
            sw = QToolButton()
            sw.setFixedSize(20, 20)
            sw.setToolTip(name)
            sw.setCursor(Qt.PointingHandCursor)
            sw.setStyleSheet(f"QToolButton{{background:{hexc}; border:1px solid #aaa;}}")
            sw.clicked.connect(lambda _=False, h=hexc, m=menu: (self.set_hex(h), m.close()))
            g.addWidget(sw, i // 6, i % 6)
        wa = QWidgetAction(menu)
        wa.setDefaultWidget(grid)
        menu.addAction(wa)
        menu.addSeparator()
        menu.addAction("Custom…").triggered.connect(self._custom)
        menu.exec(self.mapToGlobal(self.rect().bottomLeft()))

    def _custom(self):
        col = QColorDialog.getColor(QColor(self._hex), self, "Custom color")
        if col.isValid():
            self.set_hex(col.name())

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
        for c in (0, 1, 3):
            hh.setSectionResizeMode(c, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.Stretch)   # Hatch column takes the slack

        self._rows = []   # (idx, color_btn, hatch_combo, alpha_spin)
        for idx, mat in enumerate(self._materials):
            ov = mat_overrides.get(str(idx)) or {}
            name = mat.get("name") if isinstance(mat, dict) else str(mat)
            item = QTableWidgetItem(name or f"Material {idx + 1}")
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(idx, 0, item)

            color = ov.get("color") or to_hex(get_material_color(idx))
            btn = ColorButton(color)
            btn.colorChanged.connect(self._emit)
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

    # --- result / preview ------------------------------------------------
    def result(self):
        """The sparse style-delta dict: only materials that differ from defaults
        (palette color, no hatch, 0.6 opacity). Non-material keys are preserved."""
        out = copy.deepcopy(self._orig)
        mats = {}
        for idx, btn, combo, spin in self._rows:
            entry = {}
            if btn.hex.lower() != to_hex(get_material_color(idx)).lower():
                entry["color"] = btn.hex
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
