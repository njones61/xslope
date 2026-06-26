"""Per-plot-type display-options panels for the context-sensitive Display dock.

Display options vary by plot type and are *not* solve inputs, so each result view
has its own small panel of controls. The panel only holds widgets, exposes the
current values via ``options()``, and emits ``changed`` when any control moves;
the main window re-renders the associated view (no re-solve) in response.
"""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFormLayout, QSpinBox, QWidget,
)

from .dialogs import SEEP_VARIABLES

FEM_PLOT_TYPES = [
    ("shear_strain", "Shear strain"),
    ("deformation", "Deformation"),
    ("displace_vector", "Displacement vectors"),
]


class _CheckboxPanel(QWidget):
    """Base for panels that are just a column of checkboxes. Subclasses list
    ``_FIELDS`` as ``(key, label, default)`` and read values via ``options()``."""

    changed = Signal()
    _FIELDS = ()

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)
        self._boxes = {}
        for key, label, default in self._FIELDS:
            box = QCheckBox(label)
            box.setChecked(default)
            box.toggled.connect(self._emit)
            form.addRow("", box)
            self._boxes[key] = box

    def _emit(self, *_):
        self.changed.emit()

    def options(self):
        return {key: box.isChecked() for key, box in self._boxes.items()}


class InputsDisplayPanel(_CheckboxPanel):
    """Display options for the Inputs view."""
    _FIELDS = (("mat_table", "Material property table", False),)


class SolutionDisplayPanel(_CheckboxPanel):
    """Display options for an LEM solution plot."""
    _FIELDS = (
        ("slice_numbers", "Slice numbers", False),
        ("seep_contours", "Seepage contours", True),
    )


class SearchDisplayPanel(_CheckboxPanel):
    """Display options for an LEM auto-search plot."""
    _FIELDS = (("highlight_fs", "Highlight critical surface", True),)


class MeshDisplayPanel(_CheckboxPanel):
    """Display options for a mesh plot (no boundary conditions)."""
    _FIELDS = (
        ("show_nodes", "Show nodes", True),
        ("label_elements", "Element numbers", False),
        ("label_nodes", "Node numbers", False),
    )


class FeDataDisplayPanel(_CheckboxPanel):
    """Display options for a seep/FEM data plot (mesh + boundary conditions)."""
    _FIELDS = (
        ("show_bc", "Boundary conditions", True),
        ("show_nodes", "Show nodes", False),
        ("label_elements", "Element numbers", False),
        ("label_nodes", "Node numbers", False),
    )


class SeepDisplayPanel(QWidget):
    """Display options for a seepage solution plot."""

    changed = Signal()

    def __init__(self, materials, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)

        self.variable = QComboBox()
        for key, label in SEEP_VARIABLES:
            self.variable.addItem(label, key)

        self.base_mat = QComboBox()
        for i, m in enumerate(materials or []):
            self.base_mat.addItem(f"{i + 1}: {m.get('name', '')}", i + 1)
        if self.base_mat.count() == 0:
            self.base_mat.addItem("1", 1)
        self.base_mat.setToolTip("Reference material for the flow-net "
                                 "(base permeability / number of flow channels).")

        self.levels = QSpinBox()
        self.levels.setRange(2, 100)
        self.levels.setValue(20)

        self.flowlines = QCheckBox("Flow lines")
        self.flowlines.setChecked(True)
        self.vectors = QCheckBox("Velocity vectors")
        self.fill = QCheckBox("Filled contours")
        self.phreatic = QCheckBox("Phreatic surface")
        self.phreatic.setChecked(True)

        form.addRow("Variable", self.variable)
        form.addRow("Base material", self.base_mat)
        form.addRow("Contour levels", self.levels)
        form.addRow("", self.flowlines)
        form.addRow("", self.vectors)
        form.addRow("", self.fill)
        form.addRow("", self.phreatic)

        self.variable.currentIndexChanged.connect(self._on_variable)
        self.base_mat.currentIndexChanged.connect(self._emit)
        self.levels.valueChanged.connect(self._emit)
        for c in (self.flowlines, self.vectors, self.fill, self.phreatic):
            c.toggled.connect(self._emit)
        self._sync_enabled()

    def _on_variable(self, *_):
        self._sync_enabled()
        self.changed.emit()

    def _emit(self, *_):
        self.changed.emit()

    def _sync_enabled(self):
        # Flow lines / base material only apply to the head flow-net.
        head = self.variable.currentData() == "head"
        self.flowlines.setEnabled(head)
        self.base_mat.setEnabled(head)

    def options(self):
        return {
            "variable": self.variable.currentData(),
            "base_mat": self.base_mat.currentData(),
            "levels": self.levels.value(),
            "flowlines": self.flowlines.isChecked(),
            "vectors": self.vectors.isChecked(),
            "fill_contours": self.fill.isChecked(),
            "phreatic": self.phreatic.isChecked(),
        }


class FemResultsDisplayPanel(QWidget):
    """Display options for an FEM results plot (plot type + deformation scale)."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)

        self.plot_type = QComboBox()
        for key, label in FEM_PLOT_TYPES:
            self.plot_type.addItem(label, key)

        self.deform_percent = QSpinBox()
        self.deform_percent.setRange(0, 100)
        self.deform_percent.setValue(15)
        self.deform_percent.setToolTip("Deformed-shape exaggeration as a percent "
                                       "of slope height.")

        form.addRow("Plot type", self.plot_type)
        form.addRow("Deform %", self.deform_percent)

        self.plot_type.currentIndexChanged.connect(self._emit)
        self.deform_percent.valueChanged.connect(self._emit)

    def _emit(self, *_):
        self.changed.emit()

    def options(self):
        return {
            "plot_type": self.plot_type.currentData(),
            "deform_percent": self.deform_percent.value(),
        }
