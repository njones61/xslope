"""Per-plot-type display-options panels for the context-sensitive Display dock.

Display options vary by plot type and are *not* solve inputs, so each result view
has its own small panel of controls. A panel only holds widgets, exposes the
current values via ``options()``, and emits ``changed`` when any control moves;
the main window re-renders the associated view (no re-solve) in response.

Every control here maps to a keyword argument the underlying ``plot_*`` function
already accepts. Styling that those functions hardcode (colors, line widths,
fonts, fills) is not exposed yet — that arrives with the Phase 5 ``StyleConfig``.
"""

from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QSpinBox, QWidget,
)

from .dialogs import SEEP_VARIABLES

# Material-table placements accepted by plot_inputs(tab_loc=…).
TAB_LOCATIONS = [
    "top", "upper left", "upper right", "upper center", "lower left",
    "lower right", "lower center", "center left", "center right", "center",
]

# FEM Results plot types (one shown at a time). plot_fem_results also supports
# displace_mag / stress / strain / yield, but those are diagnostic and omitted here.
FEM_PLOT_TYPES = [
    ("shear_strain", "Shear strain"),
    ("deformation", "Deformation"),
    ("displace_vector", "Displacement vectors"),
]


def _dspin(lo, hi, val, step, decimals=2, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(decimals)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    return s


def _ispin(lo, hi, val, suffix=""):
    s = QSpinBox()
    s.setRange(lo, hi)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    return s


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


class SolutionDisplayPanel(_CheckboxPanel):
    """Display options for an LEM solution plot."""
    _FIELDS = (
        ("slice_numbers", "Slice numbers", False),
        ("seep_contours", "Seepage contours", True),
    )


class SearchDisplayPanel(_CheckboxPanel):
    """Display options for an LEM auto-search plot."""
    _FIELDS = (("highlight_fs", "Highlight critical surface", True),)


class InputsDisplayPanel(QWidget):
    """Display options for the Inputs view: the material-property table (and its
    placement) plus legend column layout."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)

        self.mat_table = QCheckBox("Material property table")
        self.tab_loc = QComboBox()
        self.tab_loc.addItems(TAB_LOCATIONS)
        self.legend_auto = QCheckBox("Auto legend columns")
        self.legend_auto.setChecked(True)
        self.legend_ncol = _ispin(1, 12, 3)
        self.legend_max_cols = _ispin(1, 12, 6)
        self.legend_max_rows = _ispin(1, 12, 4)

        form.addRow("", self.mat_table)
        form.addRow("Table position", self.tab_loc)
        form.addRow("", self.legend_auto)
        form.addRow("Legend columns", self.legend_ncol)
        form.addRow("Max columns", self.legend_max_cols)
        form.addRow("Max rows", self.legend_max_rows)

        self.mat_table.toggled.connect(self._sync)
        self.legend_auto.toggled.connect(self._sync)
        for w in (self.mat_table, self.legend_auto):
            w.toggled.connect(self._emit)
        self.tab_loc.currentIndexChanged.connect(self._emit)
        for s in (self.legend_ncol, self.legend_max_cols, self.legend_max_rows):
            s.valueChanged.connect(self._emit)
        self._sync()

    def _sync(self, *_):
        self.tab_loc.setEnabled(self.mat_table.isChecked())
        auto = self.legend_auto.isChecked()
        self.legend_ncol.setEnabled(not auto)         # explicit count when manual
        self.legend_max_cols.setEnabled(auto)         # caps only apply to "auto"
        self.legend_max_rows.setEnabled(auto)

    def _emit(self, *_):
        self.changed.emit()

    def options(self):
        return {
            "mat_table": self.mat_table.isChecked(),
            "tab_loc": self.tab_loc.currentText(),
            "legend_ncol": ("auto" if self.legend_auto.isChecked()
                            else self.legend_ncol.value()),
            "legend_max_cols": self.legend_max_cols.value(),
            "legend_max_rows": self.legend_max_rows.value(),
        }


class MeshDisplayPanel(QWidget):
    """Display options for a mesh plot (no boundary conditions)."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)
        self.show_nodes = QCheckBox("Show nodes")
        self.show_nodes.setChecked(True)
        self.label_elements = QCheckBox("Element numbers")
        self.label_nodes = QCheckBox("Node numbers")
        self.pad_frac = _dspin(0.0, 1.0, 0.05, 0.01)

        for c in (self.show_nodes, self.label_elements, self.label_nodes):
            form.addRow("", c)
            c.toggled.connect(self._emit)
        form.addRow("Padding", self.pad_frac)
        self.pad_frac.valueChanged.connect(self._emit)

    def _emit(self, *_):
        self.changed.emit()

    def options(self):
        return {
            "show_nodes": self.show_nodes.isChecked(),
            "label_elements": self.label_elements.isChecked(),
            "label_nodes": self.label_nodes.isChecked(),
            "pad_frac": self.pad_frac.value(),
        }


class FeDataDisplayPanel(QWidget):
    """Display options for a seep/FEM data plot (mesh + boundary conditions).
    ``include_bc_symbol`` adds the BC symbol-size control (FEM data only — the
    seepage data plot has no such parameter)."""

    changed = Signal()

    def __init__(self, include_bc_symbol=False, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)
        self.show_bc = QCheckBox("Boundary conditions")
        self.show_bc.setChecked(True)
        self.show_nodes = QCheckBox("Show nodes")
        self.label_elements = QCheckBox("Element numbers")
        self.label_nodes = QCheckBox("Node numbers")
        for c in (self.show_bc, self.show_nodes, self.label_elements, self.label_nodes):
            form.addRow("", c)
            c.toggled.connect(self._emit)

        self.alpha = _dspin(0.0, 1.0, 0.4, 0.05)
        form.addRow("Fill opacity", self.alpha)
        self.alpha.valueChanged.connect(self._emit)

        self.bc_symbol_size = None
        if include_bc_symbol:
            self.bc_symbol_size = _dspin(0.0, 1.0, 0.03, 0.01, decimals=3)
            form.addRow("BC symbol size", self.bc_symbol_size)
            self.bc_symbol_size.valueChanged.connect(self._emit)

    def _emit(self, *_):
        self.changed.emit()

    def options(self):
        d = {
            "show_bc": self.show_bc.isChecked(),
            "show_nodes": self.show_nodes.isChecked(),
            "label_elements": self.label_elements.isChecked(),
            "label_nodes": self.label_nodes.isChecked(),
            "alpha": self.alpha.value(),
        }
        if self.bc_symbol_size is not None:
            d["bc_symbol_size"] = self.bc_symbol_size.value()
        return d


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

        self.levels = _ispin(2, 100, 20)
        self.alpha = _dspin(0.0, 1.0, 0.4, 0.05)
        self.vector_scale = _dspin(0.001, 10.0, 0.05, 0.01, decimals=3)
        self.pad_frac = _dspin(0.0, 1.0, 0.05, 0.01)

        self.flowlines = QCheckBox("Flow lines")
        self.flowlines.setChecked(True)
        self.vectors = QCheckBox("Velocity vectors")
        self.fill = QCheckBox("Filled contours")
        self.phreatic = QCheckBox("Phreatic surface")
        self.phreatic.setChecked(True)

        form.addRow("Variable", self.variable)
        form.addRow("Base material", self.base_mat)
        form.addRow("Contour levels", self.levels)
        form.addRow("Fill opacity", self.alpha)
        form.addRow("Vector scale", self.vector_scale)
        form.addRow("Padding", self.pad_frac)
        form.addRow("", self.flowlines)
        form.addRow("", self.vectors)
        form.addRow("", self.fill)
        form.addRow("", self.phreatic)

        self.variable.currentIndexChanged.connect(self._on_variable)
        self.base_mat.currentIndexChanged.connect(self._emit)
        for s in (self.levels, self.alpha, self.vector_scale, self.pad_frac):
            s.valueChanged.connect(self._emit)
        for c in (self.flowlines, self.vectors, self.fill, self.phreatic):
            c.toggled.connect(self._emit)
        self.vectors.toggled.connect(self._sync_enabled)
        self._sync_enabled()

    def _on_variable(self, *_):
        self._sync_enabled()
        self.changed.emit()

    def _emit(self, *_):
        self.changed.emit()

    def _sync_enabled(self):
        # Flow lines / base material only apply to the head flow-net; vector scale
        # only matters when vectors are drawn.
        head = self.variable.currentData() == "head"
        self.flowlines.setEnabled(head)
        self.base_mat.setEnabled(head)
        self.vector_scale.setEnabled(self.vectors.isChecked())

    def options(self):
        return {
            "variable": self.variable.currentData(),
            "base_mat": self.base_mat.currentData(),
            "levels": self.levels.value(),
            "alpha": self.alpha.value(),
            "vector_scale": self.vector_scale.value(),
            "pad_frac": self.pad_frac.value(),
            "flowlines": self.flowlines.isChecked(),
            "vectors": self.vectors.isChecked(),
            "fill_contours": self.fill.isChecked(),
            "phreatic": self.phreatic.isChecked(),
        }


class FemResultsDisplayPanel(QWidget):
    """Display options for an FEM results plot: plot type, deformation scale, and
    the mesh/reinforcement/vector toggles. The displacement-vector options enable
    only when the plot type is ``displace_vector``."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)

        self.plot_type = QComboBox()
        for key, label in FEM_PLOT_TYPES:
            self.plot_type.addItem(label, key)

        self.deform_percent = _ispin(0, 100, 15, suffix=" %")
        self.deform_percent.setToolTip("Deformed-shape exaggeration as a percent "
                                       "of slope height.")

        self.show_mesh = QCheckBox("Show mesh")
        self.show_mesh.setChecked(True)
        self.show_reinforcement = QCheckBox("Reinforcement")
        self.show_reinforcement.setChecked(True)
        self.label_elements = QCheckBox("Element numbers")

        # Displacement-vector-only controls.
        self.plot_boundary = QCheckBox("Boundary edges")
        self.plot_boundary.setChecked(True)
        self.plot_nodes = QCheckBox("Node dots")
        self.plot_elements = QCheckBox("Element edges")
        self.scale_vectors = QCheckBox("Auto-scale vectors")
        self.displacement_tolerance = _dspin(0.0, 1.0, 0.5, 0.05)
        self.displacement_tolerance.setToolTip(
            "Hide vectors below this fraction of the max displacement.")

        form.addRow("Plot type", self.plot_type)
        form.addRow("Deform", self.deform_percent)
        form.addRow("", self.show_mesh)
        form.addRow("", self.show_reinforcement)
        form.addRow("", self.label_elements)
        form.addRow("", self.plot_boundary)
        form.addRow("", self.plot_nodes)
        form.addRow("", self.plot_elements)
        form.addRow("", self.scale_vectors)
        form.addRow("Vector cutoff", self.displacement_tolerance)

        self.plot_type.currentIndexChanged.connect(self._on_plot_type)
        self.deform_percent.valueChanged.connect(self._emit)
        self.displacement_tolerance.valueChanged.connect(self._emit)
        for c in (self.show_mesh, self.show_reinforcement, self.label_elements,
                  self.plot_boundary, self.plot_nodes, self.plot_elements,
                  self.scale_vectors):
            c.toggled.connect(self._emit)
        self._sync_enabled()

    @property
    def _vector_widgets(self):
        return (self.plot_boundary, self.plot_nodes, self.plot_elements,
                self.scale_vectors, self.displacement_tolerance)

    def _on_plot_type(self, *_):
        self._sync_enabled()
        self.changed.emit()

    def _emit(self, *_):
        self.changed.emit()

    def _sync_enabled(self):
        vec = self.plot_type.currentData() == "displace_vector"
        for w in self._vector_widgets:
            w.setEnabled(vec)

    def options(self):
        return {
            "plot_type": self.plot_type.currentData(),
            "deform_percent": self.deform_percent.value(),
            "show_mesh": self.show_mesh.isChecked(),
            "show_reinforcement": self.show_reinforcement.isChecked(),
            "label_elements": self.label_elements.isChecked(),
            "plot_boundary": self.plot_boundary.isChecked(),
            "plot_nodes": self.plot_nodes.isChecked(),
            "plot_elements": self.plot_elements.isChecked(),
            "scale_vectors": self.scale_vectors.isChecked(),
            "displacement_tolerance": self.displacement_tolerance.value(),
        }
