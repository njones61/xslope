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

import numpy as np
from matplotlib import colormaps as _mpl_colormaps
from PySide6.QtCore import QSize, Signal
from PySide6.QtGui import QIcon, QImage, QPixmap
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QSpinBox, QWidget,
)

from xslope.colormaps import RAMP_CHOICES
from .dialogs import SEEP_VARIABLES

# Size of the gradient preview swatch shown beside each colormap name.
_CMAP_ICON_SIZE = QSize(72, 14)

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


# Minimum width for spin boxes so the value isn't clipped by the up/down arrows
# (notably 3-decimal values like "0.050").
_SPIN_MIN_W = 78


def _dspin(lo, hi, val, step, decimals=2, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(decimals)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    s.setMinimumWidth(_SPIN_MIN_W)
    return s


def _ispin(lo, hi, val, suffix="", step=1):
    s = QSpinBox()
    s.setRange(lo, hi)
    s.setValue(val)
    s.setSingleStep(step)
    if suffix:
        s.setSuffix(suffix)
    s.setMinimumWidth(_SPIN_MIN_W)
    return s


def _add_legend_controls(panel, form):
    """Append the shared legend controls to a panel's form, grouped: a 'Legend'
    on/off toggle first, then 'Auto legend columns' + an explicit 'Legend columns'
    count + a 'Legend frame' toggle. All wired to ``panel._emit``. When Legend is
    off the layout controls disable. Stashes the widgets and exposes the values via
    ``_legend_option`` / ``_legend_frame_option`` / ``_show_legend_option``."""
    legend = QCheckBox("Legend")                 # on/off; heads the legend group
    legend.setChecked(True)
    auto = QCheckBox("Auto legend columns")
    auto.setChecked(True)
    ncol = _ispin(1, 12, 3)
    frame = QCheckBox("Legend frame")            # frameless by default
    frame.setChecked(False)
    form.addRow("", legend)
    form.addRow("", auto)
    form.addRow("Legend columns", ncol)
    form.addRow("", frame)

    def sync(*_):
        on = legend.isChecked()
        auto.setEnabled(on)
        frame.setEnabled(on)
        ncol.setEnabled(on and not auto.isChecked())

    legend.toggled.connect(sync)
    legend.toggled.connect(panel._emit)
    auto.toggled.connect(sync)
    auto.toggled.connect(panel._emit)
    ncol.valueChanged.connect(panel._emit)
    frame.toggled.connect(panel._emit)
    sync()
    panel._show_legend = legend
    panel._legend_auto = auto
    panel._legend_ncol = ncol
    panel._legend_frame = frame


def _legend_option(panel):
    """The legend_ncol value for a panel built with _add_legend_controls:
    "auto" when the Auto box is checked, else the explicit column count."""
    return "auto" if panel._legend_auto.isChecked() else panel._legend_ncol.value()


def _legend_frame_option(panel):
    """Whether the legend frame is shown (panels built with _add_legend_controls)."""
    return panel._legend_frame.isChecked()


def _add_title_toggle(panel, form):
    """Prepend a 'Title' on/off checkbox (on) to a panel's form, wired to
    ``panel._emit`` — every plot type can hide its title. Read via
    ``_show_title_option``. Also initializes ``panel._show_legend = None``; panels
    that draw a legend add their own 'Legend' toggle grouped with the legend
    controls (via ``_add_legend_controls`` or inline)."""
    title = QCheckBox("Title")
    title.setChecked(True)
    form.addRow("", title)
    title.toggled.connect(panel._emit)
    panel._show_title = title
    panel._show_legend = None


def _show_title_option(panel):
    """Whether the title is shown (panels built with _add_title_toggle)."""
    return panel._show_title.isChecked()


def _show_legend_option(panel):
    """Whether the legend is shown (True when the panel has no legend toggle)."""
    lg = getattr(panel, "_show_legend", None)
    return lg is None or lg.isChecked()


def _cmap_icon(name):
    """A horizontal gradient swatch of colormap ``name`` as a QIcon, so the
    selector previews each ramp the way the screenshot does."""
    w, h = _CMAP_ICON_SIZE.width(), _CMAP_ICON_SIZE.height()
    cmap = _mpl_colormaps[name]
    row = (cmap(np.linspace(0.0, 1.0, w)) * 255).astype(np.uint8)   # (w, 4) RGBA
    arr = np.ascontiguousarray(np.broadcast_to(row, (h, w, 4)))     # (h, w, 4)
    img = QImage(arr.data, w, h, 4 * w, QImage.Format_RGBA8888).copy()
    return QIcon(QPixmap.fromImage(img))


def _make_cmap_combo(default):
    """A QComboBox of the curated color ramps, each shown with a gradient preview,
    with ``default`` (a matplotlib colormap name) pre-selected. ``currentData()``
    returns the selected colormap name to pass through as ``cmap=…``."""
    combo = QComboBox()
    combo.setIconSize(_CMAP_ICON_SIZE)
    names = [name for name, _ in RAMP_CHOICES]
    for name, label in RAMP_CHOICES:
        combo.addItem(_cmap_icon(name), label, name)
    if default not in names:   # keep an off-list default selectable
        combo.addItem(_cmap_icon(default), default, default)
    combo.setCurrentIndex(combo.findData(default))
    return combo


class _CheckboxPanel(QWidget):
    """Base for panels that are just a column of checkboxes (plus, optionally, the
    shared legend controls). Subclasses list ``_FIELDS`` as ``(key, label,
    default)`` and set ``_WITH_LEGEND`` to add the legend column controls; values
    are read via ``options()``."""

    changed = Signal()
    _FIELDS = ()
    _WITH_LEGEND = False

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)
        _add_title_toggle(self, form)
        self._boxes = {}
        for key, label, default in self._FIELDS:
            box = QCheckBox(label)
            box.setChecked(default)
            box.toggled.connect(self._emit)
            form.addRow("", box)
            self._boxes[key] = box
        if self._WITH_LEGEND:
            _add_legend_controls(self, form)

    def _emit(self, *_):
        self.changed.emit()

    def options(self):
        d = {key: box.isChecked() for key, box in self._boxes.items()}
        d["show_title"] = _show_title_option(self)
        d["show_legend"] = _show_legend_option(self)
        if self._WITH_LEGEND:
            d["legend_ncol"] = _legend_option(self)
            d["legend_frame"] = _legend_frame_option(self)
        return d


class SolutionDisplayPanel(_CheckboxPanel):
    """Display options for an LEM solution plot."""
    _WITH_LEGEND = True
    _FIELDS = (
        ("slice_numbers", "Slice numbers", False),
        ("seep_contours", "Seepage contours", True),
    )


class SearchDisplayPanel(_CheckboxPanel):
    """Display options for an LEM auto-search plot."""
    _WITH_LEGEND = True
    _FIELDS = (("highlight_fs", "Highlight critical surface", True),)


class ReliabilityDisplayPanel(_CheckboxPanel):
    """Display options for an LEM reliability plot — just the legend controls (the
    plot itself has no other toggles)."""
    _WITH_LEGEND = True


class InputsDisplayPanel(QWidget):
    """Display options for the Inputs view: the material-property table (and its
    placement) plus legend column layout."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)
        _add_title_toggle(self, form)

        self.mat_table = QCheckBox("Material property table")
        self.tab_loc = QComboBox()
        self.tab_loc.addItems(TAB_LOCATIONS)
        form.addRow("", self.mat_table)
        form.addRow("Table position", self.tab_loc)
        self.label_coords = QCheckBox("Coordinate labels")
        self.coord_size = _ispin(4, 16, 7, suffix=" pt")
        form.addRow("", self.label_coords)
        form.addRow("Coordinate label size", self.coord_size)
        _add_legend_controls(self, form)         # 'Legend' toggle + column layout

        self.mat_table.toggled.connect(self._sync)
        self.mat_table.toggled.connect(self._emit)
        self.tab_loc.currentIndexChanged.connect(self._emit)
        self.label_coords.toggled.connect(self._sync)
        self.label_coords.toggled.connect(self._emit)
        self.coord_size.valueChanged.connect(self._emit)
        self._sync()

    def _sync(self, *_):
        self.tab_loc.setEnabled(self.mat_table.isChecked())
        self.coord_size.setEnabled(self.label_coords.isChecked())

    def _emit(self, *_):
        self.changed.emit()

    def options(self):
        return {
            "mat_table": self.mat_table.isChecked(),
            "tab_loc": self.tab_loc.currentText(),
            "label_coordinates": self.label_coords.isChecked(),
            "coord_label_size": self.coord_size.value(),
            "legend_ncol": _legend_option(self),
            "legend_frame": _legend_frame_option(self),
            "show_title": _show_title_option(self),
            "show_legend": _show_legend_option(self),
        }


class MeshDisplayPanel(QWidget):
    """Display options for a mesh plot (no boundary conditions)."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)
        _add_title_toggle(self, form)
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
        _add_legend_controls(self, form)

    def _emit(self, *_):
        self.changed.emit()

    def options(self):
        return {
            "show_nodes": self.show_nodes.isChecked(),
            "label_elements": self.label_elements.isChecked(),
            "label_nodes": self.label_nodes.isChecked(),
            "pad_frac": self.pad_frac.value(),
            "legend_ncol": _legend_option(self),
            "legend_frame": _legend_frame_option(self),
            "show_title": _show_title_option(self),
            "show_legend": _show_legend_option(self),
        }


class FeDataDisplayPanel(QWidget):
    """Display options for a seep/FEM data plot (mesh + boundary conditions).
    ``include_bc_symbol`` adds the BC symbol-size control (FEM data only — the
    seepage data plot has no such parameter)."""

    changed = Signal()

    def __init__(self, include_bc_symbol=False, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)
        _add_title_toggle(self, form)
        self.show_bc = QCheckBox("Boundary conditions")
        self.show_bc.setChecked(True)
        self.show_nodes = QCheckBox("Show nodes")
        self.label_elements = QCheckBox("Element numbers")
        self.label_nodes = QCheckBox("Node numbers")
        for c in (self.show_bc, self.show_nodes, self.label_elements, self.label_nodes):
            form.addRow("", c)
            c.toggled.connect(self._emit)

        self.alpha = _dspin(0.0, 1.0, 0.6, 0.05)   # match the inputs view's zone opacity
        form.addRow("Fill opacity", self.alpha)
        self.alpha.valueChanged.connect(self._emit)

        self.bc_symbol_size = None
        if include_bc_symbol:
            self.bc_symbol_size = _dspin(0.0, 1.0, 0.03, 0.01, decimals=3)
            form.addRow("BC symbol size", self.bc_symbol_size)
            self.bc_symbol_size.valueChanged.connect(self._emit)

        _add_legend_controls(self, form)

    def _emit(self, *_):
        self.changed.emit()

    def options(self):
        d = {
            "show_bc": self.show_bc.isChecked(),
            "show_nodes": self.show_nodes.isChecked(),
            "label_elements": self.label_elements.isChecked(),
            "label_nodes": self.label_nodes.isChecked(),
            "alpha": self.alpha.value(),
            "legend_ncol": _legend_option(self),
            "legend_frame": _legend_frame_option(self),
            "show_title": _show_title_option(self),
            "show_legend": _show_legend_option(self),
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
        _add_title_toggle(self, form)

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

        self.levels = _ispin(2, 100, 20, step=5)
        self.alpha = _dspin(0.0, 1.0, 0.85, 0.05)   # contour-fill wash opacity
        self.vector_scale = _dspin(0.001, 10.0, 0.05, 0.01, decimals=3)
        self.pad_frac = _dspin(0.0, 1.0, 0.05, 0.01)
        self.cmap = _make_cmap_combo("Spectral_r")
        self.cmap.setToolTip("Color ramp for the filled contours.")
        self.cbar_size = _dspin(0.1, 1.0, 1.0, 0.05)
        self.cbar_size.setToolTip("Length of the color-ramp legend (colorbar) as "
                                  "a fraction of the plot height.")

        self.flowlines = QCheckBox("Flow lines")
        self.flowlines.setChecked(True)
        self.vectors = QCheckBox("Velocity vectors")
        self.fill = QCheckBox("Filled contours")
        self.phreatic = QCheckBox("Phreatic surface")
        self.phreatic.setChecked(True)

        form.addRow("Variable", self.variable)
        form.addRow("", self.flowlines)
        form.addRow("", self.vectors)
        form.addRow("", self.fill)
        form.addRow("", self.phreatic)
        form.addRow("Base material", self.base_mat)
        form.addRow("Contour levels", self.levels)
        form.addRow("Color ramp", self.cmap)
        form.addRow("Colorbar size", self.cbar_size)
        form.addRow("Fill opacity", self.alpha)
        form.addRow("Vector scale", self.vector_scale)
        form.addRow("Padding", self.pad_frac)
        _add_legend_controls(self, form)

        self.variable.currentIndexChanged.connect(self._on_variable)
        self.base_mat.currentIndexChanged.connect(self._emit)
        self.cmap.currentIndexChanged.connect(self._emit)
        for s in (self.levels, self.alpha, self.vector_scale, self.pad_frac,
                  self.cbar_size):
            s.valueChanged.connect(self._emit)
        for c in (self.flowlines, self.vectors, self.fill, self.phreatic):
            c.toggled.connect(self._emit)
        self.vectors.toggled.connect(self._sync_enabled)
        self.fill.toggled.connect(self._sync_enabled)
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
        # The ramp and its colorbar only apply to the filled contours.
        self.cmap.setEnabled(self.fill.isChecked())
        self.cbar_size.setEnabled(self.fill.isChecked())

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
            "cmap": self.cmap.currentData(),
            "cbar_shrink": self.cbar_size.value(),
            "legend_ncol": _legend_option(self),
            "legend_frame": _legend_frame_option(self),
            "show_title": _show_title_option(self),
            "show_legend": _show_legend_option(self),
        }


class FemResultsDisplayPanel(QWidget):
    """Display options for an FEM results plot: plot type, deformation scale, and
    the mesh/reinforcement/vector toggles. Controls dim to the active plot type —
    color ramp + colorbar for shear strain, deform scale for deformation, and the
    boundary/node/vector controls for displacement vectors. A single "Element
    edges" toggle drives the mesh overlay in every plot type (the undeformed mesh
    on the deformation plot)."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        form = QFormLayout(self)
        _add_title_toggle(self, form)

        self.plot_type = QComboBox()
        for key, label in FEM_PLOT_TYPES:
            self.plot_type.addItem(label, key)

        self.deform_percent = _ispin(0, 100, 15, suffix=" %")
        self.deform_percent.setToolTip("Deformed-shape exaggeration as a percent "
                                       "of slope height.")

        self.cmap = _make_cmap_combo("coolwarm")
        self.cmap.setToolTip("Color ramp for the shear-strain contours.")
        # (No "Colorbar size" control: the colorbar now tracks the plot height
        # automatically via make_axes_locatable, so there's nothing to tune.)

        # Universal control, but with a per-plot-type default (and per-type memory
        # of the user's choice): off for the contour/vector plots (keep them
        # clean), on for the deformation plot (the undeformed reference).
        self._edges_state = {"shear_strain": False, "deformation": True,
                               "displace_vector": False}
        self._current_pt = self.plot_type.currentData()
        self.element_edges = QCheckBox("Element edges")
        self.element_edges.setChecked(self._edges_state.get(self._current_pt, False))
        self.element_edges.setToolTip(
            "Overlay the mesh element edges (light gray). On the deformation plot "
            "this is the original, undeformed mesh.")
        self.show_reinforcement = QCheckBox("Reinforcement")
        self.show_reinforcement.setChecked(True)
        self.label_elements = QCheckBox("Element numbers")

        # Displacement-vector-only controls.
        self.plot_boundary = QCheckBox("Boundary edges")
        self.plot_boundary.setChecked(True)
        self.plot_boundary.setToolTip("Outline the domain boundary (black).")
        self.plot_nodes = QCheckBox("Node dots")
        self.scale_vectors = QCheckBox("Auto-scale vectors")
        self.scale_vectors.setChecked(True)
        self.scale_vectors.setToolTip(
            "On: size arrows for visibility (relative magnitudes kept).\n"
            "Off: draw arrows at true displacement magnitude (may be tiny).")
        self.displacement_tolerance = _dspin(0.0, 1.0, 0.5, 0.05)
        self.displacement_tolerance.setToolTip(
            "Hide vectors below this fraction of the max displacement.")

        form.addRow("Plot type", self.plot_type)
        form.addRow("Color ramp", self.cmap)
        form.addRow("Deform", self.deform_percent)
        form.addRow("", self.element_edges)
        form.addRow("", self.show_reinforcement)
        form.addRow("", self.label_elements)
        form.addRow("", self.plot_boundary)
        form.addRow("", self.plot_nodes)
        form.addRow("", self.scale_vectors)
        form.addRow("Vector cutoff", self.displacement_tolerance)
        # No legend controls: the single-panel FEM result plots draw no legend.

        self.plot_type.currentIndexChanged.connect(self._on_plot_type)
        self.cmap.currentIndexChanged.connect(self._emit)
        self.deform_percent.valueChanged.connect(self._emit)
        self.displacement_tolerance.valueChanged.connect(self._emit)
        for c in (self.element_edges, self.show_reinforcement, self.label_elements,
                  self.plot_boundary, self.plot_nodes, self.scale_vectors):
            c.toggled.connect(self._emit)
        self._sync_enabled()

    @property
    def _vector_widgets(self):
        return (self.plot_boundary, self.plot_nodes,
                self.scale_vectors, self.displacement_tolerance)

    def _on_plot_type(self, *_):
        # "Element edges" carries a per-type default and remembers the user's
        # toggle per plot type: save the state we're leaving, load the one we're
        # entering (falling back to that type's default the first time it's seen).
        new_pt = self.plot_type.currentData()
        if new_pt != self._current_pt:
            self._edges_state[self._current_pt] = self.element_edges.isChecked()
            self.element_edges.blockSignals(True)
            self.element_edges.setChecked(self._edges_state.get(new_pt, False))
            self.element_edges.blockSignals(False)
            self._current_pt = new_pt
        self._sync_enabled()
        self.changed.emit()

    def _emit(self, *_):
        self.changed.emit()

    def _sync_enabled(self):
        pt = self.plot_type.currentData()
        # Vector-only controls.
        for w in self._vector_widgets:
            w.setEnabled(pt == "displace_vector")
        # Color ramp: shear-strain only. Deform scale: deformation only.
        self.cmap.setEnabled(pt == "shear_strain")
        self.deform_percent.setEnabled(pt == "deformation")

    def options(self):
        # One universal "Element edges" toggle drives the mesh overlay in every
        # plot type: show_mesh for the contour/deformation plots, plot_elements
        # for the vector plot.
        edges = self.element_edges.isChecked()
        return {
            "plot_type": self.plot_type.currentData(),
            "cmap": self.cmap.currentData(),
            "deform_percent": self.deform_percent.value(),
            "show_mesh": edges,
            "plot_elements": edges,
            "show_reinforcement": self.show_reinforcement.isChecked(),
            "label_elements": self.label_elements.isChecked(),
            "plot_boundary": self.plot_boundary.isChecked(),
            "plot_nodes": self.plot_nodes.isChecked(),
            "scale_vectors": self.scale_vectors.isChecked(),
            "displacement_tolerance": self.displacement_tolerance.value(),
            "show_title": _show_title_option(self),
            "show_legend": _show_legend_option(self),
        }
