"""Run-options dialogs (Phase 3).

``RunLemDialog`` collects the options for an LEM solve: method, number of slices,
analysis type (single surface / auto-search), surface type (circular /
non-circular), rapid drawdown, and a diagnostic-output toggle. Reliability is the
remaining analysis type (next increment).
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout,
    QGroupBox, QHeaderView, QLabel, QSpinBox, QTableWidget, QTableWidgetItem,
    QVBoxLayout,
)

LEM_METHODS = [
    ("oms", "Ordinary Method of Slices (OMS)"),
    ("bishop", "Bishop's Simplified"),
    ("janbu", "Janbu (Corrected)"),
    ("corps", "Corps of Engineers"),
    ("lowe", "Lowe & Karafiath"),
    ("spencer", "Spencer"),
    ("mprice", "Morgenstern-Price"),
]

ANALYSIS_TYPES = [("single_surface", "Single surface"), ("auto_search", "Auto search"),
                  ("reliability", "Reliability")]
SURFACE_TYPES = [("circular", "Circular"), ("noncircular", "Non-circular")]

MESH_ELEMENT_TYPES = [
    ("tri3", "Linear triangles (tri3)"),
    ("tri6", "Quadratic triangles (tri6)"),
    ("quad4", "Linear quads (quad4)"),
    ("quad8", "Quadratic quads (quad8)"),
    ("quad9", "Quadratic quads (quad9)"),
]


FEM_ANALYSIS_TYPES = [("single", "Single (fixed F)"), ("ssrm", "SSRM (find FS)")]
FEM_FAILURE_CRITERIA = [
    ("non_convergence", "Non-convergence"),
    ("displacement_limit", "Displacement limit"),
    ("displacement_increase", "Displacement increase"),
]


class RunFemDialog(QDialog):
    """Solve parameters for an FEM run (single trial or SSRM). Display options
    (plot type, deformation scale) live on the FEM Results view."""

    def __init__(self, parent=None, defaults=None):
        super().__init__(parent)
        self.setWindowTitle("Run FEM")
        defaults = defaults or {}

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.analysis = QComboBox()
        for key, label in FEM_ANALYSIS_TYPES:
            self.analysis.addItem(label, key)
        aidx = self.analysis.findData(defaults.get("analysis", "ssrm"))
        if aidx >= 0:
            self.analysis.setCurrentIndex(aidx)
        form.addRow("Analysis", self.analysis)

        self.F = QDoubleSpinBox()
        self.F.setRange(0.1, 10.0)
        self.F.setSingleStep(0.05)
        self.F.setValue(float(defaults.get("F", 1.0)))
        form.addRow("F (single)", self.F)

        self.F_min = QDoubleSpinBox()
        self.F_min.setRange(0.1, 10.0)
        self.F_min.setSingleStep(0.05)
        self.F_min.setValue(float(defaults.get("F_min", 1.0)))
        form.addRow("F min (SSRM)", self.F_min)

        self.F_max = QDoubleSpinBox()
        self.F_max.setRange(0.1, 20.0)
        self.F_max.setSingleStep(0.05)
        self.F_max.setValue(float(defaults.get("F_max", 2.0)))
        form.addRow("F max (SSRM)", self.F_max)

        self.tolerance = QDoubleSpinBox()
        self.tolerance.setDecimals(4)
        self.tolerance.setRange(0.0001, 1.0)
        self.tolerance.setValue(float(defaults.get("tolerance", 0.01)))
        form.addRow("Tolerance (SSRM)", self.tolerance)

        self.failure_criterion = QComboBox()
        for key, label in FEM_FAILURE_CRITERIA:
            self.failure_criterion.addItem(label, key)
        cidx = self.failure_criterion.findData(defaults.get("failure_criterion",
                                                            "non_convergence"))
        if cidx >= 0:
            self.failure_criterion.setCurrentIndex(cidx)
        form.addRow("Failure criterion", self.failure_criterion)

        layout.addLayout(form)
        note = QLabel("Plot type and deformation scale are set on the FEM Results "
                      "view after solving.")
        note.setWordWrap(True)
        layout.addWidget(note)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

        self.analysis.currentIndexChanged.connect(self._sync_enabled)
        self._sync_enabled()

    def _sync_enabled(self):
        single = self.analysis.currentData() == "single"
        self.F.setEnabled(single)
        for w in (self.F_min, self.F_max, self.tolerance, self.failure_criterion):
            w.setEnabled(not single)

    def options(self):
        return {
            "analysis": self.analysis.currentData(),
            "F": self.F.value(),
            "F_min": self.F_min.value(),
            "F_max": self.F_max.value(),
            "tolerance": self.tolerance.value(),
            "failure_criterion": self.failure_criterion.currentData(),
        }


SEEP_VARIABLES = [
    ("head", "Total head"),
    ("u", "Pore pressure"),
    ("v_mag", "Velocity magnitude"),
    ("i_mag", "Gradient magnitude"),
]


class RunSeepDialog(QDialog):
    """Solve parameters for a seepage run. Display options (variable, contours,
    flow lines, base material, …) are not here — they live on the Seep Solution
    view and re-render the cached solution without re-solving."""

    def __init__(self, parent=None, defaults=None, has_bc2=False):
        super().__init__(parent)
        self.setWindowTitle("Run Seepage")
        defaults = defaults or {}

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.bc = QComboBox()
        self.bc.addItem("Set 1", 1)
        if has_bc2:
            self.bc.addItem("Set 2 (rapid drawdown)", 2)
            self.bc.addItem("Both sets (1 & 2)", "both")
        idx = self.bc.findData(defaults.get("bc", 1))
        if idx >= 0:
            self.bc.setCurrentIndex(idx)
        form.addRow("BC set", self.bc)

        self.tol = QDoubleSpinBox()
        self.tol.setDecimals(8)
        self.tol.setRange(1e-10, 1.0)
        self.tol.setValue(float(defaults.get("tol", 1e-4)))
        form.addRow("Convergence tol", self.tol)

        layout.addLayout(form)
        note = QLabel("Display options (plotted variable, contours, flow lines, "
                      "base material) are set on the solution view after solving.")
        note.setWordWrap(True)
        layout.addWidget(note)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    def options(self):
        return {"bc": self.bc.currentData(), "tol": self.tol.value()}


class BuildMeshDialog(QDialog):
    """Options for building a finite-element mesh from the geometry.

    Target element size is either entered directly or auto-sized as
    ``(x_max - x_min) / size_divisions`` over the ground surface (see the
    main_seep / main_fem drivers)."""

    def __init__(self, parent=None, defaults=None):
        super().__init__(parent)
        self.setWindowTitle("Build mesh")
        defaults = defaults or {}

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.element_type = QComboBox()
        for key, label in MESH_ELEMENT_TYPES:
            self.element_type.addItem(label, key)
        idx = self.element_type.findData(defaults.get("element_type", "tri6"))
        if idx >= 0:
            self.element_type.setCurrentIndex(idx)
        form.addRow("Element type", self.element_type)

        self.auto_size = QCheckBox("Auto-size from geometry")
        self.auto_size.setChecked(bool(defaults.get("auto_size", True)))
        form.addRow("", self.auto_size)

        self.size_divisions = QSpinBox()
        self.size_divisions.setRange(5, 1000)
        self.size_divisions.setValue(int(defaults.get("size_divisions", 100)))
        self.size_divisions.setToolTip("Number of elements across the slope width.")
        form.addRow("Size divisions", self.size_divisions)

        self.target_size = QDoubleSpinBox()
        self.target_size.setRange(0.001, 1e6)
        self.target_size.setDecimals(3)
        self.target_size.setValue(float(defaults.get("target_size", 1.0)))
        form.addRow("Target element size", self.target_size)

        layout.addLayout(form)
        note = QLabel("Auto-size sets the target element size to the slope width "
                      "divided by the number of divisions.")
        note.setWordWrap(True)
        layout.addWidget(note)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Build")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

        self.auto_size.toggled.connect(self._sync_enabled)
        self._sync_enabled()

    def _sync_enabled(self):
        auto = self.auto_size.isChecked()
        self.size_divisions.setEnabled(auto)
        self.target_size.setEnabled(not auto)

    def options(self):
        return {
            "element_type": self.element_type.currentData(),
            "auto_size": self.auto_size.isChecked(),
            "size_divisions": self.size_divisions.value(),
            "target_size": self.target_size.value(),
        }


class RunLemDialog(QDialog):
    """Options for an LEM solve (single surface or auto-search; circular or not)."""

    def __init__(self, parent=None, defaults=None, slope_data=None):
        super().__init__(parent)
        self.setWindowTitle("Run LEM")
        defaults = defaults or {}
        slope_data = slope_data or {}

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.method = self._combo(LEM_METHODS, defaults.get("method", "bishop"))
        form.addRow("Method", self.method)

        self.analysis = self._combo(ANALYSIS_TYPES, defaults.get("analysis", "auto_search"))
        form.addRow("Analysis", self.analysis)

        # Only offer a surface-type choice when the file has both; otherwise
        # hardwire to whichever surface is present (shown as a fixed label).
        has_circ = bool(slope_data.get("circular") and slope_data.get("circles"))
        has_ncirc = bool(slope_data.get("non_circ"))
        if has_circ != has_ncirc:                      # exactly one defined
            self._fixed_surface = "circular" if has_circ else "noncircular"
            self.surface = None
            form.addRow("Surface", QLabel(dict(SURFACE_TYPES)[self._fixed_surface]))
        else:                                          # both (a real choice) or neither
            self._fixed_surface = None
            self.surface = self._combo(SURFACE_TYPES, defaults.get("surface", "circular"))
            form.addRow("Surface", self.surface)

        self.num_slices = QSpinBox()
        self.num_slices.setRange(5, 500)
        self.num_slices.setValue(int(defaults.get("num_slices", 40)))
        form.addRow("Number of slices", self.num_slices)

        self.rapid = QCheckBox("Rapid drawdown")
        self.rapid.setChecked(bool(defaults.get("rapid", False)))
        form.addRow("", self.rapid)

        self.diagnostic = QCheckBox("Diagnostic output (verbose log)")
        self.diagnostic.setChecked(bool(defaults.get("diagnostic", False)))
        form.addRow("", self.diagnostic)

        layout.addLayout(form)

        # Search-convergence tolerances — only meaningful for auto-search and
        # reliability (both drive the search loops). ``tol`` is the circular-search
        # geometric refinement tolerance and has no noncircular counterpart, so it
        # is greyed out for non-circular surfaces.
        self.tol_group = QGroupBox("Search tolerances")
        tform = QFormLayout(self.tol_group)

        self.fs_tol = QDoubleSpinBox()
        self.fs_tol.setDecimals(6)
        self.fs_tol.setRange(1e-8, 1.0)
        self.fs_tol.setSingleStep(1e-4)
        self.fs_tol.setValue(float(defaults.get("fs_tol", 5e-4)))
        self.fs_tol.setToolTip("Factor-of-safety convergence tolerance for the search.")
        tform.addRow("FS tol", self.fs_tol)

        self.tol = QDoubleSpinBox()
        self.tol.setDecimals(6)
        self.tol.setRange(1e-8, 1.0)
        self.tol.setSingleStep(1e-3)
        self.tol.setValue(float(defaults.get("tol", 1e-2)))
        self.tol.setToolTip("Geometric refinement tolerance (circular search only).")
        tform.addRow("Geometric tol", self.tol)

        self.max_iter = QSpinBox()
        self.max_iter.setRange(1, 1000)
        self.max_iter.setValue(int(defaults.get("max_iter", 50)))
        self.max_iter.setToolTip("Maximum search refinement iterations.")
        tform.addRow("Max iterations", self.max_iter)

        layout.addWidget(self.tol_group)

        note = QLabel("Single surface analyzes the first circle / the non-circular "
                      "surface as entered. Auto-search refines from there to the "
                      "critical surface.")
        note.setWordWrap(True)
        layout.addWidget(note)

        self.analysis.currentIndexChanged.connect(self._sync_tols)
        if self.surface is not None:
            self.surface.currentIndexChanged.connect(self._sync_tols)
        self._sync_tols()

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    @staticmethod
    def _combo(items, selected):
        combo = QComboBox()
        for key, label in items:
            combo.addItem(label, key)
        idx = combo.findData(selected)
        if idx >= 0:
            combo.setCurrentIndex(idx)
        return combo

    def _surface_value(self):
        return self._fixed_surface or self.surface.currentData()

    def _sync_tols(self):
        is_search = self.analysis.currentData() in ("auto_search", "reliability")
        self.tol_group.setEnabled(is_search)
        # Geometric tol applies to circular search only.
        self.tol.setEnabled(is_search and self._surface_value() == "circular")

    def options(self):
        return {
            "method": self.method.currentData(),
            "analysis": self.analysis.currentData(),
            "surface": self._surface_value(),
            "num_slices": self.num_slices.value(),
            "rapid": self.rapid.isChecked(),
            "diagnostic": self.diagnostic.isChecked(),
            "fs_tol": self.fs_tol.value(),
            "tol": self.tol.value(),
            "max_iter": self.max_iter.value(),
        }


# Import targets shown in the wizard (label, cad.DXF_TARGETS key). "Ignore" first
# so it's the obvious skip; "Material zone"/"Profile line" use the material column.
_DXF_TARGET_CHOICES = [
    ("Ignore", "ignore"),
    ("Material zone", "material_zone"),
    ("Profile line", "profile"),
    ("Piezometric line", "piezo"),
    ("Distributed load", "dload"),
    ("Reinforcement", "reinforce"),
    ("Failure circles", "circles"),
]
_MATERIAL_TARGETS = ("material_zone", "profile")


def _dxf_layer_summary(geom):
    """One-line description of a layer's geometry for the wizard."""
    bits = []
    for key, noun in (("closed", "zone"), ("open", "polyline"),
                      ("lines", "line"), ("circles", "circle"), ("points", "point")):
        n = len(geom.get(key) or [])
        if n:
            bits.append(f"{n} {noun}{'s' if n != 1 else ''}")
    return ", ".join(bits) or "(empty)"


class DxfImportDialog(QDialog):
    """Feature-aware DXF import wizard: map each layer to an xslope input feature.

    For every DXF layer the user picks a target (material zone, profile line,
    piezo line, distributed load, reinforcement, failure circles, or ignore) and,
    for material-zone / profile layers, the material it belongs to (same name →
    merge). Defaults are seeded from xslope's own export layer names and the
    geometry kind, but the user can override anything — so a DXF drawn in external
    CAD with arbitrary layer names maps just as well. ``result()`` returns
    ``{layer: {'target': key, 'material': name|None}}`` for
    ``ProjectDocument.build_from_dxf_mapping``.
    """

    def __init__(self, layers, suggest, parent=None):
        # layers: dict {name: geom} from read_dxf_layers (first-appearance order).
        # suggest: callable(name, geom) -> default target key (cad.suggest_dxf_target).
        super().__init__(parent)
        self.setWindowTitle("Import DXF — map layers to features")
        self.resize(680, 460)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "Map each DXF layer to an input feature (or Ignore).\n"
            "• Material zone / Profile line use the Material column (same name → merge).\n"
            "• Loads, reinforcement and circles import geometry only — set magnitudes "
            "and strengths in the editors afterward."))

        self.table = QTableWidget(len(layers), 4, self)
        self.table.setHorizontalHeaderLabels(["Layer", "Contents", "Import as", "Material"])
        self.table.verticalHeader().setVisible(False)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)
        hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.Stretch)

        self._rows = []   # (layer_name, target_combo, material_edit)
        for row, (name, geom) in enumerate(layers.items()):
            lyr = QTableWidgetItem(name)
            lyr.setFlags(lyr.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(row, 0, lyr)
            cont = QTableWidgetItem(_dxf_layer_summary(geom))
            cont.setFlags(cont.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(row, 1, cont)

            combo = QComboBox()
            for label, key in _DXF_TARGET_CHOICES:
                combo.addItem(label, key)
            default = suggest(name, geom)
            di = combo.findData(default)
            combo.setCurrentIndex(di if di >= 0 else 0)
            self.table.setCellWidget(row, 2, combo)

            # Material defaults to the layer name (PROFILE_ prefix stripped).
            mat_default = name[8:] if name.upper().startswith("PROFILE_") else name
            edit = QTableWidgetItem(mat_default)
            self.table.setItem(row, 3, edit)

            combo.currentIndexChanged.connect(
                lambda _i, r=row: self._sync_material_enabled(r))
            self._rows.append((name, combo, edit))
            self._sync_material_enabled(row)

        layout.addWidget(self.table, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _sync_material_enabled(self, row):
        """Grey out the Material cell unless the target uses it (zone / profile)."""
        _, combo, edit = self._rows[row]
        uses_mat = combo.currentData() in _MATERIAL_TARGETS
        flags = edit.flags()
        edit.setFlags((flags | Qt.ItemIsEditable | Qt.ItemIsEnabled) if uses_mat
                      else (flags & ~Qt.ItemIsEditable & ~Qt.ItemIsEnabled))

    def result(self):
        """{layer: {'target': key, 'material': name|None}}. Material is the column
        text for zone/profile targets (falling back to the layer name), else None."""
        out = {}
        for name, combo, edit in self._rows:
            target = combo.currentData()
            mat = None
            if target in _MATERIAL_TARGETS:
                mat = (edit.text() or "").strip() or name
            out[name] = {"target": target, "material": mat}
        return out
