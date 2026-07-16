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


FEM_ANALYSIS_TYPES = [("single", "Single (fixed F)"), ("ssrm", "SSRM (find FS)"),
                      ("reliability", "Reliability")]
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
        # Fit the dialog to its content so showing/hiding the reliability note
        # resizes the window instead of squeezing the fields.
        from PySide6.QtWidgets import QLayout
        layout.setSizeConstraint(QLayout.SetFixedSize)
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

        self.reliability_tol = QDoubleSpinBox()
        self.reliability_tol.setDecimals(4)
        self.reliability_tol.setRange(0.0001, 0.1)
        self.reliability_tol.setSingleStep(0.0005)
        self.reliability_tol.setValue(float(defaults.get("reliability_tol", 0.001)))
        self.reliability_tol.setToolTip(
            "Bisection tolerance for the reliability SSRM solves — tighter than a "
            "single run (default 0.001). TSPM amplifies factor-of-safety "
            "imprecision, so a tight tolerance keeps the reliability index stable.")
        form.addRow("Tolerance (Reliability)", self.reliability_tol)

        self.failure_criterion = QComboBox()
        for key, label in FEM_FAILURE_CRITERIA:
            self.failure_criterion.addItem(label, key)
        cidx = self.failure_criterion.findData(defaults.get("failure_criterion",
                                                            "non_convergence"))
        if cidx >= 0:
            self.failure_criterion.setCurrentIndex(cidx)
        form.addRow("Failure criterion", self.failure_criterion)

        # Optional surficial-failure filter. Off by default. A cohesionless slope's
        # true global minimum is an infinitely shallow surface-parallel "skin" slide
        # (FS = tan phi / tan beta, depth-independent), which masks the deep-seated
        # factor of safety; enabling this ignores any mechanism shallower than the
        # given depth so the search reports the deep answer instead.
        self.min_slip_on = QCheckBox("Ignore surficial (skin) failures")
        self.min_slip_on.setChecked(bool(defaults.get("min_slip_depth")))
        self.min_slip_on.setToolTip(
            "Exclude near-surface, surface-parallel failures from the SSRM search.\n\n"
            "A cohesionless slope's true critical mechanism is an infinitely shallow "
            "skin slide with FS = tan phi / tan beta — correct, but not the deep-seated "
            "answer a design usually wants. Turn this on and set a minimum depth below "
            "the ground surface; sweep the depth until FS stops rising (the plateau) to "
            "read the deep-seated factor of safety.\n\n"
            "Depth is in model length units. SSRM only.")
        form.addRow("", self.min_slip_on)

        self.min_slip_depth = QDoubleSpinBox()
        self.min_slip_depth.setDecimals(2)
        self.min_slip_depth.setRange(0.0, 1.0e6)
        self.min_slip_depth.setSingleStep(1.0)
        self.min_slip_depth.setValue(float(defaults.get("min_slip_depth") or 0.0))
        self.min_slip_depth.setToolTip("Minimum failure depth below the ground surface, "
                                       "in model length units.")
        form.addRow("Min slip depth", self.min_slip_depth)

        layout.addLayout(form)
        self._rel_note = QLabel(
            "Reliability uses the bracket above only to find the most-likely-value "
            "factor of safety; the ±σ perturbations then auto-bracket around it. It "
            "uses the material standard deviations (sigma columns) and runs 1+2N "
            "SSRM solves at Tolerance (Reliability).")
        self._rel_note.setWordWrap(True)
        # Cap the note to the form's width so it wraps (grows taller) instead of
        # widening the whole dialog under the SetFixedSize constraint.
        self._rel_note.setMaximumWidth(form.sizeHint().width())
        layout.addWidget(self._rel_note)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

        self.analysis.currentIndexChanged.connect(self._sync_enabled)
        self.min_slip_on.toggled.connect(self._sync_enabled)
        self._sync_enabled()

    def _sync_enabled(self):
        a = self.analysis.currentData()
        single = a == "single"
        # Single: only F. SSRM: F_min/F_max + the single-run tolerance. Reliability:
        # the bracket (used ONLY to find F_MLV — the perturbations auto-centre on it)
        # plus its own tighter tolerance; the single-run Tolerance does not apply.
        self.F.setEnabled(single)
        self.F_min.setEnabled(not single)
        self.F_max.setEnabled(not single)
        self.failure_criterion.setEnabled(not single)
        self.tolerance.setEnabled(a == "ssrm")
        self.reliability_tol.setEnabled(a == "reliability")
        self._rel_note.setVisible(a == "reliability")   # note only applies to reliability
        # Surficial-failure filter applies to the SSRM criterion only.
        self.min_slip_on.setEnabled(a == "ssrm")
        self.min_slip_depth.setEnabled(a == "ssrm" and self.min_slip_on.isChecked())

    def options(self):
        return {
            "analysis": self.analysis.currentData(),
            "F": self.F.value(),
            "F_min": self.F_min.value(),
            "F_max": self.F_max.value(),
            "tolerance": self.tolerance.value(),
            "reliability_tol": self.reliability_tol.value(),
            "failure_criterion": self.failure_criterion.currentData(),
            "min_slip_depth": (self.min_slip_depth.value()
                               if self.min_slip_on.isChecked()
                               and self.min_slip_depth.value() > 0 else None),
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

        # Composite surfaces: let a circle be truncated at the bottom of the model
        # and run along it. Off by default — the floor of a profile-line model is
        # max_depth, a search bound rather than a real impenetrable boundary.
        self.composite = QCheckBox("Composite surfaces (truncate circles at the base)")
        self.composite.setChecked(bool(defaults.get("composite", False)))
        self.composite.setToolTip(
            "Allow a circle deeper than the bottom of the model to be truncated at it "
            "and run along the base between the two crossings.\n\n"
            "Turn this on when the base is a real impenetrable boundary — bedrock, or a "
            "weak seam resting on it — because the critical mechanism there follows the "
            "base and no ordinary circle can reach it.\n\n"
            "Leave it off when the bottom of the model is just how deep you chose to "
            "look. Circular surfaces only.")
        form.addRow("", self.composite)

        # Grid seeding: a coarse grid-and-tangent sweep seeds the circular search
        # instead of (only) the circles sheet, protecting against the local-minimum
        # trap of a single bad starting circle. Circular auto-search only.
        self.grid_seed = QCheckBox("Grid search (auto-seed the circular search)")
        self.grid_seed.setChecked(bool(defaults.get("grid_seed", False)))
        self.grid_seed.setToolTip(
            "Sweep a grid of circle centers against a range of tangent elevations, "
            "derived from the slope geometry, and refine from the best circle of every "
            "competing family — plus your entered circles, if any.\n\n"
            "This is a GLOBAL search: it reports the most critical surface anywhere in "
            "the model. Without it the search only refines the neighborhood of your "
            "starting circles, and a single seed in the wrong place can converge to a "
            "local minimum that reads 20% or more too high, with no warning.\n\n"
            "Leave it off to interrogate a specific mechanism with your own circles.")
        form.addRow("", self.grid_seed)

        self.diagnostic = QCheckBox("Diagnostic output (verbose log)")
        self.diagnostic.setChecked(bool(defaults.get("diagnostic", False)))
        form.addRow("", self.diagnostic)

        # Optional surficial-failure filter (auto-search only): reject any trial
        # surface whose maximum depth is shallower than the given depth. Off by
        # default. Mirrors the SSRM filter so a cohesionless slope's shallow skin
        # mechanisms don't win the search over the deep-seated surface.
        self.min_slip_on = QCheckBox("Ignore surficial (skin) failures")
        self.min_slip_on.setChecked(bool(defaults.get("min_slip_depth")))
        self.min_slip_on.setToolTip(
            "Reject trial surfaces shallower than a minimum depth during the search.\n\n"
            "On a cohesionless slope the critical surface is an infinitely shallow "
            "skin slide; without this the search chases it instead of the deep-seated "
            "mechanism a design wants. Set a minimum depth below the ground surface and "
            "sweep it until the factor of safety stops rising (the plateau).\n\n"
            "Depth is in model length units. Auto-search only.")
        form.addRow("", self.min_slip_on)

        self.min_slip_depth = QDoubleSpinBox()
        self.min_slip_depth.setDecimals(2)
        self.min_slip_depth.setRange(0.0, 1.0e6)
        self.min_slip_depth.setSingleStep(1.0)
        self.min_slip_depth.setValue(float(defaults.get("min_slip_depth") or 0.0))
        self.min_slip_depth.setToolTip("Minimum failure depth below the ground surface, "
                                       "in model length units.")
        form.addRow("Min slip depth", self.min_slip_depth)

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
        self.min_slip_on.toggled.connect(self._sync_tols)
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
        circular = self._surface_value() == "circular"
        # Geometric tol applies to circular search only.
        self.tol.setEnabled(is_search and circular)
        # Grid seeding drives the circular search; meaningless for a single surface
        # or a non-circular search.
        self.grid_seed.setEnabled(is_search and circular)
        if not (is_search and circular):
            self.grid_seed.setChecked(False)
        # A composite surface is a truncated circle; the idea has no meaning for a
        # non-circular surface, whose points are placed by the user or the search.
        self.composite.setEnabled(circular)
        if not circular:
            self.composite.setChecked(False)
        # The surficial-failure filter rejects too-shallow trials, so it only has
        # meaning for the auto-search (single-surface just evaluates the given
        # surface; the reliability path does not thread it through).
        auto_search = self.analysis.currentData() == "auto_search"
        self.min_slip_on.setEnabled(auto_search)
        self.min_slip_depth.setEnabled(auto_search and self.min_slip_on.isChecked())

    def options(self):
        return {
            "method": self.method.currentData(),
            "analysis": self.analysis.currentData(),
            "surface": self._surface_value(),
            "num_slices": self.num_slices.value(),
            "rapid": self.rapid.isChecked(),
            "composite": self.composite.isChecked(),
            "grid_seed": self.grid_seed.isChecked(),
            "diagnostic": self.diagnostic.isChecked(),
            "fs_tol": self.fs_tol.value(),
            "tol": self.tol.value(),
            "max_iter": self.max_iter.value(),
            "min_slip_depth": (self.min_slip_depth.value()
                               if self.min_slip_on.isChecked()
                               and self.min_slip_depth.value() > 0 else None),
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


class GszImportDialog(QDialog):
    """GeoStudio (.gsz) import: choose which analysis to import.

    A .gsz needs no layer-mapping wizard — unlike a DXF, it already knows what its
    geometry means. The one thing it cannot decide for us is *which* analysis: a
    GeoStudio file routinely holds several over the same geometry, and they differ
    in more than the slip surface (materials are assigned per analysis, so the same
    slope can be a different soil in each). So this dialog just lists them.

    ``result()`` returns the chosen analysis ID.
    """

    def __init__(self, analyses, parent=None):
        # analyses: [{'id','name','kind','method'}, ...] from geostudio.list_analyses
        super().__init__(parent)
        self.setWindowTitle("Import GeoStudio — choose an analysis")
        self.resize(620, 360)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "This GeoStudio file contains the analyses below. Choose one to import.\n"
            "• Material zones, strengths and water conditions import automatically.\n"
            "• Analyses can use different materials over the same geometry, so the "
            "choice matters.\n"
            "• Reinforcement, loads and SLOPE/W's search are not imported — you will "
            "get a list of what was left out."))

        self.table = QTableWidget(len(analyses), 3, self)
        self.table.setHorizontalHeaderLabels(["Analysis", "Type", "Method"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)
        hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeToContents)

        self._ids = []
        for row, a in enumerate(analyses):
            self._ids.append(a["id"])
            self.table.setItem(row, 0, QTableWidgetItem(a.get("name") or ""))
            self.table.setItem(row, 1, QTableWidgetItem(a.get("kind") or ""))
            self.table.setItem(row, 2, QTableWidgetItem(a.get("method") or ""))
        self.table.selectRow(0)

        layout.addWidget(self.table, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        # Double-clicking a row is the natural "pick this one".
        self.table.doubleClicked.connect(self.accept)

    def result(self):
        """The chosen analysis ID."""
        row = self.table.currentRow()
        return self._ids[row if row >= 0 else 0]
