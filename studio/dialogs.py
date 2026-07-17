"""Run-options dialogs (Phase 3).

``RunLemDialog`` collects the options for an LEM solve: method, number of slices,
analysis type (single surface / auto-search), surface type (circular /
non-circular), rapid drawdown, and a diagnostic-output toggle. Reliability is the
remaining analysis type (next increment).
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
    QDoubleSpinBox, QFormLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel,
    QPushButton, QSpinBox, QStackedWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
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


class SensitivityDialog(QDialog):
    """Options for a sensitivity / design study (one dialog, two modes).

    Sensitivity mode sweeps several parameters, each +/- a percent about its
    current value (per-row overridable, with a one-click sigma-range preset that
    reuses the reliability sigma_* columns) — the result is a tornado diagram.
    Design mode sweeps ONE parameter across explicit From/To bounds in N steps and
    finds where FS meets a target. The sweepable set is any numeric material
    property plus the global k_seismic; geometry design stays in main_design.py.

    A thin caller of ``xslope.sensitivity``: ``options()`` returns plain
    dicts/lists the runner forwards straight to ``design()`` / ``sensitivity()``.
    """

    _GLOBAL_LABEL = "k_seismic (global)"
    # per-engine-mode output quantity: (short name, long label, design-target label)
    _OUTPUT = {
        "lem": ("FS", "Factor of Safety", "Target FS"),
        "fem": ("FS", "Factor of Safety", "Target FS"),
        "seep": ("q", "total discharge q", "Target q"),
    }

    def __init__(self, parent=None, defaults=None, slope_data=None, app_mode="lem"):
        super().__init__(parent)
        self.app_mode = app_mode if app_mode in self._OUTPUT else "lem"
        out_short, out_long, target_label = self._OUTPUT[self.app_mode]
        self._out_short = out_short
        titles = {"lem": "Sensitivity / Design study (LEM)",
                  "fem": "Sensitivity / Design study (FEM · SSRM)",
                  "seep": "Sensitivity / Design study (Seepage)"}
        self.setWindowTitle(titles[self.app_mode])
        defaults = defaults or {}
        slope_data = slope_data or {}

        from xslope.sensitivity import list_params
        self._params = list_params(slope_data, mode=self.app_mode)
        self._by_ref = {e["ref"]: e for e in self._params}
        # Group parameters for the picker's first combo: materials by name, plus
        # the pseudo-groups globals ('__global__') and seep boundary heads
        # ('__seep_bc__'). First-appearance order preserved.
        self._groups = []                 # (display, key)
        self._group_entries = {}          # key -> [entries]
        for e in self._params:
            if e["kind"] == "global":
                key, disp = "__global__", self._GLOBAL_LABEL
            elif e["kind"] == "seep_bc":
                key, disp = "__seep_bc__", "Boundary heads"
            else:
                key, disp = e["name"], e["name"]
            if key not in self._group_entries:
                self._group_entries[key] = []
                self._groups.append((disp, key))
            self._group_entries[key].append(e)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.mode = self._combo([("sensitivity", "Sensitivity (tornado)"),
                                 (f"design", f"Design ({out_short} vs one parameter)")],
                                defaults.get("mode", "sensitivity"))
        form.addRow("Mode", self.mode)

        # Engine-specific solver row(s). LEM keeps method + slices; FEM swaps in the
        # SSRM knobs (each step is a full SSRM solve); Seep takes a BC set + tol.
        self.method = self._combo(LEM_METHODS, defaults.get("method", "bishop"))
        self.num_slices = QSpinBox()
        self.num_slices.setRange(5, 500)
        self.num_slices.setValue(int(defaults.get("num_slices", 40)))
        if self.app_mode == "lem":
            form.addRow("Method", self.method)
            form.addRow("Number of slices", self.num_slices)
        elif self.app_mode == "fem":
            self._build_fem_solver_rows(form, defaults)
        elif self.app_mode == "seep":
            self._build_seep_solver_rows(form, defaults, slope_data)
        layout.addLayout(form)

        # --- parameter picker (shared by both modes) ------------------------
        picker = QGroupBox("Parameter")
        pform = QFormLayout(picker)
        self.material = QComboBox()
        for disp, key in self._groups:
            self.material.addItem(disp, key)
        pform.addRow("Material" if self.app_mode != "seep" else "Material / BC",
                     self.material)
        self.prop = QComboBox()
        pform.addRow("Property", self.prop)
        self.material.currentIndexChanged.connect(self._on_material_changed)
        self.prop.currentIndexChanged.connect(self._on_prop_changed)
        layout.addWidget(picker)

        # --- mode pages -----------------------------------------------------
        self.stack = QStackedWidget()
        self.stack.addWidget(self._build_sens_page(defaults))
        self.stack.addWidget(self._build_design_page(defaults, target_label))
        layout.addWidget(self.stack)

        # Re-search toggle is an LEM concept (FEM finds its own mechanism, seepage
        # has no failure surface) — only shown in LEM mode.
        self.search = QCheckBox("Re-search the critical surface at each step")
        self.search.setChecked(bool(defaults.get("search", True)))
        self.search.setToolTip(
            "On (recommended): re-run the search at every swept value, because the "
            "critical surface moves as the parameter changes — a fixed surface "
            "silently understates the sensitivity.\n\n"
            "Off: re-solve the entered surface only (much faster, but the answer is "
            "only right for that prescribed surface).")
        self.search.setVisible(self.app_mode == "lem")
        layout.addWidget(self.search)

        self.note = QLabel()
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self._ok = bb.button(QDialogButtonBox.Ok)
        self._ok.setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

        self.mode.currentIndexChanged.connect(self._on_mode_changed)
        self._on_material_changed()                 # populate property combo
        # Restore a remembered sensitivity table (from the previous run this session).
        for spec in defaults.get("_remember_params", []):
            e = self._by_ref.get(spec.get("ref"))
            if e is not None:
                self._add_row(e, pct=spec.get("pct", self.pct.value()),
                              use_sigma=spec.get("use_sigma", False))
        self._on_mode_changed()
        self.resize(560, 600)

    # --- engine-specific solver rows ----------------------------------------
    def _build_fem_solver_rows(self, form, defaults):
        fo = defaults.get("fem_opts", {}) or {}
        self.f_min = self._dspin(0.1, 10.0, float(fo.get("F_min", 1.0)))
        self.f_min.setSingleStep(0.05)
        self.f_max = self._dspin(0.1, 20.0, float(fo.get("F_max", 2.0)))
        self.f_max.setSingleStep(0.05)
        self.ssrm_tol = QDoubleSpinBox()
        self.ssrm_tol.setDecimals(4)
        self.ssrm_tol.setRange(0.0001, 1.0)
        self.ssrm_tol.setValue(float(fo.get("tolerance", 0.01)))
        self.failure_criterion = self._combo(FEM_FAILURE_CRITERIA,
                                              fo.get("failure_criterion",
                                                     "non_convergence"))
        form.addRow("F min (SSRM)", self.f_min)
        form.addRow("F max (SSRM)", self.f_max)
        form.addRow("Tolerance (SSRM)", self.ssrm_tol)
        form.addRow("Failure criterion", self.failure_criterion)
        note = QLabel("Each swept point is a full SSRM solve (output = FS) — minutes "
                      "per step. Runs in the background and is cancellable; keep the "
                      "step count small.")
        note.setWordWrap(True)
        form.addRow("", note)

    def _build_seep_solver_rows(self, form, defaults, slope_data):
        so = defaults.get("seep_opts", {}) or {}
        self.seep_bc = QComboBox()
        self.seep_bc.addItem("Set 1", 1)
        if (slope_data.get("seepage_bc2") or {}).get("specified_heads"):
            self.seep_bc.addItem("Set 2 (rapid drawdown)", 2)
        bidx = self.seep_bc.findData(so.get("bc", 1))
        if bidx >= 0:
            self.seep_bc.setCurrentIndex(bidx)
        self.seep_tol = QDoubleSpinBox()
        self.seep_tol.setDecimals(8)
        self.seep_tol.setRange(1e-10, 1.0)
        self.seep_tol.setValue(float(so.get("tol", 1e-4)))
        form.addRow("BC set", self.seep_bc)
        form.addRow("Convergence tol", self.seep_tol)
        note = QLabel("Output quantity = total discharge q through the section. Sweep "
                      "hydraulic conductivity / boundary heads; each step is one "
                      "seepage solve.")
        note.setWordWrap(True)
        form.addRow("", note)

    # --- page builders ------------------------------------------------------
    def _build_sens_page(self, defaults):
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(0, 0, 0, 0)

        controls = QHBoxLayout()
        controls.addWidget(QLabel("Default ±%"))
        self.pct = QDoubleSpinBox()
        self.pct.setRange(1.0, 99.0)
        self.pct.setDecimals(0)
        self.pct.setSuffix(" %")
        self.pct.setValue(float(defaults.get("default_pct", 20.0)))
        controls.addWidget(self.pct)
        controls.addSpacing(12)
        controls.addWidget(QLabel("Points"))
        self.n_points = QSpinBox()
        self.n_points.setRange(3, 31)
        self.n_points.setValue(int(defaults.get("n", 7)))
        self.n_points.setToolTip("Points per parameter's FS-vs-value curve "
                                 "(shown on click-through). The tornado uses the "
                                 "curve's two endpoints.")
        controls.addWidget(self.n_points)
        controls.addStretch(1)
        self.add_btn = QPushButton("Add parameter")
        self.add_btn.clicked.connect(self._on_add_clicked)
        controls.addWidget(self.add_btn)
        v.addLayout(controls)

        self.table = QTableWidget(0, 4, page)
        self.table.setHorizontalHeaderLabels(["Parameter", "±%", "σ range", ""])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionMode(QAbstractItemView.NoSelection)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)
        hh.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        v.addWidget(self.table, 1)
        self._rows = []            # [{ref, value, sigma, pct_spin, sigma_btn}]
        return page

    def _build_design_page(self, defaults, target_label="Target FS"):
        page = QWidget()
        form = QFormLayout(page)
        form.setContentsMargins(0, 0, 0, 0)
        self._design_echo = QLabel("—")
        f = self._design_echo.font()
        f.setBold(True)
        self._design_echo.setFont(f)
        form.addRow("Sweeping", self._design_echo)
        # Seepage sweeps hydraulic conductivity (k ~ 1e-5) and boundary heads, so
        # the bound/target spins need many decimals and a fine step or a tiny k
        # rounds to 0.000 in the display. LEM/FEM keep the compact 3-decimal spins.
        dec, step = (8, 1e-6) if self.app_mode == "seep" else (3, 1.0)
        self.d_from = self._dspin(-1e9, 1e9, float(defaults.get("low", 0.0)),
                                  decimals=dec, step=step)
        self.d_to = self._dspin(-1e9, 1e9, float(defaults.get("high", 1.0)),
                                decimals=dec, step=step)
        self.d_steps = QSpinBox()
        self.d_steps.setRange(2, 200)
        self.d_steps.setValue(int(defaults.get("steps", 11)))
        if self.app_mode == "seep":
            self.d_target = self._dspin(0.0, 1e9, float(defaults.get("target_fs", 1e-5)),
                                        decimals=8, step=1e-6)
        else:
            self.d_target = self._dspin(0.0, 100.0, float(defaults.get("target_fs", 1.5)))
        form.addRow("From", self.d_from)
        form.addRow("To", self.d_to)
        form.addRow("Steps", self.d_steps)
        form.addRow(target_label, self.d_target)
        self._design_seeded = "low" in defaults      # don't re-seed a remembered range
        return page

    # --- picker / mode logic ------------------------------------------------
    @staticmethod
    def _dspin(lo, hi, val, decimals=3, step=None):
        s = QDoubleSpinBox()
        s.setRange(lo, hi)
        s.setDecimals(decimals)                # decimals BEFORE value: a tiny value
        if step is not None:                   # must not round away on setValue
            s.setSingleStep(step)
        s.setValue(val)
        return s

    def _on_material_changed(self):
        self.prop.blockSignals(True)
        self.prop.clear()
        key = self.material.currentData()
        entries = self._group_entries.get(key, [])
        for e in entries:
            tag = "" if e["value"] is not None else "  (unset)"
            self.prop.addItem(f"{e['field']}{tag}", e["ref"])
        # a single-entry group (globals) has nothing to pick between
        self.prop.setEnabled(len(entries) > 1)
        self.prop.blockSignals(False)
        self._on_prop_changed()

    def _on_prop_changed(self):
        ref = self.prop.currentData()
        e = self._by_ref.get(ref)
        self._design_echo.setText(ref or "—")
        # Seed sensible From/To bounds around the current value the first time.
        if e is not None and not self._design_seeded:
            self._seed_design_bounds(e)

    def _seed_design_bounds(self, entry):
        val = entry.get("value")
        if val is None or val == 0:
            lo, hi = (0.0, 1.0) if entry["field"] != "k_seismic" else (0.0, 0.3)
        else:
            lo, hi = val * 0.5, val * 1.5
        self.d_from.setValue(lo)
        self.d_to.setValue(hi)

    def _on_mode_changed(self):
        design = self.mode.currentData() == "design"
        self.stack.setCurrentIndex(1 if design else 0)
        q = self._out_short
        if design:
            self.note.setText(
                f"Design: sweeps the one parameter above across [From, To] and "
                f"annotates where {q} meets the target (interpolated). If it never "
                f"crosses, the plot says which way to widen the range.")
        else:
            self.note.setText(
                f"Sensitivity: add parameters, each swept ±% about its value (or "
                f"click σ for a ±σ range). The result is a tornado; click a bar to "
                f"see that parameter's {q}-vs-value curve.")

    # --- sensitivity table --------------------------------------------------
    def _on_add_clicked(self):
        ref = self.prop.currentData()
        e = self._by_ref.get(ref)
        if e is None:
            return
        if any(r["ref"] == ref for r in self._rows):
            return                                   # already added
        self._add_row(e, pct=self.pct.value(), use_sigma=False)

    def _add_row(self, entry, pct, use_sigma):
        row = self.table.rowCount()
        self.table.insertRow(row)
        item = QTableWidgetItem(entry["ref"])
        item.setFlags(item.flags() & ~Qt.ItemIsEditable)
        self.table.setItem(row, 0, item)

        pct_spin = QDoubleSpinBox()
        pct_spin.setRange(1.0, 99.0)
        pct_spin.setDecimals(0)
        pct_spin.setSuffix(" %")
        pct_spin.setValue(float(pct))
        pct_spin.setMinimumWidth(72)         # keep the value visible past the arrows
        self.table.setCellWidget(row, 1, pct_spin)

        sigma = entry.get("sigma")
        sigma_btn = QPushButton("σ")
        sigma_btn.setCheckable(True)
        if sigma:
            sigma_btn.setToolTip(f"Use ±σ range: {entry['value'] - sigma:g} … "
                                 f"{entry['value'] + sigma:g}  (σ = {sigma:g})")
            sigma_btn.setChecked(bool(use_sigma))
        else:
            sigma_btn.setEnabled(False)
            sigma_btn.setToolTip("No standard deviation (sigma_*) for this property.")
        sigma_btn.toggled.connect(lambda on, s=pct_spin: s.setEnabled(not on))
        pct_spin.setEnabled(not sigma_btn.isChecked())
        self.table.setCellWidget(row, 2, sigma_btn)

        rm = QPushButton("✕")
        rm.setToolTip("Remove this parameter")
        rm.clicked.connect(lambda _=False, b=rm: self._remove_row(b))
        self.table.setCellWidget(row, 3, rm)

        self._rows.append({"ref": entry["ref"], "value": entry.get("value"),
                           "sigma": sigma, "pct_spin": pct_spin,
                           "sigma_btn": sigma_btn, "rm": rm})

    def _remove_row(self, btn):
        idx = next((i for i, r in enumerate(self._rows) if r["rm"] is btn), None)
        if idx is None:
            return
        self.table.removeRow(idx)
        self._rows.pop(idx)

    # --- result -------------------------------------------------------------
    def options(self):
        # 'engine_mode' selects the sweep engine (lem/fem/seep); 'mode' stays the
        # sensitivity-vs-design study type the runner already keys on.
        common = {
            "mode": self.mode.currentData(),
            "engine_mode": self.app_mode,
            "method": self.method.currentData(),
            "num_slices": self.num_slices.value(),
            "search": self.search.isChecked() if self.app_mode == "lem" else False,
        }
        if self.app_mode == "fem":
            common["fem_opts"] = {
                "F_min": self.f_min.value(), "F_max": self.f_max.value(),
                "tolerance": self.ssrm_tol.value(),
                "failure_criterion": self.failure_criterion.currentData(),
            }
        elif self.app_mode == "seep":
            common["seep_opts"] = {"bc": self.seep_bc.currentData(),
                                   "tol": self.seep_tol.value()}
        if common["mode"] == "design":
            common.update({
                "param": self.prop.currentData(),
                "low": self.d_from.value(),
                "high": self.d_to.value(),
                "steps": self.d_steps.value(),
                "target_fs": self.d_target.value(),
            })
        else:
            specs, remembered = [], []
            for r in self._rows:
                use_sigma = r["sigma_btn"].isChecked() and r["sigma"]
                if use_sigma:
                    specs.append({"ref": r["ref"], "low": r["value"] - r["sigma"],
                                  "high": r["value"] + r["sigma"]})
                else:
                    specs.append({"ref": r["ref"],
                                  "rel_range": r["pct_spin"].value() / 100.0})
                remembered.append({"ref": r["ref"], "pct": r["pct_spin"].value(),
                                   "use_sigma": bool(use_sigma)})
            common.update({"params": specs, "n": self.n_points.value(),
                           "default_pct": self.pct.value(),
                           # remembered table for the next session
                           "_remember_params": remembered})
        return common

    @staticmethod
    def _combo(items, selected):
        combo = QComboBox()
        for key, label in items:
            combo.addItem(label, key)
        idx = combo.findData(selected)
        if idx >= 0:
            combo.setCurrentIndex(idx)
        return combo


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


class Slide2ImportDialog(QDialog):
    """Rocscience Slide2 (.sli/.slim/.slmd) import: choose which scenario to
    import.

    Like a .gsz, a Slide2 model needs no layer-mapping wizard — it already knows
    what its geometry means. The one thing it cannot decide for us is *which*
    scenario: a ``.slmd`` routinely bundles several over the same geometry (a base
    case plus variants), and they can differ in more than the slip surface. So
    this dialog just lists them.

    ``result()`` returns the chosen scenario index.
    """

    def __init__(self, scenarios, parent=None):
        # scenarios: [{'index','name'}, ...] from slide2.list_scenarios
        super().__init__(parent)
        self.setWindowTitle("Import Slide2 — choose a scenario")
        self.resize(520, 360)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "This Slide2 file contains the scenarios below. Choose one to import.\n"
            "• Material zones, strengths and water conditions import automatically.\n"
            "• Scenarios can use different geometry and water conditions, so the "
            "choice matters.\n"
            "• Supports/anchors, loads and Slide2's search are not imported — you "
            "will get a list of what was left out."))

        self.table = QTableWidget(len(scenarios), 1, self)
        self.table.setHorizontalHeaderLabels(["Scenario"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)

        self._indices = []
        for row, s in enumerate(scenarios):
            self._indices.append(s["index"])
            self.table.setItem(row, 0, QTableWidgetItem(s.get("name") or ""))
        self.table.selectRow(0)

        layout.addWidget(self.table, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        # Double-clicking a row is the natural "pick this one".
        self.table.doubleClicked.connect(self.accept)

    def result(self):
        """The chosen scenario index."""
        row = self.table.currentRow()
        return self._indices[row if row >= 0 else 0]
