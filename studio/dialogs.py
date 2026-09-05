"""Run-options dialogs (Phase 3).

``RunLemDialog`` collects the options for an LEM solve: method, number of slices,
analysis type (single surface / auto-search), surface type (circular /
non-circular), rapid drawdown, and a diagnostic-output toggle. Probabilistic
reliability analysis has its own dialog (``ReliabilityDialog``), a sibling of the
Parametric study rather than an LEM/FEM run option.

Both Run dialogs also carry the transient-seepage controls, because a stability run
against a time-dependent seepage field depends on WHICH instant it reads and that
dependency should be stated where the run is started: ``SeepageTimeSelector`` names
the instant (and can store it in the model), and ``StageTimeFields`` — Run LEM only —
edits the two rapid-drawdown stage times at their point of use.
"""

from __future__ import annotations

import os

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QButtonGroup, QCheckBox, QComboBox, QDialog,
    QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QListWidget, QListWidgetItem,
    QMessageBox, QPushButton, QRadioButton, QSpinBox, QStackedWidget,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from .preflight_panel import (
    PreflightPanel, SEISMIC_NOTE_FEM, SEISMIC_NOTE_LEM, apply_capabilities,
    two_pane,
)
# The thin-zone toggle's default is declared beside the code that acts on it, so
# the checkbox and the mesh build can never disagree about what "off" means.
from .runners import REFINE_THIN_ZONES_DEFAULT

LEM_METHODS = [
    ("oms", "Ordinary Method of Slices (OMS)"),
    ("bishop", "Bishop's Simplified"),
    ("janbu", "Janbu (Corrected)"),
    ("corps", "Corps of Engineers"),
    ("lowe", "Lowe & Karafiath"),
    ("spencer", "Spencer"),
    ("mprice", "Morgenstern-Price"),
]

ANALYSIS_TYPES = [("single_surface", "Single surface"), ("auto_search", "Auto search")]
SURFACE_TYPES = [("circular", "Circular"), ("noncircular", "Non-circular")]

MESH_ELEMENT_TYPES = [
    ("tri3", "Linear triangles (tri3)"),
    ("tri6", "Quadratic triangles (tri6)"),
    ("quad4", "Linear quads (quad4)"),
    ("quad8", "Quadratic quads (quad8)"),
    ("quad9", "Quadratic quads (quad9)"),
]

# The two quadrilateral meshing styles the builder offers (mesh.build_mesh_from_
# polygons' ``quad_style``), in the order they are shown. Labels are ordinary FEA
# vocabulary; the tooltips are the whole of what a user needs to choose between
# them, including — for 'structured' — that unsweepable zones stay free-meshed.
QUAD_STYLES = [
    ("free", "Free (recommended)",
     "Robust on any section; honors the target element size; may leave a few "
     "triangles where quads don't fit."),
    ("structured", "Structured where possible",
     "Fills grid-like regions with rows and columns of quads; other regions fall "
     "back to free meshing."),
]
QUAD_STYLE_TRI_TIP = ("Quadrilateral styles apply to quad meshes only — choose a "
                      "quad element type to set one.")

# "Refine thin zones" — the whole of what a user needs to decide it. It says what
# the option guarantees (element rows across a thin zone), why it is on by default
# (the failure it prevents is silent), and that it costs nothing on a section with
# no thin zone. The mechanism differs by element family (a size field on triangles,
# a derived local size on quads) and is deliberately NOT in the tooltip: it is the
# mesher's business, not a choice the user makes.
REFINE_THIN_ZONES_TIP = (
    "Give every thin material zone about four element rows across its width. A "
    "band too thin to resolve does not fail — it returns a factor of safety that "
    "is too high — so this is on by default. A section with no thin zone meshes "
    "exactly as it would with this off.")

# "1D element size" — the model's main!D20 cell, edited where the mesh is built.
ELEMENT_SIZE_1D_TIP = (
    "Target element size along reinforcement and pile lines. Blank = the mesh "
    "element size. Refines the beam/bar elements and the soil they share nodes "
    "with.")


FEM_ANALYSIS_TYPES = [("single", "Single (fixed F)"), ("ssrm", "SSRM (find FS)")]
# v21 main!D22. 'rollers' first — it is the template's shipped value, what every
# file that declares nothing means, and the historical hardwired behaviour.
FEM_SIDE_BC = [("rollers", "Rollers (vertical movement free)"),
               ("fixed", "Fixed (both components clamped)")]
FEM_FAILURE_CRITERIA = [
    ("non_convergence", "Non-convergence"),
    ("hybrid", "Hybrid (non-convergence + displacement)"),
    ("displacement_limit", "Displacement limit"),
    ("displacement_increase", "Displacement increase"),
]


def stored_surface_family(slope_data):
    """The surface family this model states, for a deck that defines both.

    The file stores and the dialog edits: ``main!D24`` carries the answer (blank on
    every model that never had to choose), the dialog opens on it, and a change is
    written back there — so a both-family model reopens on the family its last run
    used instead of silently reverting to the circles. The legacy ``circular`` flag
    is the fallback for a model held in memory that has not been through the cell
    (an import, a document built in Studio), and circular is the fallback for a model
    that says nothing at all, because that is what the slicer does with both present.
    """
    fam = slope_data.get("surface_family")
    if fam:
        return "noncircular" if str(fam).strip().lower().startswith("non") else "circular"
    return "circular" if slope_data.get("circular", True) else "noncircular"


def _fmt_time(v):
    """A time as the file would carry it: exact enough to round-trip, clean for
    round numbers (100.0 -> '100'). Blank for None."""
    return "" if v is None else f"{float(v):.10g}"


def _optfloat(text):
    """Blank -> None, a number -> float, anything else -> the string ``'bad'`` so a
    caller can tell a typo from a deliberate blank. Mirrors the editors' optfloat
    contract, where a zero is a value and a blank is an absence."""
    text = (text or "").strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return "bad"


class SeepageTimeSelector(QGroupBox):
    """Which instant of a transient seepage solution a stability run reads.

    A transient seepage march produces a sequence of saved frames. An LEM or FEM run
    with ``u = seep`` consumes exactly ONE of them (two, for a rapid drawdown), and
    without this control the choice is invisible: it falls to wherever the results
    view's play bar happens to be sitting, or to the last saved frame. Three ways to
    name it, ordered by expected use:

    * **a saved frame** — the instants that certainly exist, so the run starts
      immediately;
    * **the frame shown in the results viewer** — the honest, named version of the
      play-bar coupling, offered only while a transient results tab is open;
    * **another time** — supported, but it is not a computed state yet, so the run
      re-marches the transient solution with that instant added to the save schedule.
      Fields are never interpolated between frames. A re-march costs a full re-solve,
      which on a long march is minutes rather than seconds.

    The choice governs THIS run. Ticking *Save as the model's stability time* also
    writes it to the model's ``tseep`` ``stability_time``, which is what makes a
    headless or scripted re-run reproduce it.

    With no transient solution in hand the whole group is disabled and says why —
    the capability-gating every optional control follows.
    """

    def __init__(self, parent=None, transient=None, current_time=None,
                 slope_data=None):
        super().__init__("Seepage time", parent)
        slope_data = slope_data if slope_data is not None else {}
        self._times = [float(t) for t in (transient or {}).get("times", [])]
        self._current_time = (None if current_time is None else float(current_time))
        ts = slope_data.get("tseep") or {}
        self._file_time = ts.get("stability_time")
        self._duration = ts.get("duration")
        unit = slope_data.get("time_unit")
        self._unit = f" {unit}" if unit else ""

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.mode = QComboBox()
        self.mode.addItem("Saved frame", "saved")
        if self._current_time is not None:
            self.mode.addItem("Frame shown in the results viewer", "current")
        self.mode.addItem("Another time (reruns the analysis)", "other")
        self.mode.setToolTip(
            "A saved frame and the viewer's frame are instant — the pore pressures "
            "already exist. Another time re-runs the transient seepage analysis with "
            "that instant added to the save schedule, because a field between two frames "
            "is not a solution of anything and is never interpolated.")
        form.addRow("Time", self.mode)

        self.frame = QComboBox()
        for t in self._times:
            self.frame.addItem(f"t = {_fmt_time(t)}{self._unit}", t)
        self.frame.setToolTip("The instants this solution actually saved.")
        form.addRow("Saved frame", self.frame)

        self.other = QLineEdit()
        self.other.setToolTip("A time in the model's declared time unit. It must lie "
                              "within the run duration.")
        form.addRow("Other time", self.other)

        layout.addLayout(form)

        self.write_back = QCheckBox("Save as the model's stability time")
        self.write_back.setToolTip(
            "Write this instant to the tseep sheet's stability_time, so a scripted or "
            "headless re-run of this model reads the same frame. Left off, the choice "
            "applies to this run only and the file is unchanged.")
        layout.addWidget(self.write_back)

        self.note = QLabel()
        self.note.setWordWrap(True)
        self.note.setMaximumWidth(420)
        layout.addWidget(self.note)

        if not self._times:
            self.setEnabled(False)
            self.note.setEnabled(True)
            self.note.setText(
                "No transient seepage solution is loaded, so there is no frame to "
                "choose. Run or reload a transient seepage analysis first; a steady "
                "run uses its single pore-pressure field.")
            return

        self._select_default()
        self.mode.currentIndexChanged.connect(self._sync)
        self.frame.currentIndexChanged.connect(self._sync)
        self.other.textChanged.connect(self._sync)
        self._sync()

    # --- state ----------------------------------------------------------
    def _select_default(self):
        """Open on the model's own answer: the file's stability_time when it declares
        one, otherwise the last saved frame — which is what a blank stability_time
        means, and the note says so."""
        if self._file_time is not None:
            idx = self._nearest_index(self._file_time)
            if idx is not None:
                self.frame.setCurrentIndex(idx)
                self.other.setText(_fmt_time(self._times[idx]))
                return
            # Declared, but this solution never saved it — open on free entry, which
            # is the path that can actually produce it.
            self.mode.setCurrentIndex(max(0, self.mode.findData("other")))
            self.other.setText(_fmt_time(self._file_time))
            return
        self.frame.setCurrentIndex(len(self._times) - 1)
        self.other.setText(_fmt_time(self._times[-1]))

    def _nearest_index(self, t):
        """Index of the saved frame at ``t``, or None when no frame is at that time."""
        if not self._times:
            return None
        i = min(range(len(self._times)), key=lambda k: abs(self._times[k] - t))
        span = (self._times[-1] - self._times[0]) or 1.0
        return i if abs(self._times[i] - t) <= 1e-6 * abs(span) + 1e-9 else None

    def _sync(self):
        mode = self.mode.currentData()
        self.frame.setEnabled(mode == "saved")
        self.other.setEnabled(mode == "other")
        t = self.time()
        if mode == "other":
            raw = _optfloat(self.other.text())
            if raw is None or raw == "bad":
                self.note.setText("Enter a time in the model's time unit.")
                return
            if self._duration is not None and not (0 < raw <= float(self._duration)):
                self.note.setText(
                    f"t = {_fmt_time(raw)}{self._unit} is outside the run "
                    f"(0 to {_fmt_time(self._duration)}{self._unit}); the run "
                    f"cannot land on it.")
                return
            if self._nearest_index(raw) is not None:
                self.note.setText(f"t = {_fmt_time(raw)}{self._unit} is already a "
                                  f"saved frame — nothing to re-solve.")
                return
            self.note.setText(
                f"t = {_fmt_time(raw)}{self._unit} is not a saved frame. The run "
                f"re-runs the transient seepage analysis with this instant added to "
                f"the save schedule — a full re-solve, which on a long run takes "
                f"minutes. It reports progress and can be cancelled.")
            return
        if mode == "current":
            self.note.setText(f"The results viewer is showing t = "
                              f"{_fmt_time(self._current_time)}{self._unit}.")
            return
        if self._file_time is None:
            last = self._times[-1]
            extra = ("" if t is None or abs(t - last) > 1e-12 else
                     "  This model sets no stability time, so this is also what a "
                     "scripted run reads: the LAST saved frame, usually the drained "
                     "end state.")
            self.note.setText(f"Reading t = {_fmt_time(t)}{self._unit}.{extra}")
        else:
            self.note.setText(
                f"Reading t = {_fmt_time(t)}{self._unit}. The model's stability time "
                f"is {_fmt_time(self._file_time)}{self._unit}.")

    # --- results --------------------------------------------------------
    def time(self):
        """The chosen instant, or ``None`` when there is nothing to choose (no
        solution) or the free entry is blank/unreadable."""
        if not self._times:
            return None
        mode = self.mode.currentData()
        if mode == "current":
            return self._current_time
        if mode == "other":
            raw = _optfloat(self.other.text())
            return None if (raw is None or raw == "bad") else raw
        return self.frame.currentData()

    def problem(self):
        """A user-facing reason this selection cannot be run, or ``None``."""
        if not self._times or self.mode.currentData() != "other":
            return None
        raw = _optfloat(self.other.text())
        if raw is None:
            return "Enter a seepage time, or pick a saved frame."
        if raw == "bad":
            return f"{self.other.text()!r} is not a number."
        if self._duration is not None and not (0 < raw <= float(self._duration)):
            return (f"The seepage time must be greater than 0 and no later than the "
                    f"run duration ({_fmt_time(self._duration)}{self._unit}).")
        return None

    def frame_selection(self):
        """What the model checks are told about this run's pore pressures.

        ``None`` when no transient solution is loaded — there is then nothing to
        stage, and a model that needs a seepage field and has none is genuinely
        incomplete. Otherwise ``{"times": [t]}``: a frame WILL be staged into the
        model before the solver starts, and this is which one. The list is empty
        while the free entry is unreadable; the run is refused for that reason (see
        :meth:`problem`), which is the accurate sentence, rather than for a seepage
        solution that is in fact right here.
        """
        if not self._times:
            return None
        t = self.time()
        return {"times": [] if t is None else [float(t)]}

    def options(self):
        """``None`` with no transient solution; otherwise
        ``{time, remarch, write_back}``. ``remarch`` is True only when the chosen
        instant is not already a saved frame."""
        if not self._times:
            return None
        t = self.time()
        if t is None:
            return None
        return {"time": float(t),
                "remarch": self._nearest_index(t) is None,
                "write_back": self.write_back.isChecked()}


class StageTimeFields(QGroupBox):
    """Rapid-drawdown stage times, edited at the point of use.

    ``stage_1`` and ``stage_2`` say which two instants of a transient seepage march
    the two drawdown stages read. They are pure EXTRACTION parameters — the drawdown
    schedule itself, how far and how fast the pool falls, lives in the boundary
    conditions — which is why they are safe to set from a Run dialog where a genuine
    physics input would not be. The file keeps storing them on the ``tseep`` sheet,
    because a headless run has to read them and because the march has to land on
    them; this is a second, point-of-use view of the same two values, and editing
    either place changes the model.

    If the times name instants the loaded solution never saved, the run re-marches
    with them added to the save schedule.
    """

    def __init__(self, parent=None, slope_data=None):
        super().__init__("Rapid-drawdown stage times", parent)
        slope_data = slope_data if slope_data is not None else {}
        ts = slope_data.get("tseep") or {}
        self._has_tseep = bool(slope_data.get("tseep"))
        self._duration = ts.get("duration")
        unit = slope_data.get("time_unit")
        self._unit = f" {unit}" if unit else ""
        self._initial = (ts.get("stage_1"), ts.get("stage_2"))

        form = QFormLayout(self)
        self.stage_1 = QLineEdit(_fmt_time(self._initial[0]))
        self.stage_1.setToolTip(
            "The instant whose pore pressures seed drawdown stage 1 (before the pool "
            "falls). Set both stage times, or neither.")
        self.stage_2 = QLineEdit(_fmt_time(self._initial[1]))
        self.stage_2.setToolTip(
            "The instant whose pore pressures drive drawdown stage 2 (after the pool "
            "has fallen). Must be later than stage 1.")
        form.addRow("Stage 1 time" + self._unit, self.stage_1)
        form.addRow("Stage 2 time" + self._unit, self.stage_2)

        self.note = QLabel()
        self.note.setWordWrap(True)
        self.note.setMaximumWidth(420)
        form.addRow(self.note)

        if not self._has_tseep:
            self.setEnabled(False)
            self.note.setEnabled(True)
            self.note.setText(
                "This model has no transient seepage sheet, so a rapid drawdown reads "
                "its two stages from the classic seep.csv / seep2.csv pair instead of "
                "from stage times.")
            return
        self.note.setText(
            "Stored on the tseep sheet, so a scripted run reads the same two "
            "instants. A stage time the loaded solution never saved re-runs the "
            "transient seepage analysis with it added to the save schedule.")

    def values(self):
        """``(stage_1, stage_2)`` as floats or ``None``; ``'bad'`` for a typo."""
        return _optfloat(self.stage_1.text()), _optfloat(self.stage_2.text())

    def problem(self):
        """A user-facing reason these stage times cannot be run, or ``None``."""
        if not self._has_tseep:
            return None
        s1, s2 = self.values()
        for name, v in (("Stage 1", s1), ("Stage 2", s2)):
            if v == "bad":
                return f"{name} time is not a number."
        if s1 is None and s2 is None:
            return None                    # no staging; the classic two-file path
        if s1 is None or s2 is None:
            return "Set BOTH stage times, or neither."
        if s1 >= s2:
            return "Stage 1 must be earlier than stage 2."
        if self._duration is not None and s2 > float(self._duration):
            return (f"Stage 2 is later than the run duration "
                    f"({_fmt_time(self._duration)}{self._unit}).")
        return None

    def changed(self):
        """True when either field differs from what the model carries."""
        s1, s2 = self.values()
        return (s1, s2) != self._initial

    def options(self):
        """``None`` when there is nothing to store; otherwise
        ``{stage_1, stage_2, changed}``."""
        if not self._has_tseep:
            return None
        s1, s2 = self.values()
        if s1 == "bad" or s2 == "bad":
            return None
        return {"stage_1": s1, "stage_2": s2, "changed": self.changed()}


class SsrExcludeDialog(QDialog):
    """Pick material zones to EXCLUDE from strength reduction (RS2's per-material
    Apply_SSR flag / "SSR Exclusion Area"). Excluded zones keep their full c and
    tan(phi) at every trial F, forcing the mechanism up into the reduced zones.

    One checkbox per zone; CHECKED = included (reduced), the default. Unchecking a
    zone excludes it. Returns the list of excluded names."""

    def __init__(self, parent=None, material_names=None, excluded=None):
        super().__init__(parent)
        self.setWindowTitle("SSR exclusions")
        material_names = list(material_names or [])
        excluded = set(excluded or [])

        layout = QVBoxLayout(self)
        note = QLabel(
            "Uncheck a material zone to hold it at FULL strength during the SSRM "
            "(an SSR exclusion). Excluded zones are not divided by the trial factor, "
            "so the failure mechanism is forced up into the reduced zones — used, for "
            "example, to keep a stiff foundation from carrying the critical surface.")
        note.setWordWrap(True)
        note.setMaximumWidth(360)
        layout.addWidget(note)

        self._checks = []
        if not material_names:
            layout.addWidget(QLabel("No material zones found."))
        for name in material_names:
            cb = QCheckBox(name)
            cb.setChecked(name not in excluded)      # checked = included (reduced)
            layout.addWidget(cb)
            self._checks.append((name, cb))

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

    def excluded(self):
        """Names of zones left UNchecked = excluded from reduction."""
        return [name for name, cb in self._checks if not cb.isChecked()]


class RunFemDialog(QDialog):
    """Solve parameters for an FEM run (single trial or SSRM). Display options
    (plot type, deformation scale) live on the FEM Results view.

    Like the other Run dialogs it carries the model checks in a column beside the
    controls (:func:`studio.preflight_panel.two_pane`), so an input the finite
    element engine would refuse — a blank Poisson's ratio, a material with no
    tensile cap, water standing on bare ground — is stated before the solve rather
    than surfacing as a bracket that will not close. When the
    seismic coefficient is nonzero the panel also carries the FEM's sign
    convention, which is not the limit-equilibrium engine's.

    When the model carries a transient seepage solution the dialog also carries the
    seepage-time selector, so the instant the pore pressures are read from is named
    at the run rather than inherited from wherever the results viewer was left."""

    def __init__(self, parent=None, defaults=None, material_names=None,
                 slope_data=None, document=None, transient=None, current_time=None):
        super().__init__(parent)
        self.setWindowTitle("Run FEM")
        defaults = defaults or {}
        self._sd = slope_data if slope_data is not None else {}
        self._material_names = list(material_names or [])
        self._ssr_exclude = [n for n in (defaults.get("ssr_exclude") or [])
                             if n in self._material_names]

        # The run controls, built into the left column of the two-pane layout the
        # dialog is assembled into at the end (see two_pane): unparented here, so
        # the layout can be moved into that column.
        layout = QVBoxLayout()
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
        self.tolerance.setToolTip(
            "How narrow the F bracket must get before the search stops — the "
            "width of [F min, F max] after bisection, not a solver convergence "
            "tolerance.")
        form.addRow("Tolerance (SSRM)", self.tolerance)

        # Per-trial viscoplastic iteration budget. Trials just below the true FS
        # need far more iterations to reach equilibrium than trials far from it;
        # too small a budget makes near-failure trials read as failed and biases
        # the reported FS low. Previously only the API exposed this
        # (solve_fem/solve_ssrm max_iterations); a Studio run was pinned at the
        # engine default with no way to reach the converged plateau.
        self.max_iterations = QSpinBox()
        self.max_iterations.setRange(500, 100000)
        self.max_iterations.setSingleStep(500)
        self.max_iterations.setValue(int(defaults.get("max_iterations") or 12000))
        self.max_iterations.setToolTip(
            "Viscoplastic iteration budget for EACH trial F (and single-F runs). "
            "It is a budget, not a hard stop: a trial that reaches it with the "
            "out-of-balance force still falling is given another budget's worth, "
            "again and again, up to the iteration ceiling below. A residual that "
            "stops improving is reported but never ends a trial on its own; a "
            "trial whose movement is clearly running away is declared failed "
            "early and does not spend the rest of its budget. "
            "Raise it if the reported FS keeps climbing when you do; it has "
            "plateaued when raising it further changes nothing.")
        form.addRow("Max iterations per trial", self.max_iterations)

        # Hard stop on the automatic budget extension. A trial that reaches THIS
        # while still improving is inconclusive - it did not fail, and it did not
        # settle - so the bisection does not count it as a failure and says so.
        self.max_iterations_ceiling = QSpinBox()
        self.max_iterations_ceiling.setRange(500, 500000)
        self.max_iterations_ceiling.setSingleStep(1000)
        self.max_iterations_ceiling.setValue(
            int(defaults.get("max_iterations_ceiling") or 50000))
        self.max_iterations_ceiling.setToolTip(
            "Hard stop on the automatic budget extension above. A trial that "
            "reaches this ceiling with its out-of-balance force still falling is "
            "reported INCONCLUSIVE - neither settled nor failed - and is not "
            "counted as a failure. The factor of safety is still the final "
            "bracket's midpoint, as on any other run; what changes is that the "
            "bracket's upper edge is an undecided trial rather than a measured "
            "failure, and the Log says so.")
        form.addRow("Iteration ceiling", self.max_iterations_ceiling)

        # Side boundary condition (v21 main!D22). Applies to both a single trial and
        # the SSRM — it is part of how the model is restrained, not part of the
        # reduction — so it is never gated by the analysis type.
        self.side_bc = QComboBox()
        for key, label in FEM_SIDE_BC:
            self.side_bc.addItem(label, key)
        sidx = self.side_bc.findData(str(defaults.get("side_bc") or "rollers").lower())
        self.side_bc.setCurrentIndex(sidx if sidx >= 0 else 0)
        self.side_bc.setToolTip(
            "What holds the left and right edges of the model.\n\n"
            "Rollers (the default, and every file that declares nothing) fixes the "
            "horizontal component and leaves the vertical free, so the truncated "
            "ground can still settle under its own weight.\n\n"
            "Fixed clamps both components, which is what RS2 does on its side "
            "boundaries. It is a vendor-parity option, not a better model: fixing the "
            "sides adds shear restraint the real ground does not have, and stiffens a "
            "domain truncated close to the slope.")
        form.addRow("Side BC", self.side_bc)

        # K0 initial stress (v19). Off = the historical gravity turn-on, where the
        # initial lateral stress is whatever plane-strain elasticity gives
        # (nu/(1-nu)·sigma_v). Seeded from the file's main!D16 when it declares one.
        self.k0_on = QCheckBox("K0 initial stress")
        self.k0_on.setChecked(defaults.get("k0") is not None)
        self.k0_on.setToolTip(
            "Start from an at-rest in-situ stress state instead of switching gravity "
            "on in one step.\n\n"
            "Off, the initial lateral stress is set by the STIFFNESS, not the soil: "
            "sigma_h = nu/(1-nu)·sigma_v, about 0.43·sigma_v at nu = 0.3. Real "
            "compacted fill and overconsolidated clay sit at K0 = 1 and above, and "
            "under-confining a thin structural column (a reinforced-soil block) "
            "lowers its factor of safety.\n\n"
            "On, sigma_v comes from the soil overburden and sigma_h = K0·sigma_v "
            "in-plane and out-of-plane; the solver then iterates that state to "
            "equilibrium under the body forces.")
        form.addRow("", self.k0_on)

        self.k0 = QDoubleSpinBox()
        self.k0.setDecimals(3)
        self.k0.setRange(0.001, 10.0)
        self.k0.setSingleStep(0.05)
        self.k0.setValue(float(defaults.get("k0") or 1.0))
        self.k0.setToolTip("At-rest lateral earth pressure coefficient K0 "
                           "(sigma_h = K0·sigma_v).")
        form.addRow("K0", self.k0)

        # Tension SRF (v19). A no-op wherever no material declares a cutoff above
        # zero — see _tensile_cap_state, which is what dims it.
        self.tension_srf = QCheckBox("Reduce the tensile cap with F (Tension SRF)")
        self.tension_srf.setChecked(bool(defaults.get("tension_srf", True)))
        self._tension_srf_tip = (
            "Divide each material's tensile-strength cap by the trial F alongside c "
            "and tan(phi), so the factor of safety is the factor on the WHOLE "
            "strength envelope — shear and tensile.\n\n"
            "This is RS2's tensilestrength_SRF = 1 (the setting behind its published "
            "verification set) and what Plaxis does. Off holds each cap at its "
            "authored value through the bisection.")
        self.tension_srf.setToolTip(self._tension_srf_tip)
        form.addRow("", self.tension_srf)

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

        # SSR exclusions: material zones held at full strength during reduction.
        # The button opens a checkbox picker; the label summarises the current set.
        self.ssr_exclude_btn = QPushButton("SSR exclusions…")
        self.ssr_exclude_btn.setToolTip(
            "Hold selected material zones at full strength during the SSRM "
            "(RS2's SSR Exclusion Area). Default: every zone reduced.")
        self.ssr_exclude_btn.clicked.connect(self._edit_ssr_exclude)
        self.ssr_exclude_label = QLabel()
        _ssr_row = QHBoxLayout()
        _ssr_row.setContentsMargins(0, 0, 0, 0)
        _ssr_row.addWidget(self.ssr_exclude_btn)
        _ssr_row.addWidget(self.ssr_exclude_label, 1)
        form.addRow("SSR exclusions", _ssr_row)
        self._refresh_ssr_label()

        # At-failure mechanism capture (see solve_ssrm's capture_failure_state):
        # once the bracket resolves, re-solve ONCE just beyond critical F with the
        # displacement cap off so the unconverged field develops the failure
        # MECHANISM the deformation/vector figures render, instead of the diffuse
        # settlement of the last converged trial. On by default (the figure path);
        # turning it off skips the extra solve. FS/the bracket are unaffected either way.
        self.capture_failure_state = QCheckBox("Capture failure-state mechanism")
        self.capture_failure_state.setChecked(
            bool(defaults.get("capture_failure_state", True)))
        self.capture_failure_state.setToolTip(
            "After the bracket resolves, re-solve once just beyond critical F "
            "(cap and early-exit off) so the unconverged viscoplastic field "
            "develops the failure mechanism for the deformation/displacement-"
            "vector figures. Off skips this extra solve; FS and the bracket are "
            "unaffected either way.")
        form.addRow("", self.capture_failure_state)

        self.capture_margin = QDoubleSpinBox()
        self.capture_margin.setDecimals(2)
        self.capture_margin.setRange(0.02, 0.5)
        self.capture_margin.setSingleStep(0.01)
        self.capture_margin.setValue(float(defaults.get("capture_margin", 0.15)))
        self.capture_margin.setToolTip(
            "How far beyond the critical F the mechanism snapshot is solved "
            "(fraction of FS).")
        form.addRow("Capture margin", self.capture_margin)

        # Iteration budget for the capture solve. Off (unchecked) mirrors
        # solve_ssrm's own default: max(max_iterations, 3000).
        self.capture_iter_on = QCheckBox("Set capture iteration budget")
        self.capture_iter_on.setChecked(
            defaults.get("capture_max_iterations") is not None)
        self.capture_iter_on.setToolTip(
            "Override the iteration ceiling for the failure-state capture solve. "
            "Left off, the budget is automatic: max(max iterations, 3000).")
        form.addRow("", self.capture_iter_on)

        self.capture_max_iterations = QSpinBox()
        self.capture_max_iterations.setRange(100, 100000)
        self.capture_max_iterations.setSingleStep(100)
        self.capture_max_iterations.setValue(
            int(defaults.get("capture_max_iterations") or 3000))
        self.capture_max_iterations.setToolTip(
            "Iteration ceiling for the failure-state capture solve.")
        form.addRow("Capture max iterations", self.capture_max_iterations)

        layout.addLayout(form)

        # Which instant of a transient seepage march this run reads. Shown only when
        # the model has a transient solution OR a tseep sheet (dimmed, with the
        # reason, in the second case) — a steady model has one field and no choice.
        self.seep_time = None
        if transient or self._sd.get("tseep"):
            self.seep_time = SeepageTimeSelector(self, transient=transient,
                                                 current_time=current_time,
                                                 slope_data=self._sd)
            layout.addWidget(self.seep_time)

        self.preflight = PreflightPanel(
            analysis=lambda: ("ssrm" if self.analysis.currentData() == "ssrm"
                              else "fem"),
            slope_data=self._sd, document=document,
            selection_fn=self._selection,
            notes=(SEISMIC_NOTE_FEM,), parent=self)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self._ok = bb.button(QDialogButtonBox.Ok)
        self._ok.setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)

        # Controls left, checks right. ``fixed`` keeps the old behaviour of fitting
        # the window to its content when the min-slip and SSR rows appear or go.
        two_pane(self, layout, self.preflight, bb, fixed=True)

        self.analysis.currentIndexChanged.connect(self._sync_enabled)
        self.k0_on.toggled.connect(self._sync_enabled)
        self.min_slip_on.toggled.connect(self._sync_enabled)
        self.capture_failure_state.toggled.connect(self._sync_enabled)
        self.capture_iter_on.toggled.connect(self._sync_enabled)
        self.preflight.changed.connect(self._sync_run)
        # The chosen instant is part of what the checks are told, so it re-evaluates
        # them rather than only the button.
        if self.seep_time is not None:
            self.seep_time.mode.currentIndexChanged.connect(self._recheck)
            self.seep_time.frame.currentIndexChanged.connect(self._recheck)
            self.seep_time.other.textChanged.connect(self._recheck)
        self._sync_enabled()
        self._sync_run()

    def _selection(self):
        """What the model checks are asked about — the transient frame this run will
        stage, which is what makes ``u = seep`` satisfiable before it is staged. See
        :meth:`RunLemDialog._seep_frame`."""
        return {"seep_frame": (self.seep_time.frame_selection()
                               if self.seep_time is not None else None)}

    def _recheck(self):
        self.preflight.refresh()

    def _sync_run(self):
        """Run refuses on a preflight ERROR, and on a seepage time it cannot read."""
        reason = self.preflight.block_reason() if self.preflight.blocked else None
        if reason is None and self.seep_time is not None:
            reason = self.seep_time.problem()
        self._ok.setEnabled(reason is None)
        self._ok.setToolTip(reason or "")

    def _tensile_cap_state(self):
        """What this model gives the Tension SRF setting to act on: ``"reducible"``
        where some material declares a cutoff above zero, ``"zero"`` where every
        declared cutoff is 0, and ``"none"`` where no material declares one.

        The rule lives in :func:`xslope.preflight.has_reducible_tensile_cap`, so
        the control this dims and the ``main.tension_srf_unset`` note beside it in
        the checks column read the same model the same way. The three states are
        here rather than there because only the tooltip needs to tell the last two
        apart: a no-tension model (every cap 0) and an unbounded one (every cap
        blank) both leave the trial factor nothing to divide, but they are
        opposite modelling statements and the reason given should say which.
        """
        from xslope.preflight import (declares_tensile_cap,
                                      has_reducible_tensile_cap)
        mats = (self._sd or {}).get("materials") or []
        if any(has_reducible_tensile_cap(m) for m in mats):
            return "reducible"
        return "zero" if any(declares_tensile_cap(m) for m in mats) else "none"

    def _sync_enabled(self):
        a = self.analysis.currentData()
        single = a == "single"
        # Single: only F. SSRM: F_min/F_max + the single-run tolerance.
        self.F.setEnabled(single)
        self.F_min.setEnabled(not single)
        self.F_max.setEnabled(not single)
        self.failure_criterion.setEnabled(not single)
        self.tolerance.setEnabled(a == "ssrm")
        # K0 initialization applies to both a single trial and the SSRM.
        self.k0.setEnabled(self.k0_on.isChecked())
        # The tensile cap is only reduced during a strength reduction, and only
        # where some material declares a cutoff ABOVE ZERO: a blank t_cut is no
        # cap at all, and a cutoff of 0 divides to 0 at every trial factor. In
        # both cases the box changes no number, so it is dimmed rather than
        # offered, with the reason in its tooltip.
        cap = self._tensile_cap_state()
        self.tension_srf.setEnabled((not single) and cap == "reducible")
        self.tension_srf.setToolTip(self._tension_srf_tip + {
            "reducible": "",
            "zero": "\n\nDimmed: every tensile cutoff on this model is 0 — a soil "
                    "that carries no tension — and 0 divided by the trial factor "
                    "is still 0, so there is nothing here to reduce.",
            "none": "\n\nDimmed: no material on this model declares a t_cut, so "
                    "there is no cap to reduce.",
        }[cap])
        # Surficial-failure filter applies to the SSRM criterion only.
        self.min_slip_on.setEnabled(a == "ssrm")
        self.min_slip_depth.setEnabled(a == "ssrm" and self.min_slip_on.isChecked())
        # SSR exclusions apply to the SSRM criterion only.
        self.ssr_exclude_btn.setEnabled(a == "ssrm" and bool(self._material_names))
        # At-failure mechanism capture is an SSRM-only extra (a single trial has
        # no bracket to capture beyond).
        self.capture_failure_state.setEnabled(not single)
        cap_on = (not single) and self.capture_failure_state.isChecked()
        self.capture_margin.setEnabled(cap_on)
        self.capture_iter_on.setEnabled(cap_on)
        self.capture_max_iterations.setEnabled(cap_on and self.capture_iter_on.isChecked())
        # A single trial and an SSRM are different analysis types to the registry
        # (the strength-reduction rules apply only to the second), so the findings
        # follow the selector.
        self.preflight.refresh()

    def _refresh_ssr_label(self):
        if self._ssr_exclude:
            self.ssr_exclude_label.setText("excluded: " + ", ".join(self._ssr_exclude))
        else:
            self.ssr_exclude_label.setText("none (all zones reduced)")

    def _edit_ssr_exclude(self):
        dlg = SsrExcludeDialog(self, material_names=self._material_names,
                               excluded=self._ssr_exclude)
        if dlg.exec():
            self._ssr_exclude = dlg.excluded()
            self._refresh_ssr_label()

    def options(self):
        return {
            "analysis": self.analysis.currentData(),
            "F": self.F.value(),
            "F_min": self.F_min.value(),
            "F_max": self.F_max.value(),
            "tolerance": self.tolerance.value(),
            "max_iterations": self.max_iterations.value(),
            "max_iterations_ceiling": self.max_iterations_ceiling.value(),
            "failure_criterion": self.failure_criterion.currentData(),
            "min_slip_depth": (self.min_slip_depth.value()
                               if self.min_slip_on.isChecked()
                               and self.min_slip_depth.value() > 0 else None),
            "ssr_exclude": list(self._ssr_exclude) or None,
            "side_bc": self.side_bc.currentData(),
            "k0": self.k0.value() if self.k0_on.isChecked() else None,
            "tension_srf": self.tension_srf.isChecked(),
            "capture_failure_state": self.capture_failure_state.isChecked(),
            "capture_margin": self.capture_margin.value(),
            "capture_max_iterations": (self.capture_max_iterations.value()
                                       if self.capture_iter_on.isChecked() else None),
            # Which instant of a transient seepage solution to read (None on a model
            # with no transient solution, which is every steady model).
            "seep_time": (self.seep_time.options()
                          if self.seep_time is not None else None),
        }


#: The nodal fields a solved seepage run can be read as, by
#: :func:`~xslope.plot_seep.plot_seep_solution`'s own ``variable`` value and the name
#: the field is given wherever it is written for a reader — the seepage
#: documentation, the report's figure captions, this control. The report prints the
#: same four in the same order (``xslope.report.SEEP_PANELS``), and the two are held
#: against each other by name as well as by key.
SEEP_VARIABLES = [
    ("head", "Total head"),
    ("u", "Pore pressure"),
    ("v_mag", "Velocity magnitude"),
    ("i_mag", "Hydraulic gradient magnitude"),
]


class RunSeepDialog(QDialog):
    """Solve parameters for a seepage run. Display options (variable, contours,
    flow lines, base material, …) are not here — they live on the Seep Solution
    view and re-render the cached solution without re-solving.

    There is no BC-set selector: a steady run solves every boundary set the file
    defines — set 1 alone for most files, sets 1 and 2 for a rapid-drawdown pair
    (each keeps its own results tab and sidecar, and the only workflow that
    defines a second set needs both solutions anyway). When the loaded file
    defines a ``tseep`` (transient seepage) sheet, a Run type selector adds a
    **Transient** choice: a time-dependent solve that produces a sequence of
    saved frames. Transient always uses BC set 1 (its series bindings).

    Transient run PARAMETERS — duration, save schedule, rapid-drawdown stage times,
    and the time-series table — are model INPUTS, edited under **Inputs → Transient**
    (the :class:`~studio.editors.TransientEditor`), not here. The runner reads the
    stage times from the document's ``tseep`` data, so this dialog carries no stage
    fields; a caption points the user to the editor."""

    def __init__(self, parent=None, defaults=None, has_bc2=False, has_tseep=False,
                 slope_data=None, document=None):
        super().__init__(parent)
        self.setWindowTitle("Run Seepage")
        defaults = defaults or {}
        self.has_tseep = bool(has_tseep)
        self._sd = slope_data if slope_data is not None else {}

        layout = QVBoxLayout()                  # the left column; see two_pane
        form = QFormLayout()

        # Run type — steady vs transient. Only offered when the file has a tseep sheet;
        # otherwise the control is absent and the run is steady (back-compat).
        self.run_type = QComboBox()
        self.run_type.addItem("Steady", "steady")
        if self.has_tseep:
            self.run_type.addItem("Transient (time-dependent)", "transient")
            want = defaults.get("mode", "transient")   # default to transient when available
            j = self.run_type.findData(want)
            if j >= 0:
                self.run_type.setCurrentIndex(j)
            form.addRow("Run type", self.run_type)

        self._has_bc2 = bool(has_bc2)

        self.tol = QDoubleSpinBox()
        self.tol.setDecimals(8)
        self.tol.setRange(1e-10, 1.0)
        self.tol.setValue(float(defaults.get("tol", 1e-4)))
        form.addRow("Convergence tol", self.tol)

        self.max_iter = QSpinBox()
        self.max_iter.setRange(50, 100000)
        self.max_iter.setValue(int(defaults.get("max_iter", 400)))
        self.max_iter.setToolTip(
            "Sweep ceiling for the unconfined iteration. A run that hits it stops "
            "and reports converged = False; steep unsaturated-conductivity curves "
            "can need more than the default.")
        form.addRow("Max iterations", self.max_iter)

        layout.addLayout(form)

        # Two boundary sets means a rapid-drawdown pair; a steady run solves both,
        # each keeping its own results tab. Stated here because there is nothing
        # else in the dialog to say so.
        if self._has_bc2:
            self.bc2_caption = QLabel(
                "This file defines two boundary sets; a steady run solves both "
                "(each keeps its own results tab).")
            self.bc2_caption.setWordWrap(True)
            layout.addWidget(self.bc2_caption)

        # Where transient inputs live now (the run dialog no longer edits them).
        if self.has_tseep:
            self.transient_caption = QLabel(
                "Transient stage times and time series are edited under "
                "Inputs → Transient.")
            self.transient_caption.setWordWrap(True)
            layout.addWidget(self.transient_caption)

        # A conductivity of zero, a boundary set that drives no flow, a stale
        # {base}_seep.csv computed on another mesh: each of those runs to completion
        # today and returns an answer. They are stated here instead.
        self.preflight = PreflightPanel(
            analysis=lambda: "tseep" if self._transient() else "seep",
            slope_data=self._sd, document=document,
            selection_fn=lambda: {}, parent=self)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self._ok = bb.button(QDialogButtonBox.Ok)
        self._ok.setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)

        two_pane(self, layout, self.preflight, bb)      # controls left, checks right

        self.run_type.currentIndexChanged.connect(self._sync_mode)
        self.preflight.changed.connect(self._sync_run)
        self._sync_mode()
        self._sync_run()

    def _sync_run(self):
        blocked = self.preflight.blocked
        self._ok.setEnabled(not blocked)
        self._ok.setToolTip(self.preflight.block_reason() if blocked else "")

    def _transient(self):
        return self.has_tseep and self.run_type.currentData() == "transient"

    def _sync_mode(self, *_):
        # Transient adds its own requirements (a declared time base, storage per
        # material) on top of every steady rule, so the findings follow the mode.
        self.preflight.refresh()
        # Both solve-parameter fields belong to the STEADY unconfined iteration;
        # the transient march runs its own step and iteration controls. Dim them
        # rather than leave live controls the run ignores.
        steady = not self._transient()
        self.tol.setEnabled(steady)
        self.max_iter.setEnabled(steady)

    def options(self):
        if self._transient():
            # Stage times are NOT carried here — the runner reads them from the
            # document's tseep data (edited under Inputs → Transient).
            return {"mode": "transient", "bc": 1, "tol": self.tol.value()}
        # Steady solves every set the file defines (see the class docstring).
        return {"mode": "steady", "bc": "both" if self._has_bc2 else 1,
                "tol": self.tol.value(), "max_iter": self.max_iter.value()}


class BuildMeshDialog(QDialog):
    """Options for building a finite-element mesh from the geometry.

    Target element size is either entered directly or auto-sized as
    ``(x_max - x_min) / size_divisions`` over the ground surface (see the
    main_seep / main_fem drivers).

    The 1D element size is the size along the reinforcement and pile lines. It IS a
    model input — ``main!D20``, ``slope_data['element_size_1d']`` — so the box opens
    on whatever the file states and an entry is written back to the model, where it
    is saved and undone like any other edit. Blank means the file states nothing and
    the lines follow the mesh element size.

    The quadrilateral style (free / structured-where-possible) is a per-run
    choice: it rides on the returned options dict to ``build_mesh_from_polygons``
    and is deliberately NOT written to the model, because the .xlsx carries no
    cell for it — a mesh built at the non-default style travels as the committed
    ``{stem}_mesh.json``, not as an input. It applies only to quadrilateral
    element types and is dimmed, with the reason as its tooltip, for the
    triangular ones."""

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

        # Element size along the constraint lines — the model's own cell, so a
        # blank box means the model states nothing and the lines follow the target
        # size. A line edit rather than a spin box for exactly that reason: a spin
        # box has no empty state, and "no value" is the default here. Same idiom as
        # the local Size field on a material zone or profile line.
        self.element_size_1d = QLineEdit()
        _e1d = defaults.get("element_size_1d")
        if _e1d is not None:
            self.element_size_1d.setText(f"{float(_e1d):g}")
        self.element_size_1d.setPlaceholderText("follows the mesh element size")
        self.element_size_1d.setToolTip(ELEMENT_SIZE_1D_TIP)
        form.addRow("1D element size", self.element_size_1d)

        # Quadrilateral style — two mutually exclusive choices, so radios rather
        # than a combo: both options and the difference between them are readable
        # without opening anything. Sits with the size controls because both
        # styles are driven by the target size above.
        self.quad_style_group = QGroupBox("Quadrilateral style")
        qs_row = QHBoxLayout(self.quad_style_group)
        self.quad_style_buttons = QButtonGroup(self)
        self.quad_style_buttons.setExclusive(True)
        self._quad_style_radios = {}
        for key, label, tip in QUAD_STYLES:
            b = QRadioButton(label)
            b.setToolTip(tip)
            self.quad_style_buttons.addButton(b)
            self._quad_style_radios[key] = b
            qs_row.addWidget(b)
        qs_row.addStretch(1)
        want = defaults.get("quad_style", "free")
        self._quad_style_radios.get(want, self._quad_style_radios["free"]).setChecked(True)
        form.addRow(self.quad_style_group)

        self.refine_near_features = QCheckBox("Refine near features")
        self.refine_near_features.setChecked(bool(defaults.get("refine_near_features", False)))
        self.refine_near_features.setToolTip(
            "Shrink elements near reinforcement/pile lines, crack tips and thin "
            "material zones, growing back to the target size away from them.")
        form.addRow("", self.refine_near_features)

        self.refine_factor = QDoubleSpinBox()
        self.refine_factor.setRange(1.1, 100.0)
        self.refine_factor.setDecimals(1)
        self.refine_factor.setSingleStep(0.5)
        self.refine_factor.setValue(float(defaults.get("refine_factor", 3.0)))
        self.refine_factor.setToolTip("Local element size = target size / factor at features.")
        form.addRow("Refinement factor", self.refine_factor)

        # Thin-zone resolution, next to the feature refinement it complements. This
        # is a GUARANTEE rather than a second refinement setting: checked, a thin
        # material zone gets enough element rows to carry a shear band whatever the
        # element family; unchecked, nothing extra happens and "Refine near
        # features" is left to do whatever it would have done. It is not gated on
        # the element type, because both families have the failure and each has its
        # own cure.
        self.refine_thin_zones = QCheckBox("Refine thin zones")
        self.refine_thin_zones.setChecked(
            bool(defaults.get("refine_thin_zones", REFINE_THIN_ZONES_DEFAULT)))
        self.refine_thin_zones.setToolTip(REFINE_THIN_ZONES_TIP)
        form.addRow("", self.refine_thin_zones)

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
        self.refine_near_features.toggled.connect(self._sync_enabled)
        self.element_type.currentIndexChanged.connect(self._sync_enabled)
        self._sync_enabled()

    def quad_style(self):
        """The selected quadrilateral style, ``'free'`` or ``'structured'``."""
        for key, b in self._quad_style_radios.items():
            if b.isChecked():
                return key
        return "free"

    def _sync_enabled(self):
        auto = self.auto_size.isChecked()
        self.size_divisions.setEnabled(auto)
        self.target_size.setEnabled(not auto)
        self.refine_factor.setEnabled(self.refine_near_features.isChecked())
        # Quad styles mean nothing to a triangular mesh: dim them and say why.
        quads = str(self.element_type.currentData() or "").startswith("quad")
        self.quad_style_group.setEnabled(quads)
        tip = "" if quads else QUAD_STYLE_TRI_TIP
        self.quad_style_group.setToolTip(tip)
        for key, label, own_tip in QUAD_STYLES:
            self._quad_style_radios[key].setToolTip(own_tip if quads else tip)

    def size_1d(self):
        """The entered 1D element size as a float, or None when the box is blank
        (the lines follow the mesh element size). Raises ``ValueError`` on text
        that is not a positive number — ``accept`` reports that to the user rather
        than building a mesh at a size nobody asked for."""
        text = self.element_size_1d.text().strip()
        if not text:
            return None
        value = float(text)
        if value <= 0:
            raise ValueError(text)
        return value

    def accept(self):
        # A typo in the 1D size must not fall through to "blank": that would build a
        # mesh at the target size and look like the entry had been accepted.
        try:
            self.size_1d()
        except ValueError:
            QMessageBox.warning(
                self, "Build mesh",
                f"'{self.element_size_1d.text().strip()}' is not a valid 1D element "
                "size. Enter a positive number, or leave the box blank to mesh the "
                "reinforcement and pile lines at the mesh element size.")
            self.element_size_1d.setFocus()
            self.element_size_1d.selectAll()
            return
        super().accept()

    def options(self):
        return {
            "element_type": self.element_type.currentData(),
            "auto_size": self.auto_size.isChecked(),
            "size_divisions": self.size_divisions.value(),
            "target_size": self.target_size.value(),
            "element_size_1d": self.size_1d(),
            "quad_style": self.quad_style(),
            "refine_near_features": self.refine_near_features.isChecked(),
            "refine_factor": self.refine_factor.value(),
            "refine_thin_zones": self.refine_thin_zones.isChecked(),
        }


class RunLemDialog(QDialog):
    """Options for an LEM solve (single surface or auto-search; circular or not).

    The dialog is also the run gate. It renders the model checks
    (:class:`studio.preflight_panel.PreflightPanel`) in a column beside the run
    controls: an error disables Run and says why in the words the run would have
    been refused with, a
    warning is shown but never blocks, and methods the selected surface family
    cannot support are dimmed with the same rule's reason as their tooltip. Change
    the surface family and the method list re-filters live.

    It is also the point of use for the transient-seepage times: the seepage-time
    selector when the model carries a transient solution, and the rapid-drawdown
    stage times, which are edited here and stored on the ``tseep`` sheet.
    """

    def __init__(self, parent=None, defaults=None, slope_data=None, document=None,
                 transient=None, current_time=None):
        super().__init__(parent)
        self.setWindowTitle("Run LEM")
        defaults = defaults or {}
        slope_data = slope_data if slope_data is not None else {}
        self._sd = slope_data

        layout = QVBoxLayout()                  # the left column; see two_pane
        form = QFormLayout()

        self.method = self._combo(LEM_METHODS, defaults.get("method", "spencer"))
        form.addRow("Method", self.method)

        self.analysis = self._combo(ANALYSIS_TYPES, defaults.get("analysis", "auto_search"))
        form.addRow("Analysis", self.analysis)

        note = QLabel("Single surface analyzes the first circle / the non-circular "
                      "surface as entered. Auto-search refines from there to the "
                      "critical surface.")
        note.setWordWrap(True)
        form.addRow(note)

        # The surface-family selector appears when the deck carries BOTH families —
        # keyed on what is present, never on the stored selection, because reading the
        # selection as presence would make a run that chose non-circular unable to
        # choose the circles back. With one family present it is a fixed label.
        has_circ = bool(slope_data.get("circles"))
        has_ncirc = bool(slope_data.get("non_circ"))
        if has_circ != has_ncirc:                      # exactly one defined
            self._fixed_surface = "circular" if has_circ else "noncircular"
            self.surface = None
            form.addRow("Surface", QLabel(dict(SURFACE_TYPES)[self._fixed_surface]))
        else:                                          # both (a real choice) or neither
            self._fixed_surface = None
            stored = (stored_surface_family(slope_data)
                      if (has_circ and has_ncirc) else "circular")
            self.surface = self._combo(SURFACE_TYPES,
                                       defaults.get("surface", stored))
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
        # Toggling it changes the run the checks are describing, so re-ask them.
        self.grid_seed.toggled.connect(lambda *_: self._recheck())
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

        # Transient seepage: which instant a single-time run reads, and — for a rapid
        # drawdown — which two instants the stages read. Both are shown only when the
        # model has something time-dependent to say; a steady model has one field.
        self._has_transient = bool((transient or {}).get("times"))
        self.seep_time = None
        if transient or slope_data.get("tseep"):
            self.seep_time = SeepageTimeSelector(self, transient=transient,
                                                 current_time=current_time,
                                                 slope_data=slope_data)
            layout.addWidget(self.seep_time)
        self.stage_times = None
        if slope_data.get("tseep") or transient:
            self.stage_times = StageTimeFields(self, slope_data=slope_data)
            layout.addWidget(self.stage_times)

        # The model checks, beside the controls: warnings inform the decision instead
        # of annotating the regret, and an error refuses before anything is started.
        self.preflight = PreflightPanel(
            analysis=lambda: "rapid" if self.rapid.isChecked() else "lem",
            slope_data=slope_data, document=document,
            selection_fn=self._selection, notes=(SEISMIC_NOTE_LEM,), parent=self)

        self.analysis.currentIndexChanged.connect(self._sync_tols)
        if self.surface is not None:
            self.surface.currentIndexChanged.connect(self._sync_tols)
        self.min_slip_on.toggled.connect(self._sync_tols)
        self.method.currentIndexChanged.connect(self._recheck)
        self.rapid.toggled.connect(self._recheck)
        self.rapid.toggled.connect(self._sync_seep_time)
        # The seepage time is part of what the checks are told (it names the frame
        # that will be staged), so a change to it re-evaluates them as well as the
        # button — otherwise the panel would go on naming the previous instant.
        if self.seep_time is not None:
            self.seep_time.mode.currentIndexChanged.connect(self._recheck)
            self.seep_time.frame.currentIndexChanged.connect(self._recheck)
            self.seep_time.other.textChanged.connect(self._recheck)
        if self.stage_times is not None:
            self.stage_times.stage_1.textChanged.connect(self._recheck)
            self.stage_times.stage_2.textChanged.connect(self._recheck)
        self._sync_tols()

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self._ok = bb.button(QDialogButtonBox.Ok)
        self._ok.setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)

        two_pane(self, layout, self.preflight, bb)      # controls left, checks right

        self.preflight.changed.connect(self._sync_run)
        self._sync_seep_time()

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

    def _selection(self):
        """What preflight and capabilities() are asked about — read live."""
        return {"method": self.method.currentData(),
                "surface": self._surface_value(),
                "search": self.analysis.currentData() == "auto_search",
                # Grid seeding changes which mechanisms a search reaches, so the
                # checks are handed the whole run description, seeding included.
                "grid_seed": (self.grid_seed.isChecked()
                              if hasattr(self, "grid_seed") else False),
                "seep_frame": self._seep_frame()}

    def _seep_frame(self):
        """The transient frame(s) this run stages, as the model checks are told it.

        With ``u = seep`` and no stored field the model is only incomplete if nothing
        is going to supply one — and the Seepage time group directly above the checks
        is exactly what supplies it: on OK, ``MainWindow._apply_transient_analysis_frame``
        writes the chosen frame into ``seep_u`` (or refuses the run, loudly) before
        the solver is built. That is the same order a scripted run follows by calling
        ``apply_transient_stability_frame`` before the entry point whose gate then
        sees the field. This says so, so the gate can agree with both.

        ``None`` when no transient solution is loaded, which is when the
        missing-field ERROR is right and must stand.
        """
        if not self._has_transient:
            return None
        if self.rapid.isChecked() and self.stage_times is not None:
            s1, s2 = self.stage_times.values()
            if not (s1 is None and s2 is None):        # staging asked for
                if self.stage_times.problem() is not None:
                    return {"times": []}
                return {"times": [float(s1), float(s2)]}
        # A single instant — including a rapid drawdown with both stage times blank,
        # which is what the staging call itself falls back to.
        return (self.seep_time.frame_selection() if self.seep_time is not None
                else {"times": []})

    def _sync_methods(self):
        """Dim every method the selected surface family cannot support.

        The reason string is the rule's own, so the tooltip on a dimmed item and
        the refusal the run would have printed are one sentence.
        """
        from xslope.preflight import capabilities
        # Read the model from the panel: a remedy applied from the findings list
        # replaces the model's contents, and the dimming must follow it.
        caps = capabilities(self.preflight.model, self._selection())
        apply_capabilities(self.method, caps.get("lem_method", {}))

    def _recheck(self):
        self.preflight.refresh()

    def _sync_seep_time(self):
        """A rapid drawdown reads its two stage times, not a single instant, so the
        single-time selector steps aside while Rapid drawdown is ticked and the stage
        fields become the live control."""
        rapid = self.rapid.isChecked()
        if self.seep_time is not None:
            self.seep_time.setVisible(not rapid)
        if self.stage_times is not None:
            self.stage_times.setVisible(rapid)
        self._sync_run()

    def _sync_run(self):
        """Run refuses on a preflight ERROR, and on times it cannot read."""
        reason = self.preflight.block_reason() if self.preflight.blocked else None
        if reason is None:
            if self.rapid.isChecked():
                if self.stage_times is not None:
                    reason = self.stage_times.problem()
            elif self.seep_time is not None:
                reason = self.seep_time.problem()
        self._ok.setEnabled(reason is None)
        self._ok.setToolTip(reason or "")

    def _sync_tols(self):
        self._sync_methods()
        is_search = self.analysis.currentData() == "auto_search"
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
        # meaning for the auto-search (single-surface just evaluates the given surface).
        self.min_slip_on.setEnabled(is_search)
        self.min_slip_depth.setEnabled(is_search and self.min_slip_on.isChecked())
        # The findings depend on the surface family and on whether this is a search
        # (a non-circular SEARCH with OMS is refused in different words from a single
        # surface), so they are re-evaluated here rather than only at open.
        self._recheck()

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
            # Transient seepage. A rapid drawdown reads the two stage times instead
            # of a single instant, so only one of these is ever live.
            "seep_time": (None if self.rapid.isChecked() or self.seep_time is None
                          else self.seep_time.options()),
            "stage_times": (self.stage_times.options()
                            if self.rapid.isChecked() and self.stage_times is not None
                            else None),
        }


class SensitivityDialog(QDialog):
    """Options for a Parametric study (one dialog, three modes).

    - **Sensitivity** sweeps several parameters, each +/- a percent about its
      current value (per-row overridable, with a one-click sigma-range preset that
      reuses the reliability sigma_* columns). A **Plot type** selector chooses the
      view: a tornado (default), scaled-sensitivity bars (elasticity / per-1% /
      per-sigma), a spider plot, a variance-contribution Pareto, or Monte Carlo
      rank-correlation bars — the last two only when the model carries sigmas.
    - **Design** sweeps ONE parameter across explicit From/To bounds in N steps and
      finds where the output meets a target (FS = 1.5 by default).
    - **Back-Analysis** is the same single-parameter sweep framed as a failure
      investigation: FS = 1.0 is known, so it back-calculates the parameter value
      consistent with the observed slide.
    - **Factor of safety vs time** sweeps the saved instants of a transient
      seepage march instead of a parameter: no input changes, and each point
      solves the same model against that instant's pore pressures. It needs a
      transient solution in hand and a material reading ``u = seep``, and says so
      when either is missing.

    The sweepable set is any numeric material property plus the global k_seismic;
    geometry design stays in main_design.py. A thin caller of ``xslope.sensitivity``:
    ``options()`` returns plain dicts/lists the runner forwards straight to
    ``design()`` / ``back_analysis()`` / ``sensitivity()`` and the plot data
    functions.
    """

    _GLOBAL_LABEL = "k_seismic (global)"
    # per-engine-mode output quantity: (short name, long label, design-target label)
    _OUTPUT = {
        "lem": ("FS", "Factor of Safety", "Target FS"),
        "fem": ("FS", "Factor of Safety", "Target FS"),
        "seep": ("q", "total discharge q", "Target q"),
    }

    def __init__(self, parent=None, defaults=None, slope_data=None, app_mode="lem",
                 document=None, transient=None):
        super().__init__(parent)
        self.app_mode = app_mode if app_mode in self._OUTPUT else "lem"
        out_short, out_long, target_label = self._OUTPUT[self.app_mode]
        self._out_short = out_short
        titles = {"lem": "Parametric study (LEM)",
                  "fem": "Parametric study (FEM · SSRM)",
                  "seep": "Parametric study (Seepage)"}
        self.setWindowTitle(titles[self.app_mode])
        defaults = defaults or {}
        slope_data = slope_data or {}

        from xslope.sensitivity import list_params
        self._params = list_params(slope_data, mode=self.app_mode)
        self._by_ref = {e["ref"]: e for e in self._params}
        # Sigma-gated plots (variance Pareto, MC rank, per-sigma scaling) are only
        # offered when the model carries at least one reliability standard deviation.
        self._has_sigma = any(e.get("sigma") for e in self._params)
        # Remembered so a Design<->Back-Analysis switch can restore the design target.
        self._design_target_default = float(
            defaults.get("target_fs", 1.5 if self.app_mode != "seep" else 1e-5))
        self._prev_mode = None
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

        layout = QVBoxLayout()                  # the left column; see two_pane
        form = QFormLayout()

        # FS-versus-time sweeps the frames of a transient march. Two things have to
        # be true for it to mean anything, and each has its own plain sentence: a
        # march to read, and a material that reads it.
        self._times = [float(t) for t in (transient or {}).get("times", [])]
        self._time_unit = str(slope_data.get("time_unit") or "").strip()
        # Grid seeding sweeps circle centers, so it has nothing to offer a model
        # whose surface is non-circular.
        self._circular = bool(slope_data.get("circular", True))
        self._fs_time_reason = self._fs_time_unavailable(slope_data)
        # The instant a drawdown curve reads stage 1 at: the file's own stage_1
        # (normally 0, full pool), or the earliest saved frame where it names none.
        s1 = (slope_data.get("tseep") or {}).get("stage_1")
        self._stage_1 = (float(s1) if s1 is not None
                         else (min(self._times) if self._times else 0.0))
        self._rapid_time_reason = self._rapid_time_unavailable(slope_data)

        ba_label = ("Back-Analysis (FS = 1)" if self.app_mode != "seep"
                    else f"Back-Analysis (target {out_short})")
        self.mode = self._combo([("sensitivity", "Sensitivity (tornado + plots)"),
                                 ("design", f"Design ({out_short} target)"),
                                 ("back_analysis", ba_label),
                                 ("fs_vs_time", "Factor of safety vs time")],
                                defaults.get("mode", "sensitivity"))
        if self._fs_time_reason:
            item = self.mode.model().item(self.mode.count() - 1)
            if item is not None:
                item.setEnabled(False)
                item.setToolTip(self._fs_time_reason)
            if self.mode.currentData() == "fs_vs_time":
                self.mode.setCurrentIndex(0)
        form.addRow("Mode", self.mode)

        # Engine-specific solver row(s). LEM keeps method + slices; FEM swaps in the
        # SSRM knobs (each step is a full SSRM solve); Seep takes a BC set + tol.
        self.method = self._combo(LEM_METHODS, defaults.get("method", "spencer"))
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
        self.picker = QGroupBox("Parameter")
        pform = QFormLayout(self.picker)
        self.material = QComboBox()
        for disp, key in self._groups:
            self.material.addItem(disp, key)
        pform.addRow("Material" if self.app_mode != "seep" else "Material / BC",
                     self.material)
        self.prop = QComboBox()
        pform.addRow("Property", self.prop)
        self.material.currentIndexChanged.connect(self._on_material_changed)
        self.prop.currentIndexChanged.connect(self._on_prop_changed)
        layout.addWidget(self.picker)

        # --- mode pages -----------------------------------------------------
        self.stack = QStackedWidget()
        self.stack.addWidget(self._build_sens_page(defaults))
        self.stack.addWidget(self._build_design_page(defaults, target_label))
        self.stack.addWidget(self._build_time_page(defaults))
        layout.addWidget(self.stack)

        # Re-search toggle is an LEM concept (FEM finds its own mechanism, seepage
        # has no failure surface) — only shown in LEM mode.
        self.search = QCheckBox("Re-search the critical surface at each step")
        self.search.setChecked(bool(defaults.get("search", True)))
        self._search_tip = (
            "On (recommended): re-run the search at every swept value, because the "
            "critical surface moves as the parameter changes — a fixed surface "
            "silently understates the sensitivity.\n\n"
            "Off: re-solve the entered surface only (much faster, but the answer is "
            "only right for that prescribed surface).")
        self.search.setToolTip(self._search_tip)
        self.search.setVisible(self.app_mode == "lem")
        self.search.toggled.connect(lambda *_: self._sync_grid_seed())
        layout.addWidget(self.search)

        self.note = QLabel()
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

        # A sweep is a composite analysis: it inherits every rule of the base
        # analysis it repeats, so the base model is checked once here rather than
        # per step (a swept value is a deliberate perturbation, not a mistake).
        self.preflight = PreflightPanel(
            analysis="sensitivity", slope_data=slope_data, document=document,
            selection_fn=self._selection,
            notes=(SEISMIC_NOTE_LEM if self.app_mode == "lem" else SEISMIC_NOTE_FEM,),
            parent=self)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self._ok = bb.button(QDialogButtonBox.Ok)
        self._ok.setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)

        two_pane(self, layout, self.preflight, bb)

        self.preflight.changed.connect(self._sync_run)
        if self.app_mode == "lem":
            self.method.currentIndexChanged.connect(self.preflight.refresh)
        self._sync_run()

        self.mode.currentIndexChanged.connect(self._on_mode_changed)
        # FS-versus-time names a frame the checks have to be told about, so a mode
        # switch re-evaluates them rather than leaving the previous mode's answer.
        self.mode.currentIndexChanged.connect(self.preflight.refresh)
        self.plot_type.currentIndexChanged.connect(self._on_plot_type_changed)
        self._on_material_changed()                 # populate property combo
        # Restore a remembered sensitivity table (from the previous run this session).
        for spec in defaults.get("_remember_params", []):
            e = self._by_ref.get(spec.get("ref"))
            if e is not None:
                self._add_row(e, pct=spec.get("pct", self.pct.value()),
                              use_sigma=spec.get("use_sigma", False))
        self._on_plot_type_changed()
        self._on_rapid_toggled()            # a remembered drawdown holds Re-search
        self._on_mode_changed()
        self.resize(1000, 620)              # two columns: controls | model checks

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
        # A sweep tracks ONE discharge per step, so with two boundary sets the set
        # is a real choice; with one there is nothing to pick and the row is noise
        # (the hidden combo still answers options() with set 1).
        if self.seep_bc.count() > 1:
            form.addRow("BC set", self.seep_bc)
        else:
            self.seep_bc.setVisible(False)
        form.addRow("Convergence tol", self.seep_tol)
        note = QLabel("Output quantity = total discharge q through the section. Sweep "
                      "hydraulic conductivity / boundary heads; each step is one "
                      "seepage solve.")
        note.setWordWrap(True)
        form.addRow("", note)

    # --- page builders ------------------------------------------------------
    def _plot_type_items(self):
        """Plot-type menu, sigma-gated: the variance Pareto and MC rank plots need
        reliability sigmas, so they appear only when the model carries one."""
        items = [("tornado", "Tornado (FS swing per parameter)"),
                 ("scaled", "Scaled-sensitivity bars"),
                 ("spider", "Spider (FS vs each parameter)")]
        # Variance Pareto and MC rank reuse the LEM Taylor-series / Monte-Carlo
        # reliability, so they are offered only for an LEM study that carries sigmas.
        if self._has_sigma and self.app_mode == "lem":
            items += [("variance", "Variance Pareto (σ)"),
                      ("rank", "Monte Carlo rank correlation (σ)")]
        return items

    def _scaling_items(self):
        items = [("elasticity", "Elasticity (∂F/∂p · p/F)"),
                 ("per_1pct", "Per 1% change")]
        if self._has_sigma:
            items.append(("per_sigma", "Per σ"))
        return items

    def _build_sens_page(self, defaults):
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(0, 0, 0, 0)

        # Plot-type selector: the view the sweep produces. Scaled bars expose a
        # scaling sub-choice; the MC rank plot exposes a sample count.
        pt = QHBoxLayout()
        pt.addWidget(QLabel("Plot type"))
        self.plot_type = self._combo(self._plot_type_items(),
                                     defaults.get("plot_type", "tornado"))
        pt.addWidget(self.plot_type)
        pt.addSpacing(10)
        self._scaling_label = QLabel("Scaling")
        pt.addWidget(self._scaling_label)
        self.scaling = self._combo(self._scaling_items(),
                                   defaults.get("scaling", "elasticity"))
        pt.addWidget(self.scaling)
        pt.addSpacing(10)
        self._mc_label = QLabel("MC samples")
        pt.addWidget(self._mc_label)
        self.mc_samples = QSpinBox()
        self.mc_samples.setRange(200, 200000)
        self.mc_samples.setSingleStep(1000)
        self.mc_samples.setValue(int(defaults.get("mc_samples", 5000)))
        pt.addWidget(self.mc_samples)
        pt.addStretch(1)
        v.addLayout(pt)

        controls = QHBoxLayout()
        self._pct_label = QLabel("Default ±%")
        controls.addWidget(self._pct_label)
        self.pct = QDoubleSpinBox()
        self.pct.setRange(1.0, 99.0)
        self.pct.setDecimals(0)
        self.pct.setSuffix(" %")
        self.pct.setValue(float(defaults.get("default_pct", 20.0)))
        controls.addWidget(self.pct)
        controls.addSpacing(12)
        self._points_label = QLabel("Points")
        controls.addWidget(self._points_label)
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

    def _fs_time_unavailable(self, slope_data):
        """Why an FS-versus-time run cannot be made on this model, or ``""``.

        Three conditions, each with the sentence a user can act on. The engine mode
        comes first because it is the one the mode strip decides: a seepage sweep
        reports discharge, and the seepage solution is this run's INPUT.
        """
        if self.app_mode == "seep":
            return ("A factor-of-safety curve needs a stability engine. Switch the "
                    "mode strip to LEM or FEM — the seepage solution is this run's "
                    "input, not its output.")
        if not self._times:
            return "Run a transient seepage analysis first."
        uses_seep = any(str(m.get("u", "")).strip().lower() == "seep"
                        for m in (slope_data.get("materials") or []))
        if not uses_seep:
            return ("No material takes its pore pressure from the seepage solution. "
                    "Set u = seep on the materials table, or the curve would be flat.")
        return ""

    @staticmethod
    def _rapid_time_unavailable(slope_data):
        """Why the instants cannot be evaluated as rapid drawdowns, or ``""``.

        The three-stage procedure is a comparison between undrained and drained
        strengths on the drawn-down section, and both come from the d / psi columns.
        With neither set anywhere, every stage reads the same strengths and the
        drawdown is a single-stage analysis wearing three names — so the option is
        offered only on a model that carries them.
        """
        for m in (slope_data.get("materials") or []):
            for field in ("d", "psi"):
                try:
                    if float(m.get(field) or 0.0) > 0.0:
                        return ""
                except (TypeError, ValueError):
                    continue
        return ("No material carries the rapid-drawdown strengths d and "
                "\N{GREEK SMALL LETTER PSI}, so all three stages would read the "
                "same strengths. Set them on the low-permeability material(s) in "
                "the materials table.")

    def _build_time_page(self, defaults):
        """The instants to evaluate: every saved frame of the march, all ticked.

        A checklist rather than a range, because the times are not a continuum — a
        curve can only be drawn through instants the march SAVED, and a time between
        two frames is never interpolated. Unticking is how a long march is sampled
        (each instant is a full stability run) or how a leading frame is dropped.
        """
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(0, 0, 0, 0)

        row = QHBoxLayout()
        row.addWidget(QLabel("Saved frames"))
        row.addStretch(1)
        self.times_all = QPushButton("All")
        self.times_all.setToolTip("Tick every saved instant.")
        self.times_all.clicked.connect(lambda: self._set_all_times(True))
        self.times_none = QPushButton("None")
        self.times_none.setToolTip("Untick every saved instant.")
        self.times_none.clicked.connect(lambda: self._set_all_times(False))
        row.addWidget(self.times_all)
        row.addWidget(self.times_none)
        v.addLayout(row)

        self.times = QListWidget(page)
        self.times.setToolTip("The instants the transient run saved. Each ticked "
                              "one is a full stability run.")
        unit = f" {self._time_unit}" if self._time_unit else ""
        remembered = defaults.get("times")
        for t in self._times:
            item = QListWidgetItem(f"t = {_fmt_time(t)}{unit}")
            item.setData(Qt.UserRole, float(t))
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            ticked = (True if remembered is None
                      else any(abs(float(r) - t) < 1e-9 for r in remembered))
            item.setCheckState(Qt.Checked if ticked else Qt.Unchecked)
            self.times.addItem(item)
        self.times.itemChanged.connect(lambda *_: self._on_times_changed())
        v.addWidget(self.times, 1)

        # Each instant as a drawdown rather than a single state. It sits with the
        # frames because it changes what a ticked instant MEANS: stage 2 of a fall
        # that started at the march's initial pool, not a state on its own.
        unit_txt = f" {self._time_unit}" if self._time_unit else ""
        self.rapid = QCheckBox("Rapid drawdown at each time")
        self.rapid.setToolTip(self._rapid_time_reason or (
            f"On: every ticked instant is a three-stage Duncan-Wright-Brandon "
            f"rapid drawdown — stage 1 the transient run's initial state at "
            f"t = {_fmt_time(self._stage_1)}{unit_txt} (full pool), stage 2 that "
            f"instant's drawn-down state, stage 3 the same section re-checked with "
            f"drained strengths. The curve reports the drawdown's own factor of "
            f"safety, the lower of stages 2 and 3.\n\n"
            f"Off: every instant is a single-stage analysis against that instant's "
            f"pore pressures."))
        if self._rapid_time_reason:
            self.rapid.setEnabled(False)
        else:
            self.rapid.setChecked(bool(defaults.get("rapid", False)))
        self.rapid.toggled.connect(self._on_rapid_toggled)
        v.addWidget(self.rapid)

        # Grid seeding, the same option Run LEM offers, and the one this mode most
        # often needs: a curve spans states whose critical mechanisms are in
        # different places, and a search refined from one starting circle stays in
        # that circle's neighborhood at every instant. Same label and tooltip as Run
        # LEM's, so the two dialogs name one option once.
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
        self.grid_seed.toggled.connect(lambda *_: self._on_grid_seed_toggled())
        v.addWidget(self.grid_seed)
        return page

    def _on_grid_seed_toggled(self, *_):
        """The note says how each point is searched and the checks read the same
        flag out of ``_selection()``, so both are re-asked when the box moves."""
        self._on_mode_changed()
        if hasattr(self, "preflight"):
            self.preflight.refresh()

    def _on_rapid_toggled(self, *_):
        """A drawdown point is always searched, so the Re-search toggle has nothing
        left to decide: it is held on and greyed with the reason, rather than
        offering a run this mode does not make."""
        on = self.rapid.isChecked()
        if hasattr(self, "search"):
            if on:
                self.search.setChecked(True)
                self.search.setEnabled(False)
                self.search.setToolTip(
                    "Held on: a rapid drawdown's critical surface is not the "
                    "drained one and moves as the drawn-down field changes, so "
                    "every instant is searched from the starting circle.")
            else:
                self.search.setEnabled(True)
                self.search.setToolTip(self._search_tip)
        self._sync_stage_1_frames(on)
        self._sync_grid_seed()
        if self.mode.currentData() == "fs_vs_time":
            self._on_mode_changed()
            self.preflight.refresh()

    def _sync_stage_1_frames(self, rapid_on):
        """With every instant a drawdown, the frames at or before stage 1 are the
        state the others fall from and cannot be a stage 2 themselves; they are
        unticked and dimmed rather than left to come back as failed rows."""
        stage_1 = getattr(self, "_stage_1", None)
        if stage_1 is None:
            return
        unit = f" {self._time_unit}" if self._time_unit else ""
        for i in range(self.times.count()):
            item = self.times.item(i)
            if float(item.data(Qt.UserRole)) > stage_1:
                continue
            if rapid_on:
                item.setCheckState(Qt.Unchecked)
                item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
                item.setToolTip(f"Stage 1 of every drawdown is read at "
                                f"t = {_fmt_time(stage_1)}{unit}; this frame is "
                                f"the state the others fall from, not a "
                                f"drawdown of its own.")
            else:
                item.setFlags(item.flags() | Qt.ItemIsEnabled)
                item.setCheckState(Qt.Checked)
                item.setToolTip("")

    def _set_all_times(self, on):
        for i in range(self.times.count()):
            item = self.times.item(i)
            if not (item.flags() & Qt.ItemIsEnabled):
                continue
            item.setCheckState(Qt.Checked if on else Qt.Unchecked)

    def selected_times(self):
        """The ticked instants, in saved order."""
        return [float(self.times.item(i).data(Qt.UserRole))
                for i in range(self.times.count())
                if self.times.item(i).checkState() == Qt.Checked]

    def _on_times_changed(self):
        if self.mode.currentData() == "fs_vs_time":
            self._on_mode_changed()
            self.preflight.refresh()

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

    def _selection(self):
        """What preflight and capabilities() are asked about — read live.

        In FS-versus-time mode the selection also names the frame that WILL be
        staged: ``fs_vs_time`` gates the base model with its first instant already
        placed, so a model whose materials read ``u = seep`` is not missing a
        pore-pressure field here, and the checks must not refuse the one run the
        transient controls exist to start.
        """
        sel = {"base": self.app_mode,
               "method": self.method.currentData(),
               "search": self.search.isChecked()}
        if self.mode.currentData() == "fs_vs_time":
            picked = self.selected_times()
            # Same key Run LEM passes, so both dialogs describe the run the
            # checks are asked about the same way.
            sel["grid_seed"] = (self.grid_seed.isChecked()
                                if hasattr(self, "grid_seed") else False)
            if self.rapid.isChecked():
                # Every point is a drawdown, so the checks are the drawdown's:
                # base = 'rapid' brings that family in on top of the LEM rules, and
                # the two instants say which route it runs (both from the march).
                sel["base"] = "rapid"
                sel["rapid"] = True
                later = [t for t in picked if t > self._stage_1]
                sel["seep_frame"] = {"times": [self._stage_1,
                                               (later[-1] if later
                                                else (picked[-1] if picked
                                                      else self._stage_1))]}
            else:
                sel["seep_frame"] = {"times": picked[:1]}
        return sel

    def _sync_run(self):
        """Run refuses on an ERROR in the base model; a warning never blocks."""
        if self.app_mode == "lem":
            from xslope.preflight import capabilities
            apply_capabilities(self.method,
                               capabilities(self.preflight.model,
                                            {"method": self.method.currentData(),
                                             "search": self.search.isChecked()}
                                            ).get("lem_method", {}))
        reason = self.preflight.block_reason() if self.preflight.blocked else None
        if reason is None and self.mode.currentData() == "fs_vs_time":
            if not self.selected_times():
                reason = "Tick at least one saved instant to evaluate."
        self._ok.setEnabled(reason is None)
        self._ok.setToolTip(reason or "")

    def _on_mode_changed(self):
        m = self.mode.currentData()
        single = m in ("design", "back_analysis")
        is_time = m == "fs_vs_time"
        self.stack.setCurrentIndex(2 if is_time else (1 if single else 0))
        # No input is substituted at any point of an FS-versus-time run — the axis
        # is time — so the parameter picker has nothing to say and steps aside.
        self.picker.setVisible(not is_time)
        if is_time:
            n = len(self.selected_times())
            unit = f" {self._time_unit}" if self._time_unit else ""
            seed_txt = ("over a grid of centers and tangent elevations"
                        if getattr(self, "grid_seed", None)
                        and self.grid_seed.isChecked()
                        else "from the starting circles")
            if self.rapid.isChecked():
                self.note.setText(
                    f"Rapid drawdown vs time: one three-stage drawdown per ticked "
                    f"instant of the transient seepage run — {n} of "
                    f"{len(self._times)} saved frames. Stage 1 is that run's "
                    f"initial state at t = {_fmt_time(self._stage_1)}{unit} and "
                    f"stage 2 is that instant, so the curve answers how safe the "
                    f"slope is if the pool falls to where it stands then. Each "
                    f"point is searched {seed_txt}, and the reported "
                    f"value is the lower of stages 2 and 3.")
            else:
                self.note.setText(
                    f"Factor of safety vs time: one stability run per ticked instant "
                    f"of the transient seepage run — {n} of {len(self._times)} "
                    f"saved frames. No input changes between the points; each "
                    f"solves the same model against that instant's pore pressures, "
                    f"and the reservoir load is re-derived from the pool as it stood "
                    f"then. Each point is searched {seed_txt}, and the circles "
                    f"sheet's search window is applied.")
            self._prev_mode = m
            self._sync_grid_seed()
            self._sync_run()
            return
        # Seed the target on an actual transition into a mode; FS = 1.0 is the
        # convention for a back-analysis, the design default otherwise. Seep sweeps
        # target a discharge q, so leave their target alone.
        if self.app_mode != "seep":
            if m == "back_analysis" and self._prev_mode != "back_analysis":
                self.d_target.setValue(1.0)
            elif m == "design" and self._prev_mode == "back_analysis":
                self.d_target.setValue(self._design_target_default)
        self._prev_mode = m
        q = self._out_short
        if m == "back_analysis":
            self.note.setText(
                "Back-Analysis: a failure occurred, so the factor of safety at "
                "failure is known to be 1.0. Sweep the one parameter above and read "
                "off the value consistent with the observed slide (the back-calculated "
                f"value at {q} = 1). Widen [From, To] if the curve never reaches it.")
        elif m == "design":
            self.note.setText(
                f"Design: sweeps the one parameter above across [From, To] and "
                f"annotates where {q} meets the target (interpolated). If it never "
                f"crosses, the plot says which way to widen the range.")
        else:
            self._on_plot_type_changed()

    def _sync_grid_seed(self):
        """Grid seeding drives a circular SEARCH, so it is live only when this run
        makes one — LEM, a circular model, and either the re-search toggle on or a
        rapid drawdown, which holds it on. Mirrors Run LEM's own rule."""
        if not hasattr(self, "grid_seed"):
            return
        searching = self.rapid.isChecked() or self.search.isChecked()
        live = (self.app_mode == "lem" and self._circular and searching
                and self.mode.currentData() == "fs_vs_time")
        self.grid_seed.setEnabled(live)
        if not live:
            self.grid_seed.setChecked(False)

    def _on_plot_type_changed(self):
        """Toggle the scaling / MC-sample controls and the sensitivity note to match
        the selected plot type (no effect while a single-parameter mode is active)."""
        pt = self.plot_type.currentData()
        scaled = pt == "scaled"
        rank = pt == "rank"
        for w in (self._scaling_label, self.scaling):
            w.setVisible(scaled)
        for w in (self._mc_label, self.mc_samples):
            w.setVisible(rank)
        # Gate the sweep controls by what the selected plot type actually reads:
        # the sigma-based plots (variance, rank) never read the table, and the
        # scaled bars read only its parameter list (the derivative is a fixed
        # ±1% central difference), so the ±% ranges and Points gray out with
        # them. Anything a selection ignores is disabled, not left live.
        mode_sens = self.mode.currentData() == "sensitivity"
        table_used = (not mode_sens) or pt not in ("variance", "rank")
        ranges_used = (not mode_sens) or pt in ("tornado", "spider")
        for w in (self.material, self.prop, self.add_btn, self.table):
            w.setEnabled(table_used)
        for w in (self._pct_label, self.pct, self._points_label, self.n_points):
            w.setEnabled(ranges_used)
        for r in self._rows:
            r["sigma_btn"].setEnabled(ranges_used and r.get("sigma") is not None)
            r["pct_spin"].setEnabled(ranges_used
                                     and not r["sigma_btn"].isChecked())
        if not mode_sens:
            return
        q = self._out_short
        notes = {
            "tornado": (f"Tornado: each table parameter swept ±% about its value "
                        f"(or ±σ); bars show the {q} swing, widest on top. Double-click "
                        f"a bar for that parameter's curve."),
            "scaled": ("Scaled-sensitivity bars: one bar per table parameter, height = "
                       "the chosen scaling of ∂F/∂p (central difference at a fixed ±1%, "
                       "so the ±% ranges and Points do not apply), color = sign. "
                       "Elasticity is unitless and comparable across parameters."),
            "spider": (f"Spider: {q} vs each table parameter over its ±% range, on one "
                       f"normalized axis (% change from base), with a base-case marker."),
            "variance": ("Variance Pareto: each uncertain parameter's share of Var(FS) "
                         "from the Taylor-series reliability, sorted tallest first. "
                         "Uses every σ-carrying material — the sweep table does not "
                         "apply and is disabled."),
            "rank": ("Monte Carlo rank correlation: Spearman correlation of each sampled "
                     "input with FS — a GLOBAL measure. Uses every σ-carrying material — "
                     "the sweep table does not apply and is disabled. Can take a while "
                     "at high sample counts."),
        }
        self.note.setText(notes.get(pt, ""))

    # --- sensitivity table --------------------------------------------------
    def _on_add_clicked(self):
        ref = self.prop.currentData()
        e = self._by_ref.get(ref)
        if e is None:
            return
        if any(r["ref"] == ref for r in self._rows):
            return                                   # already added
        self._add_row(e, pct=self.pct.value(), use_sigma=False)
        self._on_plot_type_changed()

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
        if common["mode"] == "fs_vs_time":
            common["times"] = self.selected_times()
            common["rapid"] = self.rapid.isChecked()
            common["grid_seed"] = self.grid_seed.isChecked()
            if common["rapid"]:
                common["search"] = True
        elif common["mode"] in ("design", "back_analysis"):
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
                           # plot-type selection for the sensitivity view
                           "plot_type": self.plot_type.currentData(),
                           "scaling": self.scaling.currentData(),
                           "mc_samples": self.mc_samples.value(),
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


RELIABILITY_ENGINES = [("taylor", "Taylor series (TSPM)"),
                       ("mc", "Monte Carlo"),
                       ("rs", "Response surface (RS)")]
#: The reliability engines that sample the standard deviations rather than
#: perturbing them. Both are limit-equilibrium only and share the seed and
#: distribution controls.
RELIABILITY_SAMPLING_ENGINES = ("mc", "rs")
MC_DISTRIBUTIONS = [("normal", "Normal"), ("lognormal", "Lognormal")]


SEARCH_LABEL_TAYLOR = "Search for the critical surface at each solve (MLV and MLV ± σ)"
SEARCH_LABEL_SAMPLING = "Search for the critical surface at the mean values"
#: The unticked alternative, stated under the checkbox — it is the genuinely
#: non-obvious half, and it differs by engine only in what happens around it:
#: either way the sampling engines hold ONE surface for every realization.
SEARCH_HINT_TAYLOR = ("Unticked: every solve evaluates the first circle on the "
                      "circles sheet (or the non-circular surface) — the same "
                      "surface a Single surface run solves.")
SEARCH_TIP_TAYLOR = (
    "On: the Taylor series solves the model 1+2N times — once at the "
    "most-likely values and at MLV ± σ for each uncertain parameter — and "
    "every one of those solves is a fresh search, so each perturbed model "
    "finds its own critical surface. On some models all the searches land on "
    "one circle; on layered models they separate, and that spread is part of "
    "the sensitivity.\n\n"
    "Off: every solve evaluates the first circle on the circles sheet (or the "
    "non-circular surface), exactly as a Single surface run does. Untick when "
    "the surface itself is prescribed — a published benchmark's circle, an "
    "observed failure surface, a geologically controlled plane — and the "
    "question is that surface's own reliability.")
SEARCH_TIP_SAMPLING = (
    "On: the critical surface is searched once, at the most-likely values, "
    "and every realization is evaluated on it — one surface for the whole "
    "campaign, because re-searching per realization is prohibitive at "
    "thousands of solves. This fixed-surface convention is what the "
    "commercial codes use in their probabilistic modes.\n\n"
    "Off: the fixed surface is the first circle on the circles sheet (or the "
    "non-circular surface) instead, exactly as a Single surface run solves "
    "it. Untick when the surface is prescribed — a published benchmark's "
    "circle, an observed failure surface — and the question is that "
    "surface's own reliability.")
SEARCH_HINT_SAMPLING = ("One surface serves every realization — the searched one, "
                        "or, unticked, the first circle on the circles sheet "
                        "(or the non-circular surface), as a Single surface run "
                        "solves it.")


class ReliabilityDialog(QDialog):
    """Probabilistic reliability analysis — the sibling of the Parametric study.

    Where the Parametric study answers deterministic what-ifs, this dialog turns
    the material standard deviations (the ``s(·)`` columns of the mat sheet) into a
    reliability index and a probability of failure. It offers three engines:

    - **Taylor series (TSPM)** — ``1 + 2N`` limit-equilibrium searches (or FEM SSRM
      solves in FEM mode): the mean-value factor of safety plus a ±σ perturbation
      per uncertain parameter.
    - **Monte Carlo** — ``N`` sampled realizations of every uncertain parameter,
      each evaluated on a fixed failure surface, reported as an FS histogram. Monte
      Carlo needs ~10⁴ factor-of-safety evaluations, which is affordable with a
      limit-equilibrium solve but not with the finite-element SSRM, so it is
      disabled in FEM mode (FEM reliability stays on the Taylor series).
    - **Response surface (RS)** — the same sampling, with the factor of safety taken
      from a quadratic surrogate fitted to a few dozen real solves and measured
      against held-out ones before it is used. Ten million realizations at no
      sampling cost, so the tail of the distribution is resolved. Limit-equilibrium
      only, for the same reason as Monte Carlo.

    A read-only summary lists the file's σ columns so the user can confirm what the
    analysis will vary before running.
    """

    def __init__(self, parent=None, defaults=None, slope_data=None, app_mode="lem",
                 document=None):
        super().__init__(parent)
        self.app_mode = "fem" if app_mode == "fem" else "lem"
        self.setWindowTitle("Reliability analysis"
                            + (" (FEM · SSRM)" if self.app_mode == "fem" else " (LEM)"))
        defaults = defaults or {}
        slope_data = slope_data if slope_data is not None else {}
        self._sd = slope_data

        from xslope.sensitivity import list_params
        params = list_params(slope_data, mode="lem")
        self._sigma_params = [e for e in params if e.get("sigma")]
        self._has_sigma = bool(self._sigma_params)

        layout = QVBoxLayout()                  # the left column; see two_pane
        form = QFormLayout()

        # --- engine (method) selector --------------------------------------
        self.engine = QComboBox()
        for key, label in RELIABILITY_ENGINES:
            self.engine.addItem(label, key)
        # The sampling engines are limit-equilibrium only; grey them out in FEM mode.
        if self.app_mode == "fem":
            for key in RELIABILITY_SAMPLING_ENGINES:
                item = self.engine.model().item(self.engine.findData(key))
                item.setEnabled(False)
                item.setToolTip("Sampling the standard deviations needs thousands of "
                                "solves — limit-equilibrium only.")
        want = defaults.get("engine", "taylor")
        if self.app_mode == "fem":
            want = "taylor"
        eidx = self.engine.findData(want)
        if eidx >= 0:
            self.engine.setCurrentIndex(eidx)
        form.addRow("Method", self.engine)

        self._mc_disabled_note = QLabel(
            "Monte Carlo and the response surface are disabled for FEM — ~10⁴ SSRM "
            "solves is prohibitive; FEM reliability uses the Taylor series.")
        self._mc_disabled_note.setWordWrap(True)
        self._mc_disabled_note.setVisible(self.app_mode == "fem")

        # --- LEM solver + surface (LEM mode only) --------------------------
        self.method = self._combo(LEM_METHODS, defaults.get("method", "spencer"))
        if self.app_mode == "lem":
            form.addRow("LEM method", self.method)

        # Presence, not the stored selection (see RunLemDialog): reading the selection
        # as presence would hide the circles from a file whose last run chose the
        # non-circular surface.
        has_circ = bool(slope_data.get("circles"))
        has_ncirc = bool(slope_data.get("non_circ"))
        if self.app_mode == "lem":
            if has_circ != has_ncirc:                  # exactly one defined
                self._fixed_surface = "circular" if has_circ else "noncircular"
                self.surface = None
                form.addRow("Surface", QLabel(dict(SURFACE_TYPES)[self._fixed_surface]))
            else:
                self._fixed_surface = None
                stored = (stored_surface_family(slope_data)
                          if (has_circ and has_ncirc) else "circular")
                self.surface = self._combo(SURFACE_TYPES,
                                           defaults.get("surface", stored))
                form.addRow("Surface", self.surface)
        else:
            self._fixed_surface = None
            self.surface = None

        self.num_slices = QSpinBox()
        self.num_slices.setRange(5, 500)
        self.num_slices.setValue(int(defaults.get("num_slices", 40)))
        if self.app_mode == "lem":
            form.addRow("Number of slices", self.num_slices)

        self.rapid = QCheckBox("Rapid drawdown")
        self.rapid.setChecked(bool(defaults.get("rapid", False)))
        if self.app_mode == "lem":
            form.addRow("", self.rapid)

        # The checkbox's label is DYNAMIC per engine, because the truthful
        # statement differs: the Taylor series searches at every one of its
        # 1+2N solves; the sampling engines search once at the mean values and
        # hold that surface across their realizations. One static label lies
        # for one engine or the other (owner + pushback, 2026-08-15).
        self.search = QCheckBox(SEARCH_LABEL_TAYLOR)
        self.search.setChecked(bool(defaults.get("search", True)))
        self.search.setToolTip(SEARCH_TIP_TAYLOR)
        self.search_hint = QLabel(SEARCH_HINT_TAYLOR)
        self.search_hint.setWordWrap(True)
        self.search_hint.setStyleSheet("color: gray; font-size: 11px;")
        if self.app_mode == "lem":
            form.addRow("", self.search)
            form.addRow("", self.search_hint)

        # --- Monte Carlo controls (LEM only) -------------------------------
        self.n_samples = QSpinBox()
        self.n_samples.setRange(100, 500000)
        self.n_samples.setSingleStep(1000)
        self.n_samples.setValue(int(defaults.get("n_samples", 10000)))
        self.n_samples.setToolTip(
            "Number of sampled realizations (10000 by default). Monte Carlo only: "
            "the response surface samples ten million realizations of its "
            "surrogate, a count that costs seconds of arithmetic rather than "
            "solves, so there is no trade-off to set here.")
        self.seed = QSpinBox()
        self.seed.setRange(0, 2_000_000_000)
        self.seed.setValue(int(defaults.get("rng_seed", 20240117)))
        self.seed.setToolTip("Random seed — a fixed value makes the run reproducible. "
                             "Applies to Monte Carlo and to the response surface.")
        self.distribution = self._combo(MC_DISTRIBUTIONS,
                                        defaults.get("distribution", "normal"))
        self.mc_sampling = self._combo(
            [("lhs", "Latin hypercube"), ("random", "Random")],
            defaults.get("sampling", "lhs"))
        self.mc_sampling.setToolTip(
            "How realizations are drawn. Latin hypercube (the default) "
            "stratifies each parameter into equal-probability bins (one draw "
            "per bin), which reduces sampling scatter for the same count; "
            "draws stay reproducible under the fixed seed. Random draws every "
            "value independently. The convergence stop's confidence band "
            "assumes independent draws, so under Latin hypercube it errs "
            "conservative.")
        # Statistical-convergence stopping (off by default): check the empirical
        # P_f every 100 realizations and stop when its 95% confidence half-width
        # falls inside the tolerance — the samples field becomes the cap. The
        # rule never fires before 500 realizations or 10 observed failures, so a
        # rare-event campaign runs to the cap rather than stopping on false
        # confidence.
        self.mc_converge = QCheckBox("Stop when P_f converges")
        self.mc_converge.setChecked(defaults.get("converge_rel") is not None)
        self.mc_converge.setToolTip(
            "Stop sampling once the 95% confidence half-width on the empirical "
            "probability of failure is inside the tolerance — further "
            "realizations no longer move the answer at that resolution. "
            "MC samples becomes the cap. With a fixed seed the stopped run is "
            "still exactly reproducible.")
        self.mc_converge_tol = QDoubleSpinBox()
        self.mc_converge_tol.setDecimals(0)
        self.mc_converge_tol.setRange(1.0, 50.0)
        self.mc_converge_tol.setSingleStep(5.0)
        self.mc_converge_tol.setSuffix(" % of P_f")
        self.mc_converge_tol.setValue(
            float(defaults.get("converge_rel") or 0.05) * 100)
        self.mc_converge_tol.setToolTip(
            "Stop when the probability of failure is known to this fraction of "
            "itself (95% confidence half-width ≤ tolerance × P_f). Relative, so "
            "it self-scales: a small P_f automatically demands more "
            "realizations than a large one, roughly (1−p)/p of them.")
        if self.app_mode == "lem":
            form.addRow("MC samples", self.n_samples)
            form.addRow("MC seed", self.seed)
            form.addRow("MC distribution", self.distribution)
            form.addRow("MC sampling", self.mc_sampling)
            form.addRow("", self.mc_converge)
            form.addRow("P_f tolerance (±)", self.mc_converge_tol)

        # --- FEM SSRM bracket + tolerance (FEM only) -----------------------
        self.f_min = QDoubleSpinBox()
        self.f_min.setRange(0.1, 10.0)
        self.f_min.setSingleStep(0.05)
        self.f_min.setValue(float(defaults.get("F_min", 0.7)))
        self.f_max = QDoubleSpinBox()
        self.f_max.setRange(0.1, 20.0)
        self.f_max.setSingleStep(0.05)
        self.f_max.setValue(float(defaults.get("F_max", 2.0)))
        self.reliability_tol = QDoubleSpinBox()
        self.reliability_tol.setDecimals(4)
        self.reliability_tol.setRange(0.0001, 0.1)
        self.reliability_tol.setSingleStep(0.0005)
        self.reliability_tol.setValue(float(defaults.get("reliability_tol", 0.001)))
        self.reliability_tol.setToolTip(
            "Bisection tolerance for the reliability SSRM solves — tighter than a "
            "single run (default 0.001). TSPM amplifies factor-of-safety imprecision, "
            "so a tight tolerance keeps the reliability index stable.")
        if self.app_mode == "fem":
            form.addRow("F min (SSRM)", self.f_min)
            form.addRow("F max (SSRM)", self.f_max)
            form.addRow("Tolerance (SSRM)", self.reliability_tol)

        layout.addLayout(form)
        if self.app_mode == "fem":
            layout.addWidget(self._mc_disabled_note)

        # --- read-only σ summary -------------------------------------------
        sigma_box = QGroupBox("Standard deviations in this file")
        sform = QVBoxLayout(sigma_box)
        if self._has_sigma:
            for e in self._sigma_params:
                val, sig = e.get("value"), e.get("sigma")
                cov = (f"  (COV {sig / val * 100:.0f}%)"
                       if val not in (None, 0) else "")
                lbl = QLabel(f"• {e['ref']} = {('' if val is None else f'{val:g}')} "
                             f"± {sig:g}{cov}")
                sform.addWidget(lbl)
        else:
            warn = QLabel("No standard deviations (s-columns) are set in the mat "
                          "sheet — reliability needs at least one.")
            warn.setWordWrap(True)
            sform.addWidget(warn)
        layout.addWidget(sigma_box)

        note = QLabel(
            "Reliability turns the σ columns into a reliability index β and a "
            "probability of failure. The Taylor series runs 1+2N solves; Monte Carlo "
            "samples the σ's and reports an FS histogram; the response surface fits "
            "a surrogate to a few dozen solves, checks it against held-out ones, and "
            "samples it ten million times.")
        note.setWordWrap(True)
        layout.addWidget(note)

        # A reliability run inherits every rule of the analysis it repeats, so the
        # base model is checked here — thousands of solves is a long way to travel
        # before learning that a conductivity or a Poisson's ratio was blank.
        self.preflight = PreflightPanel(
            analysis="reliability", slope_data=slope_data, document=document,
            selection_fn=lambda: {"base": self.app_mode,
                                  "method": self.method.currentData(),
                                  "surface": self._surface_value(),
                                  "search": self.search.isChecked()},
            notes=(SEISMIC_NOTE_LEM if self.app_mode == "lem" else SEISMIC_NOTE_FEM,),
            parent=self)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self._ok = bb.button(QDialogButtonBox.Ok)
        self._ok.setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)

        two_pane(self, layout, self.preflight, bb, fixed=True)

        self.engine.currentIndexChanged.connect(self._sync_enabled)
        self.mc_converge.toggled.connect(self._sync_enabled)
        self.method.currentIndexChanged.connect(self.preflight.refresh)
        if self.surface is not None:
            self.surface.currentIndexChanged.connect(self.preflight.refresh)
        self.preflight.changed.connect(self._sync_run)
        self._sync_enabled()
        self._sync_run()

    def _sync_run(self):
        """Run needs a sigma AND no error; a warning never blocks."""
        if self.app_mode == "lem":
            from xslope.preflight import capabilities
            apply_capabilities(self.method,
                               capabilities(self.preflight.model,
                                            {"method": self.method.currentData(),
                                             "surface": self._surface_value(),
                                             "search": self.search.isChecked()}
                                            ).get("lem_method", {}))
        blocked = self.preflight.blocked
        self._ok.setEnabled(self._has_sigma and not blocked)
        self._ok.setToolTip(
            self.preflight.block_reason() if blocked else
            "" if self._has_sigma else
            "No standard deviations (s-columns) are set in the mat sheet.")

    @staticmethod
    def _combo(items, selected):
        combo = QComboBox()
        for key, label in items:
            combo.addItem(label, key)
        idx = combo.findData(selected)
        if idx >= 0:
            combo.setCurrentIndex(idx)
        return combo

    def _engine_value(self):
        return self.engine.currentData()

    def _surface_value(self):
        return self._fixed_surface or (self.surface.currentData()
                                       if self.surface is not None else "circular")

    def _sync_enabled(self):
        lem = self.app_mode == "lem"
        engine = self._engine_value()
        mc = engine == "mc" and lem
        # The seed and the distribution describe the DRAWS, so they apply to both
        # sampling engines. The sample count and the convergence stop describe how
        # many realizations are solved, which is a Monte Carlo question only: the
        # response surface samples a fixed ten million realizations of its
        # surrogate and its accuracy is set by the fit, not by the count.
        sampling = engine in RELIABILITY_SAMPLING_ENGINES and lem
        for w in (self.seed, self.distribution, self.mc_sampling):
            w.setEnabled(sampling)
        self.search.setText(SEARCH_LABEL_SAMPLING if sampling
                            else SEARCH_LABEL_TAYLOR)
        self.search_hint.setText(SEARCH_HINT_SAMPLING if sampling
                                 else SEARCH_HINT_TAYLOR)
        self.search.setToolTip(SEARCH_TIP_SAMPLING if sampling
                               else SEARCH_TIP_TAYLOR)
        for w in (self.n_samples, self.mc_converge):
            w.setEnabled(mc)
        self.mc_converge_tol.setEnabled(mc and self.mc_converge.isChecked())

    def options(self):
        return {
            "engine": self._engine_value(),
            "app_mode": self.app_mode,
            "method": self.method.currentData(),
            "surface": self._surface_value(),
            "num_slices": self.num_slices.value(),
            "rapid": self.rapid.isChecked(),
            "search": self.search.isChecked(),
            "n_samples": self.n_samples.value(),
            "rng_seed": self.seed.value(),
            "distribution": self.distribution.currentData(),
            "sampling": self.mc_sampling.currentData(),
            "converge_rel": (self.mc_converge_tol.value() / 100
                             if self.mc_converge.isChecked() else None),
            "F_min": self.f_min.value(),
            "F_max": self.f_max.value(),
            "reliability_tol": self.reliability_tol.value(),
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
            "• SLOPE/W's search definition and reinforcement sets are not imported — "
            "you will get a list of what was left out."))

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


class UnpackPackageDialog(QDialog):
    """Opening a project package (.xslz): where its files go.

    A package is transport, not a place to work — so opening one extracts it to a
    folder of loose files and the project is opened from there. The default folder
    is named for the package and sits beside it, which is where the user will look
    for the workbook when they want it in Excel.

    If that folder already exists it may hold the user's own edits, so the dialog
    never quietly reuses or overwrites it: the buttons become **Open Existing**
    (leave the folder alone and open the project already in it) and **Extract
    Fresh** (unpack into a new numbered folder beside it).

    :meth:`chosen` returns ``(destination, mode)`` with mode ``"extract"`` or
    ``"existing"``. (Not ``result()`` — QDialog already has one, and its meaning is
    the accepted/rejected code every Qt caller expects.)
    """

    def __init__(self, package, parent=None):
        super().__init__(parent)
        from xslope.package import unpack_path

        self._package = str(package)
        self.setWindowTitle("Open project package")
        layout = QVBoxLayout(self)
        note = QLabel(
            f"{os.path.basename(self._package)} is a project package: a workbook and "
            f"its results zipped together to travel as one file. It will be unpacked "
            f"to the folder below and the project opened from there, so the workbook "
            f"is an ordinary .xlsx you can also open in Excel.")
        note.setWordWrap(True)
        note.setMaximumWidth(460)
        layout.addWidget(note)

        row = QHBoxLayout()
        self.dest = QLineEdit(unpack_path(self._package))
        # Wide enough for a real path in whatever font and scaling the platform is
        # using, measured rather than guessed in pixels.
        self.dest.setMinimumWidth(self.dest.fontMetrics().horizontalAdvance("n" * 48))
        change = QPushButton("Change…")
        change.clicked.connect(self._change)
        row.addWidget(self.dest, 1)
        row.addWidget(change)
        row.setContentsMargins(0, 0, 0, 0)
        holder = QWidget()
        holder.setLayout(row)
        layout.addWidget(holder)

        self.status = QLabel("")
        self.status.setWordWrap(True)
        self.status.setMaximumWidth(460)
        layout.addWidget(self.status)

        self.buttons = QDialogButtonBox(QDialogButtonBox.Cancel)
        self.btn_unpack = self.buttons.addButton("Unpack and Open",
                                                 QDialogButtonBox.AcceptRole)
        self.btn_existing = self.buttons.addButton("Open Existing",
                                                   QDialogButtonBox.AcceptRole)
        self.btn_fresh = self.buttons.addButton("Extract Fresh",
                                                QDialogButtonBox.AcceptRole)
        self.buttons.clicked.connect(self._clicked)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

        self._mode = "extract"
        self.dest.textChanged.connect(self._refresh)
        self._refresh()

    def _change(self):
        start = os.path.dirname(self.dest.text()) or os.path.dirname(self._package)
        folder = QFileDialog.getExistingDirectory(
            self, "Choose a folder to unpack into", start)
        if folder:
            # The chosen folder is the PARENT: the project still lands in a folder of
            # its own, named for the package, so its loose files never mix with
            # whatever else is in there.
            self.dest.setText(os.path.join(
                folder, os.path.splitext(os.path.basename(self._package))[0]))

    def _refresh(self):
        """Show the buttons the chosen destination allows."""
        exists = os.path.exists(self.dest.text())
        self.btn_unpack.setVisible(not exists)
        self.btn_existing.setVisible(exists)
        self.btn_fresh.setVisible(exists)
        if exists:
            self.status.setText(
                f"That folder already exists, and the project in it may hold edits of "
                f"your own. Open Existing leaves it untouched; Extract Fresh unpacks "
                f"into {os.path.basename(self._fresh_dest())} beside it.")
        else:
            self.status.setText("")

    def _fresh_dest(self):
        """The first free numbered folder beside the chosen destination."""
        base = self.dest.text().rstrip(os.sep)
        n = 2
        while os.path.exists(f"{base}-{n}"):
            n += 1
        return f"{base}-{n}"

    def _clicked(self, button):
        if button is self.btn_existing:
            self._mode = "existing"
            self.accept()
        elif button is self.btn_fresh:
            self._mode = "extract"
            self.dest.setText(self._fresh_dest())
            self.accept()
        elif button is self.btn_unpack:
            self._mode = "extract"
            self.accept()

    def chosen(self):
        """``(destination folder, "extract" | "existing")``."""
        return self.dest.text(), self._mode
