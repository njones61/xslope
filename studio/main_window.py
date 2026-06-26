"""MainWindow — the XSlope Studio shell (Phase 1: read-only viewer).

Provides the app frame: menus, a dockable Inputs summary panel, a Log pane that
captures engine stdout/stderr, a mode selector (LEM / Seep / FEM), recent-files,
and the central embedded Matplotlib canvas. File -> Open loads an Excel input
file via ProjectDocument and renders the Inputs view. Editing, running analyses,
and saving arrive in later phases.
"""

from __future__ import annotations

import os
import sys
import traceback
from pathlib import Path

from PySide6.QtCore import Qt, QObject, QSettings, QThread, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QComboBox, QDockWidget, QFileDialog, QHBoxLayout, QLabel, QMainWindow,
    QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QStackedWidget,
    QTabWidget, QToolBar, QToolButton, QTreeWidget, QTreeWidgetItem, QWidget,
)

from .canvas import MplCanvas
from .dialogs import BuildMeshDialog, RunFemDialog, RunLemDialog, RunSeepDialog
from .display_panels import (
    FeDataDisplayPanel, FemResultsDisplayPanel, InputsDisplayPanel,
    MeshDisplayPanel, SearchDisplayPanel, SeepDisplayPanel, SolutionDisplayPanel,
)
from .document import ProjectDocument
from .editors import CATEGORY_EDITORS
from .runners import FemRunner, LemRunner, MeshWorker, SeepRunner

APP_NAME = "XSlope Studio"
ORG_NAME = "XSlope"
MAX_RECENT = 8
MODES = [("LEM", "lem"), ("Seepage", "seep"), ("FEM", "fem")]
# Blank template bundled with the app; used to create files on Save As.
TEMPLATE = Path(__file__).resolve().parent / "resources" / "input_template.xlsx"
CATEGORY_ROLE = Qt.UserRole + 1


class _LogStream(QObject):
    """stdout/stderr tee: forwards to the original stream and the log pane.

    Writes are marshaled to the GUI thread via a queued signal, so engine output
    printed from a worker thread (e.g. an LEM solve) streams into the log pane
    live and safely — Qt widgets must not be touched off the GUI thread.
    """

    _emitted = Signal(str)

    def __init__(self, widget, original):
        super().__init__()
        self._original = original
        # Queued so a worker-thread write() appends on the widget's (GUI) thread.
        self._emitted.connect(widget.appendPlainText, Qt.QueuedConnection)

    def write(self, text):
        if self._original is not None:
            try:
                self._original.write(text)
            except Exception:
                pass
        text = text.rstrip("\n")
        if text:
            self._emitted.emit(text)

    def flush(self):
        if self._original is not None:
            try:
                self._original.flush()
            except Exception:
                pass


class MainWindow(QMainWindow):
    # Emitted to hand a mesh build to the persistent mesh thread (queued).
    _mesh_requested = Signal(object, object)

    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(1200, 800)
        self.settings = QSettings(ORG_NAME, APP_NAME)

        self.doc = ProjectDocument(self)
        self.doc.loaded.connect(self._on_loaded)
        self.doc.changed.connect(self._render)
        self.doc.dirty_changed.connect(lambda *_: self._update_title())

        # Central area: a tab strip of result views (plan §7). The Inputs view is
        # always present; the LEM Solution view is added after the first solve.
        self.canvas = MplCanvas(self)
        self.mesh_canvas = None
        self.search_canvas = None
        self.solution_canvas = None
        self.reliability_canvas = None
        self.seep_data_canvas = {}        # bc set -> MplCanvas
        self.seep_solution_canvas = {}    # bc set -> MplCanvas
        self.fem_data_canvas = None
        self.fem_results_canvas = None
        self.view_tabs = QTabWidget()
        self.view_tabs.addTab(self.canvas, "Inputs")
        self.view_tabs.currentChanged.connect(self._on_view_tab_changed)
        self.setCentralWidget(self.view_tabs)

        self._mode = "lem"
        self._runner = None
        self._seep_runner = None
        self._fem_runner = None
        self._mesh_busy = False
        self._run_implemented = {"lem", "seep", "fem"}   # modes whose Run is wired up
        self._last_lem_opts = {}
        self._last_mesh_opts = {}
        self._last_seep_opts = {}
        self._last_fem_opts = {}

        # gmsh must run on one consistent thread (it segfaults if re-initialized on
        # a fresh thread each build), so a single persistent worker thread handles
        # every mesh build via a queued request signal.
        self._mesh_thread = QThread(self)
        self._mesh_worker = MeshWorker()
        self._mesh_worker.moveToThread(self._mesh_thread)
        self._mesh_requested.connect(self._mesh_worker.build)
        self._mesh_worker.succeeded.connect(self._on_mesh_succeeded)
        self._mesh_worker.failed.connect(self._on_mesh_failed)
        self._mesh_thread.start()
        self._recent = [p for p in (self.settings.value("recent_files") or []) if p]

        self._display_panels = {}     # result tab widget -> its display-options panel
        self._make_inputs_dock()
        self._make_display_dock()
        self._make_log_dock()
        # Let the left dock column own the bottom-left corner so it runs the full
        # window height; the Log dock then starts at the central canvas's left edge
        # instead of spanning under the Inputs/Display docks.
        self.setCorner(Qt.BottomLeftCorner, Qt.LeftDockWidgetArea)
        self._arrange_docks()
        self._make_actions()
        self._make_menus()
        self._make_toolbar()
        self._update_recent_menu()
        self._update_title()
        self._install_log_capture()
        # A run progress bar + Cancel button live at the right of the status bar
        # (hidden when idle).
        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximumWidth(220)
        self.progress_bar.setVisible(False)
        self.statusBar().addPermanentWidget(self.progress_bar)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setVisible(False)
        self.cancel_btn.clicked.connect(self._cancel_run)
        self.statusBar().addPermanentWidget(self.cancel_btn)
        self._update_run_actions()    # initial labels / visibility (no file yet)
        self.statusBar().showMessage("Open an Excel input file to begin.")

    # --- docks -----------------------------------------------------------
    def _make_inputs_dock(self):
        self.inputs_tree = QTreeWidget()
        self.inputs_tree.setHeaderLabels(["Input", "Count"])
        self.inputs_tree.setColumnWidth(0, 180)
        self.inputs_tree.itemClicked.connect(self._on_tree_click)
        dock = QDockWidget("Inputs", self)
        dock.setObjectName("inputs_dock")
        dock.setWidget(self.inputs_tree)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.inputs_dock = dock

    def _make_display_dock(self):
        # Context-sensitive display options: a stack whose page follows the active
        # result tab. Sits under the Inputs tree.
        self.display_stack = QStackedWidget()
        self._display_placeholder = QLabel("No display options for this view.")
        self._display_placeholder.setWordWrap(True)
        self._display_placeholder.setContentsMargins(8, 8, 8, 8)
        self._display_placeholder.setAlignment(Qt.AlignTop)
        self.display_stack.addWidget(self._display_placeholder)
        dock = QDockWidget("Display", self)
        dock.setObjectName("display_dock")
        dock.setWidget(self.display_stack)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.splitDockWidget(self.inputs_dock, dock, Qt.Vertical)
        self.display_dock = dock

        # The Inputs view's display options (always present — the Inputs canvas
        # exists for the life of the window and survives _clear_result_tabs).
        self.inputs_panel = InputsDisplayPanel()
        self.inputs_panel.changed.connect(self._render)
        self.display_stack.addWidget(self.inputs_panel)
        self._display_panels[self.canvas] = self.inputs_panel

    def _make_log_dock(self):
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(5000)
        dock = QDockWidget("Log", self)
        dock.setObjectName("log_dock")
        dock.setWidget(self.log)
        # Custom title bar: the "Log" label plus a right-aligned Clear button that
        # empties the pane (like clearing a terminal). The custom widget is still
        # the dock's drag handle; visibility is toggled from the View menu.
        title = QWidget()
        row = QHBoxLayout(title)
        row.setContentsMargins(6, 2, 4, 2)
        row.addWidget(QLabel("Log"))
        row.addStretch(1)
        clear_btn = QToolButton()
        clear_btn.setText("Clear")
        clear_btn.setAutoRaise(True)
        clear_btn.setToolTip("Clear the log output")
        clear_btn.clicked.connect(self.log.clear)
        row.addWidget(clear_btn)
        dock.setTitleBarWidget(title)
        self.addDockWidget(Qt.BottomDockWidgetArea, dock)
        self.log_dock = dock

    def _arrange_docks(self):
        # A QTreeWidget gives the dock almost no width hint, so without an explicit
        # width the left column collapses to ~90px. Size it to fit the Inputs tree
        # (180px name col + Count) without being wide.
        self.resizeDocks([self.inputs_dock, self.display_dock], [290, 290],
                         Qt.Horizontal)
        # The left column spans the full height (it owns the bottom-left corner).
        # Give the Inputs tree enough height to show all categories without
        # scrolling; the Display panel takes the larger remaining share.
        self.resizeDocks([self.inputs_dock, self.display_dock], [300, 430],
                         Qt.Vertical)

    def _install_log_capture(self):
        sys.stdout = _LogStream(self.log, sys.__stdout__)
        sys.stderr = _LogStream(self.log, sys.__stderr__)

    # --- actions / menus -------------------------------------------------
    def _make_actions(self):
        self.act_new = QAction("&New", self, shortcut=QKeySequence.New,
                               triggered=self.new_project)
        self.act_open = QAction("&Open…", self, shortcut=QKeySequence.Open,
                                triggered=self.open_dialog)
        self.act_quit = QAction("&Quit", self, shortcut=QKeySequence.Quit,
                                triggered=self.close)
        self.act_undo = QAction("&Undo", self, shortcut=QKeySequence.Undo,
                                triggered=self.doc.undo, enabled=False)
        self.act_redo = QAction("&Redo", self, shortcut=QKeySequence.Redo,
                                triggered=self.doc.redo, enabled=False)
        self.act_about = QAction("&About", self, triggered=self._about)
        self.act_save = QAction("&Save", self, shortcut=QKeySequence.Save,
                                enabled=False, triggered=self.save)
        self.act_save_as = QAction("Save &As…", self, enabled=False, triggered=self.save_as)
        self.act_run = QAction("Run &LEM…", self, enabled=False, triggered=self.run_current)
        self.act_build_mesh = QAction("Build &Mesh…", self, enabled=False,
                                      triggered=self.build_mesh)

    def _make_menus(self):
        mb = self.menuBar()

        m_file = mb.addMenu("&File")
        m_file.addAction(self.act_new)
        m_file.addAction(self.act_open)
        self.recent_menu = m_file.addMenu("Open &Recent")
        m_file.addSeparator()
        m_file.addAction(self.act_save)
        m_file.addAction(self.act_save_as)
        m_file.addSeparator()
        m_file.addAction(self.act_quit)

        m_edit = mb.addMenu("&Edit")
        m_edit.addAction(self.act_undo)
        m_edit.addAction(self.act_redo)

        m_run = mb.addMenu("&Run")
        m_run.addAction(self.act_build_mesh)
        m_run.addSeparator()
        m_run.addAction(self.act_run)

        m_view = mb.addMenu("&View")
        m_view.addAction(self.inputs_dock.toggleViewAction())
        m_view.addAction(self.display_dock.toggleViewAction())
        m_view.addAction(self.log_dock.toggleViewAction())

        m_help = mb.addMenu("&Help")
        m_help.addAction(self.act_about)

    def _make_toolbar(self):
        tb = QToolBar("Main", self)
        tb.setObjectName("main_toolbar")
        self.addToolBar(tb)
        tb.addAction(self.act_new)
        tb.addAction(self.act_open)
        tb.addSeparator()
        self.mode_label = QLabel(" Mode: ")
        tb.addWidget(self.mode_label)
        self.mode_combo = QComboBox()
        for label, _ in MODES:
            self.mode_combo.addItem(label)
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        tb.addWidget(self.mode_combo)
        tb.addSeparator()
        tb.addAction(self.act_build_mesh)
        tb.addAction(self.act_run)
        # macOS's native style draws text-only toolbar buttons in the larger system
        # font and ignores setFont; a stylesheet forces the size so New/Open/Run LEM
        # match the "Mode:" label. pointSizeF() is -1 for pixel-defined fonts.
        pt = self.mode_label.font().pointSizeF()
        if pt > 0:
            tb.setStyleSheet(f"QToolButton {{ font-size: {pt:g}pt; }}")

    # --- new / open / recent ---------------------------------------------
    def new_project(self):
        if not self._confirm_discard():
            return
        self.doc.new()
        self.statusBar().showMessage(
            "New project — edit the inputs, then Save As to write an Excel file.")

    def _confirm_discard(self):
        """Prompt to save/discard if the current doc has unsaved edits.
        Returns False if the user cancels (caller should abort)."""
        if not (self.doc.is_open and self.doc.dirty):
            return True
        res = QMessageBox.question(
            self, "Unsaved changes", "Discard unsaved changes to the current project?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
        if res == QMessageBox.Cancel:
            return False
        if res == QMessageBox.Save:
            self.save()
            return not self.doc.dirty   # abort if the save failed/was cancelled
        return True

    def open_dialog(self):
        start = os.path.dirname(self._recent[0]) if self._recent else ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Open XSlope input file", start, "Excel files (*.xlsx);;All files (*)")
        if path:
            self.open_path(path)

    def open_path(self, path):
        if not self._confirm_discard():
            return
        try:
            self.doc.load(path)
        except Exception as exc:  # ValueError from the loader, or anything else
            traceback.print_exc()
            QMessageBox.critical(self, "Could not open file",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return
        self._add_recent(path)

    def _add_recent(self, path):
        path = os.path.abspath(path)
        self._recent = [path] + [p for p in self._recent if p != path]
        self._recent = self._recent[:MAX_RECENT]
        self.settings.setValue("recent_files", self._recent)
        self._update_recent_menu()

    def _update_recent_menu(self):
        self.recent_menu.clear()
        if not self._recent:
            empty = self.recent_menu.addAction("(none)")
            empty.setEnabled(False)
            return
        for path in self._recent:
            act = self.recent_menu.addAction(path)
            act.triggered.connect(lambda _=False, p=path: self.open_path(p))

    # --- document signal handlers ---------------------------------------
    def _on_loaded(self):
        self.act_save.setEnabled(True)
        self.act_save_as.setEnabled(True)
        self._clear_result_tabs()
        # Restore saved solutions first so the default mode can see whether an FEM
        # solution exists, then pick the mode that fits this file.
        self._load_solution_sidecars()
        self._mode = self._default_mode()
        self.mode_combo.blockSignals(True)    # set silently; we render explicitly below
        self.mode_combo.setCurrentIndex([m for _, m in MODES].index(self._mode))
        self.mode_combo.blockSignals(False)
        self._update_run_actions()
        self.canvas.reset_fit()       # fit the fresh file to the window
        self._render()
        self._populate_inputs_tree()
        self._update_title()
        n = len(self.doc.slope_data.get("materials", []))
        name = os.path.basename(self.doc.path) if self.doc.path else "untitled"
        self.statusBar().showMessage(
            f"Loaded {name} — {n} material(s). "
            f"Click an underlined input to edit it.")

    def _default_mode(self):
        """Pick the mode that fits the loaded file: FEM if a mesh and a restored
        FEM solution are present; Seep if the materials carry only seepage
        properties (conductivity, no strength); otherwise LEM."""
        sd = self.doc.slope_data
        if sd.get("mesh") is not None and "fem_solution" in self.doc.results:
            return "fem"
        if self._materials_seep_only(sd.get("materials", [])):
            return "seep"
        return "lem"

    @staticmethod
    def _materials_seep_only(materials):
        """True when at least one material defines seepage conductivity and none
        defines strength — i.e. a pure seepage problem. A material is usable for
        LEM only if it has unit weight / cohesion / friction angle; when gamma, c
        and phi are all blank for every material, the file cannot be analyzed by
        LEM."""
        def num(v):
            try:
                return float(v or 0)
            except (TypeError, ValueError):
                return 0.0

        def lem_capable(m):
            return num(m.get("gamma")) > 0 or num(m.get("c")) > 0 or num(m.get("phi")) > 0

        def has_seep(m):
            return num(m.get("k1")) > 0 or num(m.get("k2")) > 0

        if not materials:
            return False
        return any(has_seep(m) for m in materials) and not any(
            lem_capable(m) for m in materials)

    def _load_solution_sidecars(self):
        """Restore any saved seep / FEM solution sidecars next to the .xlsx so
        their result tabs appear immediately — no re-solve needed. The mesh is
        loaded by ``load_slope_data`` (``{stem}_mesh.json``); we rebuild the
        seep/FEM data on it and read the saved nodal/element results back in.
        Best-effort: a mismatched or unreadable sidecar is skipped, not fatal."""
        if not self.doc.path:
            return
        mesh = self.doc.slope_data.get("mesh")
        if mesh is None:                       # no mesh -> no FE solution to restore
            return
        stem = os.path.splitext(self.doc.path)[0]
        self._restore_seep_sidecar(mesh, stem)
        self._restore_fem_sidecar(mesh, stem)

    def _restore_seep_sidecar(self, mesh, stem):
        # Restore each BC set that has a sidecar: _seep.csv (BC 1) and _seep2.csv
        # (BC 2, rapid drawdown). Each lands in its own pair of tabs.
        for bc, suffix in ((1, "_seep.csv"), (2, "_seep2.csv")):
            path = f"{stem}{suffix}"
            if not os.path.exists(path):
                continue
            try:
                from xslope.seep import build_seep_data, import_seep_solution
                seep_data = build_seep_data(mesh, self.doc.slope_data, seep_bc=bc)
                solution = import_seep_solution(seep_data, path)
            except Exception:
                traceback.print_exc()   # streams to the Log pane; load still succeeds
                continue
            self.doc.results.setdefault("seep_solutions", {})[bc] = {
                "seep_data": seep_data, "solution": solution, "options": {"bc": bc}}
            self._show_seep_data(seep_data, bc)
            self._show_seep_solution(bc)
            print(f"Restored saved seepage solution (BC set {bc}) from "
                  f"{os.path.basename(path)}.")

    def _restore_fem_sidecar(self, mesh, stem):
        if not os.path.exists(f"{stem}_fem_nodes.csv"):
            return
        try:
            from xslope.fem import build_fem_data, import_fem_solution, import_fem_meta
            fem_data = build_fem_data(self.doc.slope_data, mesh)
            solution = import_fem_solution(fem_data, stem)
            meta = import_fem_meta(stem) or {}
        except Exception:
            traceback.print_exc()
            return
        self.doc.results["fem_solution"] = {
            "fem_data": fem_data, "solution": solution, "FS": meta.get("FS"),
            "analysis": meta.get("analysis") or "loaded"}
        self._show_fem_data(fem_data)
        self._show_fem_results()
        fs = meta.get("FS")
        fs_note = f" (SSRM FS = {fs:.3f})" if isinstance(fs, (int, float)) else ""
        print(f"Restored saved FEM solution from {os.path.basename(stem)}_fem_*.csv{fs_note}.")

    def _render(self):
        if not self.doc.is_open:
            self.canvas.clear()
            return
        try:
            self.canvas.render_inputs(
                self.doc.slope_data, mode=self._mode,
                opts=self.inputs_panel.options())
        except Exception:
            traceback.print_exc()
        self.act_undo.setEnabled(self.doc.can_undo())
        self.act_redo.setEnabled(self.doc.can_redo())

    def _on_mode_changed(self, index):
        self._mode = MODES[index][1]
        self._update_run_actions()
        if self.doc.is_open:
            self._render()
            self._populate_inputs_tree()

    def _update_run_actions(self):
        """Keep the single Run action's label and the Build Mesh action in sync
        with the current mode: Run text follows LEM/Seep/FEM; Build Mesh shows only
        in Seep/FEM (LEM needs no mesh); Seep/FEM Run is gated on a built mesh."""
        mode = self._mode
        self.act_run.setText({"lem": "Run &LEM…", "seep": "Run &Seep…",
                              "fem": "Run &FEM…"}.get(mode, "Run…"))
        open_ = self.doc.is_open
        busy = (self._runner is not None or self._seep_runner is not None
                or self._fem_runner is not None or self._mesh_busy)
        if mode == "lem":
            self.act_run.setEnabled(open_ and not busy)
            self.act_run.setToolTip("")
        else:
            has_mesh = open_ and self.doc.slope_data.get("mesh") is not None
            implemented = mode in self._run_implemented
            self.act_run.setEnabled(open_ and has_mesh and implemented and not busy)
            self.act_run.setToolTip(
                "Coming soon." if not implemented
                else "Build a mesh first (Build Mesh…)." if not has_mesh else "")
        # Meshing only applies to the FE workflows.
        self.act_build_mesh.setVisible(mode in ("seep", "fem"))
        self.act_build_mesh.setEnabled(open_ and mode in ("seep", "fem") and not busy)

    def run_current(self):
        """Dispatch the Run action by the current mode."""
        if self._mode == "lem":
            self.run_lem()
        elif self._mode == "seep":
            self.run_seep()
        elif self._mode == "fem":
            self.run_fem()

    def _populate_inputs_tree(self):
        d = self.doc.slope_data
        self.inputs_tree.clear()
        font_editable = None

        def add(name, value, category=None):
            item = QTreeWidgetItem(self.inputs_tree, [name, str(value)])
            if category is not None:
                item.setData(0, CATEGORY_ROLE, category)
                f = item.font(0)
                f.setUnderline(True)
                item.setFont(0, f)
                item.setToolTip(0, "Click to edit")
            return item

        sbc = d.get("seepage_bc") or {}
        profile_lines = d.get("profile_lines") or []
        polygons = d.get("polygons") or []
        add("Global parameters", "", category="global")
        add("Materials", len(d.get("materials", [])), category="materials")
        add("Profile lines", len(profile_lines),
            category="profile" if profile_lines else None)
        # Polygons are derived from profile lines for profile-based files (edit them
        # via the profile editor); only polygon-based files edit polygons directly.
        add("Polygons", len(polygons),
            category="polygons" if (polygons and not profile_lines) else None)
        add("Circles", len(d.get("circles") or []), category="circles")
        add("Non-circular pts", len(d.get("non_circ") or []), category="non_circ")
        add("Piezometric lines", len(d.get("piezo_line") or []), category="piezo")
        add("Distributed loads", len(d.get("dloads") or []), category="dloads")
        add("Reinforcement lines", len(d.get("reinforcement_lines") or []),
            category="reinforce")
        add("Piles", len(d.get("pile_lines") or []), category="piles")
        add("Seep BC", len(sbc.get("specified_heads", [])), category="seep_bc")
        self.inputs_tree.expandAll()

    # --- editing ---------------------------------------------------------
    def _on_tree_click(self, item, _column):
        category = item.data(0, CATEGORY_ROLE)
        if category:
            self.edit_category(category)

    def edit_category(self, category):
        editor = CATEGORY_EDITORS.get(category)
        if editor is None or not self.doc.is_open:
            return
        dlg = editor.build(self.doc.slope_data, self)
        if dlg.exec():
            self.doc.begin_edit()
            try:
                editor.apply(self.doc.slope_data, dlg)
            except Exception:
                traceback.print_exc()
            self.doc.commit_edit()        # -> re-render + mark dirty
            self._populate_inputs_tree()

    # --- meshing ---------------------------------------------------------
    def build_mesh(self):
        if not self.doc.is_open or self._mesh_busy:
            return
        dlg = BuildMeshDialog(self, defaults=self._last_mesh_opts)
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_mesh_opts = opts
        self._mesh_busy = True
        self._update_run_actions()    # disable Run/Build while meshing
        self.statusBar().showMessage("Building mesh …")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        # Runs on the persistent mesh thread (queued connection).
        self._mesh_requested.emit(self.doc.slope_data, opts)

    def _on_mesh_succeeded(self, mesh):
        self.doc.slope_data["mesh"] = mesh   # used by Inputs render, Seep and FEM
        self.doc.results["mesh"] = mesh
        # Persist alongside the .xlsx using the {stem}_mesh.json convention.
        if self.doc.path:
            try:
                from xslope.mesh import export_mesh_to_json
                stem = os.path.splitext(self.doc.path)[0]
                export_mesh_to_json(mesh, f"{stem}_mesh.json")
            except Exception:
                traceback.print_exc()
        self._show_mesh(mesh)
        self._render()                       # mesh now appears in the Inputs view
        n1d = len(mesh.get("elements_1d", []) or [])
        self.statusBar().showMessage(
            f"Mesh built — {len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements"
            + (f", {n1d} 1D elements." if n1d else "."))
        self._mesh_done()

    def _on_mesh_failed(self, message):
        QMessageBox.warning(self, "Build mesh failed", message)
        self.statusBar().showMessage("Build mesh failed.")
        self._mesh_done()

    def _mesh_done(self):
        self._mesh_busy = False
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self._update_run_actions()

    def _show_mesh(self, mesh):
        if self.mesh_canvas is None:
            self.mesh_canvas = MplCanvas(self)
            self.view_tabs.insertTab(1, self.mesh_canvas, "Mesh")
            panel = MeshDisplayPanel()
            panel.changed.connect(self._rerender_mesh)
            self.display_stack.addWidget(panel)
            self._display_panels[self.mesh_canvas] = panel
        self._rerender_mesh()
        self.view_tabs.setCurrentWidget(self.mesh_canvas)

    def _rerender_mesh(self):
        mesh = self.doc.results.get("mesh") or self.doc.slope_data.get("mesh")
        panel = self._display_panels.get(self.mesh_canvas)
        if mesh and panel and self.mesh_canvas is not None:
            try:
                self.mesh_canvas.render_mesh(
                    mesh, self.doc.slope_data.get("materials"), panel.options())
            except Exception:
                traceback.print_exc()

    # --- seepage ---------------------------------------------------------
    def run_seep(self):
        if not self.doc.is_open or self._seep_runner is not None:
            return
        if self.doc.slope_data.get("mesh") is None:
            QMessageBox.information(self, "Run Seepage",
                                    "Build a mesh first (Build Mesh…).")
            return
        dlg = RunSeepDialog(self, defaults=self._last_seep_opts,
                            has_bc2=bool(self.doc.slope_data.get("has_seepage_bc2")))
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_seep_opts = opts
        self.statusBar().showMessage("Running seepage …")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self._seep_runner = SeepRunner(self.doc.slope_data, opts, parent=self)
        self._seep_runner.succeeded.connect(self._on_seep_succeeded)
        self._seep_runner.failed.connect(self._on_seep_failed)
        self._seep_runner.finished.connect(self._on_seep_finished)
        self._update_run_actions()
        self._seep_runner.start()

    def _on_seep_succeeded(self, bundle):
        bc = bundle["options"].get("bc", 1)
        # Keep one solution per BC set so BC 1 and BC 2 (rapid drawdown) coexist
        # in separate tabs and can be compared side by side.
        self.doc.results.setdefault("seep_solutions", {})[bc] = bundle
        # Persist the solution next to the .xlsx ({stem}_seep.csv / _seep2.csv).
        if self.doc.path:
            try:
                from xslope.seep import export_seep_solution
                stem = os.path.splitext(self.doc.path)[0]
                suffix = "_seep.csv" if bc == 1 else f"_seep{bc}.csv"
                export_seep_solution(bundle["seep_data"], bundle["solution"], stem + suffix)
            except Exception:
                traceback.print_exc()
        self._show_seep_data(bundle["seep_data"], bc)
        self._show_seep_solution(bc)
        canvas = self.seep_solution_canvas.get(bc)
        if canvas is not None:
            self.view_tabs.setCurrentWidget(canvas)
        self.statusBar().showMessage(f"Seepage done (BC set {bc}).")

    def _on_seep_failed(self, message):
        QMessageBox.warning(self, "Seepage run failed", message)
        self.statusBar().showMessage("Seepage run failed.")

    def _on_seep_finished(self):
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        if self._seep_runner is not None:
            self._seep_runner.deleteLater()
            self._seep_runner = None
        self._update_run_actions()

    @staticmethod
    def _seep_tab_label(base, bc):
        return base if bc == 1 else f"{base} {bc}"

    def _show_seep_data(self, seep_data, bc=1):
        if bc not in self.seep_data_canvas:
            canvas = MplCanvas(self)
            self.seep_data_canvas[bc] = canvas
            self.view_tabs.addTab(canvas, self._seep_tab_label("Seep · Data", bc))
            panel = FeDataDisplayPanel()
            panel.changed.connect(lambda *_dummy, b=bc: self._rerender_seep_data(b))
            self.display_stack.addWidget(panel)
            self._display_panels[canvas] = panel
        self._rerender_seep_data(bc)

    def _rerender_seep_data(self, bc=1):
        bundle = self.doc.results.get("seep_solutions", {}).get(bc)
        canvas = self.seep_data_canvas.get(bc)
        panel = self._display_panels.get(canvas)
        if bundle and panel and canvas is not None:
            try:
                canvas.render_seep_data(bundle["seep_data"], panel.options())
            except Exception:
                traceback.print_exc()

    def _show_seep_solution(self, bc=1):
        if bc not in self.seep_solution_canvas:
            canvas = MplCanvas(self)
            self.seep_solution_canvas[bc] = canvas
            self.view_tabs.addTab(canvas, self._seep_tab_label("Seep · Solution", bc))
            panel = SeepDisplayPanel(self.doc.slope_data.get("materials"))
            panel.changed.connect(lambda *_dummy, b=bc: self._rerender_seep_solution(b))
            self.display_stack.addWidget(panel)
            self._display_panels[canvas] = panel
        self._rerender_seep_solution(bc)

    def _rerender_seep_solution(self, bc=1):
        """Re-render a cached seep solution (per BC set) with its Display options."""
        bundle = self.doc.results.get("seep_solutions", {}).get(bc)
        canvas = self.seep_solution_canvas.get(bc)
        panel = self._display_panels.get(canvas)
        if bundle and panel and canvas is not None:
            try:
                canvas.render_seep_solution(
                    bundle["seep_data"], bundle["solution"], panel.options())
            except Exception:
                traceback.print_exc()

    # --- FEM -------------------------------------------------------------
    def run_fem(self):
        if not self.doc.is_open or self._fem_runner is not None:
            return
        if self.doc.slope_data.get("mesh") is None:
            QMessageBox.information(self, "Run FEM", "Build a mesh first (Build Mesh…).")
            return
        dlg = RunFemDialog(self, defaults=self._last_fem_opts)
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_fem_opts = opts
        is_ssrm = opts["analysis"] == "ssrm"
        self.statusBar().showMessage("Running FEM …")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self._fem_runner = FemRunner(self.doc.slope_data, opts, parent=self)
        self._fem_runner.succeeded.connect(self._on_fem_succeeded)
        self._fem_runner.failed.connect(self._on_fem_failed)
        self._fem_runner.cancelled.connect(self._on_fem_cancelled)
        self._fem_runner.progress.connect(self._on_run_progress)
        self._fem_runner.finished.connect(self._on_fem_finished)
        if is_ssrm:                     # only SSRM supports cooperative cancel
            self.cancel_btn.setEnabled(True)
            self.cancel_btn.setVisible(True)
        self._update_run_actions()
        self._fem_runner.start()

    def _on_fem_succeeded(self, bundle):
        self.doc.results["fem_solution"] = bundle
        if self.doc.path:
            try:
                from xslope.fem import export_fem_solution
                export_fem_solution(bundle["fem_data"], bundle["solution"],
                                    os.path.splitext(self.doc.path)[0],
                                    meta={"FS": bundle.get("FS"),
                                          "analysis": bundle.get("analysis")})
            except Exception:
                traceback.print_exc()
        self._show_fem_data(bundle["fem_data"])
        self._show_fem_results()
        if self.fem_results_canvas is not None:
            self.view_tabs.setCurrentWidget(self.fem_results_canvas)
        if bundle.get("FS") is not None:
            self.statusBar().showMessage(f"FEM done — SSRM FS = {bundle['FS']:.3f}")
        else:
            conv = bundle["solution"].get("converged")
            self.statusBar().showMessage(f"FEM single solve done (converged={conv}).")

    def _on_fem_failed(self, message):
        QMessageBox.warning(self, "FEM run failed", message)
        self.statusBar().showMessage("FEM run failed.")

    def _on_fem_cancelled(self):
        self.statusBar().showMessage("Run cancelled.")

    def _on_fem_finished(self):
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.cancel_btn.setVisible(False)
        if self._fem_runner is not None:
            self._fem_runner.deleteLater()
            self._fem_runner = None
        self._update_run_actions()

    def _show_fem_data(self, fem_data):
        if self.fem_data_canvas is None:
            self.fem_data_canvas = MplCanvas(self)
            self.view_tabs.addTab(self.fem_data_canvas, "FEM · Data")
            panel = FeDataDisplayPanel(include_bc_symbol=True)
            panel.changed.connect(self._rerender_fem_data)
            self.display_stack.addWidget(panel)
            self._display_panels[self.fem_data_canvas] = panel
        self._rerender_fem_data()

    def _rerender_fem_data(self):
        bundle = self.doc.results.get("fem_solution")
        panel = self._display_panels.get(self.fem_data_canvas)
        if bundle and panel and self.fem_data_canvas is not None:
            try:
                self.fem_data_canvas.render_fem_data(bundle["fem_data"], panel.options())
            except Exception:
                traceback.print_exc()

    def _show_fem_results(self):
        if self.fem_results_canvas is None:
            self.fem_results_canvas = MplCanvas(self)
            self.view_tabs.addTab(self.fem_results_canvas, "FEM · Results")
            panel = FemResultsDisplayPanel()
            panel.changed.connect(self._rerender_fem_results)
            self.display_stack.addWidget(panel)
            self._display_panels[self.fem_results_canvas] = panel
        self._rerender_fem_results()

    def _rerender_fem_results(self):
        bundle = self.doc.results.get("fem_solution")
        panel = self._display_panels.get(self.fem_results_canvas)
        if bundle and panel and self.fem_results_canvas is not None:
            try:
                self.fem_results_canvas.render_fem_results(
                    bundle["fem_data"], bundle["solution"], panel.options(),
                    fs=bundle.get("FS"))
            except Exception:
                traceback.print_exc()

    # --- LEM analysis ----------------------------------------------------
    def run_lem(self):
        if not self.doc.is_open or self._runner is not None:
            return
        dlg = RunLemDialog(self, defaults=self._last_lem_opts,
                           slope_data=self.doc.slope_data)
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_lem_opts = opts
        self.act_run.setEnabled(False)
        verb = {"auto_search": "Searching", "reliability": "Running reliability"}.get(
            opts["analysis"], "Running")
        self.statusBar().showMessage(f"{verb} {opts['method']} …")
        # Show a busy bar immediately; reliability switches it to determinate.
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self.cancel_btn.setEnabled(True)
        self.cancel_btn.setVisible(True)
        self._runner = LemRunner(self.doc.slope_data, opts, parent=self)
        self._runner.succeeded.connect(self._on_lem_succeeded)
        self._runner.failed.connect(self._on_lem_failed)
        self._runner.cancelled.connect(self._on_lem_cancelled)
        self._runner.progress.connect(self._on_run_progress)
        self._runner.finished.connect(self._on_lem_finished)
        self._runner.start()

    def _cancel_run(self):
        runner = next((r for r in (self._runner, self._fem_runner)
                       if r is not None and r.isRunning()), None)
        if runner is not None:
            runner.cancel()
            self.cancel_btn.setEnabled(False)
            self.progress_bar.setRange(0, 0)   # back to busy while it winds down
            self.statusBar().showMessage("Cancelling…")

    def _on_run_progress(self, done, total, label):
        if total <= 0:                       # indeterminate
            self.progress_bar.setRange(0, 0)
        else:
            self.progress_bar.setRange(0, total)
            self.progress_bar.setValue(min(done, total))
        if label:
            self.statusBar().showMessage(f"{label}  ({done}/{total})" if total > 0 else label)

    def _on_lem_succeeded(self, bundle):
        self.doc.results["lem_solution"] = bundle
        if bundle.get("search"):
            self._show_search(bundle["search"])
        if bundle.get("reliability"):
            self._show_reliability(bundle["reliability"])
        if isinstance(bundle.get("results"), dict):
            self._show_solution(bundle)
        # Lead with the most specific result view produced by this run.
        lead = (self.reliability_canvas if bundle.get("reliability")
                else self.search_canvas if bundle.get("search")
                else self.solution_canvas)
        if lead is not None:
            self.view_tabs.setCurrentWidget(lead)
        if bundle.get("reliability"):
            r = bundle["reliability"]
            self.statusBar().showMessage(
                f"Reliability done — F_MLV = {r['F_MLV']:.3f}, "
                f"reliability = {r['reliability'] * 100:.2f}%, "
                f"Pf = {r['prob_failure'] * 100:.2f}%")
        else:
            res = bundle["results"]
            self.statusBar().showMessage(
                f"LEM done — {res.get('method')} FS = {res.get('FS'):.3f}")

    def _on_lem_failed(self, message):
        QMessageBox.warning(self, "LEM run failed", message)
        self.statusBar().showMessage("LEM run failed.")

    def _on_lem_cancelled(self):
        self.statusBar().showMessage("Run cancelled.")

    def _on_lem_finished(self):
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.cancel_btn.setVisible(False)
        if self._runner is not None:
            self._runner.deleteLater()
            self._runner = None
        self._update_run_actions()

    def _show_search(self, search):
        if self.search_canvas is None:
            self.search_canvas = MplCanvas(self)
            # Keep order Inputs → Search → Solution.
            self.view_tabs.insertTab(1, self.search_canvas, "LEM · Search")
            panel = SearchDisplayPanel()
            panel.changed.connect(self._rerender_search)
            self.display_stack.addWidget(panel)
            self._display_panels[self.search_canvas] = panel
        self._rerender_search()

    def _rerender_search(self):
        bundle = self.doc.results.get("lem_solution")
        search = bundle.get("search") if bundle else None
        panel = self._display_panels.get(self.search_canvas)
        if search and panel and self.search_canvas is not None:
            try:
                self.search_canvas.render_search(self.doc.slope_data, search, panel.options())
            except Exception:
                traceback.print_exc()

    def _show_reliability(self, reliability_data):
        if self.reliability_canvas is None:
            self.reliability_canvas = MplCanvas(self)
            self.view_tabs.insertTab(1, self.reliability_canvas, "LEM · Reliability")
        try:
            self.reliability_canvas.render_reliability(self.doc.slope_data, reliability_data)
        except Exception:
            traceback.print_exc()

    def _show_solution(self, bundle):
        if self.solution_canvas is None:
            self.solution_canvas = MplCanvas(self)
            self.view_tabs.addTab(self.solution_canvas, "LEM · Solution")
            panel = SolutionDisplayPanel()
            panel.changed.connect(self._rerender_solution)
            self.display_stack.addWidget(panel)
            self._display_panels[self.solution_canvas] = panel
        self._rerender_solution()

    def _rerender_solution(self):
        bundle = self.doc.results.get("lem_solution")
        panel = self._display_panels.get(self.solution_canvas)
        if (bundle and panel and self.solution_canvas is not None
                and isinstance(bundle.get("results"), dict)):
            try:
                self.solution_canvas.render_solution(
                    self.doc.slope_data, bundle["slice_df"],
                    bundle["failure_surface"], bundle["results"], panel.options())
            except Exception:
                traceback.print_exc()

    def _clear_result_tabs(self):
        """Drop result views (e.g. on opening another file) so they don't show
        stale results from the previous project."""
        single = ("mesh_canvas", "search_canvas", "solution_canvas",
                  "reliability_canvas", "fem_data_canvas", "fem_results_canvas")
        # The seep canvases are per-BC dicts; flatten them in with the rest.
        canvases = [getattr(self, a) for a in single]
        canvases += list(self.seep_data_canvas.values())
        canvases += list(self.seep_solution_canvas.values())
        for canvas in canvases:
            if canvas is not None:
                idx = self.view_tabs.indexOf(canvas)
                if idx >= 0:
                    self.view_tabs.removeTab(idx)
                panel = self._display_panels.pop(canvas, None)
                if panel is not None:
                    self.display_stack.removeWidget(panel)
                    panel.deleteLater()
                canvas.deleteLater()
        for a in single:
            setattr(self, a, None)
        self.seep_data_canvas = {}
        self.seep_solution_canvas = {}
        # Removing tabs may not fire currentChanged if Inputs was already active,
        # so point the Display dock at whatever tab remains current.
        self._show_display_for_tab(self.view_tabs.currentWidget())

    def _on_view_tab_changed(self, index):
        w = self.view_tabs.widget(index)
        if hasattr(w, "ensure_fitted"):
            w.ensure_fitted()
        self._show_display_for_tab(w)

    def _show_display_for_tab(self, widget):
        """Point the Display dock at the active view's options panel (or a
        placeholder when the view has none)."""
        panel = self._display_panels.get(widget)
        self.display_stack.setCurrentWidget(panel or self._display_placeholder)

    # --- save ------------------------------------------------------------
    def save(self):
        if not self.doc.path:
            return self.save_as()
        try:
            self.doc.save(template=None)   # edit in place, preserve formatting
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Save failed", str(exc))
            return
        self.statusBar().showMessage(f"Saved {os.path.basename(self.doc.path)}")

    def save_as(self):
        start = self.doc.path or (self._recent[0] if self._recent else "")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save As", start, "Excel files (*.xlsx)")
        if not path:
            return
        if not path.lower().endswith(".xlsx"):
            path += ".xlsx"
        try:
            self.doc.save(path=path, template=str(TEMPLATE))  # fresh file from template
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Save failed", str(exc))
            return
        self._add_recent(path)
        self._update_title()
        self.statusBar().showMessage(f"Saved {os.path.basename(path)}")

    # --- misc ------------------------------------------------------------
    def _update_title(self):
        if self.doc.is_open:
            name = os.path.basename(self.doc.path) if self.doc.path else "untitled"
            star = "*" if self.doc.dirty else ""
            self.setWindowTitle(f"{name}{star} — {APP_NAME}")
        else:
            self.setWindowTitle(APP_NAME)

    def _about(self):
        QMessageBox.about(
            self, f"About {APP_NAME}",
            f"<b>{APP_NAME}</b><br>Desktop GUI for the xslope slope-stability "
            f"engine.<br><br>Open an Excel input file to view its geometry and inputs.")

    def closeEvent(self, event):
        if self._runner is not None and self._runner.isRunning():
            self._runner.cancel()     # ask an in-flight run to stop, then wait briefly
            self._runner.wait(5000)
        if self._seep_runner is not None and self._seep_runner.isRunning():
            self._seep_runner.wait(10000)   # seepage has no cancel hook; let it finish
        if self._fem_runner is not None and self._fem_runner.isRunning():
            self._fem_runner.cancel()       # SSRM stops cooperatively
            self._fem_runner.wait(15000)
        # Stop the persistent mesh thread (lets an in-flight build finish first).
        self._mesh_thread.quit()
        self._mesh_thread.wait(10000)
        if self.doc.is_open and self.doc.dirty:
            res = QMessageBox.question(
                self, "Unsaved changes", "Save changes before closing?",
                QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
            if res == QMessageBox.Cancel:
                event.ignore()
                return
            if res == QMessageBox.Save:
                self.save()
                if self.doc.dirty:        # save failed or was cancelled
                    event.ignore()
                    return
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        super().closeEvent(event)
