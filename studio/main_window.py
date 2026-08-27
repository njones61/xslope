"""MainWindow — the XSLOPE Studio shell (Phase 1: read-only viewer).

Provides the app frame: menus, a dockable Inputs summary panel, a Log pane that
captures engine stdout/stderr, a mode selector (LEM / Seep / FEM), recent-files,
and the central embedded Matplotlib canvas. File -> Open loads an Excel input
file via ProjectDocument and renders the Inputs view. Editing, running analyses,
and saving arrive in later phases.
"""

from __future__ import annotations

import html
import json
import os
import re
import sys
import traceback
from urllib.parse import urlparse

from PySide6.QtCore import Qt, QObject, QSettings, QStandardPaths, QThread, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication, QButtonGroup, QDialog, QDockWidget, QFileDialog, QHBoxLayout,
    QLabel, QMainWindow, QMenu, QMessageBox, QPlainTextEdit, QProgressBar,
    QPushButton, QStackedWidget, QTabWidget, QToolBar, QToolButton, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout, QWidget,
)

from xslope.fileio import default_template_path
from xslope.package import (FEM_SOLUTION_SIDECARS as _FEM_SOLUTION_SIDECARS,
                            PACKAGE_EXT, is_package, pack, package_contents, unpack)

from . import links, urlscheme
from .canvas import MplCanvas
from .dialogs import (
    BuildMeshDialog, DxfImportDialog, GszImportDialog, ReliabilityDialog,
    RunFemDialog, RunLemDialog, RunSeepDialog, SensitivityDialog, Slide2ImportDialog,
    UnpackPackageDialog,
)
from .display_panels import (
    FeDataDisplayPanel, FemResultsDisplayPanel, InputsDisplayPanel,
    MeshDisplayPanel, ReliabilityDisplayPanel, ReliabilityMcDisplayPanel,
    SearchDisplayPanel, SeepDisplayPanel, SolutionDisplayPanel,
)
from .document import ProjectDocument
from .editors import CATEGORY_EDITORS
from .runners import (FemRunner, LemRunner, MeshWorker, ReliabilityRunner,
                      ReportRunner, SeepRunner, SensitivityRunner,
                      resume_cycle_gc, suspend_cycle_gc)
from .update_ui import UpdateController
from .welcome import WelcomeDialog

#: The app's display name — title bar, About box, menus. Matches the name the
#: installers give it on disk ("XSLOPE Studio.app", the Windows Add/Remove entry).
APP_NAME = "XSLOPE Studio"
ORG_NAME = "XSlope"
#: The QSettings identity, deliberately NOT the display name: it is the on-disk
#: settings path, and renaming it would orphan every existing user's preferences,
#: recent files and assistant configuration. studio/ai/config.py and the editors
#: open the same store with the same literal pair.
SETTINGS_APP = "XSlope Studio"
MAX_RECENT = 8
MODES = [("LEM", "lem"), ("Seepage", "seep"), ("FEM", "fem")]
# Blank template used to create files on Save As — the single copy bundled with
# the engine package (xslope/resources), so the GUI and library share one source.
TEMPLATE = default_template_path()
CATEGORY_ROLE = Qt.UserRole + 1

#: Every file ``xslope.fem.export_fem_solution`` writes beside a model, by suffix.
#: Imported from ``xslope.package``, where the whole companion convention lives, so
#: that the list this window deletes a solution by and the list the packager
#: attributes files with are the same list. (Re-exported here because the name is
#: part of this module's vocabulary and the editors import it from here.)
FEM_SOLUTION_SIDECARS = _FEM_SOLUTION_SIDECARS


# The engine prints color-emoji markers (🔁 ✅ ❌ ⚠️). Apple Color Emoji is a bitmap
# font with fixed, large metrics, so QPlainTextEdit — which sizes each line to its
# tallest glyph — stretches every emoji line, and shrinking the point size does not
# help. For the LOG PANE we map the markers to text-style symbols that render at the
# normal text height (the console keeps the color emoji), then strip any other color
# emoji as a fallback. The strip ranges skip U+2713–2717 so the ✓/✗ replacements
# survive, and stay above box-drawing/Greek/±/Δ so tables and math symbols are safe.
_LOG_EMOJI_MAP = {"🔁": "↻", "🔄": "↻", "✅": "✓", "❌": "✗", "⚠️": "!", "⚠": "!"}
_LOG_EMOJI_STRIP_RE = re.compile(
    "[\U0001F000-\U0001FAFF\U00002600-\U00002712\U00002718-\U000027BF"
    "\U00002B00-\U00002BFF\U0000FE00-\U0000FE0F]+")


def _log_sanitize(text):
    """Map the engine's color-emoji markers to text-style symbols (and strip any
    other color emoji) so the Log pane's lines stay a single text-height."""
    for k, v in _LOG_EMOJI_MAP.items():
        if k in text:
            text = text.replace(k, v)
    return _LOG_EMOJI_STRIP_RE.sub("", text)


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
                self._original.write(text)   # console: keep color emoji
            except Exception:
                pass
        text = text.rstrip("\n")
        if text:
            self._emitted.emit(_log_sanitize(text))   # log pane: text-style markers

    def flush(self):
        if self._original is not None:
            try:
                self._original.flush()
            except Exception:
                pass


class SolutionView(QWidget):
    """The LEM Solution tab: an admissibility-warning strip stacked above the
    result canvas.

    The solvers return ``results['warnings']`` — Duncan & Wright admissibility
    notes on an already-accepted solution (base tension on cohesionless slices,
    interslice tension, thrust line outside the slices; see
    ``solve._admissibility_warnings``). They reach the Log pane via the stdout tee,
    but a Studio user reading only the plot would take an inadmissible FS as a
    clean success. This strip surfaces them beside the solution: hidden when the
    list is empty, otherwise one amber line per note. It refreshes from the fresh
    results on every render, so it clears on the next (clean) solve.

    The view quacks like the ``MplCanvas`` it wraps for the two calls the main
    window makes on a result view — ``render_solution`` and ``ensure_fitted`` — so
    it slots into the existing tab / display / clear machinery with no
    special-casing beyond the one-line class swap that builds it.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.canvas = MplCanvas(self)
        self._warning = QLabel()
        self._warning.setWordWrap(True)
        self._warning.setTextInteractionFlags(Qt.TextSelectableByMouse)
        # Amber wash + deep-amber text, matching the app's existing warning hue
        # (#9a6700, used in the chat dock). Padding/radius follow the chat blocks;
        # the strip sizes to its wrapped text (no fixed height) and the canvas
        # takes the rest, so it self-adjusts to width and to the number of notes.
        self._warning.setStyleSheet(
            "QLabel {"
            " background-color: #fff4d6;"
            " color: #7a5200;"
            " border: 1px solid #e0b400;"
            " border-radius: 4px;"
            " padding: 6px 10px;"
            " }")
        self._warning.setVisible(False)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)
        layout.addWidget(self._warning)        # strip above the plot
        layout.addWidget(self.canvas, 1)       # canvas takes the remaining height

    def render_solution(self, *args, **kwargs):
        # The main window passes ``results`` as the 4th positional arg; refresh the
        # strip from it (empty -> hidden) before drawing, so it tracks each solve.
        results = args[3] if len(args) > 3 else kwargs.get("results")
        self._set_warnings(results)
        self.canvas.render_solution(*args, **kwargs)

    def ensure_fitted(self):
        self.canvas.ensure_fitted()

    def _set_warnings(self, results):
        warns = list(results.get("warnings") or []) if isinstance(results, dict) else []
        if not warns:
            self._warning.setVisible(False)
            self._warning.clear()
            return
        method = (results.get("method") or "").strip()
        plural = "s" if len(warns) > 1 else ""
        head = (f"{method.title()} — admissibility warning{plural}" if method
                else f"Admissibility warning{plural}")
        body = "<br>".join("&#8226;&nbsp;" + html.escape(w) for w in warns)
        self._warning.setText(f"<b>{html.escape(head)}</b><br>{body}")
        self._warning.setVisible(True)


class SweepCanvas(MplCanvas):
    """Result canvas for sensitivity / design sweeps.

    A thin MplCanvas subclass that renders the three sweep plots through the base
    canvas's ``_draw`` — so the Save… button, zoom/pan, and the pick machinery all
    come for free. It reuses the engine's ``plot_tornado`` / ``plot_sensitivity``
    as-is; the design view adds the target-FS crossing annotation on top (the
    engine plot draws the target line but not the interpolated crossing).
    """

    def render_tornado(self, result):
        from xslope.plot import plot_tornado
        self._draw(lambda fig: plot_tornado(result, fig=fig), dxf=False)

    def render_scaled(self, result, scaling="elasticity"):
        from xslope.plot import plot_scaled_sensitivity
        self._draw(lambda fig: plot_scaled_sensitivity(result, scaling=scaling,
                                                       fig=fig), dxf=False)

    def render_spider(self, sweeps):
        from xslope.plot import plot_spider
        self._draw(lambda fig: plot_spider(sweeps, fig=fig), dxf=False)

    def render_variance(self, result):
        from xslope.plot import plot_variance_pareto
        self._draw(lambda fig: plot_variance_pareto(result, fig=fig), dxf=False)

    def render_rank(self, result):
        from xslope.plot import plot_mc_rank_correlation
        self._draw(lambda fig: plot_mc_rank_correlation(result, fig=fig), dxf=False)

    def render_curve(self, df, target_fs=None):
        from xslope.plot import plot_sensitivity
        self._draw(lambda fig: plot_sensitivity(df, target_fs=target_fs, fig=fig),
                   dxf=False)

    def render_fs_vs_time(self, result, slope_data=None):
        from xslope.plot import plot_fs_vs_time
        self._draw(lambda fig: plot_fs_vs_time(result, slope_data=slope_data,
                                               fig=fig), dxf=False)

    def render_design(self, df, target_fs, summary):
        from xslope.plot import plot_sensitivity

        def draw(fig):
            plot_sensitivity(df, target_fs=target_fs, fig=fig)
            self._annotate_crossing(fig, target_fs, summary)

        self._draw(draw, dxf=False)

    def _annotate_crossing(self, fig, target_fs, summary):
        """Mark the interpolated output=target crossing (or an honest miss note).

        The drawing itself is ``xslope.plot.annotate_design_crossing`` — pure
        matplotlib, so anything that draws a design sweep outside Studio (the
        tutorial figure producer) puts the same crossing on the same curve
        rather than a lookalike of it.
        """
        from xslope.plot import annotate_design_crossing
        ax = self._main_axes()
        if ax is None:
            return
        annotate_design_crossing(ax, target_fs, summary)


class MainWindow(QMainWindow):
    # Emitted to hand a mesh build to the persistent mesh thread (queued).
    _mesh_requested = Signal(object, object)

    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        # Launch large — the canvas is the main feature. Cap to the screen.
        from PySide6.QtGui import QGuiApplication
        screen = QGuiApplication.primaryScreen()
        avail = screen.availableGeometry() if screen else None
        if avail is not None:
            self.resize(min(1680, avail.width() - 80), min(1040, avail.height() - 120))
        else:
            self.resize(1680, 1040)
        self.settings = QSettings(ORG_NAME, SETTINGS_APP)
        # Help → Check for Updates… and the once-a-day silent startup check. The
        # controller is created here (it needs ``settings``) but never checks
        # anything on its own: the manual item and app.py's deferred startup call
        # are the only two triggers, so constructing a MainWindow touches no network.
        self.updates = UpdateController(self)

        self.doc = ProjectDocument(self)
        self.doc.loaded.connect(self._on_loaded)
        self.doc.changed.connect(self._render)
        self.doc.dirty_changed.connect(lambda *_: self._update_title())

        # Central area: a tab strip of result views (plan §7). The Inputs view is
        # always present; the LEM Solution view is added after the first solve.
        self.canvas = MplCanvas(self)
        # Double-click an input on the Inputs canvas to edit it (plan §6/§8).
        # Only the Inputs view is wired; result-view canvases stay view-only.
        self.canvas.picked.connect(self._on_canvas_pick)
        self.canvas.set_pick_enabled(True)  # show the select cursor on the Inputs view
        self.mesh_canvas = None
        self.search_canvas = None
        self.solution_canvas = None
        self.reliability_canvas = None
        self.reliability_hist_canvas = None      # Monte Carlo FS histogram
        # Sensitivity / design study result tabs.
        self.sens_canvas = None            # tornado
        self.sens_curve_canvas = None      # click-through FS-vs-value curve
        self.design_canvas = None          # design curve + target crossing
        self.fs_time_canvas = None         # factor of safety vs time (transient)
        self.seep_data_canvas = {}        # bc set -> MplCanvas
        self.seep_solution_canvas = {}    # bc set -> MplCanvas
        self.transient_seep_view = None   # TransientSeepView (frames + play bar)
        self.fem_data_canvas = None
        self.fem_results_canvas = None
        self.fem_details_btn = None       # "1D Details…" on the FEM results toolbar
        self.fem_details_dlg = None       # the open (non-modal) details dialog
        self.view_tabs = QTabWidget()
        self.view_tabs.addTab(self.canvas, "Inputs")
        self.view_tabs.currentChanged.connect(self._on_view_tab_changed)
        self.setCentralWidget(self.view_tabs)

        self._mode = "lem"
        self._runner = None
        self._seep_runner = None
        self._fem_runner = None
        self._sens_runner = None
        self._rel_runner = None
        self._report_runner = None
        self._update_dl_runner = None    # in-flight update download (Help → Updates)
        self._mesh_busy = False
        # A stability run waiting on a transient re-march: ("lem"|"fem", options).
        # ``_pending_run_ok`` records whether that re-march actually produced frames,
        # so a failed or cancelled march never lets the analysis through.
        self._pending_run = None
        self._pending_run_ok = False
        self._run_implemented = {"lem", "seep", "fem"}   # modes whose Run is wired up
        self._last_lem_opts = {}
        self._last_sens_opts = {}          # keyed by engine mode (lem/fem/seep)
        self._last_rel_opts = {}           # keyed by engine mode (lem/fem)
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
        # Quit is not always a window close — Cmd-Q on macOS ends exec() without
        # sending one — and a QThread destroyed while it is still running aborts
        # the process. So the teardown hangs off the application's own quit as
        # well as off closeEvent; see stop_threads.
        app = QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(self.stop_threads)
        self._recent = [p for p in (self.settings.value("recent_files") or []) if p]

        self._display_panels = {}     # result tab widget -> its display-options panel
        self._make_inputs_dock()
        self._make_display_dock()
        self._make_log_dock()
        self._make_chat_dock()
        # Give both side columns the full window height (they own the bottom
        # corners), so the Log dock spans only the central canvas's width and the
        # left (Inputs/Display) and right (Assistant) columns run top to bottom.
        self.setCorner(Qt.BottomLeftCorner, Qt.LeftDockWidgetArea)
        self.setCorner(Qt.BottomRightCorner, Qt.RightDockWidgetArea)
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
        # The panel stack + a Styles button pinned at the bottom (plan §8a): styles
        # are project-global, so one button here opens the shared Styles dialog.
        display_container = QWidget()
        dv = QVBoxLayout(display_container)
        dv.setContentsMargins(0, 0, 0, 0)
        dv.addWidget(self.display_stack, 1)
        self.styles_btn = QPushButton("Styles…")
        self.styles_btn.setEnabled(False)
        self.styles_btn.setToolTip("Edit per-feature styles (color, hatch, opacity)")
        self.styles_btn.clicked.connect(self.open_styles_dialog)
        dv.addWidget(self.styles_btn)
        dock = QDockWidget("Display", self)
        dock.setObjectName("display_dock")
        dock.setWidget(display_container)
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
        # Fixed-width font so ASCII/tabulate grid tables (e.g. the reliability
        # results table) line up — a proportional font misaligns the columns.
        # Prefer an explicit terminal-style monospace (Menlo / SF Mono on macOS) so
        # the log's line spacing matches the console; the generic system fixed font
        # can render looser. Size is adjustable via the spinner in the title bar.
        from PySide6.QtGui import QFont
        self._log_font = QFont()
        self._log_font.setFamilies(["Menlo", "SF Mono", "Monaco", "DejaVu Sans Mono",
                                    "Consolas", "Courier New", "monospace"])
        self._log_font.setStyleHint(QFont.Monospace)
        self._log_font.setFixedPitch(True)
        self._log_font.setPointSize(12)
        self.log.setFont(self._log_font)
        self.log.setLineWrapMode(QPlainTextEdit.NoWrap)
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
        # Font-size control, just left of Clear.
        from PySide6.QtWidgets import QSpinBox
        font_spin = QSpinBox()
        font_spin.setRange(6, 28)
        font_spin.setValue(self._log_font.pointSize())
        font_spin.setSuffix(" pt")
        font_spin.setToolTip("Log font size")
        # Slightly smaller widget font so the "11 pt" reading isn't oversized.
        _sf = font_spin.font()
        _sf.setPointSize(11)
        font_spin.setFont(_sf)
        font_spin.valueChanged.connect(self._set_log_font_size)
        clear_btn = QToolButton()
        clear_btn.setText("Clear")
        clear_btn.setAutoRaise(True)
        clear_btn.setToolTip("Clear the log output")
        clear_btn.clicked.connect(self.log.clear)
        # Give both the same height (the taller of the two) so the spinner's native
        # up/down arrows sit centered next to the button, and vertically center them.
        _h = max(font_spin.sizeHint().height(), clear_btn.sizeHint().height())
        font_spin.setFixedHeight(_h)
        clear_btn.setFixedHeight(_h)
        row.addWidget(font_spin, 0, Qt.AlignVCenter)
        row.addWidget(clear_btn, 0, Qt.AlignVCenter)
        dock.setTitleBarWidget(title)
        self.addDockWidget(Qt.BottomDockWidgetArea, dock)
        self.log_dock = dock

    def _set_log_font_size(self, pt):
        """Resize the Log pane's (fixed-width) font from the title-bar spinner."""
        self._log_font.setPointSize(int(pt))
        self.log.setFont(self._log_font)

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
        # Keep the Assistant a relatively narrow right column so the canvas stays
        # the main feature.
        self.resizeDocks([self.chat_dock], [380], Qt.Horizontal)

    def _make_chat_dock(self):
        # AI assistant (Phase A spike) — a chat that drives the app/engine via an
        # in-process run_python tool. The anthropic dep is optional and imported
        # lazily, so the dock loads even without it (sending then reports it).
        from .ai.assistant import Assistant
        from .chat_dock import ChatDock
        self.assistant = Assistant(self)
        chat = ChatDock(self.assistant, self)
        dock = QDockWidget("Assistant", self)
        dock.setObjectName("chat_dock")
        dock.setWidget(chat)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        self.chat_dock = dock

    def refresh_inputs_view(self):
        """Re-render the Inputs canvas and the inputs tree after an external edit
        (e.g. the assistant mutated slope_data via run_python). Resyncs derived
        structures first, since the renderer reads those — e.g. reinforcement is
        plotted from the derived ``reinforce_lines``, not the ``reinforcement_lines``
        table the assistant edits, and geometry from ``polygons``/``ground_surface``."""
        if not self.doc.is_open:
            return
        sd = self.doc.slope_data
        try:
            from .editors import _resync_geometry
            _resync_geometry(sd)          # polygons/ground surface from profile_lines
        except Exception:
            traceback.print_exc()
        try:
            from xslope.fileio import build_reinforce_lines
            if sd.get("reinforcement_lines") is not None:
                sd["reinforce_lines"] = build_reinforce_lines(sd["reinforcement_lines"])
        except Exception:
            traceback.print_exc()
        self._render()
        self._populate_inputs_tree()

    def _install_log_capture(self):
        sys.stdout = _LogStream(self.log, sys.__stdout__)
        sys.stderr = _LogStream(self.log, sys.__stderr__)

    # --- actions / menus -------------------------------------------------
    def _make_actions(self):
        self.act_new = QAction("&New", self, shortcut=QKeySequence.New,
                               triggered=self.new_project)
        self.act_open = QAction("&Open…", self, shortcut=QKeySequence.Open,
                                triggered=self.open_dialog)
        self.act_import_dxf = QAction("&Import DXF…", self,
                                      triggered=self.import_dxf_dialog)
        self.act_import_gsz = QAction("Import &GeoStudio (SLOPE/W)…", self,
                                      triggered=self.import_gsz_dialog)
        self.act_import_slide2 = QAction("Import &Slide2…", self,
                                         triggered=self.import_slide2_dialog)
        self.act_import_rs2 = QAction("Import &RS2 (.fez)…", self,
                                      triggered=self.import_rs2_dialog)
        self.act_export_dxf = QAction("&Export Geometry (DXF)…", self, enabled=False,
                                      triggered=self.export_dxf_dialog)
        self.act_export_gsz = QAction("Export to GeoStudio (SLOPE/&W)…", self,
                                      enabled=False, triggered=self.export_gsz_dialog)
        self.act_export_pkg = QAction("Export &Project Package…", self, enabled=False,
                                      triggered=self.export_package_dialog)
        self.act_quit = QAction("&Quit", self, shortcut=QKeySequence.Quit,
                                triggered=self.close)
        self.act_undo = QAction("&Undo", self, shortcut=QKeySequence.Undo,
                                triggered=self._undo, enabled=False)
        self.act_redo = QAction("&Redo", self, shortcut=QKeySequence.Redo,
                                triggered=self._redo, enabled=False)
        self.act_about = QAction("&About", self, triggered=self._about)
        # The documentation is one click from anywhere in the app: it is the manual,
        # and it is also where the samples and the verification problems are listed,
        # so Studio keeps no in-app copy of either list.
        self.act_documentation = QAction("&Documentation", self,
                                         triggered=self.open_documentation)
        self.act_welcome = QAction("&Welcome", self, triggered=self.show_welcome)
        # Update items live in Help on every platform. macOS would otherwise be
        # free to re-home them by matching the text, so the role is explicit.
        self.act_check_updates = QAction("Check for &Updates…", self,
                                         triggered=self._check_updates)
        self.act_check_updates.setMenuRole(QAction.ApplicationSpecificRole)
        self.act_startup_check = QAction("Check for Updates at Startup", self,
                                         checkable=True)
        self.act_startup_check.setMenuRole(QAction.ApplicationSpecificRole)
        self.act_startup_check.setChecked(self.updates.startup_check_enabled())
        self.act_startup_check.toggled.connect(self.updates.set_startup_check_enabled)
        self.act_save = QAction("&Save", self, shortcut=QKeySequence.Save,
                                enabled=False, triggered=self.save)
        self.act_save_as = QAction("Save &As…", self, enabled=False, triggered=self.save_as)
        self.act_run = QAction("Run &LEM…", self, enabled=False, triggered=self.run_current)
        self.act_sensitivity = QAction("&Parametric…", self, enabled=False,
                                       triggered=self.run_sensitivity)
        self.act_reliability = QAction("&Reliability…", self, enabled=False,
                                       triggered=self.run_reliability)
        self.act_build_mesh = QAction("Build &Mesh…", self, enabled=False,
                                      triggered=self.build_mesh)
        self.act_report = QAction("&Generate Report…", self, enabled=False,
                                  triggered=self.generate_report)

    def _make_menus(self):
        mb = self.menuBar()

        m_file = mb.addMenu("&File")
        m_file.addAction(self.act_new)
        m_file.addAction(self.act_open)
        self.recent_menu = m_file.addMenu("Open &Recent")
        m_file.addSeparator()
        m_file.addAction(self.act_import_dxf)
        m_file.addAction(self.act_import_gsz)
        m_file.addAction(self.act_import_slide2)
        m_file.addAction(self.act_import_rs2)
        m_file.addAction(self.act_export_dxf)
        m_file.addAction(self.act_export_gsz)
        m_file.addAction(self.act_export_pkg)
        m_file.addSeparator()
        m_file.addAction(self.act_report)
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
        m_run.addAction(self.act_sensitivity)
        m_run.addAction(self.act_reliability)

        m_view = mb.addMenu("&View")
        m_view.addAction(self.inputs_dock.toggleViewAction())
        m_view.addAction(self.display_dock.toggleViewAction())
        m_view.addAction(self.log_dock.toggleViewAction())
        m_view.addAction(self.chat_dock.toggleViewAction())

        # Help reads outward first (the documentation and the welcome window),
        # then the update items, then About — which macOS moves to the application
        # menu on its own, by the action's text.
        m_help = mb.addMenu("&Help")
        m_help.addAction(self.act_documentation)
        m_help.addAction(self.act_welcome)
        m_help.addSeparator()
        m_help.addAction(self.act_check_updates)
        m_help.addAction(self.act_startup_check)
        m_help.addSeparator()
        m_help.addAction(self.act_about)

    def _make_toolbar(self):
        tb = QToolBar("Main", self)
        tb.setObjectName("main_toolbar")
        self.addToolBar(tb)
        tb.addAction(self.act_new)
        tb.addAction(self.act_open)
        tb.addSeparator()
        # Undo / Redo split-buttons: the button does a single step; the dropdown
        # lists the labeled history and jumps to a chosen point (plan §Phase 2).
        self.undo_btn = self._make_history_button(self.act_undo, redo=False)
        self.redo_btn = self._make_history_button(self.act_redo, redo=True)
        tb.addWidget(self.undo_btn)
        tb.addWidget(self.redo_btn)
        tb.addSeparator()
        self.mode_label = QLabel(" Mode: ")
        tb.addWidget(self.mode_label)
        tb.addWidget(self._make_mode_segments())
        tb.addSeparator()
        tb.addAction(self.act_build_mesh)
        tb.addAction(self.act_run)
        # The Parametric study lives on the toolbar too (Norm's ask); the action's
        # existing mode-visibility hides the button where it does not apply.
        tb.addAction(self.act_sensitivity)
        # Reliability is a sibling of Parametric: deterministic what-ifs vs the
        # probabilistic (β / probability-of-failure) study.
        tb.addAction(self.act_reliability)
        # The report is the end of the workflow, so it sits at the end of the strip.
        tb.addSeparator()
        tb.addAction(self.act_report)
        # macOS's native style draws text-only toolbar buttons in the larger system
        # font and ignores setFont; a stylesheet forces the size so New/Open/Run LEM
        # match the "Mode:" label. pointSizeF() is -1 for pixel-defined fonts.
        pt = self.mode_label.font().pointSizeF()
        css = f"QToolButton {{ font-size: {pt:g}pt; }}\n" if pt > 0 else ""
        tb.setStyleSheet(css + self._mode_segment_css())

    # --- analysis mode ----------------------------------------------------
    def _make_mode_segments(self):
        """The analysis-mode switch: one checkable segment per entry in MODES,
        joined into a single strip, exclusive, one click to switch.

        The segments live in their own container with zero spacing rather than
        being added to the toolbar one by one: a toolbar's gap between widgets is
        the platform style's to choose, and this control has to read as one strip
        on Windows and on macOS alike. Each segment also answers to Ctrl+N (⌘N on
        macOS) in MODES order, and says so in its tooltip.

        Clicks go through :meth:`set_mode_index`, which calls the same
        ``_on_mode_changed`` slot the drop-down this replaced was wired to, so
        every side effect of a mode change is unchanged."""
        box = QWidget(self)
        box.setObjectName("mode_segments")
        row = QHBoxLayout(box)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(0)
        self.mode_buttons = QButtonGroup(self)
        self.mode_buttons.setExclusive(True)
        self.mode_actions = []
        last = len(MODES) - 1
        for i, (label, key) in enumerate(MODES):
            btn = QToolButton(box)
            btn.setObjectName(f"mode_segment_{key}")
            btn.setText(label)
            btn.setCheckable(True)
            btn.setToolButtonStyle(Qt.ToolButtonTextOnly)
            btn.setFocusPolicy(Qt.NoFocus)   # the strip is a switch, not a tab stop
            btn.setProperty("segment",
                            "first" if i == 0 else "last" if i == last else "middle")
            tip = f"{label} analysis mode"
            if i < 9:                        # Ctrl+1..9; a 10th mode would need more
                seq = QKeySequence(f"Ctrl+{i + 1}")
                act = QAction(label, self)
                act.setShortcut(seq)
                act.triggered.connect(lambda _=False, n=i: self.set_mode_index(n))
                self.addAction(act)          # window-level: no menu entry needed
                self.mode_actions.append(act)
                tip += f"  ({seq.toString(QKeySequence.NativeText)})"
            btn.setToolTip(tip)
            self.mode_buttons.addButton(btn, i)
            row.addWidget(btn)
        self.mode_buttons.idClicked.connect(self.set_mode_index)
        self._show_mode(self._mode)
        return box

    def _mode_segment_css(self):
        """Stylesheet for the mode strip: flat, joined, modestly padded.

        The metrics are spelled out because the native styles disagree about a
        text-only tool button — macOS pads it generously, Windows draws its own
        hover and checked chrome — and a strip that is tight on one platform and
        loose on the other is the complaint this control answers. The colors come
        from the palette, so the strip follows a light or dark system theme."""
        pal = self.palette()
        line = pal.mid().color().name()
        text = pal.buttonText().color().name()
        hover = pal.midlight().color().name()
        sel = pal.highlight().color().name()
        sel_text = pal.highlightedText().color().name()
        return f"""
        QWidget#mode_segments QToolButton {{
            padding: 3px 12px; margin: 0px;
            border: 1px solid {line}; border-left-width: 0px; border-radius: 0px;
            background: transparent; color: {text};
        }}
        QWidget#mode_segments QToolButton[segment="first"] {{
            border-left-width: 1px;
            border-top-left-radius: 4px; border-bottom-left-radius: 4px;
        }}
        QWidget#mode_segments QToolButton[segment="last"] {{
            border-top-right-radius: 4px; border-bottom-right-radius: 4px;
        }}
        QWidget#mode_segments QToolButton:hover:!checked {{ background: {hover}; }}
        QWidget#mode_segments QToolButton:checked {{
            background: {sel}; color: {sel_text}; border-color: {sel};
        }}
        """

    def _show_mode(self, mode):
        """Highlight the segment for ``mode`` without running a mode change.

        Only ``idClicked`` reaches the handler, and setChecked does not emit it,
        so this is the silent set the drop-down needed blockSignals for."""
        keys = [m for _, m in MODES]
        if mode not in keys:
            return
        btn = self.mode_buttons.button(keys.index(mode))
        if btn is not None:
            btn.setChecked(True)

    def set_mode_index(self, index):
        """Switch to ``MODES[index]`` from a segment click or its shortcut.

        Re-selecting the mode already in force does nothing, matching the
        drop-down this replaced: it signalled on a change of index, not on every
        pick, so nothing re-rendered when the user chose what was already set."""
        if not 0 <= index < len(MODES):
            return
        mode = MODES[index][1]
        self._show_mode(mode)       # a shortcut has to move the highlight too
        if mode == self._mode:
            return
        self._on_mode_changed(index)

    # --- undo / redo history ---------------------------------------------
    def _make_history_button(self, action, redo):
        """A toolbar split-button: clicking runs ``action`` (single step); the
        dropdown lists the labeled history and jumps to a chosen point."""
        btn = QToolButton(self)
        btn.setDefaultAction(action)
        btn.setToolButtonStyle(Qt.ToolButtonTextOnly)
        btn.setPopupMode(QToolButton.MenuButtonPopup)
        menu = QMenu(btn)
        menu.aboutToShow.connect(lambda m=menu, r=redo: self._populate_history_menu(m, r))
        btn.setMenu(menu)
        return btn

    def _populate_history_menu(self, menu, redo):
        menu.clear()
        labels = self.doc.redo_labels() if redo else self.doc.undo_labels()
        verb = "Redo" if redo else "Undo"
        if not labels:
            menu.addAction(f"Nothing to {verb.lower()}").setEnabled(False)
            return
        for i, label in enumerate(labels):
            # Selecting entry i means jump past steps 0..i, i.e. i+1 steps.
            act = menu.addAction(f"{verb} {label}")
            act.triggered.connect(lambda _=False, n=i + 1, r=redo: self._history_multi(n, r))

    def _undo(self):
        self._history_step(self.doc.undo)

    def _redo(self):
        self._history_step(self.doc.redo)

    def _history_multi(self, n, redo):
        self._history_step((lambda: self.doc.redo_steps(n)) if redo
                           else (lambda: self.doc.undo_steps(n)))

    def _history_step(self, fn):
        """Run an undo/redo on the document, then reconcile derived state: the
        restored ``slope_data`` (mesh included) time-travels inside the snapshot, but
        cached LEM/seep/FEM *solutions* live outside it and are now stale, so drop
        them and re-sync the mesh/run gating to the restored geometry."""
        if not self.doc.is_open:
            return
        fn()                                  # emits changed -> _render (redraw + enable)
        self.invalidate_results(clear_mesh=False)   # drop stale solution tabs; keep mesh
        self.doc.results["mesh"] = self.doc.slope_data.get("mesh")
        self._update_run_actions()            # Seep/FEM gate follows the restored mesh
        self._populate_inputs_tree()

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
            self, "Open XSLOPE input file", start,
            "XSLOPE projects (*.xlsx *.xslz);;Excel files (*.xlsx);;"
            "Project packages (*.xslz);;All files (*)")
        if path:
            self.open_path(path)

    def open_path(self, path):
        if not self._confirm_discard():
            return
        if is_package(path):
            # A package is transport: unpack it to loose files first, then open the
            # workbook that comes out through this same path — so recent files, the
            # window title and everything downstream see an ordinary .xlsx.
            path = self._unpack_package(path)
            if not path:
                return
        try:
            self.doc.load(path)
        except Exception as exc:  # ValueError from the loader, or anything else
            traceback.print_exc()
            QMessageBox.critical(self, "Could not open file",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return
        self._add_recent(path)

    def _unpack_package(self, package):
        """Turn a .xslz package back into loose files; return the workbook to open.

        Returns None if the user cancels or the extraction fails. The destination is
        never reused or written over without being asked: the dialog offers the
        project already in an existing folder, or a fresh folder beside it.
        """
        try:
            names = package_contents(package)
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not open package",
                                 f"{os.path.basename(package)}:\n\n{exc}")
            return None
        dlg = UnpackPackageDialog(package, self)
        if not dlg.exec():
            return None
        dest, mode = dlg.chosen()
        workbook = os.path.join(dest, names[0])
        if mode == "existing":
            if not os.path.isfile(workbook):
                QMessageBox.critical(
                    self, "Could not open package",
                    f"{dest} holds no {names[0]}, so there is no project in it to "
                    f"open. Choose another folder, or Extract Fresh.")
                return None
            self.statusBar().showMessage(
                f"Opened the project already in {os.path.basename(dest)}")
            return workbook
        try:
            workbook = unpack(package, dest=dest)
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not unpack package",
                                 f"{os.path.basename(package)}:\n\n{exc}")
            return None
        self.statusBar().showMessage(
            f"Unpacked {os.path.basename(package)} to {dest}")
        return workbook

    # --- xslope:// links (the docs' "Open in Studio") ----------------------
    def download_dir(self):
        """Where a package fetched from a link is saved — the user's Downloads
        folder, or the home folder on a system that has none. Somewhere durable and
        visible, like every other destination in the packaging flow: what comes out
        of the package is a workbook the user is going to open in Excel."""
        for kind in (QStandardPaths.DownloadLocation, QStandardPaths.HomeLocation):
            path = QStandardPaths.writableLocation(kind)
            if path:
                return path
        return os.path.expanduser("~")

    def open_scheme_url(self, uri):
        """Act on an ``xslope://open?url=...`` link — the docs' Open in Studio.

        The link arrives from the operating system, which got it from a web page,
        which can be any page anywhere: the scheme is registered system-wide, not
        for one site. So the order here is the order of :mod:`studio.urlscheme`'s
        promises — understand the verb, check the host, ASK, and only then touch the
        network — and the last step is the ordinary :meth:`open_path`, so a package
        that arrives this way is unpacked and opened exactly like one that was
        double-clicked.

        Returns the downloaded package's path, or None if the link was refused, the
        user cancelled, or the download failed. Every refusal is shown to the user
        with the thing that was refused named in it.
        """
        try:
            _verb, url = urlscheme.parse_request(uri)
            urlscheme.check_url(url)
        except urlscheme.SchemeError as exc:
            QMessageBox.warning(self, "XSLOPE link", str(exc))
            return None
        except Exception as exc:
            # urlscheme's contract is that every refusal is a SchemeError, and a
            # check stands here anyway: whatever a link contains, the outcome has to
            # be a dialog. An exception reaching Qt from an OS-delivered link takes
            # the window down with the user's unsaved work in it.
            traceback.print_exc()
            QMessageBox.warning(self, "XSLOPE link",
                                f"{uri} could not be read as an XSLOPE link:\n\n{exc}")
            return None
        name = urlscheme.package_name(url)
        dest_dir = self.download_dir()
        answer = QMessageBox.question(
            self, "Open a project from the web?",
            f"XSLOPE Studio will download the project package\n\n    {name}\n\n"
            f"from {urlparse(url).hostname}, save it in {dest_dir}, and open it.\n\n"
            f"{url}\n\nDownload it?",
            QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Cancel)
        if answer != QMessageBox.Yes:
            return None
        self.statusBar().showMessage(f"Downloading {name}...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            package = urlscheme.download_package(url, dest_dir)
        except Exception as exc:
            if not isinstance(exc, urlscheme.SchemeError):
                traceback.print_exc()          # as above: never out to Qt
            QApplication.restoreOverrideCursor()
            self.statusBar().clearMessage()
            QMessageBox.critical(self, "Could not download the project", str(exc))
            return None
        finally:
            if QApplication.overrideCursor() is not None:
                QApplication.restoreOverrideCursor()
        self.statusBar().showMessage(f"Downloaded {name}")
        self.open_path(package)
        return package

    def _results_off_disk(self):
        """True if the session holds a mesh or solution with no sidecar beside the
        workbook yet — the state in which packaging would ship a project without the
        results the user is looking at. A run does not dirty the document (results
        are not edits), so the dirty flag alone does not see this."""
        if not self.doc.path:
            return True
        stem = os.path.splitext(self.doc.path)[0]
        results = self.doc.results
        if (self.doc.slope_data.get("mesh") is not None
                and not os.path.exists(f"{stem}_mesh.json")):
            return True
        for bc in (results.get("seep_solutions") or {}):
            suffix = "_seep.csv" if bc == 1 else f"_seep{bc}.csv"
            if not os.path.exists(f"{stem}{suffix}"):
                return True
        if results.get("transient_seep") and not os.path.exists(f"{stem}_tseep.csv"):
            return True
        if results.get("fem_solution") and not os.path.exists(f"{stem}_fem_nodes.csv"):
            return True
        return False

    def export_package_dialog(self):
        """Collect the project into one .xslz package — the workbook plus every
        sidecar written beside it — so it can be emailed, uploaded or handed over as
        a single file.

        The package is built from the files ON DISK, so a project with edits or with
        a mesh/solution this session has not written out yet is saved first. A
        package whose workbook disagreed with the results zipped beside it would be
        exactly the stale-sidecar trouble the format exists to prevent."""
        if not self.doc.is_open:
            return
        if self.doc.dirty or not self.doc.path or self._results_off_disk():
            res = QMessageBox.question(
                self, "Save before packaging?",
                "A package holds the project as it stands on disk, and this one has "
                "work that is not written out yet. Save it first?",
                QMessageBox.Save | QMessageBox.Cancel)
            if res != QMessageBox.Save:
                return
            self.save()
            if self.doc.dirty or not self.doc.path:
                return          # the save failed or was cancelled
        stem = os.path.splitext(self.doc.path)[0]
        # The dialog appends the extension itself (setDefaultSuffix) rather than this
        # code appending it afterwards: a name typed without one has to become
        # "name.xslz" BEFORE the dialog checks whether that file exists, or the
        # overwrite confirmation is asked about a file nobody is going to write.
        dlg = QFileDialog(self, "Export project package", stem + PACKAGE_EXT,
                          f"XSLOPE project packages (*{PACKAGE_EXT})")
        dlg.setAcceptMode(QFileDialog.AcceptSave)
        dlg.setDefaultSuffix(PACKAGE_EXT.lstrip("."))
        if not dlg.exec():
            return
        selected = dlg.selectedFiles()
        if not selected:
            return
        path = selected[0]
        try:
            pack(self.doc.path, dest=path, overwrite=True)   # the dialog confirmed
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not export package",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return
        n = len(package_contents(path))
        self.statusBar().showMessage(
            f"Exported {os.path.basename(path)} — {n} file(s)")

    def import_dxf_dialog(self):
        """Import a DXF into a fresh project (confirm discard first, like Open). A
        wizard maps each DXF layer to an input feature — material zone, profile
        line, piezo line, distributed load, reinforcement, failure circles, or
        ignore. Geometry populates the features; non-geometric properties come in
        as editable placeholders. Left unsaved so the user fills those in and
        Saves As."""
        if not self._confirm_discard():
            return
        start = os.path.dirname(self._recent[0]) if self._recent else ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Import DXF", start, "DXF drawings (*.dxf);;All files (*)")
        if not path:
            return
        try:
            layers, warnings = self.doc.read_dxf_layers(path)
        except ImportError:
            traceback.print_exc()
            QMessageBox.critical(
                self, "DXF support not installed",
                "Reading and writing DXF files needs the 'ezdxf' package, which "
                "isn't installed in this environment.\n\nInstall it with:\n\n"
                "    pip install ezdxf\n\n(or reinstall with the 'cad'/'gui' extra: "
                "pip install \"xslope[gui]\"), then restart XSLOPE Studio.")
            return
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not import DXF",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return
        from xslope.cad import suggest_dxf_target
        wizard = DxfImportDialog(layers, suggest_dxf_target, self)
        if not wizard.exec():
            return
        try:
            notes = self.doc.build_from_dxf_mapping(layers, wizard.result())
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not import DXF",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return
        d = self.doc.slope_data
        for w in warnings:                         # surface to the Log pane
            print(f"DXF import warning: {w}")
        self.statusBar().showMessage(
            f"Imported {os.path.basename(path)} — "
            f"{len(d.get('materials') or [])} material(s), "
            f"{len(d.get('profile_lines') or d.get('polygons') or [])} geometry item(s), "
            f"{len(d.get('circles') or [])} circle(s). Fill in properties, then Save As.")
        allnotes = list(notes) + list(warnings)
        if allnotes:
            QMessageBox.information(
                self, "DXF imported",
                "Imported with notes:\n\n• " + "\n• ".join(allnotes) +
                "\n\nSee the Log pane for details.")

    def import_gsz_dialog(self):
        """Import a GeoStudio SLOPE/W model (.gsz) into a fresh project (confirm
        discard first, like Open).

        A .gsz needs no mapping wizard — its geometry, materials and water conditions
        are already identified — so the only prompt is which analysis to import, since
        a file usually holds several and they can differ in materials. Whatever xslope
        cannot represent (reinforcement, loads, SLOPE/W's search definition) comes back
        as caveats and is shown, not dropped quietly. Left unsaved so the user reviews
        it and Saves As."""
        if not self._confirm_discard():
            return
        start = os.path.dirname(self._recent[0]) if self._recent else ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Import GeoStudio", start,
            "GeoStudio projects (*.gsz);;All files (*)")
        if not path:
            return
        try:
            gsz, analyses = self.doc.read_gsz_analyses(path)
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not import GeoStudio file",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return

        if len(analyses) == 1:
            analysis_id = analyses[0]["id"]        # nothing to choose
        else:
            picker = GszImportDialog(analyses, self)
            if not picker.exec():
                return
            analysis_id = picker.result()

        try:
            caveats = self.doc.build_from_gsz(gsz, analysis_id)
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not import GeoStudio file",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return

        d = self.doc.slope_data
        for c in caveats:                          # surface to the Log pane
            print(f"GeoStudio import note: {c}")
        self.statusBar().showMessage(
            f"Imported {os.path.basename(path)} — "
            f"{len(d.get('materials') or [])} material(s), "
            f"{len(d.get('polygons') or [])} zone(s), "
            f"{len(d.get('circles') or [])} circle(s). Review, then Save As.")
        if caveats:
            QMessageBox.information(
                self, "GeoStudio imported",
                "Imported with notes:\n\n• " + "\n• ".join(caveats) +
                "\n\nSee the Log pane for details.")

    def import_slide2_dialog(self):
        """Import a Rocscience Slide2 model (.sli/.slim/.slmd) into a fresh project
        (confirm discard first, like Open).

        Like a .gsz, a Slide2 file needs no mapping wizard — its geometry, materials
        and water conditions are already identified — so the only prompt is which
        scenario to import, since a .slmd usually holds several (a base case plus
        variants) and they can differ in geometry and water as well as materials.
        Whatever xslope cannot represent (supports/anchors, loads, Slide2's search
        definition) comes back as caveats and is shown, not dropped quietly. Most
        Slide2 tutorial models are search-only, so the import routinely arrives with
        no failure circle — the document still opens for editing (the caveat says
        so) and the user defines circles afterward. Left unsaved so the user reviews
        it and Saves As."""
        if not self._confirm_discard():
            return
        start = os.path.dirname(self._recent[0]) if self._recent else ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Import Slide2", start,
            "Slide2 models (*.slmd *.slim *.sli);;All files (*)")
        if not path:
            return
        try:
            d, scenarios = self.doc.read_slide2_scenarios(path)
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not import Slide2 file",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return

        if len(scenarios) == 1:
            scenario = scenarios[0]["index"]        # nothing to choose
        else:
            picker = Slide2ImportDialog(scenarios, self)
            if not picker.exec():
                return
            scenario = picker.result()

        try:
            caveats = self.doc.build_from_slide2(d, scenario)
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not import Slide2 file",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return

        d2 = self.doc.slope_data
        for c in caveats:                          # surface to the Log pane
            print(f"Slide2 import note: {c}")
        self.statusBar().showMessage(
            f"Imported {os.path.basename(path)} — "
            f"{len(d2.get('materials') or [])} material(s), "
            f"{len(d2.get('polygons') or [])} zone(s), "
            f"{len(d2.get('circles') or [])} circle(s). Review, then Save As.")
        if caveats:
            QMessageBox.information(
                self, "Slide2 imported",
                "Imported with notes:\n\n• " + "\n• ".join(caveats) +
                "\n\nSee the Log pane for details.")

    def import_rs2_dialog(self):
        """Import a Rocscience RS2 finite-element model (.fez) into a fresh project
        (confirm discard first, like Open).

        A .fez holds exactly one model — RS2 has no notion of bundling several
        scenarios/analyses in one file the way a .gsz or .slmd does — so there is no
        scenario picker: the only prompt is which file to open. Geometry, materials
        and water conditions import directly; what RS2 defines and xslope cannot
        (its Shear-Strength-Reduction settings, joints, reinforcement, loads) comes
        back as caveats and is shown, not dropped quietly. RS2's slope-stability
        result is a finite-element SSR field, not a limit-equilibrium search, so the
        import NEVER carries a failure surface — the document still opens for
        editing (the caveat says so) and the user defines circles afterward, same as
        a search-only Slide2 import. Left unsaved so the user reviews it and Saves
        As."""
        if not self._confirm_discard():
            return
        start = os.path.dirname(self._recent[0]) if self._recent else ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Import RS2", start, "RS2 models (*.fez);;All files (*)")
        if not path:
            return
        try:
            caveats = self.doc.build_from_fez(path)
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not import RS2 file",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return

        d = self.doc.slope_data
        for c in caveats:                          # surface to the Log pane
            print(f"RS2 import note: {c}")
        self.statusBar().showMessage(
            f"Imported {os.path.basename(path)} — "
            f"{len(d.get('materials') or [])} material(s), "
            f"{len(d.get('polygons') or [])} zone(s), "
            f"{len(d.get('circles') or [])} circle(s). Review, then Save As.")
        if caveats:
            QMessageBox.information(
                self, "RS2 imported",
                "Imported with notes:\n\n• " + "\n• ".join(caveats) +
                "\n\nSee the Log pane for details.")

    def export_gsz_dialog(self):
        """Export the current model to a GeoStudio SLOPE/W project (.gsz) — material
        zones become regions, materials become Mohr-Coulomb materials, and a piezo
        line becomes a piezometric surface.

        A .gsz cannot carry everything xslope models: failure surfaces, reinforcement,
        piles and loads have no mapping. Those come back as caveats and are shown, so
        the user knows what to re-create on the GeoStudio side."""
        from xslope.geostudio import export_gsz
        stem = (os.path.splitext(os.path.basename(self.doc.path))[0]
                if self.doc.path else "model")
        start = (os.path.join(os.path.dirname(self.doc.path), stem + ".gsz")
                 if self.doc.path else stem + ".gsz")
        path, _ = QFileDialog.getSaveFileName(
            self, "Export to GeoStudio", start, "GeoStudio projects (*.gsz)")
        if not path:
            return
        if not path.lower().endswith(".gsz"):
            path += ".gsz"
        try:
            caveats = export_gsz(self.doc.slope_data, path, analysis_name=stem)
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not export to GeoStudio", f"{exc}")
            return
        for c in caveats:
            print(f"GeoStudio export note: {c}")
        self.statusBar().showMessage(f"Exported {os.path.basename(path)}")
        if caveats:
            QMessageBox.information(
                self, "Exported to GeoStudio",
                "Exported with notes:\n\n• " + "\n• ".join(caveats) +
                "\n\nSee the Log pane for details.")

    def export_dxf_dialog(self):
        """Export the current model's geometry to a structured (layered) DXF via
        the engine's ``export_dxf`` — material zones on per-material layers, and
        profile lines / circles / reinforcement / dloads / piezo on their reserved
        feature layers. Unlike the per-view canvas Save→DXF (which writes the
        rendered picture), this is the clean geometry export meant for re-import."""
        if not self.doc.is_open:
            return
        stem = os.path.splitext(os.path.basename(self.doc.path))[0] if self.doc.path else "geometry"
        start = os.path.join(os.path.dirname(self.doc.path), stem + ".dxf") if self.doc.path else stem + ".dxf"
        path, _ = QFileDialog.getSaveFileName(
            self, "Export geometry (DXF)", start, "DXF drawings (*.dxf)")
        if not path:
            return
        if not path.lower().endswith(".dxf"):
            path += ".dxf"
        try:
            from xslope.cad import export_dxf
            export_dxf(self.doc.slope_data, path)
        except ImportError:
            traceback.print_exc()
            QMessageBox.critical(
                self, "DXF support not installed",
                "Writing DXF files needs the 'ezdxf' package.\n\nInstall it with:\n\n"
                "    pip install ezdxf\n\nthen restart XSLOPE Studio.")
            return
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Could not export DXF",
                                 f"{os.path.basename(path)}:\n\n{exc}")
            return
        self.statusBar().showMessage(f"Exported geometry to {os.path.basename(path)}")

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
        self.act_export_dxf.setEnabled(True)
        self.act_export_gsz.setEnabled(True)
        self.act_export_pkg.setEnabled(True)
        self.styles_btn.setEnabled(True)
        self.assistant.reset()        # new project -> fresh conversation
        self._clear_result_tabs()
        # Restore saved solutions first so the default mode can see whether an FEM
        # solution exists, then pick the mode that fits this file.
        self._load_solution_sidecars()
        self._mode = self._default_mode()
        self._show_mode(self._mode)   # set silently; we render explicitly below
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
        """Pick the mode that fits the loaded file. In priority order: FEM if a mesh
        and a restored FEM solution are present; Seep if the materials carry only
        seepage properties (conductivity, no strength); FEM if a mesh is present and
        the materials define the elastic properties FEM needs (E, nu) — the mesh is
        there for a stress analysis, not for seepage pore pressures; otherwise LEM."""
        sd = self.doc.slope_data
        has_mesh = sd.get("mesh") is not None
        materials = sd.get("materials", [])
        if has_mesh and "fem_solution" in self.doc.results:
            return "fem"
        if self._materials_seep_only(materials):
            return "seep"
        if has_mesh and self._materials_fem_ready(materials):
            return "fem"
        return "lem"

    @staticmethod
    def _materials_fem_ready(materials):
        """True when every material carries the elastic properties an FEM stress
        analysis needs — Young's modulus E > 0 (fileio defaults E to 0 when blank,
        so a non-zero E on all materials means FEM was deliberately set up). This is
        the same precondition build_fem_data enforces, so it flags a file that is
        ready to run FEM rather than one that merely has a mesh for seepage."""
        def num(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return 0.0

        if not materials:
            return False
        return all(num(m.get("E")) > 0 for m in materials)

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
        self._restore_transient_sidecar(mesh, stem)
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
            bundle = {"seep_data": seep_data, "solution": solution,
                      "options": {"bc": bc}}
            self.doc.results.setdefault("seep_solutions", {})[bc] = bundle
            # load_slope_data has normally already read this same sidecar into
            # seep_u/seep_u2. Applying it again here costs nothing and closes the
            # gap where it did not: the loader's read is best-effort and warns on
            # failure, and it reads nothing at all when no material declares
            # u = seep at load time.
            self._apply_seep_field(bundle, bc)
            self._show_seep_data(seep_data, bc)
            self._show_seep_solution(bc)
            print(f"Restored saved seepage solution (BC set {bc}) from "
                  f"{os.path.basename(path)}.")

    def _restore_transient_sidecar(self, mesh, stem):
        # Restore a saved transient run ({stem}_tseep.csv + _tseep_meta.json) into the
        # Seep · Transient tab (frames + play bar). Best-effort: a mismatched or
        # unreadable sidecar is skipped, not fatal.
        if not os.path.exists(f"{stem}_tseep.csv"):
            return
        try:
            from xslope.seep import build_seep_data, import_transient_solution
            seep_data = build_seep_data(mesh, self.doc.slope_data, seep_bc=1)
            loaded = import_transient_solution(seep_data, stem)
        except Exception:
            traceback.print_exc()   # streams to the Log pane; load still succeeds
            return
        bundle = {"mode": "transient", "seep_data": seep_data,
                  "transient": loaded, "frames": loaded["frames"],
                  "options": {"mode": "transient"}}
        self.doc.results["transient_seep"] = bundle
        self._show_transient_seep(bundle, keep_index=False)
        print(f"Restored saved transient seepage solution "
              f"({len(loaded['frames'])} frame(s)) from "
              f"{os.path.basename(stem)}_tseep.csv.")

    def _restore_fem_sidecar(self, mesh, stem):
        if not os.path.exists(f"{stem}_fem_nodes.csv"):
            return
        try:
            from xslope.fem import build_fem_data, import_fem_solution, import_fem_meta
            fem_data = build_fem_data(self.doc.slope_data, mesh)
            try:
                solution = import_fem_solution(fem_data, stem)
            except ValueError as exc:
                # Stale sidecar: saved against a different mesh than the one now on
                # disk (e.g. the mesh was rebuilt in an older build that left the
                # sidecar behind). Skip it quietly — no traceback — rather than
                # failing the load; a fresh solve + Save re-syncs it.
                print(f"Skipping stale FEM solution sidecar: {exc}")
                return
            meta = import_fem_meta(stem) or {}
            # The strength-reduction factor the saved field was solved at, where
            # the file records one. It is NOT the factor of safety: F is the last
            # trial the solve converged at, FS the bracket midpoint the run
            # reported, and the subplot titles name them differently. Falling back
            # to FS printed it as "rendered at last converged F" over a trial that
            # was never run — rs2_28a records FS = 1.606 and no F at all, and its
            # only recorded trial is the 1.847 in its failure sidecar. Absent F is
            # absent: _fs_title then titles the panel without one.
            F_saved = meta.get("F")
            if F_saved is not None:
                solution["F"] = F_saved
        except Exception:
            traceback.print_exc()
            return
        # import_fem_solution nests the reconstructed at-failure snapshot under
        # solution["failure_solution"] (absent for sidecars saved before it
        # existed, or an SSRM run with no capture). Lift it onto the bundle
        # itself — that's the shape _rerender_fem_results/_on_fem_succeeded use
        # to thread it into the canvas render opts (see render_fem_results).
        failure_solution = solution.pop("failure_solution", None)
        # What was run is "analysis" in the sidecar Studio writes and
        # "analysis_type" in the ones the benchmark figures were built with. Both
        # are read, because a restored strength reduction run that arrives as
        # "loaded" is one whose factor of safety the report will not state.
        analysis = meta.get("analysis") or meta.get("analysis_type") or "loaded"
        self.doc.results["fem_solution"] = {
            "fem_data": fem_data, "solution": solution, "FS": meta.get("FS"),
            "analysis": analysis,
            "failure_solution": failure_solution,
            # The whole record the run kept of itself, carried the way a live run
            # carries it, so a report generated from a reopened model says how the
            # answer was reached and not only what it was. This is the shape
            # xslope.report.solutions_from_sidecars builds for the same files.
            "meta": meta}
        self._show_fem_data(fem_data)
        self._show_fem_results()
        fs = meta.get("FS")
        fs_note = f" (SSRM FS = {fs:.3f})" if isinstance(fs, (int, float)) else ""
        print(f"Restored saved FEM solution from {os.path.basename(stem)}_fem_*.csv{fs_note}.")

    def _render(self):
        if not self.doc.is_open:
            self.canvas.clear()
            return
        sd = self.doc.slope_data
        if not (sd.get("profile_lines") or sd.get("polygons")):
            # Empty project (e.g. freshly created with New) — no geometry to draw
            # yet. Leave the canvas blank until the user adds inputs.
            self.canvas.clear()
        else:
            try:
                self.canvas.render_inputs(sd, mode=self._mode,
                                          opts=self.inputs_panel.options(),
                                          style=self.doc.style or None)
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
                or self._fem_runner is not None or self._sens_runner is not None
                or self._rel_runner is not None
                or self._report_runner is not None or self._mesh_busy)
        has_mesh = open_ and self.doc.slope_data.get("mesh") is not None
        # The Parametric study has a version for every mode (LEM: FS; FEM: FS via
        # SSRM; Seep: discharge q). Always visible; the FEM/Seep sweeps run on the
        # mesh, so gate those on a built mesh exactly like Run.
        self.act_sensitivity.setVisible(True)
        if mode == "lem":
            self.act_sensitivity.setEnabled(open_ and not busy)
            self.act_sensitivity.setToolTip("")
        else:
            self.act_sensitivity.setEnabled(open_ and has_mesh and not busy)
            self.act_sensitivity.setToolTip(
                "" if has_mesh else "Build a mesh first (Build Mesh…).")
        # Reliability is probabilistic; it applies to LEM (Taylor / Monte Carlo) and
        # FEM (Taylor / SSRM), but not to a seepage analysis. Hidden in Seep mode;
        # the FEM sweep runs on the mesh, so gate it on a built mesh like Run.
        self.act_reliability.setVisible(mode in ("lem", "fem"))
        if mode == "lem":
            self.act_reliability.setEnabled(open_ and not busy)
            self.act_reliability.setToolTip("")
        elif mode == "fem":
            self.act_reliability.setEnabled(open_ and has_mesh and not busy)
            self.act_reliability.setToolTip(
                "" if has_mesh else "Build a mesh first (Build Mesh…).")
        if mode == "lem":
            self.act_run.setEnabled(open_ and not busy)
            self.act_run.setToolTip("")
        else:
            implemented = mode in self._run_implemented
            self.act_run.setEnabled(open_ and has_mesh and implemented and not busy)
            self.act_run.setToolTip(
                "Coming soon." if not implemented
                else "Build a mesh first (Build Mesh…)." if not has_mesh else "")
        # Meshing only applies to the FE workflows.
        self.act_build_mesh.setVisible(mode in ("seep", "fem"))
        self.act_build_mesh.setEnabled(open_ and mode in ("seep", "fem") and not busy)
        # A report documents a solved model, so it waits for a solution and says
        # so. Any engine's solution is one: a seepage run and a strength
        # reduction run each get a section of their own.
        solved = bool(self.report_solutions()) if open_ else False
        self.act_report.setEnabled(open_ and solved and not busy)
        self.act_report.setToolTip(
            "" if solved else "Run an analysis first — a report documents results.")

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
        # A project is profile-based unless it has polygons but no profile lines.
        # An empty (new) project can start down either path, so both rows are
        # live until the first geometry lands; the first profile line makes the
        # file profile-based, the first polygon makes it polygon-based.
        profile_based = bool(profile_lines) or not polygons
        empty_geometry = not profile_lines and not polygons
        add("Profile lines", len(profile_lines),
            category="profile" if profile_based else None)
        # Polygons are derived from profile lines for profile-based files (edit them
        # via the profile editor); polygon-based files — and an empty project —
        # edit polygons directly.
        add("Polygons", len(polygons),
            category="polygons" if (not profile_based or empty_geometry) else None)
        # SSR zones are polygon-sheet rows with sentinel Mat IDs, edited in the same
        # dialog as the material zones (appended after them). Listed separately
        # because they are analysis overlays, not geometry — the Polygons count must
        # keep meaning "material zones". Shown only when the file has some.
        _ssr_zones = d.get("ssr_zones") or []
        if _ssr_zones:
            add("SSR zones", len(_ssr_zones),
                category="polygons" if not profile_based else None)
        add("Circles", len(d.get("circles") or []), category="circles")
        add("Non-circular pts", len(d.get("non_circ") or []), category="non_circ")
        add("Piezometric lines", len(d.get("piezo_line") or []), category="piezo")
        add("Distributed loads", len(d.get("dloads") or []), category="dloads")
        add("Reinforcement lines", len(d.get("reinforcement_lines") or []),
            category="reinforce")
        add("Line loads", len(d.get("line_loads") or []), category="line_loads")
        add("Piles", len(d.get("pile_lines") or []), category="piles")
        add("Seep BC", len(sbc.get("specified_heads", [])), category="seep_bc")
        # Transient seepage (the optional tseep sheet). Shown for every file — like
        # Seep BC — with an "on/off" summary; the editor carries the enable affordance.
        add("Transient", "on" if d.get("tseep") else "off", category="transient")
        self.inputs_tree.expandAll()

    # --- editing ---------------------------------------------------------
    def _on_tree_click(self, item, _column):
        category = item.data(0, CATEGORY_ROLE)
        if category:
            self.edit_category(category)

    def open_styles_dialog(self):
        """Edit per-feature styles (plan §8a). Previews live on the canvas; OK keeps
        the change (and marks dirty → written to {stem}_styles.json on Save), Cancel
        restores the prior style."""
        if not self.doc.is_open:
            return
        import copy as _copy
        from .styles_dialog import StylesDialog
        orig = _copy.deepcopy(self.doc.style)

        def preview(style):
            self.doc.style = style
            self._render()                 # Inputs
            self._rerender_styled_results()  # LEM solution/search/reliability

        dlg = StylesDialog(self.doc.slope_data.get("materials") or [],
                           self.doc.style, preview, self)
        if dlg.exec():
            self.doc.set_style(dlg.result())     # mark dirty + re-render Inputs
            self._rerender_styled_results()
        else:
            self.doc.style = orig
            self._render()
            self._rerender_styled_results()

    def _rerender_styled_results(self):
        """Re-render the styled result views (so a Styles change previews there too)."""
        self._rerender_solution()
        self._rerender_search()
        self._rerender_reliability()
        self._rerender_mesh()
        for bc in list(self.doc.results.get("seep_solutions", {})):
            self._rerender_seep_data(bc)
            self._rerender_seep_solution(bc)
        self._rerender_transient_seep()
        self._rerender_fem_data()

    def _on_canvas_pick(self, x, y, tol):
        """Open the editor for the input feature the user double-clicked on the
        Inputs canvas. The hit-test maps the click back to a slope_data object and
        returns its editor category and index (plan §6/§8)."""
        if not self.doc.is_open:
            return
        from .picking import pick_category
        hit = pick_category(self.doc.slope_data, x, y, max(tol, 1e-9), mode=self._mode)
        if hit:
            category, index = hit
            self.edit_category(category, select=index)

    def edit_category(self, category, select=None):
        editor = CATEGORY_EDITORS.get(category)
        if editor is None or not self.doc.is_open:
            return
        # Pass the picked index to editors that can pre-highlight it (profile /
        # polygon dialogs); others keep the simple build(slope_data, parent) shape.
        import inspect
        if "select" in inspect.signature(editor.build).parameters:
            dlg = editor.build(self.doc.slope_data, self, select=select)
        else:
            dlg = editor.build(self.doc.slope_data, self)
        if dlg.exec():
            mesh_before = self.mesh_signature(self.doc.slope_data)
            self.doc.begin_edit(f"Edit {self.EDIT_LABELS.get(category, category)}")
            try:
                editor.apply(self.doc.slope_data, dlg)
            except Exception:
                # A failed apply must not leave a partial edit on the undo stack —
                # restore the snapshot and bail (the document is unchanged).
                traceback.print_exc()
                self.doc.rollback_edit()
                return
            self.doc.commit_edit()        # -> re-render + mark dirty
            self._populate_inputs_tree()
            # inputs changed -> solution is stale; if the geometry changed, the
            # mesh is stale too (it embeds the profile/polygon/reinforce/pile lines).
            geom_changed = self.mesh_signature(self.doc.slope_data) != mesh_before
            self.invalidate_results(clear_mesh=geom_changed)

    # Human labels for the undo-history dropdown, keyed by CATEGORY_EDITORS key.
    EDIT_LABELS = {
        "global": "Global Parameters", "materials": "Materials", "circles": "Circles",
        "non_circ": "Non-Circular Surface", "piezo": "Piezometric Lines",
        "dloads": "Distributed Loads", "seep_bc": "Seepage BC", "piles": "Piles",
        "reinforce": "Reinforcement", "line_loads": "Line Loads",
        "profile": "Profile Lines", "polygons": "Polygons",
        "transient": "Transient Seepage",
    }

    # Source inputs whose change makes the mesh stale (the domain geometry plus
    # the reinforcement/pile lines baked in as mesh constraint lines).
    MESH_KEYS = ("profile_lines", "polygons", "max_depth", "reinforcement_lines",
                 "pile_lines")

    @staticmethod
    def mesh_signature(sd):
        """JSON signature of the mesh-affecting inputs, to detect geometry edits."""
        def clean(o):
            if hasattr(o, "wkt"):
                return o.wkt
            try:
                import numpy as np
                if isinstance(o, np.ndarray):
                    return o.tolist()
                if isinstance(o, np.generic):
                    return o.item()
            except Exception:
                pass
            if isinstance(o, dict):
                return {k: clean(v) for k, v in o.items()}
            if isinstance(o, (list, tuple)):
                return [clean(x) for x in o]
            return o
        def geometry_only(key, rows):
            # A reinforcement or pile row enters the mesh as a constraint line
            # and nothing more: only its endpoints (and how many rows there are)
            # can change the mesh. Editing S, D, Tmax, Appl or any other property
            # on a row must not throw the mesh away.
            if key not in ("reinforcement_lines", "pile_lines") or not rows:
                return clean(rows)
            try:
                return [[clean(r.get("x1")), clean(r.get("y1")),
                         clean(r.get("x2")), clean(r.get("y2"))] for r in rows]
            except Exception:
                return clean(rows)
        try:
            return json.dumps({k: geometry_only(k, sd.get(k)) for k in MainWindow.MESH_KEYS},
                              sort_keys=True, default=str)
        except Exception:
            return None

    def invalidate_results(self, clear_mesh=False):
        """Inputs changed (via an editor or the assistant), so any cached analysis
        solution is stale — drop the solution result tabs and their cached results.
        ``clear_mesh`` (set when the geometry changed) also drops the mesh, which
        is then rebuilt explicitly. Leaves the user on a valid view (Inputs), which
        is why an edit visibly refreshes."""
        if not self.doc.is_open:
            return
        single = ["search_canvas", "solution_canvas", "reliability_canvas",
                  "sens_canvas", "sens_curve_canvas", "design_canvas",
                  "fs_time_canvas",
                  "fem_data_canvas", "fem_results_canvas"]
        if clear_mesh:
            single.append("mesh_canvas")
        canvases = [getattr(self, a) for a in single]
        canvases += list(self.seep_data_canvas.values())
        canvases += list(self.seep_solution_canvas.values())
        if self.transient_seep_view is not None:
            canvases.append(self.transient_seep_view)
        removed = False
        for canvas in canvases:
            if canvas is not None:
                idx = self.view_tabs.indexOf(canvas)
                if idx >= 0:
                    self.view_tabs.removeTab(idx)
                    removed = True
                panel = self._display_panels.pop(canvas, None)
                if panel is not None:
                    self.display_stack.removeWidget(panel)
                    panel.deleteLater()
                canvas.deleteLater()
        for a in single:
            setattr(self, a, None)
        self.seep_data_canvas = {}
        self.seep_solution_canvas = {}
        self.transient_seep_view = None
        for key in ("lem_solution", "seep_solutions", "transient_seep",
                    "fem_solution", "design", "sensitivity", "fs_vs_time"):
            self.doc.results.pop(key, None)
        if clear_mesh:
            self.doc.slope_data["mesh"] = None
            self.doc.results.pop("mesh", None)
            self._update_run_actions()       # Seep/FEM Run re-gated on a built mesh
        if removed:
            self._show_display_for_tab(self.view_tabs.currentWidget())
            what = "results and mesh" if clear_mesh else "result(s)"
            self.statusBar().showMessage(f"Inputs changed — cleared the now-stale {what}.")

    def _file_defaults(self, which):
        """Dialog defaults for ``which`` ("lem" | "mesh" | "fem"), seeded from the
        run options the OPEN FILE declares (template v19, main sheet and circles
        sheet).

        Precedence is: this session's last dialog choice > the file > the dialog's
        own hardcoded default. The session cache is checked first because an
        explicit choice the user has already made in this window must never be
        overwritten by the file on the next open of the dialog; before any run this
        session the cache is empty, so the file's declaration is what greets them.
        A key the file leaves blank is simply absent, and the dialog falls through
        to its own default exactly as before.
        """
        cache = {"lem": self._last_lem_opts, "mesh": self._last_mesh_opts,
                 "fem": self._last_fem_opts}[which]
        d = dict(cache)
        sd = self.doc.slope_data or {}

        def seed(key, value, ok=lambda v: v is not None):
            if key not in d and ok(value):
                d[key] = value

        if which == "lem":
            # 'all' is a batch-report mode with no Studio equivalent — the run
            # dialog solves one method — so it seeds nothing and the default stands.
            seed("method", sd.get("lem_method"),
                 lambda v: bool(v) and v != "all")
            seed("num_slices", sd.get("num_slices"))
        elif which == "mesh":
            seed("element_type", sd.get("element_type"), bool)
            # The 1D element size is a model input the dialog writes back, so the
            # model wins over the last run's value: the box always shows what the
            # file states, including a value cleared since the last build.
            d["element_size_1d"] = sd.get("element_size_1d")
            if "target_size" not in d and sd.get("target_size") is not None:
                d["target_size"] = float(sd["target_size"])
                # An explicit size in the file means the file MEANT that size, so
                # turn auto-sizing off — otherwise the value would sit in a
                # disabled box and the mesh would come out at the auto size.
                d.setdefault("auto_size", False)
        else:                                   # fem
            seed("F_min", sd.get("ssrm_f_min"))
            seed("F_max", sd.get("ssrm_f_max"))
            seed("k0", sd.get("k0"))
            seed("tension_srf", sd.get("tension_srf"))
            seed("side_bc", sd.get("side_bc"))
        return d

    # --- meshing ---------------------------------------------------------
    def build_mesh(self):
        if not self.doc.is_open or self._mesh_busy:
            return
        dlg = BuildMeshDialog(self, defaults=self._file_defaults("mesh"))
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_mesh_opts = opts
        # The 1D element size is a model input (main!D20), not a per-run choice, so
        # an entry goes back into the model — saved with the file and undoable like
        # any other edit. Nothing is recorded when the value is unchanged, so opening
        # the dialog and building does not by itself dirty the document.
        if opts.get("element_size_1d") != self.doc.slope_data.get("element_size_1d"):
            self.doc.begin_edit("1D element size")
            self.doc.slope_data["element_size_1d"] = opts["element_size_1d"]
            self.doc.commit_edit()
        self._mesh_busy = True
        self._update_run_actions()    # disable Run/Build while meshing
        self.statusBar().showMessage("Building mesh …")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        # Runs on the persistent mesh thread (queued connection). A mesh build is a
        # background run like any other, so it takes the cyclic collector with it —
        # see runners.suspend_cycle_gc. The RunnerThread base does this for the
        # QThread runners; this one is a worker on a shared thread, so the pair is
        # explicit: taken here, returned in _mesh_done, both on the GUI thread.
        suspend_cycle_gc()
        self._mesh_requested.emit(self.doc.slope_data, opts)

    def _on_mesh_succeeded(self, mesh):
        self.doc.slope_data["mesh"] = mesh   # used by Inputs render, Seep and FEM
        self.doc.results["mesh"] = mesh
        # A new mesh invalidates any previously computed seep/FEM solution (and the
        # LEM solution, whose pore pressures come from seepage): they were built on
        # the old node/element set. Drop the stale in-memory results and their tabs
        # so they can't be re-shown or re-saved against the new mesh.
        self.invalidate_results()
        # Persist the mesh alongside the .xlsx ({stem}_mesh.json) — this write is
        # eager (before Save), so the on-disk solution sidecars, which are stale vs
        # the new mesh, must be removed now too. Otherwise the next load pairs the
        # new mesh with an old {stem}_fem_nodes.csv and import_fem_solution raises a
        # node-count mismatch.
        if self.doc.path:
            try:
                from xslope.mesh import export_mesh_to_json
                stem = os.path.splitext(self.doc.path)[0]
                export_mesh_to_json(mesh, f"{stem}_mesh.json")
                self._remove_solution_sidecars(stem)
            except Exception:
                traceback.print_exc()
        self._show_mesh(mesh)
        self._render()                       # mesh now appears in the Inputs view
        e1d = mesh.get("elements_1d")            # a numpy array when present;
        n1d = len(e1d) if e1d is not None else 0  # `array or []` raises (ambiguous truth)
        self.statusBar().showMessage(
            f"Mesh built — {len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements"
            + (f", {n1d} 1D elements." if n1d else "."))
        self._mesh_done()

    def _on_mesh_failed(self, message):
        QMessageBox.warning(self, "Build mesh failed", message)
        self.statusBar().showMessage("Build mesh failed.")
        self._mesh_done()

    @staticmethod
    def _remove_solution_sidecars(stem):
        """Delete the on-disk seep/FEM solution sidecars for ``stem``. Called when
        the mesh is rebuilt (they no longer match the new node/element set) so the
        persisted files stay self-consistent; a fresh solve re-writes them on Save."""
        for name in ("_seep.csv", "_seep2.csv") + FEM_SOLUTION_SIDECARS:
            path = f"{stem}{name}"
            try:
                if os.path.exists(path):
                    os.remove(path)
            except Exception:
                traceback.print_exc()

    def _mesh_done(self):
        self._mesh_busy = False
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        resume_cycle_gc()             # the pair of build_mesh's suspension
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
                    mesh, self.doc.slope_data.get("materials"), panel.options(),
                    style=self.doc.style or None)
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
                            has_bc2=bool(self.doc.slope_data.get("has_seepage_bc2")),
                            has_tseep=bool(self.doc.slope_data.get("tseep")),
                            slope_data=self.doc.slope_data, document=self.doc)
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_seep_opts = opts
        # Rapid-drawdown stage times are model inputs (edited under Inputs →
        # Transient); the transient runner reads them from the document's tseep data,
        # so the dialog carries none and there is nothing to persist here.
        transient = opts.get("mode") == "transient"
        self.statusBar().showMessage("Running seepage …")
        # A transient run marches over a known duration, so it reports determinate
        # progress (t / duration) and can be cancelled; a steady solve is seconds-
        # scale and stays a busy indicator with no cancel.
        self.progress_bar.setRange(0, 100 if transient else 0)
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(True)
        self._seep_runner = SeepRunner(self.doc.slope_data, opts, parent=self)
        self._seep_runner.succeeded.connect(self._on_seep_succeeded)
        self._seep_runner.failed.connect(self._on_seep_failed)
        self._seep_runner.cancelled.connect(self._on_seep_cancelled)
        self._seep_runner.progress.connect(self._on_run_progress)
        self._seep_runner.finished.connect(self._on_seep_finished)
        if transient:
            self.cancel_btn.setEnabled(True)
            self.cancel_btn.setVisible(True)
        self._update_run_actions()
        self._seep_runner.start()

    def _on_seep_succeeded(self, bundle):
        if bundle.get("mode") == "transient":
            self._on_transient_seep_succeeded(bundle)
            return
        bc = bundle["options"].get("bc", 1)
        # Keep one solution per BC set so BC 1 and BC 2 (rapid drawdown) coexist
        # in separate tabs and can be compared side by side.
        self.doc.results.setdefault("seep_solutions", {})[bc] = bundle
        # A solved field belongs to the MODEL, not only to a results tab: a stability
        # run with u = seep reads slope_data['seep_u'], and so does the gate that
        # decides whether the run may start. Applied here, right where the solve
        # lands, for the same reason the mesh is: it is a computed artifact the rest
        # of the session depends on, written directly (no undo step, no dirty flag)
        # and persisted to its own sidecar just below.
        self._apply_seep_field(bundle, bc)
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

    def _apply_seep_field(self, bundle, bc):
        """Place a steady solution's nodal pore pressures into the model (seep_u for
        BC set 1, seep_u2 for set 2).

        Shared by a fresh solve and by the restore of a saved sidecar, so a field is
        in the model on the same terms however it got here. A refusal — the field was
        computed on a different mesh — is reported and the previous field is left
        alone rather than replaced by one that does not fit the geometry.
        """
        try:
            from xslope.seep import apply_steady_stability_field
            apply_steady_stability_field(self.doc.slope_data, bundle["solution"], bc=bc)
        except Exception as e:
            traceback.print_exc()
            QMessageBox.warning(
                self, "Seepage solution",
                f"The seepage solution was computed but could not be attached to the "
                f"model, so a stability run will not read it.\n\n{e}")

    def _on_seep_failed(self, message):
        self._pending_run = None
        QMessageBox.warning(self, "Seepage run failed", message)
        self.statusBar().showMessage("Seepage run failed.")

    def _on_seep_cancelled(self):
        # A cancelled transient march stores no partial result — just reset to idle.
        # A re-march started for a stability run takes that run down with it.
        if self._pending_run is not None:
            self._pending_run = None
            self.statusBar().showMessage("Seepage re-run cancelled — the analysis was not run.")
            return
        self.statusBar().showMessage("Run cancelled.")

    def _on_seep_finished(self):
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.cancel_btn.setVisible(False)
        self.cancel_btn.setEnabled(True)
        if self._seep_runner is not None:
            self._seep_runner.deleteLater()
            self._seep_runner = None
        self._update_run_actions()
        # A re-march runs so that a stability analysis can read an instant the old
        # solution never saved; with the new frames in hand, start that analysis.
        pending, self._pending_run = self._pending_run, None
        if pending and self._pending_run_ok:
            kind, opts = pending
            self._pending_run_ok = False
            {"lem": self._start_lem, "fem": self._start_fem}[kind](opts)

    def _start_remarch(self, times, pending):
        """Re-run the transient march with ``times`` added to the save schedule, then
        start the stability run that needs them.

        A time that is not a saved frame has no pore-pressure field yet, and a field
        blended between two frames is not a solution of anything — so it is COMPUTED,
        never interpolated. That costs a full re-solve, which on a long march is
        minutes, so the cost is stated before it is incurred and the march reports
        progress and can be cancelled like any other.
        """
        if self._seep_runner is not None:
            QMessageBox.information(self, "Seepage time",
                                    "A seepage run is already in progress.")
            return
        if self.doc.slope_data.get("mesh") is None:
            QMessageBox.information(self, "Seepage time",
                                    "Build a mesh first (Build Mesh…).")
            return
        unit = self.doc.slope_data.get("time_unit")
        listed = ", ".join(f"t = {t:g}" + (f" {unit}" if unit else "") for t in times)
        if QMessageBox.question(
                self, "Re-run the transient seepage analysis",
                f"The loaded transient solution has no saved frame at {listed}.\n\n"
                f"Pore pressures are never interpolated between frames, so the "
                f"transient seepage analysis will be re-run with {'this instant' if len(times) == 1 else 'these instants'} "
                f"added to the save schedule, and the analysis will start when it "
                f"finishes. It is a full re-solve — seconds on a short run, "
                f"minutes on a long one. It can be cancelled.\n\nRe-run now?",
                QMessageBox.Yes | QMessageBox.Cancel,
                QMessageBox.Yes) != QMessageBox.Yes:
            return
        self._pending_run = pending
        self._pending_run_ok = False
        opts = {"mode": "transient", "bc": 1,
                "tol": (self._last_seep_opts or {}).get("tol", 1e-4),
                "extra_save_times": [float(t) for t in times]}
        self.statusBar().showMessage("Re-running transient seepage …")
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(True)
        self._seep_runner = SeepRunner(self.doc.slope_data, opts, parent=self)
        self._seep_runner.succeeded.connect(self._on_seep_succeeded)
        self._seep_runner.failed.connect(self._on_seep_failed)
        self._seep_runner.cancelled.connect(self._on_seep_cancelled)
        self._seep_runner.progress.connect(self._on_run_progress)
        self._seep_runner.finished.connect(self._on_seep_finished)
        self.cancel_btn.setEnabled(True)
        self.cancel_btn.setVisible(True)
        self._update_run_actions()
        self._seep_runner.start()

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
                canvas.render_seep_data(bundle["seep_data"], panel.options(),
                                        style=self.doc.style or None)
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
                # Default the flow-net base material to the solver-side pick
                # (nearest-to-squares) until the user chooses one by hand.
                from xslope.plot_seep import flownet_base_material
                panel.suggest_base_mat(flownet_base_material(
                    bundle["seep_data"], bundle["solution"]))
                canvas.render_seep_solution(
                    bundle["seep_data"], bundle["solution"], panel.options(),
                    style=self.doc.style or None)
            except Exception:
                traceback.print_exc()

    # --- transient seepage (frames + play bar) ---------------------------
    def _on_transient_seep_succeeded(self, bundle):
        self.doc.results["transient_seep"] = bundle
        # Persist the frame bundle next to the .xlsx ({stem}_tseep.csv + meta).
        if self.doc.path:
            try:
                from xslope.seep import export_transient_solution
                stem = os.path.splitext(self.doc.path)[0]
                export_transient_solution(
                    bundle["seep_data"], bundle["transient"], stem,
                    input_file=self.doc.path,
                    mesh_file=f"{os.path.basename(stem)}_mesh.json")
            except Exception:
                traceback.print_exc()
        self._show_transient_seep(bundle, keep_index=False)
        if self.transient_seep_view is not None:
            self.view_tabs.setCurrentWidget(self.transient_seep_view)
        self.statusBar().showMessage(
            f"Transient seepage done — {len(bundle['frames'])} saved frame(s).")
        # A re-march started for a stability run: the frames it was waiting on now
        # exist, so the run may go ahead once this thread has wound down.
        if self._pending_run is not None:
            self._pending_run_ok = True

    def _show_transient_seep(self, bundle, keep_index=True):
        if self.transient_seep_view is None:
            from .transient import TransientSeepView
            view = TransientSeepView(self)
            self.transient_seep_view = view
            self.view_tabs.addTab(view, "Seep · Transient")
            # transient=True: Water levels default ON (the pool drops through playback)
            # and the flow-net-only Flow lines / Base material controls are omitted (a
            # transient storage-release state has no stream function / flow net).
            panel = SeepDisplayPanel(self.doc.slope_data.get("materials"),
                                     transient=True)
            panel.changed.connect(self._rerender_transient_seep)
            self.display_stack.addWidget(panel)
            self._display_panels[view] = panel
        view = self.transient_seep_view
        panel = self._display_panels.get(view)
        view.set_frames(bundle["seep_data"], bundle["frames"],
                        opts_getter=(panel.options if panel else None),
                        style_getter=(lambda: self.doc.style or None),
                        keep_index=keep_index)

    def _rerender_transient_seep(self):
        if self.transient_seep_view is not None:
            self.transient_seep_view.rerender()

    def transient_solution(self):
        """The loaded transient seepage solution, or ``None``. The run dialogs need
        it to offer a frame to choose from."""
        bundle = self.doc.results.get("transient_seep") or {}
        return bundle.get("transient")

    def transient_current_time(self):
        """The instant the transient results viewer is showing, or ``None`` when no
        transient results tab is open."""
        view = self.transient_seep_view
        return None if view is None else view.current_time

    def _apply_transient_analysis_frame(self, rapid=False, time=None):
        """Point the pending run at one instant of the loaded transient solution.

        Places that frame's pore pressures into ``slope_data['seep_u']`` (and, for a
        rapid drawdown with stage times set, the stage_1/stage_2 frames into
        ``seep_u`` / ``seep_u2``) so an LEM or FEM run with ``u = seep`` reads the
        chosen instant. No transient result → untouched, and the steady field stands.

        ``time`` is the Run dialog's choice. Without one the model's own
        ``stability_time`` decides, and with that blank too the last saved frame does
        — the meaning of a blank stability time, stated in the log either way.

        The write is TRANSACTIONAL and failure is LOUD. Staging a rapid drawdown can
        fail for real reasons — stage times with no saved frames, one stage blank,
        stage 1 after stage 2 — and the earlier version wrote outside the document's
        edit machinery and swallowed the exception, so a failed staging left the
        previous stage-1 field in place and the run reported a factor of safety for
        pore pressures nobody asked for. Now a failure rolls the model back and
        returns False, and the caller abandons the run.

        Returns True when the run may proceed.
        """
        sol = self.transient_solution()
        if sol is None:
            return True
        from xslope.seep import apply_transient_stability_frame
        sd = self.doc.slope_data
        ts = sd.get("tseep") or {}
        staged = rapid and ts.get("stage_1") is not None and ts.get("stage_2") is not None
        self.doc.begin_edit("Seepage frame")
        try:
            apply_transient_stability_frame(sd, sol, time=time, rapid=staged)
        except Exception as e:
            traceback.print_exc()
            self.doc.rollback_edit()
            QMessageBox.warning(
                self, "Transient seepage frame",
                f"The run was not started: the pore pressures for this analysis could "
                f"not be read out of the transient seepage solution.\n\n{e}")
            self.statusBar().showMessage("Run not started — no usable seepage frame.")
            return False
        self.doc.commit_edit()
        return True

    def _remarch_times_needed(self, opts):
        """The instants this run needs that the loaded solution has not saved.

        Empty when nothing has to be recomputed — which is the dominant case, since a
        saved frame and the viewer's frame are both already computed states."""
        from xslope.seep import has_transient_frame
        sol = self.transient_solution()
        if sol is None:
            return []
        wanted = []
        st = opts.get("stage_times")
        if st and st.get("stage_1") is not None and st.get("stage_2") is not None:
            wanted += [float(st["stage_1"]), float(st["stage_2"])]
        elif opts.get("seep_time") and opts["seep_time"].get("time") is not None:
            wanted.append(float(opts["seep_time"]["time"]))
        # t = 0 is the initial condition and is always saved, and nothing past the
        # duration can be reached — neither is a re-march's business.
        dur = (self.doc.slope_data.get("tseep") or {}).get("duration")
        return [t for t in wanted
                if not has_transient_frame(sol, t)
                and t > 0 and (dur is None or t <= float(dur))]

    def _store_transient_times(self, opts):
        """Write the dialog's times back into the model where the file keeps them.

        Stage times are model inputs shown at their point of use, so an edit always
        lands. The single-run stability time is different: the dialog choice governs
        this run only, and storing it is the explicit act that makes a scripted re-run
        reproduce it — hence the checkbox rather than an automatic write.
        """
        sd = self.doc.slope_data
        st = opts.get("stage_times")
        seep_time = opts.get("seep_time")
        want = {}
        if st and st.get("changed"):
            want["stage_1"] = st.get("stage_1")
            want["stage_2"] = st.get("stage_2")
        if seep_time and seep_time.get("write_back"):
            want["stability_time"] = seep_time.get("time")
        if not want or sd.get("tseep") is None:
            return
        if all(sd["tseep"].get(k) == v for k, v in want.items()):
            return
        self.doc.begin_edit("Transient times")
        sd["tseep"] = {**sd["tseep"], **want}
        self.doc.commit_edit()

    # --- FEM -------------------------------------------------------------
    def run_fem(self):
        if not self.doc.is_open or self._fem_runner is not None:
            return
        if self.doc.slope_data.get("mesh") is None:
            QMessageBox.information(self, "Run FEM", "Build a mesh first (Build Mesh…).")
            return
        material_names = [m.get("name", f"Material {i + 1}")
                          for i, m in enumerate(self.doc.slope_data.get("materials", []))]
        dlg = RunFemDialog(self, defaults=self._file_defaults("fem"),
                           material_names=material_names,
                           slope_data=self.doc.slope_data, document=self.doc,
                           transient=self.transient_solution(),
                           current_time=self.transient_current_time())
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_fem_opts = opts
        self._store_transient_times(opts)
        need = self._remarch_times_needed(opts)
        if need:
            self._start_remarch(need, ("fem", opts))
            return
        self._start_fem(opts)

    def _start_fem(self, opts):
        """Begin the FEM run itself, once the seepage frame it reads is in hand."""
        # SSRM supports cooperative cancel (a single-trial solve is quick).
        supports_cancel = opts["analysis"] == "ssrm"
        seep_time = opts.get("seep_time") or {}
        if not self._apply_transient_analysis_frame(rapid=False,
                                                    time=seep_time.get("time")):
            return
        self.statusBar().showMessage("Running FEM …")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self._fem_runner = FemRunner(self.doc.slope_data, opts, parent=self)
        self._fem_runner.succeeded.connect(self._on_fem_succeeded)
        self._fem_runner.failed.connect(self._on_fem_failed)
        self._fem_runner.cancelled.connect(self._on_fem_cancelled)
        self._fem_runner.progress.connect(self._on_run_progress)
        self._fem_runner.finished.connect(self._on_fem_finished)
        if supports_cancel:
            self.cancel_btn.setEnabled(True)
            self.cancel_btn.setVisible(True)
        self._update_run_actions()
        self._fem_runner.start()

    @staticmethod
    def _fem_meta(bundle):
        """The meta sidecar for one finite element bundle, written by both paths
        that save one — the run that just finished, and File → Save.

        The run's own record (:func:`xslope.fem.ssrm_run_record`, put on the
        bundle by the runner) is the base, and the three facts about the bundle
        are laid over it: they are what this writer knows and the run record does
        not name. Layered rather than merged blindly so a key can only ever mean
        what this function says it means.
        """
        meta = dict(bundle.get("meta") or {})
        meta["FS"] = bundle.get("FS")
        meta["analysis"] = bundle.get("analysis")
        # The strength-reduction factor the written field was solved at
        # (solution["F"]). NOT the factor of safety, which is a different number.
        meta["F"] = (bundle.get("solution") or {}).get("F")
        return meta

    def _on_fem_succeeded(self, bundle):
        self.doc.results["fem_solution"] = bundle
        if self.doc.path:
            try:
                from xslope.fem import export_fem_solution
                export_fem_solution(bundle["fem_data"], bundle["solution"],
                                    os.path.splitext(self.doc.path)[0],
                                    # What the run knows and the arrays do not.
                                    # The solve's OWN facts — whether it closed,
                                    # its iterations, its residual, how far the
                                    # section moved — are added by
                                    # export_fem_solution off the solution being
                                    # written, so a reload can state them instead
                                    # of assuming them.
                                    meta=self._fem_meta(bundle),
                                    failure_solution=bundle.get("failure_solution"))
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
                self.fem_data_canvas.render_fem_data(bundle["fem_data"], panel.options(),
                                                     style=self.doc.style or None)
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
            self.fem_details_btn = self.fem_results_canvas.add_tool_button(
                "1D Details…",
                "Per-line profiles for reinforcement and piles.",
                self.open_fem_details)
        self._rerender_fem_results()

    def _rerender_fem_results(self):
        bundle = self.doc.results.get("fem_solution")
        panel = self._display_panels.get(self.fem_results_canvas)
        if bundle and panel and self.fem_results_canvas is not None:
            try:
                self.fem_results_canvas.render_fem_results(
                    bundle["fem_data"], bundle["solution"],
                    {**panel.options(), "fs": bundle.get("FS"),
                     "failure_solution": bundle.get("failure_solution")})
            except Exception:
                traceback.print_exc()
        self._update_fem_details_action()

    def _update_fem_details_action(self):
        """Gate the FEM results view's "1D Details…" button.

        The button stays on the toolbar and dims with its reason when there is
        nothing to detail — a model with no reinforcement lines and no piles has
        no 1D elements, and the dialog would open empty."""
        if self.fem_details_btn is None:
            return
        from xslope.fem_details import has_1d_details
        bundle = self.doc.results.get("fem_solution")
        fem_data = (bundle or {}).get("fem_data")
        ok = bool(bundle) and has_1d_details(fem_data)
        self.fem_details_btn.setEnabled(ok)
        self.fem_details_btn.setToolTip(
            "Per-line profiles for reinforcement and piles."
            if ok else
            "This model has no reinforcement lines or piles to detail.")

    def open_fem_details(self):
        """Open the non-modal reinforcement/pile details dialog.

        The dialog holds the solution it was opened with, so a second press
        replaces any open one rather than raising a stale copy — the results
        may have been re-solved since."""
        bundle = self.doc.results.get("fem_solution")
        if not bundle:
            return
        from xslope.fem_details import has_1d_details
        if not has_1d_details(bundle.get("fem_data")):
            return
        if self.fem_details_dlg is not None:
            self.fem_details_dlg.close()
        from .fem_details_dialog import FemDetailsDialog
        self.fem_details_dlg = FemDetailsDialog(
            bundle["fem_data"], bundle["solution"], self.doc.slope_data,
            model_path=self.doc.path, parent=self,
            # The at-failure snapshot sits beside the solution on the bundle
            # (both the fresh-solve and the reload path put it there) — the same
            # field the results view's own Field state switch renders.
            failure_solution=bundle.get("failure_solution"))
        self.fem_details_dlg.show()
        self.fem_details_dlg.raise_()
        return self.fem_details_dlg

    # --- LEM analysis ----------------------------------------------------
    def run_lem(self):
        if not self.doc.is_open or self._runner is not None:
            return
        dlg = RunLemDialog(self, defaults=self._file_defaults("lem"),
                           slope_data=self.doc.slope_data, document=self.doc,
                           transient=self.transient_solution(),
                           current_time=self.transient_current_time())
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_lem_opts = opts
        self._store_surface_family(opts.get("surface"))
        self._store_transient_times(opts)
        need = self._remarch_times_needed(opts)
        if need:
            self._start_remarch(need, ("lem", opts))
            return
        self._start_lem(opts)

    def _start_lem(self, opts):
        """Begin the LEM run itself, once the seepage frame(s) it reads are in hand."""
        seep_time = opts.get("seep_time") or {}
        if not self._apply_transient_analysis_frame(rapid=opts.get("rapid", False),
                                                    time=seep_time.get("time")):
            return
        self.act_run.setEnabled(False)
        verb = {"auto_search": "Searching"}.get(opts["analysis"], "Running")
        self.statusBar().showMessage(f"{verb} {opts['method']} …")
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

    def _store_surface_family(self, surface):
        """Record which surface family this run chose, when the deck carries both.

        The file stores and the dialog edits: a model defining both a circular and a
        non-circular surface has no way, in its geometry, to say which one it means,
        and the circles simply win. The dialog's choice is written back into the
        model's ``surface_family`` — the ``main`` sheet's **Surface family** cell
        (D24), so it survives save and reload — and the legacy ``circular`` flag is
        kept in step with it, because that flag is what the runners, the plots and
        the sensitivity sweep read when nobody hands them a family.
        """
        sd = self.doc.slope_data
        if not (sd.get("circles") and sd.get("non_circ")):
            return                              # only one family: nothing to choose
        want = "noncircular" if surface == "noncircular" else "circular"
        if (sd.get("surface_family") == want
                and bool(sd.get("circular")) == (want == "circular")):
            return
        self.doc.begin_edit("Surface family")
        sd["surface_family"] = want
        sd["circular"] = want == "circular"
        self.doc.commit_edit()

    def _cancel_run(self):
        runner = next((r for r in (self._runner, self._fem_runner, self._sens_runner,
                                   self._rel_runner, self._seep_runner,
                                   self._report_runner, self._update_dl_runner)
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
        if isinstance(bundle.get("results"), dict):
            self._show_solution(bundle)
        # Lead with the most specific result view produced by this run.
        lead = (self.search_canvas if bundle.get("search") else self.solution_canvas)
        if lead is not None:
            self.view_tabs.setCurrentWidget(lead)
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

    # --- report ----------------------------------------------------------
    def report_solutions(self):
        """What the report can document, in :mod:`xslope.report`'s shape.

        Every engine the session has solved, keyed as the report reads them:
        ``lem`` per method, ``seep`` per boundary condition set, ``tseep`` for a
        transient march, ``fem`` for the stress or strength reduction run. The
        runners' own bundles are the shape the report takes, so they are passed
        through as they are — the seepage ones in boundary-condition order, so the
        section that documents BC 1 and BC 2 documents them in that order, and the
        march after them, being the run between the states they hold.

        The LEM bundle alone is carried with the method it was run under: the
        runner emits the solution, and the method the user chose is the run
        options', so the two are joined here rather than guessed from the result
        dict. The rest of those options ride along too — a report asked for a
        method this run did not solve runs it, and runs it the way this run was
        run rather than inferring it from the answer.
        """
        out = {}
        bundle = self.doc.results.get("lem_solution")
        if bundle:
            opts = self._last_lem_opts or {}
            method = opts.get("method")
            # The run's own record wins; the dialog's options are the fallback
            # for a bundle that predates the runner recording them.
            out["lem"] = [dict(bundle, method=method or bundle.get("method"),
                               options=bundle.get("options") or dict(opts))]
        seep = self.doc.results.get("seep_solutions") or {}
        if seep:
            out["seep"] = [seep[bc] for bc in sorted(seep)]
        transient = self.doc.results.get("transient_seep")
        if transient:
            out["tseep"] = [transient]
        fem = self.doc.results.get("fem_solution")
        if fem:
            out["fem"] = [fem]
        return out

    def generate_report(self):
        """File → Generate Report…: compose a report, build it off the GUI
        thread, and open it.

        The build is a run like any other — a dozen plots and then Word — so it
        gets a run's treatment: a worker thread, the progress bar counting the
        figures by name, a Cancel button, and a window that stays live
        throughout. The action is disabled while it builds, so a second report
        cannot be started over the first.
        """
        from .report_dialog import ReportDialog, finalization_enabled

        if self._report_runner is not None:
            return
        solutions = self.report_solutions()
        dlg = ReportDialog(self, slope_data=self.doc.slope_data,
                           solutions=solutions, model_path=self.doc.path,
                           default_method=(self._last_lem_opts or {}).get("method"),
                           settings=self.settings)
        if dlg.exec() != QDialog.Accepted:
            return
        path, fmt, options = dlg.output_path(), dlg.output_format(), dlg.options()
        options["style"] = self.doc.style

        self.statusBar().showMessage("Generating the report…")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(True)
        self.cancel_btn.setEnabled(True)
        self.cancel_btn.setVisible(True)
        self._report_runner = ReportRunner(
            self.doc.slope_data, solutions, options, path, fmt=fmt,
            finalize=finalization_enabled(self.settings), parent=self)
        self._report_runner.succeeded.connect(self._on_report_succeeded)
        self._report_runner.failed.connect(self._on_report_failed)
        self._report_runner.cancelled.connect(self._on_report_cancelled)
        self._report_runner.progress.connect(self._on_report_progress)
        self._report_runner.finished.connect(self._on_report_finished)
        self._update_run_actions()     # the report action dims while it builds
        self._report_runner.start()

    def _on_report_progress(self, done, total, label):
        """The run's progress line, plus what Word's stretch does to Cancel.

        Word is one call that returns when it is finished, so its stretch is
        indeterminate — and it is not interruptible: a report half-updated by
        Word is worse than one whose contents page the reader refreshes. Cancel
        goes out with the determinate phase.
        """
        if total <= 0:
            self.cancel_btn.setEnabled(False)
        self._on_run_progress(done, total, label)

    def _on_report_succeeded(self, out):
        from .report_dialog import open_output

        # The Word finish already ran on the worker thread; this only shows it.
        shown = open_output(out["path"], out.get("fmt"), finalize=False)
        self.statusBar().showMessage(
            f"Report written to {os.path.basename(out['path'])} "
            f"({len(out['figures'])} figures) — opening the {shown}.")

    def _on_report_failed(self, message):
        QMessageBox.warning(self, "Generate Report", message)
        self.statusBar().showMessage("The report could not be generated.")

    def _on_report_cancelled(self):
        self.statusBar().showMessage("Report cancelled.")

    def _on_report_finished(self):
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.cancel_btn.setVisible(False)
        self.cancel_btn.setEnabled(True)
        if self._report_runner is not None:
            self._report_runner.deleteLater()
            self._report_runner = None
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
                self.search_canvas.render_search(self.doc.slope_data, search,
                                                 panel.options(), style=self.doc.style or None)
            except Exception:
                traceback.print_exc()

    def _show_reliability(self, reliability_data):
        if self.reliability_canvas is None:
            self.reliability_canvas = MplCanvas(self)
            self.view_tabs.insertTab(1, self.reliability_canvas, "LEM · Reliability")
            panel = ReliabilityDisplayPanel()
            panel.changed.connect(self._rerender_reliability)
            self.display_stack.addWidget(panel)
            self._display_panels[self.reliability_canvas] = panel
        self._rerender_reliability()

    def _rerender_reliability(self):
        """Re-render the Taylor-series LEM reliability surface plot (fs_cache of the
        MLV / F± surfaces). Only the LEM Taylor engine produces those surfaces."""
        bundle = self.doc.results.get("reliability")
        rel = bundle.get("reliability") if bundle else None
        panel = self._display_panels.get(self.reliability_canvas)
        if (rel and panel and self.reliability_canvas is not None
                and rel.get("fs_cache")):
            try:
                self.reliability_canvas.render_reliability(
                    self.doc.slope_data, rel, panel.options(), style=self.doc.style or None)
            except Exception:
                traceback.print_exc()

    def _show_reliability_histogram(self):
        if self.reliability_hist_canvas is None:
            self.reliability_hist_canvas = MplCanvas(self)
            self.view_tabs.insertTab(1, self.reliability_hist_canvas,
                                     "Reliability · MC")
            panel = ReliabilityMcDisplayPanel()
            panel.changed.connect(self._rerender_reliability_histogram)
            self.display_stack.addWidget(panel)
            self._display_panels[self.reliability_hist_canvas] = panel
        # One tab serves both sampling engines, named for the one that filled it.
        bundle = self.doc.results.get("reliability") or {}
        idx = self.view_tabs.indexOf(self.reliability_hist_canvas)
        if idx >= 0:
            self.view_tabs.setTabText(
                idx, "Reliability · RS" if bundle.get("engine") == "rs"
                else "Reliability · MC")
        self._rerender_reliability_histogram()

    def _rerender_reliability_histogram(self):
        bundle = self.doc.results.get("reliability")
        rel = bundle.get("reliability") if bundle else None
        panel = self._display_panels.get(self.reliability_hist_canvas)
        if (rel and panel and self.reliability_hist_canvas is not None
                and rel.get("fs_samples") is not None):
            try:
                self.reliability_hist_canvas.render_reliability_histogram(
                    rel, panel.options())
            except Exception:
                traceback.print_exc()

    # --- reliability run -------------------------------------------------
    def run_reliability(self):
        """Open the Reliability dialog and start a probabilistic run (the sibling of
        the Parametric study). LEM offers Taylor series / Monte Carlo; FEM offers the
        Taylor series over SSRM solves."""
        if not self.doc.is_open or self._rel_runner is not None:
            return
        if self._mode not in ("lem", "fem"):
            return
        if self._mode == "fem" and self.doc.slope_data.get("mesh") is None:
            QMessageBox.warning(self, "No mesh",
                                "Build a finite-element mesh first (Build Mesh…) — "
                                "FEM reliability runs on the mesh.")
            return
        dlg = ReliabilityDialog(self, defaults=self._last_rel_opts.get(self._mode, {}),
                                slope_data=self.doc.slope_data, app_mode=self._mode,
                                document=self.doc)
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_rel_opts[self._mode] = opts
        self.act_reliability.setEnabled(False)
        self.act_run.setEnabled(False)
        eng = {"taylor": "Taylor series", "mc": "Monte Carlo",
               "rs": "Response surface"}.get(opts["engine"], opts["engine"])
        self.statusBar().showMessage(f"Reliability ({eng}) …")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self.cancel_btn.setEnabled(True)
        self.cancel_btn.setVisible(True)
        self._rel_runner = ReliabilityRunner(self.doc.slope_data, opts, parent=self)
        self._rel_runner.succeeded.connect(self._on_rel_succeeded)
        self._rel_runner.failed.connect(self._on_rel_failed)
        self._rel_runner.cancelled.connect(self._on_rel_cancelled)
        self._rel_runner.progress.connect(self._on_run_progress)
        self._rel_runner.finished.connect(self._on_rel_finished)
        self._update_run_actions()
        self._rel_runner.start()

    def _on_rel_succeeded(self, bundle):
        self.doc.results["reliability"] = bundle
        rel = bundle["reliability"]
        engine = bundle.get("engine")
        app_mode = bundle.get("app_mode")
        if engine in ("mc", "rs"):
            # Both sampling engines report an FS distribution; the response
            # surface's histogram is drawn from a fixed-stride subsample of its
            # ten million realizations.
            self._show_reliability_histogram()
            lead = self.reliability_hist_canvas
        elif app_mode == "fem":
            # Reuse the FEM Data / Results views for the deformation at the mean values.
            self.doc.results["fem_solution"] = {
                "fem_data": bundle["fem_data"], "solution": bundle["solution"],
                "FS": bundle["FS"], "analysis": "reliability", "reliability": rel}
            self._show_fem_data(bundle["fem_data"])
            self._show_fem_results()
            lead = self.fem_results_canvas
        else:  # LEM Taylor — the MLV / F± surface plot
            self._show_reliability(rel)
            lead = self.reliability_canvas
        if lead is not None:
            self.view_tabs.setCurrentWidget(lead)
        self.statusBar().showMessage(self._reliability_status(rel, engine))
        QMessageBox.information(self, "Reliability results",
                                self._reliability_summary(rel, engine, app_mode))

    def _on_rel_failed(self, message):
        QMessageBox.warning(self, "Reliability run failed", message)
        self.statusBar().showMessage("Reliability run failed.")

    def _on_rel_cancelled(self):
        self.statusBar().showMessage("Reliability run cancelled.")

    def _on_rel_finished(self):
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.cancel_btn.setVisible(False)
        if self._rel_runner is not None:
            self._rel_runner.deleteLater()
            self._rel_runner = None
        self._update_run_actions()

    @staticmethod
    def _reliability_status(rel, engine):
        if engine in ("mc", "rs"):
            return (f"Reliability ({'MC' if engine == 'mc' else 'RS'}) — mean FS = "
                    f"{rel['mean_FS']:.3f}, β_ln = {rel['beta_ln']:.3f}, "
                    f"Pf = {rel['pf_empirical'] * 100:.2f}%")
        return (f"Reliability — F_MLV = {rel['F_MLV']:.3f}, "
                f"β_ln = {rel['beta_ln']:.3f}, "
                f"Pf = {rel['prob_failure'] * 100:.2f}%")

    @staticmethod
    def _reliability_summary(rel, engine, app_mode):
        if engine == "rs":
            return (
                f"Response surface ({rel.get('distribution', 'normal')}, "
                f"{rel.get('n_surrogate', 0):,} surrogate realizations)\n\n"
                f"mean FS = {rel['mean_FS']:.3f}\n"
                f"σ_F = {rel['sigma_F']:.3f}    COV_F = {rel['COV_F']:.3f}\n"
                f"β (normal) = {rel['beta_normal']:.3f}\n"
                f"β (lognormal) = {rel['beta_ln']:.3f}\n"
                f"Pf (empirical) = {rel['pf_empirical'] * 100:.3f}%\n"
                f"Pf (normal) = {rel['pf_normal'] * 100:.3f}%    "
                f"Pf (lognormal) = {rel['pf_lognormal'] * 100:.3f}%\n\n"
                f"Surrogate: {rel['n_design']} design solves, checked against "
                f"{rel['n_gate_valid']} held-out solves — R² = {rel['gate_r2']:.5f}, "
                f"RMS error {rel['gate_rms']:.4f} of a spread of "
                f"{rel['gate_sigma']:.3f}, {rel['gate_pf_disagree_n']} draw(s) on "
                f"the wrong side of F = 1.\n\n"
                "The FS histogram is on the Reliability · RS tab; the fit and gate "
                "summary is in the Log pane.")
        if engine == "mc":
            return (
                f"Monte Carlo ({rel.get('distribution', 'normal')}, "
                f"{rel.get('n_valid')} of {rel.get('n_samples')} valid)\n\n"
                f"mean FS = {rel['mean_FS']:.3f}\n"
                f"σ_F = {rel['sigma_F']:.3f}    COV_F = {rel['COV_F']:.3f}\n"
                f"β (normal) = {rel['beta_normal']:.3f}\n"
                f"β (lognormal) = {rel['beta_ln']:.3f}\n"
                f"Pf (empirical) = {rel['pf_empirical'] * 100:.2f}%\n"
                f"Pf (normal) = {rel['pf_normal'] * 100:.2f}%    "
                f"Pf (lognormal) = {rel['pf_lognormal'] * 100:.2f}%\n\n"
                "The FS histogram is on the Reliability · MC tab; the per-parameter "
                "table is in the Log pane.")
        tag = "FEM (SSRM)" if app_mode == "fem" else "Taylor series (TSPM)"
        extra = ("\n\nThe FEM Results view shows the deformation at the most-likely "
                 "values." if app_mode == "fem" else "")
        return (
            f"{tag}\n\n"
            f"F_MLV = {rel['F_MLV']:.3f}\n"
            f"σ_F = {rel['sigma_F']:.3f}    COV_F = {rel['COV_F']:.3f}\n"
            f"β (lognormal) = {rel['beta_ln']:.3f}\n"
            f"Reliability = {rel['reliability'] * 100:.2f}%\n"
            f"Probability of failure = {rel['prob_failure'] * 100:.2f}%"
            + extra + "\n\nThe per-parameter ΔF table is in the Log pane.")

    def _show_solution(self, bundle):
        if self.solution_canvas is None:
            # A SolutionView (warning strip + canvas) rather than a bare MplCanvas,
            # so non-empty admissibility warnings surface above the plot. It forwards
            # render_solution/ensure_fitted, so the tab machinery is unchanged.
            self.solution_canvas = SolutionView(self)
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
                    bundle["failure_surface"], bundle["results"], panel.options(),
                    style=self.doc.style or None)
            except Exception:
                traceback.print_exc()

    # --- sensitivity / design study --------------------------------------
    def run_sensitivity(self):
        if not self.doc.is_open or self._sens_runner is not None:
            return
        # FEM/Seep sweeps run on the mesh — require one, like Run.
        if self._mode != "lem" and self.doc.slope_data.get("mesh") is None:
            QMessageBox.warning(self, "No mesh",
                                "Build a finite-element mesh first (Build Mesh…) — "
                                "FEM and seepage sweeps run on the mesh.")
            return
        dlg = SensitivityDialog(self, defaults=self._last_sens_opts.get(self._mode, {}),
                                slope_data=self.doc.slope_data, app_mode=self._mode,
                                document=self.doc,
                                transient=self.transient_solution())
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_sens_opts[self._mode] = opts
        study = opts["mode"]
        if study in ("design", "back_analysis") and not opts.get("param"):
            QMessageBox.warning(self, "Nothing to sweep",
                                "Pick a material and property to sweep.")
            return
        # The variance / MC-rank plots use every σ-carrying material, so an empty
        # table is fine for them; the tornado / scaled / spider plots need a table.
        if (study == "sensitivity" and not opts.get("params")
                and opts.get("plot_type", "tornado") in ("tornado", "scaled", "spider")):
            QMessageBox.warning(self, "Nothing to sweep",
                                "Add at least one parameter to the table.")
            return
        self.act_sensitivity.setEnabled(False)
        self.act_run.setEnabled(False)
        verb = {"design": "Design sweep", "back_analysis": "Back-analysis",
                "fs_vs_time": ("Rapid drawdown vs time" if opts.get("rapid")
                               else "FS vs time")}.get(study, "Sensitivity sweep")
        self.statusBar().showMessage(f"{verb} — {opts['method']} …")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self.cancel_btn.setEnabled(True)
        self.cancel_btn.setVisible(True)
        self._sens_runner = SensitivityRunner(self.doc.slope_data, opts, parent=self,
                                             transient=self.transient_solution())
        self._sens_runner.succeeded.connect(self._on_sens_succeeded)
        self._sens_runner.failed.connect(self._on_sens_failed)
        self._sens_runner.cancelled.connect(self._on_sens_cancelled)
        self._sens_runner.progress.connect(self._on_run_progress)
        self._sens_runner.finished.connect(self._on_sens_finished)
        self._sens_runner.start()

    def _on_sens_succeeded(self, bundle):
        if bundle.get("kind") == "fs_vs_time":
            self.doc.results["fs_vs_time"] = bundle
            self._show_fs_vs_time()
            if self.fs_time_canvas is not None:
                self.view_tabs.setCurrentWidget(self.fs_time_canvas)
            unit = str(self.doc.slope_data.get("time_unit") or "").strip()
            what = ("Rapid drawdown vs time" if bundle.get("rapid")
                    else "FS vs time")
            if bundle.get("min_fs") is not None:
                self.statusBar().showMessage(
                    f"{what} — lowest {bundle['min_fs']:.3f} at t = "
                    f"{bundle['critical_time']:g}{(' ' + unit) if unit else ''} "
                    f"over {len(bundle.get('times', []))} instant(s)")
            else:
                self.statusBar().showMessage(
                    f"{what} — no instant produced a result; see the Log pane.")
        elif bundle.get("kind") == "design":
            self.doc.results["design"] = bundle
            self._show_design()
            if self.design_canvas is not None:
                self.view_tabs.setCurrentWidget(self.design_canvas)
            verb = ("Back-analysis" if bundle.get("study") == "back_analysis"
                    else "Design")
            if bundle.get("bracketed"):
                self.statusBar().showMessage(
                    f"{verb} — {bundle.get('output', 'FS')} = "
                    f"{bundle['target_fs']:g} at "
                    f"{bundle['param'].split(':')[-1]} = {bundle['crossing']:.4g}")
            else:
                self.statusBar().showMessage(bundle.get("message", f"{verb} done."))
        else:
            self.doc.results["sensitivity"] = bundle
            self._show_sensitivity()
            if self.sens_canvas is not None:
                self.view_tabs.setCurrentWidget(self.sens_canvas)
            pt = bundle.get("plot_type", "tornado")
            msg = {
                "tornado": (f"Sensitivity — {len(bundle.get('sweeps', {}))} "
                            f"parameter(s); double-click a tornado bar for its curve."),
                "scaled": (f"Scaled sensitivity — "
                           f"{len(bundle.get('scaled', {}).get('bars', []))} "
                           f"parameter(s), {bundle.get('scaling', 'elasticity')}."),
                "spider": (f"Spider — {len(bundle.get('sweeps', {}))} parameter(s)."),
                "variance": (f"Variance contribution — "
                             f"{len(bundle.get('variance', {}).get('bars', []))} "
                             f"uncertain parameter(s)."),
                "rank": (f"Monte Carlo rank correlation — "
                         f"{bundle.get('rank', {}).get('n_valid')} valid samples."),
            }.get(pt, "Sensitivity done.")
            self.statusBar().showMessage(msg)

    def _on_sens_failed(self, message):
        QMessageBox.warning(self, "Sweep failed", message)
        self.statusBar().showMessage("Sweep failed.")

    def _on_sens_cancelled(self):
        self.statusBar().showMessage("Sweep cancelled.")

    def _on_sens_finished(self):
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.cancel_btn.setVisible(False)
        if self._sens_runner is not None:
            self._sens_runner.deleteLater()
            self._sens_runner = None
        self._update_run_actions()

    def _show_fs_vs_time(self):
        if self.fs_time_canvas is None:
            self.fs_time_canvas = SweepCanvas(self)
            self.view_tabs.addTab(self.fs_time_canvas, "FS vs Time")
        # The tab says which curve it is holding: a drawdown curve and a
        # single-stage curve are different analyses of the same march, and the tab
        # is reused for whichever ran last.
        i = self.view_tabs.indexOf(self.fs_time_canvas)
        if i >= 0:
            self.view_tabs.setTabText(
                i, "Drawdown vs Time"
                if (self.doc.results.get("fs_vs_time") or {}).get("rapid")
                else "FS vs Time")
        self._rerender_fs_vs_time()

    def _rerender_fs_vs_time(self):
        bundle = self.doc.results.get("fs_vs_time")
        if bundle and self.fs_time_canvas is not None:
            try:
                self.fs_time_canvas.render_fs_vs_time(bundle,
                                                      slope_data=self.doc.slope_data)
            except Exception:
                traceback.print_exc()

    def _show_design(self):
        if self.design_canvas is None:
            self.design_canvas = SweepCanvas(self)
            self.view_tabs.addTab(self.design_canvas, "Design")
            # No display-option panel: the curve/target are set at run time. The
            # Display dock tracks the tab (shows the placeholder), like other views.
        self._rerender_design()

    def _rerender_design(self):
        bundle = self.doc.results.get("design")
        if bundle and self.design_canvas is not None:
            try:
                self.design_canvas.render_design(bundle["df"], bundle["target_fs"],
                                                 bundle)
            except Exception:
                traceback.print_exc()

    def _show_sensitivity(self):
        # A fresh sweep invalidates any prior click-through curve (different data),
        # so drop the curve tab; the user re-clicks a bar on the new tornado.
        if self.sens_curve_canvas is not None:
            idx = self.view_tabs.indexOf(self.sens_curve_canvas)
            if idx >= 0:
                self.view_tabs.removeTab(idx)
            self.sens_curve_canvas.deleteLater()
            self.sens_curve_canvas = None
        if self.sens_canvas is None:
            self.sens_canvas = SweepCanvas(self)
            self.view_tabs.addTab(self.sens_canvas, "Sensitivity")
            # Double-click a bar to open that parameter's FS-vs-value curve (tornado
            # only; enabled per-render below).
            self.sens_canvas.picked.connect(self._on_tornado_pick)
        # Click-through is a tornado affordance; disable it for the other plots.
        bundle = self.doc.results.get("sensitivity") or {}
        tornado_like = bundle.get("plot_type", "tornado") == "tornado"
        self.sens_canvas.set_pick_enabled(tornado_like)
        self.sens_canvas._hint_label.setText(
            "(double-click a bar to see its FS curve)" if tornado_like else "")
        self._rerender_sensitivity()

    def _rerender_sensitivity(self):
        bundle = self.doc.results.get("sensitivity")
        if not bundle or self.sens_canvas is None:
            return
        pt = bundle.get("plot_type", "tornado")
        try:
            if pt == "scaled":
                self.sens_canvas.render_scaled(
                    bundle["scaled"], scaling=bundle.get("scaling", "elasticity"))
            elif pt == "spider":
                self.sens_canvas.render_spider(bundle["sweeps"])
            elif pt == "variance":
                self.sens_canvas.render_variance(bundle["variance"])
            elif pt == "rank":
                self.sens_canvas.render_rank(bundle["rank"])
            else:
                self.sens_canvas.render_tornado(bundle["tornado"])
        except Exception:
            traceback.print_exc()

    def _on_tornado_pick(self, x, y, _tol):
        """Map a double-clicked tornado bar (its y row) to a parameter and show
        that parameter's FS-vs-value curve in a companion tab."""
        bundle = self.doc.results.get("sensitivity")
        if not bundle or self.sens_canvas is None:
            return
        ax = self.sens_canvas._main_axes()
        if ax is None:
            return
        labels = [t.get_text() for t in ax.get_yticklabels()]
        k = int(round(y))
        if not (0 <= k < len(labels)):
            return
        param = labels[k]
        df = bundle["sweeps"].get(param)
        if df is None:
            return
        if self.sens_curve_canvas is None:
            self.sens_curve_canvas = SweepCanvas(self)
            self.view_tabs.addTab(self.sens_curve_canvas, "Sensitivity · Curve")
        self.sens_curve_canvas.render_curve(df, target_fs=bundle.get("target_fs"))
        self.view_tabs.setCurrentWidget(self.sens_curve_canvas)

    def _clear_result_tabs(self):
        """Drop result views (e.g. on opening another file) so they don't show
        stale results from the previous project."""
        single = ("mesh_canvas", "search_canvas", "solution_canvas",
                  "reliability_canvas", "reliability_hist_canvas", "sens_canvas",
                  "sens_curve_canvas", "design_canvas", "fem_data_canvas",
                  "fem_results_canvas")
        # The seep canvases are per-BC dicts; flatten them in with the rest, along
        # with the transient view (frames + play bar) when present.
        canvases = [getattr(self, a) for a in single]
        canvases += list(self.seep_data_canvas.values())
        canvases += list(self.seep_solution_canvas.values())
        if self.transient_seep_view is not None:
            canvases.append(self.transient_seep_view)
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
        self.transient_seep_view = None
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
        # Heads-up: saving in place onto an older-format file that now uses a
        # v11-only feature (e.g. a van Genuchten material) will upgrade the file to
        # the current template format. Use the engine's own predicate so the dialog
        # appears exactly when the upgrade will happen.
        try:
            from xslope.fileio import _inplace_save_would_drop
            if _inplace_save_would_drop(self.doc.path,
                                        self.doc.slope_data.get("materials", [])):
                QMessageBox.information(
                    self, "File will be upgraded",
                    f"“{os.path.basename(self.doc.path)}” was created in an "
                    "older template version that lacks columns this model now uses "
                    "(e.g. van Genuchten unsaturated parameters). Saving will upgrade "
                    "it to the current template format.")
        except Exception:
            pass
        try:
            self.doc.save(template=None)   # edit in place, preserve formatting
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(self, "Save failed", str(exc))
            return
        self._sync_sidecars(os.path.splitext(self.doc.path)[0])
        self.statusBar().showMessage(f"Saved {os.path.basename(self.doc.path)}")

    def _sync_sidecars(self, stem):
        """Make the mesh / seep / FEM sidecars next to the saved .xlsx match the
        in-memory project. (Re)write those whose artifact is present (so Save As
        carries them to the new name), and DELETE those whose artifact was
        invalidated — e.g. a geometry edit cleared the mesh and solutions. Without
        this, a stale ``{stem}_mesh.json`` / ``_seep.csv`` / ``_fem_*.csv`` would be
        auto-loaded on the next Open and silently mismatch the edited inputs.
        Best-effort: a failure on one sidecar is logged, not fatal."""
        sd = self.doc.slope_data
        results = self.doc.results

        def remove(path):
            try:
                if os.path.exists(path):
                    os.remove(path)
            except Exception:
                traceback.print_exc()

        # Mesh ({stem}_mesh.json), auto-loaded by load_slope_data.
        mesh = sd.get("mesh")
        if mesh is not None:
            try:
                from xslope.mesh import export_mesh_to_json
                export_mesh_to_json(mesh, f"{stem}_mesh.json")
            except Exception:
                traceback.print_exc()
        else:
            remove(f"{stem}_mesh.json")

        # Seepage solutions, per BC set ({stem}_seep.csv / _seep2.csv).
        seep = results.get("seep_solutions", {})
        for bc, suffix, key in ((1, "_seep.csv", "seep_u"), (2, "_seep2.csv", "seep_u2")):
            path = stem + suffix
            bundle = seep.get(bc)
            imported = sd.get(key)
            if bundle:
                try:
                    from xslope.seep import export_seep_solution
                    export_seep_solution(bundle["seep_data"], bundle["solution"], path)
                except Exception:
                    traceback.print_exc()
            elif mesh is not None and imported is not None and len(imported):
                # A pore-pressure field xslope did not solve for -- lifted out of a solved
                # SEEP/W analysis by the GeoStudio importer. There is no solver bundle
                # behind it, and deleting the file on that basis would silently strip the
                # water out of the model: it would reload dry, with every material still
                # asking for a seepage solution that no longer existed.
                try:
                    from xslope.seep import export_seep_u
                    from xslope.units import require_gamma_water
                    export_seep_u(mesh["nodes"], imported, path,
                                  require_gamma_water(sd, "seep-field re-export"))
                except Exception:
                    traceback.print_exc()
            else:
                remove(path)

        # Transient seepage frames ({stem}_tseep.csv + _tseep_meta.json). Written when
        # a transient run is present (so Save As carries it), removed otherwise.
        tview = results.get("transient_seep")
        if tview:
            try:
                from xslope.seep import export_transient_solution
                export_transient_solution(
                    tview["seep_data"], tview["transient"], stem,
                    input_file=self.doc.path,
                    mesh_file=f"{os.path.basename(stem)}_mesh.json")
            except Exception:
                traceback.print_exc()
        else:
            remove(f"{stem}_tseep.csv")
            remove(f"{stem}_tseep_meta.json")

        # FEM solution ({stem}_fem_nodes.csv / _fem_elements.csv / _fem_meta.json).
        fem = results.get("fem_solution")
        if fem:
            try:
                from xslope.fem import export_fem_solution
                export_fem_solution(fem["fem_data"], fem["solution"], stem,
                                    meta=self._fem_meta(fem),
                                    failure_solution=fem.get("failure_solution"))
            except Exception:
                traceback.print_exc()
        else:
            for suffix in FEM_SOLUTION_SIDECARS:
                remove(f"{stem}{suffix}")

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
        self._sync_sidecars(os.path.splitext(path)[0])
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
        from xslope import __version__
        QMessageBox.about(
            self, f"About {APP_NAME}",
            f"<b>{APP_NAME}</b> {__version__}<br>Desktop GUI for the xslope "
            f"slope-stability engine.<br><br>Open an Excel input file to view its "
            f"geometry and inputs.")

    # --- documentation and welcome ---------------------------------------
    def open_documentation(self):
        """Help → Documentation — the documentation root, in the browser."""
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        return QDesktopServices.openUrl(QUrl(links.DOCS_URL))

    def show_welcome(self):
        """Help → Welcome, and the first launch (opened from ``studio.app``).

        Modal: it is a page of text with links on it, and the window behind it has
        nothing on it yet the first time this runs.
        """
        dlg = WelcomeDialog(self, settings=self.settings)
        return dlg.exec()

    # --- updates ---------------------------------------------------------
    def _check_updates(self):
        """Help → Check for Updates… — always answers with one dialog."""
        self.updates.check_manual()

    def check_updates_at_startup(self):
        """The silent once-a-day check. Called by app.py after the window is up
        (never from ``__init__``, so tests and embedded uses stay offline)."""
        return self.updates.check_at_startup()

    def stop_threads(self):
        """Bring every background thread down, and wait for it.

        Called from :meth:`closeEvent` and from the application's ``aboutToQuit``,
        because closing the window is not the only way the app ends. On macOS,
        Quit (Cmd-Q) ends ``exec()`` without ever sending the window a close
        event, and the window is then destroyed by PySide's own teardown at
        interpreter shutdown — with the mesh thread still in its event loop.
        ``QThread``'s destructor calls ``qFatal`` on a thread that is still
        running, so that exit aborts the process: the "python quit unexpectedly"
        report a user gets after quitting normally.

        Idempotent, so the two paths can both run it (a window closed with Cmd-W
        and then quit, say) and the second call finds nothing left to stop.
        """
        if self._runner is not None and self._runner.isRunning():
            self._runner.cancel()     # ask an in-flight run to stop, then wait briefly
            self._runner.wait(5000)
        if self._seep_runner is not None and self._seep_runner.isRunning():
            self._seep_runner.wait(10000)   # seepage has no cancel hook; let it finish
        if self._fem_runner is not None and self._fem_runner.isRunning():
            self._fem_runner.cancel()       # SSRM stops cooperatively
            self._fem_runner.wait(15000)
        if self._sens_runner is not None and self._sens_runner.isRunning():
            self._sens_runner.cancel()      # sweep stops at the next point
            self._sens_runner.wait(15000)
        if self._rel_runner is not None and self._rel_runner.isRunning():
            self._rel_runner.cancel()       # reliability stops at the next sample
            self._rel_runner.wait(15000)
        if self._report_runner is not None and self._report_runner.isRunning():
            self._report_runner.cancel()    # stops at the next figure
            # Long enough to outlast a Word finish, which is never cut short: the
            # document Word is holding is the user's report.
            self._report_runner.wait(90000)
        if self._update_dl_runner is not None and self._update_dl_runner.isRunning():
            self._update_dl_runner.cancel()  # abandon a part-downloaded installer
            self._update_dl_runner.wait(5000)
        # Stop the persistent mesh thread (lets an in-flight build finish first).
        if self._mesh_thread is not None and self._mesh_thread.isRunning():
            self._mesh_thread.quit()
            self._mesh_thread.wait(10000)

    def closeEvent(self, event):
        self.stop_threads()
        if self.doc.is_open and self.doc.dirty:
            box = QMessageBox(self)
            box.setIcon(QMessageBox.Question)
            box.setWindowTitle("Unsaved changes")
            box.setText("Save changes before closing?")
            save_btn = box.addButton(QMessageBox.Save)
            discard_btn = box.addButton(QMessageBox.Discard)
            cancel_btn = box.addButton(QMessageBox.Cancel)
            box.setDefaultButton(save_btn)
            box.exec()
            # Compare the actual clicked button by identity — the StandardButton
            # return value is unreliable for Discard ("Don't Save") on macOS.
            clicked = box.clickedButton()
            if clicked is cancel_btn or clicked is None:
                event.ignore()
                return
            if clicked is save_btn:
                self.save()
                if self.doc.dirty:        # save failed or was cancelled
                    event.ignore()
                    return
            # discard_btn → fall through and close without saving
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        super().closeEvent(event)
