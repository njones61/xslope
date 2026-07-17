"""MainWindow — the XSlope Studio shell (Phase 1: read-only viewer).

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

from PySide6.QtCore import Qt, QObject, QSettings, QThread, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QComboBox, QDockWidget, QFileDialog, QHBoxLayout, QLabel, QMainWindow,
    QMenu, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QStackedWidget,
    QTabWidget, QToolBar, QToolButton, QTreeWidget, QTreeWidgetItem,
    QVBoxLayout, QWidget,
)

from xslope.fileio import default_template_path

from .canvas import MplCanvas
from .dialogs import (
    BuildMeshDialog, DxfImportDialog, GszImportDialog, RunFemDialog, RunLemDialog,
    RunSeepDialog, SensitivityDialog, Slide2ImportDialog,
)
from .display_panels import (
    FeDataDisplayPanel, FemResultsDisplayPanel, InputsDisplayPanel,
    MeshDisplayPanel, ReliabilityDisplayPanel, SearchDisplayPanel,
    SeepDisplayPanel, SolutionDisplayPanel,
)
from .document import ProjectDocument
from .editors import CATEGORY_EDITORS
from .runners import FemRunner, LemRunner, MeshWorker, SeepRunner, SensitivityRunner

APP_NAME = "XSlope Studio"
ORG_NAME = "XSlope"
MAX_RECENT = 8
MODES = [("LEM", "lem"), ("Seepage", "seep"), ("FEM", "fem")]
# Blank template used to create files on Save As — the single copy bundled with
# the engine package (xslope/resources), so the GUI and library share one source.
TEMPLATE = default_template_path()
CATEGORY_ROLE = Qt.UserRole + 1


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

    def render_curve(self, df, target_fs=None):
        from xslope.plot import plot_sensitivity
        self._draw(lambda fig: plot_sensitivity(df, target_fs=target_fs, fig=fig),
                   dxf=False)

    def render_design(self, df, target_fs, summary):
        from xslope.plot import plot_sensitivity

        def draw(fig):
            plot_sensitivity(df, target_fs=target_fs, fig=fig)
            self._annotate_crossing(fig, target_fs, summary)

        self._draw(draw, dxf=False)

    def _annotate_crossing(self, fig, target_fs, summary):
        """Mark the interpolated output=target crossing (or an honest miss note)."""
        ax = self._main_axes()
        if ax is None:
            return
        param = summary.get("param", "")
        short = param.split(":")[-1] or param
        out = summary.get("output", "FS")            # 'FS' or 'q'
        if summary.get("bracketed") and summary.get("crossing") is not None:
            xc = summary["crossing"]
            ax.axvline(xc, color="#0a7d2c", linestyle="--", linewidth=1.0)
            ax.plot([xc], [target_fs], marker="D", color="#0a7d2c", ms=9, zorder=8)
            ax.annotate(f"{short} = {xc:.4g}\nfor {out} = {target_fs:g}",
                        xy=(xc, target_fs), xytext=(8, 14),
                        textcoords="offset points", color="#0a7d2c", fontsize=9,
                        fontweight="bold", zorder=9,
                        bbox=dict(boxstyle="round,pad=0.3", fc="white",
                                  ec="#0a7d2c", alpha=0.9))
        else:
            # Put the note in the empty band: just under the target line when the
            # curve sits below the target (need a higher FS), else near the bottom.
            fs_min = (summary.get("fs_range") or (None, None))[0]
            below_range = fs_min is not None and target_fs < fs_min
            y, va = (0.05, "bottom") if below_range else (0.95, "top")
            ax.text(0.5, y, summary.get("message", "Target FS not reached."),
                    transform=ax.transAxes, ha="center", va=va, fontsize=9,
                    color="#7a5200", wrap=True, zorder=9,
                    bbox=dict(boxstyle="round,pad=0.4", fc="#fff4d6",
                              ec="#e0b400", alpha=0.95))


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
        self.settings = QSettings(ORG_NAME, APP_NAME)

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
        # Sensitivity / design study result tabs.
        self.sens_canvas = None            # tornado
        self.sens_curve_canvas = None      # click-through FS-vs-value curve
        self.design_canvas = None          # design curve + target crossing
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
        self._sens_runner = None
        self._mesh_busy = False
        self._run_implemented = {"lem", "seep", "fem"}   # modes whose Run is wired up
        self._last_lem_opts = {}
        self._last_sens_opts = {}          # keyed by engine mode (lem/fem/seep)
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
        # up/down arrows sit centred next to the button, and vertically centre them.
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
        self.act_quit = QAction("&Quit", self, shortcut=QKeySequence.Quit,
                                triggered=self.close)
        self.act_undo = QAction("&Undo", self, shortcut=QKeySequence.Undo,
                                triggered=self._undo, enabled=False)
        self.act_redo = QAction("&Redo", self, shortcut=QKeySequence.Redo,
                                triggered=self._redo, enabled=False)
        self.act_about = QAction("&About", self, triggered=self._about)
        self.act_save = QAction("&Save", self, shortcut=QKeySequence.Save,
                                enabled=False, triggered=self.save)
        self.act_save_as = QAction("Save &As…", self, enabled=False, triggered=self.save_as)
        self.act_run = QAction("Run &LEM…", self, enabled=False, triggered=self.run_current)
        self.act_sensitivity = QAction("Sensitivity / &Design…", self, enabled=False,
                                       triggered=self.run_sensitivity)
        self.act_build_mesh = QAction("Build &Mesh…", self, enabled=False,
                                      triggered=self.build_mesh)

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

        m_view = mb.addMenu("&View")
        m_view.addAction(self.inputs_dock.toggleViewAction())
        m_view.addAction(self.display_dock.toggleViewAction())
        m_view.addAction(self.log_dock.toggleViewAction())
        m_view.addAction(self.chat_dock.toggleViewAction())

        m_help = mb.addMenu("&Help")
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
        self.mode_combo = QComboBox()
        for label, _ in MODES:
            self.mode_combo.addItem(label)
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        tb.addWidget(self.mode_combo)
        tb.addSeparator()
        tb.addAction(self.act_build_mesh)
        tb.addAction(self.act_run)
        # Sensitivity / Design lives on the toolbar too (Norm's ask); the action's
        # existing mode-visibility (LEM-only) hides the button in other modes.
        tb.addAction(self.act_sensitivity)
        # macOS's native style draws text-only toolbar buttons in the larger system
        # font and ignores setFont; a stylesheet forces the size so New/Open/Run LEM
        # match the "Mode:" label. pointSizeF() is -1 for pixel-defined fonts.
        pt = self.mode_label.font().pointSizeF()
        if pt > 0:
            tb.setStyleSheet(f"QToolButton {{ font-size: {pt:g}pt; }}")

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
                "pip install \"xslope[gui]\"), then restart XSlope Studio.")
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
                "    pip install ezdxf\n\nthen restart XSlope Studio.")
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
        self.styles_btn.setEnabled(True)
        self.assistant.reset()        # new project -> fresh conversation
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
            # Restore the strength-reduction factor the result plots show in their
            # subplot titles (solution["F"]); fall back to the SSRM FS.
            F_saved = meta.get("F")
            if F_saved is None:
                F_saved = meta.get("FS")
            if F_saved is not None:
                solution["F"] = F_saved
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
                or self._mesh_busy)
        has_mesh = open_ and self.doc.slope_data.get("mesh") is not None
        # Sensitivity / design has a version for every mode (LEM: FS; FEM: FS via
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
        # An empty (new) project defaults to profile-based so the user can add the
        # first profile line; polygons are then derived from it.
        profile_based = bool(profile_lines) or not polygons
        add("Profile lines", len(profile_lines),
            category="profile" if profile_based else None)
        # Polygons are derived from profile lines for profile-based files (edit them
        # via the profile editor); only polygon-based files edit polygons directly.
        add("Polygons", len(polygons),
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
        try:
            return json.dumps({k: clean(sd.get(k)) for k in MainWindow.MESH_KEYS},
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
                  "fem_data_canvas", "fem_results_canvas"]
        if clear_mesh:
            single.append("mesh_canvas")
        canvases = [getattr(self, a) for a in single]
        canvases += list(self.seep_data_canvas.values())
        canvases += list(self.seep_solution_canvas.values())
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
        for key in ("lem_solution", "seep_solutions", "fem_solution",
                    "design", "sensitivity"):
            self.doc.results.pop(key, None)
        if clear_mesh:
            self.doc.slope_data["mesh"] = None
            self.doc.results.pop("mesh", None)
            self._update_run_actions()       # Seep/FEM Run re-gated on a built mesh
        if removed:
            self._show_display_for_tab(self.view_tabs.currentWidget())
            what = "results and mesh" if clear_mesh else "result(s)"
            self.statusBar().showMessage(f"Inputs changed — cleared the now-stale {what}.")

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
        for name in ("_seep.csv", "_seep2.csv", "_fem_nodes.csv",
                     "_fem_elements.csv", "_fem_meta.json"):
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
                canvas.render_seep_solution(
                    bundle["seep_data"], bundle["solution"], panel.options(),
                    style=self.doc.style or None)
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
        # SSRM and reliability (a series of SSRM solves) both support cooperative cancel.
        supports_cancel = opts["analysis"] in ("ssrm", "reliability")
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

    def _on_fem_succeeded(self, bundle):
        self.doc.results["fem_solution"] = bundle
        if self.doc.path:
            try:
                from xslope.fem import export_fem_solution
                export_fem_solution(bundle["fem_data"], bundle["solution"],
                                    os.path.splitext(self.doc.path)[0],
                                    meta={"FS": bundle.get("FS"),
                                          "analysis": bundle.get("analysis"),
                                          # The strength-reduction factor shown in
                                          # the subplot titles (solution["F"]).
                                          "F": bundle["solution"].get("F")})
            except Exception:
                traceback.print_exc()
        self._show_fem_data(bundle["fem_data"])
        self._show_fem_results()
        if self.fem_results_canvas is not None:
            self.view_tabs.setCurrentWidget(self.fem_results_canvas)
        if bundle.get("analysis") == "reliability" and bundle.get("reliability"):
            r = bundle["reliability"]
            self.statusBar().showMessage(
                f"FEM reliability done — F_MLV = {r['F_MLV']:.3f}, "
                f"reliability = {r['reliability'] * 100:.2f}%, "
                f"Pf = {r['prob_failure'] * 100:.2f}%")
            QMessageBox.information(
                self, "FEM Reliability",
                f"F_MLV = {r['F_MLV']:.3f}\n"
                f"σ_F = {r['sigma_F']:.3f}\n"
                f"COV_F = {r['COV_F']:.3f}\n"
                f"β (lognormal) = {r['beta_ln']:.3f}\n"
                f"Reliability = {r['reliability'] * 100:.2f}%\n"
                f"Probability of failure = {r['prob_failure'] * 100:.2f}%\n\n"
                "The per-parameter ΔF table is in the Log pane; the FEM Results "
                "view shows the deformation at the most-likely values.")
        elif bundle.get("FS") is not None:
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
        self._rerender_fem_results()

    def _rerender_fem_results(self):
        bundle = self.doc.results.get("fem_solution")
        panel = self._display_panels.get(self.fem_results_canvas)
        if bundle and panel and self.fem_results_canvas is not None:
            try:
                self.fem_results_canvas.render_fem_results(
                    bundle["fem_data"], bundle["solution"], panel.options())
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
        runner = next((r for r in (self._runner, self._fem_runner, self._sens_runner)
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
        bundle = self.doc.results.get("lem_solution")
        rel = bundle.get("reliability") if bundle else None
        panel = self._display_panels.get(self.reliability_canvas)
        if rel and panel and self.reliability_canvas is not None:
            try:
                self.reliability_canvas.render_reliability(
                    self.doc.slope_data, rel, panel.options(), style=self.doc.style or None)
            except Exception:
                traceback.print_exc()

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
                                slope_data=self.doc.slope_data, app_mode=self._mode)
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_sens_opts[self._mode] = opts
        if opts["mode"] == "design" and not opts.get("param"):
            QMessageBox.warning(self, "Nothing to sweep",
                                "Pick a material and property to sweep.")
            return
        if opts["mode"] == "sensitivity" and not opts.get("params"):
            QMessageBox.warning(self, "Nothing to sweep",
                                "Add at least one parameter to the table.")
            return
        self.act_sensitivity.setEnabled(False)
        self.act_run.setEnabled(False)
        verb = "Design sweep" if opts["mode"] == "design" else "Sensitivity sweep"
        self.statusBar().showMessage(f"{verb} — {opts['method']} …")
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setVisible(True)
        self.cancel_btn.setEnabled(True)
        self.cancel_btn.setVisible(True)
        self._sens_runner = SensitivityRunner(self.doc.slope_data, opts, parent=self)
        self._sens_runner.succeeded.connect(self._on_sens_succeeded)
        self._sens_runner.failed.connect(self._on_sens_failed)
        self._sens_runner.cancelled.connect(self._on_sens_cancelled)
        self._sens_runner.progress.connect(self._on_run_progress)
        self._sens_runner.finished.connect(self._on_sens_finished)
        self._sens_runner.start()

    def _on_sens_succeeded(self, bundle):
        if bundle.get("kind") == "design":
            self.doc.results["design"] = bundle
            self._show_design()
            if self.design_canvas is not None:
                self.view_tabs.setCurrentWidget(self.design_canvas)
            if bundle.get("bracketed"):
                self.statusBar().showMessage(
                    f"Design — {bundle.get('output', 'FS')} = "
                    f"{bundle['target_fs']:g} at "
                    f"{bundle['param'].split(':')[-1]} = {bundle['crossing']:.4g}")
            else:
                self.statusBar().showMessage(bundle.get("message", "Design done."))
        else:
            self.doc.results["sensitivity"] = bundle
            self._show_sensitivity()
            if self.sens_canvas is not None:
                self.view_tabs.setCurrentWidget(self.sens_canvas)
            n = len(bundle.get("sweeps", {}))
            self.statusBar().showMessage(
                f"Sensitivity — {n} parameter(s); click a tornado bar for its curve.")

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
            # Double-click a bar to open that parameter's FS-vs-value curve.
            self.sens_canvas.set_pick_enabled(True)
            self.sens_canvas._hint_label.setText(
                "(double-click a bar to see its FS curve)")
            self.sens_canvas.picked.connect(self._on_tornado_pick)
        self._rerender_sensitivity()

    def _rerender_sensitivity(self):
        bundle = self.doc.results.get("sensitivity")
        if bundle and self.sens_canvas is not None:
            try:
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
                  "reliability_canvas", "sens_canvas", "sens_curve_canvas",
                  "design_canvas", "fem_data_canvas", "fem_results_canvas")
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
                    export_seep_u(mesh["nodes"], imported, path,
                                  sd.get("gamma_water") or 9.807)
                except Exception:
                    traceback.print_exc()
            else:
                remove(path)

        # FEM solution ({stem}_fem_nodes.csv / _fem_elements.csv / _fem_meta.json).
        fem = results.get("fem_solution")
        if fem:
            try:
                from xslope.fem import export_fem_solution
                export_fem_solution(fem["fem_data"], fem["solution"], stem,
                                    meta={"FS": fem.get("FS"),
                                          "analysis": fem.get("analysis"),
                                          "F": fem["solution"].get("F")})
            except Exception:
                traceback.print_exc()
        else:
            for f in (f"{stem}_fem_nodes.csv", f"{stem}_fem_elements.csv",
                      f"{stem}_fem_meta.json"):
                remove(f)

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
        if self._sens_runner is not None and self._sens_runner.isRunning():
            self._sens_runner.cancel()      # sweep stops at the next point
            self._sens_runner.wait(15000)
        # Stop the persistent mesh thread (lets an in-flight build finish first).
        self._mesh_thread.quit()
        self._mesh_thread.wait(10000)
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
