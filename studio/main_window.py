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

from PySide6.QtCore import Qt, QObject, QSettings, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QComboBox, QDockWidget, QFileDialog, QLabel, QMainWindow, QMessageBox,
    QPlainTextEdit, QTabWidget, QToolBar, QTreeWidget, QTreeWidgetItem, QWidget,
)

from .canvas import MplCanvas
from .dialogs import RunLemDialog
from .document import ProjectDocument
from .editors import CATEGORY_EDITORS
from .runners import LemRunner

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
        self.search_canvas = None
        self.solution_canvas = None
        self.view_tabs = QTabWidget()
        self.view_tabs.addTab(self.canvas, "Inputs")
        self.view_tabs.currentChanged.connect(self._on_view_tab_changed)
        self.setCentralWidget(self.view_tabs)

        self._mode = "lem"
        self._runner = None
        self._last_lem_opts = {}
        self._recent = [p for p in (self.settings.value("recent_files") or []) if p]

        self._make_inputs_dock()
        self._make_log_dock()
        self._make_actions()
        self._make_menus()
        self._make_toolbar()
        self._update_recent_menu()
        self._update_title()
        self._install_log_capture()
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

    def _make_log_dock(self):
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(5000)
        dock = QDockWidget("Log", self)
        dock.setObjectName("log_dock")
        dock.setWidget(self.log)
        self.addDockWidget(Qt.BottomDockWidgetArea, dock)
        self.log_dock = dock

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
        self.act_run_lem = QAction("Run &LEM…", self, enabled=False, triggered=self.run_lem)

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
        m_run.addAction(self.act_run_lem)

        m_view = mb.addMenu("&View")
        m_view.addAction(self.inputs_dock.toggleViewAction())
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
        tb.addAction(self.act_run_lem)
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
        self.mode_combo.setCurrentIndex([m for _, m in MODES].index(self._mode))
        self.act_save.setEnabled(True)
        self.act_save_as.setEnabled(True)
        self.act_run_lem.setEnabled(True)
        self._clear_result_tabs()
        self.canvas.reset_fit()       # fit the fresh file to the window
        self._render()
        self._populate_inputs_tree()
        self._update_title()
        n = len(self.doc.slope_data.get("materials", []))
        name = os.path.basename(self.doc.path) if self.doc.path else "untitled"
        self.statusBar().showMessage(
            f"Loaded {name} — {n} material(s). "
            f"Click an underlined input to edit it.")

    def _render(self):
        if not self.doc.is_open:
            self.canvas.clear()
            return
        try:
            self.canvas.render_inputs(self.doc.slope_data, mode=self._mode)
        except Exception:
            traceback.print_exc()
        self.act_undo.setEnabled(self.doc.can_undo())
        self.act_redo.setEnabled(self.doc.can_redo())

    def _on_mode_changed(self, index):
        self._mode = MODES[index][1]
        if self.doc.is_open:
            self._render()
            self._populate_inputs_tree()

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

    # --- LEM analysis ----------------------------------------------------
    def run_lem(self):
        if not self.doc.is_open or self._runner is not None:
            return
        dlg = RunLemDialog(self, defaults=self._last_lem_opts)
        if not dlg.exec():
            return
        opts = dlg.options()
        self._last_lem_opts = opts
        self.act_run_lem.setEnabled(False)
        verb = "Searching" if opts["analysis"] == "auto_search" else "Running"
        self.statusBar().showMessage(f"{verb} {opts['method']} …")
        self._runner = LemRunner(self.doc.slope_data, opts, parent=self)
        self._runner.succeeded.connect(self._on_lem_succeeded)
        self._runner.failed.connect(self._on_lem_failed)
        self._runner.finished.connect(self._on_lem_finished)
        self._runner.start()

    def _on_lem_succeeded(self, bundle):
        self.doc.results["lem_solution"] = bundle
        if bundle.get("search"):
            self._show_search(bundle["search"])
        self._show_solution(bundle)
        # For a search, lead with the Search view (all trial surfaces + critical);
        # for a single surface, the Solution view is the result.
        lead = self.search_canvas if bundle.get("search") else self.solution_canvas
        if lead is not None:
            self.view_tabs.setCurrentWidget(lead)
        res = bundle["results"]
        self.statusBar().showMessage(
            f"LEM done — {res.get('method')} FS = {res.get('FS'):.3f}")

    def _on_lem_failed(self, message):
        QMessageBox.warning(self, "LEM run failed", message)
        self.statusBar().showMessage("LEM run failed.")

    def _on_lem_finished(self):
        self.act_run_lem.setEnabled(self.doc.is_open)
        if self._runner is not None:
            self._runner.deleteLater()
            self._runner = None

    def _show_search(self, search):
        if self.search_canvas is None:
            self.search_canvas = MplCanvas(self)
            # Keep order Inputs → Search → Solution.
            self.view_tabs.insertTab(1, self.search_canvas, "LEM · Search")
        try:
            self.search_canvas.render_search(self.doc.slope_data, search)
        except Exception:
            traceback.print_exc()

    def _show_solution(self, bundle):
        if self.solution_canvas is None:
            self.solution_canvas = MplCanvas(self)
            self.view_tabs.addTab(self.solution_canvas, "LEM · Solution")
        try:
            self.solution_canvas.render_solution(
                self.doc.slope_data, bundle["slice_df"],
                bundle["failure_surface"], bundle["results"])
        except Exception:
            traceback.print_exc()

    def _clear_result_tabs(self):
        """Drop result views (e.g. on opening another file) so they don't show
        stale results from the previous project."""
        for attr in ("search_canvas", "solution_canvas"):
            canvas = getattr(self, attr)
            if canvas is not None:
                idx = self.view_tabs.indexOf(canvas)
                if idx >= 0:
                    self.view_tabs.removeTab(idx)
                canvas.deleteLater()
                setattr(self, attr, None)

    def _on_view_tab_changed(self, index):
        w = self.view_tabs.widget(index)
        if hasattr(w, "ensure_fitted"):
            w.ensure_fitted()

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
            self._runner.wait(5000)   # let an in-flight solve finish before teardown
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
