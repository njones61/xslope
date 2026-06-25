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

from PySide6.QtCore import Qt, QSettings
from PySide6.QtGui import QAction, QActionGroup, QKeySequence
from PySide6.QtWidgets import (
    QComboBox, QDockWidget, QFileDialog, QLabel, QMainWindow, QMessageBox,
    QPlainTextEdit, QToolBar, QTreeWidget, QTreeWidgetItem, QWidget,
)

from .canvas import MplCanvas
from .document import ProjectDocument

APP_NAME = "XSlope Studio"
ORG_NAME = "XSlope"
MAX_RECENT = 8
MODES = [("LEM", "lem"), ("Seepage", "seep"), ("FEM", "fem")]


class _LogStream:
    """Minimal stdout/stderr tee: forwards to the original stream and the log pane."""

    def __init__(self, widget, original):
        self._widget = widget
        self._original = original

    def write(self, text):
        if self._original is not None:
            try:
                self._original.write(text)
            except Exception:
                pass
        text = text.rstrip("\n")
        if text:
            self._widget.appendPlainText(text)

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

        self.canvas = MplCanvas(self)
        self.setCentralWidget(self.canvas)

        self._mode = "lem"
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
        self.inputs_tree.setHeaderLabels(["Input", "Count / Value"])
        self.inputs_tree.setColumnWidth(0, 180)
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
        self.act_open = QAction("&Open…", self, shortcut=QKeySequence.Open,
                                triggered=self.open_dialog)
        self.act_quit = QAction("&Quit", self, shortcut=QKeySequence.Quit,
                                triggered=self.close)
        self.act_undo = QAction("&Undo", self, shortcut=QKeySequence.Undo,
                                triggered=self.doc.undo, enabled=False)
        self.act_redo = QAction("&Redo", self, shortcut=QKeySequence.Redo,
                                triggered=self.doc.redo, enabled=False)
        self.act_about = QAction("&About", self, triggered=self._about)
        # Save actions exist but are disabled until Phase 2 wires editing.
        self.act_save = QAction("&Save", self, shortcut=QKeySequence.Save, enabled=False)
        self.act_save_as = QAction("Save &As…", self, enabled=False)

    def _make_menus(self):
        mb = self.menuBar()

        m_file = mb.addMenu("&File")
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

        m_view = mb.addMenu("&View")
        m_view.addAction(self.inputs_dock.toggleViewAction())
        m_view.addAction(self.log_dock.toggleViewAction())

        m_help = mb.addMenu("&Help")
        m_help.addAction(self.act_about)

    def _make_toolbar(self):
        tb = QToolBar("Main", self)
        tb.setObjectName("main_toolbar")
        self.addToolBar(tb)
        tb.addAction(self.act_open)
        tb.addSeparator()
        tb.addWidget(QLabel(" Mode: "))
        self.mode_combo = QComboBox()
        for label, _ in MODES:
            self.mode_combo.addItem(label)
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        tb.addWidget(self.mode_combo)

    # --- open / recent ---------------------------------------------------
    def open_dialog(self):
        start = os.path.dirname(self._recent[0]) if self._recent else ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Open XSlope input file", start, "Excel files (*.xlsx);;All files (*)")
        if path:
            self.open_path(path)

    def open_path(self, path):
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
        self._render()
        self._populate_inputs_tree()
        self._update_title()
        n = len(self.doc.slope_data.get("materials", []))
        self.statusBar().showMessage(
            f"Loaded {os.path.basename(self.doc.path)} — {n} material(s)")

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

        def add(name, value):
            QTreeWidgetItem(self.inputs_tree, [name, str(value)])

        add("Unit weight of water", d.get("gamma_water"))
        add("Tension crack depth", d.get("tcrack_depth"))
        add("Seismic coefficient k", d.get("k_seismic"))
        add("Materials", len(d.get("materials", [])))
        add("Profile lines", len(d.get("profile_lines") or []))
        add("Polygons", len(d.get("polygons") or []))
        add("Circles", len(d.get("circles") or []))
        add("Non-circular pts", len(d.get("non_circ") or []))
        add("Piezo line pts", len(d.get("piezo_line") or []))
        add("Distributed loads", len(d.get("dloads") or []))
        add("Reinforcement lines", len(d.get("reinforcement_lines") or []))
        add("Piles", len(d.get("pile_lines") or []))
        sbc = d.get("seepage_bc") or {}
        add("Seepage heads", len(sbc.get("specified_heads", [])))
        add("Mesh", "yes" if d.get("mesh") is not None else "no")
        self.inputs_tree.expandAll()

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
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        super().closeEvent(event)
