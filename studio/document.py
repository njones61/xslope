"""ProjectDocument — the single source of truth for an open project.

Holds the in-memory ``slope_data`` dict, the source file path, a dirty flag, and
a snapshot-based undo/redo stack (the dicts are small, so a deep copy per edit is
cheap and gives full undo cheaply — see plan §4/§12). The GUI mutates
``slope_data`` through this object and listens to its signals to re-render.
"""

from __future__ import annotations

import copy

from PySide6.QtCore import QObject, Signal

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx


def new_slope_data():
    """Build an empty but structurally complete in-memory ``slope_data`` for a new
    project — every category present but empty, no geometry yet. The canvas stays
    blank until the user adds inputs through the editors (start with a material and
    a profile line; the profile editor rebuilds the ground surface). Save As writes
    it back through the template.
    """
    return {
        "template_version": None,
        "gamma_water": 62.4, "tcrack_depth": 0.0, "tcrack_water": 0.0,
        "k_seismic": 0.0, "max_depth": 0.0,
        "profile_lines": [], "polygons": [],
        "domain_polygon": None, "ground_surface": None,
        "tcrack_surface": None, "materials": [],
        "piezo_line": [], "piezo_line2": [],
        "circular": False, "circles": [], "non_circ": [],
        "dloads": [], "dloads2": [],
        "reinforce_lines": [], "reinforcement_lines": [], "pile_lines": [],
        "seepage_bc": {"specified_heads": [], "exit_face": []},
        "seepage_bc2": {"specified_heads": [], "exit_face": []},
        "has_seepage_bc2": False, "mesh": None,
    }


class ProjectDocument(QObject):
    loaded = Signal()             # a project was loaded/created (full reset)
    changed = Signal()            # slope_data mutated (re-render)
    dirty_changed = Signal(bool)  # unsaved-changes flag toggled

    def __init__(self, parent=None):
        super().__init__(parent)
        self.slope_data = None
        self.path = None          # source .xlsx path (None until first save)
        self.results = {}         # cached analysis results (e.g. 'lem_solution')
        self._dirty = False
        self._undo = []           # list of slope_data snapshots (deep copies)
        self._redo = []

    # --- state -----------------------------------------------------------
    @property
    def is_open(self) -> bool:
        return self.slope_data is not None

    @property
    def dirty(self) -> bool:
        return self._dirty

    def _set_dirty(self, value: bool):
        if value != self._dirty:
            self._dirty = value
            self.dirty_changed.emit(value)

    # --- load ------------------------------------------------------------
    def load(self, path):
        """Load a project from an Excel file. Raises ValueError on bad input."""
        self.slope_data = load_slope_data(str(path))
        self.path = str(path)
        self.results.clear()
        self._undo.clear()
        self._redo.clear()
        self._dirty = False
        self.loaded.emit()
        self.dirty_changed.emit(False)

    def new(self):
        """Start a fresh, unsaved project from a minimal in-memory skeleton."""
        self.slope_data = new_slope_data()
        self.path = None          # no source file until first Save As
        self.results.clear()
        self._undo.clear()
        self._redo.clear()
        self._dirty = True        # nothing on disk yet
        self.loaded.emit()
        self.dirty_changed.emit(True)

    # --- editing / snapshot undo ----------------------------------------
    def begin_edit(self):
        """Snapshot the current state before a mutation so it can be undone."""
        if self.slope_data is not None:
            self._undo.append(copy.deepcopy(self.slope_data))
            self._redo.clear()

    def commit_edit(self):
        """Signal that an in-place mutation of slope_data is complete."""
        self._set_dirty(True)
        self.changed.emit()

    def rollback_edit(self):
        """Discard the in-place mutation made since ``begin_edit`` by restoring the
        snapshot. Used when an assistant ``run_python`` edit raised, so a partial
        (and possibly repeated) mutation from a failed snippet doesn't stick."""
        if self._undo:
            self.slope_data = self._undo.pop()
            self.changed.emit()

    def cancel_edit(self):
        """Drop the ``begin_edit`` snapshot without committing — for a code run that
        turned out not to change anything (e.g. a read-only query), so it isn't
        recorded as an edit or marked dirty."""
        if self._undo:
            self._undo.pop()

    def can_undo(self) -> bool:
        return bool(self._undo)

    def can_redo(self) -> bool:
        return bool(self._redo)

    def undo(self):
        if not self._undo:
            return
        self._redo.append(copy.deepcopy(self.slope_data))
        self.slope_data = self._undo.pop()
        self._set_dirty(True)
        self.changed.emit()

    def redo(self):
        if not self._redo:
            return
        self._undo.append(copy.deepcopy(self.slope_data))
        self.slope_data = self._redo.pop()
        self._set_dirty(True)
        self.changed.emit()

    # --- save (scaffolded; wired into the UI in Phase 2) -----------------
    def save(self, path=None, template=None):
        target = str(path) if path else self.path
        if target is None:
            raise ValueError("No path to save to (use Save As).")
        save_slope_data_to_xlsx(self.slope_data, target, template=template)
        self.path = target
        self._set_dirty(False)
