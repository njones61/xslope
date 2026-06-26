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


def _new_material(name):
    """A single material with all keys the loader produces (mirrors fileio)."""
    return {"name": name, "gamma": 125.0, "option": "mc", "c": 500.0, "phi": 0.0,
            "cp": 0.0, "r_elev": 0.0, "d": 0.0, "psi": 0.0, "u": "none",
            "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0, "sigma_cp": 0.0,
            "sigma_d": 0.0, "sigma_psi": 0.0, "k1": 0.0, "k2": 0.0, "alpha": 0.0,
            "kr0": 0.0, "h0": 0.0, "E": 0.0, "nu": 0.0}


def new_slope_data():
    """Build a minimal but valid in-memory ``slope_data`` for a new project.

    A blank template can't be ``load``ed (the loader rejects empty geometry and
    requires a failure surface), so New seeds a single-material slope with one
    profile line and one circle — the proven ``xslope_simple1`` example geometry —
    then derives the rest exactly as the loader does. The user edits everything
    from there via the input editors; Save As writes it back through the template.
    """
    from shapely.geometry import Polygon
    from xslope.mesh import build_polygons
    from xslope.fileio import build_ground_surface_from_polygons

    materials = [_new_material("soil")]
    profile_lines = [{"mat_id": 0, "coords": [(0.0, 0.0), (20.0, 20.0), (100.0, 20.0)]}]
    max_depth = 0.0
    polys = [{"polygon": Polygon(p["coords"]), "mat_id": p["mat_id"]}
             for p in build_polygons(slope_data={"profile_lines": profile_lines,
                                                  "max_depth": max_depth})]
    ground_surface, domain_polygon = build_ground_surface_from_polygons(polys)

    return {
        "template_version": None,
        "gamma_water": 62.4, "tcrack_depth": 0.0, "tcrack_water": 0.0,
        "k_seismic": 0.0, "max_depth": max_depth,
        "profile_lines": profile_lines, "polygons": polys,
        "domain_polygon": domain_polygon, "ground_surface": ground_surface,
        "tcrack_surface": None, "materials": materials,
        "piezo_line": [], "piezo_line2": [],
        # One default circle so the project is valid (the loader requires a failure
        # surface). R = Yo - Depth, matching the loader's Depth-option collapse.
        "circular": True,
        "circles": [{"Xo": 10.0, "Yo": 40.0, "Option": "Depth", "Depth": 0.0,
                     "Xi": 0.0, "Yi": 0.0, "R": 40.0}],
        "non_circ": [],
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
        self._undo.clear()
        self._redo.clear()
        self._dirty = False
        self.loaded.emit()
        self.dirty_changed.emit(False)

    def new(self):
        """Start a fresh, unsaved project from a minimal in-memory skeleton."""
        self.slope_data = new_slope_data()
        self.path = None          # no source file until first Save As
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
