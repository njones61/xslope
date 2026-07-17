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
        "line_loads": [],
        "seepage_bc": {"specified_heads": [], "exit_face": []},
        "seepage_bc2": {"specified_heads": [], "exit_face": []},
        "has_seepage_bc2": False, "mesh": None,
    }


def _blank_material(name):
    """A placeholder material carrying only a name, every property zeroed — the
    same key set ``load_slope_data`` produces. Used by DXF import, which yields
    layer names (→ material names) but no properties for the user to fill in."""
    return {
        "name": str(name), "gamma": 0.0, "gamma_sat": None, "option": "",
        "c": 0.0, "phi": 0.0, "cp": 0.0, "r_elev": 0.0, "d": 0, "psi": 0,
        "pow_a": 0.0, "pow_b": 0.0, "pow_c": 0.0, "pow_d": 0.0,
        "hb_sci": 0.0, "hb_gsi": 0.0, "hb_mi": 0.0, "hb_d": 0.0,
        "u": "none", "ru": 0.0,
        "sigma_gamma": 0.0, "sigma_c": 0.0, "sigma_phi": 0.0, "sigma_cp": 0.0,
        "sigma_d": 0.0, "sigma_psi": 0.0, "k1": 0.0, "k2": 0.0, "alpha": 0.0,
        "unsat": "lf", "kr0": 0.0, "h0": 0.0, "vg_a": 0.0, "vg_n": 0.0,
        "E": 0.0, "nu": 0.0,
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
        self.style = {}           # sparse feature-style deltas (see xslope.style; §8a)
        self._dirty = False
        # Undo/redo stacks of (label, slope_data snapshot) records. The label is a
        # short human tag for the edit ("Edit Materials", "Assistant: dloads"),
        # surfaced in the toolbar history dropdown (plan §Phase 2 undo hardening).
        self._undo = []
        self._redo = []
        # Saved-state marker: the undo position (== len(self._undo)) the document
        # was last saved at, so undo/redo back to that point clears the dirty flag
        # (Office-style). None = the saved state is unreachable (never saved, or its
        # branch was discarded by a later edit). `style` lives outside the undo
        # stack, so its dirtiness is tracked separately and OR-ed in.
        self._clean_index = 0
        self._style_dirty = False

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

    def _refresh_dirty(self):
        """Recompute dirty from the undo position vs the saved marker (plus the
        out-of-band style flag). Call after any settled history transition."""
        position_dirty = self._clean_index is None or len(self._undo) != self._clean_index
        self._set_dirty(self._style_dirty or position_dirty)

    def _snapshot(self):
        """A snapshot of ``slope_data`` for the undo stack. Everything is deep-copied
        **except** ``mesh``, whose (potentially large) object is shared by reference:
        the mesh is only ever replaced wholesale (Build Mesh / clear), never mutated
        in place, so sharing is safe and keeps snapshots cheap once a mesh exists."""
        sd = self.slope_data
        mesh = sd.get("mesh")
        if mesh is None:
            return copy.deepcopy(sd)
        sd["mesh"] = None
        try:
            snap = copy.deepcopy(sd)
        finally:
            sd["mesh"] = mesh
        snap["mesh"] = mesh
        return snap

    # --- load ------------------------------------------------------------
    def load(self, path):
        """Load a project from an Excel file. Raises ValueError on bad input."""
        self.slope_data = load_slope_data(str(path))
        self.path = str(path)
        self.results.clear()
        self.style = self._read_styles_sidecar(self.path)
        self._undo.clear()
        self._redo.clear()
        self._clean_index = 0
        self._style_dirty = False
        self._dirty = False
        self.loaded.emit()
        self.dirty_changed.emit(False)

    # --- styles (per-project sidecar; §8a) -------------------------------
    @staticmethod
    def _styles_sidecar_path(xlsx_path):
        import os
        stem = os.path.splitext(str(xlsx_path))[0]
        return stem + "_styles.json"

    def _read_styles_sidecar(self, xlsx_path):
        """Read `{stem}_styles.json` (sparse style deltas) next to the xlsx, or {}."""
        import json
        path = self._styles_sidecar_path(xlsx_path)
        try:
            with open(path) as fh:
                data = json.load(fh)
            return data if isinstance(data, dict) else {}
        except (FileNotFoundError, ValueError):
            return {}

    def _write_styles_sidecar(self, xlsx_path):
        """Write the project's style deltas to `{stem}_styles.json`; remove the file
        when there are no deltas, so unstyled projects don't litter sidecars."""
        import json
        import os
        path = self._styles_sidecar_path(xlsx_path)
        if self.style:
            with open(path, "w") as fh:
                json.dump(self.style, fh, indent=2)
        elif os.path.exists(path):
            os.remove(path)

    def set_style(self, style):
        """Replace the project's style deltas (sparse, see xslope.style), mark dirty,
        and trigger a re-render. Styles live outside the undo stack, so their
        dirtiness is tracked on its own flag (cleared on save)."""
        self.style = style or {}
        self._style_dirty = True
        self._refresh_dirty()
        self.changed.emit()

    def new(self):
        """Start a fresh, unsaved project from a minimal in-memory skeleton.

        An untouched empty project is NOT dirty — there's nothing to lose, so
        quitting/opening shouldn't prompt to save. It becomes dirty on the first
        edit (commit_edit). DXF import, which actually populates content, marks
        itself dirty.
        """
        self.slope_data = new_slope_data()
        self.path = None          # no source file until first Save As
        self.results.clear()
        self.style = {}
        self._undo.clear()
        self._redo.clear()
        self._clean_index = 0
        self._style_dirty = False
        self._dirty = False
        self.loaded.emit()
        self.dirty_changed.emit(False)

    def read_dxf_layers(self, dxf_path):
        """Read ALL geometry from a DXF, grouped/classified by layer, without
        mutating the project (the engine's ``read_dxf_layers``). The feature-aware
        import wizard uses this to show each layer, collect a target mapping, then
        call ``build_from_dxf_mapping``."""
        from xslope.cad import read_dxf_layers
        layers, warnings = read_dxf_layers(str(dxf_path))
        if not layers:
            raise ValueError("No geometry found in the DXF.")
        return layers, warnings

    def build_from_dxf_mapping(self, layers, mapping):
        """Build a fresh project from classified DXF layers under a per-layer
        mapping. ``layers`` is ``read_dxf_layers`` output; ``mapping`` maps each
        layer name to ``{'target': one of cad.DXF_TARGETS, 'material': name|None}``.

        Each layer's geometry is routed to the chosen input feature; non-geometric
        properties (load magnitudes, reinforcement strengths, material properties,
        circle depth) come in as editable placeholders. Profile-based if any layer
        maps to a profile line, else polygon-based. REPLACES the current project
        (callers confirm discard); result is unsaved. Returns caveat strings."""
        from xslope.cad import fit_circle, _stitch_lines, chain_segments
        from xslope.fileio import build_reinforce_lines
        from studio.editors import _resync_geometry

        notes = []
        sd = new_slope_data()

        names = []                        # materials used by zone/profile mappings
        for m in mapping.values():
            if m.get("target") in ("material_zone", "profile") and m.get("material"):
                if m["material"] not in names:
                    names.append(m["material"])
        name_to_idx = {n: i for i, n in enumerate(names)}
        sd["materials"] = [_blank_material(n) for n in names]

        profile_lines, polygons, piezo_polylines = [], [], []
        reinforcement_lines, dload_blocks, circles = [], [], []
        placeholder = False

        for lyr, geom in layers.items():
            t = (mapping.get(lyr) or {}).get("target", "ignore")
            if t == "ignore":
                continue
            mid = name_to_idx.get((mapping.get(lyr) or {}).get("material"), 0)
            if t == "material_zone":
                for coords in geom["closed"] + geom["open"]:   # open rings force-close
                    if len(coords) >= 3:
                        polygons.append({"coords": list(coords), "mat_id": mid})
                rings, _ = _stitch_lines(geom["lines"]) if geom["lines"] else ([], 0)
                polygons += [{"coords": r, "mat_id": mid} for r in rings if len(r) >= 3]
            elif t == "profile":
                for coords in geom["open"] + geom["closed"]:
                    if len(coords) >= 2:
                        profile_lines.append({"coords": list(coords), "mat_id": mid})
            elif t == "piezo":
                piezo_polylines += [c for c in geom["open"] + geom["closed"] if len(c) >= 2]
            elif t == "dload":
                for coords in geom["open"] + geom["closed"]:
                    if len(coords) >= 2:
                        dload_blocks.append([{"X": x, "Y": y, "Normal": 0.0} for x, y in coords])
                        placeholder = True
            elif t == "reinforce":
                # Reconnect loose LINE segments into chains (export split each line
                # into per-point segments); each chain / open polyline → one
                # reinforcement line, defined by its two extreme endpoints.
                chains = [c for c in geom["open"] if len(c) >= 2]
                chains += chain_segments(geom["lines"]) if geom["lines"] else []
                for chain in chains:
                    a, b = chain[0], chain[-1]
                    reinforcement_lines.append(
                        {"x1": a[0], "y1": a[1], "x2": b[0], "y2": b[1], "t_max": 0.0,
                         "t_res": 0.0, "lp1": 0.0, "lp2": 0.0, "E": 0.0, "area": 0.0})
                placeholder = placeholder or bool(chains)
            elif t == "circles":
                for cx, cy, r in geom["circles"]:
                    circles.append({"Xo": cx, "Yo": cy, "R": r, "Depth": 0.0})
                for pt in geom["points"]:
                    r = fit_circle(pt, geom["open"])
                    if r:
                        circles.append({"Xo": pt[0], "Yo": pt[1], "R": r, "Depth": 0.0})
                    else:
                        notes.append(f"circle center on layer '{lyr}' had no arc to "
                                     "size its radius — skipped")

        if not (profile_lines or polygons or piezo_polylines or dload_blocks
                or reinforcement_lines or circles):
            raise ValueError("Nothing selected to import (all layers set to Ignore).")

        if profile_lines:
            sd["profile_lines"] = profile_lines
            if polygons:
                notes.append("both profile lines and material zones were mapped; used "
                             "the profile lines (zones are derived from them)")
        elif polygons:
            sd["polygons"] = polygons          # coords form; _resync normalizes
        if piezo_polylines:
            sd["piezo_line"] = [tuple(p) for p in piezo_polylines[0]]
            if len(piezo_polylines) > 1:
                sd["piezo_line2"] = [tuple(p) for p in piezo_polylines[1]]
            if len(piezo_polylines) > 2:
                notes.append(f"{len(piezo_polylines) - 2} extra piezo polyline(s) ignored (max 2)")
        if dload_blocks:
            sd["dloads"] = dload_blocks
        if reinforcement_lines:
            sd["reinforcement_lines"] = reinforcement_lines
            sd["reinforce_lines"] = build_reinforce_lines(reinforcement_lines)
        if circles:
            sd["circles"] = circles
            sd["circular"] = True
        if placeholder:
            notes.append("distributed-load magnitudes and reinforcement strengths "
                         "imported as 0 — set them in the editors")

        _resync_geometry(sd)               # ground surface / domain / t-crack

        self.slope_data = sd
        self.path = None
        self.results.clear()
        self.style = {}
        self._undo.clear()
        self._redo.clear()
        self._clean_index = None      # freshly imported, never saved -> dirty
        self._style_dirty = False
        self._dirty = True
        self.loaded.emit()
        self.dirty_changed.emit(True)
        return notes

    def read_gsz_analyses(self, gsz_path):
        """Open a GeoStudio .gsz and list the analyses it holds, without mutating
        the project (the engine's ``read_gsz``). A .gsz routinely bundles several
        analyses over one geometry, and they can differ in materials as well as in
        slip surface — so the import dialog shows them and the user picks one, then
        calls ``build_from_gsz``.

        Returns ``(gsz, analyses)``; the opaque ``gsz`` is passed straight back."""
        from xslope.geostudio import read_gsz, list_analyses
        gsz = read_gsz(str(gsz_path))
        analyses = list_analyses(gsz)
        if not analyses:
            raise ValueError("This GeoStudio file defines no analyses.")
        return gsz, analyses

    def build_from_gsz(self, gsz, analysis_id):
        """Build a fresh project from one analysis of a parsed GeoStudio .gsz.

        Unlike DXF, a .gsz is semantically complete — material zones, strengths and
        water conditions arrive already identified, so there is nothing to map. What
        SLOPE/W models and xslope cannot is returned as caveats rather than dropped
        in silence. REPLACES the current project (callers confirm discard); result is
        unsaved. Returns caveat strings."""
        from xslope.geostudio import gsz_to_slope_data, gsz_style
        from studio.editors import _resync_geometry

        sd, caveats = gsz_to_slope_data(gsz, analysis_id)
        _resync_geometry(sd)               # ground surface / domain / t-crack

        self.slope_data = sd
        self.path = None
        self.results.clear()
        self.style = gsz_style(gsz, analysis_id)   # keep the user's own material colours
        self._undo.clear()
        self._redo.clear()
        self._clean_index = None      # freshly imported, never saved -> dirty
        self._style_dirty = bool(self.style)   # so the colours reach the sidecar on save
        self._dirty = True
        self.loaded.emit()
        self.dirty_changed.emit(True)
        return caveats

    def read_slide2_scenarios(self, slide2_path):
        """Open a Rocscience Slide2 model (.sli/.slim/.slmd) and list the scenarios
        it holds, without mutating the project (the engine's ``read_slide2``). A
        ``.slmd`` routinely bundles several scenarios over one geometry (a base case
        plus variants), the way a GeoStudio ``.gsz`` bundles several analyses — so
        the import dialog shows them and the user picks one, then calls
        ``build_from_slide2``. A ``.sli``/``.slim`` has exactly one.

        Returns ``(d, scenarios)``; the opaque ``d`` is passed straight back."""
        from xslope.slide2 import read_slide2, list_scenarios
        d = read_slide2(str(slide2_path))
        scenarios = list_scenarios(d)
        if not scenarios:
            raise ValueError("This Slide2 file defines no scenarios.")
        return d, scenarios

    def build_from_slide2(self, d, scenario):
        """Build a fresh project from one scenario of a parsed Slide2 model.

        Like a .gsz, a Slide2 file is semantically complete — material zones,
        strengths and water conditions arrive already identified, so there is
        nothing to map. What Slide2 models and xslope cannot is returned as
        caveats rather than dropped in silence. Most Slide2 tutorial models define
        a SEARCH rather than a specified surface, so the import routinely arrives
        with no failure circle — the document still opens for editing (the caveat
        says so) and the user defines circles afterward, same as an unsolved .gsz.
        REPLACES the current project (callers confirm discard); result is unsaved.
        Returns caveat strings."""
        from xslope.slide2 import slide2_to_slope_data
        from studio.editors import _resync_geometry

        sd, caveats = slide2_to_slope_data(d, scenario)
        _resync_geometry(sd)               # ground surface / domain / t-crack

        self.slope_data = sd
        self.path = None
        self.results.clear()
        self.style = {}                    # Slide2 carries no colour info to keep
        self._undo.clear()
        self._redo.clear()
        self._clean_index = None      # freshly imported, never saved -> dirty
        self._style_dirty = False
        self._dirty = True
        self.loaded.emit()
        self.dirty_changed.emit(True)
        return caveats

    def build_from_fez(self, fez_path):
        """Build a fresh project from a Rocscience RS2 model (.fez).

        Unlike a .gsz or a Slide2 file, a .fez holds exactly one model — RS2 has no
        concept of bundling scenarios/analyses in one file — so there is no picker
        step: read and convert happen in one call (the engine's ``read_fez`` +
        ``fez_to_slope_data``). Material zones, strengths and water conditions arrive
        already identified, so there is nothing to map; what RS2 defines and xslope
        cannot (its SSR settings, joints, reinforcement, loads) is returned as
        caveats rather than dropped in silence. RS2's slope-stability answer is a
        finite-element SSR field, not a limit-equilibrium surface, so the import
        NEVER carries a failure surface — the document still opens for editing (the
        caveat says so) and the user defines circles afterward, same as an unsolved
        .gsz or a search-only Slide2 scenario. REPLACES the current project (callers
        confirm discard); result is unsaved. Returns caveat strings."""
        from xslope.rs2 import read_fez, fez_to_slope_data
        from studio.editors import _resync_geometry

        d = read_fez(str(fez_path))
        sd, caveats = fez_to_slope_data(d)
        _resync_geometry(sd)               # ground surface / domain / t-crack

        self.slope_data = sd
        self.path = None
        self.results.clear()
        self.style = {}                    # RS2 carries no colour info to keep
        self._undo.clear()
        self._redo.clear()
        self._clean_index = None      # freshly imported, never saved -> dirty
        self._style_dirty = False
        self._dirty = True
        self.loaded.emit()
        self.dirty_changed.emit(True)
        return caveats

    # --- editing / snapshot undo ----------------------------------------
    # Cap on the undo/redo depth: each snapshot deep-copies slope_data (the mesh is
    # shared, not copied — see _snapshot), so an unbounded stack would still grow on
    # a long session.
    MAX_HISTORY = 100

    def begin_edit(self, label="Edit"):
        """Snapshot the current state before a mutation so it can be undone. ``label``
        is a short tag describing the edit, shown in the history dropdown."""
        if self.slope_data is None:
            return
        # A new edit discards the redo branch. If the saved (clean) state lived on
        # that branch, it's now unreachable -> mark it lost so we stay dirty.
        if self._clean_index is not None and self._clean_index > len(self._undo):
            self._clean_index = None
        self._undo.append((label, self._snapshot()))
        if len(self._undo) > self.MAX_HISTORY:
            self._undo.pop(0)             # drop oldest; positions shift down by 1
            if self._clean_index is not None:
                self._clean_index -= 1
                if self._clean_index < 0:
                    self._clean_index = None   # saved state scrolled off the stack
        self._redo.clear()

    def commit_edit(self):
        """Signal that an in-place mutation of slope_data is complete."""
        self._refresh_dirty()
        self.changed.emit()

    def rollback_edit(self):
        """Discard the in-place mutation made since ``begin_edit`` by restoring the
        snapshot. Used when an assistant ``run_python`` edit raised, so a partial
        (and possibly repeated) mutation from a failed snippet doesn't stick."""
        if self._undo:
            _, snap = self._undo.pop()
            self.slope_data = snap
            self._refresh_dirty()
            self.changed.emit()

    def cancel_edit(self):
        """Drop the ``begin_edit`` snapshot without committing — for a code run that
        turned out not to change anything (e.g. a read-only query), so it isn't
        recorded as an edit or marked dirty."""
        if self._undo:
            self._undo.pop()

    def relabel_last_edit(self, label):
        """Rename the most recent undo step — used by the assistant to tag a step
        with what it changed once that's known (after the snippet ran)."""
        if self._undo:
            _, snap = self._undo[-1]
            self._undo[-1] = (label, snap)

    def can_undo(self) -> bool:
        return bool(self._undo)

    def can_redo(self) -> bool:
        return bool(self._redo)

    def undo_labels(self):
        """Edit labels on the undo stack, most-recent first (for the history menu)."""
        return [label for label, _ in reversed(self._undo)]

    def redo_labels(self):
        """Edit labels on the redo stack, next-to-redo first (for the history menu)."""
        return [label for label, _ in reversed(self._redo)]

    def _undo_one(self):
        label, snap = self._undo.pop()
        self._redo.append((label, self._snapshot()))
        self.slope_data = snap

    def _redo_one(self):
        label, snap = self._redo.pop()
        self._undo.append((label, self._snapshot()))
        self.slope_data = snap

    def undo(self):
        """Undo a single step. No-arg so it can wire directly to a QAction."""
        self.undo_steps(1)

    def redo(self):
        """Redo a single step. No-arg so it can wire directly to a QAction."""
        self.redo_steps(1)

    def undo_steps(self, n):
        """Undo ``n`` steps at once (history-dropdown 'undo to here'), emitting one
        change. Clamped to the available depth."""
        n = min(int(n), len(self._undo))
        if n <= 0:
            return
        for _ in range(n):
            self._undo_one()
        self._refresh_dirty()
        self.changed.emit()

    def redo_steps(self, n):
        """Redo ``n`` steps at once, emitting one change. Clamped to the available depth."""
        n = min(int(n), len(self._redo))
        if n <= 0:
            return
        for _ in range(n):
            self._redo_one()
        self._refresh_dirty()
        self.changed.emit()

    # --- save (scaffolded; wired into the UI in Phase 2) -----------------
    def save(self, path=None, template=None):
        target = str(path) if path else self.path
        if target is None:
            raise ValueError("No path to save to (use Save As).")
        save_slope_data_to_xlsx(self.slope_data, target, template=template)
        self.path = target
        self._write_styles_sidecar(target)      # {stem}_styles.json (or remove if none)
        self._clean_index = len(self._undo)     # this undo position is now saved
        self._style_dirty = False
        self._set_dirty(False)
