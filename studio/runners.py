"""Background run controllers (Phase 3).

Long engine calls run on a ``QThread`` so the UI stays responsive. ``LemRunner``
performs an LEM analysis — single surface or auto-search, circular or
non-circular — and emits the result (or an error) back to the GUI thread. Engine
``print`` output streams to the Log pane live via the thread-safe stdout tee
installed by the main window (``_LogStream``), so the worker does not touch any Qt
widget itself.
"""

from __future__ import annotations

import threading
import traceback

from PySide6.QtCore import QObject, QThread, Signal

# The engine modules are imported inside each worker (they are heavy and some pull
# in optional dependencies), but the preflight error type is needed in the EXCEPT
# clauses themselves, which run outside those imports. xslope.preflight is pure
# stdlib and costs nothing to import here.
from xslope.preflight import PreflightError

#: Element rows a thin material zone is meshed with when *Refine thin zones* is on.
#: The mesher's own definition of "thin" is three rows across the local width
#: (``mesh._REFINE_THIN_MIN_ELEMS``); four is the sizing that reproduced the
#: historical answer on the Griffiths thin-band case, so the toggle asks for four
#: rather than for the bare minimum it is trying to clear.
THIN_ZONE_ELEMS = 4

#: The toggle's default, ON. An under-resolved thin zone does not fail — it returns
#: a factor of safety that is simply too high — and a section with no thin zone
#: pays nothing for the option, because the detector finds no zone and no size is
#: added. Measured on the Griffiths thin-band section at a 3-unit target: the band
#: carries 1.6 element rows on a default triangular mesh and 1.2 on a quad one,
#: against 4.7 and 4.1 with this on.
REFINE_THIN_ZONES_DEFAULT = True


def thin_zone_size_regions(polygons, target_size, elems=THIN_ZONE_ELEMS):
    """Local mesh-size regions giving each thin material zone ``elems`` element rows.

    This is the QUADRILATERAL half of *Refine thin zones*. The two element families
    need different mechanisms, and the reason is in the mesher: a triangular mesh
    pins a zone's boundary with transfinite curve constraints, so a local size alone
    cannot resolve across a band (measured: 2.3 rows where 4 were asked for) and the
    ``thin_zones`` refine-feature — which also unpins those boundary edges — is what
    works. A quad mesh sets no boundary constraints at all, so the refine-feature's
    size field reaches only ~1.8-3.3 rows depending on the target, while a declared
    local size delivers the requested 4.1 rows at any target. Hence: triangles take
    the refine-features path, quads take this one.

    Returns ``[{'polygon': ring, 'size': float}, ...]`` in the shape
    ``build_mesh_from_polygons(size_regions=...)`` accepts, ready to be appended to
    the model's own ``Type='refine'`` overlays. A zone whose polygon already declares
    a Size is SKIPPED: an explicit Size is the user saying how finely that zone is to
    be meshed, and an automatic option must not overrule it. So is a zone already
    resolved at the global target, whose derived size would be inert anyway.
    """
    from xslope.mesh import detect_thin_zones

    coords, declared = [], []
    for p in polygons:
        if isinstance(p, dict):
            coords.append(list(p.get("coords") or []))
            declared.append(p.get("size"))
        else:
            coords.append(list(p))
            declared.append(None)

    out = []
    for zone in detect_thin_zones(coords, target_size):
        size = float(zone["width"]) / float(elems)
        if not (size > 0.0) or size >= target_size:
            continue                    # already resolved at the global size
        if zone["kind"] == "whole":
            i = zone["poly_index"]
            if not (0 <= i < len(coords)) or declared[i] is not None:
                continue                # off the end, or the model states its own Size
            ring = list(coords[i])
        else:
            xmin, ymin, xmax, ymax = zone["bbox"]
            ring = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
        out.append({"polygon": ring, "size": size})
    return out


class MeshWorker(QObject):
    """Builds finite-element meshes on a single, long-lived thread.

    gmsh is a global C singleton with thread affinity: initializing it on a fresh
    OS thread for each build (a new ``QThread`` per call) segfaults on the second
    build. So this worker is moved onto one persistent thread and its ``build``
    slot is invoked there for every mesh, keeping gmsh on one consistent thread.
    Includes reinforcement/pile constraint lines so the mesh also serves FEM.
    """

    succeeded = Signal(object)
    failed = Signal(str)

    def build(self, slope_data, options):
        from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                                 extract_constraint_line_geometry, extract_size_regions,
                                 MeshInputError, _REFINE_FEATURES)
        try:
            sd = slope_data
            element_type = options["element_type"]
            constraint_lines, n_reinf, n_pile = extract_constraint_line_geometry(sd)
            polygons = get_material_polygons(sd, reinf_lines=constraint_lines)
            if options.get("auto_size", True):
                divisions = options.get("size_divisions", 100)
                xs = [x for x, _ in sd["ground_surface"].coords]
                target = (max(xs) - min(xs)) / divisions
                print(f"Auto element size: {target:.3f} (slope width / {divisions} divisions)")
            else:
                target = options.get("target_size", 1.0)
            extra = (f", {n_reinf} reinforcement + {n_pile} pile line(s)"
                     if (n_reinf + n_pile) else "")
            wants_quads = element_type.startswith("quad")
            size_regions = extract_size_regions(sd)
            factor = float(options.get("refine_factor", 3.0))
            features = set(_REFINE_FEATURES) if options.get(
                "refine_near_features", False) else set()
            refine_msg = (f", refine x{factor:g} near features" if features else "")
            # Refine thin zones — a guarantee that a band too thin to mesh does not
            # quietly come back with a factor of safety that is too high. It only
            # ever ADDS: unchecked, the feature refinement above is left exactly as
            # the user set it. The mechanism differs by element family; see
            # thin_zone_size_regions for the measurements behind the split.
            thin_msg = ""
            if options.get("refine_thin_zones", REFINE_THIN_ZONES_DEFAULT):
                if wants_quads:
                    derived = thin_zone_size_regions(polygons, target)
                    size_regions = list(size_regions) + derived
                    if derived:
                        thin_msg = (f", {len(derived)} thin zone(s) sized for "
                                    f"{THIN_ZONE_ELEMS} element rows")
                else:
                    features.add("thin_zones")
                    thin_msg = ", thin zones refined"
            refine = factor if features else None
            refine_features = sorted(features) or None
            # Quadrilateral style is a per-run dialog choice, not a model input; it
            # is inert on the triangular element types, so it is only announced
            # where it can act.
            quad_style = options.get("quad_style", "free")
            style_msg = (f", {quad_style} quads"
                         if element_type.startswith("quad") else "")
            print(f"Building {element_type} mesh, target size {target:.3g}"
                  f"{style_msg}{extra}{refine_msg}{thin_msg}…")
            mesh = build_mesh_from_polygons(polygons, target_size=target,
                                            element_type=element_type,
                                            lines=constraint_lines or None,
                                            refine_factor=refine,
                                            refine_features=refine_features,
                                            size_regions=size_regions,
                                            quad_style=quad_style)
            n1d = len(mesh.get("elements_1d", []))
            print(f"Mesh built: {len(mesh['nodes'])} nodes, {len(mesh['elements'])} "
                  f"elements" + (f", {n1d} 1D elements" if n1d else "") + ".")
            self.succeeded.emit(mesh)
        except MeshInputError as e:
            print(f"Mesh input error: {e}")
            self.failed.emit(str(e))
        except Exception:
            traceback.print_exc()
            self.failed.emit("Mesh build failed — see the Log pane for details.")


class SeepRunner(QThread):
    """Runs a seepage solve off the GUI thread (no gmsh, so a plain per-run
    QThread is fine). ``options['bc']`` may be 1, 2, or ``'both'``; for ``'both'``
    it solves BC set 1 then 2, emitting ``succeeded`` once per set (each bundle's
    ``options['bc']`` is the concrete set). Emits ``failed`` if a set errors, but
    a later set still runs.

    When ``options['mode'] == 'transient'`` it instead runs the time-dependent
    solver and emits ONE bundle carrying the frame sequence:
    ``{mode: 'transient', seep_data, transient, frames, options}`` where ``frames``
    are the plottable per-frame solution dicts the play bar renders. The transient
    march reports determinate progress (``progress`` — the simulated-time fraction)
    and is cancellable (``cancel`` sets a flag the solver's progress callback reads;
    a cancelled march emits ``cancelled`` and stores no partial result)."""

    succeeded = Signal(object)
    failed = Signal(str)
    cancelled = Signal()
    progress = Signal(int, int, str)   # done, total (-1 = indeterminate), label

    def __init__(self, slope_data, options, parent=None):
        super().__init__(parent)
        self._sd = slope_data
        self._options = options
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    def run(self):
        from xslope.seep import build_seep_data, run_seepage_analysis, SeepInputError
        sd = self._sd
        mesh = sd.get("mesh")
        if mesh is None:
            self.failed.emit("No mesh available — build a mesh first.")
            return
        if self._options.get("mode") == "transient":
            self._run_transient(sd, mesh)
            return
        tol = self._options.get("tol", 1e-4)
        bc_opt = self._options.get("bc", 1)
        bcs = [1, 2] if bc_opt == "both" else [bc_opt]
        errors = []
        for bc in bcs:
            label = f"BC set {bc}" if len(bcs) > 1 else "Seepage"
            try:
                print(f"Building seepage data (BC set {bc})…")
                seep_data = build_seep_data(mesh, sd, seep_bc=bc)
                print(f"Running seepage analysis (BC set {bc}, tol={tol:g})…")
                solution = run_seepage_analysis(seep_data, tol=tol)
                if solution is None:            # defensive: no solver should return None now
                    raise RuntimeError("Seepage analysis returned no solution.")
                print(f"Seepage analysis complete (BC set {bc}).")
                self.succeeded.emit({"seep_data": seep_data, "solution": solution,
                                     "options": {**self._options, "bc": bc}})
            except (SeepInputError, PreflightError) as e:
                # Expected, user-actionable input problem — show the message, no
                # scary traceback needed (still logged for the record). A
                # PreflightError names the sheet and the field itself, so it belongs
                # in the same class as the seepage builder's own refusals.
                print(f"{label}: {e}")
                errors.append(f"{label}: {e}")
            except Exception as e:
                traceback.print_exc()
                errors.append(f"{label}: {e}  (see the Log pane for details.)")
        if errors:
            self.failed.emit("\n\n".join(errors))

    def _run_transient(self, sd, mesh):
        """Transient (time-dependent) seepage: build seep + tseep data, run the
        theta-stepper, and reconstruct the plottable per-frame solution dicts the
        play bar renders (velocity/gradient derived on load, exactly like the sidecar
        reload path). Stage times are model inputs on ``sd['tseep']`` (edited under
        Inputs → Transient), so ``build_tseep_data`` picks them up here.

        The march reports determinate progress (the simulated-time fraction) and is
        cancellable: a ``progress_callback(t, duration)`` emits ``progress`` and
        returns False once ``cancel`` is requested, which stops the solver cleanly and
        returns a ``cancelled`` result — we emit ``cancelled`` and store nothing.

        ``options['extra_save_times']`` adds instants to the save schedule for this
        march only, without touching the model. That is how a stability run asks for
        a time the previous solution never saved: the instant is COMPUTED here rather
        than interpolated between frames afterwards."""
        from xslope.seep import (build_seep_data, build_tseep_data,
                                  run_transient_seepage, _transient_frame_solution,
                                  SeepInputError)
        try:
            if sd.get("tseep") is None:
                raise SeepInputError("Transient run needs a 'tseep' sheet; this file "
                                     "has none.")
            print("Building seepage data (transient, BC set 1)…")
            seep_data = build_seep_data(mesh, sd, seep_bc=1)
            tseep_data = build_tseep_data(sd)
            extra = [float(t) for t in (self._options.get("extra_save_times") or [])]
            if extra:
                tseep_data = dict(tseep_data)
                tseep_data["save_times"] = sorted(
                    set(tseep_data.get("save_times") or []) | set(extra))
                print("Re-marching with extra save time(s): "
                      + ", ".join(f"{t:g}" for t in extra))
            print("Running transient seepage analysis…")
            time_unit = sd.get("time_unit")

            def _progress(t, duration):
                # Determinate progress = simulated-time fraction; the text carries the
                # declared time unit (bare numbers when undeclared). Returning False
                # (cancel requested) tells the solver to stop cleanly.
                frac = 0.0 if duration <= 0 else max(0.0, min(1.0, t / duration))
                label = (f"Transient seepage: t = {t:g} / {duration:g}"
                         + (f" {time_unit}" if time_unit else ""))
                self.progress.emit(int(round(frac * 100)), 100, label)
                return not self._cancel.is_set()

            solution = run_transient_seepage(seep_data, tseep_data, verbose=True,
                                             progress_callback=_progress)
            if solution.get("cancelled"):
                print("Transient seepage cancelled.")
                self.cancelled.emit()
                return
            unconfined = bool(solution.get("unconfined", False))
            # Reconstruct plottable frames (head/u/phi from the run; velocity and
            # gradient derived) so each renders through the unchanged plot_seep_solution.
            frames = [
                _transient_frame_solution(
                    seep_data, fr["head"], fr["u"], fr.get("phi"),
                    fr.get("inflow"), fr.get("outflow"), unconfined, time=fr["time"])
                for fr in solution["frames"]
            ]
            print(f"Transient seepage complete — {len(frames)} saved frame(s).")
            self.succeeded.emit({"mode": "transient", "seep_data": seep_data,
                                 "transient": solution, "frames": frames,
                                 "options": self._options})
        except (SeepInputError, PreflightError) as e:
            print(f"Transient seepage: {e}")
            self.failed.emit(f"Transient seepage: {e}")
        except Exception as e:
            traceback.print_exc()
            self.failed.emit(f"Transient seepage: {e}  (see the Log pane for details.)")


class FemRunner(QThread):
    """Runs an FEM analysis (single trial or SSRM) off the GUI thread. SSRM
    supports cooperative cancellation via a cancel_check threaded into solve_ssrm.
    Emits ``succeeded`` with ``{fem_data, solution, FS, analysis}``, ``failed``,
    or ``cancelled``."""

    succeeded = Signal(object)
    failed = Signal(str)
    cancelled = Signal()
    progress = Signal(int, int, str)   # done, total (-1 = indeterminate), label

    def __init__(self, slope_data, options, parent=None):
        super().__init__(parent)
        self._sd = slope_data
        self._options = options
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    def run(self):
        from xslope.fem import build_fem_data, solve_fem, solve_ssrm
        from xslope.search import AnalysisCancelled
        try:
            sd = self._sd
            mesh = sd.get("mesh")
            if mesh is None:
                self.failed.emit("No mesh available — build a mesh first.")
                return
            opts = self._options
            # The side boundary condition (v21) is read off slope_data inside
            # build_fem_data rather than passed as a solver argument like k0 /
            # tension_srf, so the dialog's choice is layered onto a SHALLOW COPY here.
            # Copying keeps a run option from silently rewriting the open document —
            # the dialog seeds itself from the file, and the run must not edit it back.
            if opts.get("side_bc"):
                sd = dict(sd)
                sd["side_bc"] = opts["side_bc"]
            print("Building FEM data…")
            fem_data = build_fem_data(sd, mesh)
            analysis = opts.get("analysis", "ssrm")
            if analysis == "single":
                F = opts.get("F", 1.0)
                print(f"Solving FEM (single trial, F={F:g})…")

                def fem_cb(frac, label):
                    self.progress.emit(int(frac * 100), 100, str(label))

                solution = solve_fem(fem_data, F=F, debug_level=1,
                                     k0=opts.get("k0"),
                                     progress_callback=fem_cb)
                print(f"FEM solve: converged={solution.get('converged')}, "
                      f"iterations={solution.get('iterations')}")
                self.succeeded.emit({"fem_data": fem_data, "solution": solution,
                                     "FS": None, "analysis": "single"})
            else:
                print(f"Running SSRM (F in [{opts.get('F_min', 1.0):g}, "
                      f"{opts.get('F_max', 2.0):g}], {opts.get('failure_criterion')})…")
                def cb(done, total, label):
                    self.progress.emit(int(done), int(total) if total else -1, str(label))

                result = solve_ssrm(
                    fem_data, F_min=opts.get("F_min", 1.0), F_max=opts.get("F_max", 2.0),
                    tolerance=opts.get("tolerance", 0.01), debug_level=1,
                    failure_criterion=opts.get("failure_criterion", "non_convergence"),
                    min_slip_depth=opts.get("min_slip_depth"),
                    ssr_exclude=opts.get("ssr_exclude"),
                    # v19 run options. The dialog seeds both from the file, so what
                    # arrives here is the user's confirmed choice and outranks the
                    # file value solve_ssrm would otherwise fall back to.
                    k0=opts.get("k0"),
                    tension_srf=opts.get("tension_srf"),
                    capture_failure_state=opts.get("capture_failure_state", True),
                    capture_margin=opts.get("capture_margin", 0.15),
                    capture_max_iterations=opts.get("capture_max_iterations"),
                    cancel_check=self._cancel.is_set, progress_callback=cb)
                if not result.get("converged", False):
                    self.failed.emit(f"SSRM did not converge: "
                                     f"{result.get('error', 'unknown error')}")
                    return
                fs = result.get("FS")
                print(f"SSRM factor of safety = {fs:.3f}")
                self.succeeded.emit({"fem_data": fem_data,
                                     "solution": result["last_solution"],
                                     # The at-failure (unconverged) mechanism field the
                                     # deformation/vector panels render (None if absent).
                                     "failure_solution": result.get("failure_solution"),
                                     "FS": fs, "analysis": "ssrm"})
        except AnalysisCancelled:
            print("Run cancelled.")
            self.cancelled.emit()
        except PreflightError as e:
            # An input problem the registry can name. The message already says which
            # sheet and which field, so it reaches the box verbatim instead of being
            # flattened into "see the Log pane" -- the same treatment SeepInputError
            # has always had. Reachable when the model changed after the Run dialog
            # gated it, and on the assistant's run_python path, which has no dialog.
            print(str(e))
            self.failed.emit(str(e))
        except Exception:
            traceback.print_exc()
            self.failed.emit("FEM run failed — see the Log pane for details.")


class LemRunner(QThread):
    """Runs an LEM analysis off the GUI thread.

    Emits ``succeeded`` with a bundle ``{slice_df, failure_surface, results,
    search}`` (``search`` is ``None`` for single-surface, else a dict describing
    the search for the Search view), ``failed`` with an error message, or
    ``cancelled`` if the run was aborted via :meth:`cancel`.
    """

    succeeded = Signal(object)
    failed = Signal(str)
    cancelled = Signal()
    progress = Signal(int, int, str)   # done, total (-1 = indeterminate), label

    def __init__(self, slope_data, options, parent=None):
        super().__init__(parent)
        self._sd = slope_data
        self._method = options["method"]
        self._analysis = options.get("analysis", "single_surface")
        self._surface = options.get("surface", "circular")
        self._num_slices = options.get("num_slices", 40)
        self._rapid = options.get("rapid", False)
        self._composite = options.get("composite", False)
        self._seed = 'grid' if options.get("grid_seed", False) else 'circles'
        self._diagnostic = options.get("diagnostic", False)
        self._fs_tol = options.get("fs_tol")
        self._tol = options.get("tol")
        self._max_iter = options.get("max_iter")
        self._min_slip_depth = options.get("min_slip_depth")
        self._cancel = threading.Event()

    def cancel(self):
        """Request cooperative cancellation; the search loops stop at the next
        iteration boundary and the run emits ``cancelled``."""
        self._cancel.set()

    def _search_kwargs(self, circular):
        """Tolerance and search-window kwargs to forward to the search functions.
        ``tol`` and the window limits are only accepted by ``circular_search``
        (noncircular has neither a geometric tol nor entry/exit/centre limits)."""
        kw = {}
        if self._fs_tol is not None:
            kw["fs_tol"] = self._fs_tol
        if self._max_iter is not None:
            kw["max_iter"] = self._max_iter
        if self._min_slip_depth is not None:
            kw["min_slip_depth"] = self._min_slip_depth
        if circular and self._tol is not None:
            kw["tol"] = self._tol
        if circular:
            kw.update(self._file_search_window(kw))
        return kw

    def _file_search_window(self, already):
        """The circles-sheet search window (v19) as circular_search kwargs.

        Each limit is forwarded only when the file supplies BOTH ends of it — a
        half-declared range is not a window, and inventing the missing end would
        silently constrain the search in a direction nobody asked for. Anything the
        dialog already set (``already``) wins and is left alone.
        """
        sw = (self._sd.get("search_window") or {}) if isinstance(self._sd, dict) else {}
        if not sw:
            return {}
        kw = {}

        def pair(name, lo, hi):
            if name in already:
                return
            if sw.get(lo) is not None and sw.get(hi) is not None:
                kw[name] = (float(sw[lo]), float(sw[hi]))

        pair("entry_range", "entry_x_min", "entry_x_max")
        pair("exit_range", "exit_x_min", "exit_x_max")
        if ("center_box" not in already
                and all(sw.get(k) is not None for k in
                        ("center_box_x_min", "center_box_y_min",
                         "center_box_x_max", "center_box_y_max"))):
            kw["center_box"] = (float(sw["center_box_x_min"]), float(sw["center_box_y_min"]),
                                float(sw["center_box_x_max"]), float(sw["center_box_y_max"]))
        # "Max tangent depth" is a single ELEVATION: the deepest the circle's tangent
        # point may reach. circular_search takes a band, so pair it with the top of the
        # model — the only upper end that adds no constraint of its own.
        if "tangent_depth" not in already and sw.get("max_tangent_depth") is not None:
            gs = self._sd.get("ground_surface")
            try:
                y_top = max(y for _x, y in gs.coords)
            except Exception:
                y_top = None
            if y_top is not None and y_top > float(sw["max_tangent_depth"]):
                kw["tangent_depth"] = (float(sw["max_tangent_depth"]), float(y_top))
        if "min_slip_depth" not in already and sw.get("min_slip_depth") is not None:
            kw["min_slip_depth"] = float(sw["min_slip_depth"])
        if kw:
            print(f"Applying the file's search window: "
                  f"{', '.join(sorted(kw))}.")
        return kw

    def run(self):
        from xslope.search import AnalysisCancelled
        try:
            if self._analysis == "auto_search":
                self._run_search()
            else:
                self._run_single()
        except AnalysisCancelled:
            print("Run cancelled.")
            self.cancelled.emit()
        except PreflightError as e:
            # An input problem the registry can name. The message already says which
            # sheet and which field, so it reaches the box verbatim instead of being
            # flattened into "see the Log pane" -- the same treatment SeepInputError
            # has always had. Reachable when the model changed after the Run dialog
            # gated it, and on the assistant's run_python path, which has no dialog.
            print(str(e))
            self.failed.emit(str(e))
        except Exception:
            traceback.print_exc()   # streams to the Log pane via the stdout/stderr tee
            self.failed.emit("Solve failed — see the Log pane for details.")

    # --- single surface --------------------------------------------------
    def _run_single(self):
        from xslope.slice import generate_slices
        from xslope.solve import solve_selected

        sd = self._sd
        circular = self._surface == "circular"
        if circular:
            circle = sd["circles"][0] if sd.get("circular") and sd.get("circles") else None
            if circle is None:
                self.failed.emit("A circular surface is required (no circles defined).")
                return
            non_circ = None
            print(f"Running {self._method.upper()} — single circular surface "
                  f"(Xo={circle.get('Xo')}, Yo={circle.get('Yo')}, R={circle.get('R'):.3g}), "
                  f"{self._num_slices} slices{self._rapid_tag()}…")
        else:
            non_circ = sd.get("non_circ") or None
            if not non_circ:
                self.failed.emit("A non-circular surface is required "
                                 "(no non-circular points defined).")
                return
            circle = None
            print(f"Running {self._method.upper()} — single non-circular surface, "
                  f"{self._num_slices} slices{self._rapid_tag()}…")

        ok, result = generate_slices(sd, circle=circle, non_circ=non_circ,
                                     num_slices=self._num_slices,
                                     composite=self._composite)
        if not ok:
            self.failed.emit(str(result))
            return
        slice_df, failure_surface = result
        print(f"Generated {len(slice_df)} slices; solving…")
        results = solve_selected(self._method, slice_df, rapid=self._rapid)
        if not isinstance(results, dict):
            self.failed.emit(f"No solution: {results}")
            return
        self.succeeded.emit({"slice_df": slice_df, "failure_surface": failure_surface,
                             "results": results, "search": None})

    # --- auto-search -----------------------------------------------------
    def _run_search(self):
        from xslope.search import circular_search, noncircular_search

        sd = self._sd
        circular = self._surface == "circular"
        if circular:
            if self._seed != 'grid' and not (sd.get("circular") and sd.get("circles")):
                self.failed.emit("Circular search needs at least one starting circle "
                                 "(or turn on grid seeding).")
                return
            print(f"Searching for the critical circular surface with "
                  f"{self._method.upper()}{self._rapid_tag()}…")
            fs_cache, converged, search_path, circle_cache = circular_search(
                sd, self._method, rapid=self._rapid, num_slices=self._num_slices,
                diagnostic=self._diagnostic, cancel_check=self._cancel.is_set,
                composite=self._composite, seed=self._seed,
                **self._search_kwargs(circular=True))
            search = {"kind": "circular", "fs_cache": fs_cache,
                      "search_path": search_path, "circle_cache": circle_cache}
        else:
            if not sd.get("non_circ"):
                self.failed.emit("Non-circular search needs a starting non-circular surface.")
                return
            print(f"Searching for the critical non-circular surface with "
                  f"{self._method.upper()}{self._rapid_tag()}…")
            fs_cache, converged, search_path = noncircular_search(
                sd, self._method, rapid=self._rapid, num_slices=self._num_slices,
                diagnostic=self._diagnostic, cancel_check=self._cancel.is_set,
                **self._search_kwargs(circular=False))
            search = {"kind": "noncircular", "fs_cache": fs_cache,
                      "search_path": search_path, "circle_cache": None}

        if not fs_cache:
            self.failed.emit("Search produced no valid surfaces.")
            return
        critical = fs_cache[0]
        results = critical.get("solver_result")
        if not isinstance(results, dict):
            self.failed.emit("Search found no surface with a valid solution.")
            return
        tail = "" if converged else "  (search did not fully converge)"
        print(f"Critical FS = {results.get('FS'):.3f}{tail}")
        self.succeeded.emit({"slice_df": critical.get("slices"),
                             "failure_surface": critical.get("failure_surface"),
                             "results": results, "search": search})

    def _rapid_tag(self):
        return " (rapid drawdown)" if self._rapid else ""


class SensitivityRunner(QThread):
    """Runs a Parametric study (sensitivity / design / back-analysis) off the GUI
    thread.

    A sweep is N sequential solves (no multiprocessing). Progress is reported per
    solve and the run cancels cooperatively (each engine sweep is passed a
    ``cancel_check``; an in-flight solve finishes, then the sweep stops at the next
    point). A thin caller of ``xslope.sensitivity``, keyed on ``opts['mode']`` and
    (for sensitivity) ``opts['plot_type']``:

      * design         -> ``design()``        -> bundle {'kind':'design', …}
      * back_analysis  -> ``back_analysis()``  -> bundle {'kind':'design',
                          'study':'back_analysis', …}
      * sensitivity, plot_type tornado/spider -> per-parameter ``sensitivity()``
                          sweeps (+ ``tornado_from_sweeps()`` for the tornado)
      * sensitivity, plot_type scaled   -> ``scaled_sensitivity()``
      * sensitivity, plot_type variance -> ``variance_contribution()``
      * sensitivity, plot_type rank     -> ``mc_rank_correlation()``
        -> bundle {'kind':'sensitivity', 'plot_type': …, <payload>, method}.
    """

    succeeded = Signal(object)
    failed = Signal(str)
    cancelled = Signal()
    progress = Signal(int, int, str)   # done, total, label

    def __init__(self, slope_data, options, parent=None):
        super().__init__(parent)
        self._sd = slope_data
        self._opts = options
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    def run(self):
        from xslope.search import AnalysisCancelled
        try:
            if self._opts.get("mode") in ("design", "back_analysis"):
                self._run_design()
            else:
                self._run_sensitivity()
        except AnalysisCancelled:
            print("Sweep cancelled.")
            self.cancelled.emit()
        except PreflightError as e:
            # An input problem the registry can name. The message already says which
            # sheet and which field, so it reaches the box verbatim instead of being
            # flattened into "see the Log pane" -- the same treatment SeepInputError
            # has always had. Reachable when the model changed after the Run dialog
            # gated it, and on the assistant's run_python path, which has no dialog.
            print(str(e))
            self.failed.emit(str(e))
        except Exception:
            traceback.print_exc()   # streams to the Log pane via the stdout tee
            self.failed.emit("Sweep failed — see the Log pane for details.")

    def _run_design(self):
        from xslope.sensitivity import design, back_analysis
        o = self._opts
        emode = o.get("engine_mode", "lem")
        is_back = o.get("mode") == "back_analysis"
        study = "Back-analysis" if is_back else "Design sweep"
        fn = back_analysis if is_back else design
        total = int(o["steps"]) + 1        # points + base

        def cb(done, _t, label):
            self.progress.emit(int(done), total, str(label))

        out = "q" if emode == "seep" else "FS"
        engine_tag = {"lem": o.get("method", "spencer"), "fem": "SSRM",
                      "seep": f"BC {o.get('seep_opts', {}).get('bc', 1)}"}.get(emode, emode)
        print(f"{study} ({emode}): {o['param']} from {o['low']:g} to "
              f"{o['high']:g} in {o['steps']} steps, target {out} = "
              f"{o['target_fs']:g} ({engine_tag})…")
        ok, res = fn(self._sd, param=o["param"], low=o["low"], high=o["high"],
                     steps=o["steps"], target_fs=o["target_fs"], mode=emode,
                     method=o.get("method", "spencer"),
                     search=o.get("search", True),
                     num_slices=o.get("num_slices", 40),
                     fem_opts=o.get("fem_opts"), seep_opts=o.get("seep_opts"),
                     progress_callback=cb, cancel_check=self._cancel.is_set)
        if not ok:
            self.failed.emit(str(res))
            return
        print(res["message"])
        # kind stays 'design' — the result canvas renders both the same way; the
        # 'study' field (set by back_analysis) drives the title/status wording.
        self.succeeded.emit({"kind": "design",
                             "study": res.get("study", "design"), **res})

    def _run_sensitivity(self):
        o = self._opts
        plot_type = o.get("plot_type", "tornado")
        # The variance Pareto and MC rank plots are their own engine runs (they use
        # every sigma-carrying material, not the parameter table).
        if plot_type == "variance":
            self._run_variance()
            return
        if plot_type == "rank":
            self._run_rank()
            return
        self._run_local(plot_type)

    def _run_local(self, plot_type):
        """Tornado / spider (full per-parameter sweeps) or scaled bars (central
        differences at +/-1%)."""
        import numpy as np
        from xslope.sensitivity import (sensitivity, tornado_from_sweeps,
                                        scaled_sensitivity)
        o = self._opts
        emode = o.get("engine_mode", "lem")
        specs = o["params"]
        if not specs:
            self.failed.emit("Add at least one parameter to sweep.")
            return
        n = int(o["n"])
        method = o.get("method", "spencer")
        disp_method = {"lem": method, "fem": "SSRM", "seep": ""}.get(emode, method)

        if plot_type == "scaled":
            count = [0]

            def cb(done, total, label):
                count[0] += 1
                self.progress.emit(count[0], max(len(specs) * 3, 1), str(label))

            print(f"Scaled sensitivity ({emode}): {len(specs)} parameter(s), "
                  f"central differences at ±1% ({o.get('scaling', 'elasticity')})…")
            ok, res = scaled_sensitivity(
                self._sd, [s["ref"] for s in specs], method=method, mode=emode,
                search=o.get("search", True), num_slices=o.get("num_slices", 40),
                fem_opts=o.get("fem_opts"), seep_opts=o.get("seep_opts"),
                progress_callback=cb, cancel_check=self._cancel.is_set)
            if not ok:
                self.failed.emit(str(res))
                return
            print(f"Scaled sensitivity done — {len(res['bars'])} parameter(s).")
            self.succeeded.emit({"kind": "sensitivity", "plot_type": "scaled",
                                 "scaled": res, "scaling": o.get("scaling",
                                                                 "elasticity"),
                                 "method": disp_method})
            return

        # tornado / spider: full per-parameter FS-vs-value sweeps
        def n_points(spec):
            return len(spec["values"]) if spec.get("values") is not None else n

        total = sum(n_points(s) + 1 for s in specs)
        count = [0]

        def cb(_done, _t, label):
            count[0] += 1
            self.progress.emit(count[0], total, str(label))

        sweeps, base_fs, out = {}, None, "FS"
        for i, spec in enumerate(specs):
            ref = spec["ref"]
            if spec.get("low") is not None and spec.get("high") is not None:
                values = list(np.linspace(spec["low"], spec["high"], n))
                rel, tag = 0.5, f"±σ [{spec['low']:g}, {spec['high']:g}]"
            else:
                values, rel = None, spec.get("rel_range", 0.2)
                tag = f"±{rel * 100:g}%"
            print(f"[{i + 1}/{len(specs)}] Sweeping {ref} {tag} ({emode})…")
            ok, res = sensitivity(self._sd, param=ref, values=values, rel_range=rel,
                                  n=n, mode=emode, methods=(method,),
                                  search=o.get("search", True),
                                  num_slices=o.get("num_slices", 40),
                                  fem_opts=o.get("fem_opts"),
                                  seep_opts=o.get("seep_opts"),
                                  progress_callback=cb,
                                  cancel_check=self._cancel.is_set)
            if not ok:
                self.failed.emit(f"{ref}: {res}")
                return
            df = res["df"]
            out = res.get("output", "FS")
            sweeps[res["param"]] = df
            if base_fs is None:
                b = df.loc[df["is_base"] & df["success"], "fs"]
                base_fs = float(b.iloc[0]) if len(b) else None
        tornado = tornado_from_sweeps(sweeps, base_fs=base_fs, method=disp_method)
        print(f"Sensitivity done — {len(sweeps)} parameter(s), base {out} = "
              f"{base_fs:.3g}." if base_fs is not None else "Sensitivity done.")
        self.succeeded.emit({"kind": "sensitivity", "plot_type": plot_type,
                             "sweeps": sweeps, "tornado": tornado,
                             "method": disp_method})

    def _run_variance(self):
        from xslope.sensitivity import variance_contribution
        o = self._opts
        method = o.get("method", "spencer")

        def cb(done, total, label):
            self.progress.emit(int(done or 0), int(total or 0), str(label))

        print(f"Variance contribution ({method}): Taylor-series reliability over "
              f"every σ-carrying material…")
        ok, res = variance_contribution(self._sd, method=method,
                                        search=o.get("search", True),
                                        progress_callback=cb,
                                        cancel_check=self._cancel.is_set)
        if not ok:
            self.failed.emit(str(res))
            return
        print(f"Variance contribution done — {len(res['bars'])} parameter(s), "
              f"σ_F = {res.get('sigma_F'):.3g}.")
        self.succeeded.emit({"kind": "sensitivity", "plot_type": "variance",
                             "variance": res, "method": method})

    def _run_rank(self):
        from xslope.sensitivity import mc_rank_correlation
        o = self._opts
        method = o.get("method", "bishop")
        n_samples = int(o.get("mc_samples", 5000))

        def cb(done, total, label):
            self.progress.emit(int(done or 0), int(total or 0), str(label))

        print(f"Monte Carlo rank correlation ({method}): {n_samples} samples…")
        ok, res = mc_rank_correlation(self._sd, method=method, n_samples=n_samples,
                                      search=o.get("search", True),
                                      num_slices=o.get("num_slices", 40),
                                      progress_callback=cb,
                                      cancel_check=self._cancel.is_set)
        if not ok:
            self.failed.emit(str(res))
            return
        print(f"Monte Carlo rank correlation done — {res.get('n_valid')} valid "
              f"of {res.get('n_samples')} samples.")
        self.succeeded.emit({"kind": "sensitivity", "plot_type": "rank",
                             "rank": res, "method": method})


class ReliabilityRunner(QThread):
    """Runs a probabilistic reliability analysis off the GUI thread — the sibling
    of ``SensitivityRunner`` for the Reliability toolbar button.

    Dispatches on the dialog options:

      * app_mode 'lem', engine 'taylor' -> ``reliability(engine='taylor')`` (TSPM).
      * app_mode 'lem', engine 'mc'      -> ``reliability(engine='mc')`` (Monte Carlo).
      * app_mode 'fem'                   -> ``reliability_fem`` (Taylor / SSRM).

    Emits ``succeeded`` with a bundle ``{engine, app_mode, reliability, ...}``,
    ``failed`` with a message, ``cancelled``, or ``progress``. The engines stream
    their per-parameter tables to the Log pane via the stdout tee.
    """

    succeeded = Signal(object)
    failed = Signal(str)
    cancelled = Signal()
    progress = Signal(int, int, str)   # done, total (-1 = indeterminate), label

    def __init__(self, slope_data, options, parent=None):
        super().__init__(parent)
        self._sd = slope_data
        self._opts = options
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    def run(self):
        from xslope.search import AnalysisCancelled
        try:
            if self._opts.get("app_mode") == "fem":
                self._run_fem()
            elif self._opts.get("engine") == "mc":
                self._run_mc()
            else:
                self._run_taylor()
        except AnalysisCancelled:
            print("Reliability run cancelled.")
            self.cancelled.emit()
        except PreflightError as e:
            # An input problem the registry can name. The message already says which
            # sheet and which field, so it reaches the box verbatim instead of being
            # flattened into "see the Log pane" -- the same treatment SeepInputError
            # has always had. Reachable when the model changed after the Run dialog
            # gated it, and on the assistant's run_python path, which has no dialog.
            print(str(e))
            self.failed.emit(str(e))
        except Exception:
            traceback.print_exc()   # streams to the Log pane via the stdout tee
            self.failed.emit("Reliability run failed — see the Log pane for details.")

    def _progress_cb(self):
        def cb(done, total, label):
            self.progress.emit(int(done), int(total) if total else -1, str(label))
        return cb

    def _run_taylor(self):
        from xslope.reliability import reliability
        o, sd = self._opts, self._sd
        circular = o.get("surface", "circular") == "circular"
        if circular and not (sd.get("circular") and sd.get("circles")):
            self.failed.emit("Reliability (circular) needs at least one starting circle.")
            return
        if not circular and not sd.get("non_circ"):
            self.failed.emit("Reliability (non-circular) needs a starting "
                             "non-circular surface.")
            return
        print(f"Running Taylor-series reliability — {o['method'].upper()}, "
              f"{'circular' if circular else 'non-circular'}"
              f"{' (rapid drawdown)' if o.get('rapid') else ''}…")
        ok, result = reliability(
            sd, o["method"], engine="taylor", rapid=o.get("rapid", False),
            circular=circular, search=o.get("search", True), debug_level=0,
            progress_callback=self._progress_cb(), cancel_check=self._cancel.is_set)
        if not ok:
            self.failed.emit(str(result))
            return
        self.succeeded.emit({"engine": "taylor", "app_mode": "lem",
                             "reliability": result})

    def _run_mc(self):
        from xslope.reliability import reliability
        o, sd = self._opts, self._sd
        circular = o.get("surface", "circular") == "circular"
        print(f"Running Monte Carlo reliability — {o['method'].upper()}, "
              f"{o.get('n_samples')} samples, seed {o.get('rng_seed')}, "
              f"{o.get('distribution', 'normal')}…")
        ok, result = reliability(
            sd, o["method"], engine="mc", rapid=o.get("rapid", False),
            circular=circular, search=o.get("search", True),
            n_samples=int(o.get("n_samples", 10000)),
            rng_seed=int(o.get("rng_seed", 20240117)),
            distribution=o.get("distribution", "normal"),
            num_slices=int(o.get("num_slices", 40)), debug_level=0,
            progress_callback=self._progress_cb(), cancel_check=self._cancel.is_set)
        if not ok:
            self.failed.emit(str(result))
            return
        self.succeeded.emit({"engine": "mc", "app_mode": "lem",
                             "reliability": result})

    def _run_fem(self):
        from xslope.reliability import reliability_fem
        from xslope.fem import build_fem_data
        o, sd = self._opts, self._sd
        mesh = sd.get("mesh")
        if mesh is None:
            self.failed.emit("No mesh available — build a mesh first.")
            return
        print(f"Running FEM reliability (SSRM, F in "
              f"[{o.get('F_min', 0.7):g}, {o.get('F_max', 2.0):g}])…")
        ok, result = reliability_fem(
            sd, mesh=mesh, F_min=o.get("F_min", 0.7), F_max=o.get("F_max", 2.0),
            tolerance=o.get("reliability_tol", 0.001), debug_level=1,
            progress_callback=self._progress_cb(), cancel_check=self._cancel.is_set)
        if not ok:
            self.failed.emit(f"Reliability failed: {result}")
            return
        fem_data = build_fem_data(sd, mesh)
        self.succeeded.emit({
            "engine": "taylor", "app_mode": "fem", "reliability": result,
            "fem_data": fem_data,
            "solution": result["mlv_solution"]["last_solution"],
            # The at-failure (unconverged) mechanism field the deformation/vector
            # panels render (None if absent) — mirrors the normal FEM SSRM path above.
            "failure_solution": result["mlv_solution"].get("failure_solution"),
            "FS": result["F_MLV"]})


class UpdateCheckRunner(QThread):
    """Reads the release manifest off the GUI thread.

    A version check is a network round trip, so it gets the same treatment as a
    solve: its own QThread, one signal back to the GUI thread, and nothing Qt
    touched inside :meth:`run`. It has no ``failed`` signal on purpose — a check
    that could not reach the network still produces an answer
    (``UpdateInfo(status="error")``), and whether that answer is shown is the
    caller's decision, not the worker's.
    """

    checked = Signal(object)          # studio.updater.UpdateInfo

    def __init__(self, current_version=None, url=None, timeout=None, parent=None):
        super().__init__(parent)
        self._version = current_version
        self._url = url
        self._timeout = timeout

    def run(self):
        from .updater import TIMEOUT, check_for_update
        self.checked.emit(check_for_update(
            current_version=self._version, url=self._url,
            timeout=self._timeout if self._timeout is not None else TIMEOUT))


class UpdateDownloadRunner(QThread):
    """Downloads one release artifact and verifies its sha256, cancellably.

    Emits ``progress`` as a percentage (0-100) with a byte-count label, exactly
    like the transient seepage march, so the status-bar bar and Cancel button
    behave the same for an update as for a run. ``succeeded`` carries the path of
    a file whose digest MATCHED the manifest; a mismatch is a ``failed``, and the
    partial or corrupt file is gone by the time the signal arrives.
    """

    succeeded = Signal(str)           # path of the verified artifact
    failed = Signal(str)
    cancelled = Signal()
    progress = Signal(int, int, str)  # done, total (-1 = indeterminate), label

    def __init__(self, artifact, dest_dir, label="", parent=None):
        super().__init__(parent)
        self._artifact = artifact
        self._dest_dir = dest_dir
        self._label = label or artifact.get("filename", "update")
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    def run(self):
        from .updater import UpdateCancelled, UpdateError, download_artifact
        last = [-1]

        def _progress(done, total):
            if total > 0:
                pct = int(done * 100 / total)
                if pct == last[0]:
                    return                     # one emit per percent, not per chunk
                last[0] = pct
                self.progress.emit(pct, 100, f"Downloading {self._label} — "
                                             f"{done / 1e6:.1f} / {total / 1e6:.1f} MB")
            else:
                self.progress.emit(0, -1, f"Downloading {self._label} — "
                                          f"{done / 1e6:.1f} MB")

        try:
            print(f"Downloading {self._label}…")
            path = download_artifact(self._artifact, self._dest_dir,
                                     progress=_progress,
                                     cancel=self._cancel.is_set)
        except UpdateCancelled:
            print("Update download cancelled.")
            self.cancelled.emit()
            return
        except UpdateError as exc:
            print(f"Update download failed: {exc}")
            self.failed.emit(str(exc))
            return
        except Exception as exc:
            traceback.print_exc()
            self.failed.emit(f"Update download failed: {exc}")
            return
        print(f"Downloaded and checksum-verified: {path}")
        self.succeeded.emit(path)


class AssistantModelsRunner(QThread):
    """Reads a provider's model list (and the recommendations manifest) off the
    GUI thread.

    Same shape as :class:`UpdateCheckRunner`: two short network reads, one signal
    back, no Qt touched inside :meth:`run`, and no ``failed`` signal — a list that
    could not be fetched still produces an answer (``models=None`` plus the
    sentence explaining it), and the dialog falls back to its cached or shipped
    list rather than showing an error.

    Deliberately does no settings I/O: the cache and the manifest are written by
    the dialog, on the GUI thread, from the payload this emits.
    """

    listed = Signal(object)           # {"provider", "models", "error", "manifest"}

    def __init__(self, provider, api_key="", base_url=None, want_models=True,
                 want_manifest=False, models_url=None, manifest_url=None,
                 parent=None):
        super().__init__(parent)
        self._provider = provider
        self._api_key = api_key
        self._base_url = base_url
        self._want_models = want_models
        self._want_manifest = want_manifest
        self._models_url = models_url
        self._manifest_url = manifest_url

    def run(self):
        from .ai import models as model_list
        out = {"provider": self._provider, "models": None, "error": "",
               "manifest": None, "manifest_checked": bool(self._want_manifest)}
        if self._want_models:
            out["models"], out["error"] = model_list.fetch_models(
                self._provider, api_key=self._api_key, base_url=self._base_url,
                url=self._models_url)
        if self._want_manifest:
            try:
                out["manifest"] = model_list.fetch_manifest(url=self._manifest_url)
            except Exception:
                out["manifest"] = None     # silent: a manifest is advice, not data
        self.listed.emit(out)
