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
                                 extract_constraint_line_geometry)
        try:
            sd = slope_data
            element_type = options["element_type"]
            constraint_lines, n_reinf, n_pile = extract_constraint_line_geometry(sd)
            polygons = get_material_polygons(sd, reinf_lines=constraint_lines)
            if options.get("auto_size", True):
                divisions = options.get("size_divisions", 50)
                xs = [x for x, _ in sd["ground_surface"].coords]
                target = (max(xs) - min(xs)) / divisions
                print(f"Auto element size: {target:.3f} (slope width / {divisions} divisions)")
            else:
                target = options.get("target_size", 1.0)
            extra = (f", {n_reinf} reinforcement + {n_pile} pile line(s)"
                     if (n_reinf + n_pile) else "")
            print(f"Building {element_type} mesh, target size {target:.3g}{extra}…")
            mesh = build_mesh_from_polygons(polygons, target_size=target,
                                            element_type=element_type,
                                            lines=constraint_lines or None)
            n1d = len(mesh.get("elements_1d", []))
            print(f"Mesh built: {len(mesh['nodes'])} nodes, {len(mesh['elements'])} "
                  f"elements" + (f", {n1d} 1D elements" if n1d else "") + ".")
            self.succeeded.emit(mesh)
        except Exception:
            traceback.print_exc()
            self.failed.emit("Mesh build failed — see the Log pane for details.")


class SeepRunner(QThread):
    """Runs a seepage solve off the GUI thread (no gmsh, so a plain per-run
    QThread is fine). Emits ``succeeded`` with ``{seep_data, solution, options}``
    or ``failed``."""

    succeeded = Signal(object)
    failed = Signal(str)

    def __init__(self, slope_data, options, parent=None):
        super().__init__(parent)
        self._sd = slope_data
        self._options = options

    def run(self):
        from xslope.seep import build_seep_data, run_seepage_analysis
        try:
            sd = self._sd
            mesh = sd.get("mesh")
            if mesh is None:
                self.failed.emit("No mesh available — build a mesh first.")
                return
            bc = self._options.get("bc", 1)
            tol = self._options.get("tol", 1e-4)
            print(f"Building seepage data (BC set {bc})…")
            seep_data = build_seep_data(mesh, sd, seep_bc=bc)
            print(f"Running seepage analysis (tol={tol:g})…")
            solution = run_seepage_analysis(seep_data, tol=tol)
            print("Seepage analysis complete.")
            self.succeeded.emit({"seep_data": seep_data, "solution": solution,
                                 "options": self._options})
        except Exception:
            traceback.print_exc()
            self.failed.emit("Seepage run failed — see the Log pane for details.")


class FemRunner(QThread):
    """Runs an FEM analysis (single trial or SSRM) off the GUI thread. SSRM
    supports cooperative cancellation via a cancel_check threaded into solve_ssrm.
    Emits ``succeeded`` with ``{fem_data, solution, FS, analysis}``, ``failed``,
    or ``cancelled``."""

    succeeded = Signal(object)
    failed = Signal(str)
    cancelled = Signal()

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
            print("Building FEM data…")
            fem_data = build_fem_data(sd, mesh)
            opts = self._options
            if opts.get("analysis", "ssrm") == "single":
                F = opts.get("F", 1.0)
                print(f"Solving FEM (single trial, F={F:g})…")
                solution = solve_fem(fem_data, F=F, debug_level=1)
                print(f"FEM solve: converged={solution.get('converged')}, "
                      f"iterations={solution.get('iterations')}")
                self.succeeded.emit({"fem_data": fem_data, "solution": solution,
                                     "FS": None, "analysis": "single"})
            else:
                print(f"Running SSRM (F in [{opts.get('F_min', 1.0):g}, "
                      f"{opts.get('F_max', 2.0):g}], {opts.get('failure_criterion')})…")
                result = solve_ssrm(
                    fem_data, F_min=opts.get("F_min", 1.0), F_max=opts.get("F_max", 2.0),
                    tolerance=opts.get("tolerance", 0.01), debug_level=1,
                    failure_criterion=opts.get("failure_criterion", "non_convergence"),
                    cancel_check=self._cancel.is_set)
                if not result.get("converged", False):
                    self.failed.emit(f"SSRM did not converge: "
                                     f"{result.get('error', 'unknown error')}")
                    return
                fs = result.get("FS")
                print(f"SSRM factor of safety = {fs:.3f}")
                self.succeeded.emit({"fem_data": fem_data,
                                     "solution": result["last_solution"],
                                     "FS": fs, "analysis": "ssrm"})
        except AnalysisCancelled:
            print("Run cancelled.")
            self.cancelled.emit()
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
        self._diagnostic = options.get("diagnostic", False)
        self._cancel = threading.Event()

    def cancel(self):
        """Request cooperative cancellation; the search loops stop at the next
        iteration boundary and the run emits ``cancelled``."""
        self._cancel.set()

    def run(self):
        from xslope.search import AnalysisCancelled
        try:
            if self._analysis == "reliability":
                self._run_reliability()
            elif self._analysis == "auto_search":
                self._run_search()
            else:
                self._run_single()
        except AnalysisCancelled:
            print("Run cancelled.")
            self.cancelled.emit()
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
                                     num_slices=self._num_slices)
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
            if not (sd.get("circular") and sd.get("circles")):
                self.failed.emit("Circular search needs at least one starting circle.")
                return
            print(f"Searching for the critical circular surface with "
                  f"{self._method.upper()}{self._rapid_tag()}…")
            fs_cache, converged, search_path, circle_cache = circular_search(
                sd, self._method, rapid=self._rapid, num_slices=self._num_slices,
                diagnostic=self._diagnostic, cancel_check=self._cancel.is_set)
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
                diagnostic=self._diagnostic, cancel_check=self._cancel.is_set)
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

    # --- reliability -----------------------------------------------------
    def _run_reliability(self):
        from xslope.advanced import reliability

        sd = self._sd
        circular = self._surface == "circular"
        if circular and not (sd.get("circular") and sd.get("circles")):
            self.failed.emit("Reliability (circular) needs at least one starting circle.")
            return
        if not circular and not sd.get("non_circ"):
            self.failed.emit("Reliability (non-circular) needs a starting "
                             "non-circular surface.")
            return
        print(f"Running reliability — {self._method.upper()}, "
              f"{'circular' if circular else 'non-circular'}{self._rapid_tag()}…")

        def cb(done, total, label):
            self.progress.emit(int(done), int(total) if total is not None else -1, str(label))

        ok, result = reliability(sd, self._method, rapid=self._rapid, circular=circular,
                                 debug_level=1 if self._diagnostic else 0,
                                 progress_callback=cb, cancel_check=self._cancel.is_set)
        if not ok:
            self.failed.emit(str(result))
            return
        # The MLV entry carries the standard solver_result for the Solution view.
        mlv = result["fs_cache"][0]["result"] if result.get("fs_cache") else None
        solver = mlv.get("solver_result") if isinstance(mlv, dict) else None
        self.succeeded.emit({"reliability": result,
                             "slice_df": result.get("critical_slices"),
                             "failure_surface": result.get("critical_surface"),
                             "results": solver, "search": None})

    def _rapid_tag(self):
        return " (rapid drawdown)" if self._rapid else ""
