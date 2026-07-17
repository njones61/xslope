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
                                 extract_constraint_line_geometry, MeshInputError)
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
            print(f"Building {element_type} mesh, target size {target:.3g}{extra}…")
            mesh = build_mesh_from_polygons(polygons, target_size=target,
                                            element_type=element_type,
                                            lines=constraint_lines or None)
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
    a later set still runs."""

    succeeded = Signal(object)
    failed = Signal(str)

    def __init__(self, slope_data, options, parent=None):
        super().__init__(parent)
        self._sd = slope_data
        self._options = options

    def run(self):
        from xslope.seep import build_seep_data, run_seepage_analysis, SeepInputError
        sd = self._sd
        mesh = sd.get("mesh")
        if mesh is None:
            self.failed.emit("No mesh available — build a mesh first.")
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
            except SeepInputError as e:
                # Expected, user-actionable input problem — show the message, no
                # scary traceback needed (still logged for the record).
                print(f"{label}: {e}")
                errors.append(f"{label}: {e}")
            except Exception as e:
                traceback.print_exc()
                errors.append(f"{label}: {e}  (see the Log pane for details.)")
        if errors:
            self.failed.emit("\n\n".join(errors))


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
            print("Building FEM data…")
            fem_data = build_fem_data(sd, mesh)
            opts = self._options
            analysis = opts.get("analysis", "ssrm")
            if analysis == "reliability":
                from xslope.advanced import reliability_fem

                def rel_cb(done, total, label):
                    self.progress.emit(int(done), int(total) if total else -1, str(label))

                print(f"Running FEM reliability (SSRM, F in "
                      f"[{opts.get('F_min', 1.0):g}, {opts.get('F_max', 2.0):g}])…")
                # Reliability uses its own tight bisection tolerance (the dialog's
                # "Reliability tol" field, not the single-run "Tolerance"): TSPM
                # amplifies FS imprecision, so a coarse band would make
                # beta/reliability jitter between runs.
                success, result = reliability_fem(
                    sd, mesh=mesh, F_min=opts.get("F_min", 1.0),
                    F_max=opts.get("F_max", 2.0),
                    tolerance=opts.get("reliability_tol", 0.001),
                    failure_criterion=opts.get("failure_criterion", "non_convergence"),
                    debug_level=1, cancel_check=self._cancel.is_set,
                    progress_callback=rel_cb)
                if not success:
                    self.failed.emit(f"Reliability failed: {result}")
                    return
                self.succeeded.emit({
                    "fem_data": fem_data,
                    "solution": result["mlv_solution"]["last_solution"],
                    "FS": result["F_MLV"], "analysis": "reliability",
                    "reliability": result})
            elif analysis == "single":
                F = opts.get("F", 1.0)
                print(f"Solving FEM (single trial, F={F:g})…")

                def fem_cb(frac, label):
                    self.progress.emit(int(frac * 100), 100, str(label))

                solution = solve_fem(fem_data, F=F, debug_level=1,
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
                    cancel_check=self._cancel.is_set, progress_callback=cb)
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
        """Tolerance kwargs to forward to the search functions. ``tol`` is only
        accepted by ``circular_search`` (noncircular has no geometric tol)."""
        kw = {}
        if self._fs_tol is not None:
            kw["fs_tol"] = self._fs_tol
        if self._max_iter is not None:
            kw["max_iter"] = self._max_iter
        if self._min_slip_depth is not None:
            kw["min_slip_depth"] = self._min_slip_depth
        if circular and self._tol is not None:
            kw["tol"] = self._tol
        return kw

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
                                 progress_callback=cb, cancel_check=self._cancel.is_set,
                                 fs_tol=self._fs_tol, tol=self._tol,
                                 max_iter=self._max_iter, composite=self._composite,
                                 seed=self._seed)
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


class SensitivityRunner(QThread):
    """Runs a sensitivity / design study off the GUI thread.

    A sweep is N sequential LEM solves (no multiprocessing). Progress is reported
    per solve and the run cancels cooperatively (each engine sweep is passed a
    ``cancel_check``; an in-flight solve finishes, then the sweep stops at the next
    point). A thin caller of ``xslope.sensitivity``:

      * design mode  -> ``design()``  -> bundle {'kind':'design', df, crossing, …}
      * sensitivity  -> ``sensitivity()`` per parameter (full FS-vs-value curves,
        for click-through) + ``tornado_from_sweeps()`` -> bundle
        {'kind':'sensitivity', sweeps, tornado, method}.
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
            if self._opts.get("mode") == "design":
                self._run_design()
            else:
                self._run_sensitivity()
        except AnalysisCancelled:
            print("Sweep cancelled.")
            self.cancelled.emit()
        except Exception:
            traceback.print_exc()   # streams to the Log pane via the stdout tee
            self.failed.emit("Sweep failed — see the Log pane for details.")

    def _run_design(self):
        from xslope.sensitivity import design
        o = self._opts
        emode = o.get("engine_mode", "lem")
        total = int(o["steps"]) + 1        # points + base

        def cb(done, _t, label):
            self.progress.emit(int(done), total, str(label))

        out = "q" if emode == "seep" else "FS"
        engine_tag = {"lem": o.get("method", "spencer"), "fem": "SSRM",
                      "seep": f"BC {o.get('seep_opts', {}).get('bc', 1)}"}.get(emode, emode)
        print(f"Design sweep ({emode}): {o['param']} from {o['low']:g} to "
              f"{o['high']:g} in {o['steps']} steps, target {out} = "
              f"{o['target_fs']:g} ({engine_tag})…")
        ok, res = design(self._sd, param=o["param"], low=o["low"], high=o["high"],
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
        self.succeeded.emit({"kind": "design", **res})

    def _run_sensitivity(self):
        import numpy as np
        from xslope.sensitivity import sensitivity, tornado_from_sweeps
        o = self._opts
        emode = o.get("engine_mode", "lem")
        specs = o["params"]
        if not specs:
            self.failed.emit("Add at least one parameter to sweep.")
            return
        n = int(o["n"])
        method = o.get("method", "spencer")
        # display method label for the tornado axis/title (an LEM method only
        # makes sense in LEM mode; FEM/Seep carry their own quantity label)
        disp_method = {"lem": method, "fem": "SSRM", "seep": ""}.get(emode, method)

        def n_points(spec):
            if spec.get("values") is not None:
                return len(spec["values"])
            return n

        total = sum(n_points(s) + 1 for s in specs)   # +1 base per parameter
        count = [0]

        def cb(_done, _t, label):
            count[0] += 1
            self.progress.emit(count[0], total, str(label))

        sweeps, base_fs, out = {}, None, "FS"
        for i, spec in enumerate(specs):
            ref = spec["ref"]
            if spec.get("low") is not None and spec.get("high") is not None:
                values = list(np.linspace(spec["low"], spec["high"], n))
                rel = 0.5
                tag = f"±σ [{spec['low']:g}, {spec['high']:g}]"
            else:
                values = None
                rel = spec.get("rel_range", 0.2)
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
        self.succeeded.emit({"kind": "sensitivity", "sweeps": sweeps,
                             "tornado": tornado, "method": disp_method})
