"""Background run controllers (Phase 3).

Long engine calls run on a ``QThread`` so the UI stays responsive. ``LemRunner``
performs an LEM analysis — single surface or auto-search, circular or
non-circular — and emits the result (or an error) back to the GUI thread. Engine
``print`` output streams to the Log pane live via the thread-safe stdout tee
installed by the main window (``_LogStream``), so the worker does not touch any Qt
widget itself.
"""

from __future__ import annotations

import gc
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
#:
#: This is a MIRROR of ``xslope.mesh._REFINE_THIN_ELEMS``, kept here so the module
#: docstrings and log wording can name the number without importing the mesher at
#: import time. Neither element family is sized from this constant — both go through
#: ``mesh.thin_zone_plan``, which reads the mesher's own — so the two cannot drift
#: into meshing differently; ``test/refine_thin_zones_check.py`` asserts they agree.
THIN_ZONE_ELEMS = 4

#: The toggle's default, ON. An under-resolved thin zone does not fail — it returns
#: a factor of safety that is simply too high — and a section with no thin zone
#: pays nothing for the option, because the detector finds no zone and no size is
#: added. Measured on the Griffiths thin-band section at a 3-unit target, as the
#: seam's width over the median element EDGE length inside it: the band carries 1.1
#: element rows on a default triangular mesh and 1.0 on a quad one, against 4.0 on
#: both with this on.
#:
#: Rows are counted by edge length rather than by an element's root-area, which is
#: the same thing for a quad but overstates a triangle by about 1.5x — the proxy that
#: had the triangular path reading 4.7 rows while it was in fact delivering 3.
REFINE_THIN_ZONES_DEFAULT = True


def _len_unit(slope_data):
    """`` m``/`` ft`` for a model that declares a unit system, ``''`` for one that
    does not — appended to a length in the log rather than assumed."""
    try:
        from xslope.units import labels, normalize_unit_system
        lbl = labels(normalize_unit_system((slope_data or {}).get("unit_system")))
        return f" {lbl['length']}" if lbl.get("length") else ""
    except Exception:
        return ""


def thin_zone_refinement(polygons, target_size, materials=None, skip_declared=False):
    """The thin material zones *Refine thin zones* will act on, each with the local
    element size it gets and a name to report it under.

    One derivation and one mechanism for both element families. The SIZE and the set
    of zones come from ``mesh.thin_zone_plan`` — a zone's own thickness over
    ``mesh._REFINE_THIN_ELEMS`` rows, capped at the global target and at six times
    finer than it — so the two families cannot mesh a zone differently, and so the Log
    pane states the same numbers the mesher used. Zones are measured on the MATERIAL:
    polygons sharing a ``mat_id`` are unioned first, because a layer that the section
    geometry splits into pieces (a bench, a step in the ground surface) is not thin
    just because its pieces are.

    Both families reach the mesher the same way: ``thin_zones`` in ``refine_features``,
    which builds a size field over each zone. The quadrilateral path used to send a
    ``size_regions`` entry instead — a local size region — because a triangular mesh
    pinned a zone's boundary with transfinite curve constraints and a declared size
    alone could not resolve across a band, while a quad mesh had no such constraints.
    Neither family has them now, so one mechanism serves both. Either route delivers
    the four rows today (measured on the Griffiths seam: 4.02 on tri6, 4.03 on quad4,
    both ways); this is the one the mesher derives internally from its own plan, so
    the automatic refinement never has to be disguised as a Size the user declared —
    which is what a size region is, and what preflight and the input warnings read it
    as. It costs about 18 % more nodes there than the size-region route (35318 against
    30033 on tri6), because the band holds the fine size for two element widths beyond
    the zone before it starts growing back.

    ``skip_declared`` drops a zone whose polygon already declares a Size. It defaults
    off and no caller sets it: the size field composes with a declared size by taking
    the minimum rather than replacing it, so an automatic refinement cannot overwrite
    what the user asked for, and honouring a too-coarse declaration would leave the
    zone unresolved. (It existed for the size-region route, where the two WERE the same
    mechanism and the automatic one would have overwritten the declared one.)

    Returns ``[{'name', 'width', 'size', 'polygon': ring}, ...]``, ordered as the
    mesher takes them. ``'polygon'`` is the zone's ring, in the shape
    ``build_mesh_from_polygons(size_regions=...)`` accepts, so a caller that wants a
    local size region out of the same plan can build one.
    """
    from xslope.mesh import thin_zone_plan

    coords, declared, mat_ids = [], [], []
    for i, p in enumerate(polygons):
        if isinstance(p, dict):
            coords.append(list(p.get("coords") or []))
            declared.append(p.get("size"))
            mat_ids.append(p.get("mat_id"))
        else:
            coords.append(list(p))
            declared.append(None)
            mat_ids.append(None)
    group_ids = ([m if m is not None else ("_", i) for i, m in enumerate(mat_ids)]
                 if any(m is not None for m in mat_ids) else None)

    def _name(indices):
        for i in indices:
            try:
                nm = (materials or [])[int(mat_ids[i])].get("name")
            except (IndexError, TypeError, ValueError, AttributeError):
                nm = None
            if nm:
                return str(nm)
        return f"material zone #{indices[0] + 1}" if indices else "thin zone"

    out = []
    for zone in thin_zone_plan(coords, target_size, group_ids=group_ids,
                               declared_sizes=declared if skip_declared else None):
        idxs = zone.get("poly_indices") or [zone.get("poly_index")]
        idxs = [i for i in idxs if isinstance(i, int) and 0 <= i < len(coords)]
        if zone["kind"] == "whole":
            if not idxs:
                continue                            # off the end
            ring = list(zone.get("ring") or coords[idxs[0]])
        else:
            xmin, ymin, xmax, ymax = zone["bbox"]
            ring = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
        out.append({"name": _name(idxs), "width": float(zone["width"]),
                    "size": float(zone["size"]), "polygon": ring})
    return out


def thin_zone_size_regions(polygons, target_size):
    """``build_mesh_from_polygons(size_regions=...)`` entries for the thin zones —
    :func:`thin_zone_refinement` reduced to a ring and a size.

    Not what the Build-mesh dialog passes any more: both element families reach the
    mesher through ``refine_features=['thin_zones']``. This stays because it is the
    derivation stated as a plain (ring, size) list — what a caller building its own
    mesh outside the dialog wants, and what the standing check measures the sizing
    rule through.

    Takes no row count: the sizing rule is ``mesh.thin_zone_plan``'s, and a parameter
    here that no longer reached it would be a knob that silently did nothing."""
    return [{"polygon": z["polygon"], "size": z["size"]}
            for z in thin_zone_refinement(polygons, target_size)]


#: Nesting depth of :func:`suspend_cycle_gc` / :func:`resume_cycle_gc`, and whether
#: the collector was running when the outermost suspension took it.
_GC_LOCK = threading.Lock()
_GC_DEPTH = 0
_GC_WAS_ENABLED = False


def suspend_cycle_gc():
    """Stop Python's CYCLIC collector for the duration of a background run.

    A Qt object's C++ half must be destroyed on the GUI thread. Shiboken enforces
    that: when a wrapper's last Python reference is released anywhere else, it does
    not delete the C++ object there — it queues the deletion onto the main thread
    with ``Py_AddPendingCall`` (``mainThreadDeletionHandler``), to be run at the
    next Python bytecode the GUI thread executes. That queue is what turns a
    background thread into a crash on the GUI thread, and it has been observed to
    fault (PySide 6.11, macOS): the deferred destructor ran against an object that
    was already gone, from inside ``QApplicationPrivate::dispatchEnterLeave`` — the
    first Python a window sitting in ``exec()`` runs after a mouse move.

    Nothing here hands a widget to a worker (``test/thread_safety_check.py`` holds
    every runner to that), but a reference does not have to be handed over to be
    released in the wrong place. The cyclic collector runs on whichever thread
    trips its allocation threshold, and during a run that is the WORKER: a run
    starts by importing the engine — matplotlib, scipy, pandas — and then allocates
    for seconds, while a window idling in ``exec()`` executes no Python at all. It
    was measured at 33 of 114 collections on the worker in one LEM run. Any Qt
    object that is cycle garbage at that moment — a dialog and its combo boxes held
    only by the connections between them — is therefore finalized off the GUI
    thread, through no reference the run ever had.

    So the collector is suspended while a run is in flight and resumed, with one
    explicit collection, on the GUI thread when it ends. What accrues in between is
    cycle garbage only — refcounted objects (the solver's arrays and frames) are
    freed as they always were — and it is reclaimed at the end of the run, in the
    one place where reclaiming it is safe.

    Nested and thread-safe: mesh builds, seepage staging and a run can overlap, and
    the collector comes back only when the last of them is done. A caller that had
    already disabled it keeps it disabled.
    """
    global _GC_DEPTH, _GC_WAS_ENABLED
    with _GC_LOCK:
        if _GC_DEPTH == 0:
            _GC_WAS_ENABLED = gc.isenabled()
            gc.disable()
        _GC_DEPTH += 1


def resume_cycle_gc():
    """Undo one :func:`suspend_cycle_gc`. The collector comes back at depth zero,
    and the collection it was denied happens right there — but only when this runs
    on the main thread, because a collection is exactly what must not happen
    anywhere else."""
    global _GC_DEPTH
    with _GC_LOCK:
        if _GC_DEPTH == 0:
            return
        _GC_DEPTH -= 1
        if _GC_DEPTH:
            return
        if _GC_WAS_ENABLED:
            gc.enable()
    if threading.current_thread() is threading.main_thread():
        gc.collect()


def cycle_gc_suspended():
    """Whether a background run currently holds the collector. Read by the checks."""
    return _GC_DEPTH > 0


class RunnerThread(QThread):
    """The base every background run shares: a thread that suspends the cyclic
    collector for as long as it is alive.

    :meth:`start` takes the suspension on the GUI thread — before the worker exists,
    so there is no window in which it is running unguarded — and ``finished``
    returns it. ``finished`` is emitted by the worker but this runner object lives
    on the GUI thread, so the connection is a queued one and the resumption (and its
    collection) happen there. It is connected in ``__init__``, ahead of any slot the
    main window attaches, so the collector is back before the window disposes of the
    runner.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._gc_held = False
        self.finished.connect(self._release_cycle_gc)

    def start(self, *args, **kwargs):
        suspend_cycle_gc()
        self._gc_held = True
        try:
            super().start(*args, **kwargs)
        except BaseException:
            self._release_cycle_gc()
            raise

    def _release_cycle_gc(self):
        if self._gc_held:
            self._gc_held = False
            resume_cycle_gc()


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
            size_regions = extract_size_regions(sd)
            factor = float(options.get("refine_factor", 3.0))
            features = set(_REFINE_FEATURES) if options.get(
                "refine_near_features", False) else set()
            refine_msg = (f", refine x{factor:g} near features" if features else "")
            # Refine thin zones — a guarantee that a band too thin to mesh does not
            # quietly come back with a factor of safety that is too high. It only
            # ever ADDS: unchecked, the feature refinement above is left exactly as
            # the user set it. One mechanism for both element families.
            #
            # The plan is derived HERE and reported, because
            # this option changes the mesh without being asked for on the run and a
            # user who has not been told which zones it touched, and how finely,
            # cannot tell an intended refinement from an over-refinement.
            thin_msg, thin_plan = "", []
            if options.get("refine_thin_zones", REFINE_THIN_ZONES_DEFAULT):
                thin_plan = thin_zone_refinement(polygons, target,
                                                 materials=sd.get("materials"))
                if thin_plan:
                    # Only ask for the refine-features path when there is a zone for
                    # it to act on. A section with none would otherwise be handed a
                    # refine_factor that finds no feature and changes nothing, and a
                    # build that quietly carries an inert refinement request is one
                    # more thing to rule out when a mesh comes out unexpected.
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
            # Name every zone the toggle refined, the size it got, and the resolution
            # that size buys. Nothing is printed when nothing was refined, so a
            # section with no thin zone reads exactly as it did before the option
            # existed — the lines appear only where the mesh actually changed.
            for _z in thin_plan:
                print(f"  thin zone refined: {_z['name']!r} -> local size "
                      f"{_z['size']:.3g} (~{_z['width'] / _z['size']:.0f} rows across "
                      f"{_z['width']:.3g}{_len_unit(sd)})")
            if thin_plan:
                print("  ('Refine thin zones' in the Build mesh dialog disables this)")
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


class SeepRunner(RunnerThread):
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


class FemRunner(RunnerThread):
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
                # What the run chose and what its trials found. solve_ssrm returns
                # these on the RESULT, and the bundle carried only
                # result["last_solution"] — the field — so the criterion that
                # decided every trial, the interval the search ended on, the
                # trials themselves and the zones the reduction was confined to
                # were dropped at this line and could never be saved or reported.
                from xslope.fem import ssrm_run_record
                self.succeeded.emit({"fem_data": fem_data,
                                     "solution": result["last_solution"],
                                     # The at-failure (unconverged) mechanism field the
                                     # deformation/vector panels render (None if absent).
                                     "failure_solution": result.get("failure_solution"),
                                     "FS": fs, "analysis": "ssrm",
                                     "meta": ssrm_run_record(result, fem_data,
                                                             opts)})
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


class LemRunner(RunnerThread):
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
        self._opts = dict(options or {})
        self._cancel = threading.Event()

    def cancel(self):
        """Request cooperative cancellation; the search loops stop at the next
        iteration boundary and the run emits ``cancelled``."""
        self._cancel.set()

    def run(self):
        from xslope.search import AnalysisCancelled
        try:
            self._run()
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

    def _run(self):
        """One analysis, run the way every path that runs one runs it.

        The work itself is :func:`xslope.search.run_lem_analysis` — the search or
        the specified surface, the model's own search window, and the bundle the
        result views and the report both read. What is left here is the thread's
        own business: the cancel flag, and turning a run that could not be made
        into the message the Run box shows.
        """
        from xslope.search import AnalysisError, run_lem_analysis

        o = self._opts
        try:
            bundle = run_lem_analysis(
                self._sd, o["method"],
                analysis=o.get("analysis", "single_surface"),
                surface=o.get("surface", "circular"),
                num_slices=o.get("num_slices", 40),
                rapid=o.get("rapid", False),
                composite=o.get("composite", False),
                grid_seed=o.get("grid_seed", False),
                diagnostic=o.get("diagnostic", False),
                cancel_check=self._cancel.is_set,
                fs_tol=o.get("fs_tol"), tol=o.get("tol"),
                max_iter=o.get("max_iter"),
                min_slip_depth=o.get("min_slip_depth"))
        except AnalysisError as e:
            self.failed.emit(str(e))
            return
        if bundle.get("results") is None:
            # The surface was built and the method did not converge on it. The
            # run has no result to show, so it is reported as a failed run with
            # the solver's own reason.
            self.failed.emit(f"No solution: {bundle.get('failure')}")
            return
        self.succeeded.emit(bundle)


class SensitivityRunner(RunnerThread):
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


class ReliabilityRunner(RunnerThread):
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


class ReportCancelled(BaseException):
    """Raised inside the report builder to stop it at a figure boundary.

    Deliberately NOT an ``Exception``. ``xslope.report`` wraps every progress
    call in ``except Exception: pass`` — a broken progress line must never cost
    the user a report — and :func:`~xslope.report.generate_report` wraps the
    whole build the same way to turn a failure into ``(False, message)``. Both
    guards let a ``BaseException`` through, which is the only seam a cancel can
    use without the builder growing a cancel argument it has no other use for.
    The temporary figure directory is removed on the way out by the builder's own
    ``finally``, so a cancelled run leaves nothing behind.

    ``test/report_check.py`` asserts that seam directly, so a later edit that
    widened either guard to ``BaseException`` would be caught rather than
    silently turning Cancel into a button that does nothing.
    """


#: What the report's progress bar counts, beyond the figures: the document
#: itself. ``xslope.report`` announces figures and nothing else, so the sections
#: that carry no figure, and the ``.docx`` write that follows them, are one step.
REPORT_WRITE_STEPS = 1

#: The label that step carries.
REPORT_WRITE_LABEL = "writing the Word document"

#: And the finish, which is reported as an indeterminate stretch: the program
#: that lays the pages out is driven over a single call that returns when it is
#: finished and says nothing while it runs, so there is no honest number to
#: show.
REPORT_FINALIZE_LABEL = "building the page numbers…"


class ReportRunner(RunnerThread):
    """Builds an Analysis Report off the GUI thread.

    A report of five methods renders eleven figures and finishes in Word, which
    is tens of seconds with the window frozen if it runs where the window does.
    This is the same treatment every other long call gets: a QThread, progress
    back to the status bar, and nothing Qt touched inside :meth:`run`.

    The run has two phases, and :attr:`progress` reports them differently
    because they are differently knowable:

      * **The solves and figures, then the document** — determinate, over
        ``planned_steps(...) + REPORT_WRITE_STEPS`` steps. A report asked for a
        method the run never solved solves it first, and that solve is a step
        like a figure, announced by name: it is the longest thing in such a
        build, and a bar standing still through a search says the report has
        hung. Each step is announced by ``xslope.report``'s own callback, with
        the label it produces ("the force diagram — Spencer's Method"). The
        builder revises its own total upward once the methods it had to run
        exist, since what their blocks draw is not knowable before that, and the
        bar follows the revised count rather than the estimate. The last
        announcement is followed by the sections that carry no figure and by the
        ``.docx`` write, none of which the builder announces, so the runner names
        that stretch itself (:data:`REPORT_WRITE_LABEL`) as the last step is
        announced. The bar does not advance for it until the document is
        written — it is the one step that is named before it is finished, rather
        than a bar that stands still through the longest step in the run.
      * **Word** — indeterminate (``total`` of -1, the convention the main
        window's progress slot already reads), and only when the caller asked
        for the finish.

    Cancel is cooperative and bounded by the work: the flag is read at the next
    figure the builder announces, and — because the longest stretch of a report
    that has methods to run is inside one of those runs — at every iteration
    boundary of the solves themselves, through the ``cancel_check`` the builder
    is given. :class:`ReportCancelled` and
    :class:`~xslope.search.AnalysisCancelled` unwind the build from wherever it
    was. Once the last step is announced nothing checks it again — the document
    write is one call into python-docx and the finish is one call into the
    program that lays the pages out, and neither may be interrupted halfway with
    the user's report as the thing at risk.

    Emits ``succeeded`` with :func:`~xslope.report.generate_report`'s own result
    dict plus ``finalized`` (whether the page numbers were built) and ``fmt``.
    """

    succeeded = Signal(object)
    failed = Signal(str)
    cancelled = Signal()
    progress = Signal(int, int, str)   # done, total (-1 = indeterminate), label

    def __init__(self, slope_data, solutions, options, path, fmt=None,
                 finalize=False, parent=None):
        super().__init__(parent)
        self._sd = slope_data
        self._solutions = solutions or {}
        self._opts = dict(options or {})
        self._path = path
        self._fmt = fmt
        # Whether the Word finish is switched on is a QSettings question, so it
        # is answered where QSettings lives (the GUI thread) and passed in.
        self._finalize = bool(finalize)
        self._steps = 0                 # set from the builder's own count
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    def run(self):
        from xslope.search import AnalysisCancelled

        try:
            self._generate()
        except (ReportCancelled, AnalysisCancelled):
            # Two seams, one meaning: the flag was read at a figure boundary, or
            # inside a solve the report was making. Either way the user stopped
            # it (:class:`ReportCancelled`, :class:`~xslope.search.AnalysisCancelled`).
            print("Report cancelled.")
            self.cancelled.emit()
        except Exception:
            traceback.print_exc()   # streams to the Log pane via the stdout tee
            self.failed.emit("The report could not be generated — see the Log "
                             "pane for details.")

    def _generate(self):
        from xslope.report import (generate_report, methods_to_run,
                                   planned_figures, planned_steps,
                                   resolve_options)

        opts = dict(self._opts)
        resolved = resolve_options(opts)
        figures = planned_figures(self._sd, self._solutions, resolved)
        solves = len(methods_to_run(self._sd, self._solutions, resolved))
        # Revised by the builder as it goes: a method it has to run first is a
        # step of its own, and the figures that method's block draws are not
        # counted until the run exists.
        self._steps = planned_steps(self._sd, self._solutions, resolved)
        opts["progress"] = self._figure_cb()
        # Stopping a report that is searching must stop the search, not wait for
        # it: the builder polls this inside every solve it makes.
        opts["cancel_check"] = self._cancel.is_set
        # A report with every figure switched off still has a document to write,
        # and the builder will announce nothing at all: name the step up front so
        # the bar is determinate from the start either way. What the first steps
        # ARE decides what that name is — a build that has methods to run is
        # running them, and calling that "rendering the figures" names the wrong
        # minutes.
        self.progress.emit(0, self._steps + REPORT_WRITE_STEPS,
                           "running the methods the report was asked for" if solves
                           else "rendering the figures" if self._steps
                           else REPORT_WRITE_LABEL)

        # The figure count is not knowable yet where methods have still to be
        # run — the blocks they produce carry figures of their own — so the
        # opening line says what is about to happen and the count is stated when
        # it is a count of what was drawn.
        if solves:
            print(f"Generating the report — running {solves} "
                  f"method{'' if solves == 1 else 's'} it was asked for, then "
                  f"its figures…")
        else:
            print(f"Generating the report — {figures} "
                  f"figure{'' if figures == 1 else 's'}…")
        ok, out = generate_report(self._sd, self._solutions, opts, self._path,
                                  fmt=self._fmt)
        if not ok:
            self.failed.emit(str(out))
            return
        if solves:
            drawn = len(out.get("figures") or ())
            print(f"Report built — {drawn} figure{'' if drawn == 1 else 's'}.")
        total = self._steps + REPORT_WRITE_STEPS
        self.progress.emit(total, total, REPORT_WRITE_LABEL)

        finalized = False
        if self._finalize:
            from .report_dialog import document_finish
            self.progress.emit(0, -1, REPORT_FINALIZE_LABEL)
            finalized, _msg = document_finish(self._path, True)
        self.succeeded.emit(dict(out, finalized=finalized, fmt=self._fmt))

    def _figure_cb(self):
        """``xslope.report``'s progress callback: one signal per step the builder
        announces, and the cancel flag read at each one.

        The builder's own count comes with every step, and it is the one the bar
        is drawn against: a build that had to run a method before it could draw
        that method's figures knows how many there are only once it has.
        """
        def cb(done, planned, label):
            if self._cancel.is_set():
                raise ReportCancelled()
            if int(planned) > 0:
                self._steps = int(planned)
            total = self._steps + REPORT_WRITE_STEPS
            self.progress.emit(int(done), total, str(label))
            if int(done) >= self._steps:
                self.progress.emit(int(done), total, REPORT_WRITE_LABEL)
        return cb


class UpdateCheckRunner(RunnerThread):
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


class UpdateDownloadRunner(RunnerThread):
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


class AssistantModelsRunner(RunnerThread):
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
