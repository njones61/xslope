"""PythonKernel — a persistent in-process namespace the assistant drives.

``run(code)`` execs a snippet in a namespace that has the engine (``xslope``),
``numpy``/``pandas``/``matplotlib.pyplot``, and the live project handed in
(``doc``, ``slope_data``, ``results``). Variables persist across calls like a
notebook. stdout/stderr are captured (not spilled to the Log pane) and any new
pyplot figures are saved to temp PNGs so the chat can show them. Runs on the GUI
thread (the assistant marshals tool execution there), so mutating ``slope_data``
and triggering a re-render is safe.
"""

from __future__ import annotations

import os
import tempfile
import traceback
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO

# Figure formats a snippet might savefig itself — if it did, we don't also
# auto-save its open figures (which would show the same plot twice).
_FIG_EXTS = (".png", ".pdf", ".svg", ".jpg", ".jpeg")

#: Parsed corpus index, built once per process (the JSON is ~650 KB).
_CORPUS_CACHE = None

#: The inputs a geometry resync reads. Compared before and after a snippet to
#: tell an edit made on the model's OWN geometry source from one made on the
#: derived copy, which the resync silently overwrites.
_GEOM_KEYS = ("profile_lines", "polygons", "max_depth")

#: What the resync says when a geometry edit was made on the wrong source. The
#: model reads these on its own tool result, so they name the fix, not the fault.
POLYGON_EDIT_WARNING = (
    "WARNING: polygons were edited on a profile-line model and have been rebuilt "
    "from profile_lines; edit profile_lines instead (and the ground surface if it "
    "is separate), then call resync_geometry(). The polygon edit did not take.")
PROFILE_EDIT_WARNING = (
    "WARNING: profile_lines were added on a polygon-native model; the polygons are "
    "now rebuilt from those profile lines, so the model's own polygons no longer "
    "apply. Edit polygons instead on this model.")


def _geometry_native(slope_data):
    """Which geometry source this model is built on — ``"profile"`` when it has
    profile lines (:func:`studio.editors._resync_geometry` rebuilds the polygons
    from them), ``"polygon"`` when the polygons are the source, None for a model
    that has neither yet."""
    if (slope_data or {}).get("profile_lines"):
        return "profile"
    if (slope_data or {}).get("polygons"):
        return "polygon"
    return None


def _geometry_sigs(slope_data):
    """Per-key JSON signatures of the geometry inputs (shapely objects fall back
    to their WKT via ``str``), so a snippet's geometry edits can be located."""
    import json
    out = {}
    for key in _GEOM_KEYS:
        try:
            out[key] = json.dumps((slope_data or {}).get(key), sort_keys=True,
                                  default=str)
        except Exception:
            out[key] = None
    return out


def _declared_lem_method(slope_data, method=None):
    """The method a run made without one uses: the method the MODEL declares
    (main!D14 — what Studio's Run LEM dialog opens on), else spencer.

    ``'all'`` is the batch-report sweep and names no single method, so it seeds
    nothing and the fallback stands.
    """
    if method:
        return str(method).lower()
    declared = str((slope_data or {}).get("lem_method") or "").lower()
    return declared if declared and declared != "all" else "spencer"


def _surface_keys(bundle, slope_data):
    """The failure surface this run settled on, as plain numbers on the result
    dict: ``Xo``, ``Yo``, ``R``, ``Depth`` (the tangent elevation) and the trace's
    ``x_entry`` / ``x_exit``.

    Every key is always present — None on a non-circular surface — so reading the
    critical circle off a result never depends on which branch produced it. Entry
    and exit follow the search's own convention (:func:`xslope.search
    ._endpoints_in_ranges`): the ENTRY point is the crest-side (higher ground) end
    of the trace and the EXIT the toe-side one, whichever way the slope faces.
    """
    keys = {"Xo": None, "Yo": None, "R": None, "Depth": None,
            "x_entry": None, "x_exit": None}
    search = bundle.get("search") or {}
    circle = None
    if search.get("kind") == "circular":
        cache = search.get("fs_cache") or []
        if cache:
            circle = cache[0]
    elif search:
        circle = None                       # a non-circular search: no circle
    elif slope_data.get("circular", True) and slope_data.get("circles"):
        circle = slope_data["circles"][0]
    if circle is not None:
        xo, yo = circle.get("Xo"), circle.get("Yo")
        depth, r = circle.get("Depth"), circle.get("R")
        if r is None and yo is not None and depth is not None:
            r = float(yo) - float(depth)
        if depth is None and yo is not None and r is not None:
            depth = float(yo) - float(r)
        keys.update({"Xo": None if xo is None else float(xo),
                     "Yo": None if yo is None else float(yo),
                     "R": None if r is None else float(r),
                     "Depth": None if depth is None else float(depth)})
    surf = bundle.get("failure_surface")
    try:
        coords = list(surf.coords)
    except Exception:
        coords = []
    if len(coords) >= 2:
        (xl, yl), (xr, yr) = coords[0], coords[-1]
        entry, exit_ = ((xl, xr) if yl >= yr else (xr, xl))
        keys["x_entry"], keys["x_exit"] = float(entry), float(exit_)
    return keys


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _corpus_rows():
    """``(topics, rows)`` for the verification corpus.

    ``topics`` maps a topic key to its display label; ``rows`` is one dict per
    cited page — ``{topic, label, title, url}``.

    Two sources, because two installs. A repo checkout has the generated index
    (``docs/verification/corpus_index.json``, the authoritative one). A pip install
    has no docs/ tree, but it does ship the Claude Code skill, whose corpus-index
    block is the same table already rendered — so the same rows are recovered from
    it rather than shipping a second 650 KB copy of the JSON that would then need a
    sync check to stay honest.
    """
    global _CORPUS_CACHE
    if _CORPUS_CACHE is None:
        _CORPUS_CACHE = _corpus_from_json() or _corpus_from_skill() or ({}, [])
    return _CORPUS_CACHE


def _corpus_from_json():
    import json
    path = os.path.join(_repo_root(), "docs", "verification", "corpus_index.json")
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except Exception:
        return None
    labels = data.get("topic_labels") or {}
    counts = (data.get("stats") or {}).get("topics") or {}
    topics = {k: f"{labels.get(k, k)} ({counts.get(k, 0)} models)" for k in labels}
    rows = []
    for entry in data.get("models") or []:
        for topic in entry.get("topics") or []:
            for ref in entry.get("references") or []:
                if ref.get("url"):
                    rows.append({"topic": topic, "label": labels.get(topic, topic),
                                 "title": ref.get("title", ""), "url": ref["url"]})
    return (topics, rows) if rows else None


def _corpus_from_skill():
    """The same rows out of the rendered table in the shipped skill (pip install)."""
    import re
    text = ""
    try:
        from importlib import resources
        text = (resources.files("xslope") / "resources" / "xslope_skill.md").read_text(
            encoding="utf-8")
    except Exception:
        return None
    block = re.search(r"<!-- corpus-index:begin -->(.*?)<!-- corpus-index:end -->",
                      text, re.S)
    if not block:
        return None
    topics, rows = {}, []
    for line in block.group(1).splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 2 or not cells[0] or set(cells[0]) <= set(":-"):
            continue
        label = cells[0]
        if label.lower() == "topic":
            continue
        key = re.sub(r"[^a-z0-9]+", "-", label.lower()).strip("-")
        links = re.findall(r"\[([^\]]+)\]\((https?://[^)]+)\)", cells[1])
        if not links:
            continue
        topics[key] = label
        rows += [{"topic": key, "label": label, "title": t, "url": u}
                 for t, u in links]
    return (topics, rows) if rows else None


class PythonKernel:
    def __init__(self, doc, window=None):
        self._doc = doc
        # The main window, where there is one. Only ever asked for what Studio
        # itself knows and the document does not — which solutions a report may
        # document, and the app's settings — and always through ``getattr``, so
        # the kernel runs unchanged against a document with no window (the
        # guardrail checks).
        self._window = window
        self._ns = {}
        self._seeded = False
        self._outdir = None
        self._fig_seq = 0
        #: Solutions THIS snippet computed, keyed as the document stores them,
        #: with the views that put each in front of the user. Studio's
        #: stale-result sweep runs after the snippet, so a snippet that edited
        #: the model AND then ran it would otherwise lose the run it just made.
        self._fresh_results = {}
        #: Geometry warnings raised during the snippet (see
        #: :meth:`geometry_warnings`), and the state it started from.
        self._geom_warnings = []
        self._geom_watch = None

    @property
    def outdir(self):
        """Folder the assistant's generated files (plots, CSVs, …) are written to,
        so the user can open them. Snippets run with this as the working directory,
        so a naive ``plt.savefig('x.png')`` lands here. Stable for the session."""
        if self._outdir is None:
            self._outdir = os.path.join(tempfile.gettempdir(), "xslope_studio_assistant")
            os.makedirs(self._outdir, exist_ok=True)
        return self._outdir

    def reset(self):
        """Drop all variables; the engine + helpers re-seed on the next run."""
        self._ns = {}
        self._seeded = False
        self._fresh_results = {}
        self._geom_warnings = []

    # --- what a snippet produced, for the caller that has to survive it ------
    def _show_views(self, *calls):
        """Put a solved run in front of the user the way a Run does — the
        results tab for that engine — where there is a window to put it in.

        Best effort by construction: the answer is already attached to the
        model and the session, and the chat carries its plot, so a view that
        cannot be built (no window, a panel that objects) must not turn a
        successful analysis into a failed snippet.
        """
        for name, args in calls:
            try:
                method = getattr(self._window, name, None)
                if callable(method):
                    method(*args)
            except Exception:
                pass

    def _store_result(self, key, value, show=(), lead=None):
        """Attach a solved bundle where Studio attaches it, show it, and record
        that THIS snippet made it.

        ``show`` is the (method name, args) list that raises the engine's result
        tabs and ``lead`` the attribute naming the canvas the run leads with, both
        replayed by :meth:`restore_fresh_results`.
        """
        self._doc.results[key] = value
        self._fresh_results[key] = (value, tuple(show), lead)
        self._show_views(*show)
        self._lead_tab(lead)

    def _lead_tab(self, attr):
        """Bring the tab a run leads with to the front, as the dialog run does."""
        if not attr:
            return
        try:
            canvas = getattr(self._window, attr, None)
            if canvas is not None:
                self._window.view_tabs.setCurrentWidget(canvas)
        except Exception:
            pass

    def restore_fresh_results(self):
        """Re-attach (and re-show) the solutions this snippet computed.

        Studio drops every cached result when the inputs change, which is right
        for a solution that predates the edit and wrong for one made after it: a
        snippet that lays the face back and reruns computes its answer on the model
        as it now stands. That sweep runs after the snippet returns, so the run is
        put back here rather than never stored.
        """
        for key, (value, show, lead) in list(self._fresh_results.items()):
            self._doc.results[key] = value
            self._show_views(*show)
            self._lead_tab(lead)

    def geometry_warnings(self):
        """Warnings about geometry the last snippet edited on the wrong source —
        polygons on a profile-line model, or the reverse. Empty for every snippet
        whose geometry edit reached the model's own source."""
        return list(self._geom_warnings)

    def _note_geometry_source(self, slope_data=None):
        """Check, at the moment before a resync would discard it, whether the
        snippet edited the geometry the model is NOT built from.

        The resync rebuilds polygons from profile_lines on a profile-line model,
        so a snippet that rebuilt ``polygons`` there has its edit reverted with
        nothing said — the failure this exists to name. Called before every resync
        and once more when the snippet ends, since a snippet that never ran
        anything has its polygons overwritten later, by the checks' own resync.
        """
        watch = self._geom_watch
        if not watch:
            return
        native, before = watch
        sd = self._doc.slope_data if slope_data is None else slope_data
        if native is None or sd is not self._doc.slope_data:
            return                      # a copy (a sweep's own model): not the doc
        now = _geometry_sigs(sd)
        if now["max_depth"] != before["max_depth"]:
            return                      # the polygons legitimately rebuild
        if (native == "profile" and now["polygons"] != before["polygons"]
                and now["profile_lines"] == before["profile_lines"]):
            self._warn_geometry(POLYGON_EDIT_WARNING)
        elif (native == "polygon"
                and now["profile_lines"] != before["profile_lines"]):
            self._warn_geometry(PROFILE_EDIT_WARNING)

    def _warn_geometry(self, text):
        if text not in self._geom_warnings:
            self._geom_warnings.append(text)

    def _seed(self):
        import matplotlib
        matplotlib.use("Agg")          # pyplot figures render off-screen; we save them
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import xslope
        self._ns.update({"np": np, "pd": pd, "plt": plt, "xslope": xslope,
                         "OUTPUT_DIR": self.outdir})
        self._ns.update(self._helpers())
        self._seeded = True

    def _helpers(self):
        """Convenience functions seeded into the namespace so the model doesn't
        have to reconstruct the engine pipeline (a common failure mode). Seeded:

        - ``run_lem(method=None, search=False, ...)`` — one LEM solve; returns
          the result dict (with the surface it was solved on: Xo/Yo/R/Depth and
          the trace's x_entry/x_exit) and shows the solution plot. ``method``
          defaults to the one the MODEL declares, else spencer. ``search=True``
          searches for the CRITICAL surface for that method (Studio's Run LEM
          default, and what "the factor of safety of this model" means);
          ``search=False`` solves the surface the model already defines. The
          bundle is attached to the session as ``doc.results['lem_solution']``,
          which is what the report and the result tabs read.
        - ``corpus_index(query=None)`` — verification-corpus rows matching a topic
          or phrase, so a citation is looked up rather than remembered.
        - ``run_seep(bc=1, ...)`` — one steady seepage solve; builds the mesh from
          the file's declared settings if there is none, attaches the solved field
          to the model and the bundle to the session, shows the solution plot.
        - ``run_fem(analysis='ssrm', ...)`` — one finite element run (SSRM factor
          of safety, or a single trial); same mesh handling, attaches the bundle,
          shows the results plot.
        - ``suggest_elastic(material_or_soil_type=None, ...)`` — soil-type E and
          nu for a material that carries none, with the classification it came
          from. A LAST-RESORT fill, never a preference over a stated value.
        - ``generate_report(path=None, **options)`` — the Analysis Report the
          Report dialog builds, written and finished; returns the path.
        - ``resync_geometry(slope_data=None)`` — rebuild derived geometry after an
          in-snippet geometry edit (call inside sweep loops).
        - ``sensitivity(values, apply, param=..., ...)`` — callback-driven FS-vs-
          parameter sweep (you write ``apply(v)``); writes a CSV + one plot.
        - ``list_params(slope_data=None, mode='lem')`` — discover every sweepable
          parameter as a menu of canonical refs (feeds design_sweep / the engine
          sweepers). mode='seep' lists the seepage set instead (each material's
          hydraulic k fields + every specified-head boundary as a seep_bc ref).
        - ``design_sweep(param, low, high, target_fs=1.5, mode='lem', ...)`` — vary
          ONE named parameter and find where the OUTPUT meets a target; renders one
          output-vs-parameter plot with the target highlighted. ``mode`` picks the
          engine: 'lem' (output = FS, default), 'fem' (FS via SSRM — MINUTES per
          step, needs a mesh), 'seep' (total discharge q — needs a mesh; target is a
          q). The engine-driven counterpart to ``sensitivity`` — a parameter ref or
          dict spec, or (exclusive with it) the same ``modify=``/``label=`` callable
          escape hatch ``sensitivity`` takes for geometry. This is the **Design**
          mode of the Parametric study family; ``parametric_design`` is the same
          function under the umbrella name.

        The Parametric study family (all share the ``list_params`` parameter refs):

        - ``parametric_sweep(params=None, plot='scaled', ...)`` — multi-parameter
          sensitivity that renders one plot-family view: ``plot`` in
          {'tornado', 'scaled', 'spider', 'variance', 'rank'} (with ``scaling`` for
          the scaled bars). The one call that reaches every plot type.
        - ``parametric_design(...)`` — alias of ``design_sweep`` (target-FS sweep).
        - ``parametric_back_analysis(param, low, high, target_fs=1.0, ...)`` — forensic
          back-analysis: the value that makes the slope limiting (FS = 1.0), i.e.
          ``design_sweep`` framed for a failure investigation.

        The Reliability family (probabilistic; turns the σ columns into β / P_f):

        - ``reliability_taylor(method='bishop', ...)`` — Taylor Series Probability
          Method (1+2N solves); renders the MLV / F± surface plot.
        - ``reliability_mc(method='bishop', n_samples=..., ...)`` — Monte Carlo
          sampling on a fixed surface; renders the FS histogram.
        - ``reliability_rs(method='bishop', ...)`` — response-surface sampling on a
          fixed surface: a quadratic surrogate fitted to a few dozen real solves,
          gated against held-out ones, then sampled ten million times; renders the
          FS histogram.
        - ``reliability(method='bishop', engine='taylor'|'mc'|'rs')`` — the front
          door (kept under the plain name for back-compat).
        """
        doc = self._doc
        window = self._window

        def resync_geometry(slope_data=None):
            """Rebuild the derived geometry (ground_surface, polygons,
            domain_polygon, tcrack) from the current profile_lines/polygons.

            **Call this after editing geometry inside a snippet** — e.g. a
            parametric sweep that moves a profile point in a loop. The canvas's
            automatic resync runs only ONCE, after the whole snippet returns, not
            between loop iterations — so without this, every iteration analyzes the
            stale original geometry (a classic constant-result bug). `run_lem`
            calls this for you; call it yourself before `generate_slices` /
            `circular_search` / `solve_*` if you drive the pipeline directly."""
            from studio.editors import _resync_geometry
            sd = doc.slope_data if slope_data is None else slope_data
            self._note_geometry_source(sd)   # before the rebuild discards an edit
            _resync_geometry(sd)
            return sd.get("ground_surface")

        def run_lem(method=None, num_slices=40, rapid=False, plot=True,
                    slope_data=None, search=False):
            """Run an LEM analysis on the loaded project and return the result dict
            (includes 'FS'). `method` is one of oms, bishop, janbu, spencer, corps,
            lowe, mprice; left out, it is THE METHOD THE MODEL DECLARES (main!D14 —
            what Studio's Run LEM dialog opens on), or spencer where it declares
            none. Shows the solution plot when plot=True (pass plot=False in sweeps
            to avoid many figures). Rebuilds derived geometry first, so edits to
            profile_lines/polygons this snippet made are reflected.

            `search=True` runs the automated search for the CRITICAL surface for
            this method — what Studio's Run LEM does by default, and what "the
            factor of safety of this model" means. Each method searches for itself,
            so two methods legitimately settle on two surfaces. `search=False`
            (default) solves only the surface the model already defines, which is
            the right call inside a sweep or when the user named that surface.

            The run is :func:`xslope.search.run_lem_analysis` — the entry point
            Studio's Run LEM dialog drives — and its bundle is attached to the
            session (`doc.results['lem_solution']`) exactly as a dialog run attaches
            it, so the result tabs show it and the Analysis Report documents it.

            The returned dict is the solver's, plus the surface it was solved on:
            'Xo', 'Yo', 'R', 'Depth' (the circle, None on a non-circular surface)
            and 'x_entry' / 'x_exit', the crest-side and toe-side ends of the trace.
            """
            import contextlib
            import io as _io
            from xslope.search import AnalysisError, run_lem_analysis

            sd = doc.slope_data if slope_data is None else slope_data
            method = _declared_lem_method(sd, method)
            resync_geometry(sd)        # reflect any in-snippet geometry edits
            if not (sd.get("circles") or sd.get("non_circ")):
                raise RuntimeError("No failure surface defined (add a circle or "
                                   "non-circular surface first).")
            family = ("circular" if (sd.get("circular", True) and sd.get("circles"))
                      else "noncircular")
            # A converging search writes a dozen progress lines the model would
            # have to read back as tokens, and only the last of them is the answer,
            # so a SEARCH runs muted. What the search reports about ITSELF —
            # unsolved trials, admissibility notes — is kept below, because those
            # are findings about the model. A single solve prints little and what
            # it prints is the solver's, so it is left alone.
            muted = (contextlib.redirect_stdout(_io.StringIO()) if search
                     else contextlib.nullcontext())
            try:
                with muted:
                    bundle = run_lem_analysis(
                        sd, method,
                        analysis="auto_search" if search else "single_surface",
                        surface=family, num_slices=num_slices, rapid=rapid,
                        announce=False)
            except AnalysisError as exc:
                raise RuntimeError(str(exc)) from None
            result = bundle.get("results")
            if not isinstance(result, dict):
                raise RuntimeError(
                    ("The search found no solvable surface: %s" if search
                     else "No solution: %s") % (bundle.get("failure")
                                                or "no solution"))
            surface = _surface_keys(bundle, sd)
            where = ""
            if search and surface["Xo"] is not None:
                where = (f" on the circle Xo={surface['Xo']:.2f}, "
                         f"Yo={surface['Yo']:.2f}, R={surface['R']:.2f}")
            head = f"{method} (auto search, {family})" if search else method
            print(f"{head}: FS = {result['FS']:.3f}{where}")
            info = bundle.get("search") or {}
            unsolved = (info.get("unsolved") or {}) if info else {}
            if unsolved.get("sentence"):
                print("  " + unsolved["sentence"])
            for note in (result.get("warnings") or []):
                print(f"  admissibility: {note}")
            # Only a run on the OPEN model is the session's run: a sweep hands in
            # its own copy, and those answers are the sweep's, not the project's.
            if sd is doc.slope_data:
                # The options this run was made under are the session's last LEM
                # options, exactly as the dialog's are: they are what the report
                # labels the run with, and what the Run LEM dialog next opens on.
                if window is not None and hasattr(window, "_last_lem_opts"):
                    try:
                        window._last_lem_opts = dict(bundle.get("options") or {},
                                                     method=method)
                    except Exception:
                        pass
                _store_result("lem_solution", bundle,
                              show=([("_show_search", (bundle["search"],))]
                                    if bundle.get("search") else [])
                                   + [("_show_solution", (bundle,))],
                              lead=("search_canvas" if bundle.get("search")
                                    else "solution_canvas"))
            if plot:
                from xslope.plot import plot_solution
                import matplotlib.pyplot as plt
                plot_solution(sd, bundle["slice_df"], bundle["failure_surface"],
                              result, fig=plt.figure(figsize=(11, 6)))
            return dict(result, **surface)

        def _ensure_mesh(sd, quiet=False):
            """The model's finite-element mesh, building one from the FILE'S OWN
            declared settings if it carries none.

            A mesh is a computed artifact, not an input, so it is attached to the
            model directly (no undo step, no dirty flag) exactly as Studio's own
            mesh build attaches it. The element type and target size come from the
            model (main!D18 / main!D19); with no declared size the automatic one is
            used — the ground-surface width over 100 — which is the Build mesh
            dialog's own default. Reinforcement and pile lines become constraint
            lines and per-polygon Size overrides are honoured, so the mesh a run
            gets here is the mesh the dialog would have built.
            """
            if sd.get("mesh") is not None:
                return sd["mesh"]
            from xslope.mesh import (build_mesh_from_polygons,
                                     extract_constraint_line_geometry,
                                     extract_size_regions, get_material_polygons)
            resync_geometry(sd)
            lines, _n_reinf, _n_pile = extract_constraint_line_geometry(sd)
            polygons = get_material_polygons(sd, reinf_lines=lines)
            element_type = sd.get("element_type") or "tri6"
            target = sd.get("target_size")
            if not target:
                xs = [x for x, _y in sd["ground_surface"].coords]
                target = (max(xs) - min(xs)) / 100.0
            mesh = build_mesh_from_polygons(
                polygons, target_size=float(target), element_type=element_type,
                lines=lines or None,
                element_size_1d=sd.get("element_size_1d"),
                size_regions=extract_size_regions(sd))
            sd["mesh"] = mesh
            if not quiet:
                print(f"Built a {element_type} mesh (target size {float(target):.4g}): "
                      f"{len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements.")
            return mesh

        # Attach a solved bundle where a Run attaches it, raise its result tabs,
        # and record that this snippet is what made it.
        _store_result = self._store_result

        def run_seep(bc=1, tol=1e-4, max_iter=400, plot=True, slope_data=None):
            """Run a STEADY seepage analysis and return the solution dict — the
            seepage counterpart of `run_lem`.

            Builds the seepage data from the model's mesh (building the mesh from the
            file's declared settings when there is none), solves, and does what
            Studio does with the answer: the nodal pore pressures are attached to the
            model (`slope_data['seep_u']` for BC set 1, `['seep_u2']` for set 2), so a
            later stability run with a material set to `u = 'seep'` reads THIS field,
            and the bundle is stored on the session (`doc.results['seep_solutions']`)
            so the report and the results tabs can find it. Shows the head/flow-net
            plot when plot=True (pass plot=False in a sweep).

            `bc` picks the boundary-condition set (1, or 2 for the drawn-down state
            rapid drawdown reads). `tol` is the relative head-change tolerance and
            `max_iter` the cap on the unconfined iteration. Returns the solution dict
            ('head', 'u', 'velocity', 'gradient', 'phi', 'flowrate', 'converged',
            'closure_error', …). For a TRANSIENT march, drive
            `xslope.seep.run_transient_seepage` directly — this helper is the steady
            solve.
            """
            from xslope.seep import (apply_steady_stability_field, build_seep_data,
                                     run_seepage_analysis)
            sd = doc.slope_data if slope_data is None else slope_data
            mesh = _ensure_mesh(sd)
            seep_data = build_seep_data(mesh, sd, seep_bc=bc)
            solution = run_seepage_analysis(seep_data, tol=tol, max_iter=int(max_iter))
            if solution is None:
                raise RuntimeError("Seepage analysis returned no solution.")
            bundle = {"seep_data": seep_data, "solution": solution,
                      "options": {"bc": bc, "tol": tol, "max_iter": int(max_iter)}}
            doc.results.setdefault("seep_solutions", {})[bc] = bundle
            # The solved field belongs to the MODEL, not only to a results tab —
            # this is what a u='seep' stability run reads.
            apply_steady_stability_field(sd, solution, bc=bc, verbose=False)
            _store_result("seep_solutions", doc.results["seep_solutions"],
                          show=[("_show_seep_data", (seep_data, bc)),
                                ("_show_seep_solution", (bc,))])
            q = solution.get("flowrate")
            print(f"seepage (BC set {bc}): converged={solution.get('converged')}"
                  + (f", total flow q = {q:.4g}" if q is not None else ""))
            if plot:
                from xslope.plot_seep import plot_seep_solution
                import matplotlib.pyplot as plt
                plot_seep_solution(seep_data, solution, fig=plt.figure(figsize=(11, 6)))
            return solution

        def run_fem(analysis="ssrm", F=1.0, F_min=None, F_max=None, tolerance=0.01,
                    failure_criterion="non_convergence", min_slip_depth=None,
                    k0=None, tension_srf=None, plot=True, slope_data=None, **kwargs):
            """Run a finite element analysis and return the result bundle — the FEM
            counterpart of `run_lem`. MINUTES, not seconds: say so before starting one.

            `analysis='ssrm'` (default) runs the shear-strength reduction search and
            returns a bundle whose 'FS' is the factor of safety; `analysis='single'`
            runs ONE trial at strength factor `F` and returns the stress/displacement
            field with 'FS' None. Builds the mesh from the file's declared settings if
            the model carries none, attaches the bundle to the session
            (`doc.results['fem_solution']`) the way Studio does, and shows the standard
            results plot when plot=True.

            `F_min`/`F_max` bracket the SSRM search (default: the model's own
            main!D21/D22 values, else 1.0-2.0); `tolerance` closes it. `k0` and
            `tension_srf` default to the model's declared values. `failure_criterion`
            and `min_slip_depth` are the SSRM failure test. Extra keyword arguments
            pass through to `solve_ssrm` / `solve_fem`.

            Returns {'fem_data', 'solution', 'failure_solution', 'FS', 'analysis'}.
            A non-converged SSRM raises with the solver's own reason rather than
            returning a number that is not a factor of safety.
            """
            from xslope.fem import build_fem_data, solve_fem, solve_ssrm
            import matplotlib.pyplot as plt
            sd = doc.slope_data if slope_data is None else slope_data
            mesh = _ensure_mesh(sd)
            fem_data = build_fem_data(sd, mesh)
            if k0 is None:
                k0 = sd.get("k0")
            if analysis == "single":
                solution = solve_fem(fem_data, F=float(F), debug_level=1, k0=k0,
                                     **kwargs)
                bundle = {"fem_data": fem_data, "solution": solution,
                          "failure_solution": None, "FS": None, "analysis": "single"}
                print(f"FEM single trial (F = {float(F):g}): "
                      f"converged={solution.get('converged')}, "
                      f"iterations={solution.get('iterations')}")
            else:
                lo = sd.get("ssrm_f_min") if F_min is None else F_min
                hi = sd.get("ssrm_f_max") if F_max is None else F_max
                lo = 1.0 if lo is None else float(lo)
                hi = 2.0 if hi is None else float(hi)
                if tension_srf is None:
                    tension_srf = sd.get("tension_srf")
                result = solve_ssrm(fem_data, F_min=lo, F_max=hi,
                                    tolerance=float(tolerance), debug_level=1,
                                    failure_criterion=failure_criterion,
                                    min_slip_depth=min_slip_depth, k0=k0,
                                    tension_srf=tension_srf, **kwargs)
                if not result.get("converged", False):
                    raise RuntimeError(f"SSRM did not converge: "
                                       f"{result.get('error', 'unknown error')}")
                from xslope.fem import ssrm_run_record
                bundle = {"fem_data": fem_data, "solution": result["last_solution"],
                          "failure_solution": result.get("failure_solution"),
                          "FS": result.get("FS"), "analysis": "ssrm",
                          "meta": ssrm_run_record(result, fem_data, {})}
                print(f"SSRM: FS = {result['FS']:.3f}")
                if result.get("note"):
                    # The factor of safety is the bracket midpoint however the
                    # bracket closed, so an undecided upper edge is invisible in
                    # the number itself.
                    print(f"  {result['note']}")
            _store_result("fem_solution", bundle,
                          show=[("_show_fem_data", (fem_data,)),
                                ("_show_fem_results", ())])
            if plot:
                from xslope.plot_fem import plot_fem_results
                plot_fem_results(fem_data, bundle["solution"], fs=bundle["FS"],
                                 failure_solution=bundle.get("failure_solution"),
                                 fig=plt.figure(figsize=(11, 7)))
            return bundle

        def suggest_elastic(material_or_soil_type=None, unit_system=None,
                            slope_data=None):
            """Soil-type Young's modulus and Poisson's ratio for a material that
            carries none — the answer to "what E and nu should I use for this clay?".

            Thin wrapper over `xslope.units.classify_elastic`. Give it a material
            (its name, its 1-based row number, or the dict itself) and it classifies
            the material FROM ITS STRENGTH — a Hoek-Brown material is rock, a
            frictional material grades by phi, a phi = 0 material by its undrained
            shear strength, a c-phi soil on its cohesive component — and returns E in
            the MODEL'S OWN stress unit with the soil type it decided on. A soil-type
            name ('stiff clay', 'dense sand', 'soft rock') is also accepted, for a
            material that does not exist yet. With no argument it does every material
            in the model and prints the table.

            `unit_system` ('si' / 'metric' / 'imperial') overrides the model's
            declared system; without either, the unit system is inferred from the
            material unit weights.

            Returns a dict (or a list of them) with 'material', 'soil_type', 'E',
            'nu', 'unit_system' and 'reason'.

            THIS IS A LAST RESORT, never a preference. A value the problem states is
            an input: transcribe it. Use this only to fill a genuine blank — and SAY
            that you did, and which soil type it came from, because a blank E is a
            singular stiffness matrix while a wrong-magnitude E leaves the factor of
            safety untouched and corrupts every displacement the FEM reports.
            """
            from xslope import units as _units
            from xslope.units import (KPA_TO_PSF, classify_elastic,
                                      infer_unit_system, normalize_unit_system)
            sd = doc.slope_data if slope_data is None else slope_data
            materials = list(sd.get("materials") or [])

            system = normalize_unit_system(unit_system)
            if system is None:
                system = normalize_unit_system(sd.get("unit_system"))
            if system is None:
                system = infer_unit_system(sd.get("materials")) or "si"
            imperial = (system == "imperial")

            # The published per-soil-type table lives in xslope.units as one
            # constant per type; read it by value rather than restating it here, so
            # a revision to the table reaches this helper with no second edit.
            table = {}
            for name, value in vars(_units).items():
                if (name.isupper() and isinstance(value, tuple) and len(value) == 3
                        and isinstance(value[0], str)):
                    table[value[0].lower()] = value

            def _one(mat, label):
                soil, E, nu = classify_elastic(mat, imperial=imperial,
                                               declared_system=system)
                opt = str(mat.get("option", "mc") or "mc").lower()
                if opt in ("hb", "hoek", "hoek-brown"):
                    why = "a Hoek-Brown strength model means rock"
                else:
                    why = (f"c = {float(mat.get('c') or 0):g}, "
                           f"phi = {float(mat.get('phi') or 0):g}")
                return {"material": label, "soil_type": soil, "E": E, "nu": nu,
                        "unit_system": system,
                        "reason": (f"{label} classified as {soil} from {why}; E is the "
                                   f"midpoint of the published range for that soil "
                                   f"type and nu its typical value, in "
                                   f"{'psf' if imperial else 'kPa'}. Last-resort fill "
                                   f"— a stated value outranks it.")}

            if material_or_soil_type is None:
                rows = [_one(m, m.get("name") or f"Material {i + 1}")
                        for i, m in enumerate(materials)]
                for r in rows:
                    print(f"  {r['material']:<20} {r['soil_type']:<12} "
                          f"E = {r['E']:,.0f}   nu = {r['nu']:g}")
                return rows

            mat = None
            label = str(material_or_soil_type)
            if isinstance(material_or_soil_type, dict):
                mat = material_or_soil_type
                label = mat.get("name") or "material"
            elif isinstance(material_or_soil_type, int):
                mat = materials[material_or_soil_type - 1]      # 1-based, as the sheet
                label = mat.get("name") or f"Material {material_or_soil_type}"
            else:
                key = label.strip().lower()
                for m in materials:
                    if str(m.get("name", "")).strip().lower() == key:
                        mat = m
                        label = m.get("name")
                        break
                if mat is None:
                    row = table.get(key)
                    if row is None:
                        raise ValueError(
                            f"{material_or_soil_type!r} is neither a material in this "
                            f"model ({', '.join(str(m.get('name')) for m in materials) or 'none'}) "
                            f"nor a soil type ({', '.join(sorted(table))}).")
                    soil, e_kpa, nu = row
                    E = round(e_kpa * KPA_TO_PSF if imperial else e_kpa, -2)
                    out = {"material": None, "soil_type": soil, "E": E, "nu": nu,
                           "unit_system": system,
                           "reason": (f"Published range for {soil}: E is its midpoint "
                                      f"and nu its typical value, in "
                                      f"{'psf' if imperial else 'kPa'}. Last-resort "
                                      f"fill — a stated value outranks it.")}
                    print(f"  {soil}: E = {E:,.0f}, nu = {nu:g} ({system})")
                    return out
            out = _one(mat, label)
            print(f"  {out['material']}: {out['soil_type']} — E = {out['E']:,.0f}, "
                  f"nu = {out['nu']:g} ({system})")
            return out

        def generate_report(path=None, finalize=True, **options):
            """Build the Analysis Report — the same document File → Generate Report
            produces — and return the path it was written to.

            Runs `xslope.report.generate_report` over the model and everything this
            session has solved (the LEM run with the method it was run under, the
            seepage solutions in boundary-condition order, a transient march, the
            finite element run), then hands the finished .docx to whatever lays pages
            out (Word, or LibreOffice) so its contents page carries real page numbers
            — the dialog's own last step. A method the report is asked for that this
            session never ran is RUN by the builder, which is the longest thing in
            such a build: warn the user before asking for one.

            `path` defaults to `<model>_report.docx` beside the project (the home
            directory for a project never saved). `finalize=False` skips the page-
            number pass. Every other keyword is a report option, exactly as the
            Report dialog's checkboxes set them — `title`, `analyst`, and the section
            switches (`lem_slices`, `fem_piles`, …); `xslope.report.DEFAULT_OPTIONS`
            names them all. Returns the output path.
            """
            from xslope.report import generate_report as _generate
            from ..report_dialog import (default_output_path, document_finish,
                                         finalization_enabled)
            sd = doc.slope_data
            solutions = _report_solutions()
            model_path = getattr(doc, "path", None)
            out_path = path or default_output_path(model_path)
            opts = dict(options)
            # The traceability stamp names the input file and its SHA-256; without
            # this it reads "not saved to a file" even for a project opened from
            # disk. The Report dialog passes the same value.
            if model_path:
                opts.setdefault("input_path", model_path)
            style = getattr(doc, "style", None)
            if style is not None:
                opts.setdefault("style", style)
            ok, res = _generate(sd, solutions, opts, out_path)
            if not ok:
                raise RuntimeError(res)
            if finalize:
                settings = getattr(window, "settings", None)
                document_finish(out_path, finalization_enabled(settings))
            print(f"Report written to {res['path']} "
                  f"({len(res.get('figures') or ())} figures).")
            return res["path"]

        def _report_solutions():
            """What the report can document, in `xslope.report`'s shape. Studio's own
            assembly where there is a window (it also carries the method the LEM run
            was made under), and the same thing off `doc.results` where there is
            not."""
            builder = getattr(window, "report_solutions", None)
            if callable(builder):
                return builder()
            results = doc.results or {}
            out = {}
            lem = results.get("lem_solution")
            if lem:
                out["lem"] = [dict(lem, method=lem.get("method"))]
            seep = results.get("seep_solutions") or {}
            if seep:
                out["seep"] = [seep[bc] for bc in sorted(seep)]
            if results.get("transient_seep"):
                out["tseep"] = [results["transient_seep"]]
            if results.get("fem_solution"):
                out["fem"] = [results["fem_solution"]]
            return out

        def sensitivity(values, apply, param="value", method="spencer", search=True,
                        num_slices=40, rapid=False, name="sensitivity", plot=True,
                        ylabel="FS"):
            """Run a parametric sensitivity sweep and produce a summary CSV + ONE
            plot — the right tool for "FS vs <parameter>" studies.

            For each v in `values`: call apply(v) to edit slope_data, recompute the
            factor of safety, and record (v, FS) — with NO per-step solution plot.
            Writes `<name>.csv` and one `<name>.png` (FS vs param) to the output
            folder, restores the project to its original state (a sweep is analysis,
            not an edit), and returns the DataFrame.

            apply  : callable(v) that edits slope_data in place, setting the
                     parameter to an ABSOLUTE value (not an increment), e.g. to vary
                     the slope angle:
                         def set_angle(a):
                             import numpy as np
                             slope_data['profile_lines'][0]['coords'][1] = (
                                 20/np.tan(np.radians(a)), 20)
                         sensitivity(np.arange(20, 31), set_angle, param='angle (deg)')
            search : True recomputes the CRITICAL surface each step (correct but
                     slower); False uses the project's first circle / non-circ
                     surface for a quick look.
            """
            import copy
            import contextlib
            import io as _io
            import numpy as np
            import pandas as pd
            import matplotlib.pyplot as plt
            from xslope.slice import generate_slices
            from xslope.solve import solve_selected
            from xslope.search import circular_search, noncircular_search

            sd = doc.slope_data
            snapshot = copy.deepcopy(sd)
            rows = []
            try:
                for v in values:
                    apply(v)
                    resync_geometry(sd)
                    non_circ = sd.get("non_circ") or None
                    fs = float("nan")
                    with contextlib.redirect_stdout(_io.StringIO()):  # quiet the engine
                        if search:
                            if non_circ:
                                fs_cache, *_ = noncircular_search(
                                    sd, method, rapid=rapid, num_slices=num_slices,
                                    diagnostic=False)
                            else:
                                fs_cache, *_ = circular_search(
                                    sd, method, rapid=rapid, num_slices=num_slices,
                                    diagnostic=False)
                            if fs_cache:
                                sr = fs_cache[0].get("solver_result")
                                fs = sr["FS"] if isinstance(sr, dict) and "FS" in sr \
                                    else fs_cache[0].get("fs", float("nan"))
                        else:
                            circle = (sd["circles"][0] if sd.get("circular")
                                      and sd.get("circles") else None)
                            ok, res = generate_slices(sd, circle=circle,
                                                      non_circ=non_circ,
                                                      num_slices=num_slices)
                            if ok:
                                r = solve_selected(method, res[0], rapid=rapid)
                                if isinstance(r, dict):
                                    fs = r["FS"]
                    rows.append((v, fs))
                    print(f"  {param}={v}: FS={fs:.3f}")
            finally:
                sd.clear()
                sd.update(snapshot)     # leave the project exactly as it was

            df = pd.DataFrame(rows, columns=[param, ylabel])
            df.to_csv(f"{name}.csv", index=False)
            if plot:
                fig, ax = plt.subplots(figsize=(8, 5))
                ax.plot(df[param], df[ylabel], "o-")
                ax.set_xlabel(param)
                ax.set_ylabel(ylabel)
                ax.set_title(f"{ylabel} vs {param}")
                ax.grid(True, alpha=0.3)
                fig.savefig(f"{name}.png", dpi=150, bbox_inches="tight")
            print(f"Wrote {name}.csv and {name}.png ({len(df)} points).")
            return df

        def list_params(slope_data=None, mode="lem"):
            """List every sweepable parameter in the current project — the menu
            `design_sweep` and the engine sweepers draw their parameter refs from.

            Thin wrapper over `xslope.sensitivity.list_params`. Returns a list of
            dicts, one per parameter, each carrying a canonical `ref`
            ("kind:name:field", e.g. "mat:Clay:c", "global:k_seismic"), a `label`,
            the current `value`, and `sigma` (the reliability std-dev if the model
            carries one). In the default `mode='lem'` (also right for `'fem'`), covers
            every material's option-aware strength + general fields plus the global
            k_seismic. In `mode='seep'` the menu switches to the seepage-relevant set:
            each material's hydraulic fields (k1, k2, alpha, kr0, h0) plus every
            specified-head boundary value as a `seep_bc:<set>:<head_index>` ref.
            Blank/zero fields are still listed so a design sweep with explicit bounds
            can target them. Prints a compact table and returns the list so a snippet
            can pick a ref programmatically.

            Hand any `ref` straight to `design_sweep(param=...)`, or use the
            LLM-friendly dict form {'material': name_or_index, 'property': field} /
            {'global': field} / {'seep_bc': {'set': 1, 'head_index': 0}}.
            """
            from xslope.sensitivity import list_params as _list_params
            sd = doc.slope_data if slope_data is None else slope_data
            params = _list_params(sd, mode=mode)
            for p in params:
                val = "—" if p["value"] is None else f"{p['value']:g}"
                sig = f"   sigma={p['sigma']:g}" if p.get("sigma") else ""
                print(f"  {p['ref']:<30} = {val}{sig}")
            return params

        def design_sweep(param=None, low=None, high=None, steps=11, target_fs=1.5,
                         mode="lem", method="spencer", search=True, num_slices=40,
                         fem_opts=None, seep_opts=None, modify=None, label=None,
                         plot=True, slope_data=None):
            """Design sweep: vary ONE parameter from `low` to `high` and find the
            value at which the OUTPUT meets `target_fs` — the deterministic-design
            staple ("vary the undrained strength between X and Y, plot FS vs Su,
            highlight where FS = 1.5").

            Thin wrapper over `xslope.sensitivity.design`. `param` is a parameter
            reference in any form the engine accepts: a "kind:name:field" string
            (e.g. "mat:Clay:c", "global:k_seismic", "seep:Soil:k1",
            "seep_bc:1:0"), a (kind, name, field) tuple (name may be a 1-based
            material index), or the LLM-friendly dict {'material': name_or_index,
            'property': field} / {'global': field} / {'seep_bc': {'set': 1,
            'head_index': 0}}. Run `list_params(mode=...)` first to discover the refs.

            For anything that is NOT a single stored scalar — geometry above all (a
            slope angle, a berm width) — the design engine also takes the same
            `modify=` callable escape hatch as `sensitivity`: pass a Python
            `(slope_data, value) -> slope_data` function plus a `label` INSTEAD of
            `param` (exactly one of the two per call). Because a snippet defines the
            callable in Python, it crosses into this call directly (just like the
            `apply=` callable of `sensitivity`); the data-driven `param` references
            remain the natural form for a plain data-in/data-out request.

            `mode` selects the engine that evaluates each step and hence the OUTPUT
            quantity the sweep targets (`target_fs` names the target VALUE of that
            quantity, whatever it is):
              'lem'  — limit equilibrium; output = FS (default). `method` picks the
                       LEM method; `search` re-searches the critical surface per step.
              'fem'  — full SSRM solve per step (xslope.fem); output = FS. Needs a
                       finite-element mesh in slope_data['mesh'], and costs MINUTES
                       PER STEP — keep `steps` small. `fem_opts` forwards SSRM knobs
                       (F_min, F_max, tolerance, failure_criterion, min_slip_depth).
              'seep' — seepage solve per step (xslope.seep); output = total discharge
                       q (so `target_fs` is a target q, e.g. 6e-6). Needs a mesh;
                       `seep_opts` = {'bc': 1|2, 'tol': ...}. The classic study varies
                       a hydraulic k (seep:<mat>:k1) or a reservoir head
                       (seep_bc:<set>:<i>) and finds the value giving a target q.

            Runs `steps` evenly spaced solves across [low, high], re-searching the
            critical surface at each step by default in LEM (search=True — the
            critical surface MOVES as the parameter changes; a fixed-surface sweep
            understates the effect). Renders ONE output-vs-parameter plot (via
            `plot_sensitivity`, which labels the axes from the df's output_label — FS
            with an FS = 1 guide, or q — and draws the target guide line) the same way
            `run_lem` shows its plot — pass plot=False in a loop. Does NOT modify the
            project (a sweep is analysis, not an edit).

            Returns the engine's result dict:
              'crossing'  — interpolated parameter value at output = target, or None.
              'bracketed' — True iff output = target is crossed inside [low, high].
              'crossings' — every crossing found (list; usually one).
              'fs_range'  — (min, max) OUTPUT over the successful sweep points.
              'direction' — 'increasing' / 'decreasing' / 'non-monotonic'.
              'extend'    — when NOT bracketed, which way to widen the range
                            ('above {high}' / 'below {low}'). Report fs_range + this;
                            NEVER extrapolate a crossing past the swept range.
              'output' / 'output_label' — 'FS'/'q' and the axis label.
              'message'   — one-line human summary. Also: 'df', 'param', 'base_value'.
            """
            from xslope.sensitivity import design as _design
            from xslope.plot import plot_sensitivity
            import matplotlib.pyplot as plt
            sd = doc.slope_data if slope_data is None else slope_data
            if mode in ("fem", "seep") and sd.get("mesh") is None:
                raise RuntimeError(
                    f"mode='{mode}' evaluates each step with the finite-element "
                    f"{'SSRM' if mode == 'fem' else 'seepage'} solver, which needs a "
                    "mesh in slope_data['mesh'] — build a mesh first (Studio: the "
                    "Mesh toolbar/menu, or xslope.mesh.build_mesh_from_polygons).")
            ok, res = _design(sd, param=param, low=low, high=high, steps=steps,
                              target_fs=target_fs, mode=mode, method=method,
                              search=search, num_slices=num_slices, fem_opts=fem_opts,
                              seep_opts=seep_opts, modify=modify, label=label)
            if not ok:
                raise RuntimeError(res)
            out = res.get("output", "FS")
            print(res["message"])
            if not res["bracketed"]:
                lo, hi = res["fs_range"]
                print(f"  (not bracketed — {out} spans [{lo:.3g}, {hi:.3g}]; extend "
                      f"the range {res['extend']} to reach {out} = {target_fs:g})")
            if plot:
                plot_sensitivity(res["df"], target_fs=res["target_fs"],
                                 fig=plt.figure(figsize=(8, 5)))
            return res

        def parametric_back_analysis(param=None, low=None, high=None, steps=11,
                                     target_fs=1.0, mode="lem", method="spencer",
                                     search=True, num_slices=40, fem_opts=None,
                                     seep_opts=None, modify=None, label=None,
                                     plot=True, slope_data=None):
            """Forensic back-analysis: vary ONE parameter and find the value that makes
            the slope limiting (FS = 1.0 by default) — the Back-Analysis mode of the
            Parametric study family.

            A slide has occurred, so FS at failure is known to be 1.0 and the unknown
            is a strength / pore-pressure / loading parameter; this inverts the problem
            to the value CONSISTENT WITH the observed failure. Mechanically it is
            `parametric_design` with `target_fs=1.0`, so all arguments and the returned
            result dict match `design_sweep` (`crossing` is the back-calculated value).
            Thin wrapper over `xslope.sensitivity.back_analysis`. It inherits the same
            `modify=`/`label=` callable escape hatch as `design_sweep` (exclusive with
            `param`), for back-calculating a non-scalar unknown — a water-table
            elevation, say.
            """
            from xslope.sensitivity import back_analysis as _back
            from xslope.plot import plot_sensitivity
            import matplotlib.pyplot as plt
            sd = doc.slope_data if slope_data is None else slope_data
            if mode in ("fem", "seep") and sd.get("mesh") is None:
                raise RuntimeError(f"mode='{mode}' needs a mesh in slope_data['mesh'].")
            ok, res = _back(sd, param=param, low=low, high=high, steps=steps,
                            target_fs=target_fs, mode=mode, method=method,
                            search=search, num_slices=num_slices, fem_opts=fem_opts,
                            seep_opts=seep_opts, modify=modify, label=label)
            if not ok:
                raise RuntimeError(res)
            print(res["message"])
            if not res["bracketed"]:
                out = res.get("output", "FS")
                lo, hi = res["fs_range"]
                print(f"  (not bracketed — {out} spans [{lo:.3g}, {hi:.3g}]; extend "
                      f"the range {res['extend']} to reach {out} = {target_fs:g})")
            if plot:
                plot_sensitivity(res["df"], target_fs=res["target_fs"],
                                 fig=plt.figure(figsize=(8, 5)))
            return res

        def parametric_sweep(params=None, plot="scaled", scaling="elasticity",
                             mode="lem", method="spencer", rel_range=0.25, n=7,
                             num_slices=40, search=True, n_samples=None,
                             fem_opts=None, seep_opts=None, plot_fig=True,
                             slope_data=None):
            """Parametric SENSITIVITY across several parameters at once, rendering one
            of the plot-family views — the multi-parameter 'sweep call' of the
            Parametric study family (its single-parameter siblings are `design_sweep` /
            `parametric_design` and `parametric_back_analysis`).

            `params` is a list of parameter refs (any form `list_params` /
            `design_sweep` accept); None auto-selects every material strength/weight
            parameter in the model. `plot` picks the view — every plot type in the
            family is reachable from this one call:
              'tornado'  — FS swing per parameter between its ±`rel_range` bounds
                           (Duncan tornado; `plot_tornado`).
              'scaled'   — scaled-sensitivity bars (`plot_scaled_sensitivity`), with
                           `scaling` one of 'elasticity' (default, ∂F/∂p·p/F),
                           'per_1pct' (ΔFS per 1%), 'per_sigma' (ΔFS per σ; σ params
                           only). The derivative is a central difference at ±1%.
              'spider'   — FS vs each parameter over ±`rel_range`, `n` points, on a
                           normalized % x-axis (`plot_spider`).
              'variance' — variance-contribution Pareto (`plot_variance_pareto`);
                           needs σ; reuses the Taylor-series reliability machinery.
              'rank'     — Monte Carlo Spearman rank correlation with FS
                           (`plot_mc_rank_correlation`); needs σ; `n_samples` sizes the
                           campaign (default 10000). A GLOBAL measure, complementary to
                           the LOCAL 'scaled' bars.

            `mode` selects the engine ('lem' default / 'fem' / 'seep'), exactly as in
            `design_sweep`. Renders ONE figure and returns the underlying engine result
            dict (a tornado result, a scaled/variance/rank result, or `{'sweeps': ...}`
            for the spider). Does NOT modify the project.
            """
            from xslope.sensitivity import (sensitivity as _sens, tornado_from_sweeps,
                                            scaled_sensitivity, variance_contribution,
                                            mc_rank_correlation, list_params as _lp)
            from xslope import plot as _plot
            import matplotlib.pyplot as plt
            sd = doc.slope_data if slope_data is None else slope_data
            if params is None:
                params = [e["ref"] for e in _lp(sd, mode="seep" if mode == "seep"
                                                else "lem")
                          if e.get("value") is not None
                          and e["kind"] in ("mat", "seep")]
            if isinstance(params, (str, dict, tuple)):
                params = [params]
            if mode in ("fem", "seep") and sd.get("mesh") is None:
                raise RuntimeError(f"mode='{mode}' needs a mesh in slope_data['mesh'].")

            def newfig():
                return plt.figure(figsize=(8, 5)) if plot_fig else None

            if plot == "variance":
                ok, res = variance_contribution(sd, method=method, search=search)
                if not ok:
                    raise RuntimeError(res)
                if plot_fig:
                    _plot.plot_variance_pareto(res, fig=newfig())
                return res
            if plot == "rank":
                kw = {"n_samples": int(n_samples)} if n_samples else {}
                ok, res = mc_rank_correlation(sd, method=method, search=search,
                                              num_slices=num_slices, **kw)
                if not ok:
                    raise RuntimeError(res)
                if plot_fig:
                    _plot.plot_mc_rank_correlation(res, fig=newfig())
                return res
            if plot == "scaled":
                ok, res = scaled_sensitivity(sd, params, method=method, mode=mode,
                                             search=search, num_slices=num_slices,
                                             fem_opts=fem_opts, seep_opts=seep_opts)
                if not ok:
                    raise RuntimeError(res)
                if plot_fig:
                    _plot.plot_scaled_sensitivity(res, scaling=scaling, fig=newfig())
                return res
            # tornado / spider both need the full per-parameter sweeps
            sweeps = {}
            for ref in params:
                ok, res = _sens(sd, param=ref, rel_range=rel_range, n=n, mode=mode,
                                methods=(method,), search=search, num_slices=num_slices,
                                fem_opts=fem_opts, seep_opts=seep_opts)
                if not ok:
                    raise RuntimeError(f"{ref}: {res}")
                sweeps[res["param"]] = res["df"]
            if plot == "spider":
                if plot_fig:
                    _plot.plot_spider(sweeps, fig=newfig())
                return {"sweeps": sweeps}
            tor = tornado_from_sweeps(sweeps, method=method)
            if plot_fig:
                _plot.plot_tornado(tor, fig=newfig())
            return tor

        # Preferred names under the Parametric umbrella; `design_sweep` is kept as a
        # back-compatible alias (older skill recipes and saved snippets call it).
        parametric_design = design_sweep

        def reliability_taylor(method="bishop", rapid=False, circular=True,
                               search=True, plot=True, slope_data=None, **kwargs):
            """Taylor Series Probability Method (TSPM) reliability — the probabilistic
            counterpart to the deterministic Parametric study.

            Turns the material standard deviations (the s(·) columns) into a
            reliability index and probability of failure by ``1 + 2N`` limit-
            equilibrium analyses: the factor of safety at the most-likely values, then
            a ±σ perturbation of every uncertain parameter. Thin wrapper over
            `xslope.reliability.reliability_taylor`. `method` is any LEM solver; extra
            keyword arguments (e.g. `fs_tol`, `tol`, `max_iter`, `composite`, `seed`)
            pass straight through. Renders the MLV / F± surface plot when `plot=True`.
            Returns the reliability result dict (`F_MLV`, `sigma_F`, `COV_F`,
            `beta_ln`, `reliability`, `prob_failure`, `param_info`, …).
            """
            from xslope.reliability import reliability_taylor as _rt
            from xslope.plot import plot_reliability_results
            import matplotlib.pyplot as plt
            sd = doc.slope_data if slope_data is None else slope_data
            ok, res = _rt(sd, method, rapid=rapid, circular=circular, search=search,
                          **kwargs)
            if not ok:
                raise RuntimeError(res)
            if plot and res.get("fs_cache"):
                plot_reliability_results(sd, res, fig=plt.figure(figsize=(11, 6)))
            return res

        def reliability_mc(method="bishop", n_samples=None, rng_seed=None,
                           distribution="normal", rapid=False, circular=True,
                           search=True, num_slices=40, plot=True, slope_data=None,
                           **kwargs):
            """Monte Carlo reliability — draws `n_samples` realizations of every
            uncertain parameter from its standard deviation, evaluates the factor of
            safety of each on a FIXED failure surface, and reports the sample
            statistics and an FS histogram (limit-equilibrium only). Thin wrapper over
            `xslope.reliability.reliability_mc`. `rng_seed` is a fixed constant by
            default (reproducible). Renders the FS histogram when `plot=True`. Returns
            the result dict (`mean_FS`, `sigma_F`, `beta_normal`, `beta_ln`,
            `pf_empirical`, `fs_samples`, …).
            """
            from xslope.reliability import reliability_mc as _rmc
            from xslope.plot import plot_reliability_histogram
            import matplotlib.pyplot as plt
            sd = doc.slope_data if slope_data is None else slope_data
            kw = dict(kwargs)
            if n_samples is not None:
                kw["n_samples"] = int(n_samples)
            if rng_seed is not None:
                kw["rng_seed"] = int(rng_seed)
            ok, res = _rmc(sd, method, rapid=rapid, circular=circular, search=search,
                           distribution=distribution, num_slices=num_slices, **kw)
            if not ok:
                raise RuntimeError(res)
            if plot:
                plot_reliability_histogram(res, fig=plt.figure(figsize=(9, 5.5)))
            return res

        def reliability_rs(method="bishop", rng_seed=None, distribution="normal",
                           rapid=False, circular=True, search=True, num_slices=40,
                           plot=True, slope_data=None, **kwargs):
            """Response-surface reliability — fits a quadratic surrogate to a central
            composite design of real solves, measures it against held-out real solves
            and against the realizations it counts as failures, then samples it ten
            million times (limit-equilibrium only). Thin wrapper over
            `xslope.reliability.reliability_rs`. A surrogate that fails its accuracy
            gate raises with the measured errors rather than returning a degraded
            answer. Renders the FS histogram when `plot=True`. Returns the result dict
            (the Monte Carlo keys plus `gate_r2`, `gate_rms`, `n_design`, …).
            """
            from xslope.reliability import reliability_rs as _rrs
            from xslope.plot import plot_reliability_histogram
            import matplotlib.pyplot as plt
            sd = doc.slope_data if slope_data is None else slope_data
            kw = dict(kwargs)
            if rng_seed is not None:
                kw["rng_seed"] = int(rng_seed)
            ok, res = _rrs(sd, method, rapid=rapid, circular=circular, search=search,
                           distribution=distribution, num_slices=num_slices, **kw)
            if not ok:
                raise RuntimeError(res)
            if plot:
                plot_reliability_histogram(res, fig=plt.figure(figsize=(9, 5.5)))
            return res

        def reliability(method="bishop", engine="taylor", **kwargs):
            """Front door to the reliability family (mirrors
            `xslope.reliability.reliability`): `engine='taylor'` (default) runs the
            Taylor Series Probability Method, `engine='mc'` runs Monte Carlo, and
            `engine='rs'` runs the response surface. Kept under the plain name for
            back-compat; prefer `reliability_taylor` / `reliability_mc` /
            `reliability_rs` directly.
            """
            key = str(engine).lower().replace("-", "_")
            if key in ("mc", "monte_carlo", "montecarlo"):
                return reliability_mc(method=method, **kwargs)
            if key in ("rs", "rsm", "response_surface", "responsesurface"):
                return reliability_rs(method=method, **kwargs)
            return reliability_taylor(method=method, **kwargs)

        def corpus_index(query=None, limit=10):
            """Worked examples from the verification corpus, matching `query`.

            The corpus is 300-odd solved models, each carrying a published
            comparison against the source or the vendor program — the pages to cite
            when a question matches a topic. The whole table is far too big to sit
            in a system prompt on every completion, so it is asked for instead:
            `corpus_index('rapid drawdown')` returns the matching rows as
            `[{'topic', 'title', 'url'}]` (also printed), and `corpus_index()` with
            no query lists the topics.

            Matching is on topic key, topic label, and title, so a plain-English
            phrase works. `limit` caps the rows so a broad query cannot flood the
            conversation.
            """
            topics, all_rows = _corpus_rows()
            if not all_rows:
                print("The verification corpus index is not available in this "
                      "install; the pages are at "
                      "https://xslope.readthedocs.io/en/latest/verification/")
                return []
            if not query:
                print("Corpus topics (pass one to corpus_index):")
                for key in sorted(topics):
                    print(f"  {key:22s} {topics[key]}")
                return sorted(topics)
            words = [w for w in str(query).lower().split() if len(w) > 2]
            rows, seen = [], set()
            for row in all_rows:
                text = f"{row['topic']} {row['label']} {row['title']}".lower()
                if not all(w in text for w in words):
                    continue
                if row["url"] in seen:
                    continue
                seen.add(row["url"])
                rows.append({"topic": row["topic"], "title": row["title"],
                             "url": row["url"]})
                if len(rows) >= int(limit):
                    break
            if not rows:
                print(f"No corpus entry matches {query!r}. Try a topic from "
                      "corpus_index().")
                return []
            for r in rows:
                print(f"- {r['title']}  [{r['topic']}]\n  {r['url']}")
            return rows

        return {"run_lem": run_lem, "run_seep": run_seep, "run_fem": run_fem,
                "corpus_index": corpus_index,
                "suggest_elastic": suggest_elastic,
                "generate_report": generate_report,
                "resync_geometry": resync_geometry,
                "sensitivity": sensitivity, "list_params": list_params,
                "design_sweep": design_sweep,
                "parametric_sweep": parametric_sweep,
                "parametric_design": parametric_design,
                "parametric_back_analysis": parametric_back_analysis,
                "reliability": reliability,
                "reliability_taylor": reliability_taylor,
                "reliability_mc": reliability_mc,
                "reliability_rs": reliability_rs}

    @staticmethod
    def _normalize(code):
        """Salvage code from weaker models that emit literal ``\\n`` / ``\\t``
        escape sequences instead of real newlines (so it lands as one unparseable
        line). Only applied when the code won't compile as-is *and* the unescaped
        version does — well-formed code (incl. real backslashes in strings) is
        left untouched."""
        try:
            compile(code, "<assistant>", "exec")
            return code
        except SyntaxError:
            if "\\n" in code or "\\t" in code:
                fixed = code.replace("\\r\\n", "\n").replace("\\n", "\n").replace("\\t", "\t")
                try:
                    compile(fixed, "<assistant>", "exec")
                    return fixed
                except SyntaxError:
                    pass
        return code

    def run(self, code):
        """Execute ``code``; return ``(stdout_text, [output_file_paths], error_text)``.
        Output files are everything the snippet produced in ``outdir`` — pyplot
        figures it left open (auto-saved) plus any files it wrote itself (plots,
        CSVs, …). ``error_text`` is ``None`` on success."""
        if not self._seeded:
            self._seed()
        import matplotlib.pyplot as plt

        # Refresh the live references each call (the document dict is replaced on
        # New / Open, so a stale binding would point at the old project).
        self._ns["doc"] = self._doc
        self._ns["slope_data"] = self._doc.slope_data
        self._ns["results"] = self._doc.results
        self._ns["OUTPUT_DIR"] = self.outdir

        code = self._normalize(code)
        # This snippet's own bookkeeping: what it solves (so Studio's stale-result
        # sweep cannot take a run made after the edit that triggered it) and the
        # geometry it starts from (so an edit made on the source the resync
        # overwrites is named rather than silently reverted).
        self._fresh_results = {}
        self._geom_warnings = []
        sd = self._doc.slope_data
        self._geom_watch = ((_geometry_native(sd), _geometry_sigs(sd))
                            if isinstance(sd, dict) else None)
        os.makedirs(self.outdir, exist_ok=True)     # may have been cleared by the OS
        before_figs = set(plt.get_fignums())
        existing = set(os.listdir(self.outdir))     # to detect files the snippet writes
        prev_cwd = os.getcwd()
        buf = StringIO()
        error = None
        try:
            os.chdir(self.outdir)        # relative saves (savefig('x.png')) land here
            with redirect_stdout(buf), redirect_stderr(buf):
                exec(code, self._ns)
        except Exception:
            error = traceback.format_exc()
        finally:
            try:
                os.chdir(prev_cwd)
            except Exception:
                pass
        if error:
            # A snippet that raised is rolled back whole, so its geometry edit is
            # not the model's problem and its results were never the session's.
            self._geom_warnings = []
            self._fresh_results = {}
        else:
            # A snippet that ran nothing still has its polygons rebuilt — by the
            # checks' own resync, after this returns — so the last word on the
            # source of a geometry edit is said here.
            self._note_geometry_source()
        self._geom_watch = None

        # Auto-save figures the snippet left open — but only if it didn't already
        # save a figure itself (e.g. savefig('plot.png') without closing), so a
        # saved-and-left-open plot isn't shown twice.
        wrote_figure = any(f.lower().endswith(_FIG_EXTS)
                           for f in os.listdir(self.outdir) if f not in existing)
        if not wrote_figure:
            for num in plt.get_fignums():
                if num in before_figs:
                    continue
                self._fig_seq += 1
                path = os.path.join(self.outdir, f"assistant_plot_{self._fig_seq:03d}.png")
                try:
                    plt.figure(num).savefig(path, dpi=120, bbox_inches="tight")
                except Exception:
                    pass
        plt.close("all")

        # Everything newly present in outdir is this snippet's output.
        created = sorted(f for f in os.listdir(self.outdir) if f not in existing)
        outputs = [os.path.join(self.outdir, f) for f in created]
        return buf.getvalue(), outputs, error
