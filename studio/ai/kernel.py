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


class PythonKernel:
    def __init__(self, doc):
        self._doc = doc
        self._ns = {}
        self._seeded = False
        self._outdir = None
        self._fig_seq = 0

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

        - ``run_lem(method='bishop', ...)`` — one single-surface LEM solve on the
          loaded surface; returns the result dict, shows the solution plot.
        - ``resync_geometry(slope_data=None)`` — rebuild derived geometry after an
          in-snippet geometry edit (call inside sweep loops).
        - ``sensitivity(values, apply, param=..., ...)`` — callback-driven FS-vs-
          parameter sweep (you write ``apply(v)``); writes a CSV + one plot.
        - ``list_params(slope_data=None)`` — discover every sweepable parameter as
          a menu of canonical refs (feeds design_sweep / the engine sweepers).
        - ``design_sweep(param, low, high, target_fs=1.5, ...)`` — vary ONE named
          parameter and find where FS meets a target; renders one FS-vs-parameter
          plot with the target highlighted. The engine-driven counterpart to
          ``sensitivity`` — no callback needed, just a parameter ref or dict spec.
        """
        doc = self._doc

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
            _resync_geometry(sd)
            return sd.get("ground_surface")

        def run_lem(method="bishop", num_slices=40, rapid=False, plot=True,
                    slope_data=None):
            """Run a single-surface LEM analysis on the loaded project's failure
            surface and return the result dict (includes 'FS'). `method` is one of
            oms, bishop, janbu, spencer, corps, lowe, mprice. Shows the solution
            plot when plot=True (pass plot=False in sweeps to avoid many figures).
            Rebuilds derived geometry first, so edits to profile_lines/polygons
            this snippet made are reflected."""
            from xslope.slice import generate_slices
            from xslope.solve import solve_selected
            sd = doc.slope_data if slope_data is None else slope_data
            resync_geometry(sd)        # reflect any in-snippet geometry edits
            circle = (sd["circles"][0] if sd.get("circular") and sd.get("circles")
                      else None)
            non_circ = sd.get("non_circ") or None
            if circle is None and not non_circ:
                raise RuntimeError("No failure surface defined (add a circle or "
                                   "non-circular surface first).")
            ok, res = generate_slices(sd, circle=circle, non_circ=non_circ,
                                      num_slices=num_slices)
            if not ok:
                raise RuntimeError(res)
            slice_df, surface = res
            result = solve_selected(method, slice_df, rapid=rapid)
            if not isinstance(result, dict):
                raise RuntimeError(f"No solution: {result}")
            print(f"{method}: FS = {result['FS']:.3f}")
            if plot:
                from xslope.plot import plot_solution
                import matplotlib.pyplot as plt
                plot_solution(sd, slice_df, surface, result,
                              fig=plt.figure(figsize=(11, 6)))
            return result

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

        def list_params(slope_data=None):
            """List every sweepable parameter in the current project — the menu
            `design_sweep` and the engine sweepers draw their parameter refs from.

            Thin wrapper over `xslope.sensitivity.list_params`. Returns a list of
            dicts, one per parameter, each carrying a canonical `ref`
            ("kind:name:field", e.g. "mat:Clay:c", "global:k_seismic"), a `label`,
            the current `value`, and `sigma` (the reliability std-dev if the model
            carries one). Covers every material's option-aware strength + general
            fields plus the global k_seismic; blank/zero fields are still listed so a
            design sweep with explicit bounds can target them. Prints a compact table
            and returns the list so a snippet can pick a ref programmatically.

            Hand any `ref` straight to `design_sweep(param=...)`, or use the
            LLM-friendly dict form {'material': name_or_index, 'property': field} /
            {'global': field}.
            """
            from xslope.sensitivity import list_params as _list_params
            sd = doc.slope_data if slope_data is None else slope_data
            params = _list_params(sd)
            for p in params:
                val = "—" if p["value"] is None else f"{p['value']:g}"
                sig = f"   sigma={p['sigma']:g}" if p.get("sigma") else ""
                print(f"  {p['ref']:<30} = {val}{sig}")
            return params

        def design_sweep(param, low, high, steps=11, target_fs=1.5,
                         method="spencer", search=True, num_slices=40, plot=True,
                         slope_data=None):
            """Design sweep: vary ONE parameter from `low` to `high` and find the
            value at which the factor of safety meets `target_fs` — the
            deterministic-design staple ("vary the undrained strength between X and
            Y, plot FS vs Su, highlight where FS = 1.5").

            Thin wrapper over `xslope.sensitivity.design`. `param` is a parameter
            reference in any form the engine accepts: a "kind:name:field" string
            (e.g. "mat:Clay:c", "global:k_seismic"), a (kind, name, field) tuple
            (name may be a 1-based material index), or the LLM-friendly dict
            {'material': name_or_index, 'property': field} / {'global': field}. Run
            `list_params()` first to discover the refs.

            Runs `steps` evenly spaced solves across [low, high], re-searching the
            critical surface at each step by default (search=True — the critical
            surface MOVES as the parameter changes; a fixed-surface sweep understates
            the effect). Renders ONE FS-vs-parameter plot (via `plot_sensitivity`,
            with FS = 1 and FS = target_fs guide lines) the same way `run_lem` shows
            its plot — pass plot=False in a loop. Does NOT modify the project (a sweep
            is analysis, not an edit).

            Returns the engine's result dict:
              'crossing'  — interpolated parameter value at FS = target_fs, or None.
              'bracketed' — True iff FS = target_fs is crossed inside [low, high].
              'crossings' — every crossing found (list; usually one).
              'fs_range'  — (min FS, max FS) over the successful sweep points.
              'direction' — 'increasing' / 'decreasing' / 'non-monotonic'.
              'extend'    — when NOT bracketed, which way to widen the range
                            ('above {high}' / 'below {low}'). Report fs_range + this;
                            NEVER extrapolate a crossing past the swept range.
              'message'   — one-line human summary. Also: 'df', 'param', 'base_value'.
            """
            from xslope.sensitivity import design as _design
            from xslope.plot import plot_sensitivity
            import matplotlib.pyplot as plt
            sd = doc.slope_data if slope_data is None else slope_data
            ok, res = _design(sd, param, low, high, steps=steps, target_fs=target_fs,
                              method=method, search=search, num_slices=num_slices)
            if not ok:
                raise RuntimeError(res)
            print(res["message"])
            if not res["bracketed"]:
                lo, hi = res["fs_range"]
                print(f"  (not bracketed — FS spans [{lo:.3f}, {hi:.3f}]; extend "
                      f"the range {res['extend']} to reach FS = {target_fs:g})")
            if plot:
                plot_sensitivity(res["df"], target_fs=res["target_fs"],
                                 fig=plt.figure(figsize=(8, 5)))
            return res

        return {"run_lem": run_lem, "resync_geometry": resync_geometry,
                "sensitivity": sensitivity, "list_params": list_params,
                "design_sweep": design_sweep}

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
