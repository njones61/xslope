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
        have to reconstruct the engine pipeline (a common failure mode)."""
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

        return {"run_lem": run_lem, "resync_geometry": resync_geometry}

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
