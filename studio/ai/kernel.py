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

import tempfile
import traceback
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO


class PythonKernel:
    def __init__(self, doc):
        self._doc = doc
        self._ns = {}
        self._seeded = False

    def _seed(self):
        import matplotlib
        matplotlib.use("Agg")          # pyplot figures render off-screen; we save them
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import xslope
        self._ns.update({"np": np, "pd": pd, "plt": plt, "xslope": xslope})
        self._seeded = True

    def run(self, code):
        """Execute ``code``; return ``(stdout_text, [figure_png_paths], error_text)``.
        ``error_text`` is ``None`` on success."""
        if not self._seeded:
            self._seed()
        import matplotlib.pyplot as plt

        # Refresh the live references each call (the document dict is replaced on
        # New / Open, so a stale binding would point at the old project).
        self._ns["doc"] = self._doc
        self._ns["slope_data"] = self._doc.slope_data
        self._ns["results"] = self._doc.results

        before = set(plt.get_fignums())
        buf = StringIO()
        error = None
        try:
            with redirect_stdout(buf), redirect_stderr(buf):
                exec(code, self._ns)
        except Exception:
            error = traceback.format_exc()

        # Save any figures the snippet created (new since `before`).
        figures = []
        for num in plt.get_fignums():
            if num in before:
                continue
            fig = plt.figure(num)
            path = tempfile.NamedTemporaryFile(
                prefix="xslope_ai_", suffix=".png", delete=False).name
            try:
                fig.savefig(path, dpi=120, bbox_inches="tight")
                figures.append(path)
            except Exception:
                pass
        plt.close("all")

        return buf.getvalue(), figures, error
