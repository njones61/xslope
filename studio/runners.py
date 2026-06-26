"""Background run controllers (Phase 3).

Long engine calls run on a ``QThread`` so the UI stays responsive. ``LemRunner``
performs a single-surface circular LEM solve and emits the result (or an error)
back to the GUI thread. Engine stdout is captured during the run and delivered
with the result, so the GUI thread — not the worker — writes to the Log pane
(Qt widgets are not safe to touch off the GUI thread).
"""

from __future__ import annotations

import io
import traceback
from contextlib import redirect_stdout

from PySide6.QtCore import QThread, Signal


class LemRunner(QThread):
    """Runs a single circular-surface LEM solve off the GUI thread.

    Emits ``succeeded`` with ``{slice_df, failure_surface, results, log}`` or
    ``failed`` with ``(message, log)``.
    """

    succeeded = Signal(object)
    failed = Signal(str, str)

    def __init__(self, slope_data, method, num_slices, rapid=False, parent=None):
        super().__init__(parent)
        self._slope_data = slope_data
        self._method = method
        self._num_slices = num_slices
        self._rapid = rapid

    def run(self):
        from xslope.slice import generate_slices
        from xslope.solve import solve_selected

        buf = io.StringIO()
        try:
            with redirect_stdout(buf):
                sd = self._slope_data
                circle = sd["circles"][0] if sd.get("circular") and sd.get("circles") else None
                if circle is None:
                    self.failed.emit(
                        "A circular failure surface is required (no circles defined).",
                        buf.getvalue())
                    return
                ok, result = generate_slices(sd, circle=circle, non_circ=None,
                                             num_slices=self._num_slices)
                if not ok:
                    self.failed.emit(str(result), buf.getvalue())
                    return
                slice_df, failure_surface = result
                results = solve_selected(self._method, slice_df, rapid=self._rapid)
            if not isinstance(results, dict):
                self.failed.emit(f"No solution: {results}", buf.getvalue())
                return
            self.succeeded.emit({
                "slice_df": slice_df,
                "failure_surface": failure_surface,
                "results": results,
                "log": buf.getvalue(),
            })
        except Exception as exc:
            traceback.print_exc(file=buf)
            self.failed.emit(str(exc), buf.getvalue())
