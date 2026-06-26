"""Background run controllers (Phase 3).

Long engine calls run on a ``QThread`` so the UI stays responsive. ``LemRunner``
performs a single-surface circular LEM solve and emits the result (or an error)
back to the GUI thread. Engine ``print`` output streams to the Log pane live via
the thread-safe stdout tee installed by the main window (``_LogStream``), so the
worker does not touch any Qt widget itself.
"""

from __future__ import annotations

import traceback

from PySide6.QtCore import QThread, Signal


class LemRunner(QThread):
    """Runs a single circular-surface LEM solve off the GUI thread.

    Emits ``succeeded`` with ``{slice_df, failure_surface, results}`` or
    ``failed`` with the error message.
    """

    succeeded = Signal(object)
    failed = Signal(str)

    def __init__(self, slope_data, method, num_slices, rapid=False, parent=None):
        super().__init__(parent)
        self._slope_data = slope_data
        self._method = method
        self._num_slices = num_slices
        self._rapid = rapid

    def run(self):
        from xslope.slice import generate_slices
        from xslope.solve import solve_selected

        try:
            sd = self._slope_data
            circle = sd["circles"][0] if sd.get("circular") and sd.get("circles") else None
            if circle is None:
                self.failed.emit(
                    "A circular failure surface is required (no circles defined).")
                return
            ok, result = generate_slices(sd, circle=circle, non_circ=None,
                                         num_slices=self._num_slices)
            if not ok:
                self.failed.emit(str(result))
                return
            slice_df, failure_surface = result
            results = solve_selected(self._method, slice_df, rapid=self._rapid)
            if not isinstance(results, dict):
                self.failed.emit(f"No solution: {results}")
                return
            self.succeeded.emit({
                "slice_df": slice_df,
                "failure_surface": failure_surface,
                "results": results,
            })
        except Exception:
            traceback.print_exc()   # streams to the Log pane via the stdout/stderr tee
            self.failed.emit("Solve failed — see the Log pane for details.")
