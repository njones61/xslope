"""Transient-seepage results view: a seep-solution canvas plus a play bar.

A transient seepage run produces a sequence of saved frames — each a full seepage
solution (head / pore pressure / stream function, with velocities derived on load)
at one instant. This view shows one frame at a time through the SAME
``plot_seep_solution`` path the steady solution tab uses (so contours, flow lines,
velocity vectors, the phreatic surface and the colorbar all render exactly as they
do for a steady solve), and adds a play bar along the bottom of the plot window
(plan §7): transport buttons, a frame slider with a time readout that doubles as a
jump-to-time entry, a playback-speed selector, and two buttons that tag the current
time as the rapid-drawdown ``stage_1`` / ``stage_2`` time.

The per-frame ``t = {value} {unit}`` title annotation is produced by
``plot_seep_solution`` itself (it reads the frame's ``time`` key and the model's
declared time unit) — this view supplies no title text of its own.
"""

from __future__ import annotations

import traceback

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QComboBox, QHBoxLayout, QLabel, QLineEdit, QPushButton, QSlider,
    QToolButton, QVBoxLayout, QWidget,
)

from .canvas import MplCanvas

# Playback speed multipliers offered in the speed selector, and the per-frame dwell
# time at 1x (halved at 2x, etc.). Kept modest so a ~50-frame run plays in a few
# seconds without being a blur.
_SPEEDS = [("0.5x", 0.5), ("1x", 1.0), ("2x", 2.0), ("4x", 4.0)]
_BASE_INTERVAL_MS = 500


class TransientSeepView(QWidget):
    """Seep-solution tab for a transient run: one frame shown at a time with a
    play bar underneath. Register the whole widget as the tab (it forwards
    ``ensure_fitted`` so the standard tab-shown fit path works) and its inner
    :class:`MplCanvas` does the rendering."""

    # (stage number 1 | 2, current frame time) — the main window writes the time
    # back to the tseep ``stage_1`` / ``stage_2`` control (saved to the sheet on Save).
    tag_requested = Signal(int, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.canvas = MplCanvas(self)
        self._frames = []        # plottable solution dicts (each carries 'time')
        self._times = []         # aligned frame times
        self._seep_data = None
        self._idx = 0
        self._opts_getter = lambda: {}
        self._style_getter = lambda: None

        # --- transport buttons (text glyphs so they render under any font) ---
        self._btn_first = _tool("|<", "First frame")
        self._btn_prev = _tool("<", "Previous frame")
        self._btn_play = _tool("Play", "Play / pause", checkable=True)
        self._btn_next = _tool(">", "Next frame")
        self._btn_last = _tool(">|", "Last frame")

        self._slider = QSlider(Qt.Horizontal)
        self._slider.setMinimum(0)
        self._slider.setMaximum(0)
        self._slider.setPageStep(1)
        self._slider.setToolTip("Scrub through the saved frames.")

        self._readout = QLineEdit()
        self._readout.setMaximumWidth(120)
        self._readout.setToolTip("Frame time — type a time and press Enter to jump "
                                 "to the nearest saved frame.")

        self._speed = QComboBox()
        for label, mult in _SPEEDS:
            self._speed.addItem(label, mult)
        self._speed.setCurrentIndex(1)   # 1x
        self._speed.setToolTip("Playback speed.")

        self._btn_tag1 = QPushButton("Set Stage 1")
        self._btn_tag1.setToolTip("Tag the current time as rapid-drawdown stage 1 "
                                  "(saved to the tseep sheet on Save).")
        self._btn_tag2 = QPushButton("Set Stage 2")
        self._btn_tag2.setToolTip("Tag the current time as rapid-drawdown stage 2 "
                                  "(saved to the tseep sheet on Save).")

        bar = QHBoxLayout()
        for w in (self._btn_first, self._btn_prev, self._btn_play,
                  self._btn_next, self._btn_last):
            bar.addWidget(w)
        bar.addWidget(self._slider, 1)
        bar.addWidget(QLabel("t ="))
        bar.addWidget(self._readout)
        bar.addWidget(QLabel("Speed"))
        bar.addWidget(self._speed)
        bar.addSpacing(12)
        bar.addWidget(self._btn_tag1)
        bar.addWidget(self._btn_tag2)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas, 1)
        self._bar = QWidget()
        self._bar.setLayout(bar)
        layout.addWidget(self._bar)

        self._timer = QTimer(self)
        self._timer.timeout.connect(self._advance_playing)

        self._btn_first.clicked.connect(lambda: self.set_index(0))
        self._btn_prev.clicked.connect(lambda: self.set_index(self._idx - 1))
        self._btn_next.clicked.connect(lambda: self.set_index(self._idx + 1))
        self._btn_last.clicked.connect(lambda: self.set_index(len(self._frames) - 1))
        self._btn_play.toggled.connect(self._on_play_toggled)
        self._slider.valueChanged.connect(self.set_index)
        self._readout.editingFinished.connect(self._on_readout)
        self._speed.currentIndexChanged.connect(self._on_speed)
        self._btn_tag1.clicked.connect(lambda: self._emit_tag(1))
        self._btn_tag2.clicked.connect(lambda: self._emit_tag(2))

    # --- data -----------------------------------------------------------
    def set_frames(self, seep_data, frames, opts_getter=None, style_getter=None,
                   keep_index=True):
        """Load a new run's frames. ``frames`` are the plottable solution dicts
        ``plot_seep_solution`` consumes (each carrying its ``time``). ``opts_getter``
        / ``style_getter`` are called at render time to pull the current display
        options and style, so a Display-panel change re-renders the shown frame with
        the new options (via :meth:`rerender`)."""
        self._seep_data = seep_data
        self._frames = list(frames)
        self._times = [float(f["time"]) if f.get("time") is not None else float(i)
                       for i, f in enumerate(self._frames)]
        if opts_getter is not None:
            self._opts_getter = opts_getter
        if style_getter is not None:
            self._style_getter = style_getter
        n = len(self._frames)
        self._slider.blockSignals(True)
        self._slider.setMaximum(max(n - 1, 0))
        self._slider.blockSignals(False)
        idx = self._idx if keep_index else n - 1
        self._idx = min(max(idx, 0), max(n - 1, 0))
        self._sync_controls()
        self.render_current()

    def render_current(self):
        """Render the currently-selected frame through the standard seep path."""
        if not self._frames or self._seep_data is None:
            return
        sol = self._frames[self._idx]
        try:
            self.canvas.render_seep_solution(
                self._seep_data, sol, self._opts_getter() or {},
                style=self._style_getter() or None)
        except Exception:
            traceback.print_exc()

    def rerender(self):
        """Re-render the current frame — called when the Display panel changes or a
        style edit needs to preview here."""
        self.render_current()

    # --- navigation -----------------------------------------------------
    def set_index(self, i):
        if not self._frames:
            return
        i = min(max(int(i), 0), len(self._frames) - 1)
        changed = i != self._idx
        self._idx = i
        self._sync_controls()
        if changed:
            self.render_current()

    @property
    def current_time(self):
        return self._times[self._idx] if self._times else None

    def current_frame(self):
        return self._frames[self._idx] if self._frames else None

    def frame_count(self):
        return len(self._frames)

    def _sync_controls(self):
        self._slider.blockSignals(True)
        self._slider.setValue(self._idx)
        self._slider.blockSignals(False)
        t = self.current_time
        self._readout.blockSignals(True)
        self._readout.setText("" if t is None else f"{t:g}")
        self._readout.blockSignals(False)
        n = len(self._frames)
        at_start = self._idx <= 0
        at_end = self._idx >= n - 1
        self._btn_first.setEnabled(not at_start)
        self._btn_prev.setEnabled(not at_start)
        self._btn_next.setEnabled(not at_end)
        self._btn_last.setEnabled(not at_end)

    def _on_readout(self):
        txt = self._readout.text().strip()
        if not txt or not self._times:
            self._sync_controls()
            return
        try:
            t = float(txt)
        except ValueError:
            self._sync_controls()   # restore the valid value
            return
        i = int(np.argmin(np.abs(np.asarray(self._times) - t)))
        self.set_index(i)

    # --- playback -------------------------------------------------------
    def _interval_ms(self):
        mult = self._speed.currentData() or 1.0
        return max(int(_BASE_INTERVAL_MS / mult), 30)

    def _on_speed(self, *_):
        if self._timer.isActive():
            self._timer.start(self._interval_ms())

    def _on_play_toggled(self, on):
        self._btn_play.setText("Pause" if on else "Play")
        if on:
            if not self._frames:
                self._btn_play.setChecked(False)
                return
            if self._idx >= len(self._frames) - 1:   # restart from the top
                self.set_index(0)
            self._timer.start(self._interval_ms())
        else:
            self._timer.stop()

    def _advance_playing(self):
        if self._idx >= len(self._frames) - 1:
            self._btn_play.setChecked(False)   # stop at the last frame
            return
        self.set_index(self._idx + 1)

    # --- tagging --------------------------------------------------------
    def _emit_tag(self, stage_no):
        t = self.current_time
        if t is not None:
            self.tag_requested.emit(stage_no, float(t))

    # --- tab plumbing ---------------------------------------------------
    def ensure_fitted(self):
        """Forwarded so the tab-shown fit path (main window) works on this wrapper."""
        self.canvas.ensure_fitted()


def _tool(text, tip, checkable=False):
    b = QToolButton()
    b.setText(text)
    b.setToolTip(tip)
    if checkable:
        b.setCheckable(True)
    return b
