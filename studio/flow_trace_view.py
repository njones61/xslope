"""The flow trace on Studio's steady seepage view: particles that follow the
solved velocity, animated over the plot, with a play bar below it.

The plot is a pixmap in a QGraphicsView (:class:`~studio.canvas.MplCanvas`); the
streaks are a separate item drawn over it in scene coordinates, so a frame
repaints only the streaks (the plot under them is not redrawn) and zoom and pan
apply to both as one picture. The trace itself (:func:`xslope.flow_trace.
build_flow_trace`) is built once per solution and particle count, off the GUI
thread; each frame is then an interpolation along its stored paths.

The animation runs only while it can be seen: it pauses when its tab is not the
current one, when the window is minimized or hidden, and while an animation is
being written, and resumes where it left off. Unticking Flow trace stops it and
removes the streaks.

Save animation… writes one loop of the trace over the same plot (same layers,
the canvas's figure size) through :func:`xslope.plot_seep.
save_flow_trace_animation`, with the transient play bar's dialog and progress.
"""

from __future__ import annotations

import os
import threading
import traceback

import numpy as np
from PySide6.QtCore import QElapsedTimer, QEvent, QLineF, QObject, QRectF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import QGraphicsItem, QHBoxLayout, QLabel, QWidget

from .canvas import BASE_DPI
from .transient import _tool, animation_save_path, save_animation_button, write_animation

#: Live frame interval (ms). The trace's own pace is set in display seconds, so
#: the interval only sets how smooth it looks, not how fast it moves.
FRAME_MS = 33
#: Streak opacity is drawn in this many steps (one pen per step and tail piece).
_ALPHA_STEPS = 16


class _StreakItem(QGraphicsItem):
    """The streaks of one frame, drawn over the plot's pixmap and clipped to the
    plot area. Lines are grouped by pen (tail piece x opacity step), so a frame
    is a few dozen ``drawLines`` calls."""

    def __init__(self):
        super().__init__()
        self._rect = QRectF()
        self._groups = []
        self.setZValue(10)

    def set_rect(self, rect):
        self.prepareGeometryChange()
        self._rect = QRectF(rect)

    def set_groups(self, groups):
        self._groups = groups
        self.update()

    def boundingRect(self):
        return self._rect

    def paint(self, painter, option, widget=None):
        if not self._groups:
            return
        painter.setClipRect(self._rect)
        painter.setRenderHint(QPainter.Antialiasing, True)
        for pen, lines in self._groups:
            painter.setPen(pen)
            painter.drawLines(lines)


class _Built(QObject):
    """Carries a finished build from the worker thread to the GUI thread."""
    done = Signal(object, object, object)     # generation, result, error text


class FlowTraceBar(QWidget):
    """The play bar under the plot while the trace is on: play/pause, a status
    line, and Save animation…, built from the transient play bar's parts."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.btn_play = _tool("Play", "Play / pause the flow trace.", checkable=True)
        self.status = QLabel("")
        self.status.setStyleSheet("color: gray;")
        self.btn_save = save_animation_button(
            "Write one loop of the flow trace to an animated GIF, drawn over the "
            "plot as shown (same layers and size).")
        bar = QHBoxLayout(self)
        bar.addWidget(self.btn_play)
        bar.addWidget(self.status, 1)
        bar.addWidget(self.btn_save)


class FlowTracePlayer(QObject):
    """Plays a flow trace over one steady seepage canvas.

    The main window calls :meth:`start` with the shown solution whenever the
    Flow trace box is ticked or the solution or particle count changes, and
    :meth:`stop` when the box is unticked or the view goes away.
    """

    def __init__(self, canvas, suggested_name=None):
        super().__init__(canvas)
        self.canvas = canvas
        self._suggested_name = suggested_name or (lambda: "flow_trace.gif")
        self.item = _StreakItem()
        self.item.setVisible(False)
        canvas.scene.addItem(self.item)
        self.bar = FlowTraceBar()
        self.bar.setVisible(False)
        canvas.add_bottom_bar(self.bar)
        self.bar.btn_play.toggled.connect(self._on_play_toggled)
        self.bar.btn_save.clicked.connect(self.save_animation)

        self._timer = QTimer(self)
        self._timer.setTimerType(Qt.PreciseTimer)
        self._timer.setInterval(FRAME_MS)
        self._timer.timeout.connect(self._tick)
        self._clock = QElapsedTimer()
        self._t = 0.0                 # display time when the clock last started
        self._enabled = False
        self._exporting = False
        self._trace = None
        self._key = None              # (solution id, n_particles) of _trace
        self._source = None           # (seep_data, solution) of _trace
        self._gen = 0                 # build generation; a stale build is dropped
        self._building = False
        self._pending_key = None
        self._opts_getter = lambda: {}
        self._style_getter = lambda: None
        self._affine = None
        self._pens = None
        self.frames_drawn = 0         # for checks and timing

        canvas.rasterized.connect(self._on_rasterized)
        canvas.installEventFilter(self)
        self._window = None
        # Show / hide / minimize are acted on just after they land, when
        # isVisible and isMinimized already say so. A child timer, so it dies
        # with the player.
        self._recheck = QTimer(self)
        self._recheck.setSingleShot(True)
        self._recheck.setInterval(0)
        self._recheck.timeout.connect(self._update_running)

    # --- control ----------------------------------------------------------
    def start(self, seep_data, solution, n_particles, opts_getter=None,
              style_getter=None):
        """Show the trace of ``solution``; build it first (off the GUI thread)
        unless the current trace is of the same solution and particle count."""
        from xslope.plot_seep import FLOW_TRACE_DEFAULTS
        if opts_getter is not None:
            self._opts_getter = opts_getter
        if style_getter is not None:
            self._style_getter = style_getter
        self._watch_window()
        if not self._enabled:
            # Ticking the box plays at once; later calls (a display option
            # changed, the plot redrawn) leave a pause the user made alone.
            self._enabled = True
            self._set_play(True)
        self.bar.setVisible(True)
        key = (id(solution), int(n_particles))
        if key == self._key and self._trace is not None:
            self._update_running()
            return
        if self._building and self._pending_key == key:
            return
        # A new solution or count: drop what is shown and build afresh.
        self._halt()
        self._trace = self._key = self._source = None
        self.item.set_groups([])
        self._t = 0.0
        self._gen += 1
        gen = self._gen
        self._building = True
        self._pending_key = key
        self.bar.status.setText("Tracing the flow…")
        self.bar.btn_play.setEnabled(False)
        self.bar.btn_save.setEnabled(False)
        d = FLOW_TRACE_DEFAULTS
        built = _Built()
        built.done.connect(self._on_built)
        self._built = built           # keep the emitter alive until it reports

        def work():
            from xslope.flow_trace import build_flow_trace
            try:
                tr = build_flow_trace(
                    seep_data, solution, n_particles=int(n_particles),
                    inflow_share=d["inflow_share"], speed_scale=d["speed_scale"],
                    loop_seconds=d["loop_seconds"])
                result, err = (tr, key, seep_data, solution), None
            except Exception as exc:
                result, err = None, f"{exc}\n{traceback.format_exc()}"
            try:
                built.done.emit(gen, result, err)
            except RuntimeError:      # the view went away while it was building
                pass

        threading.Thread(target=work, name="flow-trace-build", daemon=True).start()

    def stop(self):
        """Stop the animation and remove the streaks and the play bar. The built
        trace is kept, so ticking the box again resumes without a rebuild."""
        self._enabled = False
        self._halt()
        self._gen += 1                # a build still running is dropped
        self._building = False
        self.item.set_groups([])
        self.item.setVisible(False)
        self.bar.setVisible(False)

    def is_running(self):
        return self._timer.isActive()

    @property
    def trace(self):
        return self._trace

    def display_time(self):
        """The trace time shown now (seconds into the loop)."""
        t = self._t
        if self._timer.isActive():
            t += self._clock.elapsed() / 1000.0
        return t % self._trace.P if self._trace is not None else t

    # --- build ------------------------------------------------------------
    def _on_built(self, gen, result, err):
        if gen != self._gen:
            return
        self._building = False
        if result is None:
            print("Flow trace: could not trace the flow:\n" + (err or ""))
            self.bar.status.setText("Could not trace the flow (see the Log pane).")
            return
        trace, key, seep_data, solution = result
        self._trace, self._key, self._source = trace, key, (seep_data, solution)
        self._pens = None
        info = trace.info
        what = ("seepage velocity" if info.get("velocity_kind") == "seepage"
                else "Darcy velocity")
        self.bar.status.setText(
            f"Particles follow the {what}; slow zones are sped up. "
            f"{trace.P:g} s loop.")
        self.bar.btn_play.setEnabled(True)
        self.bar.btn_save.setEnabled(True)
        self._set_play(True)
        self._on_rasterized()
        self._update_running()

    # --- playback ---------------------------------------------------------
    def _set_play(self, on):
        """Set the play button without running its handler (the caller
        updates the timer)."""
        self.bar.btn_play.blockSignals(True)
        self.bar.btn_play.setChecked(on)
        self.bar.btn_play.blockSignals(False)
        self.bar.btn_play.setText("Pause" if on else "Play")

    def _on_play_toggled(self, on):
        self.bar.btn_play.setText("Pause" if on else "Play")
        self._update_running()

    def _can_show(self):
        win = self.canvas.window()
        return (self.canvas.isVisible() and win is not None
                and not win.isMinimized())

    def _should_run(self):
        return (self._enabled and self._trace is not None and self._affine is not None
                and self.bar.btn_play.isChecked() and not self._exporting
                and self._can_show())

    def _halt(self):
        if self._timer.isActive():
            self._t += self._clock.elapsed() / 1000.0
            self._timer.stop()

    def _update_running(self):
        if self._should_run():
            if not self._timer.isActive():
                self._clock.start()
                self._timer.start()
                self.item.setVisible(True)
                self._draw(self.display_time())
        else:
            self._halt()

    def _tick(self):
        if not self._should_run():
            self._update_running()
            return
        self._draw(self.display_time())

    def _draw(self, t):
        if self._trace is None or self._affine is None:
            return
        seg, alpha, rank = self._trace.frame(t)
        self.item.set_groups(self._groups(seg, alpha, rank))
        self.frames_drawn += 1

    def _pen_table(self):
        """One pen per (tail piece, opacity step): the head width in scene units
        (points at BASE_DPI) tapering toward the tail, as the saved animation."""
        if self._pens is None:
            from matplotlib.colors import to_rgba
            from xslope.plot_seep import (FLOW_TRACE_DEFAULTS, _FLOW_TRACE_COLOR,
                                          flow_trace_widths)
            r, g, b, _ = to_rgba(_FLOW_TRACE_COLOR)
            nt = self._trace.n_tail
            widths = flow_trace_widths(np.arange(nt), nt,
                                       FLOW_TRACE_DEFAULTS["linewidth"]) * BASE_DPI / 72.0
            pens = []
            for k in range(nt):
                row = []
                for j in range(_ALPHA_STEPS):
                    c = QColor.fromRgbF(r, g, b, (j + 0.5) / _ALPHA_STEPS)
                    pen = QPen(c)
                    pen.setWidthF(float(widths[k]))
                    # Square caps close the joints between a streak's pieces as
                    # round ones do, at a fraction of the drawing cost.
                    pen.setCapStyle(Qt.SquareCap)
                    row.append(pen)
                pens.append(row)
            self._pens = pens
        return self._pens

    def _groups(self, seg, alpha, rank):
        if len(seg) == 0:
            return []
        A, b, _ = self._affine
        sc = seg.reshape(-1, 2) @ A.T + b
        step = np.clip((alpha * _ALPHA_STEPS).astype(int), 0, _ALPHA_STEPS - 1)
        key = rank.astype(int) * _ALPHA_STEPS + step
        order = np.argsort(key, kind="stable")
        ks = key[order]
        flat = sc.reshape(-1, 4)[order].tolist()
        cut = np.flatnonzero(np.diff(ks)) + 1
        starts = np.concatenate([[0], cut])
        ends = np.concatenate([cut, [len(ks)]])
        pens = self._pen_table()
        out = []
        for a, e in zip(starts.tolist(), ends.tolist()):
            k = int(ks[a])
            out.append((pens[k // _ALPHA_STEPS][k % _ALPHA_STEPS],
                        [QLineF(*v) for v in flat[a:e]]))
        return out

    # --- following the plot and the window ---------------------------------
    def _on_rasterized(self):
        """The plot was redrawn (new options, a resize, a zoom refine): place the
        streaks on its axes again."""
        m = self.canvas.scene_transform()
        self._affine = m
        if m is None:
            self.item.setVisible(False)
            self._halt()
            return
        self.item.set_rect(m[2])
        if self._enabled and self._trace is not None:
            self.item.setVisible(True)
            if not self._timer.isActive():
                self._draw(self.display_time())
        self._update_running()

    def _watch_window(self):
        win = self.canvas.window()
        if win is not None and win is not self._window:
            if self._window is not None:
                self._window.removeEventFilter(self)
            win.installEventFilter(self)
            self._window = win

    def eventFilter(self, obj, event):
        if event.type() in (QEvent.Show, QEvent.Hide, QEvent.WindowStateChange):
            self._recheck.start()
        return False

    # --- export -----------------------------------------------------------
    def save_animation(self, _checked=False, path=None):
        """Write one loop of the trace over the plot as shown (same layers, the
        canvas's figure size) to a GIF (MP4 when ffmpeg is installed). ``path``
        skips the file dialog (scripted use). Returns the export summary, or
        None."""
        if self._trace is None or self._source is None:
            return None
        from xslope.plot_seep import FLOW_TRACE_DEFAULTS, save_flow_trace_animation
        path = animation_save_path(self.canvas, self._suggested_name(), path)
        if path is None:
            return None
        seep_data, solution = self._source
        trace = self._trace
        fps = FLOW_TRACE_DEFAULTS["fps"]
        n = max(int(round(trace.P * fps)), 1)
        w, h = self.canvas.figure.get_size_inches()
        self._exporting = True
        self._update_running()
        try:
            return write_animation(
                self.canvas, path, n,
                lambda progress: save_flow_trace_animation(
                    seep_data, solution, trace, path,
                    options=self._opts_getter() or {},
                    style=self._style_getter() or None,
                    figsize=(float(w), float(h)), dpi=100, fps=fps,
                    progress=progress))
        finally:
            self._exporting = False
            self._update_running()


def suggested_name(doc, bc=1):
    """``{stem}_flow_trace.gif`` next to the model (``_flow_trace2`` for BC 2)."""
    tail = "_flow_trace.gif" if bc == 1 else f"_flow_trace{bc}.gif"
    path = getattr(doc, "path", None)
    if path:
        return os.path.splitext(path)[0] + tail
    return "seepage" + tail
