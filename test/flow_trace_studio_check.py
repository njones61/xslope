"""Studio's Flow trace box on the steady seepage view plays, pauses, stops and
saves as it should.

A copy of the earth dam sample (docs/seep/files/xslope_earth_dam2.xlsx, opened
with its saved steady solution, so nothing is solved) in a real Studio window,
offscreen, with Studio's settings sandboxed (tools.studio_settings).

  A. TICK. Ticking Flow trace on the Seep · Solution view builds the trace and
     starts the animation within 2 s: the timer runs, frames are drawn, and the
     streaks drawn at two moments differ. The play bar appears.
  B. HIDDEN. With another tab shown the timer stops and no frame is drawn; with
     the tab back it runs again. The same for a minimized window.
  C. REDRAWN. Turning Filled contours on (the plot is redrawn under the streaks)
     and zooming in leave it running.
  D. UNTICK. Unticking stops the timer and removes the streaks and the play bar;
     ticking again plays at once on the same trace (no rebuild).
  E. SAVE. Save animation writes a GIF of one loop (8 s at 20 fps: 160 frames of
     50 ms) at the canvas's figure size, the plot drawn with the panel's options:
     a frame matches plot_seep_solution drawn directly with them on most of its
     pixels (the streaks are the rest), and the streaks are there (blue pixels
     where the plain plot has another color).
  F. CLEARED. Opening another file stops the animation of the view it removes.
  G. TRANSIENT. The transient view's Display panel has no Flow trace box (a copy
     of xslope_earth_dam_tseep), and asks for no trace even if the hidden box is
     ticked.

  H. NO FLOW. Started on the same model with a solution whose head is one value
     everywhere (velocity round-off, 1e-13), the play bar says nothing flows,
     Play and Save stay off, and nothing is printed to the Log (no traceback).
     Before the fix it animated the round-off.
  I. PLAYBACK COST. At the Trace particles maximum, one frame (the trace's
     frame, the pen groups, and a repaint of the view) takes less than the
     animation's timer interval (FRAME_MS, 33 ms), fastest of 30 frames, so a
     machine busy with other work (the overnight suite runs 8 at once) does not
     fail it. One QLineF per segment took about 57 ms at 5000 particles.
  J. TOOLTIP. The Flow trace box's tooltip names the Darcy velocity (Studio
     passes no porosity), not the seepage velocity.

Skipped without PySide6.
"""
import contextlib
import io
import os
import shutil
import sys
import tempfile
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import matplotlib                                        # noqa: E402
matplotlib.use("Agg")

_FILES = os.path.join(_REPO, "docs", "seep", "files")
_STEADY = "xslope_earth_dam2"
_TSEEP = "xslope_earth_dam_tseep"


def _copy(stem, tmp):
    for f in os.listdir(_FILES):
        if f.startswith(stem + ".") or f.startswith(stem + "_"):
            shutil.copy(os.path.join(_FILES, f), tmp)
    return os.path.join(tmp, stem + ".xlsx")


def _spin(app, seconds, until=None):
    t = time.perf_counter()
    while time.perf_counter() - t < seconds:
        app.processEvents()
        if until is not None and until():
            return True
        time.sleep(0.002)
    return until() if until is not None else True


def _lines(player):
    """The streak coordinates the item would paint now."""
    pts = []
    for _, path in player.item._groups:
        for i in range(0, path.elementCount(), 2):
            a, b = path.elementAt(i), path.elementAt(i + 1)
            pts.append((a.x, a.y, b.x, b.y))
    return np.array(sorted(pts))


def run():
    try:
        import PySide6                                   # noqa: F401
    except Exception:
        print("  skipped: PySide6 is not installed")
        return []
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from tools import studio_settings
    return studio_settings.isolated(_checks)()


def _checks():
    from PySide6.QtWidgets import QApplication
    from studio.main_window import MainWindow
    import xslope
    app = QApplication.instance()
    failures = []
    if not os.path.abspath(xslope.__file__).startswith(_REPO):
        return [f"xslope imported from {xslope.__file__}, not {_REPO}"]

    with tempfile.TemporaryDirectory() as tmp:
        mw = MainWindow()
        mw.resize(1400, 900)
        mw.show()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                mw.open_path(_copy(_STEADY, tmp))
            canvas = mw.seep_solution_canvas.get(1)
            if canvas is None:
                return ["the steady seepage view did not open with the saved solution"]
            mw.view_tabs.setCurrentWidget(canvas)
            _spin(app, 0.3)
            panel = mw._display_panels[canvas]
            if not hasattr(panel, "flow_trace") or not panel.flow_trace.isVisibleTo(panel):
                return ["A. the steady Display panel has no Flow trace box"]

            # A. tick
            t0 = time.perf_counter()
            panel.flow_trace.setChecked(True)
            player = getattr(canvas, "flow_trace", None)
            ok = player is not None and _spin(
                app, 2.0, lambda: player.frames_drawn > 0 and player.is_running())
            t_first = time.perf_counter() - t0
            if not ok:
                return ["A. ticking Flow trace did not start the animation within 2 s"]
            a = _lines(player)
            _spin(app, 0.25)
            b = _lines(player)
            moved = a.shape != b.shape or not np.allclose(a, b)
            print(f"  A. first frame {t_first:.2f} s after the tick; {len(a)} streak pieces; "
                  f"streaks moved {moved}; play bar shown {player.bar.isVisibleTo(canvas)}")
            if not moved:
                failures.append("A. the streaks did not change between frames")
            if not player.bar.isVisibleTo(canvas):
                failures.append("A. the play bar is not shown")

            # B. hidden
            mw.view_tabs.setCurrentIndex(0)
            _spin(app, 0.1)
            n = player.frames_drawn
            _spin(app, 0.3)
            hidden_ok = not player.is_running() and player.frames_drawn == n
            mw.view_tabs.setCurrentWidget(canvas)
            back = _spin(app, 0.5, player.is_running)
            mw.showMinimized()
            _spin(app, 0.1)
            n = player.frames_drawn
            _spin(app, 0.3)
            min_ok = not player.is_running() and player.frames_drawn == n
            mw.showNormal()
            back_min = _spin(app, 0.5, player.is_running)
            print(f"  B. other tab: paused {hidden_ok}, resumed {back}; "
                  f"minimized: paused {min_ok}, resumed {back_min}")
            if not (hidden_ok and back and min_ok and back_min):
                failures.append(f"B. hidden {hidden_ok}/{back}, minimized {min_ok}/{back_min}")

            # C. redrawn plot, zoom
            panel.fill.setChecked(True)
            _spin(app, 0.5)
            canvas.zoom_in()
            n = player.frames_drawn
            running = _spin(app, 0.5, lambda: player.frames_drawn > n + 3)
            print(f"  C. after Filled contours and a zoom: running {running}")
            if not running:
                failures.append("C. the animation stopped when the plot was redrawn or zoomed")
            canvas.fit()

            # D. untick, tick again
            trace = player.trace
            panel.flow_trace.setChecked(False)
            _spin(app, 0.1)
            gone = (not player.is_running() and not player.item.isVisible()
                    and not player.item._groups and not player.bar.isVisibleTo(canvas))
            panel.flow_trace.setChecked(True)
            again = _spin(app, 0.5, player.is_running)
            same = player.trace is trace
            print(f"  D. untick: stopped and removed {gone}; tick again: playing {again}, "
                  f"same trace {same}")
            if not (gone and again and same):
                failures.append(f"D. untick removed {gone}; re-tick played {again}, same trace {same}")

            # E. save
            failures += _save_check(app, canvas, panel, player, tmp)

            # J. tooltip
            tip = panel.flow_trace.toolTip()
            print(f"  J. tooltip names the Darcy velocity {'Darcy velocity' in tip}, "
                  f"the seepage velocity {'seepage velocity' in tip}")
            if "Darcy velocity" not in tip or "seepage velocity" in tip:
                failures.append(f"J. the Flow trace tooltip reads: {tip}")

            # I. playback cost at the maximum particle count
            failures += _cost_check(app, panel, player)

            # H. no flow
            failures += _no_flow_check(app, panel, player)

            # F. another file clears the view and stops it
            with contextlib.redirect_stdout(io.StringIO()):
                mw.open_path(_copy(_TSEEP, tmp))
            _spin(app, 0.3)
            stopped = not player.is_running() if _alive(player) else True
            print(f"  F. another file opened: the old view's animation stopped {stopped}")
            if not stopped:
                failures.append("F. the animation kept running after its view was cleared")

            # G. transient
            view = mw.transient_seep_view
            tpanel = mw._display_panels.get(view) if view is not None else None
            if tpanel is None:
                failures.append("G. the transient view did not open")
            else:
                shown = tpanel.flow_trace.isVisibleTo(tpanel)
                tpanel.flow_trace.setChecked(True)
                asks = tpanel.trace_options()["on"]
                print(f"  G. transient panel: Flow trace shown {shown}, asks for a trace {asks}")
                if shown or asks:
                    failures.append(f"G. transient panel shows the box {shown}, asks {asks}")
        finally:
            mw.doc._dirty = False
            mw.close()
    return failures


def _cost_check(app, panel, player):
    from studio.flow_trace_view import FRAME_MS
    n_max = panel.trace_particles.maximum()
    panel.trace_particles.setValue(n_max)
    built = _spin(app, 120, lambda: player.trace is not None and player._key[1] == n_max)
    if not built:
        return [f"I. no trace built at {n_max} particles within 120 s"]
    player.bar.btn_play.setChecked(False)
    _spin(app, 0.1)
    vp = player.canvas.view.viewport()
    tr = player.trace
    cost, segs = [], 0
    for t in np.linspace(0, tr.P, 30, endpoint=False):
        t0 = time.perf_counter()
        seg, a, r = tr.frame(t)
        player.item.set_groups(player._groups(seg, a, r))
        vp.repaint()
        cost.append(time.perf_counter() - t0)
        segs = max(segs, len(seg))
    med = 1e3 * float(np.min(cost))
    print(f"  I. {n_max} particles (the maximum), up to {segs} streak pieces: a frame takes "
          f"{med:.1f} ms (fastest of 30; timer {FRAME_MS} ms)")
    player.bar.btn_play.setChecked(True)
    panel.trace_particles.setValue(1500)
    _spin(app, 60, lambda: player.trace is not None and player._key[1] == 1500)
    if not med < FRAME_MS:
        return [f"I. a frame at {n_max} particles takes {med:.1f} ms, over the {FRAME_MS} ms timer"]
    return []


def _no_flow_check(app, panel, player):
    seep_data, solution = player._source
    flat = dict(solution)
    flat["head"] = np.full(len(solution["head"]), 10.0)
    rng = np.random.default_rng(0)
    flat["velocity"] = rng.normal(0, 1e-13, np.shape(solution["velocity"]))
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        player.start(seep_data, flat, 1500)
        _spin(app, 30, lambda: not player._building)
    text = player.bar.status.text()
    log = out.getvalue()
    play_off = not player.bar.btn_play.isEnabled() and not player.bar.btn_save.isEnabled()
    print(f"  H. no flow: status '{text}'; Play and Save off {play_off}; "
          f"Log lines {len(log.splitlines())}")
    failures = []
    if "Nothing flows" not in text or not play_off or player.trace is not None:
        failures.append(f"H. no flow: status '{text}', buttons off {play_off}")
    if "Traceback" in log or "could not trace" in log:
        failures.append("H. no flow printed a traceback to the Log")
    player.start(seep_data, solution, 1500)
    _spin(app, 30, lambda: player.trace is not None)
    return failures


def _alive(obj):
    try:
        from shiboken6 import isValid
        return isValid(obj)
    except Exception:
        return True


def _save_check(app, canvas, panel, player, tmp):
    from PIL import Image
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from xslope.plot_seep import plot_seep_solution, seep_solution_plot_kwargs
    failures = []
    gif = os.path.join(tmp, "trace.gif")
    t = time.perf_counter()
    info = player.save_animation(path=gif)
    dt = time.perf_counter() - t
    running_after = _spin(app, 0.5, player.is_running)
    if info is None or not os.path.exists(gif):
        return ["E. Save animation wrote nothing"]
    im = Image.open(gif)
    nf, dur = im.n_frames, im.info.get("duration")
    w, h = canvas.figure.get_size_inches()
    want_size = (int(round(w * 100)), int(round(h * 100)))
    # A frame against the plot drawn directly with the panel's options.
    seep_data, solution = player._source
    fig = Figure(figsize=(w, h), dpi=100)
    FigureCanvasAgg(fig)
    plot_seep_solution(seep_data, solution, fig=fig,
                       **seep_solution_plot_kwargs(panel.options()))
    fig.canvas.draw()
    plain = np.asarray(fig.canvas.buffer_rgba())[..., :3].astype(int)
    im.seek(nf // 2)
    got = np.asarray(im.convert("RGB")).astype(int)
    same_shape = got.shape == plain.shape
    close = np.abs(got - plain).max(axis=2) <= 24 if same_shape else np.zeros(1, bool)
    match = float(close.mean())
    # The streaks: blue pixels the plain plot does not have there (the GIF's
    # shared palette moves the streak blue, so its hue is what is tested).
    blue = (got[..., 2] > got[..., 0] + 60) & (got[..., 2] > got[..., 1] + 40)
    extra = int((blue & ~close).sum()) if same_shape else 0
    print(f"  E. saved {nf} frames of {dur} ms in {dt:.1f} s, {im.size[0]} x {im.size[1]} "
          f"(want {want_size[0]} x {want_size[1]}); frame matches the plot at {match:.1%}, "
          f"{extra} streak pixels; playing again after {running_after}")
    if nf != 160 or dur != 50:
        failures.append(f"E. GIF has {nf} frames of {dur} ms (want 160 of 50)")
    if im.size != want_size:
        failures.append(f"E. GIF is {im.size}, want {want_size}")
    if not (match >= 0.9 and extra > 200):
        failures.append(f"E. frame matches the plot at {match:.1%} with {extra} streak pixels")
    if not running_after:
        failures.append("E. the animation did not resume after saving")
    return failures


if __name__ == "__main__":
    f = run()
    print("\n".join(["FAILED:"] + f) if f else "\nAll Studio flow trace checks passed.")
    sys.exit(1 if f else 0)
