"""The flow trace moves particles at the solved speed, seeds them by equal flux,
keeps them in the mesh, and loops without a seam; the animation exporter writes
what it is asked to.

The model is a 10 x 2 rectangle of linear triangles with a head of 10 on the left
face and 0 on the right, so the flow is uniform and horizontal: Darcy velocity
k * 10 / 10 = k. The rows are graded (finer toward mid-height), so equal spacing
of the release points is a property of the flux, not of the mesh. The halves
above and below the row at mid-height can carry different k.

  A. SPEED AND DIRECTION. With k = 1 and speed_scale "true", a particle's
     displacement over a stretch of display time, divided by the model time that
     stretch stands for, is the Darcy velocity (1, 0) to 1e-6; with porosity 0.25 it
     is the seepage velocity (4, 0).
  B. TRUE AND COMPRESSED PACE. With k = 1 below mid-height and 3 above, the top
     layer's particles move 3 times as fast as the bottom's on "true" and sqrt(3)
     times on "sqrt" (to 2%), along the same paths.
  C. EQUAL FLUX SEEDING. The inflow release points carry the solve's whole inflow
     (the reported flowrate, to 1e-9), sit on the inflow face, and are equally
     spaced in y on the uniform model despite the graded rows (to 1e-6); on the
     layered model three quarters of them are in the faster top half (to 1 point:
     the nodal flows spread the jump in k over the two edges beside the
     interface node, about half a point here). Ordering the inflow edges by their
     start point instead of walking them head to tail broke the spacing (0.026 to
     0.072 for 0.050), and so did sharing each node's flow equally between its two
     edges on graded rows.
  D. NOBODY LEAVES THE MESH ALIVE. Every inflow path ends on the outflow face, and
     at 40 frames over the loop every drawn point is inside the mesh or on its
     boundary (to 1e-9 of the model width).
  E. THE LOOP CLOSES. frame(P) is frame(0) exactly, and frame(t + P) is frame(t) at
     random t, for both seed kinds.
  F. THE EXPORTER. save_animation writes a GIF with the animation's frame count,
     endless looping and 50 ms frames at 20 fps. Asked for an MP4 without ffmpeg it
     raises a RuntimeError naming ffmpeg (with ffmpeg it writes a non-empty file).
  G. TRANSIENT EXPORT. earth_dam_tseep's shipped saved frames (read back, not
     re-solved), written by save_transient_animation at 0.5 s a frame: one GIF
     frame per saved frame, 500 ms each, titled with the frame times in order, and
     a written frame matches plot_seep_solution drawn directly with the transient
     view's options (97% of pixels within 8 levels; the rest are anti-aliased
     edges on the GIF's shared palette).
  H. STUDIO'S SAVE ANIMATION. A copy of the same workbook opened in a real Studio
     window (sandboxed settings, offscreen): the transient view's Save animation
     writes the view's frames in the view's order, at the play bar's speed (2x is
     250 ms a frame), drawn with the display panel's current options: with Filled
     contours turned on, a written frame matches a direct render with them on and
     not one with them off. Skipped without PySide6.
"""
import contextlib
import io
import os
import shutil
import sys
import tempfile

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import matplotlib
matplotlib.use("Agg")


def _rect(k_bottom=1.0, k_top=1.0, nx=20, L=10.0, H=2.0):
    """A confined rectangle, head 10 on the left face and 0 on the right."""
    # Graded rows, finer toward a row at mid-height (the material interface).
    ys = H * np.array([0, .15, .275, .375, .45, .5, .55, .625, .725, .85, 1.0])
    xs = np.linspace(0, L, nx + 1)
    X, Y = np.meshgrid(xs, ys)
    nodes = np.column_stack([X.ravel(), Y.ravel()])

    def nid(i, j):
        return j * (nx + 1) + i
    els, mats = [], []
    for j in range(len(ys) - 1):
        for i in range(nx):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            m = 2 if 0.5 * (ys[j] + ys[j + 1]) > H / 2 else 1
            els += [[a, b, c], [a, c, d]]
            mats += [m, m]
    n = len(nodes)
    bc_type = np.zeros(n, int)
    bc_values = np.zeros(n)
    left = np.isclose(nodes[:, 0], 0.0)
    right = np.isclose(nodes[:, 0], L)
    bc_type[left | right] = 1
    bc_values[left] = 10.0
    return {"nodes": nodes, "elements": np.array(els),
            "element_types": np.full(len(els), 3), "element_materials": np.array(mats),
            "bc_type": bc_type, "bc_values": bc_values, "flux_nodal": np.zeros(n),
            "k1_by_mat": np.array([k_bottom, k_top]), "k2_by_mat": np.array([k_bottom, k_top]),
            "angle_by_mat": np.zeros(2), "kr0_by_mat": np.full(2, 0.001),
            "h0_by_mat": np.full(2, -1.0), "unsat_by_mat": None, "vg_a_by_mat": None,
            "vg_n_by_mat": None, "vg_l_by_mat": None, "unit_weight": 9.81}


def _solve(seep):
    from xslope.seep import run_seepage_analysis
    with contextlib.redirect_stdout(io.StringIO()):
        return run_seepage_analysis(seep)


def _velocity_of(tr, k, dtau):
    """Mean model-time velocity of particle k over dtau of display time from the
    middle of its path (true scale: model time = display time * time_per_s)."""
    pid = tr.path_id[k]
    mid = 0.5 * tr.path_end[pid]
    p = tr._position(np.array([pid, pid]), np.array([mid, mid + dtau]))
    return (p[1] - p[0]) / (dtau * tr.info["time_per_display_second"])


def run():
    from xslope.flow_trace import build_flow_trace, inflow_release_points, TraceField, _nearest_edge
    failures = []

    seep = _rect()
    sol = _solve(seep)
    W = 10.0

    # A. speed and direction
    for por, want in ((None, 1.0), (0.25, 4.0)):
        tr = build_flow_trace(seep, sol, n_particles=200, inflow_share=1.0,
                              speed_scale="true", porosity=por, seed=1)
        ks = np.nonzero(tr.path_end[tr.path_id] > 0.5)[0][:50]
        v = np.array([_velocity_of(tr, k, 0.2) for k in ks])
        err = np.abs(v - [want, 0.0]).max()
        kind = tr.info["velocity_kind"]
        print(f"  A. porosity {por}: {kind} velocity, worst error {err:.2e} (want ({want:g}, 0))")
        if not err < 1e-6:
            failures.append(f"A. porosity {por}: particle velocity off by {err:.3g}")
        if kind != ("darcy" if por is None else "seepage"):
            failures.append(f"A. porosity {por}: velocity kind {kind}")

    # B. layered speeds, true and sqrt
    seep2 = _rect(1.0, 3.0)
    sol2 = _solve(seep2)
    for scale, want in (("true", 3.0), ("sqrt", np.sqrt(3.0))):
        tr = build_flow_trace(seep2, sol2, n_particles=400, inflow_share=1.0,
                              speed_scale=scale, seed=2)
        f = tr.field
        top, bot = [], []
        for k in range(len(tr.path_id)):
            pid = tr.path_id[k]
            y = tr._position(np.array([pid]), np.array([0.0]))[0, 1]
            if tr.path_end[pid] < 0.2:
                continue
            p = tr._position(np.array([pid, pid]), np.array([0.05, 0.15]))
            step = np.hypot(*(p[1] - p[0]))
            # Stay clear of the interface, where the nodal velocity is averaged.
            if y > 1.25:
                top.append(step)
            elif y < 0.75:
                bot.append(step)
        ratio = np.median(top) / np.median(bot)
        print(f"  B. {scale}: top/bottom display speed {ratio:.4f} (want {want:.4f})")
        if not abs(ratio / want - 1) < 0.02:
            failures.append(f"B. {scale}: speed ratio {ratio:.4f}, want {want:.4f}")

    # C. equal-flux seeding
    field = TraceField(seep, sol)
    n = 40
    pts, q_each, q_tot = inflow_release_points(field, n)
    ys = np.sort(pts[:, 1])
    spacing = np.diff(ys)
    print(f"  C. uniform: inflow {q_tot:.12g} vs flowrate {sol['flowrate']:.12g}; "
          f"x within {np.abs(pts[:, 0]).max():.1e} of the face; "
          f"spacing {spacing.min():.6f}..{spacing.max():.6f} (want {2.0 / n:.6f})")
    if not abs(q_tot - sol["flowrate"]) < 1e-9 * sol["flowrate"]:
        failures.append(f"C. inflow found {q_tot} != flowrate {sol['flowrate']}")
    if not np.abs(pts[:, 0]).max() < 1e-3:
        failures.append("C. release points are not on the inflow face")
    if not np.abs(spacing - 2.0 / n).max() < 1e-6:
        failures.append(f"C. release spacing {spacing.min()}..{spacing.max()}, want {2.0 / n}")
    field2 = TraceField(seep2, sol2)
    pts2, _, _ = inflow_release_points(field2, n)
    n_top = int((pts2[:, 1] > 1.0).sum())
    print(f"  C. layered k 1 / 3: {n_top} of {n} release points in the top half (want {0.75 * n:g})")
    if not abs(n_top - 0.75 * n) <= 1:
        failures.append(f"C. layered: {n_top} of {n} release points in the top half")

    # D and E on a mixed trace
    for scale in ("true", "sqrt"):
        tr = build_flow_trace(seep2, sol2, n_particles=300, inflow_share=0.5,
                              speed_scale=scale, seed=3)
        f = tr.field
        inflow = np.nonzero(tr.kind == "inflow")[0]
        ends = np.array([tr._position(np.array([p]), np.array([tr.path_end[p]]))[0]
                         for p in np.unique(tr.path_id[inflow])])
        off_face = np.abs(ends[:, 0] - 10.0).max()
        worst = 0.0
        for t in np.linspace(0, tr.P, 40, endpoint=False):
            seg, _, _ = tr.frame(t)
            p = seg.reshape(-1, 2)
            out = f.locate(p) < 0
            if out.any():
                q = p[out]
                e = _nearest_edge(f, q)
                A, B = f.bnd_pa[e], f.bnd_pb[e]
                d = B - A
                s = np.clip(((q - A) * d).sum(1) / (d * d).sum(1), 0, 1)
                worst = max(worst, float(np.hypot(*(A + d * s[:, None] - q).T).max()))
        print(f"  D. {scale}: inflow paths end within {off_face:.1e} of the outflow face; "
              f"farthest drawn point off the mesh {worst:.1e}")
        if not off_face < 1e-9 * W:
            failures.append(f"D. {scale}: an inflow path ends {off_face:.3g} from the outflow face")
        if not worst <= 1e-9 * W:
            failures.append(f"D. {scale}: a drawn point lies {worst:.3g} outside the mesh")

        s0, a0, r0 = tr.frame(0.0)
        sP, aP, rP = tr.frame(tr.P)
        same = s0.shape == sP.shape and np.array_equal(s0, sP) and np.array_equal(a0, aP)
        rng = np.random.default_rng(4)
        worst_t = 0.0
        for t in rng.random(5) * tr.P:
            s1, a1, _ = tr.frame(t)
            s2, a2, _ = tr.frame(t + tr.P)
            if s1.shape != s2.shape:
                worst_t = np.inf
                break
            worst_t = max(worst_t, float(np.abs(s1 - s2).max(initial=0)),
                          float(np.abs(a1 - a2).max(initial=0)))
        kinds = sorted(set(tr.kind[np.mod(-tr.phase, tr.P) < tr.L]))
        print(f"  E. {scale}: frame(P) == frame(0) {same}; worst frame(t+P) - frame(t) "
              f"{worst_t:.1e}; lit at t=0: {kinds}")
        if not same:
            failures.append(f"E. {scale}: frame(P) differs from frame(0)")
        if not worst_t < 1e-9:
            failures.append(f"E. {scale}: frame(t+P) differs from frame(t) by {worst_t:.3g}")

    # F. the exporter
    from PIL import Image
    from xslope.anim_export import save_animation, ffmpeg_available
    from xslope.plot_seep import animate_flow_trace
    import matplotlib.pyplot as plt
    with tempfile.TemporaryDirectory() as tmp:
        anim = animate_flow_trace(seep, sol, n_particles=100, loop_seconds=1.0, fps=20,
                                  dpi=40, figsize=(4, 2))
        gif = os.path.join(tmp, "t.gif")
        save_animation(anim, gif, fps=20)
        plt.close("all")
        im = Image.open(gif)
        nf, dur, loop = im.n_frames, im.info.get("duration"), im.info.get("loop")
        print(f"  F. GIF: {nf} frames (want 20), {dur} ms a frame (want 50), loop {loop} (want 0)")
        if (nf, dur, loop) != (20, 50, 0):
            failures.append(f"F. GIF has {nf} frames, {dur} ms, loop {loop}")
        mp4 = os.path.join(tmp, "t.mp4")
        has = ffmpeg_available()
        anim = animate_flow_trace(seep, sol, n_particles=50, loop_seconds=0.5, fps=10,
                                  dpi=40, figsize=(4, 2))
        try:
            save_animation(anim, mp4, fps=10)
            ok = has and os.path.getsize(mp4) > 0
            msg = "written" if ok else "written without ffmpeg?"
        except RuntimeError as e:
            ok = (not has) and "ffmpeg" in str(e)
            msg = f"RuntimeError: {e}"
        plt.close("all")
        print(f"  F. MP4 (ffmpeg {'found' if has else 'not found'}): {msg}")
        if not ok:
            failures.append(f"F. MP4 with ffmpeg {'found' if has else 'missing'}: {msg}")
    failures += _transient_checks()
    return failures


_TSEEP = os.path.join(_REPO, "docs/seep/files/xslope_earth_dam_tseep")


def _render(seep, frame, opts, figsize, dpi):
    """``frame`` drawn directly by plot_seep_solution, as an RGB array."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from xslope.plot_seep import plot_seep_solution, seep_solution_plot_kwargs
    fig = Figure(figsize=figsize, dpi=dpi)
    FigureCanvasAgg(fig)
    plot_seep_solution(seep, frame, fig=fig, **seep_solution_plot_kwargs(opts))
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    from PIL import Image
    return np.asarray(Image.open(buf).convert("RGB"), dtype=int)


def _gif_frame(path, k):
    from PIL import Image
    im = Image.open(path)
    im.seek(k)
    return np.asarray(im.convert("RGB"), dtype=int)


def _close_fraction(a, b):
    if a.shape != b.shape:
        return 0.0
    return float((np.abs(a - b).max(axis=2) <= 8).mean())


def _transient_checks():
    from PIL import Image
    from xslope.fileio import load_slope_data
    from xslope.seep import build_seep_data, import_transient_solution
    from xslope.plot_seep import save_transient_animation, TRANSIENT_VIEW_OPTIONS
    failures = []
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(_TSEEP + ".xlsx")
        seep = build_seep_data(sd["mesh"], sd, seep_bc=1)
        tr = import_transient_solution(seep, _TSEEP)
    frames = tr["frames"]
    times = [fr["time"] for fr in frames]

    # G. the Python export
    with tempfile.TemporaryDirectory() as tmp:
        gif = os.path.join(tmp, "t.gif")
        info = save_transient_animation(seep, tr, gif, seconds_per_frame=0.5,
                                        figsize=(10, 6), dpi=60)
        im = Image.open(gif)
        nf, dur = im.n_frames, im.info.get("duration")
        want_titles = [f"t = {t:g} " for t in times]
        in_order = all(w in t for w, t in zip(want_titles, info["titles"])) \
            and len(info["titles"]) == len(times)
        k = len(frames) // 2
        close = _close_fraction(_gif_frame(gif, k),
                                _render(seep, frames[k], TRANSIENT_VIEW_OPTIONS, (10, 6), 60))
        print(f"  G. transient GIF: {nf} frames (want {len(frames)}), {dur} ms (want 500), "
              f"titles in time order {in_order}; frame {k} matches a direct render "
              f"at {close:.1%}")
        if nf != len(frames) or dur != 500:
            failures.append(f"G. transient GIF has {nf} frames at {dur} ms")
        if not in_order:
            failures.append(f"G. frame titles {info['titles']} do not follow times {times}")
        if not close >= 0.97:
            failures.append(f"G. a written frame matches the direct render at only {close:.1%}")

    # H. Studio
    try:
        import PySide6  # noqa: F401
    except Exception:
        print("  H. PySide6 not installed: Studio check skipped")
        return failures
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from tools import studio_settings
    failures += studio_settings.isolated(_studio_check)(seep, frames)
    return failures


def _studio_check(seep, frames):
    from PySide6.QtWidgets import QApplication
    from PIL import Image
    from studio.main_window import MainWindow
    failures = []
    src = os.path.dirname(_TSEEP)
    with tempfile.TemporaryDirectory() as tmp:
        for f in os.listdir(src):
            if f.startswith(os.path.basename(_TSEEP)):
                shutil.copy(os.path.join(src, f), tmp)
        mw = MainWindow()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                mw.doc.load(os.path.join(tmp, os.path.basename(_TSEEP) + ".xlsx"))
            for _ in range(10):
                QApplication.processEvents()
            view = mw.transient_seep_view
            if view is None or view.frame_count() != len(frames):
                return [f"H. the transient view did not load the saved frames"]
            panel = mw._display_panels[view]
            panel.fill.setChecked(True)
            view._speed.setCurrentIndex(2)             # 2x
            gif = os.path.join(tmp, "studio.gif")
            info = view.save_animation(path=gif)
            im = Image.open(gif)
            nf, dur = im.n_frames, im.info.get("duration")
            same_order = info is not None and info["times"] == view._times
            w, h = view.canvas.figure.get_size_inches()
            k = len(frames) // 2
            got = _gif_frame(gif, k)
            on = _close_fraction(got, _render(seep, view._frames[k], panel.options(), (w, h), 100))
            off_opts = dict(panel.options(), fill_contours=False)
            off = _close_fraction(got, _render(seep, view._frames[k], off_opts, (w, h), 100))
            print(f"  H. Studio: {nf} frames (want {view.frame_count()}), {dur} ms at 2x "
                  f"(want 250), view's order {same_order}; frame {k} matches the panel's "
                  f"options at {on:.1%}, the same with fill off at {off:.1%}")
            if nf != view.frame_count() or dur != 250 or not same_order:
                failures.append(f"H. Studio wrote {nf} frames at {dur} ms, order kept {same_order}")
            if not (on >= 0.97 and off < 0.9):
                failures.append(f"H. Studio frame vs panel options {on:.1%}, vs fill off {off:.1%}")
        finally:
            mw.doc._dirty = False
            mw.close()
    return failures


if __name__ == "__main__":
    f = run()
    print("\n".join(["FAILED:"] + f) if f else "\nAll flow trace checks passed.")
    sys.exit(1 if f else 0)
