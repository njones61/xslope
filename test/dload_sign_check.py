"""A surface load's horizontal component must push the same way in every method.

XSLOPE hands every solution method the same surface load: one resultant ``D``
at one point, at one inclination ``beta``, carrying the force

    (D·sin beta,  −D·cos beta)

in the sliding-direction frame that ``generate_slices`` normalizes to (a
right-facing surface has its ``alpha`` and ``beta`` negated, so a load normal to
the slope face reads ``beta > 0`` whichever way the slope faces). The horizontal
part of that force is therefore a SIGNED component, not a magnitude in the
driving direction the way the seismic force ``kW`` and the tension-crack water
force ``T`` are: ``beta > 0`` pushes INTO the slope and resists sliding —
reservoir thrust on the upstream face of a dam is the everyday case — while
``beta < 0`` pushes downslope and drives it.

Line loads follow the same convention (``lload`` at ``ll_beta``).

Nothing about that is a matter of method. So this check is a two-direction
diagnostic rather than a locked number: the same model, the same load resultant
at the same point, applied once at ``+beta`` and once at ``−beta``, must give a
LOWER factor of safety at ``−beta`` than at ``+beta`` in all seven methods, on a
left-facing model and on a right-facing one. Whatever a method's absolute
answer, the load cannot be stabilising in one method and destabilising in
another.

It is run on both facings because the normalization above is the only thing that
makes the sign facing-independent, and a diagnostic that reads only one facing
would not notice if it stopped working.

The second check is the consequence at model scale: on the dam sample —
submerged upstream face, a large reservoir load — Janbu must land in the same
neighbourhood as the complete-equilibrium methods rather than a fraction of
them.

Run directly:  PYTHONPATH=. python3 test/dload_sign_check.py
"""

import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: (left-facing model, right-facing model). Both are homogeneous c-phi slopes on
#: which all seven methods converge with a surface load either way up.
MODELS = (("docs/inputs/slope/xslope_lface.xlsx", False),
          ("docs/inputs/slope/xslope_rface.xlsx", True))

DAM = "docs/inputs/slope/xslope_dam.xlsx"

#: Inclination the probe load is applied at, either side of the surface normal.
BETA_DEG = 60.0
#: Probe resultant as a fraction of the sliding mass weight. Small enough that
#: every method still converges, large enough that the sign is unmistakable.
LOAD_FRACTION = 0.02
#: The probe sits on the upper part of the surface, clear of both ends.
SPAN = (0.55, 0.85)


def _methods():
    from xslope.solve import bishop, corps, janbu, lowe, mprice, oms, spencer
    return (("oms", oms), ("bishop", bishop), ("janbu", janbu),
            ("corps", corps), ("lowe", lowe), ("spencer", spencer),
            ("mprice", mprice))


def _slices(xlsx, num_slices=40):
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices
    slope_data = load_slope_data(os.path.join(_REPO, xlsx))
    ok, out = generate_slices(slope_data, circle=slope_data["circles"][0],
                              num_slices=num_slices)
    if not ok:
        return None, out
    return out[0], None


def _solve_all(df):
    out = {}
    for name, fn in _methods():
        ok, res = fn(df.copy())
        out[name] = float(res["FS"]) if ok else None
    return out


def _with_probe(df, beta_deg, line_load=False):
    """The model with one probe load resultant on its upper surface.

    Applied straight into the solver's own columns rather than through a load
    line in the workbook, because the two directions have to differ in NOTHING
    but the sign of the inclination: the same magnitude, at the same point, on
    the same geometry. A load line steep enough to point downslope would have to
    sit on different ground.
    """
    n = len(df)
    lo, hi = int(SPAN[0] * n), int(SPAN[1] * n)
    per = LOAD_FRACTION * float(df["w"].sum()) / max(1, hi - lo)
    mag = np.zeros(n)
    ang = np.zeros(n)
    mag[lo:hi] = per
    ang[lo:hi] = beta_deg
    work = df.copy()
    if line_load:
        work["lload"] = mag
        work["ll_beta"] = ang
        work["ll_x"] = df["x_c"].values
        work["ll_y"] = df["y_ct"].values
    else:
        work["dload"] = mag
        work["beta"] = ang
        work["d_x"] = df["x_c"].values
        work["d_y"] = df["y_ct"].values
    return work


def test_two_directions():
    """A load leaning downslope lowers FS; the same load leaning into the slope
    raises it. In every method, on both facings, for distributed and line loads."""
    fails = []
    for xlsx, want_right_facing in MODELS:
        df, why = _slices(xlsx)
        if df is None:
            fails.append(f"{xlsx}: generate_slices failed: {why}")
            continue
        facing = bool(df["y_lb"].iat[0] > df["y_rb"].iat[-1])
        if facing != want_right_facing:
            fails.append(f"{xlsx}: expected right_facing={want_right_facing}, "
                         f"read {facing} — the probe is on the wrong frame")
            continue
        for kind, line_load in (("distributed load", False), ("line load", True)):
            down = _solve_all(_with_probe(df, -BETA_DEG, line_load))
            into = _solve_all(_with_probe(df, +BETA_DEG, line_load))
            for name, _fn in _methods():
                a, b = down[name], into[name]
                if a is None or b is None:
                    fails.append(f"{xlsx} {kind}: {name} did not solve "
                                 f"(downslope={a}, into-slope={b})")
                    continue
                if a >= b:
                    fails.append(
                        f"{xlsx} {kind}: {name} makes a load leaning DOWNSLOPE "
                        f"(beta=−{BETA_DEG:g}) safer than the same load leaning "
                        f"INTO the slope (beta=+{BETA_DEG:g}): FS {a:.4f} vs "
                        f"{b:.4f} — its horizontal component carries the wrong "
                        f"sign")
    return fails


def test_submerged_face():
    """Janbu on the dam sample lands with the other methods, not far below them.

    The upstream face carries a full reservoir as a perpendicular distributed
    load; its horizontal component is the largest single term in the horizontal
    balance. Read against Spencer, which satisfies both force and moment
    equilibrium. The band is wide on purpose — Janbu is an approximate method
    and is expected to differ — it is there to catch a sign, not to pin a value.
    """
    fails = []
    df, why = _slices(DAM)
    if df is None:
        return [f"{DAM}: generate_slices failed: {why}"]
    beta = np.radians(df["beta"].values)
    thrust = float(np.sum(df["dload"].values * np.sin(beta)))
    if thrust <= 0.0:
        return [f"{DAM}: the reservoir load has no net into-slope horizontal "
                f"component (sum D·sin beta = {thrust:.1f}), so the sample no "
                f"longer exercises what this check is for"]
    fs = _solve_all(df)
    if fs["janbu"] is None or fs["spencer"] is None:
        return [f"{DAM}: janbu={fs['janbu']}, spencer={fs['spencer']}"]
    ratio = fs["janbu"] / fs["spencer"]
    if not 0.7 <= ratio <= 1.3:
        fails.append(f"{DAM}: janbu FS {fs['janbu']:.4f} is {ratio:.2f}x "
                     f"spencer's {fs['spencer']:.4f} on a submerged face — the "
                     f"reservoir load is not being balanced the same way")
    return fails


CHECKS = [
    ("a surface load pushes the same way in every method", test_two_directions),
    ("Janbu on a submerged face", test_submerged_face),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    failures = []
    for name, fn in CHECKS:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:48s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("Surface-load sign:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nEvery method takes the surface load's horizontal component the same way.")


if __name__ == "__main__":
    main()
