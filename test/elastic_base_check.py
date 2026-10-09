"""A slip surface that rides along the underside of an elastic zone.

An `elastic` material is impenetrable in the LEM: a trial surface may run along
its boundary but not into it (the crossing check in slice.generate_slices). A
wall drawn as an elastic polygon sits ON its foundation, so a surface along the
wall's base rides the zone's underside, and every slice there finds the wall as
its deepest present layer. The slice must bind the strength of the material the
wall sits on, not stop the run on the wall's missing strength option.

  A. A non-circular surface along the base of Tutorial FEM-3's block wall, then
     up through the fill: generate_slices accepts it, and every slice under the
     wall takes the foundation as its base material.
  B. The same surface lifted 1 mm into the wall, inside the crossing check's
     tolerance: the same.
  C. A circular search on that wall completes with a finite factor of safety.
     It used to stop with a ValueError when a trial circle grazed the blocks,
     instead of scoring that circle as a failed trial.
"""
import math
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "docs/tutorials/files/xslope_block_wall.xlsx"
X_FACE, X_BACK = 8.0, 9.2          # the blocks' face and back, m


def _surface(y_base):
    return [{"X": X_FACE, "Y": 0.0, "Movement": "Fixed"},
            {"X": X_BACK, "Y": y_base, "Movement": "Fixed"},
            {"X": 11.0, "Y": 3.6, "Movement": "Fixed"},
            {"X": 12.0, "Y": 5.1, "Movement": "Fixed"}]


def run():
    from xslope.fileio import load_slope_data
    from xslope.search import circular_search
    from xslope.slice import generate_slices

    warnings.filterwarnings("ignore")
    failures = []
    data = load_slope_data(str(BOOK))
    names = [m["name"] for m in data["materials"]]
    if "block" not in names or data["materials"][names.index("block")]["option"] != "elastic":
        return [f"{BOOK.name}: expected an elastic 'block' material"]

    for tag, y_base in (("A", 0.0), ("B", 0.001)):
        try:
            ok, res = generate_slices(data, non_circ=_surface(y_base), num_slices=30)
        except Exception as exc:                              # noqa: BLE001
            failures.append(f"{tag}: slicing stopped: {type(exc).__name__}: {exc}")
            continue
        if not ok:
            failures.append(f"{tag}: surface along the wall base refused: {res}")
            continue
        df, _ = res
        under = df[(df["x_c"] > X_FACE) & (df["x_c"] < X_BACK)]
        got = sorted({names[int(m) - 1] for m in under["mat"]})
        if under.empty or got != ["foundation"]:
            failures.append(f"{tag}: slices under the wall bind {got}, not ['foundation']")

    data = load_slope_data(str(BOOK))
    data["circles"] = [{"Xo": 9.0, "Yo": 9.0, "Depth": -1.0, "R": 10.0}]
    try:
        out = circular_search(data, "spencer", num_slices=30)
        best = (out[0] if isinstance(out, tuple) else out)[0]
        if not math.isfinite(best["FS"]) or best["FS"] >= 9999:
            failures.append(f"C: search returned FS = {best['FS']}")
    except Exception as exc:                                  # noqa: BLE001
        failures.append(f"C: search stopped: {type(exc).__name__}: {exc}")
    return failures


if __name__ == "__main__":
    f = run()
    print("\n".join(f) if f else "elastic base check: all passed")
