# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Checks for the FEM reinforcement / pile solution-details dialog.

What is being defended:

  A. THE GATE — the "1D Details…" button sits on the FEM results view's toolbar
     and is enabled only when the solution carries 1D elements. With none it
     stays visible and dims, and its tooltip says why.

  B. THE LIST — every reinforcement line and pile appears, under its own header,
     with a utilization badge whose colour follows the reported ratio.

  C. THE ENVELOPE IS THE SOLVER'S — the capacity envelope drawn over a
     reinforcement line's force profile is the same curve the solver used to fill
     ``t_allow_by_1d_elem``, sampled at more points. Asserted at every element
     centroid, and mutation-tested: replacing the shared capacity function moves
     the drawn envelope with it, which a private copy in the plotting code could
     not do.

  D. THE RELOAD PATH — a solution written to its sidecars and read back produces
     the same profiles, series by series. The dialog must not be able to tell a
     fresh solve from a reopened file.

  E. THE EXPORT — the Export button writes a PNG and a CSV, and the CSV's
     columns are the series that were plotted.

The two sample models are solved once each and shared across the checks. Both
are solved with a reduced iteration budget: nothing here reads a converged
value, only the structure and identity of the profiles, and the capacities the
envelope check compares against are properties of the model rather than of the
solve.
"""

import os
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np                                              # noqa: E402

REINF_XLSX = os.path.join(_REPO, "docs", "fem", "files", "xslope_reinforce_fem.xlsx")
PILES_XLSX = os.path.join(_REPO, "docs", "fem", "files", "xslope_piles_fem.xlsx")

_SOLVED = {}


def _solved(path):
    """Solve one sample model once and cache it: (slope_data, fem_data, solution)."""
    if path in _SOLVED:
        return _SOLVED[path]
    from xslope.fem import build_fem_data, solve_fem
    from xslope.fileio import load_slope_data
    slope_data = load_slope_data(path)
    fem_data = build_fem_data(slope_data, slope_data.get("mesh"))
    solution = solve_fem(fem_data, F=1.0, debug_level=0, max_iterations=60)
    _SOLVED[path] = (slope_data, fem_data, solution)
    return _SOLVED[path]


def _app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


# --------------------------------------------------------------------------
# A. the gate
# --------------------------------------------------------------------------

def test_gate():
    """The toolbar button appears on the FEM results view, enabled with 1D
    elements and dimmed with a reason without them."""
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from studio.main_window import MainWindow

    slope_data, fem_data, solution = _solved(REINF_XLSX)
    mw = MainWindow()
    try:
        mw.doc.results["fem_solution"] = {
            "fem_data": fem_data, "solution": solution, "FS": None,
            "analysis": "single", "failure_solution": None}
        mw._show_fem_results()

        btn = mw.fem_details_btn
        if btn is None:
            return ["no '1D Details…' button on the FEM results toolbar"]
        if not btn.isVisibleTo(mw.fem_results_canvas):
            fails.append("the details button is not on the FEM results toolbar")
        if "Details" not in btn.text():
            fails.append(f"unexpected button text {btn.text()!r}")
        if not btn.isEnabled():
            fails.append("the details button is disabled for a model WITH 1D elements")

        # A solution with no 1D elements: still there, dimmed, and it says why.
        bare = dict(fem_data)
        bare["elements_1d"] = np.zeros((0, 3), dtype=int)
        mw.doc.results["fem_solution"] = {
            "fem_data": bare, "solution": solution, "FS": None,
            "analysis": "single", "failure_solution": None}
        mw._update_fem_details_action()
        if btn.isEnabled():
            fails.append("the details button stays enabled with no 1D elements")
        tip = btn.toolTip().lower()
        if "reinforcement" not in tip or "pile" not in tip:
            fails.append(f"the dimmed button's tooltip gives no reason: {btn.toolTip()!r}")

        # And it comes back when a 1D solution is restored — the gate is live,
        # not a one-shot decision made when the view was built.
        mw.doc.results["fem_solution"] = {
            "fem_data": fem_data, "solution": solution, "FS": None,
            "analysis": "single", "failure_solution": None}
        mw._update_fem_details_action()
        if not btn.isEnabled():
            fails.append("the details button does not re-enable when 1D results return")
    finally:
        mw.close()
    return fails


def test_gate_mutation():
    """A mutation of the gate must be caught: forcing has_1d_details to True
    should make the 'dimmed with no 1D elements' leg of test_gate fail."""
    fails = []
    import xslope.fem_details as fd
    real = fd.has_1d_details
    fd.has_1d_details = lambda fem_data: True
    try:
        if not test_gate():
            fails.append("test_gate passes with the 1D-element gate forced open — "
                         "it does not actually check the dimmed case")
    finally:
        fd.has_1d_details = real
    return fails


# --------------------------------------------------------------------------
# B. the list
# --------------------------------------------------------------------------

def test_list():
    """Both samples populate the list: a header per kind, one row per member,
    each row carrying a badge that matches its reported utilization."""
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from PySide6.QtCore import Qt
    from studio.fem_details_dialog import FemDetailsDialog
    from xslope import fem_details as fd

    for path, kind, want_n in ((REINF_XLSX, "reinforcement", 6), (PILES_XLSX, "pile", 2)):
        slope_data, fem_data, solution = _solved(path)
        dlg = FemDetailsDialog(fem_data, solution, slope_data, model_path=path)
        try:
            entries = dlg.entries()
            of_kind = [e for e in entries if e["kind"] == kind]
            if len(of_kind) != want_n:
                fails.append(f"{os.path.basename(path)}: expected {want_n} {kind} rows, "
                             f"got {len(of_kind)}")
            # One header per kind present, and it is not selectable.
            headers = [dlg.list.item(r) for r in range(dlg.list.count())
                       if dlg.list.item(r).data(Qt.UserRole) is None]
            names = {h.text() for h in headers}
            want_header = "Reinforcement" if kind == "reinforcement" else "Piles"
            if want_header not in names:
                fails.append(f"{os.path.basename(path)}: no '{want_header}' header "
                             f"(headers: {sorted(names)})")
            for h in headers:
                if h.flags() != Qt.NoItemFlags:
                    fails.append(f"header row {h.text()!r} is selectable")
            if dlg.list.count() != len(entries) + len(headers):
                fails.append(f"{os.path.basename(path)}: list has {dlg.list.count()} rows "
                             f"for {len(entries)} members and {len(headers)} headers")
            # Badge colour follows the reported ratio.
            for e in entries:
                want = ("none" if e["utilization"] is None else
                        "red" if e["utilization"] >= fd.UTIL_AT_CAPACITY else
                        "amber" if e["utilization"] >= fd.UTIL_WATCH else "green")
                if e["badge"] not in (want, "red"):   # pullout forces red
                    fails.append(f"{e['label']}: badge {e['badge']!r} for "
                                 f"utilization {e['utilization']}")
            # A member is selected and its detail is drawn.
            if dlg.current_profile() is None:
                fails.append(f"{os.path.basename(path)}: no member selected on open")
        finally:
            dlg.close()
    return fails


# --------------------------------------------------------------------------
# C. the capacity envelope
# --------------------------------------------------------------------------

def _envelope_at(profile, s):
    return float(np.interp(s, profile["env_s"], profile["env_T"]))


def test_envelope():
    """The drawn envelope equals the solver's per-element capacity, everywhere
    the solver has one — including the ramp, the plateau, and the far ramp."""
    fails = []
    from xslope import fem_details as fd

    slope_data, fem_data, solution = _solved(REINF_XLSX)
    checked = 0
    for line_id in (1, 4, 6):
        prof = fd.reinforcement_profile(fem_data, solution, line_id, slope_data)
        if prof["env_s"] is None:
            fails.append(f"line {line_id}: no analytic envelope was built")
            continue
        if len(prof["s"]) < 3:
            fails.append(f"line {line_id}: only {len(prof['s'])} probe points")
        for s, t_allow in zip(prof["s"], prof["t_allow"]):
            got = _envelope_at(prof, s)
            if abs(got - t_allow) > 1e-6 * max(1.0, abs(t_allow)):
                fails.append(f"line {line_id} at s={s:.3f}: envelope {got:.6f} != "
                             f"solver capacity {t_allow:.6f}")
            checked += 1
        # The envelope must actually have the three regions, or the equality
        # above could be passing on a flat line.
        if prof["env_T"] is not None:
            lo = _envelope_at(prof, 0.0)
            mid = _envelope_at(prof, prof["length"] / 2.0)
            hi = _envelope_at(prof, prof["length"])
            if not (lo < mid and hi < mid):
                fails.append(f"line {line_id}: envelope is not a ramp-plateau-ramp "
                             f"(ends {lo:.1f}/{hi:.1f}, middle {mid:.1f})")
    if checked < 3:
        fails.append(f"only {checked} envelope probe points were compared")
    return fails


def test_envelope_mutation():
    """The envelope is the shared capacity function, not a copy of it.

    ``xslope.fileio.reinforce_available_tension`` is replaced and the detail
    module re-imported, exactly as a fresh interpreter would take it. Both the
    solver's element capacities and the drawn envelope move to the mutated rule
    and still agree. A reimplementation inside the detail or plotting code would
    keep the original shape while the solver's capacities halved, and the
    equality would break — which is the failure this check exists to see.
    """
    fails = []
    import importlib
    import xslope.fem_details as fd
    import xslope.fileio as fileio
    from xslope.fem import build_fem_data

    slope_data, _fem_data, _solution = _solved(REINF_XLSX)
    base = fd.reinforcement_profile(_fem_data, {}, 4, slope_data)

    real = fileio.reinforce_available_tension

    def mutated(d1, d2, t_max, lp1, lp2, tend1=0.0, tend2=0.0):
        return 0.5 * real(d1, d2, t_max, lp1, lp2, tend1, tend2)

    fileio.reinforce_available_tension = mutated
    try:
        fd = importlib.reload(fd)
        fem_data = build_fem_data(slope_data, slope_data.get("mesh"))
        prof = fd.reinforcement_profile(fem_data, {}, 4, slope_data)
        if prof["env_s"] is None:
            return ["mutation run built no envelope"]
        for s, t_allow in zip(prof["s"], prof["t_allow"]):
            got = _envelope_at(prof, s)
            if abs(got - t_allow) > 1e-6 * max(1.0, abs(t_allow)):
                fails.append(f"under mutation, envelope {got:.6f} != solver "
                             f"capacity {t_allow:.6f} at s={s:.3f} — the plot has "
                             f"its own copy of the capacity rule")
        # And the mutation actually landed, so the agreement above means
        # something: the envelope is not simply unchanged.
        moved = max(abs(_envelope_at(prof, s) - _envelope_at(base, s))
                    for s in base["s"])
        if moved < 1e-6:
            fails.append("the mutation did not change the envelope — it does not "
                         "test the wire")
    finally:
        fileio.reinforce_available_tension = real
        importlib.reload(fd)
    return fails


# --------------------------------------------------------------------------
# D. the reload path
# --------------------------------------------------------------------------

_SKIP_KEYS = {"element_ids", "pile_indices", "node_ids", "mechanism"}


def _same_profile(a, b, label, fails):
    keys = sorted(set(a) | set(b))
    for k in keys:
        if k in _SKIP_KEYS:
            continue
        va, vb = a.get(k), b.get(k)
        if isinstance(va, np.ndarray) or isinstance(vb, np.ndarray):
            va, vb = np.asarray(va, dtype=float), np.asarray(vb, dtype=float)
            if va.shape != vb.shape:
                fails.append(f"{label}.{k}: shape {va.shape} != {vb.shape} after reload")
            elif not np.allclose(va, vb, rtol=1e-9, atol=1e-9, equal_nan=True):
                worst = float(np.nanmax(np.abs(va - vb)))
                fails.append(f"{label}.{k}: differs after reload (max |d| = {worst:.3g})")
        elif isinstance(va, float) and isinstance(vb, float):
            if not (np.isnan(va) and np.isnan(vb)) and abs(va - vb) > 1e-9:
                fails.append(f"{label}.{k}: {va} != {vb} after reload")
        elif va != vb:
            fails.append(f"{label}.{k}: {va!r} != {vb!r} after reload")


def test_reload():
    """Profiles from a solution written to its sidecars and read back are
    identical, series by series, to the ones from the fresh solve."""
    fails = []
    from xslope.fem import export_fem_solution, import_fem_solution
    from xslope import fem_details as fd
    from pathlib import Path

    for path, kind in ((REINF_XLSX, "reinforcement"), (PILES_XLSX, "pile")):
        slope_data, fem_data, solution = _solved(path)
        with tempfile.TemporaryDirectory() as tmp:
            stem = Path(tmp) / "reloaded"
            export_fem_solution(fem_data, solution, stem)
            back = import_fem_solution(fem_data, stem)

            if kind == "reinforcement":
                for line_id in (1, 4, 6):
                    a = fd.reinforcement_profile(fem_data, solution, line_id, slope_data)
                    b = fd.reinforcement_profile(fem_data, back, line_id, slope_data)
                    _same_profile(a, b, f"line {line_id}", fails)
            else:
                for pidx in (0, 1):
                    a = fd.pile_profile(fem_data, solution, pidx, slope_data)
                    b = fd.pile_profile(fem_data, back, pidx, slope_data)
                    _same_profile(a, b, f"pile {pidx}", fails)

            # The reinforcement sidecar now carries the cap the solve enforced,
            # and an older file without that column still reads.
            if kind == "reinforcement":
                import pandas as pd
                csv = stem.parent / f"{stem.name}_fem_reinf.csv"
                df = pd.read_csv(csv, comment="#")
                if "t_cap" not in df.columns:
                    fails.append("the reinforcement sidecar has no t_cap column")
                old = stem.parent / "old"
                df.drop(columns=["t_cap"]).to_csv(
                    stem.parent / "old_fem_reinf.csv", index=False)
                for suffix in ("nodes", "elements"):
                    (stem.parent / f"old_fem_{suffix}.csv").write_bytes(
                        (stem.parent / f"{stem.name}_fem_{suffix}.csv").read_bytes())
                try:
                    legacy = import_fem_solution(fem_data, old)
                except Exception as exc:
                    fails.append(f"a sidecar without t_cap fails to load: {exc!r}")
                else:
                    got = np.asarray(legacy.get("t_cap_1d"), dtype=float)
                    want = np.asarray(fem_data["t_allow_by_1d_elem"], dtype=float)
                    if got.shape != want.shape or not np.allclose(got, want):
                        fails.append("a sidecar without t_cap does not fall back "
                                     "to t_allow")
    return fails


def test_dialog_reload_identical():
    """The dialog itself is indifferent to where the solution came from: the
    same member list and the same drawn profile from a reloaded solution."""
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from pathlib import Path
    from studio.fem_details_dialog import FemDetailsDialog
    from xslope.fem import export_fem_solution, import_fem_solution

    slope_data, fem_data, solution = _solved(REINF_XLSX)
    with tempfile.TemporaryDirectory() as tmp:
        stem = Path(tmp) / "reloaded"
        export_fem_solution(fem_data, solution, stem)
        back = import_fem_solution(fem_data, stem)

        fresh = FemDetailsDialog(fem_data, solution, slope_data, model_path=REINF_XLSX)
        again = FemDetailsDialog(fem_data, back, slope_data, model_path=REINF_XLSX)
        try:
            if [e["label"] for e in fresh.entries()] != [e["label"] for e in again.entries()]:
                fails.append("the member list differs between a fresh and a "
                             "reloaded solution")
            for a, b in zip(fresh.entries(), again.entries()):
                if a["badge"] != b["badge"]:
                    fails.append(f"{a['label']}: badge {a['badge']} -> {b['badge']} "
                                 f"after reload")
            _same_profile(fresh.current_profile(), again.current_profile(),
                          "dialog", fails)
        finally:
            fresh.close()
            again.close()
    return fails


# --------------------------------------------------------------------------
# E. the export
# --------------------------------------------------------------------------

def test_export():
    """Export writes a PNG and a CSV beside each other, and the CSV's columns
    are the plotted series."""
    fails = []
    import csv as csvmod
    import matplotlib
    matplotlib.use("Agg")
    app = _app()
    from studio.fem_details_dialog import FemDetailsDialog
    from xslope import fem_details as fd

    for path, needed in (
            (REINF_XLSX, ("position", "axial_force", "capacity", "envelope_position",
                          "envelope_capacity", "bond_position", "bond_transfer")),
            (PILES_XLSX, ("node_depth", "lateral_displacement", "shear", "moment",
                          "soil_reaction"))):
        slope_data, fem_data, solution = _solved(path)
        dlg = FemDetailsDialog(fem_data, solution, slope_data, model_path=path)
        dlg.show()
        for _ in range(12):
            app.processEvents()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                target = os.path.join(tmp, "detail.png")
                png, csv_path = dlg.export_current(path=target)
                if not png or not os.path.exists(png) or os.path.getsize(png) < 5000:
                    fails.append(f"{os.path.basename(path)}: no usable PNG written")
                if not csv_path or not os.path.exists(csv_path):
                    fails.append(f"{os.path.basename(path)}: no CSV written")
                    continue
                with open(csv_path) as f:
                    rows = list(csvmod.reader(f))
                header = rows[0] if rows else []
                for want in needed:
                    if not any(c.split(" (")[0] == want for c in header):
                        fails.append(f"{os.path.basename(path)}: CSV has no "
                                     f"{want!r} column (got {header})")
                if len(rows) < 2:
                    fails.append(f"{os.path.basename(path)}: CSV has no data rows")
                # The CSV is the plotted data, not a re-derivation.
                cols, table = fd.profile_table(dlg.current_profile())
                if header != cols or rows[1:] != table:
                    fails.append(f"{os.path.basename(path)}: the CSV is not the "
                                 f"profile table")
                # The default name carries the model and the member.
                stem = dlg.default_export_stem()
                model = os.path.splitext(os.path.basename(path))[0]
                if not stem.startswith(model):
                    fails.append(f"default export name {stem!r} does not start "
                                 f"with the model name")
        finally:
            dlg.close()
    return fails


CHECKS = [
    ("the toolbar button and its gate", test_gate),
    ("the gate check would catch a mutation", test_gate_mutation),
    ("the member list and its badges", test_list),
    ("the envelope is the solver's capacity", test_envelope),
    ("the envelope is the shared function", test_envelope_mutation),
    ("profiles survive save + reload", test_reload),
    ("the dialog is the same on a reloaded solution", test_dialog_reload_identical),
    ("export writes the PNG and the profile CSV", test_export),
]

# Checks that need the Studio layer; skipped when PySide6 is absent.
_STUDIO_ONLY = {test_gate, test_gate_mutation, test_list,
                test_dialog_reload_identical, test_export}


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    checks = CHECKS
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("FEM 1D details: PySide6 not installed — Studio checks skipped.")
        checks = [c for c in CHECKS if c[1] not in _STUDIO_ONLY]

    failures = []
    for name, fn in checks:
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
    print("FEM 1D solution details:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll FEM 1D detail checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app_ = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
