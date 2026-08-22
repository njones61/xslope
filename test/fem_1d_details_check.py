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
     with a utilization badge whose colour follows the member's verdict, and a
     reinforcement row spells that verdict out in the one vocabulary
     ``xslope.fem_details.REINFORCEMENT_STATES`` holds.

  C. THE ENVELOPE IS THE SOLVER'S — the capacity envelope drawn over a
     reinforcement line's force profile is the same curve the solver used to fill
     ``t_allow_by_1d_elem``, sampled at more points. Asserted at every element
     centroid, and mutation-tested: replacing the shared capacity function moves
     the drawn envelope with it, which a private copy in the plotting code could
     not do.

  D. THE RELOAD PATH — a solution written to its sidecars and read back produces
     the same profiles, series by series. The dialog must not be able to tell a
     fresh solve from a reopened file.

  E. THE EXPORT — the Export button writes a PNG and a CSV, the CSV's columns
     are the series that were plotted, and both files say which field state
     they were taken at.

  F. THE FIELD STATE — the panel carries the results view's converged /
     at-failure switch. It is dimmed, and held on the converged field, when the
     run captured no at-failure snapshot; with one, it reads the snapshot's own
     member forces, which are not the converged ones. The snapshot is found
     whether the caller passes it alongside the solution (Studio's bundle) or it
     arrives nested inside a solution read back from disk.

The two sample models are solved once each and shared across the checks. Both
are solved with a reduced iteration budget: nothing here reads a converged
value, only the structure and identity of the profiles, and the capacities the
envelope check compares against are properties of the model rather than of the
solve.
"""

import contextlib
import io
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
_FAILURE = {}

# The strength reduction each sample is pushed to for its at-failure snapshot.
# Both samples stand at their own F well below this (1.49 and 1.36), so the
# trial is past the critical one: the field it leaves is the unconverged
# collapse mechanism, which is exactly what solve_ssrm captures and hands on as
# ``result['failure_solution']``.
_FAILURE_F = 2.5


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


def _failure_field(path):
    """An at-failure snapshot for one sample: the same shape solve_ssrm captures,
    a solve at a strength reduction past the model's own factor of safety."""
    if path in _FAILURE:
        return _FAILURE[path]
    from xslope.fem import solve_fem
    _slope_data, fem_data, _solution = _solved(path)
    _FAILURE[path] = solve_fem(fem_data, F=_FAILURE_F, debug_level=0,
                               max_iterations=60)
    return _FAILURE[path]


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
    each row carrying the member's verdict in words and a badge that matches
    it."""
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from PySide6.QtCore import Qt
    from studio.fem_details_dialog import FemDetailsDialog, _row_text
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
            # Badge colour follows the member's verdict, and the row carries
            # that verdict in words: two lines both standing at 100% are told
            # apart by the word, not by the dot or the percentage.
            for e in entries:
                want = fd.reinforcement_badge(e.get("status_key"),
                                              e["utilization"])
                if e["badge"] != want:
                    fails.append(f"{e['label']}: badge {e['badge']!r} for "
                                 f"{e.get('status_key')!r} at utilization "
                                 f"{e['utilization']}, expected {want!r}")
                if e["kind"] == "reinforcement":
                    if e["status_key"] not in fd.REINFORCEMENT_STATES:
                        fails.append(f"{e['label']}: unknown verdict "
                                     f"{e.get('status_key')!r}")
                    if e["status"] != fd.reinforcement_state_phrase(
                            e["status_key"]):
                        fails.append(f"{e['label']}: the row's status "
                                     f"{e['status']!r} is not the phrase for "
                                     f"{e['status_key']!r}")
                    row = _row_text(e)
                    if e["status"] not in row:
                        fails.append(f"{e['label']}: the list row {row!r} does "
                                     f"not carry the verdict {e['status']!r}")
            # A member is selected and its detail is drawn.
            if dlg.current_profile() is None:
                fails.append(f"{os.path.basename(path)}: no member selected on open")
        finally:
            dlg.close()
    return fails


def test_the_verdict_is_the_summary_s():
    """The word on a list row is the word the printed summary prints.

    The two used to be different vocabularies over the same arrays: the summary
    told a line slipping at its embedment-limited capacity near an end
    (PULLOUT) from one standing at its full tensile capacity in the middle
    (YIELDED), and the panel called both of them "at capacity" — so the panel
    could not answer the question a reader opens it to ask. One function
    decides it now, and this check compares what each surface says, line by
    line, rather than trusting that they call the same one.
    """
    fails = []
    from xslope import fem_details as fd
    from xslope.fem import print_reinforcement_summary

    slope_data, fem_data, solution = _solved(REINF_XLSX)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print_reinforcement_summary(fem_data, solution)
    printed = {}
    for raw in buf.getvalue().splitlines():
        bits = raw.split()
        if len(bits) >= 9 and bits[0].isdigit():
            printed[int(bits[0])] = " ".join(bits[8:])
    if not printed:
        return ["the summary printed no line rows to compare against"]

    rows = {e["index"]: e for e in fd.list_lines(fem_data, solution, slope_data)
            if e["kind"] == "reinforcement"}
    for line_id, word in printed.items():
        row = rows.get(line_id)
        if row is None:
            fails.append(f"the summary prints line {line_id} and the list "
                         f"carries no such row")
            continue
        if word != row["status"].upper():
            fails.append(f"line {line_id}: the summary says {word!r} and the "
                         f"panel says {row['status']!r}")
    # The distinction itself is reachable: classify a made-up line whose only
    # at-capacity element sits inside a pullout ramp, and one whose only
    # at-capacity element sits at the full Tmax, and the two must differ.
    ramp = fd.reinforcement_status([100.0, 400.0, 200.0],
                                   [200.0, 800.0, 200.0],
                                   failed=[False, False, True])
    full = fd.reinforcement_status([100.0, 800.0, 100.0],
                                   [200.0, 800.0, 200.0],
                                   failed=[False, True, False])
    if ramp[0] != "pullout" or full[0] != "yielded":
        fails.append(f"the at-capacity split does not discriminate: a ramp end "
                     f"reads {ramp[0]!r} and a full-Tmax element {full[0]!r}")
    # And a line at capacity in BOTH places is reported yielded, which is the
    # precedence the summary has always used.
    both = fd.reinforcement_status([100.0, 800.0, 200.0],
                                   [200.0, 800.0, 200.0],
                                   failed=[False, True, True])
    if both[0] != "yielded":
        fails.append(f"a line at capacity in a ramp AND at full Tmax reads "
                     f"{both[0]!r}, not yielded")
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
        elif (isinstance(va, tuple) and isinstance(vb, tuple)
              and len(va) == len(vb)):
            # A pair of numbers — a span — is two floats and is compared as
            # two floats. Exact equality is the wrong test for anything that
            # went through a file: the reload path writes decimal.
            if not all(abs(float(x) - float(y)) <= 1e-9 for x, y in zip(va, vb)):
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

def _read_export_csv(path):
    """An exported profile CSV split into its ``#`` provenance lines and its
    table (header row + data rows), the way a reader takes it."""
    import csv as csvmod
    with open(path) as f:
        text = f.read().splitlines(keepends=True)
    comments = [ln for ln in text if ln.startswith("#")]
    rows = list(csvmod.reader(ln for ln in text if not ln.startswith("#")))
    return comments, rows


def test_export():
    """Export writes a PNG and a CSV beside each other, the CSV's columns are
    the plotted series, and both files record the field state."""
    fails = []
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
                comments, rows = _read_export_csv(csv_path)
                header = rows[0] if rows else []
                # The state the profile was read at, in the file itself.
                if not any("field state" in c and "last converged" in c
                           for c in comments):
                    fails.append(f"{os.path.basename(path)}: the CSV does not "
                                 f"record its field state (comments: {comments})")
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
                # The default name carries the model, the member and the state.
                stem = dlg.default_export_stem()
                model = os.path.splitext(os.path.basename(path))[0]
                if not stem.startswith(model):
                    fails.append(f"default export name {stem!r} does not start "
                                 f"with the model name")
                if not stem.endswith("converged"):
                    fails.append(f"default export name {stem!r} does not name "
                                 f"the field state")
        finally:
            dlg.close()
    return fails


# --------------------------------------------------------------------------
# F. the field state switch
# --------------------------------------------------------------------------

def test_field_state_control():
    """The switch is on the panel with the results view's own two states, and it
    is dimmed — and held on the converged field — with no at-failure snapshot to
    switch to."""
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from studio.fem_details_dialog import FemDetailsDialog

    slope_data, fem_data, solution = _solved(REINF_XLSX)

    dlg = FemDetailsDialog(fem_data, solution, slope_data, model_path=REINF_XLSX)
    try:
        combo = getattr(dlg, "field_state", None)
        if combo is None:
            return ["the details panel has no Field state control"]
        labels = [combo.itemText(i) for i in range(combo.count())]
        keys = [combo.itemData(i) for i in range(combo.count())]
        if labels != ["At failure", "Last converged"]:
            fails.append(f"Field state offers {labels} — the results view offers "
                         f"['At failure', 'Last converged']")
        if keys != ["failure", "converged"]:
            fails.append(f"Field state keys {keys} do not match the results "
                         f"view's ('failure' / 'converged')")
        if combo.isEnabled():
            fails.append("Field state is enabled with no at-failure snapshot")
        if dlg.current_field_state() != "converged":
            fails.append(f"with no snapshot the panel reads "
                         f"{dlg.current_field_state()!r}, not the converged field")
        if dlg.current_profile().get("field_state") != "converged":
            fails.append("the drawn profile does not report the converged state")
    finally:
        dlg.close()

    # With a snapshot the control is live, and it opens on the state the results
    # view opens on.
    failure = _failure_field(REINF_XLSX)
    dlg = FemDetailsDialog(fem_data, solution, slope_data, model_path=REINF_XLSX,
                           failure_solution=failure)
    try:
        if not dlg.field_state.isEnabled():
            fails.append("Field state stays dimmed with an at-failure snapshot")
        if dlg.current_field_state() != "failure":
            fails.append(f"the panel opens on {dlg.current_field_state()!r}; the "
                         f"results view opens on 'failure'")
        if dlg.current_profile().get("field_state") != "failure":
            fails.append("the drawn profile does not report the at-failure state")
    finally:
        dlg.close()
    return fails


def test_field_state_profiles():
    """At-failure profiles are the snapshot's own: reinforcement axial force and
    pile shear both move when the state does, and the failure-band marks and the
    capacity envelope stay put."""
    fails = []
    from xslope import fem_details as fd

    slope_data, fem_data, solution = _solved(REINF_XLSX)
    failure = _failure_field(REINF_XLSX)
    banded = []
    for line_id in (1, 4, 6):
        conv = fd.reinforcement_profile(fem_data, solution, line_id, slope_data,
                                        field_state="converged",
                                        failure_solution=failure)
        fail = fd.reinforcement_profile(fem_data, solution, line_id, slope_data,
                                        field_state="failure",
                                        failure_solution=failure)
        if np.allclose(conv["T"], fail["T"]):
            fails.append(f"line {line_id}: the at-failure axial force equals the "
                         f"converged one — the switch reads the same field twice")
        if not np.allclose(conv["t_cap"], fail["t_cap"], equal_nan=True):
            fails.append(f"line {line_id}: the capacity moved with the field state")
        if not np.allclose(conv["env_T"], fail["env_T"]):
            fails.append(f"line {line_id}: the capacity envelope moved with the "
                         f"field state")
        # The band is the mechanism's, so it is the same reading in both states,
        # and the snapshot is what supplies it.
        if conv["band_lo"] != fail["band_lo"] or conv["band_hi"] != fail["band_hi"]:
            fails.append(f"line {line_id}: the failure band moved with the field "
                         f"state ({conv['band_lo']} -> {fail['band_lo']})")
        if fail["band_lo"] is not None:
            banded.append(line_id)
    # Which lines are banded is the mechanism's business (see
    # :func:`test_the_band_needs_the_mechanism`); what this check needs is that
    # the mechanism marked SOMETHING, or the comparison above compares two
    # absences.
    if not banded:
        fails.append("the snapshot marked no failure band on any line, so the "
                     "band's independence from the field state is untested")

    slope_data, fem_data, solution = _solved(PILES_XLSX)
    failure = _failure_field(PILES_XLSX)
    for pidx in (0, 1):
        conv = fd.pile_profile(fem_data, solution, pidx, slope_data,
                               field_state="converged", failure_solution=failure)
        fail = fd.pile_profile(fem_data, solution, pidx, slope_data,
                               field_state="failure", failure_solution=failure)
        if np.allclose(conv["shear"], fail["shear"]):
            fails.append(f"pile {pidx}: the at-failure shear equals the converged one")
        if np.allclose(conv["u_lateral"], fail["u_lateral"]):
            fails.append(f"pile {pidx}: the at-failure lateral displacement equals "
                         f"the converged one")
        if conv["band_depth"] != fail["band_depth"]:
            fails.append(f"pile {pidx}: the failure-band depth moved with the "
                         f"field state")
    return fails


def test_field_state_reload():
    """The snapshot is found the same way from disk: a solution exported with its
    at-failure twin and read back carries it, so a panel built on the reloaded
    solution alone has a live switch and the same at-failure profiles."""
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    _app()
    from pathlib import Path
    from studio.fem_details_dialog import FemDetailsDialog
    from xslope.fem import export_fem_solution, import_fem_solution
    from xslope import fem_details as fd

    slope_data, fem_data, solution = _solved(REINF_XLSX)
    failure = _failure_field(REINF_XLSX)
    with tempfile.TemporaryDirectory() as tmp:
        stem = Path(tmp) / "reloaded"
        export_fem_solution(fem_data, solution, stem, failure_solution=failure)
        back = import_fem_solution(fem_data, stem)

        if fd.failure_snapshot(back) is None:
            return ["a reloaded solution carries no at-failure snapshot"]
        a = fd.reinforcement_profile(fem_data, solution, 4, slope_data,
                                     field_state="failure", failure_solution=failure)
        b = fd.reinforcement_profile(fem_data, back, 4, slope_data,
                                     field_state="failure")
        _same_profile(a, b, "line 4 at failure", fails)

        dlg = FemDetailsDialog(fem_data, back, slope_data, model_path=REINF_XLSX)
        try:
            if not dlg.field_state.isEnabled():
                fails.append("the switch is dimmed for a reloaded solution that "
                             "carries an at-failure snapshot")
            if dlg.current_field_state() != "failure":
                fails.append("the reloaded panel does not open at failure")
            prof = dlg.current_profile()
            if prof.get("field_state") != "failure":
                fails.append("the reloaded panel draws the converged profile "
                             "while reading 'failure'")
            # And the switch actually moves the panel.
            dlg.field_state.setCurrentIndex(dlg.field_state.findData("converged"))
            moved = dlg.current_profile()
            if moved.get("field_state") != "converged":
                fails.append("switching to Last converged did not change the "
                             "drawn state")
            if np.allclose(prof["T"], moved["T"]):
                fails.append("switching the state did not change the drawn "
                             "axial force")
        finally:
            dlg.close()

    # A solution saved without the snapshot leaves the switch dimmed.
    with tempfile.TemporaryDirectory() as tmp:
        stem = Path(tmp) / "plain"
        export_fem_solution(fem_data, solution, stem)
        plain = import_fem_solution(fem_data, stem)
        dlg = FemDetailsDialog(fem_data, plain, slope_data, model_path=REINF_XLSX)
        try:
            if dlg.field_state.isEnabled():
                fails.append("the switch is live for a reloaded solution with no "
                             "at-failure snapshot")
            if dlg.current_field_state() != "converged":
                fails.append("a snapshot-less reloaded panel does not read the "
                             "converged field")
        finally:
            dlg.close()
    return fails


def test_field_state_export():
    """An at-failure export says so — in its default filename and in the CSV
    itself — and writes the profile that is on screen."""
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    app = _app()
    from studio.fem_details_dialog import FemDetailsDialog
    from xslope import fem_details as fd

    slope_data, fem_data, solution = _solved(REINF_XLSX)
    failure = _failure_field(REINF_XLSX)
    dlg = FemDetailsDialog(fem_data, solution, slope_data, model_path=REINF_XLSX,
                           failure_solution=failure)
    dlg.show()
    for _ in range(12):
        app.processEvents()
    try:
        stem = dlg.default_export_stem()
        if "at_failure" not in stem:
            fails.append(f"the at-failure export name {stem!r} does not name "
                         f"the state")
        with tempfile.TemporaryDirectory() as tmp:
            png, csv_path = dlg.export_current(path=os.path.join(tmp, "detail.png"))
            if not csv_path:
                return fails + ["the at-failure export wrote no CSV"]
            comments, rows = _read_export_csv(csv_path)
            if not any("field state" in c and "at failure" in c for c in comments):
                fails.append(f"the at-failure CSV does not record its state "
                             f"(comments: {comments})")
            cols, table = fd.profile_table(dlg.current_profile())
            if rows[0] != cols or rows[1:] != table:
                fails.append("the at-failure CSV is not the profile on screen")
            if not png or os.path.getsize(png) < 5000:
                fails.append("the at-failure export wrote no usable PNG")

            # And the converged export of the same member is a different file
            # with a different profile, not the same one written twice.
            dlg.field_state.setCurrentIndex(dlg.field_state.findData("converged"))
            for _ in range(12):
                app.processEvents()
            if dlg.default_export_stem() == stem:
                fails.append("both states export under the same default name")
            _png2, csv2 = dlg.export_current(path=os.path.join(tmp, "detail2.png"))
            comments2, rows2 = _read_export_csv(csv2)
            if not any("last converged" in c for c in comments2):
                fails.append(f"the converged CSV does not record its state "
                             f"(comments: {comments2})")
            if rows2[1:] == rows[1:]:
                fails.append("the two states exported identical numbers")
    finally:
        dlg.close()
    return fails


# --------------------------------------------------------------------------
# G. the failure band, and the labels around it
# --------------------------------------------------------------------------

def _label_boxes(ax, text):
    """The drawn box of every annotation on ``ax`` whose text starts with
    ``text``, in display pixels, with the panel's own frame."""
    renderer = ax.figure.canvas.get_renderer()
    boxes = [t.get_window_extent(renderer) for t in ax.texts
             if str(t.get_text()).strip().startswith(text)]
    return boxes, ax.get_window_extent(renderer)


def _drawn(profile, figsize=(9.5, 6.0)):
    """One detail figure, drawn on a real canvas so the label placement runs."""
    import matplotlib.figure as mplfig
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from xslope.plot_fem_details import plot_detail

    fig = mplfig.Figure(figsize=figsize)
    FigureCanvasAgg(fig)
    plot_detail(profile, fig=fig)
    fig.canvas.draw()
    return fig


def test_the_band_needs_the_mechanism():
    """A member gets a failure band only where the mechanism reaches it.

    The band's edges are the positions along a member where the mechanism field
    falls to BAND_FRACTION of its peak ALONG THAT MEMBER. Normalized to the
    member alone, every member has a peak and so every member has a band: on the
    reinforcement sample the lines the mechanism misses were given one-sample
    bands on their last element, drawn as a dashed rule hard against the end of
    the frame with a label running off the panel.

    The member's own peak is now measured against the mechanism's peak over the
    whole section, on the same fraction. Pinned on all six lines of the
    reinforcement sample, which has both cases, and on the drawn figures: a band
    that is one sample wide still draws inside the frame, and the label that
    names it is inside the panel in both cases.
    """
    fails = []
    from xslope import fem_details as fd

    slope_data, fem_data, solution = _solved(REINF_XLSX)
    failure = _failure_field(REINF_XLSX)
    peak = fd._mechanism_peak(fem_data, solution, failure)
    if not peak:
        return ["the sample's snapshot carries no mechanism field to band from"]

    banded, bare = [], []
    for line_id in range(1, 7):
        prof = fd.reinforcement_profile(fem_data, solution, line_id, slope_data,
                                        field_state="failure",
                                        failure_solution=failure)
        mech = prof.get("mechanism")
        if mech is None or not len(mech):
            fails.append(f"line {line_id}: the mechanism was not sampled at all")
            continue
        own = float(np.nanmax(mech))
        crossed = own >= fd.BAND_FRACTION * peak
        has_band = prof["band_lo"] is not None
        if has_band != crossed:
            fails.append(
                f"line {line_id} reaches {own / peak:.0%} of the mechanism's "
                f"peak and {'has' if has_band else 'has no'} failure band")
        (banded if has_band else bare).append((line_id, prof))

    if not banded:
        fails.append("no line is banded, so the drawn band is untested")
    if not bare:
        fails.append("every line is banded, so the gate is untested — the "
                     "sample no longer carries a member the mechanism misses")
    if fails:
        return fails

    # The drawn figures: a band inside the frame, and its label inside the
    # panel — including the narrowest band the sample has.
    narrowest = min(banded, key=lambda b: b[1]["band_hi"] - b[1]["band_lo"])
    for line_id, prof in (banded[0], narrowest):
        fig = _drawn(prof)
        ax = fig.axes[0]
        lo, hi = ax.get_xlim()
        if not (lo <= prof["band_lo"] <= hi and lo <= prof["band_hi"] <= hi):
            fails.append(f"line {line_id}: the band ({prof['band_lo']}, "
                         f"{prof['band_hi']}) falls outside the drawn frame "
                         f"({lo}, {hi})")
        boxes, frame = _label_boxes(ax, "failure band")
        if not boxes:
            fails.append(f"line {line_id}: the band is drawn and not named")
        for box in boxes:
            if not (frame.x0 <= box.x0 and box.x1 <= frame.x1
                    and frame.y0 <= box.y0 and box.y1 <= frame.y1):
                fails.append(f"line {line_id}: the band's label is cut by the "
                             f"panel edge (label {box.x0:.0f}..{box.x1:.0f}, "
                             f"panel {frame.x0:.0f}..{frame.x1:.0f})")

    # And the narrowest band there can be, at the place it was worst: one
    # sample wide, on the last element of the line. That is the band the sample
    # used to fabricate, and the label naming it ran off the right of the panel.
    line_id, prof = banded[0]
    end = float(np.asarray(prof["s"], dtype=float)[-1])
    ax = _drawn(dict(prof, band_lo=end, band_hi=end, band_peak=end)).axes[0]
    boxes, frame = _label_boxes(ax, "failure band")
    if not boxes:
        fails.append(f"line {line_id}: a band on the last element is not named")
    for box in boxes:
        if not (frame.x0 <= box.x0 and box.x1 <= frame.x1
                and frame.y0 <= box.y0 and box.y1 <= frame.y1):
            fails.append(f"line {line_id}: a band on the last element puts its "
                         f"label outside the panel ({box.x0:.0f}..{box.x1:.0f} "
                         f"against {frame.x0:.0f}..{frame.x1:.0f})")

    # A line the mechanism misses draws no band and names none.
    line_id, prof = bare[0]
    ax = _drawn(prof).axes[0]
    if _label_boxes(ax, "failure band")[0]:
        fails.append(f"line {line_id}: no band was measured and one is named")

    # Mutation: the gate removed — the member's own peak judged against itself,
    # which is what fabricated the bands this check exists to refuse.
    real = fd._band_span
    fd._band_span = lambda positions, mech, global_peak=None: real(positions, mech)
    try:
        still = 0
        for line_id, _prof in bare:
            prof = fd.reinforcement_profile(fem_data, solution, line_id,
                                            slope_data, field_state="failure",
                                            failure_solution=failure)
            still += prof["band_lo"] is not None
        if still != len(bare):
            fails.append(f"without the gate only {still} of {len(bare)} unbanded "
                         f"lines gain a band; the gate is not what removed them")
    finally:
        fd._band_span = real
    return fails


def test_the_reaction_panel_says_where_its_limit_is():
    """The pile's soil-reaction note never cites a line the panel does not show.

    The panel is scaled to the mobilized profile, because the Ito & Matsui
    limiting resistance for a pile well inside its working range would flatten
    that profile onto the axis. So the envelope sometimes lands on the panel and
    sometimes does not — and the note stating the peak mobilization was printed
    either way: "peak 20% of limit", over a panel with no limit drawn on it and
    no legend to name one.

    Both cases, off the drawn panel: with the envelope in frame it is drawn and
    named in a legend; with it off scale the note says so and says where it is.
    """
    fails = []
    from xslope import fem_details as fd
    from xslope.plot_fem_details import C_LIMIT

    slope_data, fem_data, solution = _solved(PILES_XLSX)
    profile = next(
        (p for p in (fd.pile_profile(fem_data, solution, m["index"], slope_data)
                     for m in fd.list_lines(fem_data, solution, slope_data)
                     if m["kind"] == "pile")
         if p.get("limit_p") is not None and p.get("reaction_ratio") is not None),
        None)
    if profile is None:
        fails.append("no pile carries a limiting-resistance envelope and a "
                     "mobilization to measure against it")
        return fails

    def panel(prof):
        """(legend?, limit drawn?, the note) for the soil-reaction panel."""
        fig = _drawn(prof, figsize=(11.0, 6.0))
        ax = fig.axes[3]
        drawn = any(line.get_color() == C_LIMIT for line in ax.lines)
        said = " ".join(t.get_text() for t in ax.texts)
        return ax.get_legend() is not None, drawn, said

    # In frame: drawn, named in the legend, and the note claims nothing about a
    # scale the envelope is on.
    near = dict(profile, limit_p=np.asarray(profile["limit_p"], float) * 0.05,
                reaction_ratio=min(1.0, profile["reaction_ratio"] * 20.0))
    legend, drawn, said = panel(near)
    if not drawn:
        fails.append("an envelope inside the panel's own scale is not drawn")
    if not legend:
        fails.append("the envelope is drawn beside the mobilized profile and "
                     "neither is named")
    if "off scale" in said:
        fails.append(f"the envelope is on the panel and the note reads {said!r}")

    # Off scale: not drawn, no legend of one series, and the note says where the
    # envelope went rather than citing a line that is not there.
    lo = float(np.nanmin(np.asarray(profile["limit_p"], float)))
    far = dict(profile, limit_p=np.asarray(profile["limit_p"], float) * 8.0,
               reaction_ratio=profile["reaction_ratio"] / 8.0)
    legend, drawn, said = panel(far)
    if drawn:
        fails.append("an envelope past the panel's scale is drawn on it anyway")
    if legend:
        fails.append("a panel carrying one series carries a legend")
    if "of limit" not in said:
        fails.append(f"the peak mobilization is not stated: {said!r}")
    if "off scale" not in said:
        fails.append(f"the note cites a limit the panel does not show and does "
                     f"not say where it is: {said!r}")
    if f"{lo * 8.0:,.0f}" not in said:
        fails.append(f"the note does not say how far off the limit is "
                     f"({lo * 8.0:,.0f}): {said!r}")
    return fails


def test_the_band_is_named_for_the_field_it_was_read_from():
    """The band across a member is called what it is on the run that drew it.

    The band is read from the at-failure snapshot whenever there is one, so on a
    strength reduction run it marks the mechanism and "failure band" is what it
    is. A run that converged under gravity captured no snapshot and reached no
    failure; its band is where the computed shear strain concentrates in a
    section that is standing, and calling that a failure band states a collapse
    the analysis never found.

    Both members, both runs, off the drawn figures.
    """
    fails = []
    from xslope import fem_details as fd
    from xslope.plot_fem_details import band_label

    for path, kind, read in ((REINF_XLSX, "reinforcement",
                              fd.reinforcement_profile),
                             (PILES_XLSX, "pile", fd.pile_profile)):
        slope_data, fem_data, solution = _solved(path)
        failure = _failure_field(path)
        ids = [m["index"] for m in fd.list_lines(fem_data, solution, slope_data)
               if m["kind"] == kind]
        drawn = {}
        for state, snapshot, want in (("failure", failure, "failure band"),
                                      ("converged", None, "shear strain band")):
            banded = None
            for index in ids:
                prof = read(fem_data, solution, index, slope_data,
                            field_state="converged", failure_solution=snapshot)
                if fd.band_state(solution, snapshot) != state:
                    fails.append(f"{kind}: a run {'with' if snapshot else 'without'}"
                                 f" a snapshot reads its band from "
                                 f"{fd.band_state(solution, snapshot)!r}")
                    break
                if band_label(prof) != want:
                    fails.append(f"{kind} {prof['label']}: the {state} run's "
                                 f"band is called {band_label(prof)!r}")
                has_band = (prof.get("band_depth") is not None
                            if kind == "pile" else prof.get("band_lo") is not None)
                if has_band and banded is None:
                    banded = prof
            drawn[state] = banded
            if banded is None:
                fails.append(f"{kind}: the {state} run bands no member, so the "
                             f"name on the figure is untested")
                continue
            said = " ".join(t.get_text()
                            for ax in _drawn(banded).axes for t in ax.texts)
            if want not in said:
                fails.append(f"{kind} {banded['label']}: the {state} run's "
                             f"figure does not say {want!r}: {said!r}")
            other = ("shear strain band" if want == "failure band"
                     else "failure band")
            if other in said:
                fails.append(f"{kind} {banded['label']}: the {state} run's "
                             f"figure says {other!r}")

    # Mutation: the name fixed at "failure band" whatever the run, which is what
    # the figures said before. The converged run has to catch it.
    import xslope.plot_fem_details as ppd
    real = ppd.band_label
    ppd.band_label = lambda profile: "failure band"
    try:
        slope_data, fem_data, solution = _solved(REINF_XLSX)
        prof = next((p for p in (fd.reinforcement_profile(
            fem_data, solution, i, slope_data) for i in range(1, 7))
            if p.get("band_lo") is not None), None)
        if prof is None:
            fails.append("the gravity run bands no line, so the mutation has "
                         "nothing to act on")
        else:
            said = " ".join(t.get_text() for t in _drawn(prof).axes[0].texts)
            if "shear strain band" in said:
                fails.append("the name was pinned to 'failure band' and the "
                             "figure still said 'shear strain band'")
    finally:
        ppd.band_label = real
    return fails


def _broken_tie_set(fd, fem_data, bundle, slope_data, tied_lines, fails):
    """One ``(state, line_id, profile, gaps)`` whose tie set has a hole in it.

    Takes the first line that stands at its greatest utilization over a stretch,
    copies the field it was read from, and drops the force in one bar element
    inside that stretch to half the capacity it was holding. The profile is then
    re-read through :func:`xslope.fem_details.reinforcement_profile`, so the tie
    set, the span, the gap positions and the figure's label are all computed by
    the shipping code — only the force it reads is arranged.

    Whether the shipped run produces such a hole of its own depends on the
    mechanism the model develops, which is not what this is testing.
    """
    state, line_id, prof = tied_lines[0]
    tied = np.asarray(prof["peak_indices"], dtype=int)
    element_ids = np.asarray(prof["element_ids"], dtype=int)
    hole = int(tied[len(tied) // 2])          # an interior sample of the stretch
    field = (bundle.get("failure_solution") if state == "failure"
             else bundle["solution"])
    knocked = dict(field)
    forces = np.array(field["forces_1d"], dtype=float)
    forces[element_ids[hole]] = 0.5 * forces[element_ids[hole]]
    knocked["forces_1d"] = forces
    solution = bundle["solution"]
    failure = knocked if state == "failure" else bundle.get("failure_solution")
    if state != "failure":
        solution = knocked
    out = fd.reinforcement_profile(fem_data, solution, line_id, slope_data,
                                   field_state=state, failure_solution=failure)
    gaps = [round(float(v), 6) for v in out.get("peak_gap_s", [])]
    if not gaps:
        fails.append(f"{state} line {line_id}: a bar element inside the stretch "
                     f"was dropped to half its force and the profile still "
                     f"reports an unbroken span")
    return [(state, line_id, out, gaps)]


def test_the_peak_utilization_is_tie_aware():
    """A member at its greatest utilization over a stretch reports the stretch.

    The capacity envelope is flat along the middle of a reinforcement line and
    the axial force is capped by it, so the utilization sits at one from
    wherever the bar first reaches capacity to wherever it stops. ``argmax``
    returns the first sample of that run, and the report printed it as THE point
    of greatest utilization — a position the bar does not distinguish, marked on
    the figure with a ring that says look here.

    Read off the sample's own solved run: every sample within
    :data:`UTIL_TIE_TOL` of the greatest is in the span, the span's ends are the
    first and last of them, and the figure rings all of them rather than one.

    A tie set can have a HOLE in it — a sample between two at-capacity ones that
    stands below them — and the two ends alone describe an unbroken run instead.
    The samples inside the span that do not stand with the rest are reported as
    ``peak_gap_s``, and the figure's label names them. That case is put to the
    same code path on a field with one interior bar force knocked down
    (:func:`_broken_tie_set`), since whether the shipped run develops a hole of
    its own is a property of the mechanism rather than of the reporting.
    """
    fails = []
    from xslope import fem_details as fd
    from xslope.fileio import load_slope_data
    from xslope.report import solutions_from_sidecars

    with contextlib.redirect_stdout(io.StringIO()):
        slope_data = load_slope_data(REINF_XLSX)
        solutions = solutions_from_sidecars(REINF_XLSX, slope_data, None)
    bundle = solutions.get("fem")
    if not bundle:
        return ["the reinforcement sample ships no solved run"]
    fem_data = bundle["fem_data"]

    tied_lines, single_lines, broken_lines = [], [], []
    for state in ("failure", "converged"):
        for line_id in range(1, 7):
            prof = fd.reinforcement_profile(
                fem_data, bundle["solution"], line_id, slope_data,
                field_state=state, failure_solution=bundle.get("failure_solution"))
            util = np.asarray(prof["utilization"], dtype=float)
            if not np.any(np.isfinite(util)):
                continue
            want = np.where(util >= np.nanmax(util) - fd.UTIL_TIE_TOL)[0]
            got = np.asarray(prof["peak_indices"], dtype=int)
            if not np.array_equal(got, want):
                fails.append(f"{state} line {line_id}: {list(want)} stand at the "
                             f"greatest utilization and the profile reports "
                             f"{list(got)}")
                continue
            s = np.asarray(prof["s"], dtype=float)
            if len(want) > 1:
                tied_lines.append((state, line_id, prof))
                if prof["peak_span"] != (float(s[want[0]]), float(s[want[-1]])):
                    fails.append(f"{state} line {line_id}: the span "
                                 f"{prof['peak_span']} is not the first and last "
                                 f"of the tied samples")
                T = np.asarray(prof["T"], dtype=float)[want]
                if prof["peak_T_span"] != (float(T.min()), float(T.max())):
                    fails.append(f"{state} line {line_id}: the force range "
                                 f"{prof['peak_T_span']} is not the range over "
                                 f"the span")
                # What the span leaves out: the samples between its ends that
                # are not tied. Empty for an unbroken run.
                inside = set(range(int(want[0]), int(want[-1]) + 1))
                holes = sorted(inside - set(int(i) for i in want))
                got_gaps = [round(float(v), 6)
                            for v in prof.get("peak_gap_s", [])]
                want_gaps = [round(float(s[i]), 6) for i in holes]
                if got_gaps != want_gaps:
                    fails.append(f"{state} line {line_id}: the samples inside "
                                 f"the span that stand below it are at "
                                 f"{want_gaps} and the profile reports "
                                 f"{got_gaps}")
                elif holes:
                    broken_lines.append((state, line_id, prof, want_gaps))
            else:
                single_lines.append((state, line_id, prof))
                if prof["peak_span"] is not None:
                    fails.append(f"{state} line {line_id}: one sample stands at "
                                 f"the greatest utilization and a span "
                                 f"{prof['peak_span']} is reported")

    if not tied_lines:
        fails.append("no line reaches its greatest utilization at more than one "
                     "point, so the span is untested")
    if not single_lines:
        fails.append("every line reports a span, so the single-point case is "
                     "untested")
    if fails:
        return fails

    # A tie set with a HOLE in it. Whether the shipped run leaves an interior
    # sample below the rest is a property of the mechanism and not of the
    # reporting, so the broken case is put through the same code path on a field
    # with one interior bar force knocked off capacity: everything from
    # _peak_utilization outward is the shipping code, and only the force it
    # reads is arranged.
    if not broken_lines:
        broken_lines = _broken_tie_set(fd, fem_data, bundle, slope_data,
                                       tied_lines, fails)
    if fails:
        return fails

    # The figure's label says what its thickened runs draw: a span with a hole
    # in it excepts the hole by position, and an unbroken one excepts nothing.
    for state, line_id, prof, gaps in (broken_lines[0],):
        said = " ".join(t.get_text() for t in _drawn(prof).axes[0].texts)
        if "except" not in said:
            fails.append(f"{state} line {line_id}: the label reads {said!r} and "
                         f"the line stands below capacity at {gaps} inside its "
                         f"span")
        for gap in gaps:
            if f"{gap:,.2f}" not in said:
                fails.append(f"{state} line {line_id}: {gap:,.2f} is excepted "
                             f"from the span and the label reads {said!r}")
    unbroken = next((t for t in tied_lines
                     if not len(t[2].get("peak_gap_s", []))), None)
    if unbroken is not None:
        said = " ".join(t.get_text() for t in _drawn(unbroken[2]).axes[0].texts)
        if "except" in said:
            fails.append(f"{unbroken[0]} line {unbroken[1]}: the span is "
                         f"unbroken and the label reads {said!r}")

    # The figure: as many rings as there are samples at the maximum.
    for state, line_id, prof in (tied_lines[0], single_lines[0]):
        ax = _drawn(prof).axes[0]
        rings = 0
        for line in ax.lines:
            if line.get_marker() == "o" and line.get_markerfacecolor() == "none":
                rings += len(line.get_xdata())
        want = len(prof["peak_indices"])
        if rings != want:
            fails.append(f"{state} line {line_id}: {want} sample(s) stand at the "
                         f"greatest utilization and the figure rings {rings}")

    # Mutation: the argmax restored as the whole answer, which is what the
    # profile used to report. Every span has to go, and the figure with it.
    real = fd._peak_utilization
    fd._peak_utilization = lambda util: (
        int(np.nanargmax(util)), np.array([int(np.nanargmax(util))]))
    try:
        state, line_id, _prof = tied_lines[0]
        prof = fd.reinforcement_profile(
            fem_data, bundle["solution"], line_id, slope_data,
            field_state=state, failure_solution=bundle.get("failure_solution"))
        if prof["peak_span"] is not None:
            fails.append("the argmax alone still produces a span; the span is "
                         "not read from the tie set")
        rings = sum(len(l.get_xdata()) for l in _drawn(prof).axes[0].lines
                    if l.get_marker() == "o"
                    and l.get_markerfacecolor() == "none")
        if rings != 1:
            fails.append(f"the argmax alone still rings {rings} samples")
    finally:
        fd._peak_utilization = real
    return fails


def _highlighted_ends(fig):
    """The two endpoints of the member the drawn map picked out, or None — read
    off the artist, in the colour the map highlights with."""
    from xslope.plot_fem_details import C_PEAK

    for ax in fig.axes:
        for line in ax.lines:
            if line.get_color() == C_PEAK and len(line.get_xdata()) == 2:
                xs, ys = line.get_xdata(), line.get_ydata()
                return frozenset((round(float(x), 6), round(float(y), 6))
                                 for x, y in zip(xs, ys))
    return None


def test_the_inset_follows_the_selection():
    """The details panel carries a map of the section with the SELECTED member
    picked out, and it follows the selection.

    A profile is read along a member and says nothing about where on the slope
    that member is; the list names members a user has to place before the
    profile means anything. The map is the report's own locator drawing, so the
    panel and the page cannot disagree about where a member is.
    """
    fails = []
    import matplotlib
    matplotlib.use("Agg")
    app = _app()
    from PySide6.QtCore import Qt
    from studio.fem_details_dialog import FemDetailsDialog
    from xslope.fem_details import member_lines

    for path, want_kind in ((REINF_XLSX, "reinforcement"), (PILES_XLSX, "pile")):
        slope_data, fem_data, solution = _solved(path)
        dlg = FemDetailsDialog(fem_data, solution, slope_data, model_path=path)
        dlg.show()
        try:
            for _ in range(12):
                app.processEvents()
            if not hasattr(dlg, "map_canvas"):
                fails.append(f"{want_kind}: the panel carries no map")
                continue

            geometry = {(want_kind, m["index"]): frozenset(
                (round(x, 6), round(y, 6)) for x, y in m["ends"])
                for m in member_lines(fem_data, slope_data, want_kind)}
            rows = [r for r in range(dlg.list.count())
                    if dlg.list.item(r).data(Qt.UserRole) is not None]
            if len(rows) < 2:
                fails.append(f"{want_kind}: fewer than two members to select "
                             f"between, so following the selection is untested")
                continue

            seen = set()
            for row in rows:
                dlg.list.setCurrentRow(row)
                for _ in range(6):
                    app.processEvents()
                entry = dlg.list.item(row).data(Qt.UserRole)
                key = (entry["kind"], entry["index"])
                if dlg.mapped_member() != key:
                    fails.append(f"{want_kind}: selecting {entry['label']} maps "
                                 f"{dlg.mapped_member()}")
                    continue
                dlg.map_canvas.render_now()
                ends = _highlighted_ends(dlg.map_canvas.figure)
                if ends is None:
                    fails.append(f"{want_kind}: the map picks out nothing with "
                                 f"{entry['label']} selected")
                elif ends != geometry.get(key):
                    fails.append(f"{want_kind}: {entry['label']} is picked out "
                                 f"at {sorted(ends)}, and it runs "
                                 f"{sorted(geometry.get(key) or [])}")
                seen.add(ends)
            if len(seen) < len(rows):
                fails.append(f"{want_kind}: {len(rows)} members were selected "
                             f"and the map drew {len(seen)} different "
                             f"highlights")
        finally:
            dlg.close()
    return fails


CHECKS = [
    ("the failure band needs the mechanism", test_the_band_needs_the_mechanism),
    ("the reaction panel says where its limit is",
     test_the_reaction_panel_says_where_its_limit_is),
    ("the band is named for the field it was read from",
     test_the_band_is_named_for_the_field_it_was_read_from),
    ("the inset follows the selection", test_the_inset_follows_the_selection),
    ("the peak utilization is tie-aware",
     test_the_peak_utilization_is_tie_aware),
    ("the toolbar button and its gate", test_gate),
    ("the gate check would catch a mutation", test_gate_mutation),
    ("the member list and its badges", test_list),
    ("the panel's verdict is the summary's",
     test_the_verdict_is_the_summary_s),
    ("the envelope is the solver's capacity", test_envelope),
    ("the envelope is the shared function", test_envelope_mutation),
    ("profiles survive save + reload", test_reload),
    ("the dialog is the same on a reloaded solution", test_dialog_reload_identical),
    ("export writes the PNG and the profile CSV", test_export),
    ("the field state switch and its gate", test_field_state_control),
    ("at-failure profiles are the snapshot's", test_field_state_profiles),
    ("the field state survives save + reload", test_field_state_reload),
    ("the export records the field state", test_field_state_export),
]

# Checks that need the Studio layer; skipped when PySide6 is absent.
_STUDIO_ONLY = {test_gate, test_gate_mutation, test_list,
                test_dialog_reload_identical, test_export,
                test_field_state_control, test_field_state_reload,
                test_field_state_export, test_the_inset_follows_the_selection}


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
