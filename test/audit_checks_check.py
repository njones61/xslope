"""Standing checks on the loader and search fixes from the gray-cell audit
(Norm, 2026-11): inputs the sheet accepts that the loader or the non-circular
search used to read as something else, without a word.

  A. SEEP BC BLOCKS. A block whose Value cell is blank used to END the scan, so
     that block and every block after it were dropped. Every block is read now:
     a half-filled one is kept with the value None (and the writer saves it back
     as a blank), and preflight's seep.bc_value_missing names it.
  B. UNSAT WORDS. Any word outside lf / vg / gard was read as lf, which left
     mat.unsat_vocabulary nothing to see. The word is kept; the three match in
     any case; a blank is lf.
  C. NON-CIRC MOVEMENT. Free / Horiz / Fixed match in any case (a lowercase
     "free" on an interior point used to freeze it as Fixed, and a lowercase
     "fixed" on an end point let it slide). A blank keeps its meaning: an end
     point slides along the ground surface, an interior point stays where it is.
     Any other word is refused, by preflight and by the search itself, before a
     single surface is solved.
  D. CIRCLES OPTION. Depth / Radius / Intercept match in any case, as the
     sheet's formatting and dropdown do ("depth" used to stop the file opening).
     A blank Option opens too, with the circle's size undefined, and preflight's
     surface.circle_size_missing refuses the runs that read the circles sheet.
  E. RAPID DRAWDOWN. Only an mc row's d / psi count toward the model having a
     low-permeability material: a stale pair on a cp row (whose cells the
     template grays) is read by nothing.
  F. RAPID DRAWDOWN RELIABILITY (Norm, 531a). A rapid drawdown Taylor-series run
     perturbs an mc material's d and psi by s(d) and s(psi), as it does c and phi;
     a run without drawdown leaves them out. Checked against a hand combination of
     separate deterministic rapid drawdown runs at MLV +/- sigma, on the fixed
     surface (search=False), which also checks that path now runs the drawdown.

Every check here fails on the code before the fix. No solve is run: each one
reads a workbook, a preflight report, or a refusal raised before the first trial
surface.
"""
import contextlib
import io
import os
import shutil
import sys
import tempfile
import warnings

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

#: A seepage model with two head blocks on its seep bc sheet (E/F and H/I).
SEEP_MODEL = os.path.join(_REPO, "docs", "inputs", "seep", "xslope_earth_dam1.xlsx")
#: A limit equilibrium model with circles and a one-material section.
LEM_MODEL = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_simple1.xlsx")
#: Rocscience VP96: one mc embankment with d = 1379 and psi = 18.2 under a
#: drawdown from piezometric line 1 to line 2, and one circle on its sheet.
RAPID_MODEL = os.path.join(_REPO, "docs", "verification", "files", "rocscience",
                           "vp096.xlsx")
#: The standard deviations the reliability check puts on VP96's embankment.
RAPID_SIGMAS = {"sigma_phi": 2.0, "sigma_d": 200.0, "sigma_psi": 2.0}
#: A non-circular model (Rocscience VP47): four points, all Movement = Free.
NONCIRC_MODEL = os.path.join(_REPO, "docs", "verification", "files", "rocscience",
                             "vp047.xlsx")


def _quiet(fn, *a, **kw):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        return fn(*a, **kw)


def _load(path):
    from xslope.fileio import load_slope_data
    return _quiet(load_slope_data, path)


@contextlib.contextmanager
def _edited(src, edits):
    """A temporary copy of ``src`` with ``edits`` = {sheet: {cell: value}} written
    into it through openpyxl, so the LOADER is what is tested, not the writer."""
    import openpyxl
    tmpdir = tempfile.mkdtemp(prefix="audit_checks_")
    dst = os.path.join(tmpdir, os.path.basename(src))
    try:
        shutil.copy(src, dst)
        wb = openpyxl.load_workbook(dst)
        for sheet, cells in edits.items():
            ws = wb[sheet]
            for ref, value in cells.items():
                ws[ref] = value
        wb.save(dst)
        wb.close()
        yield dst
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# --------------------------------------------------------------------------- #
# A. seep bc blocks
# --------------------------------------------------------------------------- #
def check_seep_bc_blocks():
    """Block 1 keeps its head, block 2 (points, blank Value) is kept with None,
    and block 3 after it is read; preflight names block 2; the writer saves the
    blank as a blank and the block reopens."""
    fails = []
    # E/F holds block 1 (head 18) and H/I block 2 (head 2) in the shipped model.
    # Blank block 2's value and draw a block 3 in K/L.
    edits = {"seep bc": {"I3": None, "K3": "head", "L3": 7.5,
                         "K5": 0.0, "L5": 7.5, "K6": 10.0, "L6": 7.5}}
    with _edited(SEEP_MODEL, edits) as path:
        try:
            sd = _load(path)
        except Exception as exc:
            return [f"seep bc: a block with a blank Value stopped the file "
                    f"opening ({exc!r})"]
    heads = sd["seepage_bc"]["specified_heads"]
    got = [h.get("head") for h in heads]
    if len(heads) != 3:
        fails.append(f"seep bc: expected 3 head blocks (18, blank, 7.5), read "
                     f"{len(heads)}: {got}")
    else:
        if got[0] != 18.0 or got[1] is not None or got[2] != 7.5:
            fails.append(f"seep bc: head values read {got}, expected "
                         f"[18.0, None, 7.5]")
    if len(heads) >= 2:
        from xslope.preflight import preflight
        rep = _quiet(preflight, sd, "seep", {})
        hits = [f for f in rep.errors if f.rule_id == "seep.bc_value_missing"]
        if not any("Specified head #2" in f.message for f in hits):
            fails.append("seep bc: preflight did not refuse the seepage run over "
                         "the block with a blank head")
        # ... and the writer keeps the blank a blank, so the block reopens.
        from xslope.fileio import save_slope_data_to_xlsx
        tmpdir = tempfile.mkdtemp(prefix="audit_checks_")
        try:
            out = os.path.join(tmpdir, "roundtrip.xlsx")
            try:
                _quiet(save_slope_data_to_xlsx, sd, out)
                back = _load(out)["seepage_bc"]["specified_heads"]
                if [h.get("head") for h in back] != got:
                    fails.append(f"seep bc: save and reopen read "
                                 f"{[h.get('head') for h in back]}, expected {got}")
            except Exception as exc:
                fails.append(f"seep bc: saving a block with a blank head failed "
                             f"({exc!r})")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


# --------------------------------------------------------------------------- #
# B. unsat words
# --------------------------------------------------------------------------- #
def _mat_cell(path, header, data_row=0):
    """The mat-sheet cell (e.g. 'AJ11') under the header ``header`` in the
    ``data_row``-th material row."""
    import openpyxl
    from openpyxl.utils import get_column_letter
    wb = openpyxl.load_workbook(path, read_only=True)
    try:
        ws = wb["mat"]
        for r in range(1, 21):
            row = [c.value for c in next(ws.iter_rows(min_row=r, max_row=r))]
            norm = [str(v or "").strip().replace("_", "").lower() for v in row]
            if "name" in norm and header in norm:
                return f"{get_column_letter(norm.index(header) + 1)}{r + 1 + data_row}"
    finally:
        wb.close()
    raise LookupError(f"no '{header}' header on the mat sheet of {path}")


def check_unsat_words():
    fails = []
    cell = _mat_cell(SEEP_MODEL, "unsat")
    for typed, want in (("vgn", "vgn"), ("VG", "vg"), (" Gard ", "gard"),
                        (None, "lf")):
        with _edited(SEEP_MODEL, {"mat": {cell: typed}}) as path:
            got = _load(path)["materials"][0]["unsat"]
        if got != want:
            fails.append(f"unsat: a cell holding {typed!r} loaded as {got!r}, "
                         f"expected {want!r}")
    with _edited(SEEP_MODEL, {"mat": {cell: "vgn"}}) as path:
        sd = _load(path)
    from xslope.preflight import preflight
    rep = _quiet(preflight, sd, "seep", {})
    if not any(f.rule_id == "mat.unsat_vocabulary" for f in rep.errors):
        fails.append("unsat: a saved 'vgn' reached a seepage run with no "
                     "mat.unsat_vocabulary error")
    rep = _quiet(preflight, sd, "fem", {})
    if any(f.rule_id == "mat.unsat_vocabulary" for f in rep.findings):
        fails.append("unsat: mat.unsat_vocabulary fired on a finite element run, "
                     "which never reads unsat")
    return fails


# --------------------------------------------------------------------------- #
# C. non-circ Movement
# --------------------------------------------------------------------------- #
def check_movement_loader():
    """The loader matches the three words in any case, reads a blank as None and
    keeps any other entry as typed."""
    fails = []
    edits = {"non-circ": {"C3": "fixed", "C4": "free", "C5": None, "C6": "sideways"}}
    with _edited(NONCIRC_MODEL, edits) as path:
        got = [p["Movement"] for p in _load(path)["non_circ"]]
    want = ["Fixed", "Free", None, "sideways"]
    if got != want:
        fails.append(f"Movement: the loader read {got}, expected {want}")
    return fails


def check_movement_meaning():
    """What the search does with each entry, including a blank at each end."""
    fails = []
    try:
        from xslope.search import noncirc_point_movements
    except ImportError as exc:
        return [f"Movement: no noncirc_point_movements in xslope.search ({exc})"]
    pts = lambda *words: [{"X": float(i), "Y": 0.0, "Movement": w}
                          for i, w in enumerate(words)]
    cases = [
        # blank: an end point slides along the ground, an interior point stays
        ((None, None, float("nan"), ""), ["Free", "Fixed", "Fixed", "Free"]),
        # any case
        (("fixed", "free", "HORIZ", "Fixed"), ["Fixed", "Free", "Horiz", "Fixed"]),
    ]
    for words, want in cases:
        got = noncirc_point_movements(pts(*words))
        if got != want:
            fails.append(f"Movement: {list(words)} moves as {got}, expected {want}")
    try:
        noncirc_point_movements(pts("Free", "sideways", "Free"))
        fails.append("Movement: 'sideways' was accepted by the search")
    except ValueError as exc:
        if "sideways" not in str(exc):
            fails.append(f"Movement: the refusal does not name the entry: {exc}")
    return fails


def check_movement_search_refuses():
    """noncircular_search refuses an unknown word at its gate, before it slices or
    solves a single trial surface."""
    from xslope import search
    from xslope.preflight import PreflightError
    sd = _load(NONCIRC_MODEL)
    sd["non_circ"][1]["Movement"] = "sideways"
    calls = []
    real = search.generate_slices

    def _tripwire(*a, **kw):
        calls.append(1)
        raise RuntimeError("a trial surface was sliced")
    search.generate_slices = _tripwire
    try:
        try:
            _quiet(search.noncircular_search, sd, "spencer", max_iter=1)
            return ["Movement: noncircular_search ran with Movement 'sideways'"]
        except PreflightError as exc:
            if "surface.movement_vocabulary" not in {f.rule_id for f in
                                                     exc.report.errors}:
                return [f"Movement: the search was refused, but not for the "
                        f"Movement: {exc}"]
        except RuntimeError:
            return ["Movement: noncircular_search sliced a trial surface with "
                    "Movement 'sideways' instead of refusing it"]
    finally:
        search.generate_slices = real
    return [] if not calls else ["Movement: a trial surface was sliced"]


# --------------------------------------------------------------------------- #
# D. circles Option
# --------------------------------------------------------------------------- #
def check_circle_option_case():
    fails = []
    sd = _load(LEM_MODEL)
    c0 = sd["circles"][0]
    xo, yo, depth = float(c0["Xo"]), float(c0["Yo"]), float(c0["Depth"])
    r = yo - depth
    for typed, cells in (("depth", {"E3": depth}), (" RADIUS ", {"H3": r}),
                         ("intercept", {"F3": xo, "G3": yo - r})):
        edits = {"circles": dict({"B3": xo, "C3": yo, "D3": typed}, **cells)}
        with _edited(LEM_MODEL, edits) as path:
            try:
                got = _load(path)["circles"][0]
            except Exception as exc:
                fails.append(f"circles: Option {typed!r} stopped the file "
                             f"opening ({exc})")
                continue
        if abs(float(got["R"]) - r) > 1e-9 or abs(float(got["Depth"]) - depth) > 1e-9:
            fails.append(f"circles: Option {typed!r} read R = {got['R']}, Depth = "
                         f"{got['Depth']}; expected {r}, {depth}")
    # A blank Option: the file opens, and a run that reads the circles is refused.
    with _edited(LEM_MODEL, {"circles": {"D3": None}}) as path:
        try:
            sd = _load(path)
        except Exception as exc:
            return fails + [f"circles: a blank Option stopped the file opening "
                            f"({exc})"]
    from xslope.preflight import preflight
    rep = _quiet(preflight, sd, "lem", {"surface": "circular"})
    if not any(f.rule_id == "surface.circle_size_missing"
               and "Circle 1" in f.message for f in rep.errors):
        fails.append("circles: a blank Option reached a circular run with no "
                     "surface.circle_size_missing error")
    return fails


# --------------------------------------------------------------------------- #
# E. rapid drawdown: d / psi on mc rows only
# --------------------------------------------------------------------------- #
def check_rapid_d_psi_rows():
    from xslope.advanced import validate_rapid_drawdown
    fails = []
    stale = {"materials": [
        {"name": "clay", "option": "cp", "u": "piezo", "d": 500.0, "psi": 15.0},
        {"name": "shell", "option": "mc", "u": "piezo", "d": 0.0, "psi": 0.0}],
        "water_loads": "auto", "piezo_line2": [(0.0, 0.0), (1.0, 0.0)]}
    try:
        _quiet(validate_rapid_drawdown, stale)
        fails.append("rapid: a d/psi pair on a cp row (read by nothing) passed as "
                     "the model's low-permeability material")
    except ValueError:
        pass
    good = {"materials": [dict(stale["materials"][1], d=500.0, psi=15.0)],
            "water_loads": "auto", "piezo_line2": [(0.0, 0.0), (1.0, 0.0)]}
    try:
        _quiet(validate_rapid_drawdown, good)
    except ValueError as exc:
        fails.append(f"rapid: an mc row with d and psi was refused ({exc})")
    return fails


# --------------------------------------------------------------------------- #
# F. rapid drawdown reliability varies d and psi
# --------------------------------------------------------------------------- #
def _rapid_fs(sd, circle):
    """The deterministic rapid drawdown factor of safety on ``circle``: the
    slicer and the three-stage procedure, nothing from the reliability module."""
    from xslope.slice import generate_slices
    from xslope.advanced import rapid_drawdown
    ok, res = generate_slices(sd, circle=dict(circle), num_slices=40,
                              check_inputs=False, rapid=True)
    if not ok:
        raise RuntimeError(f"slicing failed: {res}")
    ok2, r = _quiet(rapid_drawdown, res[0], "spencer", debug_level=0)
    if not ok2:
        raise RuntimeError(f"rapid drawdown failed: {r}")
    return float(r["FS"])


def rapid_reliability_by_hand(sd):
    """sigma_F and beta for VP96's embankment, combined here from separate
    deterministic rapid drawdown runs: F at the most likely values, and F+ / F-
    with one parameter at a time moved by +/- its sigma (Taylor series, Duncan
    2000), beta in the lognormal form the engine reports."""
    import copy
    import math
    circle = sd["circles"][0]
    f_mlv = _rapid_fs(sd, circle)
    m0 = sd["materials"][0]
    rows = {}
    for param, key in (("phi", "sigma_phi"), ("d", "sigma_d"), ("psi", "sigma_psi")):
        out = []
        for sign in (+1, -1):
            sd2 = copy.deepcopy(sd)
            sd2["materials"][0][param] = m0[param] + sign * m0[key]
            out.append(_rapid_fs(sd2, circle))
        rows[param] = tuple(out)
    sigma_f = math.sqrt(sum(((fp - fm) / 2.0) ** 2 for fp, fm in rows.values()))
    cov = sigma_f / f_mlv
    beta = (math.log(f_mlv / math.sqrt(1.0 + cov ** 2))
            / math.sqrt(math.log(1.0 + cov ** 2)))
    return {"F_MLV": f_mlv, "rows": rows, "sigma_F": sigma_f, "beta_ln": beta}


def check_rapid_reliability_taylor():
    from xslope.reliability import reliability_taylor
    fails = []
    sd = _load(RAPID_MODEL)
    sd["materials"][0].update(RAPID_SIGMAS)
    hand = rapid_reliability_by_hand(sd)
    ok, res = _quiet(reliability_taylor, sd, "spencer", rapid=True, search=False)
    if not ok:
        return [f"rapid reliability: the Taylor series refused a rapid drawdown "
                f"run on the fixed surface ({res})"]
    got = {p["param"]: (p["F_plus"], p["F_minus"]) for p in res["param_info"]}
    if set(got) != set(hand["rows"]):
        fails.append(f"rapid reliability: perturbed {sorted(got)}, expected "
                     f"{sorted(hand['rows'])}")
    for param, (fp, fm) in hand["rows"].items():
        if param in got and (abs(got[param][0] - fp) > 1e-9
                             or abs(got[param][1] - fm) > 1e-9):
            fails.append(f"rapid reliability: {param} F+/F- = {got[param]}, by hand "
                         f"({fp}, {fm})")
    for key in ("F_MLV", "sigma_F", "beta_ln"):
        if abs(float(res[key]) - hand[key]) > 1e-9:
            fails.append(f"rapid reliability: {key} = {res[key]!r}, by hand "
                         f"{hand[key]!r}")
    # Without drawdown nothing reads d or psi: they stay out of the set.
    ok, res = _quiet(reliability_taylor, sd, "spencer", rapid=False, search=False)
    if not ok:
        fails.append(f"rapid reliability: the run without drawdown failed ({res})")
    else:
        params = sorted(p["param"] for p in res["param_info"])
        if params != ["phi"]:
            fails.append(f"rapid reliability: a run without drawdown perturbed "
                         f"{params}, expected only phi")
    return fails


CHECKS = [
    ("seep bc: every block read", check_seep_bc_blocks),
    ("unsat: words kept, any case", check_unsat_words),
    ("non-circ Movement: loader", check_movement_loader),
    ("non-circ Movement: what each entry does", check_movement_meaning),
    ("non-circ Movement: search refuses", check_movement_search_refuses),
    ("circles Option: any case", check_circle_option_case),
    ("rapid: d/psi on mc rows only", check_rapid_d_psi_rows),
    ("rapid reliability: d and psi perturbed", check_rapid_reliability_taylor),
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
        print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("audit checks:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll audit checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
