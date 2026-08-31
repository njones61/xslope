"""Standing checks on the value fields of Studio's Seep BC editor — that a head and a
flux read a typed value by the same rule, and that Set 2 holds both of them to numbers.

A seepage boundary's value cell may name a tseep series instead of holding a number:
that is how a time-varying boundary is written, and the loader resolves the name
against the tseep sheet's series. The head field honored that; the flux field did not.
Typing ``storm`` into "Flux value:" raised ValueError on the float() and the editor
answered it by storing 0.0 — the boundary the user had just made time-varying became a
no-flow boundary, with no message, and the model saved that way. Every downstream
symptom (a transient march with no rain in it, a steady solve that matched the dry
case) points at the physics rather than at the field that dropped the name.

  A. A FLUX TAKES A SERIES NAME. Typed into Set 1's flux field, it survives _commit
     as the string it was typed as, exactly as the head field keeps one.
  B. IT SURVIVES THE FILE. The name written to the workbook comes back from the
     loader as the same string, resolved against the tseep sheet's series — so the
     editor's answer and the file format's answer are the same answer.
  C. SET 2 IS STILL CONSTANT. 'seep bc (2)' is the constant-steady rapid-drawdown
     set and fileio rejects a series there. The editor coerces the name to 0.0 in
     that set rather than writing a file that will not load — the same rule the head
     field applies, now applied to fluxes.
  D. NUMBERS ARE STILL NUMBERS. A plain number, in either notation, still parses to
     a float, an untouched display rounding still keeps the stored value, and a
     cleared field still reads 0.

Skips cleanly (exit 0) when PySide6 is not installed.
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

#: The SEEP-4 tutorial's model: one specified head, three specified-flux blocks along
#: the dam crest and shoulders, an exit face. The fluxes are the boundaries the
#: tutorial makes time-varying, so this is the file the defect was met in.
DAM = os.path.join(_REPO, "docs/tutorials/files/xslope_dam_infiltration.xlsx")

#: The series name the tutorial uses.
SERIES = "storm"


def _fail(cond, message):
    return [] if cond else [message]


def _load(path):
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()):
        return load_slope_data(path)


def _widget(bc, constant_only=False):
    """One BC set widget with no preview pane (the preview needs a full model and
    draws nothing this file has an opinion about)."""
    from studio.editors import _SeepBcSetWidget
    return _SeepBcSetWidget(bc, constant_only=constant_only)


def _type_flux(w, flux_no, text):
    """Select flux ``flux_no`` (1-based), type ``text`` into "Flux value:", and read
    the set back the way the dialog's OK does."""
    w.list.setCurrentRow(len(w._heads) + flux_no - 1)
    w.flux_edit.setText(text)
    return w.result()


def _storm_series(times):
    """A minimal tseep block defining one series named ``storm``."""
    return {"times": list(times),
            "series": {SERIES: [0.0] * len(times)},
            "dt": 1.0, "duration": float(times[-1]), "save_interval": 1.0,
            "save_times": []}


# ---------------------------------------------------------------------------
# A. A flux value takes a tseep series name
# ---------------------------------------------------------------------------
def check_flux_keeps_series_name():
    d = _load(DAM)
    w = _widget(d["seepage_bc"])
    out = _type_flux(w, 1, SERIES)
    got = out["specified_fluxes"][0]["flux"]
    failures = _fail(got == SERIES,
                     f"a series name typed into Flux value came back as {got!r}, "
                     f"not {SERIES!r}")
    # The other two flux blocks are untouched and keep their numbers.
    rest = [f["flux"] for f in out["specified_fluxes"][1:]]
    failures += _fail(all(isinstance(v, float) for v in rest),
                      f"untouched flux blocks changed type: {rest!r}")
    # The head field's behavior is the one being mirrored — assert it still holds.
    wh = _widget(d["seepage_bc"])
    wh.list.setCurrentRow(0)
    wh.head_edit.setText(SERIES)
    got_h = wh.result()["specified_heads"][0]["head"]
    failures += _fail(got_h == SERIES,
                      f"a series name typed into Head value came back as {got_h!r}")
    # Whitespace is trimmed, the same as a head's.
    w2 = _widget(d["seepage_bc"])
    got2 = _type_flux(w2, 2, f"  {SERIES}  ")["specified_fluxes"][1]["flux"]
    failures += _fail(got2 == SERIES,
                      f"a padded series name came back as {got2!r}, untrimmed")
    return failures


# ---------------------------------------------------------------------------
# B. The name survives save and load
# ---------------------------------------------------------------------------
def check_series_round_trips():
    from xslope.fileio import save_slope_data_to_xlsx

    d = _load(DAM)
    w = _widget(d["seepage_bc"])
    # All three flux blocks bound to the one series, which is what the tutorial does.
    for n in (1, 2, 3):
        bc = _type_flux(w, n, SERIES)
    d["seepage_bc"] = bc
    d["tseep"] = _storm_series([0.0, 1.0, 2.0, 3.0])
    # A tseep sheet makes the file a transient model, and the loader requires storage
    # on every material for one. Not what this check is about — just enough for the
    # file to load back.
    for m in d["materials"]:
        m["Ss"], m["Sy"] = 1e-4, 0.20

    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "seep_bc_series.xlsx")
        with contextlib.redirect_stdout(io.StringIO()):
            save_slope_data_to_xlsx(d, out)
        back = _load(out)
        got = [f["flux"] for f in back["seepage_bc"]["specified_fluxes"]]
        failures += _fail(got == [SERIES] * 3,
                          f"after save/load the flux values are {got!r}")
        failures += _fail(SERIES in (back.get("tseep") or {}).get("series", {}),
                          "the tseep sheet lost the series the fluxes name")
    return failures


# ---------------------------------------------------------------------------
# C. Set 2 stays constant
# ---------------------------------------------------------------------------
def check_set2_coerces_to_zero():
    d = _load(DAM)
    # Set 2 built from the same boundaries, as the rapid-drawdown set would be.
    w = _widget(d["seepage_bc"], constant_only=True)
    got = _type_flux(w, 1, SERIES)["specified_fluxes"][0]["flux"]
    failures = _fail(got == 0.0,
                     f"Set 2 kept a series name in a flux value ({got!r}); it must "
                     f"coerce to 0.0, since fileio refuses to load one there")
    # The head field's Set-2 rule is unchanged.
    wh = _widget(d["seepage_bc"], constant_only=True)
    wh.list.setCurrentRow(0)
    wh.head_edit.setText(SERIES)
    got_h = wh.result()["specified_heads"][0]["head"]
    failures += _fail(got_h == 0.0,
                      f"Set 2 kept a series name in a head value ({got_h!r})")
    # And the field says why, so a coerced value is not the user's first hint.
    w2 = _widget(d["seepage_bc"], constant_only=True)
    tip = w2.flux_edit.toolTip()
    failures += _fail("Set 2" in tip and "tseep" in tip,
                      f"Set 2's flux field carries no explanation (tooltip {tip!r})")
    return failures


# ---------------------------------------------------------------------------
# D. Numbers still parse
# ---------------------------------------------------------------------------
def check_numbers_still_parse():
    d = _load(DAM)
    failures = []
    for text, want in (("1e-08", 1e-08),
                       ("0.00000001", 1e-08),
                       ("-2.5", -2.5),
                       ("", 0.0)):
        w = _widget(d["seepage_bc"])
        got = _type_flux(w, 1, text)["specified_fluxes"][0]["flux"]
        failures += _fail(isinstance(got, float) and got == want,
                          f"flux text {text!r} parsed to {got!r}, not {want!r}")
    # An untouched field is a display rounding nobody typed over: the stored value
    # must survive, not the six digits the widget showed.
    stored = d["seepage_bc"]["specified_fluxes"][0]["flux"]
    w = _widget(d["seepage_bc"])
    w.list.setCurrentRow(len(w._heads))          # select it, type nothing
    got = w.result()["specified_fluxes"][0]["flux"]
    failures += _fail(got == stored,
                      f"an untouched flux moved from {stored!r} to {got!r}")
    return failures


CHECKS = [
    ("A. flux value keeps a tseep series name", check_flux_keeps_series_name),
    ("B. the series name survives save/load", check_series_round_trips),
    ("C. set 2 coerces a series name to 0", check_set2_coerces_to_zero),
    ("D. numbers still parse", check_numbers_still_parse),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("seep BC series: PySide6 not installed — skipped.")
        return []
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
    print("seep BC series:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll seep BC series checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
