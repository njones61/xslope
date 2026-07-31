"""Standing checks on the remedy panel when a fault has more than one repair.

A rule can offer several remedies (``Rule.remedies``, primary first), because a
fault can have more than one right answer and the rule cannot tell which: an empty
surface sheet takes either a starting set of trial circles, where the slope's own
geometry controls the mechanism, or a surface tracked along a weak seam, where a
layer does — and no circle can follow a seam. Studio's model-checks panel is where
that choice is put to the user, so it has to render every remedy the finding
carries rather than only the first. What is guarded here:

  A. THE FINDING CARRIES THEM — a finding from a multi-remedy rule reports every
     remedy its rule offers, with the primary first and equal to ``remedy``, and a
     single-remedy rule still reports exactly one. This is the wire between the
     registry and anything that renders from a finding.
  B. TWO BUTTONS — a model with no surface at all shows both repairs, each with its
     own label from the remedy metadata, each live. One button here means the panel
     is still rendering only the primary and the second repair is unreachable.
  C. THE AVAILABILITY CONDITION — the weak-zone surface appears only where the
     generator picks a zone on its own. On a model with two comparable candidates
     the question belongs to Studio's zone picker, so that button is absent
     entirely — not dimmed — while the circles button still stands. This is the leg
     that fails if the panel starts rendering a button per DECLARED remedy rather
     than per PROPOSAL.
  D. THE CONTRACT IS UNCHANGED — each button still goes propose -> confirm -> apply:
     the proposal states its change before anything happens, applying it lands on
     the model, and the panel re-checks so the finding that offered it is gone.
     Asserted for the SECOND remedy, which is the one that never had a button.
  E. THE DIM RULE IS UNCHANGED — a remedy that is offered but cannot be applied is
     still rendered disabled with its reason as the tooltip.

The remedies themselves — what each one changes, and the gate that decides whether
the weak-zone surface is offered at all — are checked in run_tests.py's
preflight_remedies row. Nothing here re-checks the transformation.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import io
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: A slope with one clearly weakest seam: the generator picks that zone on its own,
#: so both repairs to an empty surface sheet are live.
SEAM = os.path.join(_REPO, "docs/lem/files/xslope_acads_weak_layer.xlsx")
#: Two comparable candidate zones: the generator declines to choose, so the
#: weak-zone surface stands down and only the circles remain.
AMBIGUOUS = os.path.join(_REPO, "docs/verification/files/rocscience/vp063.xlsx")

#: The rule with two repairs, and the two remedies it offers.
RULE = "surface.none_defined"
PRIMARY = "generate_starting_circles"
SECOND = "generate_noncircular_surface"


def _quiet(fn, *args, **kwargs):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        return fn(*args, **kwargs)


def _no_surface(path):
    """The model with both surface sheets emptied — the fault the rule fires on."""
    from xslope.fileio import load_slope_data
    sd = dict(_quiet(load_slope_data, path))
    sd["circles"] = []
    sd["non_circ"] = []
    sd["circular"] = False
    return sd


def _panel(sd):
    from studio.preflight_panel import PreflightPanel
    return _quiet(PreflightPanel, "lem", slope_data=sd)


def _remedy_buttons(panel):
    """``{remedy key: QPushButton}`` for every remedy button the panel rendered."""
    from PySide6.QtWidgets import QPushButton
    out = {}
    for b in panel.findChildren(QPushButton):
        name = b.objectName()
        if name.startswith("remedy_"):
            out[name[len("remedy_"):]] = b
    return out


# ------------------------------------------------------------------ A. finding
def test_finding_carries_remedies():
    fails = []
    from xslope import preflight as pf

    found = [f for f in pf.preflight(_no_surface(SEAM), "lem") if f.rule_id == RULE]
    if not found:
        fails.append(f"{RULE} did not fire on a model with no surface at all")
        return fails
    f = found[0]
    if tuple(f.remedies) != (PRIMARY, SECOND):
        fails.append(f"the finding carries remedies={tuple(f.remedies)!r}, want "
                     f"{(PRIMARY, SECOND)!r} — primary first")
    if f.remedy != (tuple(f.remedies)[0] if f.remedies else None):
        fails.append(f"the finding's primary remedy {f.remedy!r} is not the first "
                     f"of its remedies {tuple(f.remedies)!r}")

    # A finding must report exactly what its rule declares — checked against the
    # registry, so a rule that gains or loses a remedy cannot leave this stale.
    by_id = {r.id: r for r in pf.rules()}
    for other in pf.preflight(_no_surface(SEAM), "lem"):
        rule = by_id.get(other.rule_id)
        if rule is None:
            continue
        if tuple(other.remedies) != tuple(rule.remedies):
            fails.append(f"{other.rule_id}: the finding reports "
                         f"{tuple(other.remedies)!r} but the rule declares "
                         f"{tuple(rule.remedies)!r}")
    return fails


# ------------------------------------------------------------- B. two buttons
def test_two_buttons():
    fails = []
    from xslope.preflight import REMEDIES

    panel = _panel(_no_surface(SEAM))
    buttons = _remedy_buttons(panel)
    for remedy in (PRIMARY, SECOND):
        hit = [k for k in buttons if k.startswith(remedy)]
        if not hit:
            fails.append(f"a model with no surface at all offers no {remedy!r} "
                         f"button (rendered: {sorted(buttons)})")
            continue
        btn = buttons[hit[0]]
        if not btn.isEnabled():
            fails.append(f"the {remedy!r} button is dimmed on a model where the "
                         f"remedy applies cleanly (tooltip {btn.toolTip()!r})")
        if not btn.text().strip():
            fails.append(f"the {remedy!r} button carries no label")
        if btn.text().strip().rstrip("…").strip() == remedy:
            fails.append(f"the {remedy!r} button is labelled with the remedy's own "
                         f"name rather than words from its metadata")
    if len(buttons) < 2:
        fails.append(f"the panel rendered {len(buttons)} remedy button(s) for a "
                     f"fault with two repairs — the second repair is unreachable")
    # Each label has to be distinguishable from the other, or two buttons is worse
    # than one.
    labels = {k: buttons[k].text() for k in buttons}
    if len(set(labels.values())) != len(labels):
        fails.append(f"two remedy buttons carry the same label: {labels}")
    if not REMEDIES.get(SECOND):
        fails.append(f"{SECOND!r} has no entry in REMEDIES to label a button from")
    return fails


# ------------------------------------------------------ C. availability gate
def test_availability_gate():
    fails = []
    from xslope.generators import generate_noncircular_surface

    sd = _no_surface(AMBIGUOUS)
    picked = _quiet(generate_noncircular_surface, sd, report=True)
    if picked["zone"] is not None or len(picked["candidates"]) < 2:
        fails.append("the ambiguity fixture is no longer ambiguous — the generator "
                     "picked a zone, so this leg proves nothing")
        return fails

    buttons = _remedy_buttons(_panel(sd))
    if any(k.startswith(SECOND) for k in buttons):
        fails.append(f"the weak-zone surface was offered on a model with no clearly "
                     f"weakest zone; that case belongs to the zone picker, and a "
                     f"remedy states its change rather than guessing at one "
                     f"(rendered: {sorted(buttons)})")
    if not any(k.startswith(PRIMARY) for k in buttons):
        fails.append(f"the circles button disappeared along with the weak-zone one "
                     f"(rendered: {sorted(buttons)}) — one remedy standing down must "
                     f"not take its sibling with it")
    return fails


# --------------------------------------------------------------- D. contract
def test_contract():
    fails = []
    from xslope import remedies as rem

    sd = _no_surface(SEAM)
    panel = _panel(sd)
    if not panel.blocked:
        fails.append("a model with no surface at all did not block the run")

    # propose: the change is computed and stated BEFORE anything happens.
    proposal = rem.propose(sd, SECOND)
    if not proposal.available:
        fails.append(f"the weak-zone remedy declined on the seam model: "
                     f"{proposal.reason}")
        return fails
    if not proposal.description:
        fails.append("the proposal offers no description, so a confirmation dialog "
                     "would ask the user to approve an unstated change")
    before = list(sd.get("non_circ") or [])

    # apply: through the panel, which is what a confirmed button does.
    _quiet(panel.apply_proposal, proposal)
    if not sd.get("non_circ"):
        fails.append("applying the weak-zone remedy left the model without a "
                     "non-circular surface")
    if before:
        fails.append("the fixture already carried a surface, so this leg proves "
                     "nothing")
    if panel.blocked:
        fails.append(f"the panel still blocks the run after the repair was applied "
                     f"({panel.block_reason()[:80]!r})")
    if any(f.rule_id == RULE for f in (panel.report.errors if panel.report else [])):
        fails.append("the finding that offered the remedy survived it — the panel "
                     "did not re-check the model it changed")
    return fails


# --------------------------------------------------------------- E. dim rule
def test_dim_rule():
    fails = []
    from xslope import remedies as rem

    # Any dimmed proposal will do: the rule under test is that a proposal which is
    # offered-but-unavailable renders as a disabled button carrying its reason.
    sd = _no_surface(SEAM)
    panel = _panel(sd)
    dimmed = [p for p in rem.remedy_proposals(sd) if not p.available]
    buttons = _remedy_buttons(panel)
    checked = 0
    for p in dimmed:
        btn = buttons.get(p.key)
        if btn is None:
            continue                     # its rule did not fire on this model
        checked += 1
        if btn.isEnabled():
            fails.append(f"{p.key} is unavailable but its button is live")
        if btn.toolTip() != p.reason:
            fails.append(f"{p.key} is dimmed without its reason as the tooltip "
                         f"({btn.toolTip()!r} vs {p.reason!r})")
    # ...and a live one carries its description instead, which is the other half of
    # the same rule.
    for p in rem.remedy_proposals(sd):
        btn = buttons.get(p.key)
        if btn is None or not p.available:
            continue
        checked += 1
        if btn.toolTip() != p.description:
            fails.append(f"{p.key} is live but its tooltip is not the change it "
                         f"would make ({btn.toolTip()!r})")
    if checked == 0:
        fails.append("no rendered button could be matched to a proposal, so the "
                     "dim rule was not exercised at all")
    return fails


CHECKS = [("a finding carries every remedy", test_finding_carries_remedies),
          ("two repairs, two buttons", test_two_buttons),
          ("the second button's availability gate", test_availability_gate),
          ("propose -> confirm -> apply, second remedy", test_contract),
          ("the dim rule, unchanged", test_dim_rule)]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("remedy panel: PySide6 not installed — checks skipped.")
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
    print("multi-remedy model-checks panel (Studio):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll remedy panel checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
