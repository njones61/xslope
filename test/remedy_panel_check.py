"""Standing checks on Studio's model-checks panel: its remedies and its shape.

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

The other half is the panel's SHAPE. A rule that fires once per material used to
print one full paragraph per material — same explanation, same remedy, four times
over — which made every Run dialog tall and narrow. The presentation rules that
replaced that are checked here too, because each of them is invisible to the
registry and so has nothing else standing under it:

  F. ONE RULE, ONE LINE — findings that share a rule id are one entry, whose line
     names how many there are and which, and whose detail states the explanation
     they share ONCE with the per-finding specifics under it. Four findings that
     render four rows is the regression this stands against.
  G. A REMEDY STILL APPLIES THROUGH THE AGGREGATE — a rule that fires twice and
     offers a repair for each target carries both buttons on its single entry, and
     pressing one still goes propose -> apply -> re-check.
  H. TWO PANES — the run controls and the model checks are columns beside each
     other, not one above the other, and the dialog is wider than it is tall.
  I. THE DETAIL FOLLOWS THE SELECTION — the pane opens on an error rather than on
     nothing, and selecting another line shows that line's text.
  J. THE NOTES ARE STILL COLLAPSED — infos stay behind their counted disclosure
     line and come back when it is opened.

The remedies themselves — what each one changes, and the gate that decides whether
the weak-zone surface is offered at all — are checked in run_tests.py's
preflight_remedies row. Nothing here re-checks the transformation.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import html
import io
import os
import re
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

#: A four-material model with no seepage properties entered: three separate rules
#: each fire once per material, which is the shape (one rule, several findings) the
#: aggregation exists for.
FOUR_MATERIALS = os.path.join(_REPO, "docs/inputs/slope/xslope_dam.xlsx")
#: An earth dam whose two seepage materials leave the anisotropy angle blank: two
#: findings of ONE rule, at INFO, so they sit behind the notes disclosure.
TWO_NOTES = os.path.join(_REPO, "docs/seep/files/xslope_earth_dam_tseep.xlsx")
#: A model that preflights clean, for the header with nothing to count.
CLEAN = os.path.join(_REPO, "docs/lem/files/xslope_design.xlsx")

#: What the header appends once there is a selectable line: the list-then-detail
#: interaction is not self-evident from the list alone.
DETAIL_HINT = "(select a line for details)"


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


def _panel(sd, analysis="lem"):
    from studio.preflight_panel import PreflightPanel
    return _quiet(PreflightPanel, analysis, slope_data=sd)


def _load(path):
    from xslope.fileio import load_slope_data
    return dict(_quiet(load_slope_data, path))


def _rows(panel):
    """``[(rule_id, line, hidden), ...]`` — the list as a user reads it."""
    from PySide6.QtWidgets import QLabel, QListWidget
    lst = panel.findChild(QListWidget, "preflight_list")
    out = []
    for i in range(lst.count()):
        item = lst.item(i)
        widget = lst.itemWidget(item)
        labels = [w for w in widget.findChildren(QLabel)
                  if w.objectName().startswith("finding_")]
        out.append((labels[0].objectName()[len("finding_"):] if labels else "",
                    labels[0].text() if labels else "", item.isHidden()))
    return out


def _detail(panel):
    """``(rule_id, text)`` of the detail the panel is showing right now."""
    from PySide6.QtWidgets import QLabel, QStackedWidget
    stack = panel.findChild(QStackedWidget, "preflight_detail")
    page = stack.currentWidget()
    labels = [w for w in page.findChildren(QLabel)
              if w.objectName().startswith("detail_")]
    if not labels:
        return "", ""
    text = re.sub(r"<[^>]+>", " ", labels[0].text())
    return labels[0].objectName()[len("detail_"):], html.unescape(text)


def _shared_tail(messages):
    """The sentence every one of ``messages`` ends with, read independently of the
    panel's own splitter so the check cannot agree with a broken one."""
    parts = [[s for s in m.replace("! ", ". ").split(". ") if s.strip()]
             for m in messages]
    last = {p[-1] for p in parts}
    return last.pop() if len(last) == 1 else ""


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
    if any(f.rule_id == RULE for f in (panel.report.errors if panel.report is not None else [])):
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


# ------------------------------------------------------- F. one rule, one line
def test_one_line_per_rule():
    fails = []
    from PySide6.QtWidgets import QListWidget
    from xslope import preflight as pf

    sd = _load(FOUR_MATERIALS)
    by_rule = {}
    for f in pf.preflight(sd, "seep", {}):
        by_rule.setdefault(f.rule_id, []).append(f)
    repeats = {rid: fs for rid, fs in by_rule.items() if len(fs) >= 4}
    if not repeats:
        fails.append("no rule fires four times on the fixture any more, so nothing "
                     "here exercises the aggregation")
        return fails

    names = [str(m.get("name") or "") for m in (sd.get("materials") or [])]
    shared_seen = 0            # rules whose findings DO share a closing sentence
    panel = _panel(sd, "seep")
    rows = _rows(panel)
    lst = panel.findChild(QListWidget, "preflight_list")
    for rid, fs in repeats.items():
        mine = [r for r in rows if r[0] == rid]
        if len(mine) != 1:
            fails.append(f"{rid} fired {len(fs)} times and took {len(mine)} list "
                         f"rows — findings of one rule are one entry")
            continue
        line = mine[0][1]
        if str(len(fs)) not in line:
            fails.append(f"{rid}: the line does not say how many findings stand "
                         f"behind it: {line!r}")
        if "\n" in line:
            fails.append(f"{rid}: the entry is not one line: {line!r}")
        named = [n for n in names if n and any(n in f.message for f in fs)]
        for n in named:
            if n not in line:
                fails.append(f"{rid}: the line names {len(fs)} findings but not "
                             f"which — {n!r} is missing from {line!r}")

        lst.setCurrentRow(rows.index(mine[0]))
        shown, text = _detail(panel)
        if shown != rid:
            fails.append(f"{rid}: selecting its line showed {shown!r}'s detail")
            continue
        # Where the findings of a rule end on the same sentence -- the explanation
        # and the advice that are the rule's, not the material's -- that sentence
        # is stated once. (Not every rule has one: some end on the values they
        # read, which differ per finding and belong in the sub-list.)
        tail = _shared_tail([f.message for f in fs])
        if tail:
            shared_seen += 1
            if text.count(tail) != 1:
                fails.append(f"{rid}: the sentence all {len(fs)} findings share "
                             f"appears {text.count(tail)} times in the detail, "
                             f"want once")
        for n in named:
            if n not in text:
                fails.append(f"{rid}: the detail dropped {n!r} — the specifics of "
                             f"every finding have to survive the aggregation")
    if not shared_seen:
        fails.append("no repeated rule on the fixture shares a closing sentence, so "
                     "nothing here tells one explanation from four")
    return fails


# ------------------------------------------ G. a remedy through the aggregate
def test_remedy_through_aggregate():
    fails = []
    from xslope import remedies as rem

    sd = _load(FOUR_MATERIALS)
    pts = [(float(x), float(y)) for x, y in sd["piezo_line"]]
    sd["piezo_line"] = list(reversed(pts))
    sd["piezo_line2"] = list(reversed(pts))       # one rule, two targets

    panel = _panel(sd)
    rows = [r for r in _rows(panel) if r[0] == "order.piezo_reversed"]
    if len(rows) != 1:
        fails.append(f"two reversed piezometric lines took {len(rows)} rows")
    buttons = _remedy_buttons(panel)
    for key in ("reverse_polyline:piezo1", "reverse_polyline:piezo2"):
        if key not in buttons:
            fails.append(f"the aggregated entry dropped the {key} button "
                         f"(rendered: {sorted(buttons)}) — a rule's repairs are "
                         f"offered once each, not once per finding")
        elif not buttons[key].isEnabled():
            fails.append(f"{key} is dimmed on a line that reverses cleanly")
    if fails:
        return fails

    # ...and applying one still lands: propose -> apply -> re-check.
    _quiet(panel.apply_proposal, rem.propose(sd, "reverse_polyline", "piezo1"))
    if [tuple(p) for p in panel.model["piezo_line"]] != pts:
        fails.append("the remedy offered by an aggregated entry did not reverse "
                     "the line")
    left = [r for r in _rows(panel) if r[0] == "order.piezo_reversed"]
    if len(left) != 1:
        fails.append(f"after one of two lines was fixed the rule took {len(left)} "
                     f"rows, want the one line that is still reversed")
    keys = set(_remedy_buttons(panel))
    if "reverse_polyline:piezo1" in keys:
        fails.append("the repaired line is still offering its repair")
    if "reverse_polyline:piezo2" not in keys:
        fails.append("the line that is still reversed lost its repair")
    return fails


# ------------------------------------------------------------- H. two panes
def test_two_panes():
    fails = []
    from PySide6.QtWidgets import QApplication, QGridLayout
    from studio.dialogs import RunLemDialog

    app = QApplication.instance() or QApplication([])
    sd = _load(FOUR_MATERIALS)
    pts = [(float(x), float(y)) for x, y in sd["piezo_line"]]
    sd["piezo_line"] = list(reversed(pts))
    sd["piezo_line2"] = list(reversed(pts))

    dlg = _quiet(RunLemDialog, None, slope_data=sd)
    dlg.resize(dlg.sizeHint())
    dlg.show()
    app.processEvents()
    grid = dlg.layout()
    if not isinstance(grid, QGridLayout):
        fails.append(f"the Run dialog is laid out as {type(grid).__name__}, not as "
                     f"columns")
        dlg.close()
        return fails
    right = grid.itemAtPosition(0, 1)
    if right is None:
        fails.append("the dialog has no second column — the checks are stacked "
                     "under the controls rather than beside them")
        dlg.close()
        return fails
    controls = grid.itemAtPosition(0, 0).widget()
    checks = right.widget()
    if checks is not dlg.preflight:
        fails.append("the right column is not the model checks")
    if not (controls.isVisible() and checks.isVisible()):
        fails.append(f"controls visible={controls.isVisible()}, checks "
                     f"visible={checks.isVisible()} — both are always on screen")
    if not dlg.method.isVisible():
        fails.append("the run controls are not in the left column")
    if checks.x() < controls.x() + controls.width():
        fails.append(f"the checks (x={checks.x()}) do not sit beside the controls "
                     f"(x={controls.x()}, w={controls.width()})")
    if dlg.width() <= dlg.height():
        fails.append(f"the dialog is {dlg.width()}x{dlg.height()} — taller than "
                     f"wide on an ordinary warning set")
    dlg.close()

    # The gate is unchanged by the new shape: an error still refuses the run, in
    # the words the run would have been refused with.
    blocked = _quiet(RunLemDialog, None, slope_data=_no_surface(SEAM))
    if blocked._ok.isEnabled():
        fails.append("Run stayed live on a model with no surface at all")
    if blocked._ok.toolTip() != blocked.preflight.block_reason():
        fails.append("the disabled Run button no longer carries the gate's reason")
    blocked.close()
    return fails


# --------------------------------------------- I. the detail follows selection
def test_detail_follows_selection():
    fails = []
    from PySide6.QtWidgets import QListWidget

    panel = _panel(_load(FOUR_MATERIALS), "seep")
    rows = _rows(panel)
    if len(rows) < 2:
        fails.append("the fixture reports one entry, so selection cannot be tested")
        return fails
    lst = panel.findChild(QListWidget, "preflight_list")
    if lst.currentRow() < 0:
        fails.append("the checks opened with nothing selected, so the detail pane "
                     "is blank")
        return fails
    groups = panel.groups()
    if any(g.severity == "error" for g in groups) and \
            groups[lst.currentRow()].severity != "error":
        fails.append("the pane did not open on the error that refuses the run")
    first, first_text = _detail(panel)
    if first != rows[lst.currentRow()][0]:
        fails.append(f"the detail shows {first!r} while {rows[0][0]!r} is selected")
    lst.setCurrentRow(len(rows) - 1)
    last, last_text = _detail(panel)
    if last != rows[-1][0]:
        fails.append(f"selecting the last line still shows {last!r}'s detail")
    if last_text == first_text:
        fails.append("the detail text did not change with the selection")
    return fails


# ------------------------------------------------------ J. the notes disclosure
def test_notes_collapsed():
    fails = []
    from PySide6.QtWidgets import QToolButton

    panel = _panel(_load(TWO_NOTES), "seep")
    toggle = panel.findChild(QToolButton, "preflight_infos")
    infos = [f for f in (panel.report.infos if panel.report is not None else [])]
    if not infos:
        fails.append("the fixture reports no notes, so the disclosure is not tested")
        return fails
    if toggle is None or toggle.isHidden():
        fails.append("a model with notes shows no disclosure line for them")
        return fails
    if str(len(infos)) not in toggle.text():
        fails.append(f"the disclosure does not count the notes: {toggle.text()!r}")
    rule_ids = {f.rule_id for f in infos}
    note_rows = [r for r in _rows(panel) if r[0] in rule_ids]
    if not note_rows:
        fails.append("the notes are not in the list at all")
        return fails
    if len(note_rows) != len(rule_ids):
        fails.append(f"{len(infos)} notes of {len(rule_ids)} rule(s) took "
                     f"{len(note_rows)} rows")
    if not all(hidden for _r, _t, hidden in note_rows):
        fails.append("the notes are showing before their disclosure was opened")
    toggle.setChecked(True)
    if any(hidden for _r, _t, hidden in _rows(panel) if _r in rule_ids):
        fails.append("opening the disclosure did not show the notes")
    toggle.setChecked(False)
    if not all(hidden for _r, _t, hidden in _rows(panel) if _r in rule_ids):
        fails.append("closing the disclosure did not put the notes away")
    return fails


# --------------------------------------------- K. the header counts the LIST
def test_header_counts_entries():
    """The header counts the entries on screen, not the findings behind them.

    The reported bug: "Model checks — 5 warnings" over a list of three lines. The
    header was counting raw findings while the list aggregates a rule's findings
    into one entry, and the multiplicity it was counting is already stated on the
    entry itself ("— 4 materials"). So the two numbers disagreed about the same
    screen, and the larger one was the one in the title bar.

    The fixture is chosen so the two counts CANNOT coincide: several rules fire once
    per material there, so raw > aggregated for at least one severity, and a header
    that went back to counting findings fails this leg rather than passing it by
    accident.
    """
    fails = []
    sd = _load(FOUR_MATERIALS)
    panel = _panel(sd, "seep")
    report = panel.report
    shown = [g for g in panel.groups() if g is not None and not _hidden_row(panel, g)]

    want = []
    for severity, one, many in (("error", "error", "errors"),
                                ("warning", "warning", "warnings"),
                                ("info", "note", "notes")):
        n = sum(1 for g in shown if g.severity == severity)
        if n:
            want.append(f"{n} {one if n == 1 else many}")
    if not want:
        fails.append("the fixture reports nothing, so the header is not exercised")
        return fails
    expected = "Model checks — " + " · ".join(want) + "   " + DETAIL_HINT
    if panel.title() != expected:
        fails.append(f"the header is {panel.title()!r}, want {expected!r}")

    # The mutation this leg exists to catch: counting the findings instead.
    raw = {"error": len(report.errors), "warning": len(report.warnings)}
    agg = {s: sum(1 for g in shown if g.severity == s) for s in raw}
    if not any(raw[s] > agg[s] for s in raw if agg[s]):
        fails.append("no severity on this fixture aggregates several findings into "
                     "one line, so the header could count either way and pass")
    for severity in raw:
        if agg[severity] and raw[severity] > agg[severity]:
            if f"{raw[severity]} " in panel.title():
                fails.append(f"the header carries the raw {severity} count "
                             f"({raw[severity]}) over a list of {agg[severity]} "
                             f"line(s): {panel.title()!r}")

    # Nothing to select: no counts, and no hint pointing at a line that does not open.
    clean = _panel(_load(CLEAN))
    if clean.report is not None and clean.report.findings:
        fails.append("the clean fixture is no longer clean, so the all-clear header "
                     "is not exercised")
    elif clean.title() != "Model checks":
        fails.append(f"the all-clear header is {clean.title()!r}, want plain "
                     f"'Model checks' -- its one line is not selectable")
    return fails


def _hidden_row(panel, group):
    """True when ``group``'s row is behind the notes disclosure."""
    from PySide6.QtWidgets import QListWidget
    lst = panel.findChild(QListWidget, "preflight_list")
    groups = panel.groups()
    try:
        row = groups.index(group)
    except ValueError:
        return False
    item = lst.item(row)
    return bool(item is not None and item.isHidden())


CHECKS = [("a finding carries every remedy", test_finding_carries_remedies),
          ("two repairs, two buttons", test_two_buttons),
          ("the second button's availability gate", test_availability_gate),
          ("propose -> confirm -> apply, second remedy", test_contract),
          ("the dim rule, unchanged", test_dim_rule),
          ("one rule, one line", test_one_line_per_rule),
          ("a remedy through the aggregate", test_remedy_through_aggregate),
          ("two panes: controls | checks", test_two_panes),
          ("the detail follows the selection", test_detail_follows_selection),
          ("the notes stay collapsed", test_notes_collapsed),
          ("the header counts what the list shows", test_header_counts_entries)]


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
    print("model-checks panel (Studio):")
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
