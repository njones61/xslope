"""Standing checks on the assistant's guardrails — the rules it is told, and the
check it cannot skip.

Three live sessions established the failure mode this file exists to prevent: a
modeling rule stated in the prompt is a rule the model follows MOST of the time.
It was measured losing to guidance elsewhere in the prompt stack, to a
confident-sounding tool error that prescribed the violation, and to the
edit-cascade — a one-field correction whose dependents go stale in silence (max
depth is fixed, and the circles that were tangent to the old base are now under
the new one, with nothing saying so).

Prompt text cannot be unit-tested: what a model does with a sentence is measured
in sessions, not asserted here. What CAN be asserted is everything the harness
does regardless of the model, and that is what this file covers:

  A. THE IRON RULES REACH EVERY PROVIDER — the block lives in STUDIO_SYSTEM, the
     short prompt assembled for every tier, so it is present exactly once whether
     the full skill is loaded (Anthropic and friends) or not (a local model), and
     no rule of it is duplicated by the compact rulebook that only the skill-less
     path gets. Duplication is not harmless: two copies drift, and the drifted
     one is the one that gets read.
  B. THE MODELING BRIEF AGREES WITH THEM — the non-skill providers' rulebook
     carries the circle floor rule and the rigid-base corollary, since its extents
     paragraph used to say the opposite of the corollary and that contradiction is
     the (a) in the evidence chain.
  C. THE CHECKS ARE INJECTED — every snippet that CHANGES the model comes back
     with a MODEL CHECKS block appended to its own output, carrying the preflight
     findings (errors AND warnings) or the explicit clean line. This is the
     mechanism that does not depend on the model's cooperation.
  D. THE THREE FAILURE SHAPES — replayed as scripted edit sequences through the
     real snippet path: ground left running at the base elevation past the toe,
     the edit-cascade that strands a circle under a raised base, and a sound build
     that must come back clean. Each asserts the finding the block has to carry.
  E. THE COST GUARD — a read-only snippet gets NO block (the checks run only on a
     real edit, decided by the same signature comparison that decides whether the
     run becomes an undo step), and a snippet that raised is rolled back and gets
     none either.
  F. THE MUTATION — with the injection disabled, C and D must fail. A test that
     passes either way measures nothing.

Everything runs offscreen against the real Assistant, the real PythonKernel and
the real document; no provider is contacted and no network is touched (the agent
loop that would call one is never started — tool calls are delivered to the
handler directly, exactly as the worker thread delivers them).

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import contextlib
import io
import math
import os
import sys
import tempfile
import threading

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: One Mohr-Coulomb material, every field the schema names, so the model the
#: snippets build is complete enough for the checks to have something to say.
MAT = {'name': 'Clay', 'gamma': 125.0, 'option': 'mc', 'c': 600.0, 'phi': 20.0,
       't_cut': None, 'phi_b': None, 's_cap': None, 'u': 'none',
       'cp': 0.0, 'r_elev': 0.0, 'd': 0, 'psi': 0, 'sigma_gamma': 0.0,
       'sigma_c': 0.0, 'sigma_phi': 0.0, 'sigma_cp': 0.0, 'sigma_d': 0.0,
       'sigma_psi': 0.0, 'k1': 0.0, 'k2': 0.0, 'alpha': 0.0, 'kr0': 0.0,
       'h0': 0.0, 'E': 0.0, 'nu': 0.0}

#: The section every scripted session builds: a 15-unit slope at 2:1 whose ground
#: ENDS at the toe (0, 5) — the rigid-base shape, so raising the base to the toe
#: elevation is a legal model rather than a degenerate one. The deep base starts
#: at -20 and the circles are built against it.
TOE = (0.0, 5.0)
GROUND = [(0.0, 5.0), (30.0, 20.0), (90.0, 20.0)]
DEEP_BASE = -20.0
CENTER = (15.0, 50.0)

_R_TOE = math.hypot(CENTER[0] - TOE[0], CENTER[1] - TOE[1])

#: Snippet 1 of a build: materials only. The model is not checkable yet.
SNIPPET_MATERIALS = f"""
slope_data['materials'] = [{MAT!r}]
print('material added')
"""

#: Snippet 2: geometry, the deep base, and two circles built against it — a sound
#: model, and the one the cascade later strands.
SNIPPET_GEOMETRY = f"""
slope_data['profile_lines'] = [{{'mat_id': 0, 'coords': {GROUND!r}}}]
slope_data['max_depth'] = {DEEP_BASE!r}
slope_data['circular'] = True
slope_data['circles'] = [
    {{'Xo': {CENTER[0]!r}, 'Yo': {CENTER[1]!r}, 'Depth': {DEEP_BASE!r},
      'R': {CENTER[1] - DEEP_BASE!r}}},
    {{'Xo': {CENTER[0]!r}, 'Yo': {CENTER[1]!r}, 'Depth': {CENTER[1] - _R_TOE!r},
      'R': {_R_TOE!r}}},
]
print('geometry and circles added')
"""

#: The cascade: the base is corrected up to the toe elevation the problem states,
#: and NOTHING in the snippet touches the circles that were tangent to the old one.
SNIPPET_RAISE_BASE = """
slope_data['max_depth'] = 5.0
print('max_depth set to the stated rigid base at the toe elevation')
"""

#: The other cure the sessions reached for: keep the base where it is and run the
#: ground on at that elevation past the toe. Zero-thickness soil, degenerate domain.
SNIPPET_GROUND_PAST_TOE = """
slope_data['max_depth'] = 5.0
slope_data['profile_lines'][0]['coords'] = (
    [(-40.0, 5.0)] + list(slope_data['profile_lines'][0]['coords']))
print('ground extended left of the toe at the base elevation')
"""

SNIPPET_READ_ONLY = """
print('materials:', len(slope_data['materials']),
      'max_depth:', slope_data['max_depth'])
"""

SNIPPET_RAISES = """
slope_data['max_depth'] = 999.0
raise RuntimeError('boom')
"""


# --- harness ---------------------------------------------------------------

def _quiet(fn, *args, **kwargs):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        return fn(*args, **kwargs)


def _session():
    """A live Assistant on a real ProjectDocument, with no main window.

    The window stands in for exactly what the assistant asks of it — the
    document, a settings store, the mesh signature and two refresh hooks — so the
    checks exercise the assistant's own code and nothing of the GUI. The settings
    live in a throwaway ini file, so running this can neither read nor write the
    user's own Studio settings (or their API keys).
    """
    from PySide6.QtCore import QObject, QSettings
    from studio.ai.assistant import Assistant
    from studio.document import ProjectDocument
    from studio.main_window import MainWindow

    tmp = tempfile.NamedTemporaryFile(suffix=".ini", delete=False)
    tmp.close()

    class _Window(QObject):
        def __init__(self):
            super().__init__()
            self.doc = ProjectDocument(self)
            self.settings = QSettings(tmp.name, QSettings.IniFormat)
            self.refreshed = 0
            self.invalidated = 0

        mesh_signature = staticmethod(MainWindow.mesh_signature)

        def invalidate_results(self, clear_mesh=False):
            self.invalidated += 1

        def refresh_inputs_view(self):
            self.refreshed += 1

    mw = _Window()
    mw.doc.new()
    asst = Assistant(mw)
    asst.config.set_confirm_before_run(False)     # no modal dialog in a test
    return mw, asst


def _run(asst, code):
    """Deliver one snippet the way the worker thread delivers it, and return the
    tool result the model would receive."""
    holder = {"event": threading.Event(), "content": ""}
    _quiet(asst._on_tool_call,
           {"id": "x", "name": "run_python", "input": {"code": code},
            "holder": holder})
    if not holder["event"].is_set():
        raise AssertionError("the tool call never completed")
    return holder["content"]


def _block(text):
    """The MODEL CHECKS block out of a tool result, or None."""
    from studio.ai import assistant as A
    if A.MODEL_CHECKS_CLEAN in text:
        return A.MODEL_CHECKS_CLEAN
    if A.MODEL_CHECKS_OPEN not in text:
        return None
    start = text.index(A.MODEL_CHECKS_OPEN)
    end = text.index(A.MODEL_CHECKS_END) + len(A.MODEL_CHECKS_END)
    return text[start:end]


# --- A. the iron rules reach every provider --------------------------------

#: One anchor sentence per iron rule. Each must appear exactly once in each
#: assembled prompt — present, and present once.
IRON_ANCHORS = (
    "Max depth is the lowest elevation the problem DESCRIBES.",
    "No circle bottom below the base.",
    "On a rigid base at the toe elevation, the ground ENDS at the toe.",
    "Building is not running.",
    "MODEL CHECKS are unresolved business.",
)


def _prompts():
    """The system prompt as each tier assembles it: with the skill, and without."""
    from studio.ai.assistant import (MODELING_BRIEF, SCHEMA_BRIEF, STUDIO_SYSTEM,
                                     _load_skill_text)
    skill = _load_skill_text()
    if not skill:
        raise AssertionError("the skill body could not be loaded")
    with_skill = STUDIO_SYSTEM + skill + SCHEMA_BRIEF
    without = STUDIO_SYSTEM + SCHEMA_BRIEF + MODELING_BRIEF
    return {"skill": with_skill, "compact": without}


def check_iron_rules_once():
    out = []
    for tier, prompt in _prompts().items():
        if prompt.count("IRON RULES.") != 1:
            out.append(f"{tier}: the iron-rules block appears "
                       f"{prompt.count('IRON RULES.')} time(s), expected 1")
        for anchor in IRON_ANCHORS:
            n = prompt.count(anchor)
            if n != 1:
                out.append(f"{tier}: {anchor!r} appears {n} time(s), expected 1")
    return out


def check_assembled_prompt_is_the_real_one():
    """The prompts checked above are the prompts the two tiers actually send.

    _prompts() reassembles them from the parts; this asserts that reassembly is
    faithful by asking the live config for each tier's real system prompt and
    requiring the same iron-rule content — so a future change that stops
    including STUDIO_SYSTEM on one path fails here rather than passing on a
    string this file built for itself.
    """
    out = []
    mw, asst = _session()
    for provider, model, tier in (("anthropic", "claude-opus-5", "skill"),
                                  ("ollama", "llama3.1", "compact")):
        asst.config.set_selection(provider, model)
        want_skill = asst.config.wants_skill()
        if want_skill != (tier == "skill"):
            out.append(f"{provider}: wants_skill={want_skill}, expected "
                       f"{tier == 'skill'}")
        prompt = asst._system()
        if prompt.count("IRON RULES.") != 1:
            out.append(f"{provider}: the live system prompt carries the iron-rules "
                       f"block {prompt.count('IRON RULES.')} time(s), expected 1")
        for anchor in IRON_ANCHORS:
            if prompt.count(anchor) != 1:
                out.append(f"{provider}: live prompt carries {anchor!r} "
                           f"{prompt.count(anchor)} time(s), expected 1")
    mw.deleteLater()
    return out


def check_modeling_brief_agrees():
    """B — the compact rulebook carries the floor rule and the rigid-base
    corollary, and does not restate the iron rules a second time."""
    from studio.ai.assistant import MODELING_BRIEF
    out = []
    for phrase in ("FLOOR RULE", "Depth >= max_depth",
                   "RIGID-BASE COROLLARY", "ENDS at the toe"):
        if phrase not in MODELING_BRIEF:
            out.append(f"MODELING_BRIEF is missing {phrase!r}")
    if "never to lower max_depth" not in MODELING_BRIEF:
        out.append("MODELING_BRIEF does not rule out lowering max_depth as the cure")
    for anchor in IRON_ANCHORS:
        if anchor in MODELING_BRIEF:
            out.append(f"MODELING_BRIEF duplicates the iron rule {anchor!r}")
    return out


# --- C/D. the checks are injected, on the three failure shapes --------------

def check_clean_build():
    """A sound build, in two snippets: the incomplete model says so, the finished
    one comes back clean."""
    from studio.ai import assistant as A
    out = []
    mw, asst = _session()

    first = _run(asst, SNIPPET_MATERIALS)
    blk = _block(first)
    if blk is None:
        out.append("materials-only snippet: no MODEL CHECKS block at all")
    elif "not run" not in blk or "geometry" not in blk:
        out.append(f"materials-only snippet: expected the incomplete-model note, "
                   f"got {blk[:120]!r}")

    second = _run(asst, SNIPPET_GEOMETRY)
    if _block(second) != A.MODEL_CHECKS_CLEAN:
        out.append(f"finished build: expected the clean line, got "
                   f"{_block(second)!r}")
    if not mw.doc.slope_data.get("domain_polygon"):
        out.append("finished build: the derived geometry was not rebuilt")
    mw.deleteLater()
    return out


def check_edit_cascade():
    """The cascade: a snippet that touches ONLY max_depth must surface the circles
    it just stranded under the new base — the finding no one asked for, on the
    edit that caused it."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)

    text = _run(asst, SNIPPET_RAISE_BASE)
    blk = _block(text)
    if blk is None:
        out.append("raising max_depth produced no MODEL CHECKS block")
        mw.deleteLater()
        return out
    if "surface.circle_below_domain_floor" not in blk:
        out.append(f"the stranded circle was not reported: {blk[:200]!r}")
    if "ERROR" not in blk:
        out.append(f"the stranded first circle was not an ERROR: {blk[:200]!r}")
    if "Circle 1" not in blk:
        out.append(f"the finding does not name the circle: {blk[:200]!r}")
    # The snippet's own output is still there — the block is appended, not a
    # replacement for the result the model asked for.
    if "max_depth set to the stated rigid base" not in text:
        out.append("the block displaced the snippet's own stdout")
    if not text.rstrip().endswith("=== END MODEL CHECKS ==="):
        out.append("the block is not the last thing in the tool result")
    mw.deleteLater()
    return out


def check_ground_past_the_toe():
    """The other cure the sessions reached for — ground carried on at the base
    elevation past the toe — must come back as the degenerate domain it is."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)

    blk = _block(_run(asst, SNIPPET_GROUND_PAST_TOE))
    if blk is None:
        out.append("extending the ground past the toe produced no block")
    else:
        if "domain.degenerate_ring" not in blk:
            out.append(f"the degenerate domain was not reported: {blk[:200]!r}")
        if "ERROR" not in blk:
            out.append(f"the degenerate domain was not an ERROR: {blk[:200]!r}")
    mw.deleteLater()
    return out


def check_surfaceless_models():
    """A limit-equilibrium model with no failure surface yet is reported as one;
    a finite-element model, whose failure surface is an OUTPUT, is not.

    The exclusion that separates them names a rule by id, and a misspelled id
    excludes nothing at all while looking exactly like this — so both halves are
    asserted, and the id is checked against the live registry.
    """
    from xslope.preflight import rules
    out = []
    known = {r.id for r in rules()}
    if "surface.none_defined" not in known:
        out.append("surface.none_defined is not a rule id — the exclusion is dead")

    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    blk = _block(_run(asst, f"""
slope_data['profile_lines'] = [{{'mat_id': 0, 'coords': {GROUND!r}}}]
slope_data['max_depth'] = {DEEP_BASE!r}
print('geometry only')
"""))
    if blk is None or "surface.none_defined" not in blk:
        out.append(f"a model with no failure surface was not reported: {blk!r}")

    # Attach a mesh: the same model is now a finite-element model in progress.
    # The mesh is not one of the tracked source inputs, so the snippet also edits
    # one that is -- otherwise nothing would be checked and the assertion below
    # would pass on an absent block rather than on a correct one.
    blk = _block(_run(asst, """
slope_data['mesh'] = {'nodes': [], 'elements': [], 'element_types': [],
                      'element_materials': []}
slope_data['materials'][0]['E'] = 500000.0
print('mesh attached')
"""))
    if blk is None:
        out.append("the finite-element leg produced no block to assert against")
    elif "surface.none_defined" in blk:
        out.append("a finite-element model was told it has no failure surface")
    mw.deleteLater()
    return out


def check_staged_by_run_is_annotated():
    """A finding an EDIT cannot answer must not arrive as a bare error to fix.

    A material reading u = seep with no solved field is the case: iron rule 4
    forbids the assistant running the seepage analysis unasked, so a bare ERROR
    leaves one apparent way out — change the pore-pressure option, which does not
    supply the field, it deletes the requirement and analyses a different problem.
    That is failure shape (b) rebuilt inside the guardrail, so the block has to
    name these for what they are and route them to the run.
    """
    from xslope.preflight import STAGED_BY_RUN
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)

    blk = _block(_run(asst, """
slope_data['materials'][0]['u'] = 'seep'
print('material now reads a seepage field')
"""))
    if blk is None:
        out.append("the u = seep edit produced no block")
        mw.deleteLater()
        return out
    if "seep_field.missing" not in blk:
        out.append(f"the missing seepage field was not reported at all: {blk[:200]!r}")
    if "STAGED BY A RUN" not in blk:
        out.append(f"a staged finding arrived unannotated: {blk[:300]!r}")
    if "resolved by running the seepage analysis" not in blk:
        out.append(f"the block does not name the run that resolves it: {blk[:300]!r}")
    if "offer to run it" not in blk:
        out.append("the block does not route the model to offering the run")
    # It must NOT be counted among the problems to fix by editing.
    if "found 1 problem(s)" in blk and "seep_field.missing" in blk.split(
            "STAGED BY A RUN")[0]:
        out.append("the staged finding was also listed as a problem to fix")
    mw.deleteLater()

    # The declaration must cover every rule whose cure is an upstream run. The
    # criterion is the message naming one, so the module's own source is read for
    # that phrasing: a sibling added later cannot slip past as a bare error.
    import re
    src = open(os.path.join(_REPO, "xslope", "preflight.py"), encoding="utf-8").read()
    for block in re.split(r'\n@rule\(', src)[1:]:
        rid = re.match(r'"([^"]+)"', block)
        body = block.split("\n@rule(")[0]
        if rid and re.search(r'"?[Rr]e-?[Rr]un the \w+ analysis', body):
            if rid.group(1) not in STAGED_BY_RUN:
                out.append(f"{rid.group(1)}: its cure is an upstream run but it is "
                           f"not declared in preflight.STAGED_BY_RUN")
    return out


def check_surface_family_flip_is_an_edit():
    """Flipping the stored surface family changes which surface the run analyses,
    so it is an edit: an undo step, a dirty document, and a checks block. It was
    none of the three, because `circular` was not among the tracked source keys."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    steps = len(mw.doc.undo_labels())

    blk = _block(_run(asst, """
slope_data['non_circ'] = [{'X': -10.0, 'Y': 5.0, 'Movement': 'Free'},
                          {'X': 20.0, 'Y': -2.0, 'Movement': 'Free'},
                          {'X': 50.0, 'Y': 20.0, 'Movement': 'Free'}]
slope_data['circular'] = False
print('switched to the non-circular surface')
"""))
    if blk is None:
        out.append("adding a surface and flipping the family produced no block")
    # Now flip it BACK, alone: nothing else in the model changes.
    before = len(mw.doc.undo_labels())
    blk = _block(_run(asst, """
slope_data['circular'] = True
print('switched back to the circular surface')
"""))
    if blk is None:
        out.append("flipping the surface family alone produced no checks block")
    if len(mw.doc.undo_labels()) != before + 1:
        out.append("flipping the surface family alone left no undo step")
    elif not mw.doc.undo_labels()[0].startswith("Assistant:"):
        out.append(f"the undo step is unlabelled: {mw.doc.undo_labels()[0]!r}")
    elif "surface family" not in mw.doc.undo_labels()[0]:
        out.append(f"the undo step does not name what changed: "
                   f"{mw.doc.undo_labels()[0]!r}")
    if not mw.doc.dirty:
        out.append("flipping the surface family left the document clean")
    if len(mw.doc.undo_labels()) <= steps:
        out.append("the surface-family edits produced no history at all")
    mw.deleteLater()
    return out


def check_followup_command_carries_the_selection():
    """The '+N more' command the block hands over must reproduce the block's own
    findings — so it carries the same selection. Without it the ambiguity warning
    the block suppressed comes back, and the model is told about a finding the
    block never made."""
    from studio.ai import assistant as A
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    # A deck carrying BOTH families (so the selection is what keeps
    # surface.family_ambiguous off the block) and one real fault (so there is
    # something for the overflow line to be counting).
    _run(asst, SNIPPET_RAISE_BASE)
    _run(asst, """
slope_data['non_circ'] = [{'X': -10.0, 'Y': 5.0, 'Movement': 'Free'},
                          {'X': 20.0, 'Y': -2.0, 'Movement': 'Free'},
                          {'X': 50.0, 'Y': 20.0, 'Movement': 'Free'}]
print('both families present')
""")
    sd = mw.doc.slope_data
    sel = A._checks_selection(sd)
    if not sel:
        out.append("a both-family deck produced no surface selection")
    from xslope.preflight import preflight
    with_sel = {f.rule_id for f in preflight(sd, "lem", sel).findings}
    without = {f.rule_id for f in preflight(sd, "lem").findings}
    if with_sel == without:
        out.append("the selection changes nothing here — the leg proves nothing")
    elif "surface.family_ambiguous" not in (without - with_sel):
        out.append(f"unexpected difference: {without - with_sel}")

    # Force the overflow line and read the command out of it.
    real_max = A.MAX_CHECK_FINDINGS
    try:
        A.MAX_CHECK_FINDINGS = 0
        text = A.model_checks_text(sd)
    finally:
        A.MAX_CHECK_FINDINGS = real_max
    line = [ln for ln in text.splitlines() if "more —" in ln]
    if not line:
        out.append("no overflow line was produced to check")
    elif repr(sel) not in line[0]:
        out.append(f"the follow-up command drops the selection: {line[0]!r}")
    mw.deleteLater()
    return out


def check_warning_is_reported_too():
    """A WARNING is unresolved business as much as an ERROR — the sessions failed
    on warnings passed over in silence, so the block must carry them."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    # A piezometric line entered right to left: WARNING, nothing more.
    blk = _block(_run(asst, """
slope_data['piezo_line'] = [(90.0, 12.0), (30.0, 12.0), (0.0, 8.0)]
slope_data['materials'][0]['u'] = 'piezo'
print('piezo line added')
"""))
    if blk is None or "WARNING" not in (blk or ""):
        out.append(f"a warning-only model reported no warning: {blk!r}")
    elif "order.piezo_reversed" not in blk:
        out.append(f"the reversed piezo line was not reported: {blk[:200]!r}")
    mw.deleteLater()
    return out


# --- E. the cost guard -----------------------------------------------------

def check_read_only_costs_nothing():
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)

    text = _run(asst, SNIPPET_READ_ONLY)
    if _block(text) is not None:
        out.append("a read-only snippet ran the checks anyway")
    if "materials: 1" not in text:
        out.append(f"the read-only snippet lost its own output: {text[:120]!r}")

    before = dict(mw.doc.slope_data)
    failed = _run(asst, SNIPPET_RAISES)
    if _block(failed) is not None:
        out.append("a snippet that raised ran the checks on a rolled-back model")
    if "ERROR" not in failed:
        out.append("a snippet that raised did not report its traceback")
    if mw.doc.slope_data.get("max_depth") != before.get("max_depth"):
        out.append("the failed snippet's edit was not rolled back")
    mw.deleteLater()
    return out


# --- F. the mutation -------------------------------------------------------

def check_mutation_disables_the_checks():
    """Disable the injection and the three shapes must stop being detected. A
    check that passes with the mechanism removed is measuring nothing."""
    from studio.ai import assistant as A
    out = []
    original = A.model_checks_text
    A.model_checks_text = lambda sd: ""
    try:
        survivors = [name for name, fn in (
            ("clean build", check_clean_build),
            ("edit cascade", check_edit_cascade),
            ("ground past the toe", check_ground_past_the_toe),
            ("warning reported", check_warning_is_reported_too),
            ("staged by a run", check_staged_by_run_is_annotated),
            ("surface-family flip", check_surface_family_flip_is_an_edit),
        ) if not fn()]
    finally:
        A.model_checks_text = original
    for name in survivors:
        out.append(f"mutation: '{name}' still passed with the injection disabled")
    # ... and the guards must come back green once it is restored.
    if check_clean_build():
        out.append("mutation: the injection did not come back")
    return out


CHECKS = [
    ("A. iron rules, once per prompt tier", check_iron_rules_once),
    ("A. the live prompts carry them", check_assembled_prompt_is_the_real_one),
    ("B. modeling brief agrees", check_modeling_brief_agrees),
    ("C/D. sound build comes back clean", check_clean_build),
    ("D. the edit cascade surfaces", check_edit_cascade),
    ("D. ground past the toe surfaces", check_ground_past_the_toe),
    ("D. warnings are reported too", check_warning_is_reported_too),
    ("D. surface-less LEM vs FEM models", check_surfaceless_models),
    ("G. staged-by-a-run findings annotated", check_staged_by_run_is_annotated),
    ("G. a surface-family flip is an edit", check_surface_family_flip_is_an_edit),
    ("G. the follow-up carries the selection",
     check_followup_command_carries_the_selection),
    ("E. read-only and failed runs cost nothing", check_read_only_costs_nothing),
    ("F. mutation: injection disabled", check_mutation_disables_the_checks),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("assistant guardrails: PySide6 not installed — checks skipped.")
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
    print("assistant guardrails (Studio):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll assistant guardrail checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
