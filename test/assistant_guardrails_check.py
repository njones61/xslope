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
  H. THE DELTA — a finding is quoted in full once. On later edits the block quotes
     what is new or changed and names the rest by rule key on one line, so a build
     on a deck with standing faults does not re-read the same paragraphs after
     every edit. The honesty of that line is what these legs are for: it says
     "reported in full earlier", so a finding may collapse only once it has been
     QUOTED — never merely counted on the overflow line the quote cap produces —
     and two faults that differ only in a digit belonging to a field name (Lp1
     against Lp2) must not read as one finding. What else must survive: an
     edit-answerable ERROR is quoted every time, a changed finding un-collapses,
     the staged-by-a-run label survives the collapse, nothing is dropped without
     being named, and a new session starts from silence. Mutations for each.
  L. A RUN IS THE SESSION'S RUN — a recorded session was told its report covered
     a run the report did not contain, because `run_lem` computed an answer and
     stored it nowhere. So: the bundle lands where `report_solutions` reads it, a
     report built afterwards carries the run and a traceability stamp that names
     the input file, an unqualified run runs the method the MODEL declares rather
     than a hardcoded one, the result says which surface it is the answer for, a
     geometry edit made on the source the resync overwrites is named rather than
     silently reverted, and the transcript renders the assistant's markdown
     instead of printing its tables as rows of pipes — math included, which the
     dialect has no notion of: the brief says write it as plain text, and LaTeX
     that arrives anyway is converted rather than shown as its own source.
  M. WHAT GETS RUN, AND WHERE A FILE LANDS — code is executed because the model
     ASKED for it to be, never because it happens to compile: a provider with tool
     calling has a protocol for asking and so is never read for code in prose, and
     the fallback the tool-less local models keep runs only a fence marked
     ```python run (the convention their prompt, and only their prompt, states).
     A signature written to show the user the shape of a call was executed once,
     and the NameError explained back to the user as their own. And a file the
     assistant says can be opened has to be openable: a report generated with no
     path lands in the output folder the dock's Files button opens, and is named
     among the files the snippet reports.
  N. TRANSIENT SEEPAGE, AND THE SNIPPET THAT NEVER RETURNS — a corpus sweep met a
     transient model with no helper for it, and the assistant spent the turn
     reverse-engineering engine signatures instead of answering; another sweep sat
     forty minutes on a snippet that never returned, because the snippet runs on
     the GUI thread and so does the timer that was supposed to end the turn. So:
     the transient march is a helper, it marches a real sample into stored frames
     and hands the second caller the march it already has, the route to one instant
     is in the brief the model reads, and a snippet that overruns its limit is
     stopped, rolled back and reported to the model as its own result — with the
     limit off, the same snippet finishes.
  I. THE THINGS THE HARNESS OWNS WHATEVER THE MODEL DOES — every helper the
     prompt tells the model to call is in the namespace with the arguments the
     prompt names (a helper described and absent is a NameError followed by the
     hand-rolled pipeline the helpers exist to prevent); every provider and every
     model the settings dialog offers can accept an image, since the headline
     request is a cross section handed over as a picture; the Kimi key is stored
     under the name the ruling fixed; and the token meter reads each completion's
     usage, accumulates it per turn and per session, and survives a response that
     carries no usage at all.

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
#: The deep circle of the sound build: it bottoms below the toe (so raising the
#: base to the toe elevation strands it) and still daylights on the slope and on
#: the crest. A circle reaching the -20 base from CENTER would swallow the toe
#: and produce no failure surface at all.
DEEP_CENTER = (30.0, 50.0)
DEEP_DEPTH = 0.0

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
    {{'Xo': {DEEP_CENTER[0]!r}, 'Yo': {DEEP_CENTER[1]!r}, 'Depth': {DEEP_DEPTH!r},
      'R': {DEEP_CENTER[1] - DEEP_DEPTH!r}}},
    {{'Xo': {CENTER[0]!r}, 'Yo': {CENTER[1]!r}, 'Depth': {CENTER[1] - _R_TOE!r},
      'R': {_R_TOE!r}}},
]
print('geometry and circles added')
"""

#: The cascade: the base is corrected up to the toe elevation the problem states,
#: and NOTHING in the snippet touches the circles that were tangent to the old one.
#: The cascade edit: the crest is raised to el. 60 behind the slope, so both
#: circles now meet the ground ABOVE their own centers and neither can be built
#: as a failure surface. (Raising the base instead is not a cascade at all -- a
#: circle whose nadir falls below a raised floor is truncated along it and still
#: slices.)
SNIPPET_RAISE_CREST = """
slope_data['profile_lines'] = [{'mat_id': 0, 'coords': [(0.0, 5.0), (30.0, 20.0), (60.0, 60.0), (90.0, 60.0)]}]
print('crest raised to el. 60 behind the slope')
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

#: A piezometric line entered right to left, and a material that reads it: two
#: standing warnings, the fault a later edit leaves untouched.
SNIPPET_PIEZO_REVERSED = """
slope_data['piezo_line'] = [(90.0, 12.0), (30.0, 12.0), (0.0, 8.0)]
slope_data['materials'][0]['u'] = 'piezo'
print('piezo line added')
"""

#: An edit that touches a tracked input and changes no finding — the shape of the
#: turn the delta exists for.
SNIPPET_NEUTRAL_EDIT = """
slope_data['materials'][0]['gamma'] = 126.0
print('unit weight nudged')
"""

#: A third circle, far below the base: a NEW finding, arriving on a model that
#: already carries standing ones.
SNIPPET_ADD_DEAD_CIRCLE = """
slope_data['circles'].append({'Xo': 15.0, 'Yo': 50.0, 'Depth': -300.0, 'R': 350.0})
print('added a circle far below the base')
"""

#: The same fault, on a different row: circle 3 is repaired and circle 2 is broken
#: in its place. Same rule, same wording, same numbers — only the row differs, and
#: that makes it a different finding.
SNIPPET_MOVE_DEAD_CIRCLE = """
slope_data['circles'][2] = {'Xo': 15.0, 'Yo': 50.0, 'Depth': 6.0, 'R': 44.0}
slope_data['circles'][1] = {'Xo': 15.0, 'Yo': 50.0, 'Depth': -300.0, 'R': 350.0}
print('moved the dead circle from row 3 to row 2')
"""

#: Eight more circles, every one of them under the base: more findings than one
#: block quotes, which is what makes the quote cap decide anything.
SNIPPET_MANY_DEAD_CIRCLES = """
for i in range(8):
    slope_data['circles'].append(
        {'Xo': 15.0, 'Yo': 50.0, 'Depth': -100.0 - 10.0 * i, 'R': 150.0 + 10.0 * i})
print('added 8 circles below the base')
"""

#: A material that reads a seepage field there is no run to supply: the staged
#: finding.
SNIPPET_SEEP_MATERIAL = """
slope_data['materials'][0]['u'] = 'seep'
print('material now reads a seepage field')
"""

#: Drop the deep base circle and keep the toe circle, so the surface the model
#: DEFINES is one a single-surface run can solve.
SNIPPET_ONE_CIRCLE = """
slope_data['circles'] = [slope_data['circles'][1]]
print('kept the toe circle')
"""

#: The session-2 failure, replayed: the face is laid back by rebuilding
#: `polygons` on a model whose geometry source is `profile_lines`. The resync
#: rebuilds them from the profile lines and the edit is gone.
SNIPPET_POLYGON_EDIT_ON_PROFILE_MODEL = """
slope_data['polygons'] = [{'coords': [(0.0, 5.0), (45.0, 20.0), (90.0, 20.0),
                                      (90.0, -20.0), (0.0, -20.0)], 'mat_id': 0}]
print('polygons rebuilt for a 3:1 face')
"""

#: The same discarded edit, with the READ-BACK that made it dangerous: the
#: snippet swaps the zone's material and prints the polygons back, and its print
#: runs BEFORE the rebuild, so its own output shows the edit as applied. Two
#: benchmark sessions printed exactly this and reported the model repaired.
SNIPPET_POLYGON_EDIT_READ_BACK = """
slope_data['polygons'] = [{'coords': [(0.0, 5.0), (45.0, 20.0), (90.0, 20.0),
                                      (90.0, -20.0), (0.0, -20.0)], 'mat_id': 1}]
for i, p in enumerate(slope_data['polygons']):
    print('polygon', i, 'mat_id', p['mat_id'], 'coords', p['coords'])
print('READ BACK: the zone is now material', slope_data['polygons'][0]['mat_id'])
"""

#: The same section built the OTHER way — polygons and no profile lines — and
#: then edited on its own source, which must not be warned about.
SNIPPET_POLYGON_NATIVE_BUILD = """
slope_data['profile_lines'] = []
slope_data['max_depth'] = -20.0
slope_data['polygons'] = [{'coords': [(0.0, 5.0), (30.0, 20.0), (90.0, 20.0),
                                      (90.0, -20.0), (0.0, -20.0)], 'mat_id': 0}]
slope_data['circular'] = True
slope_data['circles'] = [{'Xo': 15.0, 'Yo': 50.0, 'Depth': 2.57, 'R': 47.43}]
print('polygon-native model built')
"""

SNIPPET_POLYGON_NATIVE_EDIT = """
slope_data['polygons'][0]['coords'] = [(0.0, 5.0), (45.0, 20.0), (90.0, 20.0),
                                       (90.0, -20.0), (0.0, -20.0)]
print('face laid back on the polygon the model is built from')
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
            # What MainWindow.report_solutions labels the LEM run with. Empty
            # until something runs, exactly as it is in a fresh window.
            self._last_lem_opts = {}

        mesh_signature = staticmethod(MainWindow.mesh_signature)

        def invalidate_results(self, clear_mesh=False):
            self.invalidated += 1
            # What the real sweep does to the document: every cached solution
            # goes, because the inputs it was solved on have changed. Counting
            # the calls is not enough — a run made AFTER the edit that triggers
            # this has to survive it, and that is only visible if the keys
            # actually go.
            for key in ("lem_solution", "seep_solutions", "transient_seep",
                        "fem_solution", "design", "sensitivity", "fs_vs_time"):
                self.doc.results.pop(key, None)

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
                                     _load_brief_text)
    brief = _load_brief_text()
    if not brief:
        raise AssertionError("the Studio brief could not be loaded")
    with_brief = STUDIO_SYSTEM + brief + SCHEMA_BRIEF
    without = STUDIO_SYSTEM + SCHEMA_BRIEF + MODELING_BRIEF
    return {"brief": with_brief, "compact": without}


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
    for provider, model, tier in (("anthropic", "claude-opus-5", "brief"),
                                  ("ollama", "llama3.1", "compact")):
        asst.config.set_selection(provider, model)
        want_skill = asst.config.wants_skill()
        if want_skill != (tier == "brief"):
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
    """The cascade: a snippet that touches ONLY the profile must surface the circles
    it just stranded under the new base — the finding no one asked for, on the
    edit that caused it."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)

    text = _run(asst, SNIPPET_RAISE_CREST)
    blk = _block(text)
    if blk is None:
        out.append("raising the crest produced no MODEL CHECKS block")
        mw.deleteLater()
        return out
    if "surface.circle_daylights_above_center" not in blk:
        out.append(f"the stranded circle was not reported: {blk[:200]!r}")
    if "ERROR" not in blk:
        out.append(f"the stranded first circle was not an ERROR: {blk[:200]!r}")
    if "Circle 1" not in blk:
        out.append(f"the finding does not name the circle: {blk[:200]!r}")
    # The snippet's own output is still there — the block is appended, not a
    # replacement for the result the model asked for.
    if "crest raised to el. 60" not in text:
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
    _run(asst, SNIPPET_RAISE_CREST)
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
    line = [ln for ln in text.splitlines() if "not quoted here" in ln]
    if not line:
        out.append("no overflow line was produced to check")
    elif "print(preflight(" not in line[0]:
        out.append(f"the overflow line lost the command that reads them: {line[0]!r}")
    elif repr(sel) not in line[0]:
        out.append(f"the follow-up command drops the selection: {line[0]!r}")
    mw.deleteLater()
    return out


def check_secondary_circle_wording():
    """What a dead circle costs is what the run that reads it loses — so the
    finding must not tell a single-surface run that its "other circles" are fine.
    That run reads circles[0] and nothing else; the circle reported is one it
    never looks at, and the loss belongs to a search seeded from the same sheet.
    """
    import copy
    from xslope.preflight import preflight
    from studio.editors import _resync_geometry
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    _run(asst, """
slope_data['circles'].append({'Xo': 15.0, 'Yo': 50.0, 'Depth': -300.0, 'R': 350.0})
print('added a circle far below the base')
""")
    sd = copy.deepcopy(mw.doc.slope_data)
    _resync_geometry(sd)
    rid = "surface.circle_below_domain_floor"

    single = [f for f in preflight(sd, "lem", {"surface": "circular"}).findings
              if f.rule_id == rid]
    search = [f for f in preflight(
        sd, "lem", {"surface": "circular", "search": True}).findings
        if f.rule_id == rid]
    if len(single) != 1 or len(search) != 1:
        out.append(f"expected one finding each; got {len(single)} / {len(search)}")
        mw.deleteLater()
        return out
    if single[0].severity != "warning":
        out.append(f"a secondary dead circle is not a warning: {single[0].severity}")
    if "other circles are unaffected" in single[0].message:
        out.append("a single-surface run was told its other circles are unaffected")
    if "does not read this one" not in single[0].message:
        out.append(f"the single-surface finding does not say the run skips it: "
                   f"{single[0].message[:160]!r}")
    if "a search seeded from this sheet would lose it" not in single[0].message:
        out.append("the single-surface finding does not say who does lose it")
    if "search's other circles are unaffected" not in search[0].message:
        out.append(f"the search finding lost its own wording: "
                   f"{search[0].message[:160]!r}")
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


# --- H. the delta ----------------------------------------------------------

#: A fragment of each standing finding's own paragraph. Present = quoted in full.
PIEZO_FULL = "the points run right to left"
CIRCLE_FULL = "above its own center"          # the cascade ERROR (crest raised)
DEAD_CIRCLE_FULL = "below the bottom of the model domain"   # a dead circle added by hand
SEEP_FULL = "takes pore pressure from a seepage solution"


def _collapsed_lines(blk):
    """Every one-line summary of findings the block reported earlier."""
    return [ln for ln in (blk or "").splitlines() if "earlier finding" in ln]


def _collapsed(blk):
    """The block's one-line summary of findings already reported, or None."""
    lines = _collapsed_lines(blk)
    return lines[0] if lines else None


def check_delta_identity_is_the_row():
    """What makes two findings the same finding.

    Keyed on the message text, a fault whose quoted values drift under an
    unrelated edit reads as new and is quoted again — which is the cost this
    feature exists to remove. Keyed on the rule alone, the same fault moving to a
    different row reads as old and is never quoted at all. So identity is the rule
    plus the row, and a change is the severity plus the message with its numbers
    blanked.
    """
    from xslope.preflight import Finding
    from studio.ai.assistant import _finding_key, _finding_sig, _finding_subject
    out = []
    rid = "surface.circle_below_domain_floor"

    def f(msg, severity="warning", rule=rid):
        return Finding(rule_id=rule, severity=severity, message=msg)

    a = f("Circle 3 has Depth = -300, below the bottom of the model domain at y = -20.")
    b = f("Circle 3 has Depth = -412.5, below the bottom of the model domain at y = -20.")
    c = f("Circle 2 has Depth = -300, below the bottom of the model domain at y = -20.")
    if _finding_key(a) != _finding_key(b) or _finding_sig(a) != _finding_sig(b):
        out.append("a finding whose quoted values moved reads as a different one")
    if _finding_key(a) == _finding_key(c):
        out.append("the same fault on a different circle reads as the same finding")
    if _finding_sig(a) == _finding_sig(f(a.message, severity="error")):
        out.append("a severity change does not register as a change")
    if _finding_key(a) == _finding_key(f(a.message, rule="other.rule")):
        out.append("two rules about one row read as the same finding")

    # A row label is a row label, not any number a sentence happens to contain.
    if _finding_subject("Material 3 ('Core') has no unit weight.") != "Material 3 ('Core')":
        out.append("a named material row is not recognised as the row")
    for msg in ("None of the 5 reinforcement line(s) crosses any failure surface.",
                "This model defines both a circular and a non-circular surface.",
                "5 reinforcement line(s) leave Type blank, starting at line 1."):
        if _finding_subject(msg) is not None:
            out.append(f"a row was read out of prose that names none: "
                       f"{_finding_subject(msg)!r}")
    # ... and with no row to key on, two findings of one rule stay distinct.
    m1 = f("Rules read left to right.", rule="x.y")
    m2 = f("Rules read right to left.", rule="x.y")
    if _finding_key(m1) == _finding_key(m2):
        out.append("two unlabelled findings of one rule collapse into one key")
    return out


def check_new_finding_is_reported_in_full():
    """A finding the session has not seen arrives in full, whatever else stands."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    _run(asst, SNIPPET_PIEZO_REVERSED)
    _run(asst, SNIPPET_NEUTRAL_EDIT)

    blk = _block(_run(asst, SNIPPET_ADD_DEAD_CIRCLE))
    if blk is None:
        out.append("adding a dead circle produced no block")
        mw.deleteLater()
        return out
    if DEAD_CIRCLE_FULL not in blk:
        out.append(f"the new finding was not quoted in full: {blk[:200]!r}")
    if "surface.circle_below_domain_floor" not in blk:
        out.append("the new finding is not named by its rule key either")
    # The standing ones are still accounted for, on one line.
    line = _collapsed(blk)
    if line is None or "order.piezo_reversed" not in line:
        out.append(f"the standing findings were not carried: {blk[:200]!r}")
    if PIEZO_FULL in blk:
        out.append("a standing finding was quoted in full alongside the new one")
    mw.deleteLater()
    return out


def check_unchanged_finding_collapses():
    """The second report of an unchanged warning is one line naming its rule."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    first = _block(_run(asst, SNIPPET_PIEZO_REVERSED))
    if first is None or PIEZO_FULL not in first:
        out.append(f"the first report did not quote the finding in full: {first!r}")

    second = _block(_run(asst, SNIPPET_NEUTRAL_EDIT))
    if second is None:
        out.append("the next edit produced no block")
        mw.deleteLater()
        return out
    line = _collapsed(second)
    if line is None:
        out.append(f"the unchanged findings did not collapse: {second!r}")
    else:
        if "order.piezo_reversed" not in line:
            out.append(f"the collapsed line does not name the rule: {line!r}")
        if "2 earlier findings" not in line:
            out.append(f"the collapsed line does not count them: {line!r}")
        if "unresolved" not in line:
            out.append(f"the collapsed line does not say they still stand: {line!r}")
    if PIEZO_FULL in second:
        out.append("an unchanged finding was quoted in full a second time")
    if len(second) >= len(first or ""):
        out.append("the second report is no shorter than the first")
    mw.deleteLater()
    return out


def check_error_never_collapses():
    """A run-blocking ERROR refuses the analysis until it is answered, so it is
    live business on every edit — quoted in full in every block, however many
    edits it survives."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    blk = _block(_run(asst, SNIPPET_RAISE_CREST))
    if blk is None or "ERROR" not in blk:
        out.append(f"the stranded circle did not arrive as an error: {blk!r}")

    for i, code in enumerate((SNIPPET_NEUTRAL_EDIT, SNIPPET_PIEZO_REVERSED,
                              SNIPPET_NEUTRAL_EDIT.replace("126.0", "127.0"))):
        blk = _block(_run(asst, code))
        if blk is None:
            out.append(f"edit {i + 2} produced no block")
            continue
        if "ERROR [surface.circle_daylights_above_center]" not in blk:
            out.append(f"edit {i + 2}: the error stopped being quoted: {blk[:200]!r}")
        elif CIRCLE_FULL not in blk:
            out.append(f"edit {i + 2}: the error lost its message: {blk[:200]!r}")
    mw.deleteLater()
    return out


def check_error_survives_the_quote_cap():
    """The error exemption outranks the rotation.

    Both rules are right on their own and they compete for the same six slots:
    the rotation puts never-quoted findings first so nothing is counted as told
    without being read, and an error already quoted has been read. Order by that
    alone and a batch of new warnings pushes the error off the quoted list
    entirely — the block the model reads carries no error line while the run is
    still refused. So the error takes the quota first and the rotation orders
    what is left.
    """
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    blk = _block(_run(asst, SNIPPET_RAISE_CREST))
    if blk is None or "ERROR [surface.circle_daylights_above_center]" not in blk:
        out.append(f"the stranded circle did not arrive as an error: {blk!r}")

    # Now bury it: more never-quoted warnings arrive at once than a block quotes.
    blocks = [_block(_run(asst, SNIPPET_MANY_DEAD_CIRCLES))]
    blocks.append(_block(_run(asst, SNIPPET_NEUTRAL_EDIT)))
    if "not quoted here" not in (blocks[0] or ""):
        out.append("the warnings did not overflow the quote cap — the leg proves "
                   "nothing")
    for n, blk in enumerate(blocks, 2):
        quoted = [ln for ln in (blk or "").splitlines()
                  if ln.startswith("  ERROR [surface.circle_daylights_above_center]")]
        if not quoted:
            out.append(f"block {n}: the error is not quoted at all: "
                       f"{(blk or '')[:300]!r}")
        elif CIRCLE_FULL not in quoted[0]:
            out.append(f"block {n}: the error line lost its message: "
                       f"{quoted[0][:120]!r}")
    mw.deleteLater()
    return out


def check_changed_finding_uncollapses():
    """The same rule, the same wording, the same numbers — on a different row. It
    is a different fault about a different circle, and it is quoted in full."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    first = _block(_run(asst, SNIPPET_ADD_DEAD_CIRCLE))
    if first is None or "Circle 3" not in first:
        out.append(f"the dead circle was not reported on row 3: {first!r}")

    blk = _block(_run(asst, SNIPPET_MOVE_DEAD_CIRCLE))
    if blk is None:
        out.append("moving the dead circle produced no block")
        mw.deleteLater()
        return out
    if "Circle 2" not in blk or DEAD_CIRCLE_FULL not in blk:
        out.append(f"the fault's new row was not quoted in full: {blk[:200]!r}")
    if _collapsed(blk) is not None:
        out.append(f"the moved fault was written off as unchanged: "
                   f"{_collapsed(blk)!r}")
    if "Circle 3" in blk:
        out.append("the repaired row is still being reported")
    mw.deleteLater()
    return out


def check_fresh_session_reports_in_full():
    """Tracking belongs to the session. A new conversation, and a second
    assistant on the same document, both start from silence — so a model that
    never saw the first report is not handed a rule key and nothing else."""
    from studio.ai.assistant import Assistant
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    _run(asst, SNIPPET_PIEZO_REVERSED)
    if _collapsed(_block(_run(asst, SNIPPET_NEUTRAL_EDIT))) is None:
        out.append("the finding never collapsed, so this leg proves nothing")

    asst.reset()                       # New chat: nothing has been said yet.
    blk = _block(_run(asst, SNIPPET_NEUTRAL_EDIT.replace("126.0", "127.0")))
    if blk is None or PIEZO_FULL not in blk:
        out.append(f"a fresh conversation was not given the findings in full: "
                   f"{blk!r}")
    if _collapsed(blk) is not None:
        out.append(f"a fresh conversation still collapsed a finding: "
                   f"{_collapsed(blk)!r}")

    other = Assistant(mw)              # A second session over the same document.
    other.config.set_confirm_before_run(False)
    blk = _block(_run(other, SNIPPET_NEUTRAL_EDIT.replace("126.0", "128.0")))
    if blk is None or PIEZO_FULL not in blk:
        out.append(f"a second session inherited the first one's report: {blk!r}")
    mw.deleteLater()
    return out


def check_staged_collapse_keeps_the_label():
    """A staged finding collapses like a warning — no edit answers it, so it is
    not edit-blocking — but the collapsed line still says it is staged by a run.
    A line that stopped saying so would read as one more edit outstanding, which
    is the mistake the annotation exists to prevent."""
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    first = _block(_run(asst, SNIPPET_SEEP_MATERIAL))
    if first is None or SEEP_FULL not in first:
        out.append(f"the staged finding was not quoted in full first: {first!r}")

    blk = _block(_run(asst, SNIPPET_NEUTRAL_EDIT))
    if blk is None:
        out.append("the next edit produced no block")
        mw.deleteLater()
        return out
    line = _collapsed(blk)
    if line is None:
        out.append(f"the staged finding did not collapse: {blk!r}")
    else:
        if "STAGED BY A RUN" not in line:
            out.append(f"the collapsed line dropped the staged label: {line!r}")
        if "not by an edit" not in line:
            out.append(f"the collapsed line no longer says an edit cannot answer "
                       f"it: {line!r}")
        if "seep_field.missing" not in line:
            out.append(f"the collapsed line does not name the rule: {line!r}")
    if SEEP_FULL in blk:
        out.append("the staged finding was quoted in full a second time")
    mw.deleteLater()
    return out


def check_every_finding_is_quoted_once():
    """A finding may collapse only once it has been QUOTED — never merely once it
    has been counted.

    The quote cap is what makes this a live risk: a model carrying more findings
    than one block quotes names the rest on the overflow line, and if that counted
    as telling the model, the next block would file them under "reported in full
    earlier" — a claim about text the model never received. So the cap rotates:
    what was named and not quoted is quoted in a later block, and until then the
    line naming it carries the command that prints it.
    """
    from studio.ai import assistant as A
    from xslope.preflight import preflight
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)

    blocks = [_block(_run(asst, SNIPPET_MANY_DEAD_CIRCLES))]
    for gamma in (126.0, 127.0, 128.0, 129.0):
        blocks.append(_block(_run(asst, f"slope_data['materials'][0]['gamma'] = "
                                        f"{gamma}\nprint('nudged')")))
    sd = mw.doc.slope_data
    standing = [f for f in preflight(sd, "lem", A._checks_selection(sd)).findings
                if f.severity in ("error", "warning")]
    if len(standing) <= A.MAX_CHECK_FINDINGS:
        out.append(f"{len(standing)} findings do not overflow the cap of "
                   f"{A.MAX_CHECK_FINDINGS} — this leg proves nothing")

    quoted = set()
    for n, blk in enumerate(blocks, 1):
        blk = blk or ""
        for line in _collapsed_lines(blk):
            claimed = int(line.split()[0])
            if claimed > len(quoted):
                out.append(f"block {n}: {claimed} finding(s) filed as reported in "
                           f"full earlier, but only {len(quoted)} have been quoted")
        if "not quoted here" in blk and "print(preflight(" not in blk:
            out.append(f"block {n}: named findings it did not quote without the "
                       f"command that reads them")
        quoted |= {f.message for f in standing if f.message in blk}
    never = [f.rule_id for f in standing if f.message not in quoted]
    if never:
        out.append(f"{len(never)} finding(s) were never quoted in full across "
                   f"{len(blocks)} blocks: {never}")
    mw.deleteLater()
    return out


def check_two_faults_on_one_row():
    """Two different faults about one row are two findings.

    vp047 states one pullout length per end, and both exceed the line: the same
    rule fires twice on one reinforcement line, once for Lp1 and once for Lp2. The
    messages differ only in a digit that belongs to a FIELD NAME, so blanking it
    with the values made the pair one finding — the second was never reported, and
    repairing the first left the model reading the survivor as unchanged text it
    had already been given. Digit-suffixed names run through the whole schema
    (x1/y1, x2/y2, lp1/lp2, k1/k2), so this is a class, not a case.
    """
    from studio.ai import assistant as A
    from xslope.preflight import preflight
    out = []
    deck = os.path.join(_REPO, "docs", "verification", "files", "rocscience",
                        "vp047.xlsx")
    if not os.path.exists(deck):
        return [f"the regression deck is missing: {deck}"]
    mw, asst = _session()
    mw.doc.load(deck)
    sd = mw.doc.slope_data
    pair = [f for f in preflight(sd, "lem", A._checks_selection(sd)).findings
            if f.rule_id == "reinforce.envelope_inconsistent"
            and f.message.startswith("Reinforcement line 1 ")]
    if len(pair) != 2 or not ("Lp1 =" in pair[0].message
                              and "Lp2 =" in pair[1].message):
        mw.deleteLater()
        return [f"vp047 no longer carries the Lp1/Lp2 pair on one row: "
                f"{[f.message[:60] for f in pair]}"]
    # Precisely: the two share a _finding_key, because identity is the rule plus
    # the row and both are about row 1. What separates them is the occurrence
    # counter ChecksMemo.keys appends — which makes them two entries rather than
    # one — and the signature, which is what decides whether the survivor of the
    # pair is a finding the model has been given or a new one.
    if A._finding_sig(pair[0]) == A._finding_sig(pair[1]):
        out.append("the Lp1 and Lp2 faults carry the same signature, so the "
                   "survivor of the pair reads as one already reported")
    if len(set(A.ChecksMemo.keys(pair))) != 2:
        out.append("the two faults on one row are tracked as one entry")

    # ... and through the real path: both quoted, and repairing Lp1 leaves the
    # Lp2 fault quoted in full rather than filed under the text of Lp1.
    blk = _block(_run(asst, "slope_data['materials'][0]['gamma'] += 1.0\n"
                            "print('unit weight nudged')"))
    for f in pair:
        if f.message not in (blk or ""):
            out.append(f"a fault on the pair was not quoted: {f.message[:80]!r}")

    blk = _block(_run(asst, """
for line in slope_data['reinforcement_lines']:
    line['lp1'] = 1.0
print('first pullout length repaired on every line')
"""))
    left = [f for f in preflight(mw.doc.slope_data, "lem",
                                 A._checks_selection(mw.doc.slope_data)).findings
            if f.rule_id == "reinforce.envelope_inconsistent"]
    if not left or any("Lp1 =" in f.message for f in left):
        out.append(f"the repair did not leave the Lp2 faults standing alone: "
                   f"{[f.message[:40] for f in left]}")
    for f in left:
        if f.message not in (blk or ""):
            out.append(f"the surviving fault was written off as already reported: "
                       f"{f.message[:80]!r}")
    mw.deleteLater()
    return out


def check_nothing_is_dropped_unnamed():
    """The cap on quoted findings is a cap on paragraphs, not on accountability:
    what it leaves out is still named by rule key, so iron rule 5 stays literally
    satisfiable — nothing reaches the model as silence."""
    from studio.ai import assistant as A
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    _run(asst, SNIPPET_PIEZO_REVERSED)
    sd = mw.doc.slope_data
    real_max = A.MAX_CHECK_FINDINGS
    try:
        A.MAX_CHECK_FINDINGS = 1
        text = A.model_checks_text(sd)          # no memo: a first report
    finally:
        A.MAX_CHECK_FINDINGS = real_max
    from xslope.preflight import preflight
    # Errors and warnings only: INFO rows are not what the block reports, so
    # requiring them here would fail on a model that trips one.
    ids = {f.rule_id for f in preflight(sd, "lem", A._checks_selection(sd)).findings
           if f.severity in ("error", "warning")}
    if len(ids) < 2:
        out.append("this model carries too few findings to overflow the cap")
    for rid in ids:
        if rid not in text:
            out.append(f"{rid} was dropped from the block without being named")
    mw.deleteLater()
    return out


def check_delta_mutations():
    """Seven mutations, one leg each. Every one of them is a plausible reading of
    'collapse repeats', and each breaks something the block has to keep."""
    from studio.ai import assistant as A
    out = []

    def _all_fresh(self, findings):
        keys = self.keys(findings)
        return keys, set(keys)

    def _mark_all_seen(self, keys, findings, quoted):
        self._seen = {k: A._finding_sig(f) for k, f in zip(keys, findings)}

    mutations = (
        ("nothing is exempt: errors collapse too", A, "_never_collapses",
         lambda f: False, check_error_never_collapses,
         "an ERROR is quoted every time"),
        ("nothing ever collapses", A.ChecksMemo, "delta", _all_fresh,
         check_unchanged_finding_collapses, "an unchanged finding collapses"),
        ("tracking never resets", A.ChecksMemo, "reset", lambda self: None,
         check_fresh_session_reports_in_full, "a fresh session reports in full"),
        ("identity is the rule alone", A, "_finding_key", lambda f: (f.rule_id,),
         check_changed_finding_uncollapses, "a changed finding un-collapses"),
        # Counted is not quoted: mark everything the block SAW as told and the
        # findings the quote cap left out are collapsed having never been read.
        ("counted is treated as quoted", A.ChecksMemo, "record", _mark_all_seen,
         check_every_finding_is_quoted_once,
         "every finding is quoted in full at least once"),
        # The digits inside a field name are part of its name: blank them and the
        # Lp1 and Lp2 faults on one row become one finding.
        ("field-name digits are values", A, "_NUMBER_RE",
         __import__("re").compile(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?"),
         check_two_faults_on_one_row, "two faults on one row stay two"),
        # Rotation without the error tier: an error already quoted sorts behind
        # every new warning and is pushed off the quoted list by a batch of them.
        ("the rotation outranks the error exemption", A, "_quote_order",
         lambda rows, keys, memo: (lambda i: memo.quoted_before(keys[i])),
         check_error_survives_the_quote_cap,
         "an ERROR is quoted even when the cap is full"),
    )
    for name, target, attr, mutant, leg, legname in mutations:
        real = getattr(target, attr)
        setattr(target, attr, mutant)
        try:
            noticed = bool(leg())
        except Exception:
            noticed = True          # a leg that blows up has noticed
        finally:
            setattr(target, attr, real)
        if not noticed:
            out.append(f"mutation ({name}): '{legname}' still passed")
        if leg():
            out.append(f"mutation ({name}): '{legname}' did not come back green")
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
    A.model_checks_text = lambda sd, memo=None, extra_errors=(): ""
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


# --- I. the helpers, the providers, the meter ------------------------------
# Three things the harness owns whatever the model does: which functions the
# snippet namespace can actually call, which providers the dialog can offer, and
# what the turn cost.

#: Every helper the system prompt tells the model to reach for, with the
#: arguments the prompt names. A helper described in the prompt and absent from
#: the namespace is the worst shape of this failure: the model writes the call
#: the prompt taught it, gets a NameError, and reconstructs the pipeline by hand
#: — which is the failure the helpers exist to prevent.
HELPERS = {
    "run_lem": ("method", "num_slices", "plot", "slope_data", "search",
                "seep_time"),
    "corpus_index": ("query",),
    "run_seep": ("bc", "tol", "max_iter", "plot", "slope_data", "transient"),
    "run_tseep": ("times", "plot", "slope_data", "rerun", "save"),
    "fs_vs_time": ("times", "method", "mode", "search", "rapid", "plot"),
    "transient_solution": ("slope_data",),
    "run_fem": ("analysis", "F", "F_min", "F_max", "tolerance", "plot",
                "slope_data"),
    "suggest_elastic": ("material_or_soil_type", "unit_system", "slope_data"),
    "generate_report": ("path", "finalize"),
    "resync_geometry": ("slope_data",),
    "sensitivity": ("values", "apply", "param"),
    "list_params": ("slope_data", "mode"),
    "design_sweep": ("param", "low", "high", "target_fs", "mode"),
    "parametric_sweep": ("params", "plot"),
    "parametric_design": ("param", "low", "high"),
    "parametric_back_analysis": ("param", "low", "high"),
    "reliability": ("method", "engine"),
    "reliability_taylor": ("method",),
    "reliability_mc": ("method", "n_samples"),
    "reliability_rs": ("method",),
}


def check_helpers_are_callable():
    """Every helper the prompt names is in the namespace a snippet runs in, takes
    the arguments the prompt gives it, and is described where the model reads."""
    import inspect
    from studio.ai.assistant import SCHEMA_BRIEF, STUDIO_SYSTEM, _load_brief_text
    out = []
    mw, asst = _session()
    kernel = asst._kernel
    _quiet(kernel.run, "pass")                     # force the seed
    ns = kernel._ns
    for name, args in HELPERS.items():
        fn = ns.get(name)
        if not callable(fn):
            out.append(f"{name}() is not in the kernel namespace")
            continue
        params = inspect.signature(fn).parameters
        for arg in args:
            if arg not in params:
                out.append(f"{name}() takes no {arg!r} argument")
    prompt = STUDIO_SYSTEM + SCHEMA_BRIEF + (_load_brief_text() or "")
    for name in HELPERS:
        if name not in prompt:
            out.append(f"{name}() is never mentioned to the model")
    mw.deleteLater()
    return out


def check_every_offered_model_can_see():
    """The provider list is exactly the providers that can be handed a picture.

    The assistant's headline request is a cross section given as a sketch or a
    photograph. A provider whose models cannot accept one turns that request into
    a conversation about what the picture shows, so it is not offered at a lower
    tier — it is not offered, and neither is a text-only model inside a mixed
    catalogue.
    """
    from studio.ai.config import PROVIDERS, model_is_vision
    from studio.ai.models import LISTING, is_offerable
    out = []
    if "kimi" not in PROVIDERS:
        out.append("Kimi (Moonshot) is not offered")
    if "deepseek" in PROVIDERS:
        out.append("DeepSeek is still offered — its API takes no images")
    if "deepseek" in LISTING:
        out.append("DeepSeek still has a list-models endpoint")
    for provider, spec in PROVIDERS.items():
        for model in spec.get("models") or []:
            if model_is_vision(provider, model) is not True:
                out.append(f"{provider} offers {model!r}, which cannot read an image")
            if not is_offerable(provider, model):
                out.append(f"{provider}'s own list drops {model!r}")
        if spec.get("vision") is False:
            out.append(f"{provider} is declared text-only and still listed")
    # The mixed catalogues are filtered, not merely labelled.
    for provider, text_only in (("kimi", "moonshot-v1-8k"), ("zai", "glm-5.2"),
                                ("ollama", "llama3.1:8b")):
        if is_offerable(provider, text_only):
            out.append(f"{provider} would still list the text-only {text_only!r}")
    # Kimi reaches litellm as its own route, with its key under its own name.
    kimi = PROVIDERS.get("kimi", {})
    if kimi.get("prefix") != "moonshot/":
        out.append(f"Kimi's litellm prefix is {kimi.get('prefix')!r}")
    if not str(kimi.get("base", "")).startswith("https://api.moonshot.ai"):
        out.append(f"Kimi's base URL is {kimi.get('base')!r}")
    return out


def check_the_keychain_name():
    """Kimi's key is stored under ``kimi_api_key`` — the name the ruling fixed,
    and the one an existing install would already carry."""
    out = []
    stored = {}

    class _Keyring:
        @staticmethod
        def set_password(service, user, key):
            stored[(service, user)] = key

        @staticmethod
        def get_password(service, user):
            return stored.get((service, user))

        @staticmethod
        def delete_password(service, user):
            stored.pop((service, user), None)

    from studio.ai import config as C
    original = C._keyring
    C._keyring = lambda: _Keyring
    try:
        mw, asst = _session()
        asst.config.set_api_key("kimi", "sk-kimi-test")
        if (C.KEYRING_SERVICE, "kimi_api_key") not in stored:
            out.append(f"the key went to {sorted(stored)} instead of kimi_api_key")
        if asst.config.api_key("kimi") != "sk-kimi-test":
            out.append("the stored key did not read back")
        mw.deleteLater()
    finally:
        C._keyring = original
    return out


class _FakeUsage:
    """One provider's usage object: attributes, and a nested details object for
    the cached-prompt count."""

    class _Details:
        cached_tokens = 900

    def __init__(self, prompt=1000, completion=250, cached=True):
        self.prompt_tokens = prompt
        self.completion_tokens = completion
        self.prompt_tokens_details = self._Details() if cached else None


class _FakeResponse:
    def __init__(self, usage):
        self.usage = usage


def check_usage_accumulates():
    """Tokens are read off each completion, accumulate per turn and per session,
    and survive a response that carries no usage at all."""
    from studio.ai.assistant import format_usage, usage_from_response
    out = []
    got = usage_from_response(_FakeResponse(_FakeUsage()))
    if got != {"input": 1000, "cached_input": 900, "output": 250}:
        out.append(f"a completion's usage read as {got}")
    # A dict-shaped usage, and Anthropic's own cache-read field name.
    got = usage_from_response(_FakeResponse(
        {"prompt_tokens": 10, "completion_tokens": 2,
         "cache_read_input_tokens": 8}))
    if got != {"input": 10, "cached_input": 8, "output": 2}:
        out.append(f"a dict usage read as {got}")
    for empty in (_FakeResponse(None), object(), None):
        if usage_from_response(empty) != {"input": 0, "cached_input": 0, "output": 0}:
            out.append(f"a response with no usage did not read as zero: {empty!r}")

    mw, asst = _session()
    seen = []
    asst.usage_changed.connect(seen.append)
    if asst.usage["session"] != {"input": 0, "cached_input": 0, "output": 0}:
        out.append("a fresh session did not start at zero")
    for _ in range(3):
        asst.add_usage(usage_from_response(_FakeResponse(_FakeUsage())))
    if asst.usage["session"] != {"input": 3000, "cached_input": 2700, "output": 750}:
        out.append(f"three completions accumulated to {asst.usage['session']}")
    if asst.usage["turn"] != asst.usage["session"]:
        out.append("the turn and the session diverged inside one turn")
    if len(seen) != 3:
        out.append(f"the dock was told {len(seen)} time(s), expected 3")
    # A new turn counts from zero while the session keeps running.
    asst._usage["turn"] = {"input": 0, "cached_input": 0, "output": 0}
    asst.add_usage(usage_from_response(_FakeResponse(_FakeUsage(prompt=5, completion=1,
                                                               cached=False))))
    if asst.usage["turn"] != {"input": 5, "cached_input": 0, "output": 1}:
        out.append(f"the new turn read as {asst.usage['turn']}")
    if asst.usage["session"]["input"] != 3005:
        out.append("the session was reset by a new turn")
    line = format_usage(asst.usage["turn"], asst.usage["session"])
    for fragment in ("this turn:", "session:", " in", " out", "3,005"):
        if fragment not in line:
            out.append(f"the readout {line!r} is missing {fragment!r}")
    # A new conversation starts the meter over.
    _quiet(asst.reset)
    if asst.usage["session"] != {"input": 0, "cached_input": 0, "output": 0}:
        out.append("a new chat kept the old session's tokens")
    mw.deleteLater()
    return out


# --- J. the brief is what ships, and the skill body is not ------------------
# The in-app assistant used to be handed the whole Claude Code skill — ~34k tokens
# of file-first workflow, corpus tables and template layout, re-read on every
# completion of every turn, six of them in a measured single-question turn. It is
# replaced by a Studio brief that ships as package data. What is asserted here is
# the shipping and the size, plus the negative that matters: no part of the skill
# body is in the prompt any more, whatever the skill says next month.

#: Ceiling on the brief, in the offline counter :func:`count_tokens` uses (a check
#: must not make a network call to know whether it passes). That counter reads
#: about 1.48x low against Anthropic's own count_tokens endpoint — the brief
#: measures 4,347 here and 6,432 there — so this bound is the ~9k of billed tokens
#: the brief is allowed, expressed in the units the check can measure.
BRIEF_TOKEN_BOUND = 6500  # offline counter; ~1.48x low against the billed count, so ~9.6k billed


def check_the_brief_ships():
    """The brief is package data, is what the loader returns, and fits the budget."""
    from importlib import resources
    from studio.ai.assistant import BRIEF_RESOURCE, _load_brief_text, count_tokens
    out = []
    brief = _load_brief_text()
    if not brief.strip():
        return ["the Studio brief could not be loaded"]
    try:
        shipped = (resources.files("xslope") / "resources"
                   / BRIEF_RESOURCE).read_text(encoding="utf-8")
    except Exception as exc:
        return [f"the brief is not package data under xslope/resources: {exc!r}"]
    if shipped != brief:
        out.append("the loader does not return the shipped brief")
    n = count_tokens(brief)
    if n is None:
        print("      (litellm unavailable — brief token bound not measured)")
    elif n > BRIEF_TOKEN_BOUND:
        out.append(f"the brief measures {n:,} tokens, over the "
                   f"{BRIEF_TOKEN_BOUND:,} bound")
    for anchor in ("run_lem(", "corpus_index(", "gamma_water", "MODEL CHECKS"):
        if anchor not in brief:
            out.append(f"the brief never mentions {anchor!r}")
    return out


def check_skill_body_is_not_in_the_prompt():
    """No part of the Claude Code skill reaches the Studio system prompt.

    Sampled rather than spot-checked on a phrase: a named sentence stops proving
    anything the moment the skill is edited, whereas a spread of slices of whatever
    the skill says today cannot go stale.
    """
    from importlib import resources
    out = []
    try:
        skill = (resources.files("xslope") / "resources"
                 / "xslope_skill.md").read_text(encoding="utf-8")
    except Exception:
        return ["the shipped skill could not be read, so its absence is unproven"]
    slices = [skill[i:i + 120]
              for i in range(0, max(len(skill) - 120, 1), max(len(skill) // 40, 1))]
    for tier, prompt in _prompts().items():
        hits = [sl for sl in slices if sl and sl in prompt]
        if hits:
            out.append(f"{tier}: the skill body is still in the prompt "
                       f"({len(hits)} of {len(slices)} sampled slices matched)")
    return out


# --- K. the model summary ---------------------------------------------------
# The other half of the cost: the model opened a turn by discovering what it was
# looking at — dump slope_data, dump the materials, help(run_lem) — three or four
# completions before the one that did the work. The summary is given instead, once,
# and again whenever an edit has made it wrong.

def check_model_summary_is_injected():
    """First turn carries it, the next does not, and an edit brings it back."""
    from studio.ai.assistant import MODEL_SUMMARY_END, MODEL_SUMMARY_OPEN
    out = []
    mw, asst = _session()

    def carries(text):
        return MODEL_SUMMARY_OPEN in text and MODEL_SUMMARY_END in text

    first = asst._user_content("what is this model?")
    if not carries(first):
        out.append("the first turn carries no model summary")
    if "what is this model?" not in first:
        out.append("the user's own words were lost from the first turn")
    if carries(asst._user_content("and now?")):
        out.append("an unchanged model is described twice")

    _run(asst, "slope_data['max_depth'] = -12.0")
    if not carries(asst._user_content("after the edit")):
        out.append("an edit does not refresh the model summary")
    if carries(asst._user_content("and after that")):
        out.append("the refreshed summary repeats on the following turn")

    # An image turn is the multimodal content shape; the summary must ride in its
    # text part rather than being dropped on the floor.
    _run(asst, "slope_data['max_depth'] = -14.0")
    content = asst._user_content("look", images=["data:image/png;base64,AAAA"])
    if not isinstance(content, list):
        out.append("an image turn no longer produces multimodal content")
    elif not carries(content[0].get("text", "")):
        out.append("an image turn drops the model summary")
    mw.deleteLater()
    return out


def check_model_summary_says_what_the_model_is():
    """It names the things a first snippet used to be spent discovering."""
    from studio.ai.assistant import model_summary_text
    sd = {"unit_system": "imperial", "gamma_water": 62.4, "max_depth": -10.0,
          "circular": True, "template_version": 24,
          "materials": [{"name": "shell", "gamma": 130.0, "option": "mc",
                         "c": 300.0, "phi": 37.0, "u": "none"},
                        {"name": "clay", "gamma": 120.0, "option": "mc",
                         "c": 0.0, "phi": 0.0, "u": "piezo"}],
          "profile_lines": [{"mat_id": 0, "coords": [(0, 0), (30, 24)]}],
          "circles": [{"Xo": 0.0, "Yo": 40.0, "Depth": 0.0, "R": 40.0}],
          "non_circ": [], "piezo_line": [(0, 5), (100, 5)],
          "dloads": [[{"X": 30, "Y": 24, "Normal": 240}]],
          "reinforcement_lines": [{}] * 6, "pile_lines": [], "line_loads": []}
    text = model_summary_text(sd, {"lem_solution": object()},
                              name="thing.xlsx")
    out = []
    for want in ("thing.xlsx", "imperial", "max_depth=-10", "shell", "clay",
                 "phi=0", "gamma_water=62.4", "piezo line", "1 circle(s)",
                 "6 reinforcement line(s)", "an LEM solution"):
        if want not in text:
            out.append(f"the model summary never says {want!r}")
    if "help()" not in text:
        out.append("the model summary does not tell the model the helpers are "
                   "preloaded")
    # A model with nothing open must still produce a block, not an exception.
    if "No project is open" not in model_summary_text(None):
        out.append("an empty session produces no summary")
    return out


def check_corpus_index_returns_rows():
    """The corpus pointer resolves: a topic query comes back with real pages."""
    out = []
    mw, asst = _session()
    kernel = asst._kernel
    _quiet(kernel.run, "pass")
    corpus_index = kernel._ns.get("corpus_index")
    if not callable(corpus_index):
        mw.deleteLater()
        return ["corpus_index() is not in the kernel namespace"]
    topics = _quiet(corpus_index)
    if not topics:
        out.append("corpus_index() lists no topics")
    for query in ("rapid drawdown", "piles", "reinforcement"):
        rows = _quiet(corpus_index, query)
        if not rows:
            out.append(f"corpus_index({query!r}) returned no rows")
            continue
        for row in rows:
            if not str(row.get("url", "")).startswith("https://"):
                out.append(f"corpus_index({query!r}) returned a row with no URL")
                break
            if not row.get("title"):
                out.append(f"corpus_index({query!r}) returned an untitled row")
                break
    if _quiet(corpus_index, "zzzz no such topic zzzz"):
        out.append("corpus_index() invents rows for a query that matches nothing")
    mw.deleteLater()
    return out


# --- L. a run is the session's run -----------------------------------------
# The eighth recorded session asked for a report on a run it had just made, and
# got a three-page document with no results section in it -- because `run_lem`
# computed an answer and told nobody. What follows is that whole chain, asserted:
# the bundle lands where the report reads it, the report then documents it, the
# method is the model's own rather than a hardcoded one, and the answer says
# which surface it is the answer for.

def _solved_session():
    """A live session on the standard section, with one solvable circle and one
    single-surface LEM run already made. Returns ``(window, assistant, result)``
    where ``result`` is the tool output of the run."""
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    _run(asst, SNIPPET_ONE_CIRCLE)
    out = _run(asst, "res = run_lem(search=False, plot=False)\n"
                     "print('CIRCLE', {k: res.get(k) for k in "
                     "('Xo', 'Yo', 'R', 'Depth', 'x_entry', 'x_exit')})")
    return mw, asst, out


def _docx_text(path):
    """The visible text of a .docx, tags stripped."""
    import re
    import zipfile
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8", "replace")
    return re.sub(r"<[^>]+>", "", xml)


def check_run_lem_is_the_sessions_run():
    """`run_lem` stores its bundle where the report reads it, and the report
    built afterwards carries the run."""
    from studio.main_window import MainWindow
    out = []
    mw, asst, _ = _solved_session()
    bundle = mw.doc.results.get("lem_solution")
    if not isinstance(bundle, dict):
        mw.deleteLater()
        return ["run_lem() left no bundle on doc.results['lem_solution']"]
    for key in ("slice_df", "failure_surface", "results", "search", "method",
                "options"):
        if key not in bundle:
            out.append(f"the stored LEM bundle has no {key!r}")
    solutions = MainWindow.report_solutions(mw)
    if not solutions.get("lem"):
        out.append("MainWindow.report_solutions() finds no LEM run to document")
    elif solutions["lem"][0].get("method") != bundle.get("method"):
        out.append("the report would label the run with a method it was not run "
                   f"under ({solutions['lem'][0].get('method')!r} vs "
                   f"{bundle.get('method')!r})")
    # "Lay the face back and rerun" is ONE snippet, and the stale-result sweep
    # the edit triggers runs after it — so the run the user was just shown has to
    # survive a sweep aimed at the solution that predates the edit.
    _run(asst, "slope_data['profile_lines'][0]['coords'] = "
               "[(0.0, 5.0), (45.0, 20.0), (90.0, 20.0)]\n"
               "print(run_lem(search=False, plot=False)['FS'])")
    if not mw.doc.results.get("lem_solution"):
        out.append("a snippet that edited the model and then ran it kept no run: "
                   "the stale-result sweep took the run made after the edit")
    # And the document that comes out of it: the results section, and the
    # traceability stamp that names the file the model was opened from.
    bundle = mw.doc.results.get("lem_solution") or bundle
    fs = (bundle.get("results") or {}).get("FS")
    tmpdir = tempfile.mkdtemp(prefix="xslope_assistant_report_")
    model_path = os.path.join(tmpdir, "small_slope.xlsx")
    with open(model_path, "wb") as fh:
        fh.write(b"not a workbook -- only the traceability digest reads it")
    mw.doc.path = model_path
    report_path = os.path.join(tmpdir, "report.docx")
    result = _run(asst, f"print('REPORT', generate_report(path={report_path!r}, "
                        f"finalize=False))")
    if not os.path.exists(report_path):
        out.append(f"generate_report() wrote no document: {result[-400:]}")
    else:
        text = _docx_text(report_path)
        if "Limit Equilibrium" not in text:
            out.append("the report has no limit equilibrium section — the run it "
                       "was asked to document is missing from it")
        if fs is not None and f"{fs:.3f}" not in text:
            out.append(f"the report never states the factor of safety ({fs:.3f})")
        if "not saved to a file" in text:
            out.append("the traceability stamp says the model is not saved to a "
                       "file, for a project opened from one")
        if os.path.basename(model_path) not in text:
            out.append("the traceability stamp does not name the input file")
        if "SHA-256" not in text:
            out.append("the traceability stamp carries no SHA-256 of the input")
    mw.deleteLater()
    return out


def check_run_lem_runs_the_models_method():
    """A run made with no method named runs the method the MODEL declares.

    Two recorded sessions asked to "rerun" a Spencer model and got Bishop, from a
    hardcoded default. The model's declared method is what the Run LEM dialog
    opens on, so it is what an unqualified run means.
    """
    from studio.ai.kernel import _declared_lem_method
    out = []
    cases = [({}, "spencer"), ({"lem_method": "janbu"}, "janbu"),
             ({"lem_method": "all"}, "spencer"), ({"lem_method": None}, "spencer"),
             ({"lem_method": "bishop"}, "bishop")]
    for sd, want in cases:
        got = _declared_lem_method(sd)
        if got != want:
            out.append(f"a model declaring {sd.get('lem_method')!r} defaults to "
                       f"{got!r}, expected {want!r}")
    if _declared_lem_method({"lem_method": "janbu"}, "oms") != "oms":
        out.append("a named method does not win over the model's declaration")
    # And live: the declared method is the one that actually runs.
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    _run(asst, SNIPPET_ONE_CIRCLE)
    _run(asst, "slope_data['lem_method'] = 'janbu'\nprint('method declared')")
    _run(asst, "run_lem(search=False, plot=False)")
    ran = (mw.doc.results.get("lem_solution") or {}).get("method")
    if ran != "janbu":
        out.append(f"a model declaring janbu was run as {ran!r}")
    mw.deleteLater()
    return out


def check_run_lem_returns_the_surface():
    """The result dict says which surface it is the answer for.

    A session crashed with ``KeyError: 'Xo'`` reading the critical circle off a
    run's result, which carried no geometry at all.
    """
    out = []
    mw, asst, result = _solved_session()
    for key in ("Xo", "Yo", "R", "Depth", "x_entry", "x_exit"):
        if f"'{key}': None" in result or f"'{key}'" not in result:
            out.append(f"the result dict carries no usable {key!r}")
    # The circle it reports is the circle the model defines.
    circle = mw.doc.slope_data["circles"][0]
    if f"'Xo': {float(circle['Xo'])}" not in result:
        out.append("the result's Xo is not the solved circle's center")
    mw.deleteLater()
    return out


def check_polygon_edit_on_a_profile_model_is_named():
    """A geometry edit made on the source the resync OVERWRITES is reported.

    On a profile-line model the polygons are rebuilt from the profile lines after
    every snippet, so a face laid back by rewriting `polygons` is reverted — and
    the input checks call the result clean, because after the rebuild the geometry
    is valid. The warning is the only thing that says the edit did not take.
    """
    from studio.ai.kernel import POLYGON_EDIT_WARNING
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    result = _run(asst, SNIPPET_POLYGON_EDIT_ON_PROFILE_MODEL)
    if POLYGON_EDIT_WARNING not in result:
        out.append("a polygon edit on a profile-line model is not reported: "
                   f"{result[:400]!r}")
    if "profile_lines" not in result:
        out.append("the warning does not name the source to edit instead")
    # The premise of the warning: the edit really was discarded.
    face = list(mw.doc.slope_data["profile_lines"][0]["coords"])
    if face != [tuple(c) for c in GROUND]:
        out.append("the profile lines moved, so the polygon edit was not the "
                   "silent revert this warning is about")
    mw.deleteLater()

    # The same edit on a model whose polygons ARE the source is just an edit.
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_POLYGON_NATIVE_BUILD)
    result = _run(asst, SNIPPET_POLYGON_NATIVE_EDIT)
    if "WARNING:" in result:
        out.append("editing polygons on a polygon-native model is warned about: "
                   f"{result[:400]!r}")
    mw.deleteLater()
    return out


def check_discarded_polygon_edit_cannot_read_as_applied():
    """A discarded polygon edit is stated BEFORE the read-back that contradicts it,
    and the checks block calls it an error.

    The snippet's prints run before the resync, so a session that edits `polygons`
    and prints them back sees its own edit in its own output: two benchmark
    sessions did that, read the warning underneath, and still delivered the model
    as repaired. So the discard leads the tool result — ahead of every printed
    value it invalidates — and rides in the MODEL CHECKS block as an ERROR, which
    is the block iron rule 5 forbids reporting a model ready over. The input
    checks cannot produce it themselves: after the rebuild the model is valid, and
    merely not the one the snippet wrote.
    """
    from studio.ai.kernel import (POLYGON_EDIT_DISCARDED,
                                  POLYGON_EDIT_DISCARDED_FINDING)
    rule = POLYGON_EDIT_DISCARDED_FINDING[0]
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    result = _run(asst, SNIPPET_POLYGON_EDIT_READ_BACK)

    if not result.startswith(POLYGON_EDIT_DISCARDED):
        out.append("the tool result does not OPEN on the discard: "
                   f"{result[:300]!r}")
    if "READ BACK" in result and result.index("READ BACK") < len(
            POLYGON_EDIT_DISCARDED):
        out.append("the snippet's read-back is printed before the discard is said")
    for phrase in ("discarded", "profile_lines", "mat_id"):
        if phrase not in POLYGON_EDIT_DISCARDED.lower():
            out.append(f"the discard line does not say {phrase!r}")
    block = _block(result)
    if block is None:
        out.append("a discarded polygon edit produced no MODEL CHECKS block")
    elif f"ERROR [{rule}]" not in block:
        out.append("the MODEL CHECKS block carries no error for the discard: "
                   f"{block[:400]!r}")
    if result.count(POLYGON_EDIT_DISCARDED_FINDING[1]) > 1:
        out.append("the discard finding is quoted more than once")
    # The premise: the model really is unchanged by that snippet.
    zones = mw.doc.slope_data.get("polygons") or []
    if any(z.get("mat_id") == 1 for z in zones):
        out.append("the polygon edit survived, so this is not the discard case")
    mw.deleteLater()

    # The control: the same model, edited on the source it IS built from. A
    # profile_lines edit is an ordinary edit — no discard line, no error.
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)
    _run(asst, SNIPPET_GEOMETRY)
    result = _run(asst, SNIPPET_RAISE_CREST)
    if POLYGON_EDIT_DISCARDED in result:
        out.append("a profile_lines edit was reported as a discarded polygon edit")
    if rule in result:
        out.append("a profile_lines edit raised the discard error: "
                   f"{result[:300]!r}")
    mw.deleteLater()
    return out


def check_chat_renders_markdown():
    """The transcript renders the assistant's markdown.

    A reply is markdown — headings, tables, fenced code — and the dock used to
    escape it and break the lines, so a table of results reached the user as rows
    of pipes and a heading as a line of hashes.
    """
    from studio.chat_dock import ChatDock
    out = []
    mw, asst = _session()
    dock = ChatDock(asst)
    reply = ("Here is what the sweep found.\n\n"
             "## Results\n\n"
             "| Tmax | FS |\n| ---: | ---: |\n| 500 | 1.425 |\n| 3000 | 1.681 |\n\n"
             "Run it yourself with `run_lem(search=True)`:\n\n"
             "```python\nres = run_lem(search=True)\n```\n\n"
             "- the plateau is the circle migrating\n"
             "- FS < 1.5 below 800 lb/ft\n")
    dock._on_assistant_text(reply)
    text = dock.transcript.toPlainText()
    if "|" in text:
        out.append("the table reached the transcript as raw pipes")
    if "##" in text:
        out.append("the heading reached the transcript as raw hashes")
    if "```" in text:
        out.append("the code fence reached the transcript as raw backticks")
    for wanted in ("Results", "1.425", "res = run_lem(search=True)",
                   "the plateau is the circle migrating", "FS < 1.5"):
        if wanted not in text:
            out.append(f"the rendered reply lost {wanted!r}")
    html = dock.transcript.toHtml()
    if "<table" not in html:
        out.append("the table is not rendered as a table")
    if "monospace" not in html:
        out.append("the code block is not monospaced")
    # The user's own text is still verbatim, markdown or not.
    dock._add_block("You", "use | pipes | and ## hashes", "#1a5fb4")
    if "| pipes |" not in dock.transcript.toPlainText():
        out.append("the user's own message is no longer shown as typed")
    dock.deleteLater()
    mw.deleteLater()
    return out


def check_chat_degrades_latex():
    """Math reaches the reader as text, whichever way the model writes it.

    The dialect Qt renders has no math, so the moment the transcript began
    rendering markdown the assistant reached for LaTeX and a display equation
    arrived as its own source — literal ``$$\\frac{\\tan\\varphi}{...}$$`` in the
    reply. Both ends are checked: the brief tells the model to write math as
    plain text, and the dock converts LaTeX that arrives anyway.
    """
    from studio.ai.assistant import _load_brief_text
    from studio.chat_dock import ChatDock, strip_latex
    out = []

    brief = _load_brief_text()
    if "LaTeX" not in brief:
        out.append("the brief never tells the model to avoid LaTeX")
    if "tan φ / tan β" not in brief:
        out.append("the brief gives no plain-text example of an equation")

    mw, asst = _session()
    dock = ChatDock(asst)
    reply = ("The reinforcement is what carries it.\n\n"
             "$$\\frac{\\tan\\varphi}{\\tan\\beta} = "
             "\\frac{\\tan 3°}{\\tan 38.7°} = 0.066$$\n\n"
             "So $\\varphi$ contributes almost nothing at $\\beta$ = 38.7°, and "
             "$T \\le 3000$ lb/ft with $\\sigma_{n}$ in units of $x^{2}$. "
             "The answer is $FS = 1.587$ on that circle.\n")
    dock._on_assistant_text(reply)
    text = dock.transcript.toPlainText()
    if "$" in text:
        out.append(f"a dollar delimiter reached the transcript: {text!r}")
    if "\\" in text:
        out.append(f"a backslash command reached the transcript: {text!r}")
    for wanted in ("tan φ / tan β = tan 3° / tan 38.7° = 0.066",
                   "φ contributes almost nothing at β = 38.7°",
                   "T ≤ 3000 lb/ft", "σₙ", "x²", "FS = 1.587"):
        if wanted not in text:
            out.append(f"the converted math lost {wanted!r}: {text!r}")

    # What must NOT be touched: a price, a shell variable, and any code span —
    # a snippet means its dollar signs and backslashes literally.
    kept = ("At $5.00 / MTok input and $0.50 cached, `echo $HOME` costs "
            "`\\frac{a}{b}` nothing.")
    if strip_latex(kept) != kept:
        out.append(f"prose prices or code spans were rewritten as math: "
                   f"{strip_latex(kept)!r}")
    dock.deleteLater()
    mw.deleteLater()
    return out


# --- M. what gets run, and where a file lands -------------------------------
# Two more live sessions, two harness faults — both of them the harness telling
# the user something that was not so.
#
# One wrote a documentation signature into its answer to show the user the shape
# of a call — `sigma_from_range(hcv, lcv, n=None)`, in a ```python fence. The
# no-tool-call fallback scraped that fence out of the PROSE and ran it; the
# NameError came back, and the assistant explained it to the user as the user's
# own error. A call that compiles is not a call that was meant to run, and no
# amount of reading the code can tell the two apart. So: a provider that HAS tool
# calling never reads code out of prose at all, and the fallback kept for the
# providers without one runs only a fence marked ```python run.
#
# The other generated a report and told the user to open it from the Files button.
# The document was real and the sentence was false: with no path it went beside
# the project, which is neither the folder that button opens nor among the files
# the snippet reports back. So the default lands in OUTPUT_DIR.

#: The message that was executed: a signature written to be READ, which happens
#: to be a syntactically valid call.
PROSE_FENCE_MESSAGE = (
    "The helper's shape is:\n\n"
    "```python\n"
    "sigma_from_range(hcv, lcv, n=None)\n"
    "```\n\n"
    "Pass it the high and low conceivable values and it returns the standard "
    "deviation.")

#: The same message written the one way that ASKS for execution.
MARKED_FENCE_MESSAGE = (
    "Reading the base elevation:\n\n"
    "```python run\n"
    "print(slope_data['max_depth'])\n"
    "```\n")


def _drive_worker(content, text_fallback):
    """Run the real agent loop over one scripted assistant message.

    LiteLLM is replaced for the duration by a stub that hands back ``content``
    and then a plain sentence, so the loop — the actual decision under test —
    runs with no provider and no network. Returns ``(code_that_ran, failures)``.
    """
    import types
    from studio.ai import assistant as A

    ran = []

    class _Worker(A._AgentWorker):
        def _run(self, name, args):     # the GUI-thread hop, short-circuited
            ran.append(args.get("code", ""))
            return "(no output)"

    scripted = [content]

    def completion(**kwargs):
        text = scripted.pop(0) if scripted else "Done."
        msg = types.SimpleNamespace(content=text, tool_calls=None)
        return types.SimpleNamespace(
            choices=[types.SimpleNamespace(message=msg)], usage=None)

    stub = types.ModuleType("litellm")
    stub.drop_params = False
    stub.suppress_debug_info = False
    stub.completion = completion

    saved = sys.modules.get("litellm")
    sys.modules["litellm"] = stub
    try:
        worker = _Worker({}, "system", [], [A.RUN_PYTHON_TOOL],
                         text_fallback=text_fallback)
        failures = []
        worker.failed.connect(failures.append)
        worker.run()
    finally:
        if saved is None:
            sys.modules.pop("litellm", None)
        else:
            sys.modules["litellm"] = saved
    return ran, failures


def check_prose_fence_is_not_code():
    """`_extract_code` reads only a fence marked for execution."""
    from studio.ai import assistant as A
    out = []
    code, pure = A._extract_code(PROSE_FENCE_MESSAGE)
    if code is not None:
        out.append(f"a signature quoted in prose is still read as code: {code!r}")
    code, pure = A._extract_code(MARKED_FENCE_MESSAGE)
    if code != "print(slope_data['max_depth'])":
        out.append(f"a ```python run fence was not recovered: {code!r}")
    if pure:
        out.append("a marked fence suppressed the prose around it")
    # An unmarked fence stays prose whatever is in it — including a snippet that
    # would have been perfectly reasonable to run.
    if A._extract_code("Try this:\n\n```python\nprint(1 + 1)\n```\n")[0] is not None:
        out.append("an unmarked fence is executed when its code looks runnable")
    # The other two recovery shapes are unambiguous and stay: a message that is
    # only Python, and a text-encoded tool call.
    if A._extract_code("print(slope_data['max_depth'])")[0] is None:
        out.append("a message that is only Python is no longer recovered")
    json_call = ('{"name": "run_python", "arguments": {"code": "print(1)"}}')
    if A._extract_code(json_call)[0] != "print(1)":
        out.append("a text-encoded tool call is no longer recovered")
    return out


def check_the_fallback_is_off_where_tools_exist():
    """Only a provider without tool calling reads code out of a message — and it
    is the only one told the convention for saying so."""
    from studio.ai import assistant as A
    from studio.ai.config import PROVIDERS
    out = []
    for tools, want in ((True, False), (False, True), (None, True)):
        got = A.text_fallback_enabled({"tools": tools})
        if got is not want:
            out.append(f"a provider with tools={tools!r} has the text fallback "
                       f"{'on' if got else 'off'}")
    for name, spec in PROVIDERS.items():
        enabled = A.text_fallback_enabled({"tools": spec.get("tools")})
        if spec.get("tools") is True and enabled:
            out.append(f"{name} does tool calling but still scrapes prose")
    if not A.text_fallback_enabled({"tools": PROVIDERS["ollama"].get("tools")}):
        out.append("a local model with no tool calling has no way to run code")

    # And the note that states the convention rides with the fallback, not with
    # every prompt: a model that has tool calling is never told about a channel
    # it does not use.
    mw, asst = _session()
    caps = dict(asst.config.capabilities())
    asst.config.capabilities = lambda: dict(caps, tools=True)
    if A.RUN_FENCE_NOTE in asst._system():
        out.append("a tool-calling provider is told the fence convention")
    asst.config.capabilities = lambda: dict(caps, tools=None)
    system = asst._system()
    if A.RUN_FENCE_NOTE not in system:
        out.append("a model with no tool calling is never told how to run code")
    if "```python run" not in system:
        out.append("the prompt never shows the marker the fallback honours")
    mw.deleteLater()
    return out


def check_the_loop_honours_the_gate():
    """The agent loop itself, driven over a scripted message with no provider."""
    out = []
    ran, failures = _drive_worker(PROSE_FENCE_MESSAGE, text_fallback=False)
    if ran:
        out.append(f"a tool-calling provider executed a fence out of prose: {ran!r}")
    ran, failures = _drive_worker(PROSE_FENCE_MESSAGE, text_fallback=True)
    if ran:
        out.append(f"the fallback executed an unmarked prose fence: {ran!r}")
    ran, failures = _drive_worker(MARKED_FENCE_MESSAGE, text_fallback=True)
    if ran != ["print(slope_data['max_depth'])"]:
        out.append(f"the fallback did not run a marked fence: {ran!r}")
    ran, failures = _drive_worker(MARKED_FENCE_MESSAGE, text_fallback=False)
    if ran:
        out.append(f"a tool-calling provider ran a marked fence: {ran!r}")
    if failures:
        out.append(f"the loop failed: {failures!r}")
    return out


def check_report_lands_in_the_files_folder():
    """A report generated with no path is one the user can actually open.

    It goes to the assistant's output folder — what the dock's Files button opens
    — and the snippet's own result names it, so the sentence the assistant writes
    about opening it is true.
    """
    out = []
    mw, asst, _ = _solved_session()
    outdir = asst.output_dir()
    tmpdir = tempfile.mkdtemp(prefix="xslope_assistant_report_")
    model_path = os.path.join(tmpdir, "small_slope.xlsx")
    with open(model_path, "wb") as fh:
        fh.write(b"not a workbook -- only the traceability digest reads it")
    mw.doc.path = model_path
    beside_the_project = os.path.join(tmpdir, "small_slope_report.docx")
    in_the_folder = os.path.join(outdir, "small_slope_report.docx")
    if os.path.exists(in_the_folder):
        os.remove(in_the_folder)        # a copy from an earlier run proves nothing

    result = _run(asst, "print('REPORT', generate_report(finalize=False))")

    if os.path.exists(beside_the_project):
        out.append("a report asked for with no path was written beside the "
                   "project, outside the folder the Files button opens")
    if not os.path.exists(in_the_folder):
        out.append(f"no report at {in_the_folder}: {result[-400:]}")
    if f"REPORT {in_the_folder}" not in result:
        out.append("generate_report() returned a path other than the one in the "
                   f"output folder: {result[-400:]}")
    if "small_slope_report.docx" not in result.split("[saved", 1)[-1]:
        out.append("the report is not among the files the snippet reports the "
                   "user can open")
    # An explicit path still goes exactly where it was asked to go.
    asked_for = os.path.join(tmpdir, "named.docx")
    _run(asst, f"generate_report(path={asked_for!r}, finalize=False)")
    if not os.path.exists(asked_for):
        out.append("an explicit report path was overridden by the default")
    mw.deleteLater()
    return out


# --- N. transient seepage, and the snippet that never returns ---------------
# Two failures measured in the same corpus sweep. A transient model had no helper
# at all, so the assistant went looking for one — `inspect`, `getsource`, a
# reverse-engineered `build_seep_data()` call missing its first argument — and
# spent the turn's whole completion budget on the search. And a snippet that never
# returned took the GUI thread with it, which is also the thread the per-turn
# timeout runs on, so nothing ended it: one measured run sat there for forty
# minutes. What is asserted: the transient helper exists, is named where the model
# reads, and marches a real sample into frames; and a snippet that overruns the
# limit is stopped, rolled back, and reported to the model as its own result.

#: A real transient model — the earth dam, whose march is a few seconds once the
#: duration is cut to two saved steps. Cutting it keeps the leg about the helper
#: (the pipeline, the frames, where they are stored) rather than about the solver,
#: which has its own verification.
TSEEP_MODEL = os.path.join(_REPO, "docs", "seep", "files",
                           "xslope_earth_dam_tseep.xlsx")


def check_transient_helper_marches():
    """`run_seep(transient=True)` / `run_tseep()` runs a real transient model and
    leaves its frames where Studio leaves them."""
    if not os.path.exists(TSEEP_MODEL):
        return ["the transient sample is missing: %s" % TSEEP_MODEL]
    out = []
    mw, asst = _session()
    _quiet(mw.doc.load, TSEEP_MODEL)
    sd = mw.doc.slope_data
    if not sd.get("tseep"):
        mw.deleteLater()
        return ["the transient sample carries no tseep sheet"]
    # Two saved steps is a march; the file's own 360-day schedule is a minute of
    # solver time this check has no reason to spend.
    sd["tseep"]["duration"] = 20.0
    sd["tseep"]["save_interval"] = 10.0

    result = _run(asst, "b = run_seep(transient=True, plot=False)\n"
                        "print('MODE', b['mode'])\n"
                        "print('FRAMES', [round(f['time'], 3) for f in b['frames']])")
    if "Traceback" in result:
        mw.deleteLater()
        return ["the transient run raised:\n" + result[-600:]]
    if "MODE transient" not in result:
        out.append("the bundle is not a transient bundle: %s" % result[:200])
    if "FRAMES [" not in result or "FRAMES []" in result:
        out.append("the march produced no frames: %s" % result[:200])
    bundle = mw.doc.results.get("transient_seep")
    if not isinstance(bundle, dict):
        out.append("the march left nothing on doc.results['transient_seep'] — "
                   "the Seep · Transient tab and the report both read it there")
    else:
        frames = bundle.get("frames") or []
        if len(frames) < 2:
            out.append("the stored bundle carries %d frame(s)" % len(frames))
        elif any("u" not in fr or "head" not in fr for fr in frames):
            out.append("a stored frame carries no head/u field")
        if bundle.get("transient") is None:
            out.append("the stored bundle carries no march to read an instant from")
    # The march is not paid for twice: a second call hands back what the session
    # already holds (which is also what an opened file's sidecar leaves there).
    again = _run(asst, "print('SAME', run_tseep(plot=False) is "
                       "results['transient_seep'])")
    if "SAME True" not in again:
        out.append("a second run_tseep() re-marched instead of returning the "
                   "solution the session holds: %s" % again[:200])
    # And one instant of it reaches a stability run as pore pressure.
    staged = _run(asst, "info = xslope.seep.apply_transient_stability_frame("
                        "slope_data, transient_solution(), time=10.0)\n"
                        "print('U', slope_data.get('seep_u') is not None)")
    if "U True" not in staged:
        out.append("an instant of the march does not reach slope_data['seep_u']: "
                   "%s" % staged[:300])
    mw.deleteLater()
    return out


def check_transient_helper_is_told_to_the_model():
    """The transient route is in the brief the model reads — the helper, the way
    to pick an instant, and the rule that a refused helper is a call to fix."""
    from studio.ai.assistant import SCHEMA_BRIEF, STUDIO_SYSTEM, _load_brief_text
    out = []
    prompt = STUDIO_SYSTEM + SCHEMA_BRIEF + (_load_brief_text() or "")
    for anchor in ("run_tseep(", "run_seep(transient=True)", "seep_time",
                   "fs_vs_time("):
        if anchor not in prompt:
            out.append(f"the prompt never mentions {anchor!r}")
    brief = _load_brief_text() or ""
    if "tseep" not in brief:
        out.append("the brief never names the tseep sheet, so a transient model "
                   "is not recognizable from it")
    # The rule the reverse-engineering failure needs: read the error, fix the
    # call — do not go reading the package.
    lowered = brief.lower()
    if not any(w in lowered for w in ("getsource", "inspect")):
        out.append("the brief never rules out reading the package when a helper "
                   "refuses, which is what the measured failure did")
    return out


def check_a_runaway_snippet_is_stopped():
    """A snippet that overruns the limit is stopped, and stopped where the model
    can read it: as its own tool result, with the edit rolled back."""
    import time
    from studio.ai import kernel as K
    out = []
    mw, asst = _session()
    _run(asst, SNIPPET_MATERIALS)          # something for the stopped edit to touch
    mw.settings.setValue(K.RUN_TIMEOUT_KEY, 0.5)
    if abs(asst._kernel.run_timeout() - 0.5) > 1e-9:
        out.append("the kernel does not read the ai/run_timeout setting "
                   "(got %r)" % asst._kernel.run_timeout())
    started = time.time()
    result = _run(asst, "import time\n"
                        "slope_data['materials'][0]['c'] = 999.0\n"
                        "time.sleep(30)")
    elapsed = time.time() - started
    if elapsed > 20:
        out.append("the snippet ran %.1fs against a 0.5s limit" % elapsed)
    if "TimeoutError" not in result:
        out.append("the model is not told the snippet was stopped: %s"
                   % result[-300:])
    if float(mw.doc.slope_data["materials"][0].get("c") or 0) == 999.0:
        out.append("the stopped snippet's edit was kept; a snippet that raises "
                   "is rolled back whole")

    # The other direction: with the limit off, the same shape of snippet runs to
    # the end. A limit that fires on everything measures nothing.
    mw.settings.setValue(K.RUN_TIMEOUT_KEY, 0)
    ran = _run(asst, "import time\n"
                     "time.sleep(1.2)\n"
                     "print('FINISHED')")
    if "FINISHED" not in ran:
        out.append("with the limit off, a 1.2s snippet did not finish: %s"
                   % ran[-200:])
    mw.settings.remove(K.RUN_TIMEOUT_KEY)

    # The default, and the per-process override the producer harness uses.
    if asst._kernel.run_timeout() != K.DEFAULT_RUN_TIMEOUT_S:
        out.append("with nothing set, the limit is %r rather than the %gs default"
                   % (asst._kernel.run_timeout(), K.DEFAULT_RUN_TIMEOUT_S))
    os.environ[K.RUN_TIMEOUT_ENV] = "12"
    try:
        if asst._kernel.run_timeout() != 12.0:
            out.append("%s does not override the setting" % K.RUN_TIMEOUT_ENV)
    finally:
        os.environ.pop(K.RUN_TIMEOUT_ENV, None)
    mw.deleteLater()
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
    ("G. a secondary dead circle says whose loss it is",
     check_secondary_circle_wording),
    ("E. read-only and failed runs cost nothing", check_read_only_costs_nothing),
    ("F. mutation: injection disabled", check_mutation_disables_the_checks),
    ("H. a finding is the rule plus the row", check_delta_identity_is_the_row),
    ("H. a new finding arrives in full", check_new_finding_is_reported_in_full),
    ("H. an unchanged finding collapses", check_unchanged_finding_collapses),
    ("H. an ERROR never collapses", check_error_never_collapses),
    ("H. an ERROR outranks the quote cap", check_error_survives_the_quote_cap),
    ("H. a changed finding un-collapses", check_changed_finding_uncollapses),
    ("H. a fresh session reports in full", check_fresh_session_reports_in_full),
    ("H. the staged label survives collapse", check_staged_collapse_keeps_the_label),
    ("H. nothing is dropped unnamed", check_nothing_is_dropped_unnamed),
    ("H. every finding is quoted, not just counted",
     check_every_finding_is_quoted_once),
    ("H. two faults on one row stay two", check_two_faults_on_one_row),
    ("H. mutation: the delta's seven seams", check_delta_mutations),
    ("I. every helper the prompt names is callable", check_helpers_are_callable),
    ("I. every offered model can read an image",
     check_every_offered_model_can_see),
    ("I. Kimi's key is stored as kimi_api_key", check_the_keychain_name),
    ("I. the token meter accumulates", check_usage_accumulates),
    ("J. the Studio brief ships and fits", check_the_brief_ships),
    ("J. the skill body is out of the prompt",
     check_skill_body_is_not_in_the_prompt),
    ("K. the model summary is injected", check_model_summary_is_injected),
    ("K. the model summary says what the model is",
     check_model_summary_says_what_the_model_is),
    ("K. the corpus pointer resolves", check_corpus_index_returns_rows),
    ("L. the LEM run is the session's run", check_run_lem_is_the_sessions_run),
    ("L. an unqualified run runs the model's method",
     check_run_lem_runs_the_models_method),
    ("L. the result names its own surface", check_run_lem_returns_the_surface),
    ("L. a geometry edit on the wrong source is named",
     check_polygon_edit_on_a_profile_model_is_named),
    ("L. a discarded polygon edit cannot read as applied",
     check_discarded_polygon_edit_cannot_read_as_applied),
    ("L. the chat renders markdown", check_chat_renders_markdown),
    ("L. the chat degrades LaTeX to text", check_chat_degrades_latex),
    ("M. a prose fence is not code", check_prose_fence_is_not_code),
    ("M. the fallback is off where tools exist",
     check_the_fallback_is_off_where_tools_exist),
    ("M. the loop honours the gate", check_the_loop_honours_the_gate),
    ("M. a report lands in the Files folder",
     check_report_lands_in_the_files_folder),
    ("N. the transient march is a helper", check_transient_helper_marches),
    ("N. the transient route is in the brief",
     check_transient_helper_is_told_to_the_model),
    ("N. a runaway snippet is stopped", check_a_runaway_snippet_is_stopped),
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
