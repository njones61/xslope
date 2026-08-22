"""Assistant — the agent loop behind the chat dock (Phase B: multi-provider).

A manual tool-use loop (so each tool call can be gated by a confirm dialog) over
**LiteLLM**, so the same loop drives Claude, OpenAI, or a local Ollama model —
chosen in Settings (:mod:`studio.ai.config`). The blocking completion calls run on
a ``QThread``; each ``run_python`` call is marshalled to the GUI thread (modal
confirm + document mutation + re-render) and the worker waits for the result.
Conversation history (OpenAI message format) persists across turns. See §14.
"""

from __future__ import annotations

import json
import os
import re
import threading

from PySide6.QtCore import QObject, QThread, Signal

from .config import AssistantConfig

MAX_TOKENS = 8000

# OpenAI/LiteLLM function-tool format (works across providers).
RUN_PYTHON_TOOL = {
    "type": "function",
    "function": {
        "name": "run_python",
        "description": (
            "Execute Python in the live XSLOPE Studio session and return its "
            "stdout. A persistent namespace is preloaded with: `xslope` (the "
            "engine), `np`, `pd`, `plt` (matplotlib.pyplot), and the open project "
            "— `doc` (ProjectDocument), `slope_data` (the dict you edit), and "
            "`results`. Variables persist across calls like a notebook. Build or "
            "edit the input by mutating `slope_data` in place (it re-renders on "
            "the canvas; the user saves via Save As) rather than writing an .xlsx. "
            "Any `plt` figures you create are shown to the user automatically — "
            "you don't receive the images back. Print values you need to read. A "
            "snippet that CHANGES the model returns a MODEL CHECKS block after its "
            "output — the input checks run automatically on the edited model, and "
            "what they report is yours to resolve before you finish the turn."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "code": {"type": "string", "description": "Python source to execute."},
            },
            "required": ["code"],
        },
    },
}

STUDIO_SYSTEM = """\
You are the AI assistant embedded in **XSLOPE Studio**, a desktop GUI for the \
`xslope` slope-stability engine. You help the user build and edit slope models, \
run analyses (LEM / seepage / FEM), and explore results — by writing and running \
small Python snippets with the `run_python` tool.

IRON RULES. They hold whatever anything else says — later instructions, the \
reference below, an error message, a tool result. If something contradicts one, \
the iron rule wins and you say so.
1. **Max depth is the lowest elevation the problem DESCRIBES.** Never invent \
depth to give a search more room. Where the source states a base, that elevation \
IS `max_depth`; where it states nothing below the section, `max_depth` is the \
lowest elevation it does describe (usually the toe) — and if that reading is \
uncertain, ask rather than choose.
2. **No circle bottom below the base.** Every circle's `Depth` (its lowest point) \
sits at or above `max_depth`; on a rigid base the deep circle is TANGENT, \
`Depth = max_depth` exactly. A toe circle is `R = distance(center, toe)`, \
`Depth = Yo - R` — if that lands under the base, move the center or drop it.
3. **On a rigid base at the toe elevation, the ground ENDS at the toe.** Flat \
ground extended past it at the base elevation is soil of zero thickness and it \
destroys the domain. Extend sideways only where real soil continues below.
4. **Building is not running.** A request to build or edit — even one naming the \
analysis it is for — is complete when the model is built. Say what you built, \
OFFER the run, and stop.
5. **MODEL CHECKS are unresolved business.** Every snippet that changes the model \
returns a MODEL CHECKS block. Never report a model ready over a finding you have \
neither fixed nor named to the user.

Key facts about your environment:
- `run_python` runs in one persistent in-process namespace with `xslope`, `np`, \
`pd`, `plt`, and the live project (`doc`, `slope_data`, `results`) preloaded.
- To build or change inputs, mutate `slope_data` in place. The canvas re-renders \
automatically and the user persists the result with Save As — do NOT write .xlsx \
files for the build case.
- Figures you make with `plt` are shown to the user; you won't see them back, so \
print any numbers you need to reason about. Files you save (plots, CSVs, exported \
.xlsx, …) go to the working directory — the assistant output folder (also \
`OUTPUT_DIR`) — and the user can open them from the chat or the dock's Files \
button, so save them with clear, descriptive filenames when the user asks for a \
file or a saved plot.
- Prefer one focused snippet per step. Keep code short and readable. Print a \
brief result. Don't reformat or refactor the user's data beyond the request.
- **Only do what was asked.** If asked to edit inputs, just edit `slope_data` — \
do NOT run an analysis afterward unless the user also asked to run/solve.
- **Edits apply immediately and persist on success.** Put a mutating edit in its \
own snippet, separate from analysis. NEVER re-run a mutating snippet — re-running \
`x += 5` doubles the effect. A relative change (`+=`, append, scale) must run \
exactly once. If a snippet fails, fix that snippet; do not repeat an edit that \
already ran. If unsure whether an edit applied, print the current value first.
- When unsure of a key or signature, inspect at runtime instead of guessing: \
`print(sorted(slope_data))`, `print(slope_data['materials'][0])`, \
`import xslope; print([n for n in dir(xslope.solve)])`, `help(fn)`.
- To run a limit-equilibrium analysis, prefer the preloaded helper \
`run_lem(method='bishop')` (methods: oms, bishop, janbu, spencer, corps, lowe, \
mprice). It handles the loaded project's failure surface, returns the result \
dict (with 'FS'), and shows the solution plot — don't rebuild that pipeline by \
hand.

For slope-stability modeling rules (geometry extents, starting circles, ponded \
water, diagram reading, etc.) FOLLOW the xslope skill reference — it is the \
authoritative source. The compact reminders below restate its key rules; the \
skill has the full detail.
"""

# Appended only for Anthropic, where prompt caching makes the large skill body
# cheap. Local/other models get the compact prompt above and introspect at
# runtime (the skill body measures ~33k tokens today; sending that every turn
# makes a local model crawl).
_SKILL_HEADER = ("\n\n---\n\nReference — the `slope_data` schema and engine API "
                 "(ground truth for keys and signatures; it builds the same "
                 "in-memory `slope_data` dict you edit here, so reuse its schema "
                 "and geometry knowledge — only skip its final file-save step):\n\n")

# Appended AFTER the skill body so it is the last thing the model reads. The
# standalone skill builds the in-memory `slope_data` dict (exactly what we want)
# and then persists it with save_slope_data_to_xlsx. Inside Studio the project is
# already open and the user saves via Save As, so we override only that final
# file-save — the dict-building core carries over unchanged.
_SKILL_TRAILER = """\

---

CRITICAL — you are inside XSLOPE Studio, NOT the standalone file-based skill. The \
reference above builds the in-memory `slope_data` dict (which is exactly right) and \
then SAVES it to an .xlsx. Reuse everything about constructing the dict; override \
only the save:
- There is ALREADY an open project. Build and edit it by mutating the in-memory \
`slope_data` dict in `run_python`. The canvas re-renders automatically.
- NEVER create, open, write, or save an .xlsx (or any) file. Do not call \
`save_slope_data_to_xlsx`, `pd.ExcelWriter`, `openpyxl`, `load_slope_data`, or \
read/write the filesystem. The user persists the model themselves via Save As.
- When the reference (or your own narration) says "build the input file," that \
means populate `slope_data` in memory — say "I'll build the model" and mutate the \
dict, never a file.
- Construct geometry as the in-memory schema expects (e.g. `profile_lines` are \
dicts with `mat_id` + shapely-ready coords, `materials` are dicts). Inspect an \
existing record first when one is present; for an empty project, print \
`sorted(slope_data)` and build the lists directly.
"""

# Authoritative in-memory record schemas. These ARE the ground truth for field
# names, so the model never loads reference .xlsx files to reverse-engineer them
# (slow, and reads the filesystem needlessly). Included on every path — compact
# enough for local models, and it overrides the file-oriented skill.
SCHEMA_BRIEF = """\

Authoritative `slope_data` record schemas — these key names are the ground truth.
Do NOT load or read .xlsx files (or call `load_slope_data`) to discover the schema;
the records below ARE the schema. A project is ALWAYS open and `slope_data` is a
live dict (every list is `[]` for a fresh one) — just build the lists directly.
Do NOT call `doc.new()` or reassign `doc.slope_data`; your edits are detected and
the canvas re-renders automatically.

- materials[i]: name, gamma, option ('mc' Mohr-Coulomb | 'cp' | 'elastic' — infinite strength,
  cannot fail; ignores every strength key below, uses only gamma/gsat/E/nu (+ seepage); LEM
  treats it as impenetrable), c, phi, cp, r_elev, d, psi, t_cut (Rankine tensile-strength
  cutoff, stress units; None/blank = no cutoff (pre-v16 default), 0 = no tension; FEM only,
  LEM ignores it), phi_b, s_cap (v17 matric-suction apparent cohesion, LEM & FEM — Fredlund
  extended MC (in the FEM/SSRM the suction term is reduced by F alongside tan phi');
  phi_b = unsaturated friction angle phi^b in degrees, None/blank = no suction
  credit (default); s_cap = cap on the credited suction, stress units, None/blank = uncapped.
  Active only for option in {mc,pow,hb} with u in {piezo,seep}; leave both None for cp/elastic
  or u in {none,ru}. CAUTION: with u='piezo' the hydrostatic suction above the line is
  unbounded, so s_cap is essential; with u='seep' the FE field self-bounds it), u
  ('none'|'piezo'|'seep'), sigma_gamma/c/phi/cp/d/psi (stddevs, 0 if unused), k1, k2, alpha,
  kr0, h0 (seepage), E, nu (FEM; also the sole mechanical properties when option='elastic').
  Minimal MC example:
  {'name':'Clay','gamma':130.0,'option':'mc','c':400.0,'phi':0.0,'t_cut':None,
   'phi_b':None,'s_cap':None,'u':'none',
   'cp':0.0,'r_elev':0.0,'d':0,'psi':0,'sigma_gamma':0.0,'sigma_c':0.0,
   'sigma_phi':0.0,'sigma_cp':0.0,'sigma_d':0.0,'sigma_psi':0.0,'k1':0.0,'k2':0.0,
   'alpha':0.0,'kr0':0.0,'h0':0.0,'E':0.0,'nu':0.0}
- profile_lines[i]: {'coords':[(x,y),...], 'mat_id':0}  # mat_id = 0-based index into
  materials; lines are layer-top boundaries, ordered top to bottom.
- polygons[i]: {'polygon': <shapely Polygon>, 'mat_id':0}  # zoned geometry (dams,
  lenses); use INSTEAD of profile_lines, not both. Build the Polygon from coords
  (`from shapely.geometry import Polygon`). The canvas rebuilds ground_surface and
  domain_polygon from the polygons automatically — do NOT call plot_inputs,
  build_ground_surface_*, or set ground_surface/domain_polygon yourself.
- circles[i]: {'Xo':20.0,'Yo':40.0,'Depth':-10.0,'R':50.0}  # Depth = elevation of
  the circle's lowest point, R = Yo - Depth. In-memory there is no intercept key,
  so a TOE circle (one passing THROUGH the toe point) is R = distance(center, toe),
  Depth = Yo - R. Do NOT use Depth = toe_elevation for the toe circle — that is
  merely tangent to the toe LEVEL (lowest point below the center), not through the
  toe point, and is a different circle. That construction is also where iron rule 2
  bites: if the resulting Depth falls below max_depth, the circle is not usable —
  move the center (or drop the candidate), never lower max_depth to admit it.
- non_circ[i]: {'X':-10.0,'Y':0.0,'Movement':'Free'}
- piezo_line / piezo_line2: list of (x, y) tuples.
- dloads / dloads2: list of blocks; each block is a list of {'X','Y','Normal'} pts.
- reinforcement_lines[i]: {'x1','y1','x2','y2','t_max','t_res','lp1','lp2','area','E',
  'tend1','tend2','spacing','type','dir','appl','adhesion','delta'}
  ('adhesion'+'delta' both set = overburden pullout law; 'lp1'/'lp2' then unused)
  # EDIT THIS one; reinforce_lines (capitalized X/Y/T/Tres) is derived from it.
- pile_lines[i]: {'x1','y1','x2','y2','D_pile','S','E','I','area','M_cap','V_cap',
  'theta_p','fixity','label','H'}
- scalars: gamma_water, max_depth (elevation of the hard base — the lowest
  elevation the PROBLEM describes, never one you chose; see iron rule 1),
  k_seismic, tcrack_depth, tcrack_water, circular (bool).

Editing the source lists re-renders the canvas — derived geometry (ground_surface,
polygons, domain_polygon) is rebuilt automatically AFTER the snippet returns, so a
one-off edit needs no plot_inputs call. But that auto-rebuild does NOT run between
iterations WITHIN a snippet: if you vary geometry in a loop, call
`resync_geometry()` after each edit before analyzing — or use `run_lem(...)`, which
resyncs for you. Otherwise every iteration analyzes the STALE original geometry and
you get a constant/flat result; if a sweep looks suspiciously constant, suspect
stale derived geometry FIRST before reaching for a physics explanation.

For a sensitivity study ("FS vs <parameter>" — slope angle, a strength, a water
level, …) use the preloaded `sensitivity(values, apply, param=...)` helper: it
sweeps the parameter with NO per-step plot, computes the critical FS each step,
writes a summary CSV + ONE plot to the output folder, and restores the project.
`apply(v)` is your callback that edits slope_data to an absolute value (e.g. moves
the crest). Don't hand-roll the loop with a solution plot per step. For a single
analysis, use `run_lem(method=...)` (pass plot=False to suppress its plot).

For a DESIGN question ("what c gives FS = 1.5?", "vary Su between X and Y and
show where FS hits the target") use the preloaded `design_sweep(param, low,
high, steps=..., target_fs=..., mode='lem', method=...)`: param takes
`{'material': name_or_index, 'property': field}` or `{'global': 'k_seismic'}`;
it sweeps, plots the OUTPUT-vs-value curve with the target line (axes auto-
labelled), and returns the result dict — read `crossing` (the interpolated
required value) only when `bracketed` is True; on a miss report `fs_range` and
extend the way `extend` says, never extrapolate. Discover sweepable parameters
with `list_params(mode=...)`. For material properties prefer these over the
callback-style `sensitivity()`.
`mode` picks the engine and the OUTPUT quantity `target_fs` names: 'lem'
(output = FS, default), 'fem' (output = FS from a full SSRM solve — needs a
finite-element mesh in slope_data['mesh'] and costs MINUTES PER STEP, so tell
the user to expect a wait, keep `steps` tiny (2-3), and pass `fem_opts` for SSRM
knobs), and 'seep' (output = total discharge q, so `target_fs` is a target q
like 6e-6 — also needs a mesh; `seep_opts={'bc':1}` picks the BC set). A seep
design study varies a hydraulic conductivity (`seep:<mat>:k1`) or a reservoir-
head boundary (`seep_bc:<set>:<head_index>`, or the dict
`{'seep_bc': {'set': 1, 'head_index': 0}}`) and finds the value giving the
target q — the classic "reservoir level vs seepage discharge" study. Both mesh
modes error with "build a mesh first" if none is present.
"""

# Compact modeling rules — appended ONLY when the full skill is NOT loaded (i.e.
# for local/non-Anthropic models). On the Anthropic path the skill already carries
# these in full, so mirroring them here would just duplicate cached tokens.
MODELING_BRIEF = """\

Modeling rules (slope-stability physics):
- Starting circles: put Xo at the toe-crest midpoint and Yo = toe_elev + 2x slope
  height; ALWAYS include one circle through the toe and one tangent to the base of
  EACH material layer (including the hard base at max_depth). A toe circle PASSES
  THROUGH the toe point: R = distance(center, toe), Depth = Yo - R — NOT Depth =
  toe_elevation (that is only tangent to the toe level, a different circle). Then
  apply the FLOOR RULE, the same one the generator applies to its own candidates:
  a circle whose Depth falls below max_depth is DROPPED (or its center moved and
  R recomputed until Depth >= max_depth) — never kept and never cured by lowering
  max_depth. On a rigid base the deepest circle is tangent to it, Depth =
  max_depth exactly. `generate_starting_circles(slope_data)` (from
  xslope.search / xslope.generators) implements the whole strategy including this
  filter — prefer it over hand-building.
- Extent: extend the flat ground far enough on BOTH sides that every trial circle
  daylights on the ground INSIDE the model, never at a vertical edge (>= ~2x slope
  height beyond toe and crest, more for deep base circles). Do NOT copy the source
  diagram's width — it is usually cropped. The hard base spans the full width of
  the ground that exists. RIGID-BASE COROLLARY: where the base sits AT the toe
  elevation, the ground surface ENDS at the toe — extending flat ground beyond it
  at the base elevation is soil of zero thickness, which makes the domain
  degenerate and every analysis fail on a geometry error. The cure is to end the
  profile at the toe, never to lower max_depth (that invents depth the problem
  does not describe). No search room is lost: a circle tangent to the base
  daylights at or above the toe. Extend sideways only where real soil continues
  below the ground being extended.
- Ponded/standing/reservoir water above the ground surface is ALWAYS a hydrostatic
  dloads load: Normal = gamma_water * (water_surface_elev - y_ground) at each ground
  point. Apply it over the ENTIRE submerged surface — flat foundation/bench areas
  AND sloping faces, following the ground profile — not just the slope face. It is
  SEPARATE from any phreatic/piezo line (pore pressure): a reservoir on a dam needs
  BOTH the upstream surface-water load and the internal phreatic line. Apply even
  for total-stress phi=0 (`u='none'`); never skip it. A water table is the inverted-
  triangle (▽) symbol (may sit above the crest = submerged); ask if its level/extent
  is unclear rather than guessing.
- Reading sketches: attribute every dimension arrow to the right feature — a dimension
  near a water-table line usually measures the LAYER thickness, not the WT depth
  (cross-check that thicknesses sum to the section depth); if the WT elevation stays
  ambiguous, ask. Reinforcement drawn as "N layers spaced s vertically" means bottom
  layer AT the toe/base elevation (y = 0, s, ..., (N-1)s), each line starting ON the
  slope face with length = the labeled dimension from the face; ask if unclear.
- Method expectations: for phi=0 soils, OMS = Bishop ~= Spencer (identical FS is
  normal, not a bug). OMS is unreliable on submerged/high-pore-pressure slopes and
  most prone to a different search minimum — trust Spencer/Bishop there.
- FEM/SSRM domains follow the same extent rule as LEM (flats >= ~2x slope height,
  plus foundation depth below the toe); a domain cropped at the toe inflates FS.
"""


TOOL_NAMES = {"run_python"}
MAX_STEPS = 25                  # safety cap on agentic iterations per turn


def _parse_text_tool_call(content, names):
    """Weaker (often local) models emit a tool call as plain-text JSON instead of
    a structured tool_call. Recover ``(name, args)`` from such content, or None."""
    import re
    if not content:
        return None
    s = content.strip()
    candidates = []
    fence = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", s, re.S)
    if fence:
        candidates.append(fence.group(1))
    candidates.append(s)
    brace = re.search(r"\{.*\}", s, re.S)
    if brace:
        candidates.append(brace.group(0))
    for cand in candidates:
        try:
            obj = json.loads(cand)
        except Exception:
            continue
        if not isinstance(obj, dict):
            continue
        name = obj.get("name") or obj.get("tool")
        args = obj.get("arguments", obj.get("parameters", obj.get("input", {})))
        if isinstance(args, str):
            try:
                args = json.loads(args)
            except Exception:
                args = {"code": args}
        if name in names and isinstance(args, dict):
            return name, args
    return None


def _is_runnable(code):
    """True if ``code`` compiles (trying the literal-escape repair too)."""
    if not code:
        return False
    for variant in (code, code.replace("\\n", "\n").replace("\\t", "\t")):
        try:
            compile(variant, "<assistant>", "exec")
            return True
        except SyntaxError:
            pass
    return False


def _extract_code(content):
    """Recover ``run_python`` code from an assistant turn that made no structured
    tool call — a JSON tool-call, a fenced ```python block, or whole-content
    Python (weaker models often just *write* the code in chat). Returns
    ``(code, pure)`` where ``pure`` means the message is only the code (so it is
    suppressed from the transcript); ``(None, False)`` if nothing runnable."""
    import re
    if not content:
        return None, False
    s = content.strip()
    call = _parse_text_tool_call(s, TOOL_NAMES)
    if call:
        return call[1].get("code", ""), True
    fence = re.search(r"```(?:python|py)?\s*\n?(.*?)```", s, re.S)
    if fence and _is_runnable(fence.group(1).strip()):
        return fence.group(1).strip(), False        # keep any surrounding prose
    if _is_runnable(s) and re.search(
            r"\b(slope_data|doc|results|run_lem|xslope|plt|np|pd|print)\b", s):
        return s, True
    return None, False


_SOURCE_KEYS = (
    # `circular` is the stored surface FAMILY -- which of the two surfaces on the
    # deck a run analyses. Flipping it changes the answer as surely as moving a
    # circle does, so it belongs here: without it the flip left no undo step, no
    # dirty flag, and (once the checks were keyed off this same comparison) no
    # input check either.
    "gamma_water", "tcrack_depth", "tcrack_water", "k_seismic", "max_depth",
    "circular",
    "materials", "profile_lines", "polygons", "circles", "non_circ", "piezo_line",
    "piezo_line2", "dloads", "dloads2", "reinforcement_lines", "pile_lines",
    "seepage_bc", "seepage_bc2",
)

# Short tags for the undo-history dropdown, keyed by source input. Several keys map
# to the same tag (e.g. both piezo lines -> "piezo") so a multi-key edit reads cleanly.
_KEY_LABELS = {
    "gamma_water": "global", "tcrack_depth": "global", "tcrack_water": "global",
    "k_seismic": "global", "max_depth": "global", "circular": "surface family",
    "materials": "materials", "profile_lines": "profile", "polygons": "polygons",
    "circles": "circles", "non_circ": "non-circular", "piezo_line": "piezo",
    "piezo_line2": "piezo", "dloads": "dloads", "dloads2": "dloads",
    "reinforcement_lines": "reinforcement", "pile_lines": "piles",
    "seepage_bc": "seep BC", "seepage_bc2": "seep BC",
}


# --- harness-enforced validation --------------------------------------------
# A rule the model is told to follow is a rule it follows MOST of the time, and a
# rule stated three levels deep in a prompt loses to a confident-sounding error
# message. So the rules that matter are not left to the prompt alone: every
# snippet that CHANGED the model gets its derived geometry rebuilt and the input
# checks run, and the findings come back appended to that snippet's own output —
# in the conversation, as the model's own tool result, where they cannot be missed
# and must be addressed before the turn ends (STUDIO_SYSTEM iron rule 5).
#
# It also closes the edit-cascade, which no prompt can: a one-field fix leaves its
# dependents stale, and the dependent is silent. Correct Max depth upward to the
# base the problem states, and the circles that were tangent to the old, deeper
# base now sit underneath the new one; the checks say so on the very snippet that
# moved it, instead of surfacing as a failed run several turns later.
MODEL_CHECKS_OPEN = "=== MODEL CHECKS ==="
MODEL_CHECKS_END = "=== END MODEL CHECKS ==="
MODEL_CHECKS_CLEAN = "=== MODEL CHECKS: clean ==="
#: Cap on findings quoted back. Preflight messages are whole paragraphs, so an
#: uncapped block on a half-built model could outweigh the snippet's own output.
MAX_CHECK_FINDINGS = 6


# --- what has already been said ---------------------------------------------
# A finding is quoted in full the first time a block has room for it. After that
# it is a standing fault, and a long build re-reads the same paragraphs after
# every edit — attention spent on findings that have not moved, and the new one
# arrives in the middle of the pile. So each block quotes what is NEW or CHANGED
# and names the rest by rule key on one line.
#
# What may collapse is bounded by one rule: a finding collapses only after it has
# been QUOTED, never merely after it has been counted. The block quotes at most
# MAX_CHECK_FINDINGS paragraphs, so on a model carrying more than that the rest
# are named on the overflow line with the command that prints them — and they
# stay unrecorded, so the next block quotes them instead. A model with twelve
# findings receives all twelve in full over two or three blocks rather than the
# same six forever. Findings never yet quoted take the quota first, which is what
# makes that rotation terminate.
#
# An edit-answerable ERROR is also never collapsed: it refuses the run until it
# is answered, so it is live business on every edit rather than history. A
# finding STAGED BY A RUN is an error too, and it does block the run, but no edit
# can answer it — so it collapses like a warning, onto a line that keeps saying
# it is staged.
#
# Keying is the whole difficulty. A finding is "the same one" when it is the same
# rule about the same row — NOT when its message text matches, because a message
# carries values ("Circle 3 has Depth = -20.4") that drift under unrelated edits
# and would un-collapse a fault that never changed. So identity is the rule id
# plus the row the message leads with (preflight's own convention: "Circle 2",
# "Material 3 ('Core')"), and what counts as a CHANGE is the severity plus the
# message with its numbers blanked — a reworded value collapses, a materially
# different message or a new severity does not.

#: Words that can sit in front of a number without the number being a row index.
#: "None of the 5 reinforcement lines" is not a subject called "None of the 5".
_NOT_A_SUBJECT = {
    "a", "all", "an", "and", "any", "are", "at", "both", "by", "for", "in", "is",
    "no", "not", "of", "on", "only", "or", "than", "the", "to", "up", "was",
    "were", "with",
}

#: A leading row label: up to three words then an index, optionally followed by
#: the row's own name in parentheses. The index must be a whole number (a decimal
#: is a value, not a row).
_SUBJECT_RE = re.compile(
    r"^([A-Z][\w'’-]*(?: [\w'’-]+){0,2} \d+)(?![\d.])( \('[^']*'\))?")

#: Values, for blanking. A digit glued to letters is part of a FIELD NAME, not a
#: value: "Lp1 = 7.9" and "Lp2 = 7.9" are two different faults on one row, and
#: blanking both digits made them one finding — the second was never reported and
#: the model was pointed at the text of the first. Digit-suffixed field names are
#: everywhere in the schema (x1/y1, x2/y2, lp1/lp2, k1/k2), so the lookbehind is
#: load-bearing, not tidiness.
_NUMBER_RE = re.compile(r"(?<![A-Za-z])-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")


def _finding_subject(message):
    """The row a finding is about (``"Circle 2"``, ``"Material 3 ('Core')"``), or
    None when the message does not lead with one.

    One row per label, deliberately: a message that opens on two rows ("Circle 2
    and 3 have …") matches no label and falls back to the message skeleton, which
    is a correct key for a finding about a pair. No rule writes one today; this is
    what happens if one starts.
    """
    m = _SUBJECT_RE.match(message or "")
    if not m:
        return None
    words = m.group(1).split()
    if words[-2].lower() in _NOT_A_SUBJECT:
        return None
    return m.group(1) + (m.group(2) or "")


def _finding_skeleton(message):
    """The message with every number blanked and whitespace normalised — what is
    left when the values a finding quotes drift but the finding does not."""
    return _NUMBER_RE.sub("#", " ".join((message or "").split()))


def _finding_key(f):
    """Identity: the rule, and the row it is about. Falls back to the message
    skeleton when a message names no row, so two findings of one rule stay
    distinct without keying on values that drift."""
    subject = _finding_subject(f.message)
    return f.rule_id, subject if subject is not None else _finding_skeleton(f.message)


def _finding_sig(f):
    """What a re-report has to differ in to count as changed."""
    return f.severity, _finding_skeleton(f.message)


def _never_collapses(f):
    """True for a finding that is quoted in full in every block that has room.

    The test is whether an EDIT can answer it, not whether it blocks a run. An
    ordinary ERROR is answerable by the edit being made, and it refuses the run
    until it is, so it stays live business on every edit rather than becoming
    history. A finding STAGED BY A RUN is an error and refuses the run too, but
    no edit resolves it — its cure is the analysis it names (see
    :data:`xslope.preflight.STAGED_BY_RUN`) — so it collapses like a warning
    once it has been quoted, onto a line that keeps the staged label.
    """
    from xslope.preflight import STAGED_BY_RUN
    return f.severity == "error" and f.rule_id not in STAGED_BY_RUN


class ChecksMemo:
    """What the last MODEL CHECKS block reported, so the next one can report the
    delta rather than the list again.

    One per session: :meth:`Assistant.reset` clears it for a new conversation and
    the document's ``loaded`` signal clears it for a new (or reopened) project, so
    tracking never carries findings across into a model they were not about.
    """

    def __init__(self):
        self._seen = {}

    def reset(self):
        self._seen = {}

    @staticmethod
    def keys(findings):
        """One key per finding, in order. Repeats of a single key (one rule, one
        row, twice) are numbered so they stay distinct."""
        out, counts = [], {}
        for f in findings:
            key = _finding_key(f)
            counts[key] = counts.get(key, 0) + 1
            out.append(key + (counts[key],))
        return out

    def delta(self, findings):
        """``(keys, fresh)`` — every finding's key, and the keys of the ones this
        session has not been given in full, or has been given differently.
        Records nothing: what a block actually quotes is decided after this."""
        keys = self.keys(findings)
        fresh = {k for k, f in zip(keys, findings)
                 if self._seen.get(k) != _finding_sig(f)}
        return keys, fresh

    def record(self, keys, findings, quoted):
        """Remember the block that was just written.

        Only the findings in ``quoted`` (a set of indices) become history — a
        finding the block merely counted was never given to the model, so it must
        stay eligible to be quoted in a later block. A finding that was quoted in
        an earlier block and only counted in this one keeps the entry it already
        had, and one that has disappeared from the model loses its entry.
        """
        now = {}
        for i, (k, f) in enumerate(zip(keys, findings)):
            if i in quoted:
                now[k] = _finding_sig(f)
            elif k in self._seen:
                now[k] = self._seen[k]
        self._seen = now

    def quoted_before(self, key):
        """True when this session has already been given that finding in full."""
        return key in self._seen


def _quote_order(rows, keys, memo):
    """How a block spends its quote cap: a sort key over row indices, two tiers.

    An edit-answerable error takes the quota first, whatever else is waiting. It
    refuses the run, and a block that pushed it onto the overflow line to make
    room for a batch of new warnings would be the one block the model reads while
    the run stays refused.

    Below that, findings this session has never been given go before findings it
    has. That is what makes a model carrying more findings than the cap work
    through them over successive blocks instead of re-quoting the same few and
    counting the tail as "reported in full earlier".
    """
    return lambda i: (not _never_collapses(rows[i]), memo.quoted_before(keys[i]))


def _rule_key_list(findings):
    """``"mat.option_missing ×2, main.seismic_missing"`` — every finding named by
    its rule key, so nothing is dropped without being named."""
    counts = {}
    for f in findings:
        counts[f.rule_id] = counts.get(f.rule_id, 0) + 1
    return ", ".join(k + (f" ×{n}" if n > 1 else "") for k, n in counts.items())


def _unchanged_line(findings, standing):
    """The one line a block spends on findings it has already reported in full."""
    n = len(findings)
    noun = "finding" if n == 1 else "findings"
    return (f"{n} earlier {noun} unchanged, {standing} (reported in full "
            f"earlier): {_rule_key_list(findings)}")


def _checks_selection(sd):
    """The preflight ``selection`` for the checks: which surface family the model
    would run as. Stating it suppresses the ambiguity warning on a deck that
    carries both, which is a question for a run dialog and not for an edit."""
    has_c, has_n = bool(sd.get("circles")), bool(sd.get("non_circ"))
    if has_c and has_n:
        return {"surface": "circular" if sd.get("circular", True) else "noncircular"}
    if has_n:
        return {"surface": "noncircular"}
    if has_c:
        return {"surface": "circular"}
    return {}


def model_checks_text(slope_data, memo=None):
    """The MODEL CHECKS block for a model that a snippet just changed.

    Rebuilds the derived geometry first (so the checks read the domain the edit
    actually produced, not the one before it) and then runs :func:`xslope.preflight
    .preflight` for a limit-equilibrium run — the analysis every Studio model has
    to be valid as, and the one whose rules cover geometry, materials, water and
    surfaces alike. Errors AND warnings are reported: a warning here is the shape
    the three live failures took, not a formality.

    ``memo`` is the session's :class:`ChecksMemo`. With one, a finding already
    reported in full and unchanged since collapses to a named entry on one line
    (errors excepted — they never collapse); without one, every finding is quoted
    in full, which is what a first report is anyway.

    Returns the block as text — never raises, and never an empty string, because
    silence would read as "clean" to the model.
    """
    sd = slope_data
    if not isinstance(sd, dict):
        return ""
    if memo is None:
        memo = ChecksMemo()         # no session tracking: everything is new
    try:
        from studio.editors import _resync_geometry
        _resync_geometry(sd)
    except Exception:
        pass                    # a half-built model may have nothing to resync yet
    # Mid-build the model is legitimately incomplete, and reporting every rule it
    # does not satisfy yet would train the reader to skim the block. Wait until
    # there is a model to check.
    missing = [w for w, got in (("materials", sd.get("materials")),
                                ("geometry", sd.get("profile_lines")
                                 or sd.get("polygons"))) if not got]
    if missing:
        # Nothing was reported, so nothing is history: the next block that does
        # report starts from a clean sheet and quotes everything in full.
        memo.reset()
        return (f"{MODEL_CHECKS_OPEN}\nnot run: the model has no "
                f"{' or '.join(missing)} yet. The checks run by themselves as soon "
                f"as it does.\n{MODEL_CHECKS_END}")
    selection = _checks_selection(sd)
    try:
        from xslope.preflight import STAGED_BY_RUN, preflight
        skip = None
        if sd.get("mesh") is not None and not (sd.get("circles") or sd.get("non_circ")):
            # A finite-element model's failure surface is an OUTPUT, so "no failure
            # surface" is not a defect in one.
            skip = ["surface.none_defined"]
        report = preflight(sd, "lem", selection, skip=skip)
        rows = list(report.errors) + list(report.warnings)
    except Exception as exc:
        memo.reset()
        return (f"{MODEL_CHECKS_OPEN}\ncould not be evaluated on this model "
                f"({type(exc).__name__}: {exc}). Check the model yourself before "
                f"reporting it ready.\n{MODEL_CHECKS_END}")
    if not rows:
        memo.reset()
        return MODEL_CHECKS_CLEAN
    # Split off the findings an EDIT cannot answer (xslope.preflight.STAGED_BY_RUN).
    # Reporting "a material takes pore pressure from a seepage solution and there
    # isn't one" as a plain error to fix, on a turn where the iron rules also forbid
    # running anything, leaves exactly one apparent way out: change the material's
    # pore-pressure option — which does not supply the missing field, it deletes the
    # requirement and silently analyses a different problem. So they are labelled
    # for what they are and routed to the only honest response: offer the run.
    # What this block wants to quote in full: whatever is new, changed, or not yet
    # quoted, plus every edit-answerable error. Kept as ROW INDICES rather than
    # findings, so two findings that read alike are still two.
    keys, fresh = memo.delta(rows)
    idx = range(len(rows))
    wanted = [i for i in idx if keys[i] in fresh or _never_collapses(rows[i])]
    wanted.sort(key=_quote_order(rows, keys, memo))
    staged_wanted = [i for i in wanted if rows[i].rule_id in STAGED_BY_RUN]
    fault_wanted = [i for i in wanted if rows[i].rule_id not in STAGED_BY_RUN]
    shown_i = sorted(fault_wanted[:MAX_CHECK_FINDINGS])
    over_i = sorted(fault_wanted[MAX_CHECK_FINDINGS:])
    memo.record(keys, rows, set(shown_i) | set(staged_wanted))

    staged = [i for i in idx if rows[i].rule_id in STAGED_BY_RUN]
    faults = [i for i in idx if rows[i].rule_id not in STAGED_BY_RUN]
    sel_arg = f", {selection!r}" if selection else ""
    parts = [MODEL_CHECKS_OPEN]
    if faults:
        old_faults = [rows[i] for i in faults if i not in fault_wanted]
        if shown_i:
            parts.append(f"The input checks found {len(faults)} problem(s) in the "
                         f"model as it now stands. Fix them, or put them to the user "
                         f"with your reason for leaving them — do not report this "
                         f"model ready over them.")
            parts += [f"  {rows[i].severity.upper()} [{rows[i].rule_id}] "
                      f"{rows[i].message}" for i in shown_i]
        if over_i:
            # Named, not quoted — so the command that prints them rides on this
            # line every time it appears, and none of them is recorded as told.
            rest = [rows[i] for i in over_i]
            parts.append(f"  (+{len(rest)} more, not quoted here — "
                         f"{_rule_key_list(rest)} — read them with "
                         f"`from xslope.preflight import preflight; "
                         f"print(preflight(slope_data, 'lem'{sel_arg}).format())`)")
        if old_faults:
            parts.append(_unchanged_line(old_faults, "still unresolved"))
    if staged:
        new_staged = [rows[i] for i in staged_wanted]
        old_staged = [rows[i] for i in staged if i not in staged_wanted]
        if new_staged:
            parts.append(
                f"STAGED BY A RUN, not by an edit ({len(staged)}). These name an "
                f"analysis that has not been run yet. Do NOT edit the model to "
                f"silence one — changing the input they ask about changes the "
                f"physics being analysed. Tell the user what is outstanding and "
                f"offer to run it.")
        for f in new_staged:
            # The rule's own message is written for a user at the sheet, where
            # picking a different input IS one of their options. It is not one of
            # yours, so say whose it is.
            parts.append(f"  [{f.rule_id}] {f.message} — resolved by running "
                         f"{STAGED_BY_RUN[f.rule_id]}, not by editing. Where that "
                         f"message offers a different input instead, that is the "
                         f"user's call to make, not yours.")
        if old_staged:
            # The label survives the collapse: what makes these different from a
            # fault to fix is the whole point of reporting them separately, and a
            # one-line reminder that is not still saying so would read as a list
            # of edits outstanding.
            parts.append(_unchanged_line(
                old_staged, "still STAGED BY A RUN, not by an edit"))
    parts.append(MODEL_CHECKS_END)
    return "\n".join(parts)


def _clean(o):
    if hasattr(o, "wkt"):
        return o.wkt
    try:
        import numpy as np
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.generic):
            return o.item()
    except Exception:
        pass
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(x) for x in o]
    return o


def _source_key_sigs(sd):
    """Per-key JSON signatures of the editable source inputs, to tell whether (and
    which part of) a code run changed the inputs vs a read-only query. Returns None
    if it can't be built (then the caller treats the run as a possible edit)."""
    try:
        return {k: json.dumps(_clean(sd.get(k)), sort_keys=True, default=str)
                for k in _SOURCE_KEYS}
    except Exception:
        return None


def _assistant_edit_label(changed_keys):
    """A short 'Assistant: …' label from the source keys an edit touched."""
    tags = []
    for k in changed_keys:
        tag = _KEY_LABELS.get(k, k)
        if tag not in tags:
            tags.append(tag)
    if not tags:
        return "Assistant edit"
    shown = ", ".join(tags[:3]) + ("…" if len(tags) > 3 else "")
    return f"Assistant: {shown}"


def _load_skill_text():
    """The /xslope skill body (schema + API knowledge), best-effort. The docs file
    is the editable master (fresh in a repo checkout); the copy shipped as package
    data (xslope/resources/xslope_skill.md) covers pip installs where docs/ is
    absent. A run_tests.py sync check keeps the two identical."""
    import os
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    path = os.path.join(here, "docs", "usage", "claude", "xslope.md")
    try:
        with open(path, encoding="utf-8") as f:
            return f.read()
    except Exception:
        pass
    try:
        from importlib import resources
        return (resources.files("xslope") / "resources" / "xslope_skill.md").read_text(
            encoding="utf-8")
    except Exception:
        return ""


class _AgentWorker(QThread):
    """Runs the LiteLLM tool-use loop off the GUI thread."""

    text = Signal(str)              # an assistant text block
    tool_call = Signal(object)      # {id, name, input, holder} — handled on GUI thread
    failed = Signal(str)
    done = Signal()

    def __init__(self, kwargs, system, messages, tools, cache_system=False, parent=None):
        super().__init__(parent)
        self._kwargs = kwargs       # provider/model kwargs for litellm.completion
        self._system = system
        self._messages = messages   # shared list (OpenAI format), mutated in place
        self._tools = tools
        self._cache_system = cache_system   # Anthropic prompt caching of the system
        self._cancel = threading.Event()

    def _system_message(self):
        # Anthropic prompt caching: a content block carrying cache_control. The
        # large skill prompt then bills at cache-read rates on repeat turns.
        # Plain string for other providers (a list-content system can 400 there).
        if self._cache_system:
            return {"role": "system", "content": [
                {"type": "text", "text": self._system,
                 "cache_control": {"type": "ephemeral"}}]}
        return {"role": "system", "content": self._system}

    def cancel(self):
        self._cancel.set()

    def run(self):
        try:
            import litellm
            litellm.drop_params = True          # ignore params a provider doesn't support
            litellm.suppress_debug_info = True
        except Exception:
            self.failed.emit("The 'litellm' package is not installed. "
                             "Install it with: pip install \"xslope[ai]\"")
            return
        try:
            for _ in range(MAX_STEPS):
                if self._cancel.is_set():
                    break
                resp = litellm.completion(
                    messages=[self._system_message()] + self._messages,
                    tools=self._tools, max_tokens=MAX_TOKENS, **self._kwargs)
                msg = resp.choices[0].message
                tool_calls = getattr(msg, "tool_calls", None) or []

                # Fallback for models that don't make a structured tool call but
                # write the code (as JSON, a fenced block, or bare Python).
                code_call, pure = ((None, False) if tool_calls
                                   else _extract_code(msg.content))

                assistant_msg = {"role": "assistant", "content": msg.content or ""}
                if tool_calls:
                    assistant_msg["tool_calls"] = [
                        {"id": tc.id, "type": "function",
                         "function": {"name": tc.function.name,
                                      "arguments": tc.function.arguments}}
                        for tc in tool_calls]
                self._messages.append(assistant_msg)

                # Show assistant text, but suppress a message that is only code.
                if msg.content and msg.content.strip() and not (code_call and pure):
                    self.text.emit(msg.content)

                if tool_calls:
                    produced = False
                    for tc in tool_calls:
                        if self._cancel.is_set():
                            break
                        try:
                            args = json.loads(tc.function.arguments or "{}")
                        except Exception:
                            args = {}
                        self._messages.append({
                            "role": "tool", "tool_call_id": tc.id,
                            "content": self._run(tc.function.name, args)})
                        produced = True
                    if not produced:
                        break
                elif code_call:
                    # No tool-call protocol — feed the result back as a plain user
                    # turn so any provider can continue.
                    result = self._run("run_python", {"code": code_call})
                    self._messages.append({
                        "role": "user", "content": f"[run_python output]\n{result}"})
                else:
                    break
            self.done.emit()
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}")

    def _run(self, name, args):
        """Marshal one tool call to the GUI thread and wait for its result text."""
        holder = {"event": threading.Event(), "content": ""}
        self.tool_call.emit({"id": "x", "name": name, "input": args, "holder": holder})
        holder["event"].wait()
        return holder["content"]


class Assistant(QObject):
    """Owns the conversation + kernel; drives a worker per user turn and executes
    tools (with a confirm gate) on the GUI thread."""

    assistant_text = Signal(str)
    tool_ran = Signal(str, str, object)   # code, output_text, [figure_paths]
    tool_declined = Signal(str)           # code
    failed = Signal(str)
    finished = Signal()

    def __init__(self, main_window):
        super().__init__(main_window)
        from .kernel import PythonKernel
        self._mw = main_window
        self._kernel = PythonKernel(main_window.doc)
        self.config = AssistantConfig(getattr(main_window, "settings", None))
        self._messages = []
        self._worker = None
        # What the MODEL CHECKS blocks have already said, so a later block reports
        # the delta. A new or reopened project is a different model, and its
        # findings are different findings, so the document clears it too.
        self._checks_memo = ChecksMemo()
        try:
            main_window.doc.loaded.connect(self._checks_memo.reset)
        except AttributeError:
            pass

    # --- lifecycle -------------------------------------------------------
    def is_busy(self):
        return self._worker is not None and self._worker.isRunning()

    def reset(self):
        self._messages = []
        self._kernel.reset()          # fresh kernel — variables cleared, re-seeds
        self._checks_memo.reset()     # nothing has been reported to this session

    def _system(self):
        # Full skill body only for Anthropic (prompt-cached, so cheap). Local /
        # other models get the compact prompt and introspect at runtime — the skill
        # body measures ~33k tokens today, and sending that every turn makes a
        # local model crawl.
        if self.config.wants_skill():
            skill = _load_skill_text()
            if skill:
                # The skill already carries the modeling rules in full, so only the
                # record schemas are appended (LAST, so they win on recency — no
                # .xlsx diving to learn the schema). No MODELING_BRIEF here: it would
                # just duplicate the skill.
                return (STUDIO_SYSTEM + _SKILL_HEADER + skill + _SKILL_TRAILER
                        + SCHEMA_BRIEF)
        # No skill loaded (local/other model): the modeling rules live only here.
        return STUDIO_SYSTEM + SCHEMA_BRIEF + MODELING_BRIEF

    def send(self, user_text, images=None):
        if self.is_busy():
            return
        if not self.config.is_ready():
            self.failed.emit("No API key set for "
                             f"{self.config.display_name()}. Open the assistant "
                             "Settings to add one (or switch to a local Ollama model).")
            return
        if images:
            # Multimodal user turn (OpenAI/LiteLLM format; translated to the
            # Anthropic image format for Claude). ``images`` are data: URLs.
            content = [{"type": "text", "text": user_text or "(see image)"}]
            content += [{"type": "image_url", "image_url": {"url": u}} for u in images]
            self._messages.append({"role": "user", "content": content})
        else:
            self._messages.append({"role": "user", "content": user_text})
        self._worker = _AgentWorker(self.config.completion_kwargs(), self._system(),
                                    self._messages, [RUN_PYTHON_TOOL],
                                    cache_system=self.config.supports_prompt_cache(),
                                    parent=self)
        self._worker.text.connect(self.assistant_text)
        self._worker.tool_call.connect(self._on_tool_call)   # queued -> GUI thread
        self._worker.failed.connect(self._on_failed)
        self._worker.done.connect(self._on_done)
        self._worker.start()

    def cancel(self):
        if self._worker is not None:
            self._worker.cancel()

    # --- tool execution (GUI thread) ------------------------------------
    def _on_tool_call(self, req):
        holder = req["holder"]
        if req["name"] != "run_python":
            holder["content"] = f"Unknown tool: {req['name']}"
            holder["event"].set()
            return
        code = (req["input"] or {}).get("code", "")

        if self.config.confirm_before_run() and not self._confirm_run(code):
            self.tool_declined.emit(code)
            holder["content"] = "The user declined to run this code."
            holder["event"].set()
            return

        stdout, outputs, error, checks = self._run_python(code)
        parts = []
        if stdout.strip():
            parts.append(stdout.rstrip())
        if outputs:
            names = ", ".join(os.path.basename(p) for p in outputs)
            parts.append(f"[saved {len(outputs)} file(s) the user can open: {names}]")
        if error:
            parts.append("ERROR:\n" + error)
        if not parts:
            parts.append("(no output)")
        if checks:
            parts.append(checks)        # LAST: the final thing the model reads
        result_text = "\n".join(parts)

        holder["content"] = result_text
        holder["event"].set()
        self.tool_ran.emit(code, result_text, outputs)

    def _run_python(self, code):
        """Run one snippet and return ``(stdout, outputs, error, checks)``.

        ``checks`` is the MODEL CHECKS block, and it is produced ONLY when the
        snippet actually changed the model — the same per-key signature comparison
        that decides whether the run becomes an undo step. A read-only query costs
        nothing extra, and an edit pays one preflight pass. The session's
        :class:`ChecksMemo` carries what earlier blocks already said, so a repeat
        of a standing finding costs one line instead of its paragraph.
        """
        doc = self._mw.doc
        if doc.slope_data is None:
            # The assistant builds into a live document; if nothing is open, start
            # an empty project so `slope_data` is a real dict from the first snippet
            # (otherwise the model wastes turns on doc.new() + re-fetch).
            doc.new()
        before = _source_key_sigs(doc.slope_data)
        mesh_before = self._mw.mesh_signature(doc.slope_data)
        doc.begin_edit("Assistant edit")    # snapshot for undo / rollback
        stdout, outputs, error = self._kernel.run(code)
        edited = geom_changed = False
        if error:
            # Transactional: a snippet that raised leaves no partial edit, so
            # trial-and-error retries can't compound (e.g. a +5 applied 5 times).
            doc.rollback_edit()
        else:
            after = _source_key_sigs(doc.slope_data)
            if before is None or after is None:
                changed = None              # couldn't compare -> assume an edit
            else:
                changed = [k for k in _SOURCE_KEYS if before.get(k) != after.get(k)]
            if changed is None or changed:
                doc.commit_edit()           # real edit -> re-render + mark dirty
                edited = True
                if changed:                 # tag the undo step by what it touched
                    doc.relabel_last_edit(_assistant_edit_label(changed))
                geom_changed = self._mw.mesh_signature(doc.slope_data) != mesh_before
            else:
                doc.cancel_edit()           # read-only run -> no dirty, no undo step
        try:
            if edited:
                # Inputs changed: any cached solution is now stale (and the mesh
                # too, if the geometry changed).
                self._mw.invalidate_results(clear_mesh=geom_changed)
            self._mw.refresh_inputs_view()
        except Exception:
            pass
        checks = (model_checks_text(doc.slope_data, self._checks_memo)
                  if edited else "")
        return stdout, outputs, error, checks

    def output_dir(self):
        """Folder where the assistant writes generated files (plots, CSVs, …)."""
        return self._kernel.outdir

    def _confirm_run(self, code):
        from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QLabel,
                                       QPlainTextEdit, QVBoxLayout)
        dlg = QDialog(self._mw)
        dlg.setWindowTitle("Run code?")
        lay = QVBoxLayout(dlg)
        lay.addWidget(QLabel("The assistant wants to run this Python:"))
        view = QPlainTextEdit(code)
        view.setReadOnly(True)
        view.setMinimumSize(560, 280)
        lay.addWidget(view)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Run")
        bb.button(QDialogButtonBox.Cancel).setText("Don't run")
        bb.accepted.connect(dlg.accept)
        bb.rejected.connect(dlg.reject)
        lay.addWidget(bb)
        return dlg.exec() == QDialog.Accepted

    # --- worker signals --------------------------------------------------
    def _on_failed(self, message):
        self.failed.emit(self._friendly(message))
        self._cleanup()

    def _friendly(self, message):
        """Turn a raw provider/litellm error into a short, actionable line."""
        import re
        m = message.lower()
        prov = self.config.provider()
        if "credit balance is too low" in m or "plans & billing" in m:
            return ("Your Anthropic API account is out of credits. Add credits at "
                    "console.anthropic.com → Plans & Billing, then try again.")
        if ("invalid x-api-key" in m or "authentication" in m
                or "invalid_api_key" in m or "no api key" in m or " 401" in m):
            return "Invalid or missing API key — open Settings and check the key."
        if ("connection error" in m or "connection refused" in m
                or "max retries" in m or "failed to establish" in m
                or "all connection attempts failed" in m or "timed out" in m):
            if prov == "ollama":
                return (f"Can't reach Ollama at {self.config.ollama_base()} — is it "
                        "running? Start the Ollama app, or run 'ollama serve'.")
            return "Network error reaching the model — check your connection."
        if "not found" in m and prov == "ollama":
            model = self.config.model()
            return (f"Model '{model}' isn't available in Ollama. Pull it with "
                    f"'ollama pull {model}', or fix the name in Settings.")
        if "rate limit" in m or " 429" in m or "overloaded" in m:
            return "Rate limited or overloaded — wait a moment and try again."
        # Fallback: pull the human message out of the provider JSON, else trim.
        inner = re.search(r'"message"\s*:\s*"([^"]+)"', message)
        if inner:
            return inner.group(1)
        return message if len(message) < 200 else message[:200] + "…"

    def _on_done(self):
        self.finished.emit()
        self._cleanup()

    def _cleanup(self):
        if self._worker is not None:
            self._worker.deleteLater()
            self._worker = None
