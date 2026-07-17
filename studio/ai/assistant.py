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
            "Execute Python in the live XSlope Studio session and return its "
            "stdout. A persistent namespace is preloaded with: `xslope` (the "
            "engine), `np`, `pd`, `plt` (matplotlib.pyplot), and the open project "
            "— `doc` (ProjectDocument), `slope_data` (the dict you edit), and "
            "`results`. Variables persist across calls like a notebook. Build or "
            "edit the input by mutating `slope_data` in place (it re-renders on "
            "the canvas; the user saves via Save As) rather than writing an .xlsx. "
            "Any `plt` figures you create are shown to the user automatically — "
            "you don't receive the images back. Print values you need to read."
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
You are the AI assistant embedded in **XSlope Studio**, a desktop GUI for the \
`xslope` slope-stability engine. You help the user build and edit slope models, \
run analyses (LEM / seepage / FEM), and explore results — by writing and running \
small Python snippets with the `run_python` tool.

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
# runtime (sending ~13k tokens of skill every turn makes a local model crawl).
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

CRITICAL — you are inside XSlope Studio, NOT the standalone file-based skill. The \
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

- materials[i]: name, gamma, option ('mc' Mohr-Coulomb | 'cp'), c, phi, cp, r_elev,
  d, psi, u ('none'|'piezo'|'seep'), sigma_gamma/c/phi/cp/d/psi (stddevs, 0 if
  unused), k1, k2, alpha, kr0, h0 (seepage), E, nu (FEM). Minimal MC example:
  {'name':'Clay','gamma':130.0,'option':'mc','c':400.0,'phi':0.0,'u':'none',
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
  toe point, and is a different circle.
- non_circ[i]: {'X':-10.0,'Y':0.0,'Movement':'Free'}
- piezo_line / piezo_line2: list of (x, y) tuples.
- dloads / dloads2: list of blocks; each block is a list of {'X','Y','Normal'} pts.
- reinforcement_lines[i]: {'x1','y1','x2','y2','t_max','t_res','lp1','lp2','area','E'}
  # EDIT THIS one; reinforce_lines (capitalized X/Y/T/Tres) is derived from it.
- pile_lines[i]: {'x1','y1','x2','y2','D_pile','S','E','I','area','M_cap','V_cap',
  'theta_p','fixity','label','H'}
- scalars: gamma_water, max_depth (hard-base elevation), k_seismic, tcrack_depth,
  tcrack_water, circular (bool).

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
  toe_elevation (that is only tangent to the toe level, a different circle).
- Extent: extend the flat ground far enough on BOTH sides that every trial circle
  daylights on the ground INSIDE the model, never at a vertical edge (>= ~2x slope
  height beyond toe and crest, more for deep base circles). Do NOT copy the source
  diagram's width — it is usually cropped. The hard base spans the full width.
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
    "gamma_water", "tcrack_depth", "tcrack_water", "k_seismic", "max_depth",
    "materials", "profile_lines", "polygons", "circles", "non_circ", "piezo_line",
    "piezo_line2", "dloads", "dloads2", "reinforcement_lines", "pile_lines",
    "seepage_bc", "seepage_bc2",
)

# Short tags for the undo-history dropdown, keyed by source input. Several keys map
# to the same tag (e.g. both piezo lines -> "piezo") so a multi-key edit reads cleanly.
_KEY_LABELS = {
    "gamma_water": "global", "tcrack_depth": "global", "tcrack_water": "global",
    "k_seismic": "global", "max_depth": "global",
    "materials": "materials", "profile_lines": "profile", "polygons": "polygons",
    "circles": "circles", "non_circ": "non-circular", "piezo_line": "piezo",
    "piezo_line2": "piezo", "dloads": "dloads", "dloads2": "dloads",
    "reinforcement_lines": "reinforcement", "pile_lines": "piles",
    "seepage_bc": "seep BC", "seepage_bc2": "seep BC",
}


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

    # --- lifecycle -------------------------------------------------------
    def is_busy(self):
        return self._worker is not None and self._worker.isRunning()

    def reset(self):
        self._messages = []
        self._kernel.reset()          # fresh kernel — variables cleared, re-seeds

    def _system(self):
        # Full skill body only for Anthropic (prompt-cached, so cheap). Local /
        # other models get the compact prompt and introspect at runtime — sending
        # ~13k skill tokens every turn makes a local model crawl.
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

        stdout, outputs, error = self._run_python(code)
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
        result_text = "\n".join(parts)

        holder["content"] = result_text
        holder["event"].set()
        self.tool_ran.emit(code, result_text, outputs)

    def _run_python(self, code):
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
        return stdout, outputs, error

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
