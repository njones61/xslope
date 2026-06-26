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
print any numbers you need to reason about.
- Prefer one focused snippet per step. Keep code short and readable. Print a \
brief result. Don't reformat or refactor the user's data beyond the request.
- When unsure of a key or signature, inspect at runtime instead of guessing: \
`print(sorted(slope_data))`, `print(slope_data['materials'][0])`, \
`import xslope; print([n for n in dir(xslope.solve)])`, `help(fn)`.
"""

# Appended only for Anthropic, where prompt caching makes the large skill body
# cheap. Local/other models get the compact prompt above and introspect at
# runtime (sending ~13k tokens of skill every turn makes a local model crawl).
_SKILL_HEADER = ("\n\n---\n\nReference — the `slope_data` schema and engine API "
                 "(ground truth for keys and signatures; it was written for a "
                 "file-based workflow, so favor the in-memory document over its "
                 ".xlsx-writing patterns):\n\n")

# A compact schema cheat-sheet for non-Anthropic models (the full skill is too
# large to send every turn to a local model). Keeps them from inventing keys.
SCHEMA_BRIEF = """\

`slope_data` top-level keys you edit in place (real names — don't invent others):
- materials: list of dicts (name, gamma, c, phi, option, cp, k1, k2, E, nu, …).
- profile_lines: list of {mat_id, coords:[(x,y),…]} — geometry (profile files).
- polygons: list of {polygon (shapely), mat_id} — geometry (polygon files).
- circles: list of {Xo, Yo, R, Option, Depth} — circular surfaces.
- non_circ: list of {X, Y, Movement} — non-circular surface points.
- piezo_line / piezo_line2: list of (x, y).
- dloads / dloads2: distributed-load blocks.
- reinforce_lines / reinforcement_lines: reinforcement; pile_lines: piles.
- seepage_bc / seepage_bc2, mesh, max_depth, gamma_water, k_seismic, tcrack_depth.

ALWAYS print one example record of a category before editing it — e.g.
`print(slope_data['reinforcement_lines'][0])` — so you use its real field names.
Then run analyses via the engine, e.g.
`from xslope.slice import generate_slices; from xslope.solve import solve_selected`.
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


def _load_skill_text():
    """The /xslope skill body (schema + API knowledge), best-effort. Repo-bound
    for now; packaging it travels with §14.5."""
    import os
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    path = os.path.join(here, "docs", "usage", "claude", "xslope.md")
    try:
        with open(path, encoding="utf-8") as f:
            return f.read()
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

                # Fallback for models that emit the call as plain-text JSON.
                text_call = (_parse_text_tool_call(msg.content, TOOL_NAMES)
                             if not tool_calls else None)

                assistant_msg = {"role": "assistant", "content": msg.content or ""}
                if tool_calls:
                    assistant_msg["tool_calls"] = [
                        {"id": tc.id, "type": "function",
                         "function": {"name": tc.function.name,
                                      "arguments": tc.function.arguments}}
                        for tc in tool_calls]
                self._messages.append(assistant_msg)

                # Show assistant text, but suppress a raw tool-call JSON.
                if msg.content and msg.content.strip() and not text_call:
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
                elif text_call:
                    name, args = text_call
                    # No tool-call protocol here — feed the result back as a
                    # plain user turn so any provider can continue.
                    result = self._run(name, args)
                    self._messages.append({
                        "role": "user", "content": f"[{name} output]\n{result}"})
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

    def _system(self):
        # Full skill body only for Anthropic (prompt-cached, so cheap). Local /
        # other models get the compact prompt and introspect at runtime — sending
        # ~13k skill tokens every turn makes a local model crawl.
        if self.config.supports_prompt_cache():
            skill = _load_skill_text()
            if skill:
                return STUDIO_SYSTEM + _SKILL_HEADER + skill
        return STUDIO_SYSTEM + SCHEMA_BRIEF

    def send(self, user_text):
        if self.is_busy():
            return
        if not self.config.is_ready():
            self.failed.emit("No API key set for "
                             f"{self.config.display_name()}. Open the assistant "
                             "Settings to add one (or switch to a local Ollama model).")
            return
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

        stdout, figures, error = self._run_python(code)
        parts = []
        if stdout.strip():
            parts.append(stdout.rstrip())
        if figures:
            parts.append(f"[{len(figures)} figure(s) shown to the user]")
        if error:
            parts.append("ERROR:\n" + error)
        if not parts:
            parts.append("(no output)")
        result_text = "\n".join(parts)

        holder["content"] = result_text
        holder["event"].set()
        self.tool_ran.emit(code, result_text, figures)

    def _run_python(self, code):
        doc = self._mw.doc
        doc.begin_edit()                    # snapshot for undo
        stdout, figures, error = self._kernel.run(code)
        doc.commit_edit()                   # re-render + mark dirty
        try:
            self._mw.refresh_inputs_view()
        except Exception:
            pass
        return stdout, figures, error

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
        self.failed.emit(message)
        self._cleanup()

    def _on_done(self):
        self.finished.emit()
        self._cleanup()

    def _cleanup(self):
        if self._worker is not None:
            self._worker.deleteLater()
            self._worker = None
