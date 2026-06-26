"""Assistant — the agent loop behind the chat dock (Phase A spike).

A manual Claude tool-use loop (so tool execution can be gated by a confirm
dialog) running on a ``QThread``: the blocking ``messages.create`` calls happen
off the GUI thread, while each ``run_python`` tool call is marshalled back to the
GUI thread (modal confirm + document mutation + re-render) and the worker waits
for the result. Conversation history persists across turns on the controller.

Claude-only for now via the official ``anthropic`` SDK (model ``claude-opus-4-8``,
adaptive thinking). The system prompt is the Studio framing plus the existing
``/xslope`` skill knowledge. See plan_gui.md §14.
"""

from __future__ import annotations

import json
import os
import threading

from PySide6.QtCore import QObject, QThread, Signal

MODEL = "claude-opus-4-8"
MAX_TOKENS = 16000

RUN_PYTHON_TOOL = {
    "name": "run_python",
    "description": (
        "Execute Python in the live XSlope Studio session and return its stdout. "
        "A persistent namespace is preloaded with: `xslope` (the engine), `np`, "
        "`pd`, `plt` (matplotlib.pyplot), and the open project — `doc` "
        "(ProjectDocument), `slope_data` (the dict you edit), and `results`. "
        "Variables persist across calls like a notebook. Build or edit the input "
        "by mutating `slope_data` in place (it re-renders on the canvas; the user "
        "saves via Save As) rather than writing an .xlsx. Any `plt` figures you "
        "create are shown to the user automatically, so just create them — you "
        "don't receive the images back. Print values you need to read."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "code": {"type": "string", "description": "Python source to execute."},
        },
        "required": ["code"],
        "additionalProperties": False,
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

The reference below documents the `slope_data` schema and the engine API \
(`load_slope_data`, `generate_slices`, the LEM solvers, `circular_search`, the \
seep/FEM builders, the `plot_*` helpers). Use it as ground truth for keys and \
function signatures; it was written for a file-based workflow, so favor the \
in-memory `slope_data` document over its .xlsx-writing patterns.
"""


def _load_skill_text():
    """The /xslope skill body (schema + API knowledge), best-effort. Repo-bound
    for now; packaging it travels with §14.5."""
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    path = os.path.join(here, "docs", "usage", "claude", "xslope.md")
    try:
        with open(path, encoding="utf-8") as f:
            return f.read()
    except Exception:
        return ""


class _AgentWorker(QThread):
    """Runs the Claude tool-use loop off the GUI thread."""

    text = Signal(str)              # an assistant text block
    tool_call = Signal(object)      # {id, name, input, holder} — handled on GUI thread
    failed = Signal(str)
    done = Signal()

    def __init__(self, client, system, messages, tools, parent=None):
        super().__init__(parent)
        self._client = client
        self._system = system
        self._messages = messages   # shared list, mutated in place across turns
        self._tools = tools
        self._cancel = threading.Event()

    def cancel(self):
        self._cancel.set()

    def run(self):
        try:
            while not self._cancel.is_set():
                resp = self._client.messages.create(
                    model=MODEL, max_tokens=MAX_TOKENS,
                    system=self._system, messages=self._messages,
                    tools=self._tools, thinking={"type": "adaptive"})
                # Append the full content (incl. thinking blocks) for correct replay.
                self._messages.append({"role": "assistant", "content": resp.content})
                for block in resp.content:
                    if block.type == "text" and block.text.strip():
                        self.text.emit(block.text)

                if resp.stop_reason != "tool_use":
                    break

                results = []
                for block in resp.content:
                    if block.type != "tool_use":
                        continue
                    if self._cancel.is_set():
                        break
                    holder = {"event": threading.Event(), "content": "", "is_error": False}
                    self.tool_call.emit({"id": block.id, "name": block.name,
                                         "input": block.input, "holder": holder})
                    holder["event"].wait()
                    results.append({"type": "tool_result", "tool_use_id": block.id,
                                    "content": holder["content"],
                                    "is_error": holder["is_error"]})
                if not results:
                    break
                self._messages.append({"role": "user", "content": results})
            self.done.emit()
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}")


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
        self._messages = []
        self._worker = None
        self._client = None
        self.confirm = True               # autonomy: confirm before running code

    # --- lifecycle -------------------------------------------------------
    def is_busy(self):
        return self._worker is not None and self._worker.isRunning()

    def reset(self):
        """Clear the conversation (e.g. on opening a different project)."""
        self._messages = []

    def _ensure_client(self):
        if self._client is None:
            import anthropic              # imported lazily — optional dependency
            self._client = anthropic.Anthropic()
        return self._client

    def _system(self):
        skill = _load_skill_text()
        text = STUDIO_SYSTEM + ("\n\n---\n\n" + skill if skill else "")
        return [{"type": "text", "text": text, "cache_control": {"type": "ephemeral"}}]

    def send(self, user_text):
        if self.is_busy():
            return
        try:
            client = self._ensure_client()
        except ImportError:
            self.failed.emit("The 'anthropic' package is not installed. "
                             "Install it with: pip install \"xslope[ai]\"")
            return
        except Exception as exc:
            self.failed.emit(f"Could not start the assistant: {exc}")
            return
        self._messages.append({"role": "user", "content": user_text})
        self._worker = _AgentWorker(client, self._system(), self._messages,
                                    [RUN_PYTHON_TOOL], parent=self)
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
            holder["is_error"] = True
            holder["event"].set()
            return
        code = (req["input"] or {}).get("code", "")

        if self.confirm and not self._confirm_run(code):
            self.tool_declined.emit(code)
            holder["content"] = "The user declined to run this code."
            holder["is_error"] = True
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
        holder["is_error"] = bool(error)
        holder["event"].set()
        self.tool_ran.emit(code, result_text, figures)

    def _run_python(self, code):
        doc = self._mw.doc
        # Snapshot for undo, then re-render + refresh the inputs tree after.
        doc.begin_edit()
        stdout, figures, error = self._kernel.run(code)
        doc.commit_edit()
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
