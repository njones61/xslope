"""ChatDock — the AI assistant panel.

A narrow right-side dock: a transcript (the assistant's markdown rendered — tables,
headings, fenced code, and math reduced to plain text since the dialect renders
none — with inline figures and "ran code" blocks), a multi-line
input (Ctrl/Cmd+Enter sends), Send/Stop, and a Settings button (provider/model,
key, confirm-before-running). It wires to
``studio.ai.assistant.Assistant`` and renders its signals. The header is kept
minimal (a wrapping model label + Settings) so the dock can be made narrow.
"""

from __future__ import annotations

import html
import os
import re
import subprocess
import sys

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, Qt, QUrl, Signal
from PySide6.QtGui import (
    QDesktopServices, QImage, QTextDocument, QTextOption,
)
from PySide6.QtWidgets import (
    QHBoxLayout, QLabel, QPlainTextEdit, QPushButton, QTextBrowser,
    QVBoxLayout, QWidget,
)

_IMAGE_EXTS = (".png", ".jpg", ".jpeg", ".gif", ".bmp", ".webp")

#: Fenced blocks and inline code spans, kept whole so their contents are never
#: HTML-escaped — inside a code span ``&lt;`` is four characters, not a ``<``.
#: The unterminated fence is matched too, so a reply cut off mid-block still
#: renders as code rather than as a wall of entities.
_CODE_SPANS = re.compile(r"(```.*?```|```[\s\S]*$|``.*?``|`[^`\n]*`)", re.S)

#: Greek letters a slope-stability answer actually uses, as the letters
#: themselves. ``\phi`` renders as the same phi the rest of the program shows
#: (φ, not the variant glyph) so a reply and the material table agree.
_LATEX_LETTERS = {
    "alpha": "α", "beta": "β", "gamma": "γ", "Gamma": "Γ", "delta": "δ",
    "Delta": "Δ", "epsilon": "ε", "varepsilon": "ε", "eta": "η", "theta": "θ",
    "Theta": "Θ", "kappa": "κ", "lambda": "λ", "Lambda": "Λ", "mu": "μ",
    "nu": "ν", "xi": "ξ", "pi": "π", "rho": "ρ", "sigma": "σ", "Sigma": "Σ",
    "tau": "τ", "phi": "φ", "varphi": "φ", "Phi": "Φ", "chi": "χ", "psi": "ψ",
    "Psi": "Ψ", "omega": "ω", "Omega": "Ω",
}

#: Operators and relations, as the Unicode character.
_LATEX_SYMBOLS = {
    "times": "×", "cdot": "·", "div": "÷", "pm": "±", "mp": "∓",
    "le": "≤", "leq": "≤", "ge": "≥", "geq": "≥", "ne": "≠", "neq": "≠",
    "approx": "≈", "sim": "~", "equiv": "≡", "propto": "∝", "infty": "∞",
    "partial": "∂", "nabla": "∇", "int": "∫", "sum": "Σ", "prod": "Π",
    "rightarrow": "→", "to": "→", "leftarrow": "←", "Rightarrow": "⇒",
    "circ": "°", "degree": "°", "deg": "°", "ldots": "…", "dots": "…",
    "quad": " ", "qquad": " ", "%": "%", "$": "$", "&": "&", "#": "#",
    "{": "{", "}": "}", "_": "_", "^": "^",
}

#: Functions, kept as their plain names with the space LaTeX implies.
_LATEX_FUNCS = ("arctan", "arcsin", "arccos", "tan", "sin", "cos", "tanh",
                "sinh", "cosh", "ln", "log", "exp", "sec", "csc", "cot",
                "max", "min")

#: Macros handled structurally (they take arguments), so the letter/symbol pass
#: must leave them alone.
_LATEX_STRUCTURAL = {"frac", "dfrac", "tfrac", "sqrt", "text", "textrm",
                     "mathrm", "mathit", "mathbf", "operatorname",
                     "left", "right", "big", "Big", "bigg", "Bigg"}

_SUP_OK = "0123456789+-()ni"
_SUB_OK = "0123456789+-()aehijklmnoprstuvx"
_SUPERSCRIPTS = str.maketrans(_SUP_OK, "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁽⁾ⁿⁱ")
_SUBSCRIPTS = str.maketrans(_SUB_OK, "₀₁₂₃₄₅₆₇₈₉₊₋₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓ")
_SUP_OK, _SUB_OK = set(_SUP_OK), set(_SUB_OK)

#: ``$$…$$`` and ``\[…\]`` are always math; so is ``\(…\)``.
_MATH_DISPLAY = re.compile(r"\$\$(.+?)\$\$", re.S)
_MATH_BRACKET = re.compile(r"\\\[(.+?)\\\]", re.S)
_MATH_PAREN = re.compile(r"\\\((.+?)\\\)", re.S)
#: ``$…$`` is the ambiguous one — a reply that quotes a price ("$5.00 input,
#: $0.50 cached") has the same shape — so it is math only where it carries a
#: LaTeX signal (a macro, a brace, a script) or reads as a stated equation whose
#: first character is not a digit.
_MATH_INLINE = re.compile(r"(?<![\\$])\$(?!\s)([^$\n]+?)(?<!\s)\$(?!\$)")
_MATH_SIGNAL = re.compile(r"[\\{^_]")


def _latex_body_to_text(body):
    """One math expression, rendered as the plain text a reader can read.

    Qt's markdown is the GitHub dialect and has no math, so an equation arriving
    as LaTeX reaches the transcript as its own source. Every command this
    converts is one seen in a real reply; anything else loses its backslash and
    keeps its name, which is still readable, where the raw macro is not.
    """
    text = body

    # Functions first, and letters with them: the fraction pass below decides
    # its spacing from whether the converted numerator reads as one token
    # ("1/2") or as words ("tan φ / tan β").
    for name in _LATEX_FUNCS:
        text = re.sub(r"\\" + name + r"(?![A-Za-z])", name + " ", text)

    def _one(m):
        name = m.group(1)
        if name in _LATEX_STRUCTURAL:
            return m.group(0)
        return _LATEX_LETTERS.get(name, _LATEX_SYMBOLS.get(name, name))

    text = re.sub(r"\\([A-Za-z]+|[%${}&#_^])", _one, text)
    text = re.sub(r"\\[ ,;:!]", " ", text)               # spacing macros

    # A script is lifted only when EVERY character of it has a Unicode form —
    # "sigma_{n}" becomes σₙ, while "x^{a+b}" keeps its caret rather than
    # becoming half-raised nonsense. It runs before the argument-taking macros
    # so their braces are the only ones left for those to match.
    def _script(table, allowed):
        def sub(m):
            group = m.group(1) if m.group(1) is not None else m.group(2)
            if group and set(group) <= allowed:
                return group.translate(table)
            return m.group(0)
        return sub

    text = re.sub(r"\^\{([^{}]*)\}|\^(\w)", _script(_SUPERSCRIPTS, _SUP_OK), text)
    text = re.sub(r"_\{([^{}]*)\}|_(\w)", _script(_SUBSCRIPTS, _SUB_OK), text)

    def _fraction(m):
        # Parenthesized where the operand is a sum, since a fraction bar groups
        # what a slash does not: (c' + σ' tan φ') / τ, never c' + σ' tan φ' / τ.
        def side(s):
            s = s.strip()
            return (f"({s})" if re.search(r"(?<=\S)\s*[+\-±]\s*(?=\S)", s)
                    else s)
        top, bottom = side(m.group(1)), side(m.group(2))
        joint = " / " if " " in top or " " in bottom else "/"
        return f"{top}{joint}{bottom}"

    # Argument-taking macros, innermost group first so nesting resolves.
    for _ in range(6):
        before = text
        text = re.sub(r"\\(?:text|textrm|mathrm|mathit|mathbf|operatorname)"
                      r"\s*\{([^{}]*)\}", r"\1", text)
        text = re.sub(r"\\sqrt\s*\{([^{}]*)\}", r"sqrt(\1)", text)
        text = re.sub(r"\\(?:d|t)?frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}",
                      _fraction, text)
        if text == before:
            break
    text = re.sub(r"\\(?:left|right|[Bb]igg?)(?![A-Za-z])", "", text)

    text = text.replace("\\\\", " ").replace("{", "").replace("}", "")
    text = re.sub(r"[ \t]{2,}", " ", text)
    return text.strip()


def strip_latex(text):
    """A reply's LaTeX math as plain text, leaving code and prose alone.

    The assistant is told to write math as plain text; this is what happens when
    it writes LaTeX anyway. Delimiters (``$$…$$``, ``$…$``, ``\\[…\\]``,
    ``\\(…\\)``) are removed and the expression inside converted, so the reader
    gets ``tan φ / tan β = 0.066`` rather than the macro that produced it. Code
    spans are untouched: a snippet that contains a dollar sign or a backslash
    means them literally.
    """
    def inline(m):
        body = m.group(1)
        if _MATH_SIGNAL.search(body) or ("=" in body and not body[:1].isdigit()):
            return _latex_body_to_text(body)
        return m.group(0)

    def convert(part):
        for pattern in (_MATH_DISPLAY, _MATH_BRACKET, _MATH_PAREN):
            part = pattern.sub(lambda m: _latex_body_to_text(m.group(1)), part)
        return _MATH_INLINE.sub(inline, part)

    return "".join(part if part.startswith("`") else convert(part)
                   for part in _CODE_SPANS.split(text or ""))


def _escape_outside_code(text):
    """HTML-escape prose while leaving code spans alone.

    The assistant writes markdown, and markdown passes raw HTML through: an
    unescaped ``<not html>`` in a sentence is swallowed by the parser and the
    reader never sees it. Escaping first fixes that everywhere the escape is
    decoded again — which is everywhere except inside code, where markdown is
    literal by definition.
    """
    return "".join(part if part.startswith("`") else html.escape(part, quote=False)
                   for part in _CODE_SPANS.split(text or ""))


def markdown_to_html(text):
    """One markdown message as an HTML fragment for the transcript.

    Qt's own parser (GitHub dialect: tables, headings, lists, fenced code) does
    the work, and the fragment is the body of what it produces — the document
    wrapper and its font declaration are dropped so the appended block inherits
    the transcript's. A parse that fails falls back to escaped plain text, since
    an unrendered reply must still be a readable one.

    The dialect has no math, so any LaTeX the reply carries is converted to plain
    text first (:func:`strip_latex`) rather than shown as its own source.
    """
    try:
        doc = QTextDocument()
        doc.setMarkdown(_escape_outside_code(strip_latex(text)),
                        QTextDocument.MarkdownDialectGitHub)
        rendered = doc.toHtml()
        open_tag = rendered.find("<body")
        start = rendered.find(">", open_tag)
        end = rendered.rfind("</body>")
        if open_tag != -1 and start != -1 and end > start:
            return rendered[start + 1:end]
        return rendered
    except Exception:
        return html.escape(text or "").replace("\n", "<br>")


def qimage_to_data_url(img):
    """Encode a QImage as a base64 PNG ``data:`` URL (the OpenAI/LiteLLM image
    format Claude and GPT both accept). Down-scales to a 1568px long edge —
    Claude's recommended max — to keep the token cost sane."""
    import base64
    if img.width() > 1568 or img.height() > 1568:
        img = img.scaled(1568, 1568, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    buf = QBuffer()
    buf.open(QIODevice.WriteOnly)
    img.save(buf, "PNG")
    b64 = base64.b64encode(bytes(buf.data())).decode("ascii")
    return f"data:image/png;base64,{b64}"


class _ChatInput(QPlainTextEdit):
    """Input box where Enter sends, Shift+Enter inserts a newline, and a pasted
    or dropped image is captured as an attachment instead of garbled text."""

    submit = Signal()
    image_added = Signal(QImage)

    def keyPressEvent(self, event):
        if (event.key() in (Qt.Key_Return, Qt.Key_Enter)
                and not (event.modifiers() & Qt.ShiftModifier)):
            self.submit.emit()
            return
        super().keyPressEvent(event)

    def canInsertFromMimeData(self, source):
        if source.hasImage() or self._image_urls(source):
            return True
        return super().canInsertFromMimeData(source)

    def insertFromMimeData(self, source):
        captured = False
        if source.hasImage():
            img = QImage(source.imageData())
            if not img.isNull():
                self.image_added.emit(img)
                captured = True
        for path in self._image_urls(source):
            img = QImage(path)
            if not img.isNull():
                self.image_added.emit(img)
                captured = True
        if captured:
            return
        super().insertFromMimeData(source)

    @staticmethod
    def _image_urls(source):
        if not source.hasUrls():
            return []
        return [u.toLocalFile() for u in source.urls()
                if u.isLocalFile() and u.toLocalFile().lower().endswith(_IMAGE_EXTS)]


class ChatDock(QWidget):
    def __init__(self, assistant, parent=None):
        super().__init__(parent)
        self._assistant = assistant
        self._img_seq = 0
        self._pending = []        # list[QImage] attached to the next message
        self._paths = {}          # link-id -> output file path (open / reveal)

        self.transcript = QTextBrowser()
        # A long unbroken token (a URL, a path, a JSON blob) wraps rather than
        # widening the dock — what the per-block ``word-wrap`` style used to do,
        # set once here because the markdown blocks carry Qt's own styles.
        self.transcript.setWordWrapMode(QTextOption.WrapAtWordBoundaryOrAnywhere)
        # A fenced code block in a reply reads as the "Ran code" blocks do.
        self.transcript.document().setDefaultStyleSheet(
            "pre { background-color: #f4f4f4; }")
        self.transcript.setOpenExternalLinks(False)
        self.transcript.setOpenLinks(False)            # we handle clicks ourselves
        self.transcript.anchorClicked.connect(self._on_anchor)

        self.attach_label = QLabel()
        self.attach_label.setStyleSheet("color:#555;")
        self.clear_attach_btn = QPushButton("Remove")
        self.clear_attach_btn.clicked.connect(self._clear_attachments)
        self.attach_row = QWidget()
        ar = QHBoxLayout(self.attach_row)
        ar.setContentsMargins(0, 0, 0, 0)
        ar.addWidget(self.attach_label, 1)
        ar.addWidget(self.clear_attach_btn, 0)
        self.attach_row.hide()

        self.input = _ChatInput()
        self.input.setPlaceholderText("Ask the assistant…  (Enter to send, "
                                      "Shift+Enter for a new line; paste or drop "
                                      "an image to attach it)")
        self.input.setFixedHeight(72)

        self.send_btn = QPushButton("Send")
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.clear_btn = QPushButton("New chat")
        self.clear_btn.setToolTip("Start a fresh conversation (clears history and "
                                  "the assistant's Python variables; your project "
                                  "is unaffected).")
        self.files_btn = QPushButton("Files…")
        self.files_btn.setToolTip("Open the folder where the assistant saves "
                                  "generated plots and files.")
        self.settings_btn = QPushButton("Settings…")

        self.model_label = QLabel()
        self.model_label.setWordWrap(True)        # wraps so the dock can be narrow
        self.model_label.setStyleSheet("color:#444;")
        # Buttons on their own row; the model label sits just below them so it has
        # the full dock width instead of being squeezed beside the buttons.
        btn_row = QHBoxLayout()
        btn_row.addWidget(self.files_btn)
        btn_row.addWidget(self.clear_btn)
        btn_row.addWidget(self.settings_btn)
        btn_row.addStretch(1)
        top = QVBoxLayout()
        top.addLayout(btn_row)
        top.addWidget(self.model_label)

        self.status_label = QLabel()
        self.status_label.setStyleSheet("color:#666; font-style:italic;")

        row = QHBoxLayout()
        row.addWidget(self.send_btn)
        row.addWidget(self.stop_btn)

        # What the conversation has cost, in tokens: this turn, and the session
        # behind it. Tokens only — a price would need a per-model rate table, and
        # those change faster than a release, so a stale one quoting dollars would
        # be worse than no dollars at all. It sits under the input, where the user
        # decides whether to send the next one.
        self.usage_label = QLabel()
        self.usage_label.setObjectName("assistant_usage")
        self.usage_label.setWordWrap(True)
        self.usage_label.setStyleSheet("color:#777; font-size:11px;")
        self.usage_label.setToolTip(
            "Tokens read and written by the model. 'cached' is the part of the "
            "input the provider served from its prompt cache.")

        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.transcript, 1)
        layout.addWidget(self.status_label)
        layout.addWidget(self.attach_row)
        layout.addWidget(self.input)
        layout.addLayout(row)
        layout.addWidget(self.usage_label)
        self.setMinimumWidth(220)

        self.input.image_added.connect(self._add_attachment)
        self.input.submit.connect(self._send)         # Enter sends
        self.send_btn.clicked.connect(self._send)
        self.stop_btn.clicked.connect(self._assistant.cancel)
        self.clear_btn.clicked.connect(self._clear)
        self.files_btn.clicked.connect(self._open_files_folder)
        self.settings_btn.clicked.connect(self._open_settings)

        self._assistant.assistant_text.connect(self._on_assistant_text)
        self._assistant.usage_changed.connect(self._on_usage)
        self._assistant.tool_ran.connect(self._on_tool_ran)
        self._assistant.tool_declined.connect(self._on_tool_declined)
        self._assistant.failed.connect(self._on_failed)
        self._assistant.finished.connect(self._on_finished)
        self._refresh_model_label()

    # --- actions ---------------------------------------------------------
    def _refresh_model_label(self):
        self.model_label.setText(self._assistant.config.display_name())

    def _open_settings(self):
        from .ai.settings_dialog import AssistantSettingsDialog
        if AssistantSettingsDialog(self._assistant.config, self).exec():
            self._refresh_model_label()

    def _clear(self):
        """Start a fresh conversation: wipe the transcript and the model's message
        history. The Python kernel (variables/state) is left intact."""
        if self._assistant.is_busy():
            return
        self._assistant.reset()
        self.transcript.clear()
        self._img_seq = 0
        self._paths.clear()
        self._clear_attachments()

    def _open_files_folder(self):
        d = self._assistant.output_dir()
        os.makedirs(d, exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(d))

    # --- image attachments ----------------------------------------------
    def _add_attachment(self, img):
        if img.isNull():
            return
        if not self._assistant.config.capabilities().get("vision", True):
            self.transcript.append(
                '<div style="color:#9a6700;margin-top:8px;">The current model '
                'cannot read images. Every model the Settings dialog lists can — '
                'pick one from the list rather than a typed-in id.</div>')
            return
        self._pending.append(img)
        self._refresh_attachments()

    def _clear_attachments(self):
        self._pending = []
        self._refresh_attachments()

    def _refresh_attachments(self):
        n = len(self._pending)
        if n:
            self.attach_label.setText(f"📎 {n} image{'s' if n != 1 else ''} attached")
            self.attach_row.show()
        else:
            self.attach_row.hide()

    def _send(self):
        text = self.input.toPlainText().strip()
        images = list(self._pending)
        if (not text and not images) or self._assistant.is_busy():
            return
        self.input.clear()
        self._add_block("You", text or "(image)", "#1a5fb4")
        for img in images:
            self._append_qimage(img)
        self._clear_attachments()
        self.send_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.status_label.setText(
            "Working… (a local model may take a minute to load on first use)")
        self._assistant.send(text, images=[qimage_to_data_url(i) for i in images])

    # --- assistant signals ----------------------------------------------
    def _on_usage(self, usage):
        """Show the running token count. Blank until the first completion, so a
        fresh dock (and a new chat) carries no line at all rather than zeros."""
        from .ai.assistant import format_usage
        session = (usage or {}).get("session") or {}
        if not any(session.get(k) for k in ("input", "cached_input", "output")):
            self.usage_label.clear()
            return
        self.usage_label.setText(format_usage(usage["turn"], session))

    def _on_assistant_text(self, text):
        self._add_markdown_block("Assistant", text, "#2e7d32")

    def _on_tool_ran(self, code, output, outputs):
        pre = ("background:#f4f4f4;padding:6px;border-radius:4px;"
               "white-space:pre-wrap;word-wrap:break-word;")
        frag = (f'<div style="margin-top:8px;"><b>Ran code</b>'
                f'<pre style="{pre}">{html.escape(code)}</pre>')
        if output and output != "(no output)":
            frag += (f'<pre style="{pre.replace("#f4f4f4", "#eef3ee")}">'
                     f'{html.escape(output)}</pre>')
        frag += "</div>"
        self.transcript.append(frag)            # one complete block -> own paragraph
        for path in outputs or []:              # plots inline, other files as links
            if path.lower().endswith(_IMAGE_EXTS):
                self._append_image(path)
            else:
                self._append_file(path)

    def _on_tool_declined(self, code):
        self.transcript.append('<div style="color:#9a6700;margin-top:8px;">'
                               'Declined to run the code.</div>')

    def _on_failed(self, message):
        self.transcript.append(f'<div style="color:#b00020;margin-top:8px;">'
                               f'<b>Error:</b> {html.escape(message)}</div>')
        self._idle()

    def _on_finished(self):
        self._idle()

    # --- rendering helpers ----------------------------------------------
    def _idle(self):
        self.send_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.status_label.clear()

    def _add_block(self, who, text, color):
        """A block of VERBATIM text (what the user typed), escaped and line-broken
        as written."""
        body = html.escape(text).replace("\n", "<br>")
        self.transcript.append(
            f'<div style="margin-top:8px;">'
            f'<b style="color:{color};">{who}:</b> {body}</div>')

    def _add_markdown_block(self, who, text, color):
        """A block of MARKDOWN (what the assistant wrote), rendered: tables as
        tables, ``##`` as headings, fenced code as monospaced blocks.

        The speaker's label leads on its own line rather than inline, because a
        reply that opens on a heading or a table has no first line to sit in.
        """
        self.transcript.append(f'<div style="margin-top:8px;">'
                               f'<b style="color:{color};">{who}:</b></div>')
        body = markdown_to_html(text)
        if body.strip():
            self.transcript.append(body)

    # --- generated outputs (clickable to open / reveal) -----------------
    def _register(self, path):
        """Map a file path to a short link id used by the open/reveal anchors."""
        ident = str(len(self._paths) + 1)
        self._paths[ident] = path
        return ident

    def _append_image(self, path):
        """An assistant-generated image: shown inline, clickable to open, with an
        open / show-in-folder caption."""
        img = QImage(path)
        if img.isNull():
            self._append_file(path)         # not loadable as image -> link it
            return
        ident = self._register(path)
        url = self._image_resource(img)
        name = html.escape(os.path.basename(path))
        self.transcript.append(
            f'<a href="xopen:{ident}"><img src="{url}"></a>'
            f'<div style="font-size:11px;color:#666;">{name} — '
            f'<a href="xopen:{ident}">open</a> &middot; '
            f'<a href="xreveal:{ident}">show in folder</a></div>')

    def _append_file(self, path):
        """A non-image generated file (CSV, xlsx, …): shown as a link row."""
        ident = self._register(path)
        name = html.escape(os.path.basename(path))
        self.transcript.append(
            f'<div style="margin-top:6px;">📄 <a href="xopen:{ident}">{name}</a> '
            f'<span style="font-size:11px;color:#666;">'
            f'(<a href="xreveal:{ident}">show in folder</a>)</span></div>')

    def _append_qimage(self, img):
        """Inline a QImage with no backing file (e.g. a pasted user attachment)."""
        url = self._image_resource(img)
        if url:
            self.transcript.append(f'<img src="{url}">')

    def _image_resource(self, img):
        if img.isNull():
            return None
        cap = max(220, self.transcript.viewport().width() - 24)
        if img.width() > cap:
            img = img.scaledToWidth(cap, Qt.SmoothTransformation)
        self._img_seq += 1
        url = QUrl(f"xslope-fig://{self._img_seq}")
        self.transcript.document().addResource(QTextDocument.ImageResource, url, img)
        return url.toString()

    def _on_anchor(self, url):
        """Open a web link in the browser, or open/reveal a generated file."""
        s = url.toString()
        scheme, _, ident = s.partition(":")
        if scheme in ("http", "https"):
            # The assistant cites docs pages by URL; hand them to the browser.
            QDesktopServices.openUrl(url)
            return
        path = self._paths.get(ident)
        if not path or not os.path.exists(path):
            return
        if scheme == "xreveal":
            self._reveal(path)
        else:
            QDesktopServices.openUrl(QUrl.fromLocalFile(path))

    @staticmethod
    def _reveal(path):
        """Reveal a file in the OS file manager (select it where supported)."""
        try:
            if sys.platform == "darwin":
                subprocess.run(["open", "-R", path], check=False)
            elif sys.platform.startswith("win"):
                subprocess.run(["explorer", f"/select,{path}"], check=False)
            else:
                QDesktopServices.openUrl(QUrl.fromLocalFile(os.path.dirname(path)))
        except Exception:
            pass
