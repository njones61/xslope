"""ChatDock — the AI assistant panel.

A narrow right-side dock: a transcript (markdown-ish HTML with inline figures and
"ran code" blocks), a multi-line input (Ctrl/Cmd+Enter sends), Send/Stop, and a
Settings button (provider/model, key, confirm-before-running). It wires to
``studio.ai.assistant.Assistant`` and renders its signals. The header is kept
minimal (a wrapping model label + Settings) so the dock can be made narrow.
"""

from __future__ import annotations

import html
import os
import subprocess
import sys

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QImage, QTextDocument
from PySide6.QtWidgets import (
    QHBoxLayout, QLabel, QPlainTextEdit, QPushButton, QTextBrowser,
    QVBoxLayout, QWidget,
)

_IMAGE_EXTS = (".png", ".jpg", ".jpeg", ".gif", ".bmp", ".webp")


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

        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.transcript, 1)
        layout.addWidget(self.status_label)
        layout.addWidget(self.attach_row)
        layout.addWidget(self.input)
        layout.addLayout(row)
        self.setMinimumWidth(220)

        self.input.image_added.connect(self._add_attachment)
        self.input.submit.connect(self._send)         # Enter sends
        self.send_btn.clicked.connect(self._send)
        self.stop_btn.clicked.connect(self._assistant.cancel)
        self.clear_btn.clicked.connect(self._clear)
        self.files_btn.clicked.connect(self._open_files_folder)
        self.settings_btn.clicked.connect(self._open_settings)

        self._assistant.assistant_text.connect(self._on_assistant_text)
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
                'has no vision support — switch to a Claude/OpenAI (or vision-'
                'capable Ollama) model in Settings to send images.</div>')
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
    def _on_assistant_text(self, text):
        self._add_block("Assistant", text, "#2e7d32")

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
        body = html.escape(text).replace("\n", "<br>")
        # word-wrap so a long unbroken token (URL/JSON) still wraps in the box.
        self.transcript.append(
            f'<div style="margin-top:8px;word-wrap:break-word;">'
            f'<b style="color:{color};">{who}:</b> {body}</div>')

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
        """Open or reveal a generated file when its link is clicked."""
        s = url.toString()
        scheme, _, ident = s.partition(":")
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
