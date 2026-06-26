"""ChatDock — the AI assistant panel.

A narrow right-side dock: a transcript (markdown-ish HTML with inline figures and
"ran code" blocks), a multi-line input (Ctrl/Cmd+Enter sends), Send/Stop, and a
Settings button (provider/model, key, confirm-before-running). It wires to
``studio.ai.assistant.Assistant`` and renders its signals. The header is kept
minimal (a wrapping model label + Settings) so the dock can be made narrow.
"""

from __future__ import annotations

import html

from PySide6.QtCore import QBuffer, QByteArray, QIODevice, Qt, QUrl, Signal
from PySide6.QtGui import QImage, QTextDocument
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

        self.transcript = QTextBrowser()
        self.transcript.setOpenExternalLinks(False)

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
        self.settings_btn = QPushButton("Settings…")

        self.model_label = QLabel()
        self.model_label.setWordWrap(True)        # wraps so the dock can be narrow
        top = QHBoxLayout()
        top.addWidget(self.model_label, 1)
        top.addWidget(self.clear_btn, 0, Qt.AlignTop)
        top.addWidget(self.settings_btn, 0, Qt.AlignTop)

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
        self._clear_attachments()

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

    def _on_tool_ran(self, code, output, figures):
        pre = ("background:#f4f4f4;padding:6px;border-radius:4px;"
               "white-space:pre-wrap;word-wrap:break-word;")
        frag = (f'<div style="margin-top:8px;"><b>Ran code</b>'
                f'<pre style="{pre}">{html.escape(code)}</pre>')
        if output and output != "(no output)":
            frag += (f'<pre style="{pre.replace("#f4f4f4", "#eef3ee")}">'
                     f'{html.escape(output)}</pre>')
        frag += "</div>"
        self.transcript.append(frag)            # one complete block -> own paragraph
        for path in figures or []:
            self._append_image(path)

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

    def _append_image(self, path):
        self._append_qimage(QImage(path))

    def _append_qimage(self, img):
        if img.isNull():
            return
        cap = max(220, self.transcript.viewport().width() - 24)
        if img.width() > cap:
            img = img.scaledToWidth(cap, Qt.SmoothTransformation)
        self._img_seq += 1
        url = QUrl(f"xslope-fig://{self._img_seq}")
        self.transcript.document().addResource(QTextDocument.ImageResource, url, img)
        self.transcript.append(f'<img src="{url.toString()}">')
