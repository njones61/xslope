"""ChatDock — the AI assistant panel.

A narrow right-side dock: a transcript (markdown-ish HTML with inline figures and
"ran code" blocks), a multi-line input (Ctrl/Cmd+Enter sends), Send/Stop, and a
Settings button (provider/model, key, confirm-before-running). It wires to
``studio.ai.assistant.Assistant`` and renders its signals. The header is kept
minimal (a wrapping model label + Settings) so the dock can be made narrow.
"""

from __future__ import annotations

import html

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import (
    QImage, QKeySequence, QShortcut, QTextCursor, QTextDocument,
)
from PySide6.QtWidgets import (
    QHBoxLayout, QLabel, QPlainTextEdit, QPushButton, QTextBrowser,
    QVBoxLayout, QWidget,
)


class ChatDock(QWidget):
    def __init__(self, assistant, parent=None):
        super().__init__(parent)
        self._assistant = assistant
        self._img_seq = 0

        self.transcript = QTextBrowser()
        self.transcript.setOpenExternalLinks(False)

        self.input = QPlainTextEdit()
        self.input.setPlaceholderText("Ask the assistant…  (Ctrl/Cmd+Enter to send)")
        self.input.setFixedHeight(72)

        self.send_btn = QPushButton("Send")
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.settings_btn = QPushButton("Settings…")

        self.model_label = QLabel()
        self.model_label.setWordWrap(True)        # wraps so the dock can be narrow
        top = QHBoxLayout()
        top.addWidget(self.model_label, 1)
        top.addWidget(self.settings_btn, 0, Qt.AlignTop)

        row = QHBoxLayout()
        row.addWidget(self.send_btn)
        row.addWidget(self.stop_btn)

        layout = QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.transcript, 1)
        layout.addWidget(self.input)
        layout.addLayout(row)
        self.setMinimumWidth(220)

        self.send_btn.clicked.connect(self._send)
        self.stop_btn.clicked.connect(self._assistant.cancel)
        self.settings_btn.clicked.connect(self._open_settings)
        for seq in ("Ctrl+Return", "Ctrl+Enter", "Meta+Return", "Meta+Enter"):
            QShortcut(QKeySequence(seq), self.input, activated=self._send)

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

    def _send(self):
        text = self.input.toPlainText().strip()
        if not text or self._assistant.is_busy():
            return
        self.input.clear()
        self._add_block("You", text, "#1a5fb4")
        self.send_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self._assistant.send(text)

    # --- assistant signals ----------------------------------------------
    def _on_assistant_text(self, text):
        self._add_block("Assistant", text, "#2e7d32")

    def _on_tool_ran(self, code, output, figures):
        self._append_html(
            f'<div style="margin:6px 0;"><b>Ran code</b>'
            f'<pre style="background:#f4f4f4;padding:6px;border-radius:4px;'
            f'white-space:pre-wrap;">{html.escape(code)}</pre>')
        if output and output != "(no output)":
            self._append_html(
                f'<pre style="background:#eef3ee;padding:6px;border-radius:4px;'
                f'white-space:pre-wrap;">{html.escape(output)}</pre>')
        for path in figures or []:
            self._append_image(path)
        self._append_html("</div>")

    def _on_tool_declined(self, code):
        self._append_html('<div style="color:#9a6700;">Declined to run the code.</div>')

    def _on_failed(self, message):
        self._append_html(f'<div style="color:#b00020;"><b>Error:</b> '
                          f'{html.escape(message)}</div>')
        self._idle()

    def _on_finished(self):
        self._idle()

    # --- rendering helpers ----------------------------------------------
    def _idle(self):
        self.send_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _add_block(self, who, text, color):
        self._append_html(
            f'<div style="margin:6px 0;"><b style="color:{color};">{who}:</b> '
            f'{html.escape(text).replace(chr(10), "<br>")}</div>')

    def _append_image(self, path):
        img = QImage(path)
        if img.isNull():
            return
        cap = max(220, self.transcript.viewport().width() - 24)
        if img.width() > cap:
            img = img.scaledToWidth(cap, Qt.SmoothTransformation)
        self._img_seq += 1
        url = QUrl(f"xslope-fig://{self._img_seq}")
        self.transcript.document().addResource(QTextDocument.ImageResource, url, img)
        self._append_html(f'<img src="{url.toString()}"><br>')

    def _append_html(self, fragment):
        self.transcript.moveCursor(QTextCursor.End)
        self.transcript.insertHtml(fragment)
        self.transcript.moveCursor(QTextCursor.End)
