"""AssistantSettingsDialog — pick provider/model and store the API key.

Cross-platform: the key goes to the OS keychain via ``keyring`` (no environment
variables). Ollama needs no key (local, free) but exposes a base-URL field.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QFormLayout, QLabel, QLineEdit,
    QVBoxLayout,
)

from .config import PROVIDERS


class AssistantSettingsDialog(QDialog):
    def __init__(self, config, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Assistant settings")
        self._config = config

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.provider = QComboBox()
        for key, spec in PROVIDERS.items():
            self.provider.addItem(spec["label"], key)
        idx = self.provider.findData(config.provider())
        if idx >= 0:
            self.provider.setCurrentIndex(idx)
        form.addRow("Provider", self.provider)

        self.model = QComboBox()
        form.addRow("Model", self.model)

        self.api_key = QLineEdit()
        self.api_key.setEchoMode(QLineEdit.Password)
        self.api_key.setPlaceholderText("sk-…")
        form.addRow("API key", self.api_key)

        self.ollama_base = QLineEdit(config.ollama_base())
        form.addRow("Ollama URL", self.ollama_base)

        layout.addLayout(form)
        self.note = QLabel()
        self.note.setWordWrap(True)
        layout.addWidget(self.note)
        if not config.keyring_available():
            warn = QLabel("⚠ No system keychain found — the key will be stored in "
                          "app settings (less secure). Install 'keyring' for secure "
                          "storage.")
            warn.setWordWrap(True)
            layout.addWidget(warn)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self._save)
        bb.rejected.connect(self.reject)
        layout.addWidget(bb)

        self.provider.currentIndexChanged.connect(self._sync)
        self._sync()

    def _spec(self):
        return PROVIDERS[self.provider.currentData()]

    def _sync(self):
        prov = self.provider.currentData()
        spec = self._spec()
        # Repopulate the model list for this provider; allow free text for Ollama.
        self.model.blockSignals(True)
        self.model.clear()
        self.model.addItems(spec["models"])
        self.model.setEditable(bool(spec.get("editable_model")))
        cur = self._config._s.value(f"ai/model/{prov}", spec["models"][0])
        i = self.model.findText(cur)
        if i >= 0:
            self.model.setCurrentIndex(i)
        elif spec.get("editable_model"):
            self.model.setCurrentText(cur)
        self.model.blockSignals(False)

        needs_key = spec["needs_key"]
        self.api_key.setEnabled(needs_key)
        self.api_key.setText(self._config.api_key(prov) if needs_key else "")
        self.ollama_base.setEnabled(bool(spec.get("needs_base")))
        self.note.setText(
            "Local model — no key needed; runs on your machine, free." if not needs_key
            else "Hosted model — billed per token to your API account.")

    def _save(self):
        prov = self.provider.currentData()
        model = self.model.currentText().strip()
        self._config.set_selection(prov, model, self.ollama_base.text().strip())
        if self._spec()["needs_key"]:
            self._config.set_api_key(prov, self.api_key.text().strip())
        self.accept()
