"""AssistantSettingsDialog — pick provider/model and store the API key.

Cross-platform: the key goes to the OS keychain via ``keyring`` (no environment
variables). Ollama needs no key (local, free) but exposes a base-URL field.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFormLayout, QLabel,
    QLineEdit, QMessageBox, QVBoxLayout,
)

from .config import PROVIDERS


class AssistantSettingsDialog(QDialog):
    def __init__(self, config, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Assistant settings")
        self._config = config
        self.setMinimumWidth(460)

        layout = QVBoxLayout(self)
        form = QFormLayout()
        # Let the field column grow to the dialog width (so the URL isn't clipped)
        # and give the rows breathing room.
        form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(10)

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

        self.base_edit = QLineEdit()
        self.base_edit.setPlaceholderText("API base URL")
        form.addRow("Base URL", self.base_edit)

        self.confirm = QCheckBox("Confirm before running code")
        self.confirm.setChecked(config.confirm_before_run())
        self.confirm.setToolTip("Show the code and ask before each run. Uncheck to "
                                "let the assistant run code automatically.")
        form.addRow("", self.confirm)

        layout.addLayout(form)
        self.note = QLabel()
        self.note.setWordWrap(True)
        layout.addWidget(self.note)
        # Capability warning (shown only when the model can't / may not run code).
        self.caps_warn = QLabel()
        self.caps_warn.setWordWrap(True)
        self.caps_warn.setStyleSheet("color:#9a6700;")
        layout.addWidget(self.caps_warn)
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
        self.model.currentTextChanged.connect(self._refresh_caps_note)  # vision is per-model
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
        has_base = spec.get("base") is not None
        self.base_edit.setEnabled(has_base)
        self.base_edit.setText(self._config.base_url(prov) or "")
        self._refresh_caps_note()

    def _refresh_caps_note(self):
        """Set the cost / tool / vision note for the current provider+model. Vision
        can depend on the chosen model (e.g. Z.ai's GLM-V models), so this also runs
        when the model changes, not just the provider."""
        spec = self._spec()
        needs_key = spec["needs_key"]
        cost = ("Local model — no key needed; runs on your machine, free."
                if not needs_key else
                "Hosted model — billed per token to your API account.")
        tools = spec.get("tools")
        tool_note = ("Tool use (run code): supported." if tools is True else
                     "Tool use (run code): not supported — chat only." if tools is False
                     else "Tool use (run code): depends on the local model you pick.")
        # Vision can be per-model (GLM-V etc.). The selection isn't saved yet, so
        # resolve it from the spec + the currently shown model name.
        vis = spec.get("vision")
        pat = spec.get("vision_match")
        if pat is not None:
            import re
            vis = bool(re.search(pat, self.model.currentText().strip().lower()))
        vision_note = ("  Vision (images): supported." if vis is True else
                       "  Vision (images): not on this model." if vis is False else "")
        cache = " Prompt caching reduces repeat-turn cost." if spec.get("prompt_cache") else ""
        self.note.setText(f"{cost}\n{tool_note}{vision_note}{cache}")

        # The relocated capability warning (was a caption in the dock).
        if tools is False:
            warn = "This model has no tool support — it can chat but can't run code."
        elif tools is None:
            warn = ("Running code needs a tool-calling model; some local models "
                    "don't support it, so run-code may not work.")
        else:
            warn = ""
        self.caps_warn.setText(warn)
        self.caps_warn.setVisible(bool(warn))

    def _save(self):
        prov = self.provider.currentData()
        if self._spec()["needs_key"]:
            key = self.api_key.text().strip()
            if key and any(ch.isspace() for ch in key):
                QMessageBox.warning(self, "Invalid API key",
                                    "That doesn't look like an API key — it contains "
                                    "spaces or line breaks. Paste only the key "
                                    "(e.g. starts with 'sk-').")
                return                     # keep the dialog open
            self._config.set_api_key(prov, key)
        model = self.model.currentText().strip()
        base = self.base_edit.text().strip() if self._spec().get("base") else None
        self._config.set_selection(prov, model, base)
        self._config.set_confirm_before_run(self.confirm.isChecked())
        self.accept()
