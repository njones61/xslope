"""AssistantSettingsDialog — pick provider/model and store the API key.

Cross-platform: the key goes to the OS keychain via ``keyring`` (no environment
variables). Ollama needs no key (local, free) but exposes a base-URL field.

**The model box is built, not baked.** Opening the dialog (or switching provider,
or pressing Refresh) asks the selected provider for its own model list in the
background, falling back to the last list it saw and then to the list this
version shipped with; a curated manifest marks the recommended model and any
superseded ones. The box is editable throughout, so any model id can be typed —
including one published after this build. See :mod:`studio.ai.models`.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFormLayout, QHBoxLayout,
    QLabel, QLineEdit, QMessageBox, QPushButton, QVBoxLayout,
)

from . import models as model_list
from .config import PROVIDERS


class AssistantSettingsDialog(QDialog):
    def __init__(self, config, parent=None, auto_refresh=True):
        super().__init__(parent)
        self.setWindowTitle("Assistant settings")
        self._config = config
        self._auto_refresh = auto_refresh
        self._runner = None
        self._source = "static"
        self._listed = 0
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

        # The model box is ALWAYS editable: enumeration can fail, and a model can
        # be a day old. Its item text carries the manifest's marking while the id
        # rides in the item data, so what gets saved is the id either way.
        self.model = QComboBox()
        self.model.setEditable(True)
        self.model.setInsertPolicy(QComboBox.NoInsert)
        self.model.setToolTip("Pick a model, or type any model id the provider "
                              "accepts.")
        self.refresh_btn = QPushButton("Refresh")
        self.refresh_btn.setObjectName("refresh_models")
        self.refresh_btn.setToolTip("Ask the provider for its current model list.")
        self.refresh_btn.clicked.connect(lambda: self._start_refresh(manual=True))
        model_row = QHBoxLayout()
        model_row.setContentsMargins(0, 0, 0, 0)
        model_row.addWidget(self.model, 1)
        model_row.addWidget(self.refresh_btn)
        form.addRow("Model", model_row)

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
        # Where the model list came from (live / last known / shipped) and why.
        self.models_note = QLabel()
        self.models_note.setWordWrap(True)
        self.models_note.setObjectName("models_note")
        self.models_note.setStyleSheet("color:#666;")
        layout.addWidget(self.models_note)
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

    # --- the model list --------------------------------------------------
    def _spec(self):
        return PROVIDERS[self.provider.currentData()]

    def _settings(self):
        return getattr(self._config, "_s", None)

    def current_model(self):
        """The model id the dialog would save — the selected row's id, or exactly
        what the user typed when it matches no row."""
        text = self.model.currentText().strip()
        idx = self.model.findText(text)
        if idx >= 0:
            data = self.model.itemData(idx)
            if data:
                return str(data)
        return text

    def _stored_model(self, provider):
        """This provider's remembered model, or — on a fresh install — the
        manifest's recommendation, falling back to the shipped list."""
        settings = self._settings()
        stored = settings.value(f"ai/model/{provider}", "") if settings else ""
        stored = str(stored or "").strip()
        return stored or model_list.default_model(provider, settings)

    def _populate(self, models, source, keep=""):
        """Fill the model box from ``models``, overlaid with the manifest."""
        manifest = model_list.load_manifest(self._settings())
        provider = self.provider.currentData()
        choices = model_list.model_choices(provider, models, manifest)
        self._source = source
        self._listed = len([m for m in (models or []) if str(m).strip()])
        self.model.blockSignals(True)
        self.model.clear()
        for choice in choices:
            self.model.addItem(choice.text, choice.id)
            if choice.note:
                self.model.setItemData(self.model.count() - 1, choice.note,
                                       Qt.ToolTipRole)
        wanted = (keep or self._stored_model(provider)).strip()
        idx = self.model.findData(wanted)
        if idx >= 0:
            self.model.setCurrentIndex(idx)
        else:
            self.model.setCurrentText(wanted)     # free text survives, always
        self.model.blockSignals(False)
        self._refresh_caps_note()

    def _sync(self):
        prov = self.provider.currentData()
        spec = self._spec()
        cached = model_list.load_cached(self._settings(), prov)
        models, source = model_list.resolve_models(prov, cached=cached)
        self._populate(models, source)

        needs_key = spec["needs_key"]
        self.api_key.setEnabled(needs_key)
        self.api_key.setText(self._config.api_key(prov) if needs_key else "")
        has_base = spec.get("base") is not None
        self.base_edit.setEnabled(has_base)
        self.base_edit.setText(self._config.base_url(prov) or "")
        self._set_models_note()
        self._refresh_caps_note()
        if self._auto_refresh:
            self._start_refresh(manual=False)

    def _set_models_note(self, error=""):
        label = self._spec()["label"]
        if self._source == "live":
            text = f"{self._listed} models listed by {label}."
        elif self._source == "cache":
            text = f"Last list seen from {label}."
        else:
            text = "Showing the models this version shipped with."
        if error:
            text += f" {error}"
        elif self._source != "live" and self._needs_key_first():
            text += " Enter the API key to list the current models."
        self.models_note.setText(text + " Any model id can be typed in.")

    def _needs_key_first(self):
        return self._spec()["needs_key"] and not self.api_key.text().strip()

    def _start_refresh(self, manual=False):
        """Enumerate the provider's models in the background, and — at most once
        a day — refresh the recommendations manifest with it."""
        if self._runner is not None and self._runner.isRunning():
            return
        prov = self.provider.currentData()
        settings = self._settings()
        want_manifest = manual or model_list.manifest_due(settings)
        want_models = not self._needs_key_first()
        if not (want_models or want_manifest):
            return
        from ..runners import AssistantModelsRunner
        self._runner = AssistantModelsRunner(
            prov, api_key=self.api_key.text().strip(),
            base_url=(self.base_edit.text().strip() or None),
            want_models=want_models, want_manifest=want_manifest, parent=self)
        self._runner.listed.connect(self._on_models_listed)
        self._runner.finished.connect(self._on_refresh_finished)
        if manual:
            self.refresh_btn.setEnabled(False)
        self._runner.start()

    def _on_models_listed(self, payload):
        """Apply one enumeration result. Runs on the GUI thread, so this is also
        where the cache and the manifest are written."""
        if not isinstance(payload, dict):
            return
        provider = payload.get("provider")
        if provider != self.provider.currentData():
            return                       # the user moved on; the answer is stale
        settings = self._settings()
        if payload.get("manifest") is not None or payload.get("manifest_checked"):
            model_list.store_manifest(settings, payload.get("manifest"))
        models = payload.get("models")
        keep = self.current_model()
        if models:
            model_list.save_cached(settings, provider, models)
            self._populate(models, "live", keep=keep)
            self._set_models_note()
        else:
            # Nothing live: keep whatever is showing, but re-apply the overlay in
            # case the manifest is what changed.
            cached = model_list.load_cached(settings, provider)
            shown, source = model_list.resolve_models(provider, cached=cached)
            self._populate(shown, source, keep=keep)
            self._set_models_note(error=payload.get("error", ""))

    def _on_refresh_finished(self):
        self.refresh_btn.setEnabled(True)
        if self._runner is not None:
            self._runner.deleteLater()
            self._runner = None

    # --- the capability note ---------------------------------------------
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
            vis = bool(re.search(pat, self.current_model().lower()))
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
        model = self.current_model()
        base = self.base_edit.text().strip() if self._spec().get("base") else None
        self._config.set_selection(prov, model, base)
        self._config.set_confirm_before_run(self.confirm.isChecked())
        self.accept()
