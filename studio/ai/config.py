"""AssistantConfig — provider/model selection and credential storage.

Holds the user's choice of provider (Claude / OpenAI / local Ollama) and model in
``QSettings``, and API keys in the OS keychain via ``keyring`` (falling back to
QSettings if keyring is unavailable). Produces the kwargs for a ``litellm``
completion call, so the agent loop stays provider-agnostic.
"""

from __future__ import annotations

from PySide6.QtCore import QSettings

KEYRING_SERVICE = "XSlope Studio"

# Provider catalog. `prefix` is the LiteLLM model-name prefix; `models` are
# suggestions (the Ollama list is editable since the user picks their own).
# `tools`/`vision` are capability defaults (True/False, or None = depends on the
# chosen model — used for the hosted presets; local models vary).
PROVIDERS = {
    "anthropic": {
        "label": "Claude (Anthropic)", "prefix": "anthropic/", "needs_key": True,
        "models": ["claude-opus-4-8", "claude-sonnet-4-6", "claude-haiku-4-5"],
        "tools": True, "vision": True, "prompt_cache": True,
    },
    "openai": {
        "label": "OpenAI", "prefix": "openai/", "needs_key": True,
        "models": ["gpt-4o", "gpt-4o-mini", "o4-mini"],
        "tools": True, "vision": True,
    },
    "ollama": {
        "label": "Ollama (local, free)", "prefix": "ollama_chat/", "needs_key": False,
        "needs_base": True, "editable_model": True,
        "models": ["llama3.1", "qwen2.5-coder", "mistral"],
        "tools": None, "vision": None,   # depends on the local model
    },
}
DEFAULT_PROVIDER = "anthropic"
DEFAULT_OLLAMA_BASE = "http://localhost:11434"


def _keyring():
    try:
        import keyring
        return keyring
    except Exception:
        return None


class AssistantConfig:
    def __init__(self, settings=None):
        self._s = settings or QSettings("XSlope", "XSlope Studio")

    # --- selection -------------------------------------------------------
    def provider(self):
        p = self._s.value("ai/provider", DEFAULT_PROVIDER)
        return p if p in PROVIDERS else DEFAULT_PROVIDER

    def model(self):
        prov = self.provider()
        return self._s.value(f"ai/model/{prov}", PROVIDERS[prov]["models"][0])

    def ollama_base(self):
        return self._s.value("ai/ollama_base", DEFAULT_OLLAMA_BASE)

    def set_selection(self, provider, model, ollama_base=None):
        self._s.setValue("ai/provider", provider)
        self._s.setValue(f"ai/model/{provider}", model)
        if ollama_base is not None:
            self._s.setValue("ai/ollama_base", ollama_base)

    def confirm_before_run(self):
        v = self._s.value("ai/confirm", True)
        return v if isinstance(v, bool) else str(v).lower() in ("true", "1")

    def set_confirm_before_run(self, value):
        self._s.setValue("ai/confirm", bool(value))

    # --- credentials -----------------------------------------------------
    def api_key(self, provider):
        kr = _keyring()
        if kr is not None:
            try:
                key = kr.get_password(KEYRING_SERVICE, f"{provider}_api_key")
                if key:
                    return key
            except Exception:
                pass
        return self._s.value(f"ai/key_fallback/{provider}", "") or ""

    def set_api_key(self, provider, key):
        kr = _keyring()
        if kr is not None:
            try:
                if key:
                    kr.set_password(KEYRING_SERVICE, f"{provider}_api_key", key)
                else:
                    kr.delete_password(KEYRING_SERVICE, f"{provider}_api_key")
                self._s.remove(f"ai/key_fallback/{provider}")   # don't keep a stale copy
                return
            except Exception:
                pass
        # No keyring backend — fall back to QSettings (less secure).
        self._s.setValue(f"ai/key_fallback/{provider}", key)

    def keyring_available(self):
        return _keyring() is not None

    # --- derived ---------------------------------------------------------
    def is_ready(self):
        spec = PROVIDERS[self.provider()]
        if spec["needs_key"] and not self.api_key(self.provider()):
            return False
        return True

    def display_name(self):
        return f"{PROVIDERS[self.provider()]['label']} · {self.model()}"

    def capabilities(self):
        """Capability hints for the current provider: ``{tools, vision}`` each
        True / False / None (None = depends on the chosen local model)."""
        spec = PROVIDERS[self.provider()]
        return {"tools": spec.get("tools"), "vision": spec.get("vision")}

    def supports_prompt_cache(self):
        """Whether to mark the system prompt cacheable (Anthropic only)."""
        return bool(PROVIDERS[self.provider()].get("prompt_cache"))

    def completion_kwargs(self):
        """The provider-specific kwargs for ``litellm.completion(...)``."""
        prov = self.provider()
        spec = PROVIDERS[prov]
        kw = {"model": spec["prefix"] + self.model()}
        if spec["needs_key"]:
            kw["api_key"] = self.api_key(prov)
        if spec.get("needs_base"):
            kw["api_base"] = self.ollama_base()
        return kw
