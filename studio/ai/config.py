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
        # Gets the full skill; OpenAI caches the prompt prefix server-side
        # automatically, so repeat turns are cheap.
        "tools": True, "vision": True, "skill": True,
    },
    "deepseek": {
        "label": "DeepSeek", "prefix": "deepseek/", "needs_key": True,
        "models": ["deepseek-chat"],   # V3 chat: function-calling capable
        # Text-only (no vision). Cheap + caches the prompt prefix server-side, so
        # it gets the full skill (`skill`) without Anthropic cache_control blocks.
        "tools": True, "vision": False, "skill": True,
    },
    "zai": {
        "label": "Z.ai (GLM)", "prefix": "openai/", "needs_key": True,
        # OpenAI-compatible endpoint (editable: GLM coding-plan keys use a
        # different base — .../api/coding/paas/v4). Gets the full skill.
        "base": "https://api.z.ai/api/paas/v4",
        # Suggestions only — model is editable so you can type any current GLM id
        # (the lineup changes fast). Text models first, then the V (vision) models.
        # GLM-OCR is omitted: OCR-only, no function calling, so it can't drive the
        # assistant (which needs tool calls).
        "editable_model": True,
        "models": ["glm-4.6", "glm-4.7", "glm-5.2", "glm-5.1", "glm-5",
                   "glm-5-turbo", "glm-4.7-flashx", "glm-4.7-flash", "glm-4.5",
                   "glm-4.5-x", "glm-4.5-air", "glm-4.5-airx", "glm-4.5-flash",
                   "glm-4-32b-0414-128k",
                   "glm-5v-turbo", "glm-4.6v", "glm-4.6v-flashx", "glm-4.6v-flash",
                   "glm-4.5v"],
        # vision depends on the chosen model — GLM vision models carry a "V" right
        # after the version (glm-4.6v, glm-5v-turbo). tools across all.
        "tools": True, "vision": None, "skill": True,
        "vision_match": r"\dv",
    },
    "ollama": {
        "label": "Ollama (local, free)", "prefix": "ollama_chat/", "needs_key": False,
        "base": "http://localhost:11434", "editable_model": True,
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

    def base_url(self, provider=None):
        """The API base URL for an OpenAI-compatible / local provider (those with a
        `base` in PROVIDERS), or None. Per-provider override falls back to the spec
        default (and, for Ollama, to the legacy single-key setting)."""
        prov = provider or self.provider()
        default = PROVIDERS[prov].get("base")
        if default is None:
            return None
        if prov == "ollama":
            default = self._s.value("ai/ollama_base", default)   # migrate old key
        return self._s.value(f"ai/base/{prov}", default)

    def set_base_url(self, provider, url):
        self._s.setValue(f"ai/base/{provider}", url)

    def ollama_base(self):
        return self.base_url("ollama")

    def set_selection(self, provider, model, base=None):
        self._s.setValue("ai/provider", provider)
        self._s.setValue(f"ai/model/{provider}", model)
        if base is not None and PROVIDERS[provider].get("base") is not None:
            self.set_base_url(provider, base)

    def confirm_before_run(self):
        v = self._s.value("ai/confirm", True)
        return v if isinstance(v, bool) else str(v).lower() in ("true", "1")

    def set_confirm_before_run(self, value):
        self._s.setValue("ai/confirm", bool(value))

    # --- credentials -----------------------------------------------------
    def api_key(self, provider):
        kr = _keyring()
        raw = None
        if kr is not None:
            try:
                raw = kr.get_password(KEYRING_SERVICE, f"{provider}_api_key")
            except Exception:
                pass
        if not raw:
            raw = self._s.value(f"ai/key_fallback/{provider}", "") or ""
        key = (raw or "").strip()
        # A real API key has no internal whitespace/newlines; reject a corrupt
        # value (e.g. text pasted by mistake) so it can't become an illegal HTTP
        # header — treat it as unset instead.
        if key and any(ch.isspace() for ch in key):
            return ""
        return key

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
        """Capability hints for the current provider/model: ``{tools, vision}`` each
        True / False / None (None = depends on the chosen local model). When the
        provider lists ``vision_models`` (e.g. Z.ai's GLM-V models), vision is
        resolved from the selected model rather than the provider default."""
        spec = PROVIDERS[self.provider()]
        vision = spec.get("vision")
        pat = spec.get("vision_match")
        if pat is not None:
            import re
            vision = bool(re.search(pat, (self.model() or "").lower()))
        return {"tools": spec.get("tools"), "vision": vision}

    def supports_prompt_cache(self):
        """Whether to mark the system prompt cacheable with cache_control blocks
        (Anthropic only)."""
        return bool(PROVIDERS[self.provider()].get("prompt_cache"))

    def wants_skill(self):
        """Whether to send the full xslope skill body in the system prompt. True
        for providers that can afford it — Anthropic (prompt-cached) and any
        provider flagged ``skill`` (e.g. DeepSeek: cheap + server-side prefix
        caching). Others get the compact prompt to keep per-turn cost/latency low."""
        spec = PROVIDERS[self.provider()]
        return bool(spec.get("prompt_cache") or spec.get("skill"))

    def completion_kwargs(self):
        """The provider-specific kwargs for ``litellm.completion(...)``."""
        prov = self.provider()
        spec = PROVIDERS[prov]
        kw = {"model": spec["prefix"] + self.model()}
        if spec["needs_key"]:
            kw["api_key"] = self.api_key(prov)
        if spec.get("base"):
            kw["api_base"] = self.base_url(prov)
        return kw
