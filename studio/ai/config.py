"""AssistantConfig — provider/model selection and credential storage.

Holds the user's choice of provider (Claude / OpenAI / Kimi / Z.ai / local
Ollama) and model in ``QSettings``, and API keys in the OS keychain via
``keyring`` (falling back to QSettings if keyring is unavailable). Produces the
kwargs for a ``litellm`` completion call, so the agent loop stays
provider-agnostic.

**Every provider offered here can read an image.** The assistant's headline case
is a cross section handed over as a sketch or a photograph, and a model that
cannot see one turns that request into a conversation about what the picture
shows. So a text-only provider is not offered at a lower tier — it is not
offered. Where a provider's catalogue is mixed (Z.ai's GLM-V models beside its
text GLMs, a local Ollama library that is mostly text), ``vision_only`` filters
the list down to the models that can, and the settings dialog says so.
"""

from __future__ import annotations

import re

from PySide6.QtCore import QSettings

# Storage identity, not a display name: this string is the key under which the
# OS keychain already holds every user's API keys, exactly as the QSettings pair
# below is their settings path. Both stay spelled this way even though the app
# displays "XSLOPE Studio" — renaming either one would orphan existing data.
KEYRING_SERVICE = "XSlope Studio"

# Provider catalog. `prefix` is the LiteLLM model-name prefix; `models` is the
# NEVER-CONNECTED FALLBACK list — the settings dialog normally builds its list
# from the provider's own list-models endpoint (or the last one it saw), and only
# falls back to these; see studio/ai/models.py. Every model box is editable, so
# any id can be typed whatever the list says.
# `tools`/`vision` are capability defaults (True/False, or None = depends on the
# chosen model, resolved from `vision_match` against the model id).
# `vision_only` means the dialog LISTS only the models whose id matches
# `vision_match` — the filter that keeps a mixed catalogue offering only models
# this assistant can hand a sketch to.
PROVIDERS = {
    "anthropic": {
        "label": "Claude (Anthropic)", "prefix": "anthropic/", "needs_key": True,
        "models": ["claude-opus-5", "claude-sonnet-5", "claude-haiku-4-5"],
        "tools": True, "vision": True, "prompt_cache": True,
    },
    "openai": {
        "label": "OpenAI", "prefix": "openai/", "needs_key": True,
        "models": ["gpt-5.5", "gpt-5.4", "gpt-5.4-mini", "gpt-5.4-nano"],
        # Gets the Studio brief; OpenAI caches the prompt prefix server-side
        # automatically, so repeat turns are cheap.
        "tools": True, "vision": None, "skill": True,
        # OpenAI's catalogue mixes vision chat models with text-only ones
        # (gpt-3.5, the dated gpt-4 snapshots, o1-mini / o3-mini, the search
        # previews, bare chat-latest). Image input: every gpt-4o / gpt-4.1 /
        # gpt-4-turbo / gpt-5.x model, and the o1 / o3 / o4-mini reasoners.
        "vision_only": True,
        "vision_match": r"^(gpt-4o(?!.*search)|gpt-4\.1|gpt-4-turbo|gpt-5(?!.*search)|"
                        r"o1(?!-mini)|o3(?!-mini)|o4-mini|chatgpt-4o)",
    },
    "kimi": {
        # Moonshot AI's OpenAI-compatible endpoint. litellm routes it natively
        # under the `moonshot/` prefix (its default base is this same URL), so
        # the base stays editable only for the .cn host and for proxies.
        "label": "Kimi (Moonshot AI)", "prefix": "moonshot/", "needs_key": True,
        "base": "https://api.moonshot.ai/v1", "editable_model": True,
        # The K-series models that accept images, newest first, then the
        # moonshot-v1 vision previews. The text-only ids Moonshot also lists
        # (moonshot-v1-8k, kimi-k2-0711-preview, kimi-k2-thinking) are filtered
        # out by `vision_only` wherever the list comes from.
        "models": ["kimi-k2.6", "kimi-k2.5", "kimi-k3", "kimi-latest",
                   "kimi-latest-128k", "moonshot-v1-128k-vision-preview",
                   "moonshot-v1-32k-vision-preview"],
        # Moonshot prices a cache-hit input token separately, i.e. it caches the
        # prompt prefix server-side like OpenAI and needs no cache_control
        # blocks — so it can afford the Studio brief.
        "tools": True, "vision": None, "skill": True,
        "vision_only": True,
        "vision_match": r"(kimi-latest|kimi-k2\.[5-9]|kimi-k[3-9]|vision)",
    },
    "zai": {
        "label": "Z.ai (GLM)", "prefix": "openai/", "needs_key": True,
        # OpenAI-compatible endpoint (editable: GLM coding-plan keys use a
        # different base — .../api/coding/paas/v4). Gets the Studio brief.
        "base": "https://api.z.ai/api/paas/v4",
        # Suggestions only — model is editable so you can type any current GLM id
        # (the lineup changes fast). The V (vision) models ONLY: the text GLMs
        # cannot read a cross section, and GLM-OCR — the other model that takes an
        # image — has no function calling, so it cannot drive the assistant.
        "editable_model": True,
        "models": ["glm-5v-turbo", "glm-4.6v", "glm-4.6v-flashx",
                   "glm-4.6v-flash", "glm-4.5v"],
        # A GLM vision model carries a "V" right after the version number
        # (glm-4.6v, glm-5v-turbo); that is what both the list filter and the
        # per-model capability read.
        "tools": True, "vision": None, "skill": True,
        "vision_only": True,
        "vision_match": r"\dv",
    },
    "ollama": {
        "label": "Ollama (local, free)", "prefix": "ollama_chat/", "needs_key": False,
        "base": "http://localhost:11434", "editable_model": True,
        # Local vision models, by family. An Ollama library is whatever the user
        # pulled, so this is a starting list and the tag list is filtered to the
        # same families — a text-only local model would leave the assistant
        # unable to read the sketch it is being shown.
        "models": ["llava", "llama3.2-vision", "gemma3", "qwen2.5vl",
                   "minicpm-v"],
        # Tool use still depends on the local model (many vision builds have no
        # function calling); vision is resolved from the tag name.
        "tools": None, "vision": None,
        "vision_only": True,
        "vision_match": (r"(llava|bakllava|vision|moondream|minicpm-v|gemma3|"
                         r"qwen2\.?5-?vl|qwen3-?vl|internvl|llama4|"
                         r"granite3\.\d+-vision|mistral-small3)"),
    },
}
DEFAULT_PROVIDER = "anthropic"
DEFAULT_OLLAMA_BASE = "http://localhost:11434"


def model_is_vision(provider, model):
    """Whether ``model`` can accept images, for ``provider``.

    The provider default when it has one, otherwise the ``vision_match`` pattern
    read against the model id — which is what makes a typed-in id answer the
    question too, rather than inheriting the capability of whatever was listed.
    """
    spec = PROVIDERS.get(provider) or {}
    pat = spec.get("vision_match")
    if pat is None:
        return spec.get("vision")
    return bool(re.search(pat, str(model or "").lower()))


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
        """The chosen model — or, on a fresh install, the recommended one.

        With nothing stored yet the answer comes from the hosted recommendations
        manifest if one has been fetched, and from the shipped fallback list if
        not (``models.default_model``). Imported lazily: models.py reads this
        module's PROVIDERS table.
        """
        prov = self.provider()
        stored = str(self._s.value(f"ai/model/{prov}", "") or "").strip()
        if stored:
            return stored
        from .models import default_model
        return default_model(prov, self._s)

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
        v = self._s.value("ai/confirm", False)
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
        """Capability flags for the current provider/model.

        ``{tools, vision, prompt_cache, skill}``. ``tools`` and ``vision`` are
        True / False / None (None = depends on the chosen local model); where the
        provider carries a ``vision_match`` pattern (Kimi, Z.ai's GLM-V models, a
        local Ollama tag) vision is resolved from the SELECTED model rather than
        the provider default, so a typed-in id answers for itself.
        ``prompt_cache`` is whether the system prompt is sent in a cache_control
        block (Anthropic), and ``skill`` whether the Studio reference brief is sent
        at all — the two questions :meth:`supports_prompt_cache` and
        :meth:`wants_skill` answer, reported here so one call describes the
        selection.
        """
        spec = PROVIDERS[self.provider()]
        return {"tools": spec.get("tools"),
                "vision": model_is_vision(self.provider(), self.model()),
                "prompt_cache": self.supports_prompt_cache(),
                "skill": self.wants_skill()}

    def supports_prompt_cache(self):
        """Whether to mark the system prompt cacheable with cache_control blocks
        (Anthropic only)."""
        return bool(PROVIDERS[self.provider()].get("prompt_cache"))

    def wants_skill(self):
        """Whether to send the Studio reference brief in the system prompt.

        True for providers that can afford a reference body — Anthropic
        (prompt-cached) and any provider flagged ``skill`` (Kimi and Z.ai:
        server-side prefix caching, so a repeat turn re-reads it at cache rates).
        Others get the compact prompt, to keep per-turn cost and latency low on a
        local model.

        The name is historical: what this used to gate was the full Claude Code
        skill body, ~34k tokens re-read on every completion of every turn. It now
        gates ``studio_assistant_brief.md``, which says the same things Studio
        actually needs in a fraction of that. The tiers are unchanged — the
        providers that got the skill get the brief.
        """
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
