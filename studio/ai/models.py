"""Where the assistant's model list comes from — enumeration, cache, manifest.

A model list baked into a release goes stale the week after it ships: providers
retire ids, rename aliases, and publish new flagships far faster than Studio is
rebuilt. So the list the settings dialog offers is assembled at open time, from
three sources in a strict order:

1. **Live** — the selected provider's own list-models endpoint, read with the
   user's stored key over the same ``urllib`` path the updater uses (Anthropic
   ``GET /v1/models``, OpenAI ``GET /v1/models``, the OpenAI-compatible
   ``{base}/models`` for Kimi and Z.ai, Ollama's ``GET /api/tags``).
2. **Cache** — whatever the last successful enumeration returned for that
   provider, kept in the app's own settings store, so an offline laptop still
   shows a current list.
3. **Static** — the ``models`` list in :mod:`studio.ai.config`, which is now only
   the never-connected fallback.

On top of whichever list is showing, a curated **recommendations manifest** is
overlaid: ``assistant_models.json``, published as a release asset beside
``latest.json`` and fetched at most once a day. It marks one recommended model
per provider (the default for a fresh install), sorts a handful of good choices
to the top with friendly labels, and marks superseded ids as deprecated —
without hiding anything, and without an app release.

**Nothing here raises at the user, and nothing here is required.** Every fetch
has a short timeout and a silent-failure wrapper; a dead network, a bad key, a
404 and a corrupt manifest all end at the next source in the chain. The dialog's
model box is editable in every case, so a model published an hour ago can always
be typed in by hand.

``XSLOPE_NO_UPDATE_CHECK`` suppresses the manifest fetch along with the version
check. ``XSLOPE_ASSISTANT_MODELS_URL`` points it somewhere else (a beta channel,
or a ``file://`` fixture) without touching code.
"""

from __future__ import annotations

import json
import os
import re
import urllib.request

from .config import PROVIDERS, model_is_vision
from ..updater import REPO, TIMEOUT, due_for_check, utc_stamp

__all__ = [
    "LISTING", "MANIFEST_NAME", "MANIFEST_URL", "SUPPORTED_SCHEMA", "TIMEOUT",
    "ModelChoice", "listing_url", "auth_headers", "parse_models", "is_chat_model",
    "is_offerable", "list_models", "fetch_models", "manifest_url",
    "fetch_manifest", "valid_manifest", "load_manifest", "store_manifest", "manifest_due",
    "load_cached", "save_cached", "resolve_models", "model_choices",
    "recommended_model", "default_model",
]

#: Per-provider list-models endpoint. ``base`` is the default host; providers that
#: carry a user-editable base in :data:`studio.ai.config.PROVIDERS` (Z.ai, Ollama)
#: use that instead, which is what makes a custom OpenAI-compatible host work.
#: ``style`` names the response shape :func:`parse_models` decodes.
LISTING = {
    "anthropic": {"base": "https://api.anthropic.com", "path": "/v1/models",
                  "style": "openai"},          # {"data": [{"id": ...}, ...]}
    "openai": {"base": "https://api.openai.com", "path": "/v1/models",
               "style": "openai"},
    "kimi": {"base": "https://api.moonshot.ai/v1", "path": "/models",
             "style": "openai"},
    "zai": {"base": None, "path": "/models", "style": "openai"},
    "ollama": {"base": "http://localhost:11434", "path": "/api/tags",
               "style": "ollama"},             # {"models": [{"name": ...}, ...]}
}

#: The Anthropic API version header. Its own docs make this mandatory on every
#: request, including the models list.
ANTHROPIC_VERSION = "2023-06-01"

#: The curated manifest, published as a release asset beside ``latest.json`` and
#: read through the same ``latest/download`` redirect.
MANIFEST_NAME = "assistant_models.json"
MANIFEST_URL = f"https://github.com/{REPO}/releases/latest/download/{MANIFEST_NAME}"
#: The only ``assistant_models.json`` schema this build knows how to read.
SUPPORTED_SCHEMA = 1

# Settings keys (the app's existing QSettings store; new keys only — the
# provider/model/credential identities in config.py are untouched).
SETTING_CACHE = "ai/model_cache/{provider}"
SETTING_MANIFEST = "ai/models_manifest"
SETTING_LAST_MANIFEST = "ai/models_manifest_check"

#: Ids that are not chat models, by provider, matched case-insensitively against
#: the id. Deliberately conservative: a provider list is only worth filtering
#: where it is noisy, and anything not recognised is KEPT — a hidden model the
#: user wanted is a worse failure than one extra row in a combo box.
_NOT_CHAT = {
    "openai": re.compile(
        r"(embedding|whisper|^tts-|-tts|dall-e|gpt-image|sora|moderation|"
        r"transcribe|^babbage|^davinci|realtime)", re.I),
    # Z.ai's catalogue carries image/video/OCR/rerank families alongside the chat
    # GLMs. GLM-OCR has no function calling, so it cannot drive the assistant.
    "zai": re.compile(r"(embedding|rerank|ocr|cogview|cogvideo|vidu|asr|tts|audio)", re.I),
    # An Ollama tag list is whatever the user pulled; only embedding models are
    # unusable as a chat model.
    "ollama": re.compile(r"(embed|bge-|nomic-)", re.I),
    # Moonshot lists chat models only, and kimi-thinking-preview — the one id with
    # no tool support at all, so it cannot drive the assistant.
    "kimi": re.compile(r"(kimi-thinking-preview|embedding|moderation)", re.I),
    # Anthropic lists chat models only — nothing to filter.
}


# ---------------------------------------------------------------------------
# The endpoint
# ---------------------------------------------------------------------------

def listing_url(provider, base_url=None):
    """The list-models URL for ``provider``.

    ``base_url`` (the user's own host for the OpenAI-compatible and local
    providers) wins over the spec default, so a GLM coding-plan base or a remote
    Ollama is enumerated from where it actually lives.
    """
    spec = LISTING.get(provider)
    if spec is None:
        return None
    base = (base_url or spec["base"] or "").strip()
    if not base:
        return None
    return base.rstrip("/") + spec["path"]


def auth_headers(provider, api_key=None):
    """The headers this provider's list endpoint needs.

    Anthropic authenticates with ``x-api-key`` plus its version header; the
    OpenAI-compatible hosts use ``Authorization: Bearer``; Ollama is local and
    needs nothing.
    """
    headers = {"User-Agent": _user_agent(), "Accept": "application/json"}
    key = (api_key or "").strip()
    if provider == "anthropic":
        headers["anthropic-version"] = ANTHROPIC_VERSION
        if key:
            headers["x-api-key"] = key
    elif provider != "ollama" and key:
        headers["Authorization"] = f"Bearer {key}"
    return headers


def _user_agent():
    try:
        from xslope import __version__ as v
    except Exception:
        v = "unknown"
    import sys
    return f"XSLOPE-Studio/{v} ({sys.platform})"


# ---------------------------------------------------------------------------
# The response
# ---------------------------------------------------------------------------

def parse_models(provider, payload):
    """Model ids out of one provider's list response. Never raises.

    Anything that is not the shape this provider documents — a list where an
    object was promised, an entry with no id — is skipped rather than failing the
    whole enumeration.
    """
    style = (LISTING.get(provider) or {}).get("style", "openai")
    if not isinstance(payload, dict):
        return []
    ids = []
    if style == "ollama":
        rows = payload.get("models")
        key_order = ("model", "name")
    else:
        rows = payload.get("data")
        key_order = ("id",)
    if not isinstance(rows, list):
        return []
    for row in rows:
        value = None
        if isinstance(row, dict):
            for key in key_order:
                if isinstance(row.get(key), str) and row[key].strip():
                    value = row[key].strip()
                    break
        elif isinstance(row, str):
            value = row.strip()
        if value and value not in ids:
            ids.append(value)
    return ids


def is_chat_model(provider, model_id):
    """Whether ``model_id`` looks like something the assistant could talk to.

    Conservative by construction: only ids matching a documented non-chat naming
    pattern for that provider are rejected (embeddings, speech, image and video
    families). Everything else — including ids this build has never seen — is
    kept.
    """
    pattern = _NOT_CHAT.get(provider)
    if pattern is None:
        return True
    return not pattern.search(str(model_id or ""))


def is_offerable(provider, model_id):
    """Whether ``model_id`` belongs in this provider's list at all.

    A chat model (:func:`is_chat_model`) and — for the providers whose catalogue
    is mixed (``vision_only`` in :data:`studio.ai.config.PROVIDERS`: Kimi, Z.ai,
    Ollama) — one that can read an image. The assistant's headline request is a
    cross section handed over as a picture, and a text-only model turns that into
    a conversation about what the picture shows, so those catalogues are offered
    as the part of themselves that can do the job rather than in full.

    The filter is a NAME test, so it is as conservative as the naming convention
    it reads: a provider that renames its vision family drops out of the list
    rather than being mislabelled, and the model box stays editable, so an id
    published this morning can always be typed in.
    """
    if not is_chat_model(provider, model_id):
        return False
    if not (PROVIDERS.get(provider) or {}).get("vision_only"):
        return True
    return bool(model_is_vision(provider, model_id))


def list_models(provider, api_key=None, base_url=None, timeout=TIMEOUT, url=None):
    """GET, parse and filter one provider's model list. **Raises** on failure.

    See :func:`fetch_models` for the version that never does.
    """
    target = url or listing_url(provider, base_url)
    if not target:
        raise ValueError(f"no list-models endpoint for provider {provider!r}")
    req = urllib.request.Request(target, headers=auth_headers(provider, api_key))
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read(1 << 20)          # a model list is a few KB; cap the read
    payload = json.loads(raw.decode("utf-8"))
    return [m for m in parse_models(provider, payload) if is_offerable(provider, m)]


def fetch_models(provider, api_key=None, base_url=None, timeout=TIMEOUT, url=None):
    """One enumeration attempt, as ``(models, error)``. NEVER raises.

    ``models`` is ``None`` when nothing could be read — offline, a rejected key, a
    host that answered something other than JSON — and ``error`` is the short
    sentence the dialog puts under the combo box. An empty list is also treated
    as a failure: a provider that returned no models cannot replace a cached list
    that has some.
    """
    try:
        models = list_models(provider, api_key=api_key, base_url=base_url,
                             timeout=timeout, url=url)
    except Exception as exc:
        return None, _reason(exc)
    if not models:
        return None, "The provider returned no models."
    return models, ""


def _reason(exc):
    status = getattr(exc, "code", None)
    if status in (401, 403):
        return "The provider rejected that API key."
    if status == 404:
        return "The provider has no model list at that address."
    return f"Could not reach the provider ({exc.__class__.__name__})."


# ---------------------------------------------------------------------------
# The cache (last known list, per provider)
# ---------------------------------------------------------------------------

def load_cached(settings, provider):
    """The last successful enumeration for ``provider``, or ``None``.

    ``settings`` is anything with ``value``/``setValue`` — the app's ``QSettings``
    in Studio, a stub in the checks.
    """
    if settings is None:
        return None
    raw = settings.value(SETTING_CACHE.format(provider=provider), "")
    models = _decode_list(raw)
    return models or None


def save_cached(settings, provider, models):
    """Remember ``models`` as this provider's last known list."""
    if settings is None or not models:
        return
    settings.setValue(SETTING_CACHE.format(provider=provider),
                      json.dumps([str(m) for m in models]))


def _decode_list(raw):
    if isinstance(raw, list):
        return [str(m) for m in raw if str(m).strip()]
    if not isinstance(raw, str) or not raw.strip():
        return []
    try:
        data = json.loads(raw)
    except Exception:
        return []
    return [str(m) for m in data if str(m).strip()] if isinstance(data, list) else []


# ---------------------------------------------------------------------------
# The recommendations manifest
# ---------------------------------------------------------------------------

def manifest_url():
    """The manifest to read. ``XSLOPE_ASSISTANT_MODELS_URL`` overrides the
    release default."""
    return os.environ.get("XSLOPE_ASSISTANT_MODELS_URL") or MANIFEST_URL


def fetch_manifest(url=None, timeout=TIMEOUT):
    """GET and parse the recommendations manifest. Raises on any failure."""
    req = urllib.request.Request(url or manifest_url(),
                                 headers={"User-Agent": _user_agent()})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read(1 << 20)
    return json.loads(raw.decode("utf-8"))


def valid_manifest(obj):
    """``obj`` if it is a manifest this build understands, else ``None``.

    A manifest is advice, so anything doubtful is simply not applied: a wrong
    schema number, a missing ``providers`` map, junk where a list belongs. The
    dialog then shows exactly the list it would have shown with no manifest at
    all, which is the whole point of keeping the overlay separate from the list.
    """
    if not isinstance(obj, dict) or obj.get("schema") != SUPPORTED_SCHEMA:
        return None
    providers = obj.get("providers")
    if not isinstance(providers, dict):
        return None
    return obj


def manifest_due(settings, now=None):
    """Is the once-a-day manifest fetch allowed to run?

    ``XSLOPE_NO_UPDATE_CHECK`` suppresses it exactly as it suppresses the version
    check — one switch turns off everything this app reads from the network.
    """
    if os.environ.get("XSLOPE_NO_UPDATE_CHECK"):
        return False
    last = settings.value(SETTING_LAST_MANIFEST, "") if settings is not None else ""
    return due_for_check(last, now=now)


def load_manifest(settings):
    """The stored manifest, or ``{}``. No network."""
    if settings is None:
        return {}
    raw = settings.value(SETTING_MANIFEST, "")
    if isinstance(raw, dict):
        return valid_manifest(raw) or {}
    if not isinstance(raw, str) or not raw.strip():
        return {}
    try:
        return valid_manifest(json.loads(raw)) or {}
    except Exception:
        return {}


def store_manifest(settings, manifest, now=None):
    """Stamp the check and keep ``manifest`` if it is readable.

    The stamp is written whatever the outcome — an unreachable or malformed
    manifest must not make the app retry on every dialog open — and a bad
    manifest never replaces a good one already on disk.
    """
    if settings is None:
        return False
    settings.setValue(SETTING_LAST_MANIFEST, utc_stamp(now))
    good = valid_manifest(manifest)
    if good is None:
        return False
    settings.setValue(SETTING_MANIFEST, json.dumps(good))
    return True


def provider_advice(manifest, provider):
    """The manifest's entry for one provider, normalised. Unknown provider,
    unknown keys and wrong types all come back as empty advice."""
    empty = {"recommended": "", "good": [], "deprecated": []}
    if not isinstance(manifest, dict):
        return empty
    providers = manifest.get("providers")
    entry = providers.get(provider) if isinstance(providers, dict) else None
    if not isinstance(entry, dict):
        return empty
    rec = entry.get("recommended")
    good = []
    rows = entry.get("good_choices")
    # A string is iterable and a dict yields its keys — a field of the wrong type
    # must contribute nothing rather than a model id per character.
    for row in rows if isinstance(rows, list) else []:
        if isinstance(row, dict) and isinstance(row.get("id"), str):
            good.append({"id": row["id"].strip(),
                         "label": str(row.get("label") or "").strip(),
                         "note": str(row.get("note") or "").strip()})
        elif isinstance(row, str) and row.strip():
            good.append({"id": row.strip(), "label": "", "note": ""})
    dropped = entry.get("deprecated")
    dep = [d.strip() for d in (dropped if isinstance(dropped, list) else [])
           if isinstance(d, str) and d.strip()]
    return {"recommended": rec.strip() if isinstance(rec, str) else "",
            "good": good, "deprecated": dep}


def recommended_model(provider, manifest):
    """The manifest's recommended id for ``provider``, or ``""``."""
    return provider_advice(manifest, provider)["recommended"]


def default_model(provider, settings=None):
    """What a fresh install should select for ``provider``.

    The manifest's recommendation when one has been fetched, otherwise the first
    entry of the shipped static list. No network.
    """
    rec = recommended_model(provider, load_manifest(settings))
    if rec:
        return rec
    static = PROVIDERS.get(provider, {}).get("models") or [""]
    return static[0]


# ---------------------------------------------------------------------------
# The fallback chain, and what the combo box shows
# ---------------------------------------------------------------------------

def resolve_models(provider, live=None, cached=None, static=None):
    """The list to show, and where it came from: ``(models, source)``.

    Strictly ordered — live, then cache, then the shipped static list — with
    ``source`` one of ``"live"``, ``"cache"``, ``"static"``, ``"none"``. An empty
    or missing source is skipped, never shown: an empty combo box is the one
    outcome that leaves a user with nothing to click.
    """
    if static is None:
        static = PROVIDERS.get(provider, {}).get("models") or []
    for models, source in ((live, "live"), (cached, "cache"), (static, "static")):
        if models:
            out, seen = [], set()
            for m in models:
                m = str(m).strip()
                if m and m not in seen:
                    seen.add(m)
                    out.append(m)
            if out:
                return out, source
    return [], "none"


class ModelChoice:
    """One row of the model box: the id that gets saved, and how it is shown."""

    def __init__(self, model_id, label="", note="", recommended=False,
                 deprecated=False):
        self.id = model_id
        self.label = label
        self.note = note
        self.recommended = recommended
        self.deprecated = deprecated

    @property
    def text(self):
        """The combo-box text. The id always leads, so the list stays scannable
        and typing an id still finds its row; the marking follows it."""
        parts = [self.id]
        tail = []
        if self.label:
            tail.append(self.label)
        if self.recommended:
            tail.append("recommended")
        elif self.deprecated:
            tail.append("superseded")
        if tail:
            parts.append("— " + ", ".join(tail))
        return " ".join(parts)

    def __repr__(self):
        return f"ModelChoice({self.id!r}, recommended={self.recommended}, " \
               f"deprecated={self.deprecated})"


def model_choices(provider, models, manifest=None):
    """Apply the manifest overlay to ``models`` and order the result.

    The recommended model first, then the manifest's other good choices in the
    order it lists them, then everything else the provider offers, then the
    deprecated ids last — still present, still selectable, and marked. Curated
    ids the list does not contain are added rather than dropped, so a fresh
    install can pick the recommended model before it has ever enumerated
    anything.
    """
    advice = provider_advice(manifest, provider)
    labels = {row["id"]: row for row in advice["good"]}
    order = [row["id"] for row in advice["good"]]
    rec = advice["recommended"]
    deprecated = set(advice["deprecated"])

    known = []
    for m in models or []:
        m = str(m).strip()
        if m and m not in known:
            known.append(m)
    for curated in ([rec] if rec else []) + order:
        if curated and curated not in known:
            known.append(curated)

    def rank(model_id):
        if model_id == rec:
            return (0, 0)
        if model_id in deprecated:
            return (3, 0)
        if model_id in order:
            return (1, order.index(model_id))
        return (2, 0)

    out = []
    for model_id in known:
        row = labels.get(model_id, {})
        out.append(ModelChoice(model_id, label=row.get("label", ""),
                               note=row.get("note", ""),
                               recommended=bool(rec) and model_id == rec,
                               deprecated=model_id in deprecated))
    out.sort(key=lambda c: rank(c.id) + (known.index(c.id),))
    return out
