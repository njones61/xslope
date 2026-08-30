"""Standing checks on where the assistant's model list comes from — the chain
between a provider's catalogue and the box the user picks a model out of.

Nothing here touches the network. Every provider response and the
recommendations manifest are files in a temporary directory served over
``file://``, which is the same ``urllib`` path the real endpoints take, so the
request, the parse and the filter are genuinely exercised — only the host is
local. The Qt leg builds the real dialog with its background refresh turned off,
so opening it in a test can never reach a provider even on a machine that has a
key in its keychain.

What is guarded:

  A. THE ENDPOINT — each provider is asked at the address its own documentation
     gives (Anthropic and OpenAI ``/v1/models``, the OpenAI-compatible
     ``{base}/models``, Ollama ``/api/tags``), with that provider's auth header,
     and a user-supplied base URL overrides the default host — which is the only
     reason a GLM coding-plan account or a remote Ollama can be enumerated at all.
  B. THE PARSE — one decoder per documented response shape, and a response that
     is not that shape yields no models rather than an exception.
  C. THE FILTER — embeddings, speech and image families are dropped where a
     provider's list is noisy, and an id this build has never seen is KEPT. A
     filter that hides a model the user wanted is worse than one extra row.
  D. ENUMERATION, END TO END — a real GET (over ``file://``) per provider shape.
  E. THE FALLBACK CHAIN (mutation) — live, then the cached list, then the list
     the build shipped with. Killing each source in turn must move the answer
     exactly one step down, and an empty combo box must never be the result.
  F. THE MANIFEST OVERLAY — the recommended model is marked and sorts first,
     good choices follow in the manifest's own order with their labels,
     deprecated ids are marked and sink to the bottom but stay selectable, a
     provider the manifest says nothing about is unchanged, and a recommendation
     the list does not contain is ADDED rather than dropped.
  G. THE MALFORMED-MANIFEST GUARD (mutation) — every way a manifest can be wrong
     (bad schema, no providers map, junk types, not JSON at all) must be ignored
     silently, must not replace a good manifest already stored, and must leave
     the model list exactly as it would have been with no manifest.
  H. THE THROTTLE — the manifest is fetched at most once a day, and
     ``XSLOPE_NO_UPDATE_CHECK`` suppresses it entirely.
  I. THE BACKGROUND RUNNER — the QThread the dialog uses does both reads off the
     GUI thread and returns one payload, and a manifest that could not be read
     never takes the model list down with it.
  J. THE DIALOG — the model box is editable, a typed id survives being saved and
     re-read, what gets saved is the model ID and never the "(recommended)"
     decoration, a Refresh button exists, and a live result replaces the fallback
     list without losing what the user had selected.

Skips its Studio leg cleanly (exit 0) when PySide6 is not installed; A-H run
either way.
"""
import json
import os
import shutil
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from studio.ai import models as ml
from studio.ai.config import PROVIDERS


# --------------------------------------------------------------------------
# Fixtures: provider responses and a manifest, on disk, addressed as file:// URLs
# --------------------------------------------------------------------------

def _file_url(path):
    from urllib.request import pathname2url
    return "file://" + pathname2url(os.path.abspath(path))


class Fixture:
    """A temp directory of JSON documents served over ``file://``."""

    def __init__(self):
        self.dir = tempfile.mkdtemp(prefix="xslope-models-test-")

    def write(self, name, payload):
        path = os.path.join(self.dir, name)
        with open(path, "w", encoding="utf-8") as fh:
            if isinstance(payload, str):
                fh.write(payload)
            else:
                json.dump(payload, fh, indent=2)
        return _file_url(path)

    def close(self):
        shutil.rmtree(self.dir, ignore_errors=True)


class Settings:
    """The two QSettings methods this feature uses, without Qt."""

    def __init__(self, initial=None):
        self._d = dict(initial or {})

    def value(self, key, default=None):
        return self._d.get(key, default)

    def setValue(self, key, value):
        self._d[key] = value


#: A manifest with one of everything the dialog reads.
MANIFEST = {
    "schema": 1,
    "updated": "2026-07-31",
    "providers": {
        "anthropic": {
            "recommended": "claude-opus-5",
            "good_choices": [
                {"id": "claude-opus-5", "label": "flagship", "note": "the note"},
                {"id": "claude-haiku-4-5", "label": "economical"},
            ],
            "deprecated": ["claude-opus-4-1", "claude-3-opus-20240229"],
        },
    },
}


# ------------------------------------------------------------- A. the endpoint

def test_endpoints():
    fails = []
    want = {
        "anthropic": "https://api.anthropic.com/v1/models",
        "openai": "https://api.openai.com/v1/models",
        "kimi": "https://api.moonshot.ai/v1/models",
        "ollama": "http://localhost:11434/api/tags",
    }
    for provider, url in want.items():
        got = ml.listing_url(provider)
        if got != url:
            fails.append(f"listing_url({provider!r}) = {got!r}, expected {url!r}")

    # Z.ai has no default host of its own: the base URL the user configured (the
    # one litellm calls) is what gets enumerated.
    if ml.listing_url("zai") is not None:
        fails.append("zai must have no endpoint without a base URL")
    zai = ml.listing_url("zai", base_url=PROVIDERS["zai"]["base"])
    if zai != "https://api.z.ai/api/paas/v4/models":
        fails.append(f"zai base gave {zai!r}")
    # A user-supplied base overrides the default host, trailing slash and all.
    if ml.listing_url("ollama", base_url="http://box.local:11434/") != \
            "http://box.local:11434/api/tags":
        fails.append("a custom Ollama base URL is not honoured")
    if ml.listing_url("openai", base_url="https://proxy.example/v1") != \
            "https://proxy.example/v1/v1/models":
        fails.append("an OpenAI base override is not honoured")
    if ml.listing_url("nope") is not None:
        fails.append("an unknown provider produced an endpoint")

    # Each provider's own auth scheme, and no key leaked into the wrong header.
    h = ml.auth_headers("anthropic", "sk-ant-test")
    if h.get("x-api-key") != "sk-ant-test":
        fails.append(f"Anthropic must authenticate with x-api-key: {sorted(h)}")
    if h.get("anthropic-version") != ml.ANTHROPIC_VERSION:
        fails.append("the Anthropic version header is missing")
    if "Authorization" in h:
        fails.append("Anthropic must not send a bearer token")
    for provider in ("openai", "kimi", "zai"):
        h = ml.auth_headers(provider, "sk-test")
        if h.get("Authorization") != "Bearer sk-test":
            fails.append(f"{provider} must authenticate with a bearer token")
        if "x-api-key" in h:
            fails.append(f"{provider} sent an x-api-key header")
    h = ml.auth_headers("ollama", "")
    if "Authorization" in h or "x-api-key" in h:
        fails.append("the local provider must send no credential at all")
    if not ml.auth_headers("openai", "")["User-Agent"].startswith("XSLOPE-Studio/"):
        fails.append("the request does not identify itself")
    if "Authorization" in ml.auth_headers("openai", ""):
        fails.append("an empty key produced an Authorization header")
    return fails


# ---------------------------------------------------------------- B. the parse

def test_parse():
    fails = []
    # Anthropic / OpenAI / OpenAI-compatible: {"data": [{"id": ...}]}
    payload = {"data": [{"type": "model", "id": "claude-opus-5"},
                        {"type": "model", "id": "claude-haiku-4-5"}],
               "has_more": False}
    got = ml.parse_models("anthropic", payload)
    if got != ["claude-opus-5", "claude-haiku-4-5"]:
        fails.append(f"the Anthropic list parsed as {got}")
    if ml.parse_models("openai", {"object": "list",
                                  "data": [{"id": "gpt-5.1", "object": "model"}]}) \
            != ["gpt-5.1"]:
        fails.append("the OpenAI list did not parse")
    if ml.parse_models("zai", {"data": [{"id": "glm-5.2"}]}) != ["glm-5.2"]:
        fails.append("the OpenAI-compatible list did not parse")

    # Ollama: {"models": [{"name": ...}]} — its own shape, not OpenAI's.
    ollama = {"models": [{"name": "llama3.1:8b", "model": "llama3.1:8b"},
                         {"name": "qwen2.5-coder:7b", "model": "qwen2.5-coder:7b"}]}
    got = ml.parse_models("ollama", ollama)
    if got != ["llama3.1:8b", "qwen2.5-coder:7b"]:
        fails.append(f"the Ollama tag list parsed as {got}")

    # Junk in, nothing out — never an exception.
    for bad in (None, [], "", {"data": "nope"}, {"models": {}}, {},
                {"data": [{}, {"id": ""}, 7, None]}):
        try:
            got = ml.parse_models("openai", bad)
        except Exception as exc:
            fails.append(f"parse_models raised on {bad!r}: {exc!r}")
            continue
        if got != []:
            fails.append(f"parse_models({bad!r}) = {got}, expected []")
    # Duplicates collapse (a paginated list can repeat an id).
    if ml.parse_models("openai", {"data": [{"id": "a"}, {"id": "a"}]}) != ["a"]:
        fails.append("a repeated id was listed twice")
    return fails


# --------------------------------------------------------------- C. the filter

def test_filter():
    fails = []
    drop = [
        ("openai", "text-embedding-3-large"), ("openai", "whisper-1"),
        ("openai", "tts-1-hd"), ("openai", "dall-e-3"),
        ("openai", "gpt-image-1"), ("openai", "omni-moderation-latest"),
        ("openai", "davinci-002"), ("openai", "gpt-4o-transcribe"),
        ("zai", "embedding-3"), ("zai", "glm-ocr"), ("zai", "cogview-4"),
        ("kimi", "kimi-thinking-preview"),
        ("ollama", "nomic-embed-text"), ("ollama", "mxbai-embed-large"),
    ]
    keep = [
        ("openai", "gpt-5.1"), ("openai", "gpt-5-mini"), ("openai", "o4-mini"),
        ("openai", "some-model-invented-next-year"),
        ("anthropic", "claude-opus-5"), ("anthropic", "claude-anything-at-all"),
        ("kimi", "kimi-k2.6"), ("kimi", "some-kimi-published-tomorrow"),
        ("zai", "glm-5.2"), ("zai", "glm-4.6v"),
        ("ollama", "llama3.1:8b"), ("ollama", "qwen2.5-coder:7b"),
    ]
    for provider, model in drop:
        if ml.is_chat_model(provider, model):
            fails.append(f"{provider}: {model!r} is not a chat model but was kept")
    for provider, model in keep:
        if not ml.is_chat_model(provider, model):
            fails.append(f"{provider}: {model!r} was filtered out — when unsure, keep")
    return fails


# ------------------------------------------------------- D. enumeration, live

def test_enumeration():
    fails = []
    fx = Fixture()
    try:
        url = fx.write("anthropic.json", {"data": [
            {"id": "claude-opus-5"}, {"id": "claude-sonnet-5"}]})
        models, error = ml.fetch_models("anthropic", api_key="sk-ant-x", url=url)
        if models != ["claude-opus-5", "claude-sonnet-5"] or error:
            fails.append(f"the Anthropic enumeration gave {models!r} / {error!r}")

        url = fx.write("openai.json", {"data": [
            {"id": "gpt-5.1"}, {"id": "text-embedding-3-small"}, {"id": "gpt-5-mini"}]})
        models, _ = ml.fetch_models("openai", api_key="sk-x", url=url)
        if models != ["gpt-5.1", "gpt-5-mini"]:
            fails.append(f"the OpenAI enumeration did not filter: {models!r}")

        # A local library is mostly text models; only the ones that can read an
        # image are offered, since that is what the assistant is shown.
        url = fx.write("ollama.json", {"models": [
            {"name": "llava:13b"}, {"name": "llama3.1:8b"},
            {"name": "nomic-embed-text:latest"}]})
        models, _ = ml.fetch_models("ollama", url=url)
        if models != ["llava:13b"]:
            fails.append(f"the Ollama enumeration gave {models!r}")

        # Moonshot lists its text K2 models beside the vision ones; the text ones
        # do not reach the box.
        url = fx.write("kimi.json", {"data": [
            {"id": "kimi-k2.6"}, {"id": "moonshot-v1-8k"},
            {"id": "kimi-k2-thinking"}, {"id": "moonshot-v1-32k-vision-preview"}]})
        models, _ = ml.fetch_models("kimi", api_key="sk-x", url=url)
        if models != ["kimi-k2.6", "moonshot-v1-32k-vision-preview"]:
            fails.append(f"the Kimi enumeration gave {models!r}")

        # Failures are answers, not exceptions.
        models, error = ml.fetch_models(
            "openai", url=_file_url(os.path.join(fx.dir, "missing.json")))
        if models is not None or not error:
            fails.append(f"an unreachable endpoint gave {models!r} / {error!r}")
        models, error = ml.fetch_models("openai", url=fx.write("bad.json", "{nope"))
        if models is not None or not error:
            fails.append("a non-JSON response was not reported as a failure")
        # An empty catalogue must not out-rank a cached list that has models.
        models, error = ml.fetch_models("openai", url=fx.write("empty.json",
                                                              {"data": []}))
        if models is not None:
            fails.append("an empty provider list was treated as a live list")
    finally:
        fx.close()
    return fails


# --------------------------------------------- E. the fallback chain (mutation)

def test_fallback_chain():
    fails = []
    static = PROVIDERS["anthropic"]["models"]
    live = ["live-a", "live-b"]
    cached = ["cached-a"]

    # All three present: live wins.
    models, source = ml.resolve_models("anthropic", live=live, cached=cached)
    if (models, source) != (live, "live"):
        fails.append(f"with a live list the chain chose {source!r}: {models}")

    # MUTATION 1: kill the live list. The cache must serve.
    models, source = ml.resolve_models("anthropic", live=None, cached=cached)
    if (models, source) != (cached, "cache"):
        fails.append(f"with no live list the chain chose {source!r}: {models} — "
                     f"the cache did not serve")
    for empty in (None, [], ()):
        models, source = ml.resolve_models("anthropic", live=empty, cached=cached)
        if source != "cache":
            fails.append(f"live={empty!r} did not fall through to the cache")

    # MUTATION 2: kill the cache too. The shipped list must serve.
    models, source = ml.resolve_models("anthropic", live=None, cached=None)
    if source != "static" or models != list(static):
        fails.append(f"with no live list and no cache the chain chose {source!r}: "
                     f"{models} — the shipped list did not serve")

    # MUTATION 3: kill all three. Nothing is offered, and nothing raises.
    models, source = ml.resolve_models("anthropic", live=None, cached=None, static=[])
    if models != [] or source != "none":
        fails.append(f"an empty chain gave {source!r}: {models}")

    # The cache is a real round trip through the settings store.
    settings = Settings()
    if ml.load_cached(settings, "openai") is not None:
        fails.append("an empty settings store returned a cached list")
    ml.save_cached(settings, "openai", ["gpt-5.1", "gpt-5-mini"])
    if ml.load_cached(settings, "openai") != ["gpt-5.1", "gpt-5-mini"]:
        fails.append("the cached list did not survive the settings round trip")
    if ml.load_cached(settings, "anthropic") is not None:
        fails.append("one provider's cache leaked into another's")
    # A corrupt cache entry is a missing one, not a crash.
    settings.setValue(ml.SETTING_CACHE.format(provider="openai"), "{not json")
    if ml.load_cached(settings, "openai") is not None:
        fails.append("a corrupt cache entry was returned as a list")
    _, source = ml.resolve_models("openai", cached=ml.load_cached(settings, "openai"))
    if source != "static":
        fails.append("a corrupt cache did not fall through to the shipped list")

    # Order and duplicates: what the provider listed, in its order, once each.
    models, _ = ml.resolve_models("openai", live=["b", "a", "b", " ", "a"])
    if models != ["b", "a"]:
        fails.append(f"the live list was reordered or duplicated: {models}")
    return fails


# ------------------------------------------------------- F. the manifest overlay

def test_manifest_overlay():
    fails = []
    listed = ["claude-sonnet-5", "claude-opus-4-1", "claude-haiku-4-5",
              "claude-opus-5", "claude-3-opus-20240229"]
    choices = ml.model_choices("anthropic", listed, MANIFEST)
    ids = [c.id for c in choices]

    if ids[0] != "claude-opus-5":
        fails.append(f"the recommended model does not sort first: {ids}")
    rec = choices[0]
    if not rec.recommended or "recommended" not in rec.text:
        fails.append(f"the recommended model is not marked: {rec.text!r}")
    if rec.text.split()[0] != "claude-opus-5":
        fails.append(f"the id must lead the row text: {rec.text!r}")
    if rec.note != "the note":
        fails.append("the manifest's note did not reach the choice")

    # Good choices next, in the manifest's own order, with their labels.
    if ids[1] != "claude-haiku-4-5":
        fails.append(f"the other good choice does not follow: {ids}")
    if "economical" not in choices[1].text:
        fails.append(f"the friendly label is missing: {choices[1].text!r}")

    # Deprecated last, marked, and STILL PRESENT.
    for dep in ("claude-opus-4-1", "claude-3-opus-20240229"):
        if dep not in ids:
            fails.append(f"deprecated {dep!r} was removed instead of marked")
    tail = ids[-2:]
    if sorted(tail) != sorted(["claude-opus-4-1", "claude-3-opus-20240229"]):
        fails.append(f"deprecated ids did not sink to the bottom: {ids}")
    for choice in choices:
        if choice.id in ("claude-opus-4-1", "claude-3-opus-20240229"):
            if not choice.deprecated or "superseded" not in choice.text:
                fails.append(f"{choice.id} is not marked superseded: {choice.text!r}")
    # An ordinary model keeps its bare id as its row text.
    plain = [c for c in choices if c.id == "claude-sonnet-5"][0]
    if plain.text != "claude-sonnet-5" or plain.recommended or plain.deprecated:
        fails.append(f"an unmarked model was decorated: {plain.text!r}")

    # A recommendation the provider list does not contain is ADDED, so a fresh
    # install can select it before it has ever enumerated anything.
    choices = ml.model_choices("anthropic", ["claude-sonnet-5"], MANIFEST)
    if [c.id for c in choices][0] != "claude-opus-5":
        fails.append("a recommendation absent from the list was dropped")

    # A provider the manifest says nothing about passes through untouched.
    listed_ollama = ["llama3.1:8b", "mistral:latest"]
    choices = ml.model_choices("ollama", listed_ollama, MANIFEST)
    if [c.id for c in choices] != listed_ollama:
        fails.append(f"an unknown provider's list was reordered: {choices}")
    if any(c.recommended or c.deprecated or c.text != c.id for c in choices):
        fails.append("an unknown provider's models were marked")

    # No manifest at all: the list, unchanged.
    for nothing in (None, {}, {"schema": 1, "providers": {}}):
        choices = ml.model_choices("anthropic", listed, nothing)
        if [c.id for c in choices] != listed:
            fails.append(f"manifest {nothing!r} altered the list")

    # The recommendation feeds the fresh-install default; a stored choice wins
    # over it (that part is config's, but the default itself lives here).
    settings = Settings()
    if ml.default_model("anthropic", settings) != PROVIDERS["anthropic"]["models"][0]:
        fails.append("with no manifest the default is not the shipped first entry")
    ml.store_manifest(settings, MANIFEST)
    if ml.default_model("anthropic", settings) != "claude-opus-5":
        fails.append("the manifest's recommendation is not the fresh-install default")
    if ml.default_model("openai", settings) != PROVIDERS["openai"]["models"][0]:
        fails.append("a provider the manifest omits lost its shipped default")
    return fails


# --------------------------------------- G. the malformed-manifest guard (mutation)

def test_malformed_manifest():
    fails = []
    fx = Fixture()
    listed = ["claude-opus-5", "claude-sonnet-5"]
    try:
        # Baseline: the good manifest is fetched, validated and stored.
        settings = Settings()
        url = fx.write("good.json", MANIFEST)
        fetched = ml.fetch_manifest(url=url)
        if not ml.store_manifest(settings, fetched):
            fails.append("a valid manifest was refused")
        if ml.recommended_model("anthropic", ml.load_manifest(settings)) != \
                "claude-opus-5":
            fails.append("the stored manifest did not survive the round trip")

        # MUTATIONS: every way one can be wrong. None may be applied, none may
        # replace the good one already stored, and none may raise.
        mutations = [
            ("wrong schema", {"schema": 2, "providers": {"anthropic": {}}}),
            ("no schema", {"providers": {"anthropic": {}}}),
            ("no providers", {"schema": 1}),
            ("providers is a list", {"schema": 1, "providers": []}),
            ("providers is a string", {"schema": 1, "providers": "anthropic"}),
            ("not an object", ["schema", 1]),
            ("null", None),
            ("a bare number", 7),
        ]
        for why, mutant in mutations:
            if ml.valid_manifest(mutant) is not None:
                fails.append(f"a manifest with {why} passed validation")
            kept = ml.store_manifest(settings, mutant)
            if kept:
                fails.append(f"a manifest with {why} was stored")
            if ml.recommended_model("anthropic", ml.load_manifest(settings)) != \
                    "claude-opus-5":
                fails.append(f"a manifest with {why} replaced the good one")
            try:
                choices = ml.model_choices("anthropic", listed, mutant)
            except Exception as exc:
                fails.append(f"model_choices raised on {why}: {exc!r}")
                continue
            if [c.id for c in choices] != listed:
                fails.append(f"a manifest with {why} still altered the list")

        # Junk types INSIDE a valid envelope are ignored field by field: the
        # manifest is still readable, the unusable advice simply does not apply.
        messy = {"schema": 1, "providers": {"anthropic": {
            "recommended": 42, "good_choices": "gpt-5.1",
            "deprecated": {"a": 1}}}}
        if ml.valid_manifest(messy) is None:
            fails.append("a manifest with one bad field was thrown away whole")
        choices = ml.model_choices("anthropic", listed, messy)
        if [c.id for c in choices] != listed:
            fails.append("junk advice altered the list")
        if any(c.recommended or c.deprecated for c in choices):
            fails.append("junk advice produced a marking")

        # A manifest that is not JSON at all: fetch_manifest raises, the caller
        # stamps the check and keeps what it had.
        bad_url = fx.write("bad.json", "<html>404</html>")
        try:
            ml.fetch_manifest(url=bad_url)
            fails.append("fetch_manifest accepted an HTML body")
        except Exception:
            pass
        # A corrupt STORED manifest reads back as no manifest, not as a crash.
        settings.setValue(ml.SETTING_MANIFEST, "{not json")
        if ml.load_manifest(settings) != {}:
            fails.append("a corrupt stored manifest was returned")
        if ml.model_choices("anthropic", listed,
                            ml.load_manifest(settings))[0].recommended:
            fails.append("a corrupt stored manifest still marked a model")
    finally:
        fx.close()
    return fails


# ------------------------------------------------------------- H. the throttle

def test_throttle():
    fails = []
    from datetime import datetime, timedelta, timezone
    from studio import updater

    now = datetime(2026, 7, 31, 12, 0, 0, tzinfo=timezone.utc)
    settings = Settings()
    had = os.environ.pop("XSLOPE_NO_UPDATE_CHECK", None)
    try:
        if not ml.manifest_due(settings, now=now):
            fails.append("the first manifest check was skipped")
        # Storing stamps the check, whatever the outcome — including a refusal.
        ml.store_manifest(settings, MANIFEST, now=now)
        if ml.manifest_due(settings, now=now):
            fails.append("a manifest fetched now is immediately due again")
        if ml.manifest_due(settings, now=now + timedelta(hours=23)):
            fails.append("the manifest was re-fetched within the day")
        if not ml.manifest_due(settings, now=now + timedelta(hours=24)):
            fails.append("the manifest was never re-fetched after a day")
        # A failed fetch must also stamp, or an offline machine retries forever.
        settings2 = Settings()
        ml.store_manifest(settings2, {"schema": 99}, now=now)
        if ml.manifest_due(settings2, now=now):
            fails.append("a failed manifest fetch did not stamp the check")

        # The offline switch is the updater's, and it covers this too.
        os.environ["XSLOPE_NO_UPDATE_CHECK"] = "1"
        if ml.manifest_due(Settings(), now=now):
            fails.append("XSLOPE_NO_UPDATE_CHECK did not suppress the manifest fetch")
        os.environ.pop("XSLOPE_NO_UPDATE_CHECK")

        # And the URL is redirectable without touching code.
        os.environ["XSLOPE_ASSISTANT_MODELS_URL"] = "file:///tmp/x.json"
        if ml.manifest_url() != "file:///tmp/x.json":
            fails.append("XSLOPE_ASSISTANT_MODELS_URL is not honoured")
        os.environ.pop("XSLOPE_ASSISTANT_MODELS_URL")
        if not ml.manifest_url().endswith("/releases/latest/download/"
                                          "assistant_models.json"):
            fails.append(f"the default manifest URL is {ml.manifest_url()!r} — it "
                         f"must ride latest.json's release-asset convention")
        if updater.MANIFEST_URL.rsplit("/", 1)[0] != \
                ml.MANIFEST_URL.rsplit("/", 1)[0]:
            fails.append("the two manifests are not published side by side")
    finally:
        os.environ.pop("XSLOPE_ASSISTANT_MODELS_URL", None)
        if had is not None:
            os.environ["XSLOPE_NO_UPDATE_CHECK"] = had
    return fails


# --------------------------------------------------------------- I. the dialog

def test_dialog():
    """The real settings dialog, with its background refresh turned off — so this
    leg cannot reach a provider even on a machine that has a key stored."""
    fails = []
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication, QComboBox, QPushButton
    from studio.ai.config import AssistantConfig
    from studio.ai.settings_dialog import AssistantSettingsDialog

    QApplication.instance() or QApplication([])
    tmp = tempfile.mkdtemp(prefix="xslope-ai-settings-")
    settings = QSettings(os.path.join(tmp, "ai.ini"), QSettings.IniFormat)
    config = AssistantConfig(settings)
    had = os.environ.get("XSLOPE_NO_UPDATE_CHECK")
    os.environ["XSLOPE_NO_UPDATE_CHECK"] = "1"       # belt and braces: no network
    try:
        # --- offline, never connected: the shipped list, and it is editable ---
        dlg = AssistantSettingsDialog(config, auto_refresh=False)
        if not dlg.model.isEditable():
            fails.append("the model box is not editable — an unlisted model could "
                         "never be used")
        if dlg.refresh_btn.findChild(QComboBox) is not None:
            fails.append("unexpected widget nesting")
        if not isinstance(dlg.refresh_btn, QPushButton) or \
                not dlg.refresh_btn.isEnabled():
            fails.append("there is no working Refresh button")
        shipped = PROVIDERS["anthropic"]["models"]
        shown = [dlg.model.itemData(i) for i in range(dlg.model.count())]
        for model in shipped:
            if model not in shown:
                fails.append(f"the shipped fallback list is missing {model!r}")
        if "shipped" not in dlg.models_note.text():
            fails.append(f"the note does not say where the list came from: "
                         f"{dlg.models_note.text()!r}")

        # --- free text: type an id nothing has ever listed, and save it -------
        typed = "claude-model-published-this-morning"
        dlg.model.setCurrentText(typed)
        if dlg.current_model() != typed:
            fails.append(f"a typed model id read back as {dlg.current_model()!r}")
        dlg._save()
        if settings.value("ai/model/anthropic") != typed:
            fails.append(f"the typed model was not saved: "
                         f"{settings.value('ai/model/anthropic')!r}")
        if config.model() != typed:
            fails.append("the assistant would not use the typed model")

        # It survives the dialog being closed and reopened.
        dlg = AssistantSettingsDialog(config, auto_refresh=False)
        if dlg.current_model() != typed:
            fails.append(f"the typed model did not survive a reopen: "
                         f"{dlg.current_model()!r}")

        # --- a live result arrives: the list is replaced and cached ----------
        ml.store_manifest(settings, MANIFEST)
        live = ["claude-opus-5", "claude-sonnet-5", "claude-opus-4-1"]
        dlg._on_models_listed({"provider": "anthropic", "models": live,
                               "error": "", "manifest": None})
        ids = [dlg.model.itemData(i) for i in range(dlg.model.count())]
        if ids[0] != "claude-opus-5":
            fails.append(f"the live list did not reach the box in manifest order: "
                         f"{ids}")
        if "recommended" not in dlg.model.itemText(0):
            fails.append(f"the recommended model is not marked in the dialog: "
                         f"{dlg.model.itemText(0)!r}")
        if ml.load_cached(settings, "anthropic") != live:
            fails.append("a live list was not cached for the next offline open")
        if dlg.current_model() != typed:
            fails.append("a live refresh discarded what the user had typed")

        # --- the decoration must never reach the saved value -----------------
        dlg.model.setCurrentIndex(0)
        if dlg.current_model() != "claude-opus-5":
            fails.append(f"selecting the marked row yields "
                         f"{dlg.current_model()!r} — the marking leaked into the id")
        dlg._save()
        if settings.value("ai/model/anthropic") != "claude-opus-5":
            fails.append(f"the saved model carries its decoration: "
                         f"{settings.value('ai/model/anthropic')!r}")

        # --- a result for another provider is ignored ------------------------
        before = [dlg.model.itemData(i) for i in range(dlg.model.count())]
        dlg._on_models_listed({"provider": "openai", "models": ["gpt-5.1"],
                               "error": "", "manifest": None})
        after = [dlg.model.itemData(i) for i in range(dlg.model.count())]
        if before != after:
            fails.append("a stale answer for another provider overwrote the list")

        # --- a fresh install selects the recommended model -------------------
        fresh = QSettings(os.path.join(tmp, "fresh.ini"), QSettings.IniFormat)
        ml.store_manifest(fresh, MANIFEST)
        fresh_config = AssistantConfig(fresh)
        if fresh_config.model() != "claude-opus-5":
            fails.append(f"a fresh install does not default to the recommended "
                         f"model: {fresh_config.model()!r}")
        dlg2 = AssistantSettingsDialog(fresh_config, auto_refresh=False)
        if dlg2.current_model() != "claude-opus-5":
            fails.append(f"a fresh dialog does not preselect the recommended "
                         f"model: {dlg2.current_model()!r}")
        dlg2.deleteLater()
        dlg.deleteLater()
    finally:
        if had is None:
            os.environ.pop("XSLOPE_NO_UPDATE_CHECK", None)
        else:
            os.environ["XSLOPE_NO_UPDATE_CHECK"] = had
        shutil.rmtree(tmp, ignore_errors=True)
    return fails


# -------------------------------------------------------- the background runner

def test_runner():
    """The QThread the dialog actually uses: both reads happen off the GUI thread
    and come back as ONE payload, and a failed manifest never fails the models."""
    fails = []
    import time
    from PySide6.QtCore import QCoreApplication
    from PySide6.QtWidgets import QApplication
    from studio.runners import AssistantModelsRunner

    QApplication.instance() or QApplication([])
    fx = Fixture()

    def _run(**kw):
        got = []
        runner = AssistantModelsRunner(**kw)
        runner.listed.connect(got.append)
        runner.start()
        t0 = time.time()
        while not got and time.time() - t0 < 20:
            QCoreApplication.processEvents()
            time.sleep(0.01)
        runner.wait()
        return got[0] if got else None

    try:
        models_url = fx.write("openai.json", {"data": [
            {"id": "gpt-5.1"}, {"id": "text-embedding-3-small"}]})
        manifest_url = fx.write("manifest.json", MANIFEST)

        payload = _run(provider="openai", api_key="sk-x", want_models=True,
                       want_manifest=True, models_url=models_url,
                       manifest_url=manifest_url)
        if payload is None:
            fails.append("the runner never emitted a result")
            return fails
        if payload.get("provider") != "openai":
            fails.append(f"the payload names {payload.get('provider')!r}")
        if payload.get("models") != ["gpt-5.1"]:
            fails.append(f"the runner returned {payload.get('models')!r}")
        if ml.valid_manifest(payload.get("manifest")) is None:
            fails.append("the runner did not bring back the manifest")

        # An unreachable manifest is silent — the models still arrive.
        payload = _run(provider="openai", api_key="sk-x", want_models=True,
                       want_manifest=True, models_url=models_url,
                       manifest_url=_file_url(os.path.join(fx.dir, "gone.json")))
        if payload.get("manifest") is not None:
            fails.append("a failed manifest fetch was not swallowed")
        if payload.get("models") != ["gpt-5.1"]:
            fails.append("a failed manifest fetch took the model list with it")
        if not payload.get("manifest_checked"):
            fails.append("a failed manifest fetch did not report being attempted")

        # And with the models switched off (no key yet) only the manifest runs.
        payload = _run(provider="openai", want_models=False, want_manifest=True,
                       manifest_url=manifest_url)
        if payload.get("models") is not None or payload.get("error"):
            fails.append("the runner enumerated models it was told to skip")
        if ml.valid_manifest(payload.get("manifest")) is None:
            fails.append("the manifest-only run brought nothing back")
    finally:
        fx.close()
    return fails


# ------------------------------------------------------- the shipped manifest

def test_shipped_manifest():
    """The curated file in packaging/ must be one this build can read — it is
    published as a release asset, so a mistake there is a mistake on every
    install."""
    fails = []
    path = os.path.join(_REPO, "packaging", "assistant_models.json")
    if not os.path.exists(path):
        return [f"missing {path}"]
    sys.path.insert(0, os.path.join(_REPO, "packaging"))
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "make_assistant_manifest",
            os.path.join(_REPO, "packaging", "make_assistant_manifest.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        with open(path, encoding="utf-8") as fh:
            manifest = json.load(fh)
        for problem in mod.validate(manifest):
            fails.append(f"packaging/assistant_models.json: {problem}")
        if ml.valid_manifest(manifest) is None:
            fails.append("the shipped manifest would be ignored by Studio")
        # Every provider it advises must be one Studio offers, and its
        # recommendation must be usable as a model id.
        for name, entry in manifest.get("providers", {}).items():
            if name not in PROVIDERS:
                fails.append(f"the manifest advises on unknown provider {name!r}")
            rec = entry.get("recommended", "")
            if rec and rec != rec.strip():
                fails.append(f"{name}: recommended id has stray whitespace")
            choices = ml.model_choices(name, [], manifest)
            if rec and (not choices or choices[0].id != rec):
                fails.append(f"{name}: the recommendation does not sort first")
    finally:
        sys.path.remove(os.path.join(_REPO, "packaging"))
    return fails


CHECKS = [("the list-models endpoint per provider", test_endpoints),
          ("the response parse per provider", test_parse),
          ("the conservative non-chat filter", test_filter),
          ("enumeration end to end", test_enumeration),
          ("live -> cache -> shipped, under mutation", test_fallback_chain),
          ("the recommendations overlay", test_manifest_overlay),
          ("a malformed manifest is ignored", test_malformed_manifest),
          ("the once-a-day manifest throttle", test_throttle),
          ("the shipped packaging/assistant_models.json", test_shipped_manifest),
          ("the background runner", test_runner),
          ("the settings dialog", test_dialog)]

#: The legs that need Qt; everything else runs on an engine-only install.
_QT_CHECKS = (test_runner, test_dialog)


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
        checks = CHECKS
    except Exception:
        print("assistant models: PySide6 not installed — dialog check skipped.")
        checks = [c for c in CHECKS if c[1] not in _QT_CHECKS]
    failures = []
    for name, fn in checks:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("Assistant model list:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll assistant model-list checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
