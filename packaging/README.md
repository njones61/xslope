# Native installers for XSLOPE Studio

XSLOPE ships through two channels. `pip install "xslope[gui]"` serves anyone who
already has Python; the installers built here serve everyone who does not — a
double-clickable **XSLOPE Studio** with the whole scientific stack, gmsh, the AI
assistant and the compiled fast kernel inside it.

| File | What it is |
| --- | --- |
| `entry_studio.py` | The frozen app's entry point: launches Studio, or answers `--version` / `--self-test`. |
| `studio.spec` | PyInstaller spec — what goes in the bundle and what stays out. |
| `build_app.py` | One command: verify → compile kernel → freeze → smoke-test → `.dmg`. |
| `make_latest_json.py` | Writes the `latest.json` update manifest from the built artifacts. |
| `assistant_models.json` | Curated model recommendations for the AI assistant — hand-edited, published as a release asset. |
| `make_assistant_manifest.py` | Validates `assistant_models.json` and stages it beside `latest.json`. |
| `windows/installer.nsi` | NSIS script: per-user install, Start-menu shortcut, uninstaller, silent `/S` upgrade. |
| `macos/entitlements.plist` | Hardened-runtime entitlements used when a signed build is notarized. |

Build outputs go to `packaging/dist` (artifacts) and `packaging/build` (scratch).
Both are git-ignored; no build output is ever committed.

## Building locally

```bash
python -m venv .venv && source .venv/bin/activate
pip install ".[gui,ai,fem,cad]" pyinstaller Cython pillow
python packaging/build_app.py            # macOS: also assembles the .dmg
python packaging/build_app.py --clean    # remove every build output afterwards
```

The build refuses to finish unless the frozen app passes its own smoke test: the
shipped binary is launched offscreen (`QT_QPA_PLATFORM=offscreen`) with
`--self-test`, and must load a bundled sample project, find the packaged input
template and skill prompt, build a **tri6 mesh with the bundled gmsh**, run a
Bishop solve, import the compiled fast kernel, and open the real Studio main
window. Every step prints an `[ok]` line, and `build_app.py` asserts on those
lines rather than on the exit code alone.

On Windows a windowed build has no console, so the app also writes that
transcript to the file named by `XSLOPE_SELFTEST_LOG`, which is what the build
script reads back.

## What the bundle contains

* **gmsh** — the pip package is a single `gmsh.py` plus an 80 MB shared library
  installed *outside* site-packages (`<prefix>/lib/libgmsh.4.15.dylib`,
  `<prefix>/Lib/gmsh-4.15.dll`). PyInstaller's module graph never sees the
  library, so the spec asks the imported `gmsh` module where its library
  actually is (`gmsh.lib._name`, the path ctypes resolved) and collects it to
  the bundle root — which is exactly where `gmsh.py` looks first, since for a
  frozen module its own directory is `sys._MEIPASS`.
* **The compiled fast kernel** — `build_app.py` runs
  `setup_kernel.py build_ext --inplace` before freezing, and the spec then names
  `xslope._fem_kernel` as a hidden import so PyInstaller places the extension
  module inside the frozen `xslope` package. PyPI wheels stay pure-Python; only
  the installers carry the compiled kernel.
* **`xslope/resources/`** — the packaged input template and the AI skill prompt,
  kept at their normal package-relative path so `default_template_path()` and
  `importlib.resources.files("xslope")` work unchanged when frozen.
* **`samples/`** — every sample project in `docs/inputs/slope` plus the blank
  master template, so a first-time user has something to open.
* **Trimmed** — no tkinter/Tcl, one Qt binding only, Qt modules Studio never
  touches (WebEngine, Quick/QML, 3D, Charts, Multimedia, …) excluded, matplotlib
  restricted to the Qt and Agg backends, and no notebook/test/docs stacks.

Onedir, not onefile: onefile would unpack ~400 MB to a temp directory on every
launch, and the updater wants a stable install tree it can patch in place.

## How versioning flows

`xslope/_version.py` is the single source of truth. Everything downstream reads
it — nothing is typed twice:

```
xslope/_version.py  __version__ = "0.1.58"
   |
   +-- pyproject.toml  [tool.setuptools.dynamic] version = { attr = ... }   -> the wheel
   |
   +-- studio.spec     VERSION -> BUNDLE(version=...) and Info.plist
   |                     CFBundleShortVersionString / CFBundleVersion
   |
   +-- build_app.py    -> the .dmg filename  XSLOPE-Studio-0.1.58-macos-arm64.dmg
   |                   -> the Windows VSVersionInfo resource on XSLOPE Studio.exe
   |
   +-- CI              -> makensis /DVERSION=0.1.58, which stamps
   |                        VIProductVersion, the Add/Remove Programs
   |                        DisplayVersion, and HKCU\Software\XSLOPEStudio\Version
   |
   +-- latest.json     "version": "0.1.58"
```

Three cheap version probes exist, in increasing cost:

1. **Metadata, no launch** — `Info.plist` → `CFBundleShortVersionString` on
   macOS; the `.exe` file-version resource, or
   `HKCU\Software\XSLOPEStudio\Version`, on Windows.
2. **`XSLOPE Studio --version`** — prints the version and exits before Qt starts.
3. `xslope.__version__` from inside the running app.

A release tag `vX.Y.Z` must match `__version__`; CI derives the manifest version
by stripping the leading `v`.

## `latest.json` — the update manifest

Attached to every tagged release. The in-app updater (Stage B) fetches it from
the release's `latest` redirect and never scrapes HTML. Platform keys are
`f"{sys.platform}-{platform.machine().lower()}"` as the *running app* sees them,
so the lookup is one line.

```json
{
  "schema": 1,
  "version": "0.1.58",
  "tag": "v0.1.58",
  "released": "2026-07-31T19:04:11Z",
  "notes_url": "https://github.com/njones61/xslope/releases/tag/v0.1.58",
  "minimum_version": "0.0.0",
  "artifacts": {
    "darwin-arm64": {
      "filename": "XSLOPE-Studio-0.1.58-macos-arm64.dmg",
      "url": "https://github.com/njones61/xslope/releases/download/v0.1.58/XSLOPE-Studio-0.1.58-macos-arm64.dmg",
      "sha256": "1f0c…",
      "size": 177900000
    },
    "win32-amd64": {
      "filename": "XSLOPE-Studio-0.1.58-windows-x64-setup.exe",
      "url": "https://github.com/njones61/xslope/releases/download/v0.1.58/XSLOPE-Studio-0.1.58-windows-x64-setup.exe",
      "sha256": "9ab3…",
      "size": 210300000
    }
  }
}
```

Fields: `schema` is bumped only on an incompatible change; `minimum_version` is
the oldest install that may update straight to this one (an older install is
told to download a fresh installer instead); `sha256` must be verified before
anything downloaded is run. A universal2 dmg would simply be listed under both
`darwin-arm64` and `darwin-x86_64`.

## `assistant_models.json` — the model recommendations

The second release asset, read from the same `releases/latest/download/`
redirect and at most once a day. It carries no build output — it is the
hand-curated `packaging/assistant_models.json`, which
`make_assistant_manifest.py` validates and stages — so the assistant's model
advice can be corrected between releases.

```json
{
  "schema": 1,
  "updated": "2026-07-31",
  "providers": {
    "anthropic": {
      "recommended": "claude-opus-5",
      "good_choices": [
        {"id": "claude-opus-5", "label": "flagship", "note": "…"}
      ],
      "deprecated": ["claude-opus-4-1"]
    }
  }
}
```

The manifest only marks and orders what a provider itself lists: `recommended`
is what a fresh install selects and carries a "recommended" marking,
`good_choices` sort to the top with their labels, and `deprecated` ids stay
selectable but are marked as superseded. A provider with no entry (Ollama, whose
list is whatever the user has pulled) is shown unchanged, and a manifest this
build cannot read is ignored — the dialog then shows exactly the list it would
have shown without one. Edit the file, tag a release, and every install picks the
new advice up within a day.

## CI

`.github/workflows/release-installers.yml` builds macOS arm64 and Windows x64
(an Intel-Mac leg is present but commented out). `workflow_dispatch` produces
test artifacts; pushing a `v*` tag additionally creates the GitHub Release and
attaches the installers, their `.sha256` files, `latest.json` and
`assistant_models.json`.

**Signing is optional and gated on secrets.** Each signing step runs only when
its secret is present, so the workflow is green — and the artifacts work — with
no certificate at all. Today's builds are unsigned: macOS users right-click →
Open on first launch, Windows users click through SmartScreen. Adding the
secrets below turns signing on with no change to the workflow.

| Secret | Used for |
| --- | --- |
| `MACOS_CERT_P12`, `MACOS_CERT_PASSWORD`, `MACOS_SIGN_IDENTITY` | `codesign` the `.app` with a Developer ID |
| `APPLE_ID`, `APPLE_TEAM_ID`, `APPLE_APP_PASSWORD` | `notarytool` submit + `stapler staple` the `.dmg` |
| `WINDOWS_CERT_PFX`, `WINDOWS_CERT_PASSWORD` | `signtool` on the app `.exe` and the installer |

## Installing, per platform

* **macOS** — open the `.dmg`, drag **XSLOPE Studio** to Applications. Until the
  build is signed, first launch needs right-click → Open.
* **Windows** — run the setup `.exe`. It installs per user under
  `%LOCALAPPDATA%\Programs\XSLOPE Studio` (no administrator rights), adds a
  Start-menu shortcut and an Add/Remove Programs entry, and supports
  `setup.exe /S` for a silent in-place upgrade — the path the updater uses.
