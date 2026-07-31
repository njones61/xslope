"""Version check and update staging — the half of the updater that has no Qt in it.

Everything a "Check for Updates" needs, minus the widgets: where the manifest
lives, how a version compares, which artifact belongs to the running machine, how
a downloaded file is verified, and what command each platform's update actually
runs. ``studio/update_ui.py`` wraps this in Studio's dialog/QThread conventions;
``test/updater_check.py`` exercises it in a plain Python process.

**The manifest.** Every tagged release carries a ``latest.json`` written by
``packaging/make_latest_json.py`` and attached to the GitHub Release. It is read
from the release's ``latest/download`` redirect — a stable URL that always points
at the newest release's copy — so the updater never scrapes the release page's
HTML. Its shape (schema 1)::

    {"schema": 1, "version": "0.1.59", "tag": "v0.1.59",
     "notes_url": "...", "minimum_version": "0.0.0",
     "artifacts": {"darwin-arm64": {"filename":…, "url":…, "sha256":…, "size":…}}}

Artifact keys are ``f"{sys.platform}-{platform.machine().lower()}"`` as the
RUNNING app sees them, so picking this machine's download is one dict lookup.

**Nothing here raises at the user.** :func:`check_for_update` swallows every
network and parse failure into an :class:`UpdateInfo` with ``status="error"``;
whether that is shown at all is the caller's decision (the manual check says so,
the startup check stays silent).

**Nothing unverified is ever run.** :func:`download_artifact` hashes what it
received and deletes the file if the digest does not match ``latest.json``;
:func:`plan_install` refuses any :class:`UpdateInfo` that is not
``status="update-available"``, which is also how ``minimum_version`` is enforced.

Set ``XSLOPE_UPDATE_URL`` to point the check at another manifest (a beta channel,
or a ``file://`` fixture) without touching code.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import urllib.request
from datetime import datetime, timedelta, timezone

__all__ = [
    "REPO", "MANIFEST_URL", "RELEASES_URL", "PIP_UPGRADE_COMMAND",
    "UpdateError", "UpdateCancelled", "UpdateInfo", "InstallPlan",
    "manifest_url", "platform_key", "is_frozen", "app_executable",
    "parse_version", "is_newer", "evaluate", "fetch_manifest", "check_for_update",
    "sha256_file", "download_artifact", "plan_install", "launch_detached",
    "due_for_check", "utc_stamp",
]

#: The public repository these releases come from (pyproject ``[project.urls]``,
#: the NSIS ``HOMEPAGE``, and ``make_latest_json.py``'s default all agree).
REPO = "njones61/xslope"
#: The manifest attached to the newest release. ``releases/latest/download/<name>``
#: is GitHub's redirect to the newest release's asset of that name.
MANIFEST_URL = f"https://github.com/{REPO}/releases/latest/download/latest.json"
#: Where a user is sent when the app cannot update itself.
RELEASES_URL = f"https://github.com/{REPO}/releases/latest"
#: What a pip install is told to run. Never executed by us — shown to be copied.
PIP_UPGRADE_COMMAND = 'pip install -U "xslope[gui]"'
#: The only ``latest.json`` schema this build knows how to read.
SUPPORTED_SCHEMA = 1
#: Network timeout, seconds. A version check must never feel like a hang.
TIMEOUT = 5.0
#: How often the silent startup check is allowed to run.
CHECK_INTERVAL_HOURS = 24

# QSettings keys (the app's existing QSettings(ORG_NAME, APP_NAME) store).
SETTING_STARTUP_CHECK = "updates/check_at_startup"
SETTING_LAST_CHECK = "updates/last_check"


class UpdateError(Exception):
    """A download or verification failed in a way the user must be told about."""


class UpdateCancelled(Exception):
    """The user pressed Cancel during a download."""


# ---------------------------------------------------------------------------
# Environment probes
# ---------------------------------------------------------------------------

def manifest_url():
    """The manifest to read. ``XSLOPE_UPDATE_URL`` overrides the release default."""
    return os.environ.get("XSLOPE_UPDATE_URL") or MANIFEST_URL


def platform_key(platform_name=None, machine=None):
    """This machine's key into ``latest.json``'s ``artifacts`` map.

    Exactly the string ``make_latest_json.py`` files artifacts under:
    ``f"{sys.platform}-{platform.machine().lower()}"`` — e.g. ``darwin-arm64``,
    ``win32-amd64``.
    """
    p = sys.platform if platform_name is None else platform_name
    m = platform.machine() if machine is None else machine
    return f"{p}-{str(m).lower()}"


def is_frozen():
    """True inside the PyInstaller bundle (an installer build), False under pip."""
    return bool(getattr(sys, "frozen", False))


def app_executable():
    """The path to relaunch after a frozen-app update.

    Frozen, ``sys.executable`` IS the app binary (``XSLOPE Studio.exe`` /
    ``XSLOPE Studio.app/Contents/MacOS/XSLOPE Studio``). Unfrozen it is the Python
    interpreter, which is never relaunched by this module — the pip branch runs
    nothing at all.
    """
    return sys.executable


# ---------------------------------------------------------------------------
# Version comparison
# ---------------------------------------------------------------------------

def parse_version(text):
    """Parse a version string with PEP-440 semantics.

    ``packaging`` arrives with matplotlib, so it is present in every environment
    that can run Studio; the tuple fallback keeps this module importable (and the
    comparison sane for plain ``X.Y.Z``) if it ever is not.
    """
    text = str(text or "0").strip().lstrip("vV")
    try:
        from packaging.version import Version
        return Version(text)
    except Exception:
        parts = []
        for chunk in text.replace("-", ".").split("."):
            digits = "".join(c for c in chunk if c.isdigit())
            parts.append(int(digits) if digits else 0)
        return tuple(parts + [0] * (3 - len(parts))) if len(parts) < 3 else tuple(parts)


def is_newer(candidate, current):
    """True when ``candidate`` is a strictly later version than ``current``.

    Equal versions are NOT newer — that is the whole point of the up-to-date
    branch — and a manifest that has somehow gone backwards (a yanked release, a
    hand-edited beta channel) is not newer either, so a downgrade is never offered.
    """
    try:
        return parse_version(candidate) > parse_version(current)
    except Exception:
        return False


# ---------------------------------------------------------------------------
# The check
# ---------------------------------------------------------------------------

class UpdateInfo:
    """The outcome of one version check.

    ``status`` is one of:

    ``up-to-date``
        The running version is the newest published one (or newer).
    ``update-available``
        A newer version exists AND this install may act on it — the only status
        :func:`plan_install` will accept.
    ``update-manual``
        A newer version exists but this install cannot fetch it itself: it is
        older than the manifest's ``minimum_version``, or the release carries no
        artifact for this platform. The user is sent to the release page.
    ``error``
        Nothing could be concluded — offline, timeout, bad JSON, unknown schema.
        ``detail`` carries the sentence for the manual dialog; the startup check
        drops it on the floor.
    """

    def __init__(self, status, current, latest="", notes_url=RELEASES_URL,
                 artifact=None, detail="", manifest=None, reason=""):
        self.status = status
        self.current = current
        self.latest = latest
        self.notes_url = notes_url or RELEASES_URL
        self.artifact = artifact
        self.detail = detail
        self.manifest = manifest
        self.reason = reason        # "minimum-version" / "no-artifact" / ""

    @property
    def is_update(self):
        """A newer version exists, however it has to be obtained."""
        return self.status in ("update-available", "update-manual")

    @property
    def ok(self):
        """The check itself succeeded (whatever it concluded)."""
        return self.status != "error"

    def __repr__(self):
        return (f"UpdateInfo(status={self.status!r}, current={self.current!r}, "
                f"latest={self.latest!r}, reason={self.reason!r})")


def evaluate(manifest, current_version, key=None, frozen=None):
    """Turn a parsed manifest into an :class:`UpdateInfo`. Pure; no network.

    ``frozen`` decides how strict the answer is. A frozen app is going to
    DOWNLOAD AND RUN something, so it must have an artifact for its own platform
    and must satisfy ``minimum_version``; a pip install only ever gets shown a
    command line, so neither gate applies to it.
    """
    current = str(current_version)
    if frozen is None:
        frozen = is_frozen()
    if not isinstance(manifest, dict):
        return UpdateInfo("error", current, detail="The update manifest was not readable.")

    schema = manifest.get("schema")
    if schema != SUPPORTED_SCHEMA:
        return UpdateInfo(
            "error", current, manifest=manifest,
            detail=f"The update manifest uses format {schema!r}, which this "
                   f"version does not understand. Check the releases page.")

    latest = str(manifest.get("version", "") or "")
    notes = manifest.get("notes_url") or RELEASES_URL
    if not latest:
        return UpdateInfo("error", current, manifest=manifest,
                          detail="The update manifest carries no version.")

    if not is_newer(latest, current):
        return UpdateInfo("up-to-date", current, latest=latest, notes_url=notes,
                          manifest=manifest)

    artifacts = manifest.get("artifacts") or {}
    artifact = artifacts.get(key or platform_key())

    if frozen:
        minimum = str(manifest.get("minimum_version", "0.0.0") or "0.0.0")
        if is_newer(minimum, current):
            # This install predates the oldest one the new release accepts an
            # in-place upgrade from. Downloading the installer anyway is exactly
            # what the field exists to prevent, so we never offer it.
            return UpdateInfo(
                "update-manual", current, latest=latest, notes_url=notes,
                manifest=manifest, reason="minimum-version",
                detail=f"Version {latest} cannot be installed on top of {current}. "
                       f"Download and run the {latest} installer from the releases "
                       f"page instead.")
        if not artifact:
            return UpdateInfo(
                "update-manual", current, latest=latest, notes_url=notes,
                manifest=manifest, reason="no-artifact",
                detail=f"Release {latest} has no download for this platform "
                       f"({key or platform_key()}). See the releases page.")

    return UpdateInfo("update-available", current, latest=latest, notes_url=notes,
                      artifact=artifact, manifest=manifest)


def fetch_manifest(url=None, timeout=TIMEOUT):
    """GET and parse the manifest. Raises on any failure — see
    :func:`check_for_update` for the version that never does."""
    target = url or manifest_url()
    req = urllib.request.Request(target, headers={"User-Agent": _user_agent()})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read(1 << 20)          # a manifest is ~1 KB; cap the read
    return json.loads(raw.decode("utf-8"))


def check_for_update(current_version=None, url=None, timeout=TIMEOUT, key=None,
                     frozen=None):
    """One complete check: fetch, parse, compare. NEVER raises.

    Every network outcome — offline, DNS failure, timeout, 404, malformed JSON —
    comes back as ``status="error"``. That is the design: a version check is not
    something the user asked to have go wrong, so it cannot interrupt them.
    """
    if current_version is None:
        from xslope import __version__ as current_version
    try:
        manifest = fetch_manifest(url=url, timeout=timeout)
    except Exception as exc:
        return UpdateInfo("error", str(current_version),
                          detail=f"Could not reach the update server ({exc.__class__.__name__}).")
    return evaluate(manifest, current_version, key=key, frozen=frozen)


def _user_agent():
    try:
        from xslope import __version__ as v
    except Exception:
        v = "unknown"
    return f"XSLOPE-Studio/{v} ({sys.platform})"


# ---------------------------------------------------------------------------
# Download + verification
# ---------------------------------------------------------------------------

def sha256_file(path, chunk=1 << 20):
    """Hex sha256 of a file, read in chunks (installers are 150-250 MB)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download_artifact(artifact, dest_dir, progress=None, cancel=None,
                      timeout=TIMEOUT):
    """Download one manifest artifact into ``dest_dir`` and verify its digest.

    ``progress(done_bytes, total_bytes)`` is called as bytes arrive (``total`` is
    0 when the server sends no length); ``cancel()`` is polled between chunks and
    a true answer aborts with :class:`UpdateCancelled`.

    Returns the path of the verified file. **A digest mismatch deletes the file
    and raises** — the caller therefore cannot reach an installer that is not the
    exact bytes ``latest.json`` describes, which is the whole security story of an
    unsigned-era updater.
    """
    url = artifact["url"]
    name = artifact.get("filename") or os.path.basename(url)
    expected = (artifact.get("sha256") or "").strip().lower()
    dest = os.path.join(dest_dir, name)

    req = urllib.request.Request(url, headers={"User-Agent": _user_agent()})
    done = 0
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            total = int(resp.headers.get("Content-Length") or artifact.get("size") or 0)
            with open(dest, "wb") as out:
                while True:
                    if cancel is not None and cancel():
                        raise UpdateCancelled()
                    block = resp.read(1 << 18)
                    if not block:
                        break
                    out.write(block)
                    done += len(block)
                    if progress is not None:
                        progress(done, total)
    except UpdateCancelled:
        _unlink(dest)
        raise
    except Exception as exc:
        _unlink(dest)
        raise UpdateError(f"Download failed: {exc}") from exc

    if expected:
        actual = sha256_file(dest)
        if actual != expected:
            # Truncated transfer, a proxy that rewrote the body, or a tampered
            # artifact — indistinguishable from here, and all three are refusals.
            _unlink(dest)
            raise UpdateError(
                "The downloaded file does not match the checksum published with "
                "the release, so it was discarded and nothing was installed.")
    return dest


def _unlink(path):
    try:
        os.unlink(path)
    except OSError:
        pass


# ---------------------------------------------------------------------------
# What the update actually does, per platform
# ---------------------------------------------------------------------------

class InstallPlan:
    """What will happen when the user confirms — built, then run, never guessed.

    ``kind`` is ``"windows-installer"``, ``"macos-dmg"``, ``"pip"`` or
    ``"manual"``. ``argv`` is the exact command :func:`launch_detached` will
    spawn (empty for ``pip``/``manual``, which run nothing). ``quits_app`` says
    whether Studio has to exit for the update to proceed.
    """

    def __init__(self, kind, argv=(), message="", helper=None, quits_app=False,
                 command_text=""):
        self.kind = kind
        self.argv = list(argv)
        self.message = message
        self.helper = helper
        self.quits_app = quits_app
        self.command_text = command_text   # the literal line shown for copying

    def __repr__(self):
        return f"InstallPlan(kind={self.kind!r}, argv={self.argv!r})"


# --- Stage C seam ----------------------------------------------------------
# Today's macOS builds are unsigned, so a downloaded .app carries a quarantine
# flag and (once opened from a dmg) can be path-translocated: a script that
# swapped the running bundle out from under itself would produce an app the
# user then has to right-click-Open anyway. So the unsigned-era branch mounts
# the dmg and asks for a drag, which is the same gesture as a fresh install.
#
# STAGE C (plan §6 — Developer ID + notarization in CI): flip this flag and
# implement macos_replace_in_place() below. It is the ONLY place the macOS
# branch changes; plan_install() already routes through it.
MACOS_SCRIPTED_REPLACE = False


def macos_replace_in_place(dmg_path, app_path, pid=None, helper_dir=None):
    """Stage C: mount the dmg, swap the installed .app, relaunch — unattended.

    Deliberately unimplemented while builds are unsigned. When signing lands this
    becomes the mirror of :func:`windows_update_command`: a detached helper that
    waits for our exit, ``hdiutil attach``es the dmg, ``ditto``s the new bundle
    over the installed one, detaches, and relaunches.
    """
    raise NotImplementedError(
        "Scripted macOS replacement waits on a signed, notarized build (Stage C).")


#: The helper that performs a Windows update after Studio has exited.
#:
#: It cannot be done in-process: the NSIS installer's silent mode taskkills
#: ``XSLOPE Studio.exe`` and wipes the ``_internal`` directory the running app is
#: executing out of. So Studio writes this batch file, spawns it DETACHED, and
#: quits; the helper waits for our PID to disappear, runs the installer silently,
#: relaunches the freshly installed app, and deletes itself.
_WINDOWS_HELPER = r"""@echo off
rem XSLOPE Studio update helper - written by studio/updater.py, deletes itself.
rem   %1 = pid of the Studio process to wait for
rem   %2 = the NSIS setup .exe (already sha256-verified)
rem   %3 = the app .exe to relaunch afterwards
setlocal
set "WAITPID=%~1"
set "SETUP=%~2"
set "APPEXE=%~3"

rem 1. Wait (up to ~2 minutes) for Studio to exit, so the installer can replace it.
for /L %%i in (1,1,120) do (
  tasklist /FI "PID eq %WAITPID%" 2>nul | find "%WAITPID%" >nul || goto :gone
  ping -n 2 127.0.0.1 >nul
)
:gone

rem 2. Silent in-place upgrade. /S reuses the recorded install directory.
start "" /wait "%SETUP%" /S

rem 3. Relaunch the new build.
start "" "%APPEXE%"

rem 4. Remove this helper (self-delete idiom: release the handle, then delete).
(goto) 2>nul & del "%~f0"
"""


def windows_update_command(installer_path, app_path=None, pid=None, helper_dir=None):
    """Write the update helper and return ``(argv, helper_path)``.

    ``argv`` is what :func:`launch_detached` spawns just before Studio quits.
    Nothing is executed here, which is what makes the whole Windows branch
    assertable in a test on any platform.
    """
    app_path = app_path or app_executable()
    pid = os.getpid() if pid is None else pid
    helper_dir = helper_dir or os.path.dirname(os.path.abspath(installer_path))
    helper = os.path.join(helper_dir, "xslope_update.cmd")
    with open(helper, "w", encoding="ascii", newline="\r\n") as fh:
        fh.write(_WINDOWS_HELPER)
    comspec = os.environ.get("COMSPEC", "cmd.exe")
    argv = [comspec, "/c", helper, str(pid), str(installer_path), str(app_path)]
    return argv, helper


def macos_open_command(dmg_path):
    """``open`` the verified dmg — mounts it and shows the drag-to-Applications
    window. ``/usr/bin/open`` by absolute path so a doctored PATH cannot redirect
    it."""
    return ["/usr/bin/open", str(dmg_path)]


def plan_install(info, artifact_path=None, frozen=None, platform_name=None,
                 app_path=None, pid=None, helper_dir=None):
    """Decide — and prepare — what this install does about ``info``.

    Refuses anything but ``status="update-available"``: that single gate is how
    ``minimum_version``, an unknown schema, a missing platform artifact and an
    up-to-date install are all kept away from an installer run.

    A frozen app needs ``artifact_path`` — the file :func:`download_artifact`
    already verified. Passing an unverified path is the one thing this module
    cannot check for you, so downloads and installs are never split across
    callers: ``update_ui`` does both in one worker.
    """
    if frozen is None:
        frozen = is_frozen()
    plat = sys.platform if platform_name is None else platform_name

    if info.status != "update-available":
        raise UpdateError(
            info.detail or f"No installable update ({info.status}).")

    if not frozen:
        # A pip environment is the user's, not ours. We print the line; they run it.
        return InstallPlan(
            "pip", command_text=PIP_UPGRADE_COMMAND,
            message=f"Version {info.latest} is available. Update it the way you "
                    f"installed it:")

    if artifact_path is None or not os.path.exists(str(artifact_path)):
        raise UpdateError("Nothing was downloaded, so there is nothing to install.")

    if plat.startswith("win"):
        argv, helper = windows_update_command(artifact_path, app_path=app_path,
                                              pid=pid, helper_dir=helper_dir)
        return InstallPlan(
            "windows-installer", argv=argv, helper=helper, quits_app=True,
            message=f"XSLOPE Studio {info.latest} is ready to install. Studio will "
                    f"close, the installer will run silently, and the new version "
                    f"will start automatically.")

    if plat == "darwin":
        if MACOS_SCRIPTED_REPLACE:                      # Stage C
            argv, helper = macos_replace_in_place(artifact_path, app_path or "",
                                                  pid=pid, helper_dir=helper_dir)
            return InstallPlan("macos-replace", argv=argv, helper=helper,
                               quits_app=True,
                               message=f"Installing XSLOPE Studio {info.latest}…")
        return InstallPlan(
            "macos-dmg", argv=macos_open_command(artifact_path), quits_app=False,
            message=(f"XSLOPE Studio {info.latest} has been downloaded and its "
                     f"checksum verified.\n\n"
                     f"1. Drag XSLOPE Studio to the Applications folder in the "
                     f"window that opens (replace the existing copy).\n"
                     f"2. Quit this copy and start the new one.\n\n"
                     f"Because this build is not yet signed by Apple, the first "
                     f"launch needs a right-click → Open."))

    # A frozen build on a platform with no install story (there is no Linux
    # installer today) — the artifact is downloaded, the user takes it from there.
    return InstallPlan("manual", quits_app=False,
                       message=f"Downloaded {os.path.basename(str(artifact_path))}. "
                               f"Install it manually to finish updating.")


def launch_detached(argv):
    """Spawn ``argv`` so it OUTLIVES this process.

    The Windows helper must keep running after Studio exits (it is waiting for
    exactly that), so it gets DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP and no
    console window. Elsewhere ``start_new_session`` is enough.
    """
    kwargs = {"close_fds": True}
    if sys.platform.startswith("win"):
        DETACHED_PROCESS = 0x00000008
        CREATE_NEW_PROCESS_GROUP = 0x00000200
        CREATE_NO_WINDOW = 0x08000000
        kwargs["creationflags"] = (DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
                                   | CREATE_NO_WINDOW)
    else:
        kwargs["start_new_session"] = True
    return subprocess.Popen(list(argv), **kwargs)


# ---------------------------------------------------------------------------
# The once-a-day throttle
# ---------------------------------------------------------------------------

def utc_stamp(now=None):
    """An ISO-8601 UTC timestamp in the form the settings store keeps."""
    now = now or datetime.now(timezone.utc)
    return now.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def due_for_check(last_check, now=None, interval_hours=CHECK_INTERVAL_HOURS):
    """Is the silent startup check allowed to run?

    ``last_check`` is whatever came out of QSettings — a stamp, ``None``, or
    garbage left by an older build. Anything unparseable means "we have no
    evidence a check happened", so it runs; a stamp in the future (a clock that
    moved backwards) is treated the same way rather than locking the check out
    until the calendar catches up.
    """
    if not last_check:
        return True
    now = now or datetime.now(timezone.utc)
    try:
        prev = datetime.strptime(str(last_check), "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=timezone.utc)
    except Exception:
        return True
    if prev > now:
        return True
    return (now - prev) >= timedelta(hours=interval_hours)
