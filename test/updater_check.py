"""Standing checks on Studio's in-app updater — the decisions that stand between
a version number on a web server and an installer running on the user's machine.

Nothing here touches the network. The manifest and the "installer" are files in a
temporary directory served over ``file://``, which is the same code path
``urllib`` takes for the real release URL, so the fetch, the parse, the download
loop and the digest check are all genuinely exercised — only the host is local.
No installer is ever run: the platform branches are asserted on the COMMAND they
would spawn.

What is guarded:

  A. THE COMPARISON — a newer manifest version is an update, an equal one is not,
     and an OLDER one is not either (a yanked release must never offer a
     downgrade). ``v``-prefixed and pre-release strings compare by PEP 440, not
     by string order, which is the difference between 0.1.9 and 0.1.10.
  B. THE PLATFORM KEY — the key is ``f"{sys.platform}-{platform.machine().lower()}"``,
     byte-identical to the one ``packaging/make_latest_json.py`` files artifacts
     under, and a release with no artifact for the running machine does not
     become an offer a frozen app can act on.
  C. MINIMUM_VERSION (mutation) — the same manifest, with ``minimum_version``
     raised above the running version, must stop being installable: the status
     changes and ``plan_install`` REFUSES. This is the field's only job, and it
     is invisible unless something asserts on it.
  D. THE CHECKSUM (mutation) — a download whose sha256 matches is returned; the
     same bytes against a corrupted digest raise, LEAVE NO FILE BEHIND, and
     never reach an install plan. An updater that runs an unverified binary is
     the one bug in this feature that matters.
  E. FROZEN vs PIP — the branch is taken on ``sys.frozen`` alone. A frozen app
     gets its platform's install command; a pip install gets the literal
     ``pip install -U "xslope[gui]"`` line and an argv that spawns NOTHING, because
     self-modifying someone's Python environment is not this app's business.
  F. THE THROTTLE — the startup check runs at most once a day, and an absent,
     malformed or future-dated stamp all mean "check now" rather than "never
     check again".
  G. THE MENU AND THE TWO PATHS — Help carries the item and the preference; a
     manual check always ends in exactly ONE dialog (up to date / new version /
     could-not-check) and a startup check ends in one only for a new version.
     A dead server produces a dialog on the manual path and silence on the
     automatic one.
  H. END TO END — check → new version → download a (small, fake) artifact →
     verify → the Windows branch reaches its hand-off with the exact argv and
     helper script it would run; the macOS branch, for an installed bundle it can
     replace, reaches its own helper and quits the app, and for one it cannot (not
     an app, translocated) falls back to ``open <dmg>`` with drag instructions
     that say nothing about signing.
  I. THE DOWNLOAD, THROUGH THE APP — the same happy path driven by the real
     controller and its background runner (progress bar up, plan produced, bar
     down), and the same run with a corrupted digest, which must land on the
     failure dialog and produce NO install plan at all.

  J. THE macOS SWAP, RUN — on a Mac, the helper is run for real against a small
     dmg and a stand-in installed bundle: the bundle is replaced, nothing is left
     beside it, the dmg and the helper are deleted. Two mutations: a dmg with no
     app in it, and a folder the swap cannot write to, must both leave the
     installed bundle exactly as it was and report failure.

Skips its Studio legs cleanly (exit 0) when PySide6 is not installed, and J off
macOS; A-F run either way.
"""
import hashlib
import json
import os
import shutil
import sys
import tempfile
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from studio import updater

CURRENT = "0.1.58"
NEWER = "0.1.59"


# --------------------------------------------------------------------------
# Fixtures: a manifest and an "installer", on disk, addressed as file:// URLs
# --------------------------------------------------------------------------

def _file_url(path):
    from urllib.request import pathname2url
    return "file://" + pathname2url(os.path.abspath(path))


class Fixture:
    """A temp directory holding a fake artifact and a latest.json that describes
    it — the same JSON shape ``make_latest_json.py`` writes, with real sha256s."""

    def __init__(self, version=NEWER, minimum="0.0.0", keys=None, schema=1,
                 payload=b"not really an installer, but it hashes like one\n"):
        self.dir = tempfile.mkdtemp(prefix="xslope-updater-test-")
        self.payload = payload
        self.name = f"XSLOPE-Studio-{version}-windows-x64-setup.exe"
        self.artifact_path = os.path.join(self.dir, self.name)
        with open(self.artifact_path, "wb") as fh:
            fh.write(payload)
        self.sha = hashlib.sha256(payload).hexdigest()
        entry = {"filename": self.name, "url": _file_url(self.artifact_path),
                 "sha256": self.sha, "size": len(payload)}
        if keys is None:
            keys = [updater.platform_key(), "darwin-arm64", "win32-amd64"]
        self.manifest = {
            "schema": schema,
            "version": version,
            "tag": f"v{version}",
            "released": "2026-07-31T19:04:11Z",
            "notes_url": f"https://github.com/{updater.REPO}/releases/tag/v{version}",
            "minimum_version": minimum,
            "artifacts": {k: dict(entry) for k in keys},
        }
        self.manifest_path = os.path.join(self.dir, "latest.json")
        self.write()

    def write(self):
        with open(self.manifest_path, "w", encoding="utf-8") as fh:
            json.dump(self.manifest, fh, indent=2)

    @property
    def url(self):
        return _file_url(self.manifest_path)

    def close(self):
        shutil.rmtree(self.dir, ignore_errors=True)


# --------------------------------------------------------- A. the comparison

def test_version_comparison():
    fails = []
    cases = [
        ("0.1.59", "0.1.58", True,  "a later patch is newer"),
        ("0.1.58", "0.1.58", False, "the same version is not newer"),
        ("0.1.57", "0.1.58", False, "an older version is never offered"),
        ("0.2.0",  "0.1.99", True,  "minor beats patch"),
        ("0.1.10", "0.1.9",  True,  "0.1.10 > 0.1.9 (numeric, not string, order)"),
        ("v0.1.59", "0.1.58", True, "a v-prefixed tag compares as its version"),
        ("0.1.59",  "0.1.59rc1", True, "a release is newer than its own rc"),
        ("0.1.59rc1", "0.1.59", False, "an rc is not newer than the release"),
    ]
    for cand, cur, want, why in cases:
        got = updater.is_newer(cand, cur)
        if got != want:
            fails.append(f"is_newer({cand!r}, {cur!r}) = {got}, expected {want} — {why}")

    # The three statuses the comparison alone produces, through evaluate().
    fx = Fixture(version=CURRENT)
    info = updater.evaluate(fx.manifest, CURRENT, frozen=False)
    if info.status != "up-to-date":
        fails.append(f"same version gave {info.status!r}, expected 'up-to-date'")
    fx.manifest["version"] = NEWER
    info = updater.evaluate(fx.manifest, CURRENT, frozen=False)
    if info.status != "update-available" or info.latest != NEWER:
        fails.append(f"newer version gave {info!r}")
    if not info.is_update or not info.ok:
        fails.append("an available update must report is_update and ok")
    fx.manifest["version"] = "0.1.0"
    info = updater.evaluate(fx.manifest, CURRENT, frozen=False)
    if info.status != "up-to-date":
        fails.append(f"an OLDER manifest gave {info.status!r} — a downgrade was offered")

    # An unreadable manifest is an error, never an offer.
    for bad in (None, {}, {"schema": 99, "version": NEWER}, {"schema": 1}):
        info = updater.evaluate(bad, CURRENT, frozen=True)
        if info.status != "error":
            fails.append(f"manifest {bad!r} gave {info.status!r}, expected 'error'")
        if info.ok:
            fails.append(f"manifest {bad!r} reported ok")
    fx.close()
    return fails


# -------------------------------------------------------- B. the platform key

def test_platform_key():
    fails = []
    import platform as _p
    want = f"{sys.platform}-{_p.machine().lower()}"
    if updater.platform_key() != want:
        fails.append(f"platform_key() = {updater.platform_key()!r}, expected {want!r}")
    # The exact strings make_latest_json.py writes.
    if updater.platform_key("darwin", "arm64") != "darwin-arm64":
        fails.append("darwin/arm64 must key as 'darwin-arm64'")
    if updater.platform_key("win32", "AMD64") != "win32-amd64":
        fails.append("machine() is lower-cased: win32/AMD64 -> 'win32-amd64'")

    fx = Fixture(keys=["darwin-arm64", "win32-amd64"])
    # The running machine's entry is the one selected.
    info = updater.evaluate(fx.manifest, CURRENT, key="win32-amd64", frozen=True)
    if info.artifact is None or info.artifact["filename"] != fx.name:
        fails.append("the win32-amd64 key did not select its own artifact")
    info = updater.evaluate(fx.manifest, CURRENT, key="darwin-arm64", frozen=True)
    if info.artifact is None:
        fails.append("the darwin-arm64 key did not select its own artifact")

    # A release with nothing for this machine is not an offer a frozen app acts on.
    info = updater.evaluate(fx.manifest, CURRENT, key="linux-x86_64", frozen=True)
    if info.status != "update-manual" or info.reason != "no-artifact":
        fails.append(f"a missing platform artifact gave {info!r}, expected "
                     f"update-manual/no-artifact")
    try:
        updater.plan_install(info, artifact_path=fx.artifact_path, frozen=True,
                             platform_name="linux")
        fails.append("plan_install accepted an update-manual info")
    except updater.UpdateError:
        pass
    # A pip install has no artifact to miss, so the same manifest is still an offer.
    info = updater.evaluate(fx.manifest, CURRENT, key="linux-x86_64", frozen=False)
    if info.status != "update-available":
        fails.append(f"pip install with no artifact gave {info.status!r}")
    fx.close()
    return fails


# ------------------------------------------------ C. minimum_version (mutation)

def test_minimum_version():
    fails = []
    key = "win32-amd64"

    # Baseline: minimum_version below us -> installable.
    fx = Fixture(minimum="0.0.0")
    info = updater.evaluate(fx.manifest, CURRENT, key=key, frozen=True)
    if info.status != "update-available":
        fails.append(f"baseline manifest gave {info.status!r}, expected installable")
    plan = updater.plan_install(info, fx.artifact_path, frozen=True,
                                platform_name="win32", pid=4242,
                                app_path=r"C:\App\XSLOPE Studio.exe")
    if plan.kind != "windows-installer":
        fails.append(f"baseline plan was {plan.kind!r}")

    # MUTATION: raise minimum_version above the running version. Nothing else
    # changes — same version, same artifact, same digest.
    fx.manifest["minimum_version"] = "0.1.99"
    info = updater.evaluate(fx.manifest, CURRENT, key=key, frozen=True)
    if info.status != "update-manual" or info.reason != "minimum-version":
        fails.append(f"minimum_version above the running version gave {info!r} — "
                     f"the in-app update was NOT blocked")
    if not info.is_update:
        fails.append("a blocked update must still tell the user a newer version exists")
    if "installer" not in info.detail.lower():
        fails.append(f"the blocked message must send the user to the installer: "
                     f"{info.detail!r}")
    try:
        updater.plan_install(info, fx.artifact_path, frozen=True,
                            platform_name="win32")
        fails.append("plan_install ran an update below minimum_version")
    except updater.UpdateError:
        pass

    # Exactly at the minimum is allowed (it is a floor, not an exclusion).
    fx.manifest["minimum_version"] = CURRENT
    info = updater.evaluate(fx.manifest, CURRENT, key=key, frozen=True)
    if info.status != "update-available":
        fails.append(f"minimum_version == current gave {info.status!r}, expected "
                     f"'update-available' (the floor is inclusive)")

    # minimum_version is about installer upgrade paths; pip is unaffected.
    fx.manifest["minimum_version"] = "0.1.99"
    info = updater.evaluate(fx.manifest, CURRENT, key=key, frozen=False)
    if info.status != "update-available":
        fails.append(f"pip install blocked by minimum_version ({info.status!r})")
    fx.close()
    return fails


# ---------------------------------------------------- D. the checksum (mutation)

def test_checksum_refusal():
    fails = []
    fx = Fixture()
    dest = tempfile.mkdtemp(prefix="xslope-dl-")

    # The honest download: bytes arrive, digest matches, path comes back.
    seen = []
    path = updater.download_artifact(fx.manifest["artifacts"][updater.platform_key()],
                                     dest, progress=lambda d, t: seen.append((d, t)))
    if not os.path.exists(path):
        fails.append("a verified download produced no file")
    elif open(path, "rb").read() != fx.payload:
        fails.append("the downloaded bytes differ from the artifact")
    if updater.sha256_file(path) != fx.sha:
        fails.append("sha256_file disagrees with hashlib on the same bytes")
    if not seen:
        fails.append("no progress was reported during a download")
    os.unlink(path)

    # MUTATION: corrupt the published digest by one character. Same bytes, same
    # URL, same size — only the manifest's claim about them changes.
    entry = dict(fx.manifest["artifacts"][updater.platform_key()])
    bad = "0" if entry["sha256"][0] != "0" else "1"
    entry["sha256"] = bad + entry["sha256"][1:]
    try:
        updater.download_artifact(entry, dest)
        fails.append("a checksum MISMATCH was accepted — an unverified file would "
                     "have been installed")
    except updater.UpdateError as exc:
        if "checksum" not in str(exc).lower():
            fails.append(f"the refusal must say why: {exc!r}")
    leftovers = os.listdir(dest)
    if leftovers:
        fails.append(f"the rejected download was left on disk: {leftovers}")

    # And a truncated body fails the same way (the digest is the only witness).
    entry2 = dict(fx.manifest["artifacts"][updater.platform_key()])
    short = os.path.join(fx.dir, "short.exe")
    with open(short, "wb") as fh:
        fh.write(fx.payload[:-5])
    entry2["url"] = _file_url(short)
    try:
        updater.download_artifact(entry2, dest)
        fails.append("a truncated download passed verification")
    except updater.UpdateError:
        pass
    if os.listdir(dest):
        fails.append("the truncated download was left on disk")

    shutil.rmtree(dest, ignore_errors=True)
    fx.close()
    return fails


# ------------------------------------------------------- E. frozen vs pip

def test_frozen_vs_pip():
    fails = []
    fx = Fixture()
    info = updater.evaluate(fx.manifest, CURRENT, key="win32-amd64", frozen=True)
    had_frozen = hasattr(sys, "frozen")
    original = getattr(sys, "frozen", None)
    try:
        # --- pip install: sys.frozen absent -----------------------------
        if had_frozen:
            del sys.frozen
        if updater.is_frozen():
            fails.append("is_frozen() true with no sys.frozen")
        plan = updater.plan_install(info, fx.artifact_path)     # frozen inferred
        if plan.kind != "pip":
            fails.append(f"an unfrozen install planned {plan.kind!r}, expected 'pip'")
        if plan.command_text != 'pip install -U "xslope[gui]"':
            fails.append(f"the pip line must be literal: {plan.command_text!r}")
        if plan.argv:
            fails.append(f"the pip branch would SPAWN {plan.argv!r} — it must run "
                         f"nothing at all")
        if plan.quits_app:
            fails.append("the pip branch must not quit the app")

        # --- frozen app: sys.frozen set by PyInstaller ------------------
        sys.frozen = True
        if not updater.is_frozen():
            fails.append("is_frozen() false with sys.frozen set")
        plan = updater.plan_install(info, fx.artifact_path, platform_name="win32",
                                    pid=1, app_path="X.exe")
        if plan.kind != "windows-installer":
            fails.append(f"a frozen win32 app planned {plan.kind!r}")
        plan = updater.plan_install(info, fx.artifact_path, platform_name="darwin")
        if plan.kind != "macos-dmg":
            fails.append(f"a frozen darwin app planned {plan.kind!r}")
        # A frozen app with nothing downloaded never gets a plan.
        try:
            updater.plan_install(info, None, platform_name="win32")
            fails.append("a frozen plan was made with no downloaded artifact")
        except updater.UpdateError:
            pass
    finally:
        if had_frozen:
            sys.frozen = original
        elif hasattr(sys, "frozen"):
            del sys.frozen
    fx.close()
    return fails


# ---------------------------------------------------------- F. the throttle

def test_once_per_day():
    fails = []
    from datetime import datetime, timedelta, timezone
    now = datetime(2026, 7, 31, 12, 0, 0, tzinfo=timezone.utc)

    def stamp(delta):
        return updater.utc_stamp(now + delta)

    cases = [
        (None, True, "no stamp at all -> check"),
        ("", True, "an empty stamp -> check"),
        ("not a date", True, "an unparseable stamp -> check"),
        (stamp(timedelta(minutes=-1)), False, "a minute ago -> skip"),
        (stamp(timedelta(hours=-23)), False, "23 h ago -> skip"),
        (stamp(timedelta(hours=-24)), True, "exactly 24 h ago -> check"),
        (stamp(timedelta(hours=-30)), True, "30 h ago -> check"),
        (stamp(timedelta(hours=+5)), True, "a future stamp (clock moved back) -> check"),
    ]
    for last, want, why in cases:
        got = updater.due_for_check(last, now=now)
        if got != want:
            fails.append(f"due_for_check({last!r}) = {got}, expected {want} — {why}")
    if updater.utc_stamp(now) != "2026-07-31T12:00:00Z":
        fails.append(f"utc_stamp format changed: {updater.utc_stamp(now)!r}")
    # A stamp written now is not due again immediately — the round trip is the
    # property that actually throttles anything.
    if updater.due_for_check(updater.utc_stamp()):
        fails.append("a check stamped now is immediately due again")
    return fails


# ------------------------------------- G. the menu, and the two dialog paths

def _pump(predicate, timeout=20.0):
    from PySide6.QtCore import QCoreApplication
    t0 = time.time()
    while not predicate() and time.time() - t0 < timeout:
        QCoreApplication.processEvents()
        time.sleep(0.01)
    return predicate()


def _isolated(mw):
    """Point the window's controller at a throwaway settings file so the checks
    never write to the real user's preferences."""
    from PySide6.QtCore import QSettings
    tmp = tempfile.mkdtemp(prefix="xslope-settings-")
    s = QSettings(os.path.join(tmp, "u.ini"), QSettings.IniFormat)
    mw.updates.settings = s
    return s, tmp


def test_menu_and_paths():
    fails = []
    from PySide6.QtWidgets import QApplication, QMessageBox
    from studio import update_ui
    from studio.main_window import MainWindow

    QApplication.instance() or QApplication([])
    mw = MainWindow()
    settings, tmpdir = _isolated(mw)

    # --- the menu item exists, in Help, and is wired ---------------------
    # Hold a reference: a QMenu fetched from an action inside a generator is
    # collected the moment the generator is done with it.
    help_menu = None
    for act in mw.menuBar().actions():
        if "Help" in act.text():
            help_menu = act.menu()
            break
    if help_menu is None:
        fails.append("no Help menu")
    else:
        texts = [a.text() for a in help_menu.actions()]
        if not any("Check for" in t and "Updates" in t and t.endswith("…")
                   for t in texts):
            fails.append(f"Help has no 'Check for Updates…' item: {texts}")
        if not any("Startup" in t for t in texts):
            fails.append(f"Help has no startup-check preference: {texts}")
        if mw.act_check_updates not in help_menu.actions():
            fails.append("act_check_updates is not in the Help menu")
        if mw.act_startup_check not in help_menu.actions():
            fails.append("act_startup_check is not in the Help menu")
    if not mw.act_startup_check.isCheckable():
        fails.append("the startup preference is not a checkable item")

    # --- the preference: default ON, and it persists ---------------------
    if not mw.updates.startup_check_enabled():
        fails.append("the startup check must default to ON")
    mw.act_startup_check.setChecked(False)
    if settings.value(updater.SETTING_STARTUP_CHECK) not in (False, "false"):
        fails.append("unchecking the preference did not reach QSettings")
    if mw.updates.startup_check_enabled():
        fails.append("the controller still reports the startup check enabled")
    if mw.updates.check_at_startup() is not False:
        fails.append("a disabled startup check still ran")
    mw.act_startup_check.setChecked(True)
    if not mw.updates.startup_check_enabled():
        fails.append("re-checking the preference did not persist")

    # --- record every dialog instead of showing one ----------------------
    shown = []
    dialogs = []
    orig_info = QMessageBox.information
    QMessageBox.information = staticmethod(
        lambda *a, **k: shown.append(a[1:3]))
    mw.updates.show_update_dialog = lambda info, modal=True: dialogs.append(
        (info, modal))

    def run_check(manual, url):
        del shown[:], dialogs[:]
        os.environ["XSLOPE_UPDATE_URL"] = url
        settings.setValue(updater.SETTING_LAST_CHECK, "")
        done = []
        if manual:
            mw.act_check_updates.trigger()
        else:
            mw.updates.check_at_startup()
        runner = mw.updates._runner
        if runner is None:
            return False
        runner.finished.connect(lambda: done.append(True))
        _pump(lambda: bool(done) or mw.updates._runner is None)
        _pump(lambda: bool(shown) or bool(dialogs), timeout=1.0)
        return True

    try:
        up_to_date = Fixture(version=CURRENT)
        newer = Fixture(version=NEWER)
        broken = os.path.join(up_to_date.dir, "does-not-exist.json")

        from xslope import __version__ as running
        # The fixtures are written around CURRENT; if the package has moved on,
        # rebuild them around the version actually running.
        if running != CURRENT:
            up_to_date.manifest["version"] = running
            up_to_date.write()
            newer.manifest["version"] = f"{running}.1"
            newer.write()

        # 1. manual + up to date -> exactly one informational dialog, no offer
        run_check(True, up_to_date.url)
        if len(shown) != 1:
            fails.append(f"a manual up-to-date check showed {len(shown)} dialogs, "
                         f"expected exactly 1")
        elif "up to date" not in " ".join(str(x) for x in shown[0]).lower():
            fails.append(f"the up-to-date dialog does not say so: {shown[0]}")
        if dialogs:
            fails.append("an up-to-date check offered an update")

        # 2. manual + newer -> the update dialog, modal
        run_check(True, newer.url)
        if len(dialogs) != 1:
            fails.append(f"a manual check on a newer release raised {len(dialogs)} "
                         f"update dialogs")
        elif dialogs[0][1] is not True:
            fails.append("the manual update dialog must be modal")
        if shown:
            fails.append(f"the newer-version path also showed {shown}")

        # 3. manual + dead server -> exactly one could-not-check dialog
        run_check(True, _file_url(broken))
        if len(shown) != 1:
            fails.append(f"a failed manual check showed {len(shown)} dialogs, "
                         f"expected exactly 1")
        if dialogs:
            fails.append("a failed check offered an update")

        # 4. startup + dead server -> SILENCE
        run_check(False, _file_url(broken))
        if shown or dialogs:
            fails.append(f"the startup check surfaced a network failure: "
                         f"{shown} {dialogs}")

        # 5. startup + up to date -> silence
        run_check(False, up_to_date.url)
        if shown or dialogs:
            fails.append(f"the startup check spoke while up to date: {shown} {dialogs}")

        # 6. startup + newer -> one NON-modal notification
        run_check(False, newer.url)
        if len(dialogs) != 1:
            fails.append(f"the startup check raised {len(dialogs)} notifications "
                         f"for a new version")
        elif dialogs[0][1] is not False:
            fails.append("the startup notification must be non-modal")

        # 7. the throttle, through the controller: a check just ran, so the next
        #    startup check must not start one.
        settings.setValue(updater.SETTING_LAST_CHECK, updater.utc_stamp())
        if mw.updates.check_at_startup() is not False:
            fails.append("a second startup check ran within the same day")

        up_to_date.close()
        newer.close()
    finally:
        QMessageBox.information = orig_info
        os.environ.pop("XSLOPE_UPDATE_URL", None)

    # --- the dialog itself: frozen and pip look different -----------------
    info = updater.evaluate(Fixture(version=NEWER).manifest, CURRENT,
                            key="win32-amd64", frozen=True)
    dlg = update_ui.UpdateAvailableDialog(info, frozen=True, parent=mw)
    if dlg.install_btn is None:
        fails.append("a frozen app's update dialog has no Download & Install button")
    if dlg.findChild(type(dlg), "pip_command") is not None:
        fails.append("a frozen app's dialog shows a pip command")
    dlg.deleteLater()

    dlg = update_ui.UpdateAvailableDialog(info, frozen=False, parent=mw)
    if dlg.install_btn is not None:
        fails.append("a pip install's dialog offers to install for the user")
    if dlg.cmd_edit.text() != 'pip install -U "xslope[gui]"':
        fails.append(f"the pip dialog shows {dlg.cmd_edit.text()!r}")
    if not dlg.cmd_edit.isReadOnly():
        fails.append("the pip command line is editable")
    dlg._copy_command()
    from PySide6.QtWidgets import QApplication as _QA
    cb = _QA.clipboard()
    if cb is not None and cb.text() != 'pip install -U "xslope[gui]"':
        fails.append(f"Copy put {cb.text()!r} on the clipboard")
    dlg.deleteLater()

    blocked = updater.UpdateInfo("update-manual", CURRENT, latest=NEWER,
                                 reason="minimum-version", detail="too old")
    dlg = update_ui.UpdateAvailableDialog(blocked, frozen=True, parent=mw)
    if dlg.install_btn is not None:
        fails.append("a blocked update still offered Download & Install")
    if dlg.releases_btn is None:
        fails.append("a blocked update does not offer the releases page")
    dlg.deleteLater()

    mw.close()
    shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


# ------------------------------------------------------------ H. end to end

def test_end_to_end():
    """check -> new version -> download -> verify -> the hand-off command.

    The installer is never run: both platform branches are asserted on the argv
    and the helper script they WOULD hand to the OS."""
    fails = []
    fx = Fixture(version=NEWER)
    dest = tempfile.mkdtemp(prefix="xslope-e2e-")
    try:
        # 1. the check, over the wire (file://), exactly as the runner calls it
        info = updater.check_for_update(current_version=CURRENT, url=fx.url,
                                        key="win32-amd64", frozen=True)
        if info.status != "update-available":
            fails.append(f"the end-to-end check gave {info!r}")
            return fails
        if info.latest != NEWER:
            fails.append(f"latest = {info.latest!r}")
        if not info.notes_url.startswith("https://github.com/njones61/xslope"):
            fails.append(f"notes_url = {info.notes_url!r}")

        # 2. download + verify
        path = updater.download_artifact(info.artifact, dest)
        if updater.sha256_file(path) != fx.sha:
            fails.append("the downloaded artifact does not match the manifest")

        # 3a. WINDOWS: the hand-off. Assert the command, run nothing.
        plan = updater.plan_install(info, path, frozen=True, platform_name="win32",
                                    app_path=r"C:\Users\n\XSLOPE Studio.exe",
                                    pid=31337, helper_dir=dest)
        if plan.kind != "windows-installer":
            fails.append(f"windows planned {plan.kind!r}")
        if not plan.quits_app:
            fails.append("the windows update must quit the app first")
        comspec = os.environ.get("COMSPEC", "cmd.exe")
        want = [comspec, "/c", os.path.join(dest, "xslope_update.cmd"), "31337",
                path, r"C:\Users\n\XSLOPE Studio.exe"]
        if plan.argv != want:
            fails.append(f"windows argv:\n  got  {plan.argv}\n  want {want}")
        if not os.path.exists(plan.helper):
            fails.append("the windows helper script was not written")
        else:
            # newline="" so the CRLF the batch file needs survives the read.
            helper = open(plan.helper, encoding="ascii", newline="").read()
            for needed, why in [
                ('"%SETUP%" /S', "the silent NSIS upgrade"),
                ('start "" "%APPEXE%"', "the relaunch after installing"),
                ('tasklist /FI "PID eq %WAITPID%"', "waiting for our own exit"),
                ('del "%~f0"', "the helper deleting itself"),
            ]:
                if needed not in helper:
                    fails.append(f"the windows helper is missing {why}: {needed!r}")
            if "\r\n" not in helper:
                fails.append("the windows helper needs CRLF line endings")

        # 3b. macOS, a bundle it can replace: the helper, then quit. Run nothing.
        bundle = os.path.join(dest, "XSLOPE Studio.app")
        os.makedirs(os.path.join(bundle, "Contents", "MacOS"))
        exe = os.path.join(bundle, "Contents", "MacOS", "XSLOPE Studio")
        plan = updater.plan_install(info, path, frozen=True, platform_name="darwin",
                                    app_path=exe, pid=31337, helper_dir=dest)
        if plan.kind != "macos-replace":
            fails.append(f"darwin with a writable bundle planned {plan.kind!r}")
        if not plan.quits_app:
            fails.append("the macOS replace must quit the app first")
        want = ["/bin/sh", os.path.join(dest, "xslope_update.sh"), "31337", path,
                bundle, "interactive"]
        if plan.argv != want:
            fails.append(f"darwin argv:\n  got  {plan.argv}\n  want {want}")
        if not os.path.exists(plan.helper or ""):
            fails.append("the macOS helper script was not written")
        else:
            helper = open(plan.helper, encoding="utf-8").read()
            for needed, why in [
                ('kill -0 "$WAITPID"', "waiting for our own exit"),
                ('/usr/bin/hdiutil attach -nobrowse', "mounting the dmg quietly"),
                ('/usr/bin/ditto "$SRC" "$NEW"', "copying the new bundle with ditto"),
                ('/bin/mv "$NEW" "$APP"', "swapping the new bundle in"),
                ('/bin/mv "$OLD" "$APP"', "putting the old bundle back on a failure"),
                ('/usr/bin/open "$APP"', "the relaunch"),
                ('/bin/rm -f "$0"', "the helper deleting itself"),
            ]:
                if needed not in helper:
                    fails.append(f"the macOS helper is missing {why}: {needed!r}")
            if not os.access(plan.helper, os.X_OK):
                fails.append("the macOS helper is not executable")
        if "right-click" in plan.message.lower() or "signed" in plan.message.lower():
            fails.append(f"the macOS message still talks about signing: {plan.message!r}")

        # 3c. macOS, a bundle it cannot replace: open the dmg and instruct.
        trans = os.path.join(dest, "AppTranslocation", "X", "XSLOPE Studio.app")
        os.makedirs(trans)
        for label, where in (("not an app", os.path.join(dest, "python3")),
                             ("translocated", trans)):
            plan = updater.plan_install(info, path, frozen=True,
                                        platform_name="darwin", app_path=where)
            if plan.kind != "macos-dmg":
                fails.append(f"darwin ({label}) planned {plan.kind!r}")
                continue
            if plan.argv != ["/usr/bin/open", path]:
                fails.append(f"darwin ({label}) argv = {plan.argv!r}, expected open <dmg>")
            if plan.quits_app:
                fails.append(f"darwin ({label}) must not quit the app")
            for phrase in ("Applications", "checksum"):
                if phrase.lower() not in plan.message.lower():
                    fails.append(f"the macOS instructions ({label}) never mention "
                                 f"{phrase!r}")
            if "right-click" in plan.message.lower() or "signed" in plan.message.lower():
                fails.append(f"the macOS fallback still talks about signing ({label})")
    finally:
        shutil.rmtree(dest, ignore_errors=True)
        fx.close()
    return fails


def test_download_path():
    """The frozen app's download, driven through the real controller and runner:
    Download & Install → progress → sha256 → the hand-off. And the same run with a
    corrupted digest, which must reach the failure dialog and never a plan."""
    fails = []
    from PySide6.QtWidgets import QApplication, QMessageBox
    from studio.main_window import MainWindow

    QApplication.instance() or QApplication([])
    mw = MainWindow()
    _isolated(mw)
    fx = Fixture(version=NEWER)
    plans, errors = [], []
    mw.updates.run_plan = plans.append
    orig_critical = QMessageBox.critical
    QMessageBox.critical = staticmethod(lambda *a, **k: errors.append(a[1:3]))
    had_frozen = hasattr(sys, "frozen")
    original = getattr(sys, "frozen", None)
    sys.frozen = True                     # take the frozen branch on this machine
    try:
        info = updater.evaluate(fx.manifest, CURRENT, key=updater.platform_key(),
                                frozen=True)
        # --- the happy path ---------------------------------------------
        mw.updates.start_download(info)
        if not mw.progress_bar.isVisible() and mw.isVisible():
            fails.append("the download shows no progress bar")
        if mw._update_dl_runner is None:
            fails.append("no download runner was started")
        _pump(lambda: bool(plans) or bool(errors), timeout=30.0)
        if errors:
            fails.append(f"a good download reported a failure: {errors}")
        if len(plans) != 1:
            fails.append(f"the download produced {len(plans)} install plans")
        else:
            plan = plans[0]
            want = {"win32": "windows-installer", "darwin": "macos-dmg"}.get(
                sys.platform, "manual")
            if plan.kind != want:
                fails.append(f"on {sys.platform} the plan was {plan.kind!r}, "
                             f"expected {want!r}")
        _pump(lambda: mw._update_dl_runner is None, timeout=5.0)
        if mw.progress_bar.isVisible():
            fails.append("the progress bar was left up after the download")

        # --- MUTATION: the manifest's digest no longer describes the bytes
        del plans[:], errors[:]
        entry = info.artifact
        entry["sha256"] = ("0" if entry["sha256"][0] != "0" else "1") + entry["sha256"][1:]
        mw.updates.start_download(info)
        _pump(lambda: bool(plans) or bool(errors), timeout=30.0)
        if plans:
            fails.append("a checksum mismatch still produced an install plan — "
                         "Studio would have run an unverified file")
        if len(errors) != 1:
            fails.append(f"a checksum mismatch raised {len(errors)} error dialogs")
    finally:
        QMessageBox.critical = orig_critical
        if had_frozen:
            sys.frozen = original
        elif hasattr(sys, "frozen"):
            del sys.frozen
        _pump(lambda: mw._update_dl_runner is None, timeout=5.0)
        mw.close()
        fx.close()
    return fails


def _bundle(path, version):
    os.makedirs(os.path.join(path, "Contents", "MacOS"), exist_ok=True)
    with open(os.path.join(path, "Contents", "version.txt"), "w") as fh:
        fh.write(version)


def _bundle_version(path):
    try:
        return open(os.path.join(path, "Contents", "version.txt")).read()
    except OSError:
        return None


def test_macos_swap():
    """The macOS helper, run for real: a dmg holding a new bundle replaces the
    installed one; a dmg with no app, and a folder it cannot write to, leave the
    installed bundle as it was."""
    if sys.platform != "darwin":
        print("  (macOS swap: not on a Mac, skipped)")
        return []
    import subprocess
    fails = []
    root = tempfile.mkdtemp(prefix="xslope-swap-")
    try:
        dead = subprocess.Popen(["/usr/bin/true"])
        dead.wait()

        def make_dmg(name, with_app):
            src = os.path.join(root, name + "-src")
            os.makedirs(src)
            if with_app:
                _bundle(os.path.join(src, "XSLOPE Studio.app"), "new")
            else:
                open(os.path.join(src, "README.txt"), "w").close()
            dmg = os.path.join(root, name + ".dmg")
            subprocess.run(["/usr/bin/hdiutil", "create", "-quiet", "-volname", name,
                            "-srcfolder", src, "-ov", "-format", "UDZO", dmg],
                           check=True)
            return dmg

        def run(dmg, bundle):
            argv, helper = updater.macos_replace_in_place(
                dmg, bundle, pid=dead.pid, helper_dir=root, mode="quiet")
            done = subprocess.run(argv, capture_output=True, timeout=120)
            return done.returncode, helper

        # the swap
        install = os.path.join(root, "Applications")
        bundle = os.path.join(install, "XSLOPE Studio.app")
        _bundle(bundle, "old")
        dmg = make_dmg("good", True)
        code, helper = run(dmg, bundle)
        if code != 0:
            fails.append(f"the swap reported failure (exit {code})")
        if _bundle_version(bundle) != "new":
            fails.append(f"after the swap the installed bundle reads "
                         f"{_bundle_version(bundle)!r}, not the new one")
        leftovers = [n for n in os.listdir(install) if n != "XSLOPE Studio.app"]
        if leftovers:
            fails.append(f"the swap left {leftovers} beside the bundle")
        if os.path.exists(dmg):
            fails.append("the dmg was not deleted after a good swap")
        if os.path.exists(helper):
            fails.append("the helper did not delete itself")

        # MUTATION: a dmg with no app in it
        _bundle(bundle, "old")
        dmg = make_dmg("empty", False)
        code, helper = run(dmg, bundle)
        if code == 0:
            fails.append("a dmg with no app in it reported success")
        if _bundle_version(bundle) != "old":
            fails.append("a dmg with no app in it changed the installed bundle")
        if not os.path.exists(dmg):
            fails.append("a failed swap deleted the dmg the user needs to finish by hand")

        # MUTATION: a folder the swap cannot write to
        dmg = make_dmg("good2", True)
        os.chmod(install, 0o555)
        try:
            code, helper = run(dmg, bundle)
        finally:
            os.chmod(install, 0o755)
        if code == 0:
            fails.append("a read-only folder reported success")
        if _bundle_version(bundle) != "old":
            fails.append("a read-only folder changed the installed bundle")
        if sorted(os.listdir(install)) != ["XSLOPE Studio.app"]:
            fails.append(f"a failed swap left {os.listdir(install)} in the folder")
        if updater.macos_replace_blocker(bundle) is not None:
            fails.append(f"a writable bundle was refused: "
                         f"{updater.macos_replace_blocker(bundle)}")
        os.chmod(install, 0o555)
        try:
            if updater.macos_replace_blocker(bundle) is None:
                fails.append("a bundle in a read-only folder was not refused")
        finally:
            os.chmod(install, 0o755)
    finally:
        for dirpath, dirnames, _files in os.walk(root):
            for d in dirnames:
                try:
                    os.chmod(os.path.join(dirpath, d), 0o755)
                except OSError:
                    pass
        shutil.rmtree(root, ignore_errors=True)
    return fails


CHECKS = [("version comparison", test_version_comparison),
          ("the platform artifact key", test_platform_key),
          ("minimum_version blocks the install", test_minimum_version),
          ("a bad checksum refuses to install", test_checksum_refusal),
          ("frozen vs pip branch", test_frozen_vs_pip),
          ("the once-per-day throttle", test_once_per_day),
          ("the Help menu and both check paths", test_menu_and_paths),
          ("end to end: check -> verify -> hand-off", test_end_to_end),
          ("the macOS swap, run", test_macos_swap),
          ("Download & Install through the controller", test_download_path)]

#: The legs that need Qt; everything else runs on an engine-only install.
_QT_CHECKS = (test_menu_and_paths, test_download_path)


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
        checks = CHECKS
    except Exception:
        print("updater: PySide6 not installed — Studio checks skipped.")
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
    print("Studio updater:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll updater checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
