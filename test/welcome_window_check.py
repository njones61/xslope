"""Standing checks on the welcome window and the Help menu's way out of the app.

The first launch is the one moment Studio gets to say where everything else is, and
every way that can rot is a quiet one — a link that still opens, at a page that was
renamed; a preference that is written but never read, so the window either never
comes back or never goes away; a menu item that lost its slot. What is guarded:

  A. THE WINDOW SAYS WHAT IT WAS RULED TO SAY. Four documentation links, the
     version it is running, and the support address as a mail link — and NO
     bug-report link, which was ruled out of this window deliberately (it is a Help
     menu item's job). The hrefs are read out of the rendered text, so a label that
     says "Verification" over the samples URL fails here; the bug-report refusal is
     read off every widget and action on the window, because a Report a Bug button
     is a button and a guard that reads labels alone would pass one.
  B. THE PREFERENCE ROUND-TRIPS, AND IT IS THE SAME KEY BOTH WAYS. An unset key is
     a first launch: the window opens, with the tick already in the box, so the
     welcome is a first-launch window and not a greeting at every start. Closing
     writes the answer by whatever route the window was dismissed — a preference
     that only the Close button writes is lost by pressing Escape — and unticking
     brings the launch greeting back. Checked against a throwaway settings file:
     nothing here reads or writes the preferences of whoever runs the suite.
  C. THE LAUNCH ASKS, AND OBEYS THE ANSWER. The window is raised by
     ``studio.app.main`` rather than by the MainWindow constructor, so that
     building a window — which every Studio test does — never raises a modal
     dialog. Both answers are exercised, and so is the ``XSLOPE_NO_WELCOME``
     suppression, because a check that only proves the yes-path cannot tell a gate
     from an unconditional call. A launch that arrived with work in hand opens the
     work and no greeting over it, by either arrival route — a document on the
     command line, and the event macOS sends instead — and leaves the preference
     alone, so the next bare launch is still greeted.
  D. THE HELP MENU CARRIES BOTH ITEMS, AND DOCUMENTATION OPENS THE DOCS. Order
     included: the outward items first, then the update items, then About.
     ``QDesktopServices.openUrl`` is intercepted, so the check reads the URL the
     menu would open and no browser is launched.
  E. THE ADDRESSES ARE REAL. Every link is on the documentation site — the same
     host ``mkdocs.yml`` publishes and the same one the ``xslope://`` handler
     allowlists — and each one names a page that exists in ``docs/``. A page
     renamed in the documentation tree fails here rather than in a user's browser.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import os
import re
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

#: A real project to launch WITH, so the "arrived with work in hand" case is a
#: launch that actually opened something rather than one that ignored a bad path.
SAMPLE = os.path.join(_REPO, "docs", "inputs", "slope", "xslope_dam.xlsx")

#: The image the documentation page embeds, and the page that embeds it.
DOCS_PAGE = os.path.join(_REPO, "docs", "studio", "index.md")
DOCS_IMAGE = os.path.join(_REPO, "docs", "studio", "images",
                          "overview_welcome_window.png")


def _hrefs(label):
    """Every anchor target in a rich-text label, in the order they are shown."""
    return re.findall(r'href="([^"]+)"', label.text())


#: The ways a Qt widget carries text a person can read: labels and buttons say
#: ``text()``, group boxes ``title()``, an input its ``placeholderText()``, and
#: anything at all a ``toolTip()``.
_TEXT_GETTERS = ("text", "title", "placeholderText", "toolTip", "whatsThis")


def _visible_text(widget):
    """Everything readable on ``widget`` — its title, every child widget's text,
    and every action on it.

    Deliberately not "every QLabel": what a window offers is offered by whatever
    widget it is offered with, so a guard that reads one widget class only proves
    something about that class. Getters that want arguments (``QComboBox.text`` is
    not one, ``QTableWidget.item`` is) are skipped rather than guessed at.
    """
    from PySide6.QtGui import QAction
    from PySide6.QtWidgets import QWidget

    seen = [widget.windowTitle()]
    for child in [widget] + widget.findChildren(QWidget):
        for name in _TEXT_GETTERS:
            getter = getattr(child, name, None)
            if not callable(getter):
                continue
            try:
                value = getter()
            except TypeError:              # a getter that takes an index or a role
                continue
            if isinstance(value, str) and value:
                seen.append(value)
    for action in list(widget.actions()) + widget.findChildren(QAction):
        seen.append(action.text())
        seen.append(action.toolTip())
    return " ".join(s for s in seen if s)


def _settings(tmp, name="w.ini"):
    from PySide6.QtCore import QSettings
    return QSettings(os.path.join(tmp, name), QSettings.IniFormat)


# ------------------------------------------------------------------ A. contents
def test_contents():
    fails = []
    import xslope
    from PySide6.QtWidgets import QDialogButtonBox
    from studio import links
    from studio.welcome import WelcomeDialog

    with tempfile.TemporaryDirectory() as tmp:
        dlg = WelcomeDialog(settings=_settings(tmp))

        want = [links.DOCS_URL, links.GETTING_STARTED_URL, links.SAMPLES_URL,
                links.VERIFICATION_URL]
        got = _hrefs(dlg.link_label)
        if got != want:
            fails.append(f"the window's links are {got}, want {want}")
        if not dlg.link_label.openExternalLinks():
            fails.append("the links are not handed to the browser — clicking one "
                         "does nothing")

        # Each link is named, and named for what it opens: an href with no visible
        # text beside it is a link nobody can choose between.
        text = dlg.link_label.text()
        for name in ("Documentation", "Getting started", "Sample problems",
                     "Verification"):
            if f">{name}</a>" not in text:
                fails.append(f"no link labelled {name!r}")

        if _hrefs(dlg.support_label) != [links.SUPPORT_MAILTO]:
            fails.append(f"the support address is "
                         f"{_hrefs(dlg.support_label)}, want {links.SUPPORT_MAILTO}")
        if links.SUPPORT_EMAIL not in dlg.support_label.text():
            fails.append("the support address is linked but not shown, so it "
                         "cannot be read off the screen or copied")

        if xslope.__version__ not in dlg.version_label.text():
            fails.append(f"the version note does not name the running version "
                         f"({xslope.__version__})")
        if "Updates" not in dlg.version_label.text():
            fails.append("the version note does not say how a newer release is "
                         "found")

        # Ruled out of this window (plan_docs_studio_links.md §6): reporting a bug
        # is a Help menu item, not something the first launch asks for. Read off
        # EVERY widget on the window and every action on it, not the labels alone —
        # a Report a Bug button is a button, and a guard that only reads labels is a
        # guard that would have let one through.
        whole = _visible_text(dlg)
        if re.search(r"report (a )?(bug|problem)|issues/new|github\.com|"
                     r"file an issue|bug report", whole, re.I):
            fails.append("the welcome window carries a bug-report link or button, "
                         "which was ruled out of it")

        buttons = dlg.findChildren(QDialogButtonBox)
        if not buttons or not buttons[0].button(QDialogButtonBox.Close):
            fails.append("the window has no Close button")
        dlg.deleteLater()
    return fails


# --------------------------------------------------------------- B. the setting
def test_setting():
    fails = []
    from studio.welcome import (SETTING_SHOW_AT_LAUNCH, WelcomeDialog,
                                show_at_launch)

    with tempfile.TemporaryDirectory() as tmp:
        s = _settings(tmp)
        if s.value(SETTING_SHOW_AT_LAUNCH) is not None:
            fails.append("a fresh settings file already carries the preference")
        if not show_at_launch(s):
            fails.append("an unset preference does not show the window — a first "
                         "launch would never see it")

        dlg = WelcomeDialog(settings=s)
        if not dlg.dont_show.isChecked():
            fails.append("the first launch opens with 'Don't show this again' "
                         "unticked, so the welcome would greet every launch")
        dlg.show()
        dlg.close()                      # the title-bar route, not the button
        if show_at_launch(s):
            fails.append("closing with the box ticked left the window still due "
                         "at launch")
        if s.value(SETTING_SHOW_AT_LAUNCH) is None:
            fails.append("closing the window wrote no preference at all")

        # Reopened from Help → Welcome: the box shows the preference as it stands,
        # and unticking it is how the launch greeting comes back.
        dlg = WelcomeDialog(settings=s)
        if not dlg.dont_show.isChecked():
            fails.append("a reopened window does not show the stored preference")
        dlg.show()
        dlg.dont_show.setChecked(False)
        dlg.reject()                     # Escape
        if not show_at_launch(s):
            fails.append("unticking the box did not bring the launch greeting back")

        # And back off again, so the key is not one-way.
        dlg = WelcomeDialog(settings=s)
        if dlg.dont_show.isChecked():
            fails.append("the box is ticked although the preference says to show "
                         "the window at launch")
        dlg.show()
        dlg.dont_show.setChecked(True)
        dlg.accept()                     # the Close button
        if show_at_launch(s):
            fails.append("ticking the box again did not switch the greeting off")

        # An unreadable stored value is a shown window, not a crash: the store is a
        # text file a user can edit.
        s.setValue(SETTING_SHOW_AT_LAUNCH, "perhaps")
        if not show_at_launch(s):
            fails.append("an unreadable stored value suppressed the window")
    return fails


# ---------------------------------------------------------------- C. the launch
def test_launch_gate():
    fails = []
    from PySide6.QtWidgets import QApplication
    from studio import app as studio_app
    from studio import welcome
    from studio.main_window import MainWindow

    real_gate = welcome.show_at_launch
    real_show = MainWindow.show_welcome
    real_exec = QApplication.exec
    calls = []

    def _run(answer, suppress, args=()):
        """One ``studio.app.main`` launch, with the gate answering ``answer``."""
        del calls[:]
        welcome.show_at_launch = lambda settings=None: answer
        MainWindow.show_welcome = lambda self: calls.append(self)
        # exec() pumps the deferred welcome and returns instead of running a GUI.
        QApplication.exec = lambda self=None: (
            [QApplication.processEvents() for _ in range(20)] and 0 or 0)
        os.environ["XSLOPE_NO_UPDATE_CHECK"] = "1"
        if suppress:
            os.environ["XSLOPE_NO_WELCOME"] = "1"
        else:
            os.environ.pop("XSLOPE_NO_WELCOME", None)
        try:
            studio_app.main(["xslope-studio", *args])
        finally:
            os.environ.pop("XSLOPE_NO_WELCOME", None)
            for w in QApplication.topLevelWidgets():
                if isinstance(w, MainWindow):
                    w.stop_threads()
                    w.close()
        return len(calls)

    try:
        if _run(answer=True, suppress=False) != 1:
            fails.append("a launch whose preference asks for the welcome window "
                         "did not open one")
        if _run(answer=False, suppress=False) != 0:
            fails.append("a launch opened the welcome window although the "
                         "preference says not to")
        if _run(answer=True, suppress=True) != 0:
            fails.append("XSLOPE_NO_WELCOME did not suppress the window")
        # A launch that arrived with work in hand — a document on the command line,
        # which is how a double-clicked file reaches Studio on Windows and Linux —
        # opens that document and no greeting over it. The preference is untouched,
        # so the next bare launch still shows the window.
        if os.path.exists(SAMPLE):
            if _run(answer=True, suppress=False, args=(SAMPLE,)) != 0:
                fails.append("a launch that was handed a file still opened the "
                             "welcome window over it")
        else:
            fails.append(f"missing fixture {os.path.relpath(SAMPLE, _REPO)}")
        if _run(answer=True, suppress=False) != 1:
            fails.append("a bare launch after a file launch showed no welcome "
                         "window, so opening a file switched it off for good")
    finally:
        welcome.show_at_launch = real_gate
        MainWindow.show_welcome = real_show
        QApplication.exec = real_exec
        os.environ.pop("XSLOPE_NO_UPDATE_CHECK", None)

    # What "work in hand" is decided by: a request that was ACTED ON marks the
    # window, in ``open_request`` — the one place both arrival routes meet — so a
    # command-line argument and the event macOS sends instead are answered the same
    # way, and an argument that named nothing openable is not work.
    class _Window:
        def __init__(self):
            self.opened = []

        def open_path(self, path):
            self.opened.append(path)

        def open_scheme_url(self, url):
            self.opened.append(url)

    win = _Window()
    studio_app.open_request(win, SAMPLE)
    if win.opened != [SAMPLE]:
        fails.append(f"a file argument reached {win.opened}, not the normal open")
    if not getattr(win, "opened_by_request", False):
        fails.append("opening a file from the command line did not mark the "
                     "launch as one that arrived with work in hand")

    ignored = _Window()
    studio_app.open_request(ignored, os.path.join(_REPO, "no_such_project.xlsx"))
    if getattr(ignored, "opened_by_request", False):
        fails.append("an argument that named nothing openable counted as work in "
                     "hand, so a mistyped path would suppress the welcome window")

    # macOS delivers a double-clicked document as an event rather than on the
    # command line, and when the double-click is what launched Studio the request
    # arrives before there is a window to give it to. Held and replayed, it must
    # reach the same mark. Qt refuses a second application object in one process,
    # so the real methods are called against a stand-in holding the same state.
    class _Stub:
        _deliver = studio_app.StudioApplication._deliver

    sa = _Stub()
    sa._window, sa._pending = None, []
    studio_app.StudioApplication._deliver(sa, SAMPLE)          # arrives too early
    if sa._pending != [SAMPLE]:
        fails.append("a request that arrived before the window was dropped rather "
                     "than held")
    held = _Window()
    studio_app.StudioApplication.set_window(sa, held)
    if not getattr(held, "opened_by_request", False):
        fails.append("a request held until the window existed did not mark the "
                     "launch — the welcome would open over a double-clicked "
                     "document on macOS")

    # Building a window must never raise a dialog on its own: the tests, the
    # captures and the report pipeline all construct one.
    MainWindow.show_welcome = lambda self: calls.append(self)
    del calls[:]
    try:
        mw = MainWindow()
        if calls:
            fails.append("constructing a MainWindow opened the welcome window")
        mw.stop_threads()
        mw.close()
    finally:
        MainWindow.show_welcome = real_show
    return fails


# ------------------------------------------------------------- D. the Help menu
def test_help_menu():
    fails = []
    from PySide6.QtGui import QDesktopServices
    from PySide6.QtWidgets import QDialog
    from studio import links
    from studio.main_window import MainWindow
    from studio.welcome import WelcomeDialog

    mw = MainWindow()
    try:
        help_menu = None
        for act in mw.menuBar().actions():
            if "Help" in act.text():
                help_menu = act.menu()       # held: a QMenu from a loose action
                break                        # is collected the moment it is loose
        if help_menu is None:
            return ["no Help menu"]

        actions = help_menu.actions()
        texts = [a.text() for a in actions]
        for want, action in (("Documentation", mw.act_documentation),
                             ("Welcome", mw.act_welcome)):
            if action not in actions:
                fails.append(f"the {want} action is not in the Help menu: {texts}")
            if not any(want in t for t in texts):
                fails.append(f"Help has no {want} item: {texts}")

        def _at(action):
            return actions.index(action) if action in actions else -1

        if not (_at(mw.act_documentation) < _at(mw.act_welcome)
                < _at(mw.act_check_updates) < _at(mw.act_about)):
            fails.append(f"the Help menu reads {texts}, want Documentation, "
                         "Welcome, the update items, then About")

        # --- Documentation opens the docs root, and nothing else -------------
        opened = []
        real_open = QDesktopServices.openUrl
        QDesktopServices.openUrl = staticmethod(
            lambda url: opened.append(url.toString()) or True)
        try:
            mw.act_documentation.trigger()
        finally:
            QDesktopServices.openUrl = real_open
        if opened != [links.DOCS_URL]:
            fails.append(f"Help → Documentation opened {opened}, "
                         f"want [{links.DOCS_URL!r}]")

        # --- Welcome opens the welcome window, on the app's own settings -----
        shown = []
        real_exec = QDialog.exec
        QDialog.exec = lambda self: shown.append(self) or 0
        try:
            mw.act_welcome.trigger()
        finally:
            QDialog.exec = real_exec
        if len(shown) != 1 or not isinstance(shown[0], WelcomeDialog):
            fails.append(f"Help → Welcome showed {shown}, want one WelcomeDialog")
        elif shown[0]._settings is not mw.settings:
            fails.append("the welcome window opened on a settings store of its "
                         "own, so its checkbox would not match the app's")
    finally:
        mw.stop_threads()
        mw.close()
    return fails


# --------------------------------------------------------------- E. the targets
def test_targets():
    fails = []
    from urllib.parse import urlparse

    from studio import links, urlscheme

    urls = {"Documentation": links.DOCS_URL,
            "Getting started": links.GETTING_STARTED_URL,
            "Sample problems": links.SAMPLES_URL,
            "Verification": links.VERIFICATION_URL}

    # The documentation root is the site mkdocs publishes.
    with open(os.path.join(_REPO, "mkdocs.yml"), encoding="utf-8") as fh:
        mkdocs = fh.read()
    site_url = re.search(r"^site_url:\s*(\S+)\s*$", mkdocs, re.M)
    if site_url is None:
        fails.append("mkdocs.yml declares no site_url to check the links against")
    elif links.DOCS_URL.rstrip("/") != site_url.group(1).rstrip("/"):
        fails.append(f"DOCS_URL is {links.DOCS_URL}, but mkdocs publishes "
                     f"{site_url.group(1)}")

    host = urlparse(links.DOCS_URL).hostname
    if host not in urlscheme.ALLOWED_HOSTS:
        fails.append(f"the documentation host {host} is not the host the "
                     f"xslope:// handler fetches from ({urlscheme.ALLOWED_HOSTS})")

    for name, url in urls.items():
        if not url.startswith(links.DOCS_URL):
            fails.append(f"the {name} link ({url}) is not on the documentation site")
            continue
        if urlparse(url).scheme != "https":
            fails.append(f"the {name} link is not https: {url}")
        rel = url[len(links.DOCS_URL):]
        # mkdocs serves directory URLs, so ``a/b/`` is ``docs/a/b.md`` or
        # ``docs/a/b/index.md``; the root is ``docs/index.md``.
        stem = rel.strip("/")
        candidates = [os.path.join(_REPO, "docs", *(stem.split("/") if stem else []),
                                   "index.md")]
        if stem:
            candidates.append(os.path.join(_REPO, "docs", *stem.split("/")) + ".md")
        if not any(os.path.exists(p) for p in candidates):
            fails.append(f"the {name} link points at {url}, which is no page in "
                         f"docs/ (looked for "
                         f"{[os.path.relpath(p, _REPO) for p in candidates]})")

    if links.SUPPORT_MAILTO != "mailto:" + links.SUPPORT_EMAIL:
        fails.append("the support mailto: and the address shown disagree")
    if "@" not in links.SUPPORT_EMAIL:
        fails.append(f"the support address is not an address: {links.SUPPORT_EMAIL}")

    # The documentation section that shows the window, and the image it embeds.
    if not os.path.exists(DOCS_IMAGE):
        fails.append(f"the welcome-window screenshot is missing: "
                     f"{os.path.relpath(DOCS_IMAGE, _REPO)} "
                     "(tools/capture_studio_screenshots.py writes it)")
    with open(DOCS_PAGE, encoding="utf-8") as fh:
        page = fh.read()
    if os.path.basename(DOCS_IMAGE) not in page:
        fails.append(f"{os.path.relpath(DOCS_PAGE, _REPO)} does not show the "
                     "welcome window")
    return fails


CHECKS = [("what the window says", test_contents),
          ("the show-at-launch preference", test_setting),
          ("the launch gate (studio.app)", test_launch_gate),
          ("the Help menu + Documentation", test_help_menu),
          ("the link targets", test_targets)]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
    except Exception:
        print("welcome window: PySide6 not installed — checks skipped.")
        return []
    failures = []
    for name, fn in CHECKS:
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
    print("welcome window + Help menu (Studio):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll welcome window checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication, QMessageBox
        _app = QApplication.instance() or QApplication([])
        QMessageBox.warning = staticmethod(lambda *a, **k: QMessageBox.Ok)
        QMessageBox.information = staticmethod(lambda *a, **k: QMessageBox.Ok)
        QMessageBox.critical = staticmethod(lambda *a, **k: QMessageBox.Ok)
    except Exception:
        pass
    main()
