"""Standing check on Studio's Help → Report a Problem… and the bug-report form it opens.

The menu item opens GitHub's bug-report form with the version and the platform
already entered. GitHub fills a form's fields from query parameters named by the
fields' ids, and fills a dropdown only from an option spelled exactly as listed,
so the menu item and ``.github/ISSUE_TEMPLATE/bug_report.yml`` must agree: the
ids ``version`` and ``platform`` must exist, and every label Studio can send must
be one of the dropdown's options. A rename on either side breaks nothing visibly;
the form just opens empty.

The window is built offscreen and the browser call is intercepted, so nothing
opens.

Skips cleanly (exit 0) when PySide6 is not installed.
"""
import os
import sys
from urllib.parse import parse_qs, urlparse

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

FORM = os.path.join(_REPO, ".github", "ISSUE_TEMPLATE", "bug_report.yml")
#: Every label studio_platform_label() can return, whatever this machine is.
LABELS = ("XSLOPE Studio on macOS", "XSLOPE Studio on Windows", "Other")


def _form_fields():
    import yaml
    with open(FORM, encoding="utf-8") as fh:
        form = yaml.safe_load(fh)
    return {item.get("id"): item for item in form.get("body", []) if item.get("id")}


def run():
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("  PySide6 not installed; skipped")
        return []
    QApplication.instance() or QApplication([])
    from PySide6.QtGui import QDesktopServices

    from studio import main_window
    from xslope import __version__

    fails = []

    # 1. The form: the two ids exist, and the dropdown lists every label.
    fields = _form_fields()
    for fid in ("version", "platform"):
        if fid not in fields:
            fails.append(f"bug_report.yml has no field with id '{fid}', so the "
                         f"menu item's {fid} is dropped")
    options = (fields.get("platform", {}).get("attributes", {}).get("options") or [])
    for label in LABELS:
        if label not in options:
            fails.append(f"the form's platform dropdown has no option {label!r}")
    if main_window.studio_platform_label() not in LABELS:
        fails.append(f"Studio sends an unlisted platform label "
                     f"{main_window.studio_platform_label()!r}")

    # 2. The menu: the item is in Help, and it opens the form with both filled.
    opened = []
    real = QDesktopServices.openUrl
    QDesktopServices.openUrl = staticmethod(lambda url: opened.append(url.toString()) or True)
    mw = main_window.MainWindow()
    try:
        from PySide6.QtWidgets import QMenu
        help_menu = next((m for m in mw.menuBar().findChildren(QMenu)
                          if m.title().replace("&", "") == "Help"), None)
        texts = [a.text().replace("&", "") for a in help_menu.actions()] if help_menu else []
        if "Report a Problem…" not in texts:
            fails.append(f"Help has no Report a Problem… item (has {texts})")
        mw.act_report_problem.trigger()
    finally:
        QDesktopServices.openUrl = real
        mw.doc._dirty = False
        mw.close()
    if not opened:
        fails.append("Report a Problem… opened nothing")
    else:
        url = urlparse(opened[-1])
        query = {k: v[0] for k, v in parse_qs(url.query).items()}
        if not url.path.endswith("/issues/new"):
            fails.append(f"it opened {opened[-1]}, not the new-issue page")
        if query.get("template") != "bug_report.yml":
            fails.append("it does not open the bug-report form")
        if query.get("version") != __version__:
            fails.append(f"version {query.get('version')!r}, expected {__version__!r}")
        if query.get("platform") != main_window.studio_platform_label():
            fails.append(f"platform {query.get('platform')!r} was not sent")
    return fails


def main():
    print("Help → Report a Problem…:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("  the menu item opens the bug-report form with the version and platform")
    print("\nAll report-a-problem checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    main()
