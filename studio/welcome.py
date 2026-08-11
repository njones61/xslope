"""The welcome window — what a first launch says, and where it points.

Studio opens on an empty canvas, so the first launch is the one moment the
application has to say what it is and where everything else lives. This window is
that: a sentence about the application, four links into the documentation, the
version it is running, and the address to write to. It is short on purpose — the
documentation is the manual, and a welcome screen that tries to be one is a screen
people close without reading.

**The links go to documentation pages, not to anything in the application.** The
documentation organizes every sample and verification problem already, so there is
no in-app list to keep in step with it (``plan_docs_studio_links.md`` §3). The
addresses themselves live in :mod:`studio.links`.

**When it appears.** On the first launch, and after that only if the reader asks
for it — ``welcome/show_at_launch`` in the app's settings store carries the answer,
and it is unset until this window has been closed once. **Don't show this again**
is therefore already ticked when the window opens with nothing stored: reading it
once is what it is for. Closing writes the preference by any route (the button,
Escape, the title bar), because a window dismissed with Escape is dismissed.
**Help → Welcome** reopens it whatever the preference says, which is what makes
turning it off a safe thing to do — and unticking the box there is how a launch
gets it back.
"""

from __future__ import annotations

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QLabel, QVBoxLayout,
)

from . import links

#: The preference, in the app's own ``QSettings(ORG_NAME, SETTINGS_APP)`` store —
#: the same store the update check's ``updates/…`` keys live in.
SETTING_SHOW_AT_LAUNCH = "welcome/show_at_launch"

#: How long a line of the window's text is, in characters. The pixel width is
#: measured from the platform's own font (as :class:`studio.dialogs.
#: UnpackPackageDialog` measures its path field), so the paragraphs wrap at the
#: same place whatever font and scaling the machine is using.
TEXT_COLUMNS = 68


def _store(settings=None):
    """The settings store to read and write. The literal organization/application
    pair is the one :mod:`studio.main_window` opens (and ``studio.editors``, and
    the report dialog); it is the on-disk settings path, so it is never derived."""
    return settings if settings is not None else QSettings("XSlope", "XSlope Studio")


def _as_bool(value, default):
    """QSettings hands back "true"/"false" strings on some platforms and real bools
    on others — and ``None`` when the key has never been written."""
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in ("true", "1", "yes"):
        return True
    if text in ("false", "0", "no"):
        return False
    return default


def show_at_launch(settings=None):
    """Whether a launch should open the welcome window.

    An unset key is a launch that has never seen this window — the first one, which
    is the launch it exists for — so an unset key is a yes. After that the answer is
    whatever the reader left the checkbox saying.
    """
    return _as_bool(_store(settings).value(SETTING_SHOW_AT_LAUNCH, None), True)


def proposed_dont_show(settings=None):
    """What **Don't show this again** is ticked to when the window opens.

    The same key, read with the opposite default, and deliberately so: an unset key
    means the reader is looking at this window for the first time, and having read
    it once is the whole of what it is for — so the tick is already in the box and
    a first launch is the only launch that opens it. Unticking is how someone who
    wants it back gets it back, and Help → Welcome opens it either way.
    """
    stored = _store(settings).value(SETTING_SHOW_AT_LAUNCH, None)
    return not _as_bool(stored, False)


def set_show_at_launch(on, settings=None):
    """Write the preference. Called when the window closes, by any route."""
    _store(settings).setValue(SETTING_SHOW_AT_LAUNCH, bool(on))


class WelcomeDialog(QDialog):
    """The window itself: text, links, version, address, and the preference.

    Everything on it is a label. The links are rich text with
    ``setOpenExternalLinks``, so the browser (or the mail client) is opened by Qt
    from the anchor's own href and there is no slot in between that could point
    somewhere else than the text says.
    """

    def __init__(self, parent=None, settings=None):
        super().__init__(parent)
        # The window's own store when it has one, so a test — or a capture — can
        # hand this dialog a throwaway file and never touch the user's settings.
        if settings is None:
            settings = getattr(parent, "settings", None)
        self._settings = _store(settings)

        self.setWindowTitle("Welcome to XSLOPE Studio")
        lay = QVBoxLayout(self)

        heading = QLabel("<b>Welcome to XSLOPE Studio.</b>")
        heading.setTextFormat(Qt.RichText)
        lay.addWidget(heading)

        lay.addWidget(self._body(
            "Studio analyzes slope stability, seepage and finite element models "
            "from an Excel input file. Use <b>File → Open</b> to load a "
            "problem, or <b>File → New</b> to start an empty one and build it "
            "up with the input editors or the AI assistant."))

        lay.addWidget(self._body(
            "The documentation is the manual, and it is where the sample and "
            "verification problems live — each with the input file it was "
            "solved from:"))

        self.link_label = self._body(
            "<ul>"
            f'<li><a href="{links.DOCS_URL}">Documentation</a> — '
            "the whole site: inputs, methods, and both engines.</li>"
            f'<li><a href="{links.GETTING_STARTED_URL}">Getting started</a> — '
            "installing Studio, and what a first session looks like.</li>"
            f'<li><a href="{links.SAMPLES_URL}">Sample problems</a> — '
            "worked examples, each with its input file to download.</li>"
            f'<li><a href="{links.VERIFICATION_URL}">Verification</a> — '
            "XSLOPE's results against published benchmarks and other software.</li>"
            "</ul>")
        self.link_label.setOpenExternalLinks(True)
        lay.addWidget(self.link_label)

        self.version_label = self._body(self._version_note())
        lay.addWidget(self.version_label)

        self.support_label = self._body(
            "Questions, problems and suggestions: "
            f'<a href="{links.SUPPORT_MAILTO}">{links.SUPPORT_EMAIL}</a>.')
        self.support_label.setOpenExternalLinks(True)
        lay.addWidget(self.support_label)

        self.dont_show = QCheckBox("Don't show this again "
                                   "(Help → Welcome reopens it)")
        self.dont_show.setObjectName("dont_show_again")
        self.dont_show.setChecked(proposed_dont_show(self._settings))
        lay.addWidget(self.dont_show)

        # One button, and it carries Qt's reject role — so the Close button, the
        # Escape key and the title bar are one path, the one :meth:`done` covers.
        buttons = QDialogButtonBox(QDialogButtonBox.Close, self)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

    def _body(self, text):
        """A wrapped rich-text paragraph, one text column wide."""
        label = QLabel(text)
        label.setTextFormat(Qt.RichText)
        label.setWordWrap(True)
        width = label.fontMetrics().horizontalAdvance("n" * TEXT_COLUMNS)
        # Both bounds, so the column is the same width whether the paragraph is a
        # line long or five: a maximum alone lets the shortest paragraph decide how
        # wide the window opens.
        label.setMinimumWidth(width)
        label.setMaximumWidth(width)
        return label

    def _version_note(self):
        """The version this is, and how it learns about a newer one."""
        from xslope import __version__

        return (f"You are running version <b>{__version__}</b>. "
                "<b>Help → Check for Updates…</b> reports a newer "
                "release, and Studio makes that check quietly once a day unless "
                "<b>Help → Check for Updates at Startup</b> is switched off.")

    def done(self, result):
        """Persist the preference however the window was closed.

        ``done`` is the one path every close goes through — the Close button,
        Escape, and the title-bar close — so the preference cannot be lost by
        dismissing the window in a way nobody wired up.
        """
        set_show_at_launch(not self.dont_show.isChecked(), self._settings)
        super().done(result)
