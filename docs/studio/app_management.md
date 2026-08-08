# App Management

This page covers keeping Studio up to date. Installing it, where it puts its
files, and uninstalling it are on [Install](../getting_started/install.md).

## Checking for updates

**Help → Check for Updates…** asks whether a newer XSLOPE Studio has been
released. The check reads a small version manifest attached to the newest
[GitHub release](https://github.com/njones61/xslope/releases) — it never scrapes a
web page — and answers with one of three dialogs: you are up to date, version *N*
is available (with a link to its release notes), or the check could not reach the
network.

Studio also checks on its own, quietly. **Help → Check for Updates at Startup** is
on by default and runs at most once a day, a few seconds after the window opens.
It says nothing at all unless a newer version exists, in which case it raises the
same notification the menu item does — a window you can dismiss and keep working
behind. A failed automatic check is silent; only a check you asked for reports
that it failed. Turn the preference off and Studio never contacts the network on
its own.

What the update offers depends on how Studio was installed:

| Installed with | What the update dialog offers |
| --- | --- |
| **pip** | The literal `pip install -U "xslope[gui]"` line, with a **Copy** button. Studio never modifies your Python environment itself — you run the upgrade where you installed it, and restart Studio. |
| **A native installer** (`.dmg` / `.exe`) | **Download & Install**. Studio downloads the installer for your platform with a progress bar you can cancel, and verifies its SHA-256 checksum against the one published with the release. A file that does not match is discarded and nothing is installed. |

On **Windows**, a verified installer runs silently: Studio closes, the installer
upgrades the existing per-user installation in place, and the new version starts
by itself. On **macOS**, Studio mounts the downloaded disk image and asks you to
drag **XSLOPE Studio** to Applications, replacing the old copy, then quit and
relaunch.

If a release states that it cannot be installed on top of your version — or
publishes no download for your platform — the dialog says so and sends you to the
releases page to install it fresh instead.

## Managing the pip install

The `pip`-installed Studio is managed like any Python package:

```bash
pip install --upgrade "xslope[gui]"     # update
pip uninstall xslope                    # remove
```

Settings and credentials are stored by the OS (Qt `QSettings` and the system
keychain) and persist across upgrades — the same store an installed Studio uses,
so the two installations share their preferences on one machine.

## Managing a native install

Where the application and its settings live, and how to remove them on each
platform, are covered under
[Install → What gets installed](../getting_started/install.md#what-gets-installed)
and [Install → Uninstalling](../getting_started/install.md#uninstalling).
