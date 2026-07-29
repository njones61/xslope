# App Management

This page covers installing, updating, and uninstalling Studio as a native desktop
application (`.dmg` / `.msi`). Those installers are not yet available, so the page is
a stub for now; install and manage Studio with `pip` instead — see
[Installation](index.md#installation).

## Planned coverage

When the native installers are released, this page will cover:

- **Installing** — running the macOS `.dmg` (drag to Applications) and the Windows
  `.msi` installer, and the first-launch security prompts (Gatekeeper /
  SmartScreen) for a signed app.
- **Updating** — how to move to a newer release, and whether updates are manual or
  checked automatically.
- **Uninstalling** — removing the app cleanly on each platform.
- **Where files live** — the location of settings, recent-files history, API keys
  (OS keychain), and any cached data, on macOS and Windows.
- **System requirements** — supported OS versions and disk footprint (the bundle
  includes Python, Qt, and the scientific stack).

## Managing the pip install

In the meantime, the `pip`-installed Studio is managed like any Python package:

```bash
pip install --upgrade "xslope[gui]"     # update
pip uninstall xslope                    # remove
```

Settings and credentials for the pip install are stored by the OS (Qt `QSettings`
and the system keychain) and persist across upgrades.
