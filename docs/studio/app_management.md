# App Management

!!! warning "Placeholder — pending native installers (Phase 7)"
    This page covers installing, updating, and uninstalling Studio as a **native
    desktop application** (`.dmg` / `.msi`). Those installers are not yet
    available, so this page is a stub. It will be written when the installers
    ship. For now, install and manage Studio with `pip` — see
    [Installation](index.md#installation).

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

In the meantime, the `pip`-installed Studio is managed like any Python package:

```bash
pip install --upgrade "xslope[gui]"     # update
pip uninstall xslope                    # remove
```

Settings and credentials for the pip install are stored by the OS (Qt `QSettings`
and the system keychain) and persist across upgrades.
