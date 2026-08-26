---
title: "Download XSLOPE Studio for macOS and Windows"
description: "Download and install XSLOPE Studio, the free slope stability and seepage analysis application for macOS and Windows — system requirements, first launch, updates, and the pip alternative for Python users."
---

# Install

**XSLOPE Studio** is a desktop application for slope stability and seepage
analysis. Download it, open it, and start working — there is nothing else to
install and no Python to set up.

<div class="download-buttons" markdown="1">

[Download for macOS](https://github.com/njones61/xslope/releases/latest/download/XSLOPE-Studio-macos-arm64.dmg){ .btn .btn-neutral .download-btn }

[Download for Windows](https://github.com/njones61/xslope/releases/latest/download/XSLOPE-Studio-windows-x64-setup.exe){ .btn .btn-neutral .download-btn }

</div>

Each button downloads the newest release for that platform directly — the macOS
disk image (`.dmg`, Apple silicon) or the Windows installer (`.exe`, 64-bit).
The [releases page](https://github.com/njones61/xslope/releases/latest){target="blank"}
lists the same files under **Assets** along with the `.sha256` checksum beside
each one and `latest.json`, the manifest Studio's own updater reads; neither of
those needs to be downloaded by hand.

![The XSLOPE Studio main window](../studio/images/studio_main_window.png){width="1200"}

## System requirements

| Platform | Requirement |
| --- | --- |
| **macOS** | macOS 11 (Big Sur) or later, on Apple silicon (M1 or newer). Intel Macs are not covered by the installer — install the [Python package](#for-python-users) instead. |
| **Windows** | Windows 10 or 11, 64-bit (x64). |

Each download is a few hundred megabytes, and the installed application is larger
still: it carries a complete Python runtime, the scientific stack, the **gmsh**
mesher, and the compiled fast kernel, so nothing else has to be installed. A
network connection is needed to download and to check for updates; analyses
themselves run entirely offline, and only the
[AI assistant](../studio/assistant.md) contacts a service while you work.

## First launch

### macOS

Open the downloaded `.dmg` and drag **XSLOPE Studio** onto the **Applications**
shortcut in the same window. Eject the disk image, then launch Studio from
Applications, Launchpad, or Spotlight. The first launch takes a few seconds
longer than later ones while macOS verifies the application.

### Windows

Run the downloaded `…-setup.exe` and follow the installer. It installs for the
current user only, so Windows never asks for administrator rights, and the last
page offers to launch Studio immediately. Afterwards, Studio is on the Start menu
under **XSLOPE Studio**.

<details markdown="1">
<summary>If the download is not code-signed</summary>

A build that is not code-signed behaves like any other unsigned download, and the
operating system warns you once before it runs:

- **Windows** — SmartScreen shows a blue *Windows protected your PC* dialog.
  Click **More info**, then **Run anyway** to start the installer.
- **macOS** — Gatekeeper refuses a double-click. Right-click (or Control-click)
  **XSLOPE Studio** in Applications, choose **Open**, and confirm once. Later
  launches open normally.

The `.sha256` file published beside each installer lets you confirm that what you
downloaded is what was built.

</details>

Once the window is open, use **File → Open** to load an Excel problem, or
**File → New** to start an empty project and build it up with the editors or the
[AI assistant](../studio/assistant.md). The [XSLOPE Studio](../studio/index.md)
section is the tour: the [interface](../studio/interface.md),
[editing inputs](../studio/editing.md), and
[running analyses](../studio/analysis.md).

![The File menu](../studio/images/overview_file_menu.png)

## What gets installed

Studio is self-contained: the Python runtime, the analysis engine, gmsh, the
compiled fast kernel, the Excel input template, and a set of sample projects all
live inside the application. Nothing is added to a system Python, and an existing
`pip` or `conda` installation of `xslope` is untouched.

| Platform | Where it lands |
| --- | --- |
| **macOS** | `XSLOPE Studio.app` in `/Applications` — or wherever you drag it. |
| **Windows** | `%LOCALAPPDATA%\Programs\XSLOPE Studio`, plus a Start-menu shortcut and an entry in **Apps & features**. Everything is per user: no administrator rights and no system-wide changes. |

Your own settings live outside the application, so they survive updates and
reinstalls: window layout, recent files, and assistant preferences in the OS
settings store (`HKCU\Software\XSlope` on Windows,
`~/Library/Preferences/com.XSlope.XSlope Studio.plist` on macOS), and assistant
API keys in the system keychain under the service **XSlope Studio**.

## Keeping it up to date

**Help → Check for Updates…** asks whether a newer release exists, and
**Help → Check for Updates at Startup** — on by default — repeats that check
quietly at most once a day, saying nothing unless there is something to say.

When an update is available, Studio offers **Download & Install**: it downloads
the installer for your platform, verifies its SHA-256 checksum against the one
published with the release, and discards anything that does not match. On Windows
the verified installer then runs silently and the new version restarts itself; on
macOS Studio opens the disk image and asks you to drag the new copy into
Applications. See [App management](../studio/app_management.md) for the full
behavior, including what the dialogs say and how to turn the startup check off.

## Uninstalling

- **macOS** — drag `XSLOPE Studio.app` from Applications to the Trash.
- **Windows** — **Settings → Apps → XSLOPE Studio → Uninstall**, or the
  **Uninstall XSLOPE Studio** shortcut in the Start-menu folder. Program files,
  the shortcut, and the Apps & features entry are removed.

Uninstalling leaves your projects alone, and it leaves the settings and keychain
entries listed above in place, so a later reinstall starts where you left off.
Remove them by hand if you want a clean slate.

## Python package

XSLOPE is also a Python package. The engine Studio runs is the `xslope` package
documented under [API](../api/solve.md), and it reads and writes the same Excel
format, so files move freely between Studio, scripts, and notebooks. It installs
from PyPI:

```bash
pip install xslope                # limit equilibrium only
pip install "xslope[fem]"         # add gmsh, for seepage and FEM
pip install "xslope[gui]"         # add Studio itself
pip install "xslope[gui,fem,ai]"  # everything, including the AI assistant
```

The `gui` extra registers a console command that opens the same window the
installers do:

```bash
xslope-studio
```

Installed this way, Studio reports updates but never changes your environment:
the update dialog shows the `pip install -U "xslope[gui]"` line for you to run
yourself. This is also how Studio runs on Linux and on Intel Macs, which the
installers do not cover.

**Linux and Google Colab.** On Debian/Ubuntu (including Colab), gmsh needs system
OpenGL libraries that are not installed by default. Run this once before
installing the `fem` extra:

```bash
apt-get update && apt-get install -y libgl1 libglu1-mesa
```

macOS and Windows need no extra step — gmsh ships its own libraries. To run XSLOPE
in a browser with nothing installed at all, see
[Colab Notebooks](../usage/notebooks.md).

**Using the package.** After installing, import what you need:

```python
from xslope.fileio import load_slope_data
from xslope.slice import generate_slices
from xslope.solve import spencer
from xslope.plot import plot_inputs, plot_solution
```

**Installing from source.** To work with the code itself, bypass PyPI and clone
the repository at
[github.com/njones61/xslope](https://github.com/njones61/xslope/tree/main):

```bash
git clone https://github.com/njones61/xslope.git
```
