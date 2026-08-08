---
title: "Install the xslope Python package — XSLOPE"
description: "Install the xslope Python package for slope stability and seepage analysis: pip install with the gui, fem and ai extras, prerequisites, and installing from source."
---

# Python Package Install

This page covers the `xslope` **Python package** — for scripts, notebooks, and the
API. For the desktop application, download the installer for your platform from
[Install](../getting_started/install.md); it carries its own Python and needs
nothing on this page.

The package installs from PyPI (with `pip`). You can also install directly from
the source repository.

If you would rather not install anything, XSLOPE also runs entirely in your browser
through the pre-configured Google Colab notebooks — see
[Colab Notebooks](notebooks.md).

## Installing with pip

XSLOPE is available on PyPI. If you intend to only run limit equilibrium problems, you can install xslope as follows:

```bash
pip install xslope
```

If you intend to perform either seepage analysis or slope stability using the finite element method, install the
`fem` extra, which adds the **gmsh** meshing package:

```bash
pip install xslope[fem]
```

The graphical desktop application, [XSLOPE Studio](../studio/index.md), is also in
this package, behind the `gui` extra (`pip install "xslope[gui]"`) — the route to
use when you want Studio inside an environment you already manage. Most users
should take the [native installer](../getting_started/install.md) instead. See
[XSLOPE Studio → Installation](../studio/index.md#installation) for the full set of
extras (`gui`, `fem`, `ai`) and for how to launch it.

### Linux and Google Colab

On Debian/Ubuntu Linux (including Google Colab), gmsh needs system OpenGL libraries that
aren't installed by default. Run this once before installing the `fem` extra:

```bash
apt-get update && apt-get install -y libgl1 libglu1-mesa
```

macOS and Windows need no extra step — gmsh ships its own libraries.

## Accessing the functions

After installing xslope, you can access the functions as follows:

```python
import xslope as xslope

from xslope.fileio import load_slope_data, print_dictionary
from xslope.mesh import build_polygons, build_mesh_from_polygons, export_mesh_to_json, import_mesh_from_json
from xslope.plot import plot_inputs, plot_mesh, plot_polygons, plot_polygons_separately
from xslope.plot_seep import plot_seep_data, plot_seep_solution
from xslope.seep import build_seep_data, run_seepage_analysis, save_seep_data_to_json, export_seep_solution
```

## Installing from source

If you wish to access all of the underlying code in xslope, you can bypass PyPI and install directly from the 
repository:

[https://github.com/njones61/xslope/tree/main](https://github.com/njones61/xslope/tree/main)

To clone the repository, click on the green Code button and click the clipboard icon to copy the clone URL.

![clone_repo.png](images/clone_repo.png)

Then, in your terminal, navigate to the directory where you want to clone the repository and run the following 
command (you can either type the URL or paste it from the clipboard):

```bash
git clone https://github.com/njones61/xslope.git
```

You can also clone the repo directly into your IDE.

