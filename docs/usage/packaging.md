# Project Packaging

## A project is a set of files

An XSLOPE project starts as one Excel workbook, but it does not stay that way. Every
analysis that produces something worth keeping writes it to a **sidecar** file beside
the workbook, named after it:

| File | Written by | Holds |
|------|------------|-------|
| `{base}.xlsx` | you, or Studio | the model: geometry, materials, water, loads |
| `{base}_mesh.json` | Build Mesh | the finite element mesh |
| `{base}_seep.csv` | a steady seepage run | the pore-pressure field |
| `{base}_seep2.csv` | a rapid-drawdown run | the drawn-down field (stage 2) |
| `{base}_tseep.csv`, `{base}_tseep_meta.json` | a transient seepage run | every frame of the march |
| `{base}_fem_*.csv`, `{base}_fem_meta.json` | an FEM or SSRM run | nodal and element results, member forces |

The workbook is loaded by name and the sidecars follow it automatically — opening
`slope1.xlsx` picks up `slope1_mesh.json`, and the pore-pressure field from
`slope1_seep.csv` where the materials ask for one, in Studio and in Python alike. That
is why they are named the way they are, and why a project moved by copying only the
workbook arrives without its results.

## What a package is

A **project package** is a `.xslz` file: a plain zip holding the workbook and all of its
sidecars, and nothing else. It exists so a project can be emailed, uploaded, or handed
to a colleague as one file.

A package is a way for a project to travel, not a place to work in. Nothing edits or
saves inside one: opening a package unpacks it back to loose files first, and the
workbook that comes out is an ordinary `.xlsx` you can also open in Excel. Day-to-day
work never involves a package at all — Save writes loose files, and recent files point
at the loose workbook.

Because the workbook and the results are zipped at one moment, everything in a package
agrees with everything else. A project sent as loose attachments can arrive with a
workbook from Tuesday and a mesh from Monday; a package cannot.

## Creating a package

### Studio

**File → Export Project Package…** collects the current project and writes
`{base}.xslz` beside it. The whole set goes in: the workbook and every sidecar next to
it.

The package is built from the files on disk, so if the project has edits or a solution
this session has not written out yet, Studio offers to save first. Cancelling that
prompt cancels the export — a package whose workbook disagreed with the results zipped
beside it would defeat the point.

### Python

```python
import xslope

xslope.pack("slope1.xlsx")                    # writes slope1.xslz beside it
xslope.pack("slope1.xlsx", dest="~/outbox")   # or into a folder of your choosing
```

`pack` returns the path it wrote. It refuses to replace a package that already exists
unless you pass `overwrite=True`.

To see what would go in without writing anything:

```python
xslope.project_files("slope1.xlsx")
# ['/work/slope1.xlsx', '/work/slope1_mesh.json', '/work/slope1_seep.csv']
```

Sidecars are collected by the naming convention rather than a list of known file types,
so anything named after the project travels — including results a future analysis
writes, and any `{base}_*` file of your own that nothing in xslope recognizes.

Where two workbooks share a folder and one name starts with the other, a file goes to
the workbook that would actually **read** it. With `slope1.xlsx` and
`slope1_drained.xlsx` side by side, `slope1_drained_mesh.json` is the mesh
`slope1_drained` loads, so it packs with `slope1_drained`; but `slope1_fem_nodes.csv`
is the FEM result `slope1` loads, so it packs with `slope1` — even though
`slope1_fem.xlsx` may also be sitting there. (That workbook's own FEM results would be
`slope1_fem_fem_nodes.csv`.) A workbook never travels inside another project's
package.

## Opening a package

Opening a package always unpacks it first. The default destination is a folder named
for the package, beside it — `slope1.xslz` unpacks to `slope1/` — and never a temporary
directory, because the workbook that comes out is one you may want to open in Excel
later.

### Studio

**File → Open** accepts `.xslz` alongside `.xlsx`. Choosing a package shows where its
files will go:

![The Open project package dialog](../studio/images/usage_unpack_package_dialog.png)

**Change…** picks a different folder to unpack into. Once the files are out, Studio
opens the extracted workbook through its normal open path, so the window title, the
recent files list and everything else refer to the loose workbook — the package is not
mentioned again.

If the destination folder already exists, the dialog asks instead of guessing, because
the copy already there may hold your own edits:

- **Open Existing** leaves the folder untouched and opens the project already in it.
- **Extract Fresh** unpacks into a new folder beside it (`slope1-2/`).

### Python

`load_slope_data` accepts a package and does the whole thing — unpack, then load the
workbook that came out:

```python
from xslope.fileio import load_slope_data

slope_data = load_slope_data("slope1.xslz")   # unpacks to slope1/, loads slope1.xlsx
```

To extract without loading:

```python
import xslope

workbook = xslope.unpack("slope1.xslz")       # returns .../slope1/slope1.xlsx
```

A script cannot be asked a question, so the already-exists case is a refusal rather
than a guess. If the destination folder is there, both calls raise, and the message
names the three ways forward:

```
FileExistsError: /work/slope1 already exists, and unpacking slope1.xslz would write
over what is in it. Either open the workbook already there (/work/slope1/slope1.xlsx),
or pass dest= to extract to a different folder, or pass overwrite=True to replace the
files in this one.
```

```python
load_slope_data("slope1.xslz", dest="/work/run2")     # somewhere else
load_slope_data("slope1.xslz", overwrite=True)        # replace what is there
```

`overwrite=True` replaces the files the package holds and leaves anything else in the
folder alone.

There is no way to read a project out of a package without extracting it, by design:
the solvers write their results beside the workbook, and inside a zip there is nowhere
for them to land.

## Samples on this site

Every sample and verification problem in this documentation is published as a project
package, so what you get is the whole model — the workbook plus whatever mesh and
results were saved with it — rather than a workbook that has to be re-solved. The
packages are rebuilt each time the site is built, from the same files the results on
these pages were computed from.

Each input file is offered twice:

- **Download** saves the `.xslz` package. Open it afterwards with File → Open, or
  `load_slope_data`, exactly as above.
- **Open in Studio** hands it straight to XSLOPE Studio. Studio shows you what it is
  about to fetch and from where, downloads the package to your Downloads folder when
  you agree, and then unpacks and opens it the same way File → Open would.

Open in Studio works when Studio was installed from an
[installer](../getting_started/install.md), which is what registers `xslope://` links
with your operating system; it also registers `.xslz`, so a downloaded package opens on
a double-click. A `pip install` registers neither — those links will do nothing, and
Download is the way in.

Studio only fetches from the sites XSLOPE publishes from, refuses anything else by
name, and never opens a file on your computer through a link: a page can put an
`xslope://` link in front of you, but it cannot make Studio act on it without your
say-so.
