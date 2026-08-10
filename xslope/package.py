"""Project packaging — collecting a project into one file, and taking it apart again.

A project is a *set* of files, not one file: the workbook ``{base}.xlsx`` plus the
sidecars that basename-key to it (``{base}_mesh.json``, ``{base}_seep.csv``,
``{base}_seep2.csv``, the ``{base}_fem_*`` result tables, and anything else written
beside the workbook under the same convention). The loose folder is the working
format — the workbook is deliberately plain Excel so anyone can open and edit it
with no special tools, and every solver writes its results next to it.

A ``.xslz`` package is a *transport* format: a plain zip of that set, flat, with the
workbook first and no manifest. It exists so a project can be emailed, uploaded or
linked as one artifact. Nothing lives in it. Opening a package extracts it back to
loose files and then opens those; neither Studio nor the library ever edits inside a
package or saves back into one.

The naming convention is defined once, here, in :func:`project_files`: a project's
sidecars are its ``{base}_*`` siblings. There is no list of known suffixes to keep in
step with the solvers — a sidecar a future solver invents travels automatically.
"""

import os
import zipfile

#: The package extension. Studio's file dialogs and the docs both spell it from here.
PACKAGE_EXT = ".xslz"


def is_package(path):
    """True if ``path`` names a project package (by extension, case-insensitively)."""
    return os.path.splitext(str(path))[1].lower() == PACKAGE_EXT


def _abs(path):
    """An absolute path, with ``~`` expanded — every path in and out of this module
    goes through here, so ``dest="~/outbox"`` means what it says."""
    return os.path.abspath(os.path.expanduser(str(path)))


def project_files(xlsx_path):
    """Return every file that belongs to the project at ``xlsx_path``.

    The workbook comes first, then its sidecars in name order. A sidecar is any
    ``{base}_*`` file sitting beside the workbook — the same basename convention
    ``load_slope_data`` uses to find ``{base}_mesh.json`` and ``{base}_seep.csv``.

    Two kinds of ``{base}_*`` sibling are NOT this project's:

    * another workbook (``{base}_something.xlsx``) — folders of samples routinely
      hold ``xslope_earth_dam1.xlsx`` beside ``xslope_earth_dam1_vg.xlsx``, and
    * that other workbook's own sidecars (``{base}_something_mesh.json``).

    A sidecar therefore belongs to the LONGEST workbook basename that prefixes it,
    which is exactly the attribution the loader makes when it goes looking for a
    companion of one workbook.
    """
    xlsx_path = _abs(xlsx_path)
    if not os.path.isfile(xlsx_path):
        raise FileNotFoundError(f"No workbook at {xlsx_path}")
    folder = os.path.dirname(xlsx_path)
    base = os.path.splitext(os.path.basename(xlsx_path))[0]

    names = sorted(os.listdir(folder))
    # The other workbooks whose names extend this one; each owns its own sidecars.
    others = [os.path.splitext(n)[0] for n in names
              if n.lower().endswith(".xlsx") and not n.startswith("~$")
              and os.path.splitext(n)[0].startswith(base + "_")]

    files = [xlsx_path]
    for name in names:
        if not name.startswith(base + "_"):
            continue
        path = os.path.join(folder, name)
        if not os.path.isfile(path):
            continue
        if name.lower().endswith(".xlsx"):
            continue                        # another project's workbook
        stem = os.path.splitext(name)[0]
        if any(stem.startswith(other + "_") for other in others):
            continue                        # that project's sidecar
        files.append(path)
    return files


def package_path(xlsx_path, dest=None):
    """Where :func:`pack` would write the package for ``xlsx_path``.

    ``dest=None`` puts it beside the workbook, named for it; a directory ``dest``
    puts it in that directory under the same name; anything else is taken as the
    package path itself.
    """
    xlsx_path = _abs(xlsx_path)
    name = os.path.splitext(os.path.basename(xlsx_path))[0] + PACKAGE_EXT
    if dest is None:
        return os.path.join(os.path.dirname(xlsx_path), name)
    dest = _abs(dest)
    if os.path.isdir(dest):
        return os.path.join(dest, name)
    return dest


def pack(xlsx_path, dest=None, overwrite=False):
    """Collect a project into a ``.xslz`` package and return the package path.

    ``xlsx_path`` is the project's workbook; its sidecars come along by the
    :func:`project_files` convention. The archive is flat — no folder inside — with
    the workbook as its first entry, so the workbook is never buried.

    ``dest`` follows :func:`package_path`. An existing package is left alone unless
    ``overwrite=True``.
    """
    files = project_files(xlsx_path)
    out = package_path(xlsx_path, dest)
    if os.path.exists(out) and not overwrite:
        raise FileExistsError(
            f"A package already exists at {out}. Pass a different dest=, or "
            f"overwrite=True to replace it.")
    folder = os.path.dirname(out)
    if folder:
        os.makedirs(folder, exist_ok=True)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in files:
            zf.write(path, arcname=os.path.basename(path))
    return out


def package_contents(package):
    """Return the flat member names in a package, workbook first.

    Raises ValueError if the archive is not a project package: entries must be plain
    file names (no folders, no paths that would escape the destination), and exactly
    one of them must be the ``.xlsx`` workbook.
    """
    package = _abs(package)
    with zipfile.ZipFile(package) as zf:
        names = zf.namelist()
    for name in names:
        if name.endswith("/") or "/" in name or "\\" in name or os.path.isabs(name) \
                or name in (".", "..") or name.startswith(".."):
            raise ValueError(
                f"{os.path.basename(package)} is not a project package: it contains "
                f"{name!r}, and a package holds only plain files.")
    books = [n for n in names if n.lower().endswith(".xlsx")]
    if len(books) != 1:
        raise ValueError(
            f"{os.path.basename(package)} is not a project package: it holds "
            f"{len(books)} .xlsx workbooks, and a package holds exactly one.")
    return books + [n for n in names if n not in books]


def unpack_path(package, dest=None):
    """The folder :func:`unpack` would extract ``package`` into.

    The default is a folder named for the package, beside it (``slope1.xslz`` ->
    ``slope1/``). Never a temp directory: the workbook that comes out is one the user
    is going to open in Excel, so it has to land somewhere durable and visible.
    """
    package = _abs(package)
    if dest is not None:
        return _abs(dest)
    return os.path.join(os.path.dirname(package),
                        os.path.splitext(os.path.basename(package))[0])


def unpack(package, dest=None, overwrite=False):
    """Extract a project package to loose files and return the workbook's path.

    ``dest`` defaults to :func:`unpack_path` — a folder named for the package,
    beside it. If that folder already exists this raises, because the copy already
    there may hold the user's edits; ``overwrite=True`` writes the package's files
    over the ones of the same name (anything else in the folder is left alone).

    The package is never read in place: solvers write their results beside the
    workbook, and inside a zip there is nowhere for them to land.
    """
    package = _abs(package)
    names = package_contents(package)
    out = unpack_path(package, dest)
    if os.path.exists(out) and not overwrite:
        raise FileExistsError(
            f"{out} already exists, and unpacking {os.path.basename(package)} "
            f"would write over what is in it. Either open the workbook already "
            f"there ({os.path.join(out, names[0])}), or pass dest= to extract to a "
            f"different folder, or pass overwrite=True to replace the files in "
            f"this one.")
    os.makedirs(out, exist_ok=True)
    with zipfile.ZipFile(package) as zf:
        for name in names:
            with zf.open(name) as src, open(os.path.join(out, name), "wb") as dst:
                dst.write(src.read())
    return os.path.join(out, names[0])
