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
sidecars are its ``{base}_*`` siblings. Collection is a glob, not a list of known
suffixes, so a sidecar a future solver invents travels automatically.
"""

import os
import zipfile

#: The package extension. Studio's file dialogs and the docs both spell it from here.
PACKAGE_EXT = ".xslz"

# === THE SIDECAR SUFFIXES, IN ONE PLACE =====================================
# Every name a solver writes beside a workbook, and every name a loader goes looking
# for. This is the ONE definition: xslope.fileio reads MESH_SIDECAR / SEEP_SIDECARS
# to find its companions, studio.main_window reads FEM_SOLUTION_SIDECARS to write and
# delete an FEM solution, and project_files below uses the whole set to decide which
# of two workbooks in one folder a file belongs to.
#
# The set does NOT decide what travels — an unrecognized {base}_* sibling still packs
# (see project_files). It decides only ATTRIBUTION, and it has to be exact: a
# nearly-right list here would hand one project's results to another.

#: The mesh, read by ``load_slope_data`` and written by Build Mesh.
MESH_SIDECAR = "_mesh.json"
#: Steady pore-pressure fields, per boundary-condition set (2 = rapid drawdown).
SEEP_SIDECARS = ("_seep.csv", "_seep2.csv")
#: A transient march: the frames and the run metadata (written by xslope.seep).
TSEEP_SIDECARS = ("_tseep.csv", "_tseep_meta.json")
#: Every file ``xslope.fem.export_fem_solution`` writes beside a model, by suffix —
#: the converged fields, the at-failure snapshot, the run metadata, and the per-member
#: force tables for reinforcement and piles, converged and at-failure. Deleting a
#: solution means deleting ALL of them: the member tables are read back by name, so a
#: list naming only the nodal files leaves the last run's bar forces on disk to be
#: grafted onto whatever solution is imported next.
FEM_SOLUTION_SIDECARS = (
    "_fem_nodes.csv", "_fem_elements.csv", "_fem_meta.json",
    "_fem_reinf.csv", "_fem_piles.csv",
    "_fem_failure_nodes.csv", "_fem_failure_elements.csv",
    "_fem_failure_meta.json",
    "_fem_failure_reinf.csv", "_fem_failure_piles.csv",
)
#: Studio's per-project display styles (studio.document).
STYLES_SIDECAR = "_styles.json"

#: The union: every suffix that makes a file some workbook's own companion.
SIDECAR_SUFFIXES = ((MESH_SIDECAR,) + SEEP_SIDECARS + TSEEP_SIDECARS
                    + FEM_SOLUTION_SIDECARS + (STYLES_SIDECAR,))


def is_package(path):
    """True if ``path`` names a project package (by extension, case-insensitively)."""
    return os.path.splitext(str(path))[1].lower() == PACKAGE_EXT


def _abs(path):
    """An absolute path, with ``~`` expanded — every path in and out of this module
    goes through here, so ``dest="~/outbox"`` means what it says."""
    return os.path.abspath(os.path.expanduser(str(path)))


def companion_of(name, stems):
    """Which of ``stems`` claims the file ``name``, or None.

    Claimed means what the LOADERS mean by it: the file is exactly some workbook's
    stem plus one of the suffixes those loaders read and those writers write
    (:data:`SIDECAR_SUFFIXES`). Nothing is inferred from the shape of the name, and
    a file no suffix matches is claimed by nobody.

    No suffix is a tail of another at an underscore boundary, so at most one stem can
    match and the answer does not depend on the order of ``stems``.
    """
    for stem in stems:
        if not name.startswith(stem):
            continue
        for suffix in SIDECAR_SUFFIXES:
            if name == stem + suffix:
                return stem
    return None


def project_files(xlsx_path):
    """Return every file that belongs to the project at ``xlsx_path``.

    The workbook comes first, then its sidecars in name order. A sidecar is any
    ``{base}_*`` file sitting beside the workbook — the same basename convention
    ``load_slope_data`` uses to find ``{base}_mesh.json`` and ``{base}_seep.csv``.
    A ``{base}_*`` sibling nothing recognizes travels too: the package is transport,
    and a file named after the project is part of it until some other project claims
    it.

    Two kinds of ``{base}_*`` sibling are NOT this project's:

    * another workbook (``{base}_something.xlsx``), and
    * any other workbook's own companions — the files ITS loaders would read, which
      means its stem plus one of :data:`SIDECAR_SUFFIXES`, and nothing else.

    That second test is the loaders' test, and it cannot be replaced by "the longest
    workbook name that prefixes the file". ``vp091.xlsx`` and ``vp091_fem.xlsx`` sit
    in one corpus folder, and ``vp091``'s own FEM results are ``vp091_fem_nodes.csv``,
    ``vp091_fem_meta.json`` and their at-failure twins — every one of which the
    longest-prefix reading hands to ``vp091_fem``, whose real companions are
    ``vp091_fem_fem_nodes.csv`` and friends. Under that reading ``vp091`` packs as a
    bare workbook and its entire solved FEM set is dropped on the floor.

    The claim is checked against every other workbook in the folder, not only those
    whose names extend this one: read the same pair the other way round and it is
    ``vp091_fem`` that would otherwise walk off with ``vp091``'s results.
    """
    xlsx_path = _abs(xlsx_path)
    if not os.path.isfile(xlsx_path):
        raise FileNotFoundError(f"No workbook at {xlsx_path}")
    folder = os.path.dirname(xlsx_path)
    base = os.path.splitext(os.path.basename(xlsx_path))[0]

    names = sorted(os.listdir(folder))
    others = [os.path.splitext(n)[0] for n in names
              if n.lower().endswith(".xlsx") and not n.startswith("~$")
              and os.path.splitext(n)[0] != base]

    files = [xlsx_path]
    for name in names:
        if not name.startswith(base + "_"):
            continue
        path = os.path.join(folder, name)
        if not os.path.isfile(path):
            continue
        if name.lower().endswith(".xlsx"):
            continue                        # another project's workbook
        # This project's own companion wins any tie by construction — its loader is
        # going to read the file by that exact name.
        if companion_of(name, [base]) is None and companion_of(name, others):
            continue                        # another project's companion
        files.append(path)
    return files


def package_path(xlsx_path, dest=None):
    """Where :func:`pack` would write the package for ``xlsx_path``.

    ``dest=None`` puts it beside the workbook, named for it. A ``dest`` ending in
    ``.xslz`` is the package path itself; anything else is a folder to put it in,
    under the same name — whether or not that folder exists yet, so a destination
    that is meant as a folder cannot silently become an extensionless file.
    """
    xlsx_path = _abs(xlsx_path)
    name = os.path.splitext(os.path.basename(xlsx_path))[0] + PACKAGE_EXT
    if dest is None:
        return os.path.join(os.path.dirname(xlsx_path), name)
    dest = _abs(dest)
    if is_package(dest):
        return dest
    return os.path.join(dest, name)


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
