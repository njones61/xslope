"""W-2, Importing Models: the three imports in Studio.

Walked: Part 1 — File → Import DXF… on the page's drawing, the wizard's layer
mapping and the three material names typed in, the Inputs tree counts and the
note the import gives, the material rows pasted from the page's table, the
surcharge set in the distributed loads editor, and the Spencer search. Part 2 —
File → Import GeoStudio (SLOPE/W)… with the Spencer analysis picked, the notes,
and the single-surface run on SLOPE/W's circle. Part 3 — File → Import Slide2…
with No Seismic picked, the starting circle added in the circles editor, and the
search. Parts 2 and 3 open vendor files the page does not redistribute; where
they are not on this machine (``~/python_projects/vendor_files``) those parts
are noted as not walked rather than failed.
"""
import os
import re

from tutorial_walks.harness import (StepFailed, column, editable, page_rows, paste,
                                    press, repo_path, rows, set_field, set_value,
                                    show_parameters_for)

PAGE = "w02_import.md"
TIMEOUT = 1800
VENDOR = os.path.join(os.path.expanduser("~"), "python_projects", "vendor_files")
GSZ = os.path.join(VENDOR, "gsz_corpus", "Baker and Leschinsky - Earth Dam.gsz")
SLMD_NAME = "Tutorial 28 Seismic.slmd"


def _find_slmd():
    for root, _dirs, files in os.walk(os.path.join(VENDOR, "rocscience_downloads")):
        if SLMD_NAME in files:
            return os.path.join(root, SLMD_NAME)
    return None


def _choose(path):
    """Answer the next File → Open dialog with ``path`` (a copy in the work folder)."""
    from PySide6.QtWidgets import QFileDialog
    real = QFileDialog.getOpenFileName

    def pick(*a, **k):
        QFileDialog.getOpenFileName = real
        return path, ""
    QFileDialog.getOpenFileName = staticmethod(pick)


def _spencer(analysis):
    def act(dlg):
        set_field(dlg, "Method", "Spencer")
        set_field(dlg, "Analysis", analysis)
    return act


def _notes(st, mark):
    return " ".join((st.log_since(mark) + " " + " ".join(t for _k, _ti, t in st.dialogs.boxes)
                     + " " + st.status()).split())


def walk(w):
    p = w.page
    st = w.studio()
    import shutil

    with w.step("Part 1: File → Import DXF"):
        dxf = os.path.join(w.workdir, "w02_section.dxf")
        shutil.copy2(repo_path(w.page_link("w02_section.dxf")), dxf)
        names = {"EMBANKMENT": "embankment", "SILTY_CLAY": "silty clay", "DENSE_SAND": "dense sand"}

        def wizard(dlg):
            t = dlg.table
            hdr = [t.horizontalHeaderItem(c).text() for c in range(t.columnCount())]
            for want in ("Layer", "Contents", "Import as", "Material"):
                w.require(any(h.startswith(want) for h in hdr), "the wizard has no %s column (%s)" % (want, hdr))
            c_mat = next(c for c, h in enumerate(hdr) if h.startswith("Material"))
            seen = set()
            for r in range(t.rowCount()):
                layer = t.item(r, 0).text()
                if layer in names:
                    t.item(r, c_mat).setText(names[layer])
                    seen.add(layer)
            w.require(seen == set(names), "the wizard does not list the three material layers: %s" % seen)
        _choose(dxf)
        mark = st.mark()
        st.trigger("File", "Import DXF", wizard)
        counts = {k: st.tree_count(k) for k in ("Materials", "Polygons", "Circles",
                                                "Piezometric lines", "Distributed loads")}
        want = {"Materials": "3", "Polygons": "3", "Circles": "2", "Piezometric lines": "3",
                "Distributed loads": "1"}
        w.check(counts == want, "the Inputs tree counts %s, not the page's %s" % (counts, want))
        note = p.find(r"one note comes back: \*(.+?)\*\. There is no reinforcement").group(1)
        w.check(" ".join(note.split()) in _notes(st, mark), "the import does not note %r" % note)
        w.check([m["name"] for m in st.data["materials"]] == list(names.values()),
                "the materials came in as %s" % [m["name"] for m in st.data["materials"]])

    with w.step("Part 1: the material rows and the surcharge"):
        tbl = p.table(r"the columns none of the three soils use stay blank:")

        def mats(dlg):
            press(dlg, "Table view")
            show_parameters_for(dlg, {"LEM", "FEM"})
            t = editable(dlg)
            paste(t, page_rows(tbl), 0, column(t, "name"))
        st.edit("Materials", mats)
        for m, r in zip(st.data["materials"], tbl[1:]):
            got = (m["name"], float(m["gamma"]), float(m["gamma_sat"]), m["option"], float(m["c"]),
                   float(m["phi"]), m["u"])
            want = (r[0], float(r[1]), float(r[2]), r[3], float(r[4]), float(r[5]), r[13])
            w.require(got == want, "a material saved as %s, not %s" % (got, want))
        q = p.number(r"set the surcharge to (\d+) psf at both")[0]

        def dl(dlg):
            t = editable(dlg)
            c = next(i for i, h in enumerate(_hdr(t)) if h.startswith("Normal"))
            for r in range(t.table.rowCount()):
                t.table.item(r, _cols(t)[c]).setText(q)
        st.edit("Distributed loads", dl)
        loads = st.data.get("dloads") or []
        w.require(loads and all(abs(float(pt["Normal"]) - float(q)) < 1e-9 for pt in loads[0]),
                  "the surcharge saved as %s" % loads)

    with w.step("Part 1: the Spencer search"):
        st.run_lem(_spencer("Auto search"))
        b = st.results["lem_solution"]
        w.expect("DXF model FS", b["results"]["FS"], r"The search returns \*\*FS = (1\.\d+)\*\* on a circle tangent")
        crit = b["search"]["fs_cache"][0]
        top_sand = -15.0
        w.check(abs(crit["Depth"] - top_sand) < 0.05,
                "the circle is not tangent to the top of the dense sand (lowest point %.2f)" % crit["Depth"])

    if not os.path.exists(GSZ):
        w.note("Part 2 not walked: %s is not on this machine" % GSZ)
    else:
        with w.step("Part 2: File → Import GeoStudio"):
            gsz = os.path.join(w.workdir, os.path.basename(GSZ))
            shutil.copy2(GSZ, gsz)

            def pick(dlg):
                t = dlg.table
                names_ = [t.item(r, 0).text() for r in range(t.rowCount())]
                w.require(names_ and names_[0].startswith("Spencer"), "the first analysis is %s" % names_[:1])
                t.selectRow(0)
            _choose(gsz)
            mark = st.mark()
            st.trigger("File", "Import GeoStudio (SLOPE/W)", pick)
            notes = _notes(st, mark)
            slw = p.number(r"SLOPE/W's factor of safety on it, (\d\.\d+) by Spencer")[0]
            w.check(slw in notes, "the notes do not give SLOPE/W's %s" % slw)
            w.check(abs(float(st.data.get("tcrack_depth") or 0) - 5.0) < 1e-9,
                    "the tension crack imported at %s" % st.data.get("tcrack_depth"))
            w.check(not (st.data.get("dloads") or []), "the loads sheet is not empty")
            st.run_lem(_spencer("Single surface"))
            w.expect("GeoStudio surface FS", st.results["lem_solution"]["results"]["FS"],
                     r"XSLOPE returns \*\*FS = (1\.\d+)\*\* against SLOPE/W's")

    slmd = _find_slmd()
    if slmd is None:
        w.note("Part 3 not walked: %s is not on this machine" % SLMD_NAME)
        return
    with w.step("Part 3: File → Import Slide2"):
        src = os.path.join(w.workdir, SLMD_NAME)
        shutil.copy2(slmd, src)
        seen = {}

        def pick(dlg):
            t = dlg.table
            seen["open"] = t.item(t.currentRow(), 0).text() if t.currentRow() >= 0 else None
            r = next(i for i in range(t.rowCount()) if t.item(i, 0).text() == "No Seismic")
            t.selectRow(r)
        _choose(src)
        st.trigger("File", "Import Slide2", pick)
        w.check(seen["open"] == "Master Scenario", "the picker opens on %r" % seen["open"])
        w.check(not st.data.get("circles"), "a failure surface was imported")
        m = p.find(r"\*\*Xo\*\* = (\d+), \*\*Yo\*\* = (\d+), \*\*Option\*\* = `(\w+)`, "
                   r"\*\*Xi\*\* = (\d+), \*\*Yi\*\* = (\d+)")
        vals = dict(zip(("Xo", "Yo", "Option", "Xi", "Yi"), m.groups()))

        def circ(dlg):
            press(dlg, "Add row")
            t = editable(dlg)
            r = t.table.rowCount() - 1
            for k, v in vals.items():
                c = column(t, k)
                cw = t.table.cellWidget(r, c)
                if cw is not None:
                    from PySide6.QtWidgets import QComboBox
                    set_value(cw if isinstance(cw, QComboBox) else cw.findChild(QComboBox), v)
                else:
                    t.table.item(r, c).setText(v)
        st.edit("Circles", circ)
        st.run_lem(_spencer("Auto search"))
        w.expect("Slide2 model FS", st.results["lem_solution"]["results"]["FS"],
                 r"The search returns \*\*FS = (1\.\d+)\*\* against Slide2's")


def _hdr(t):
    from tutorial_walks.harness import headers
    return headers(t)[0]


def _cols(t):
    from tutorial_walks.harness import headers
    return headers(t)[1]
