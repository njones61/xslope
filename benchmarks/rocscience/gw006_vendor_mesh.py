"""Same-mesh cross-check for GW6 case 4: solve the vendor's own model in XSLOPE.

The corpus file ``gw006d.xlsx`` states the problem the way XSLOPE takes it —
geometry, material, boundary polylines — and XSLOPE meshes it afresh.  That
leaves two things varying at once against the vendor's published answer: the
model transcription AND the discretization.  This script removes the second.  It
reads the vendor's RS2 model, rebuilds ITS mesh and ITS boundary cards node for
node, and solves that with XSLOPE's own seepage solver, so the remaining
difference is the physics.

The vendor archive is LOCAL and is never committed; neither is anything this
script reads.  What is committed is the script and whatever finding it supports.

    ~/python_projects/vendor_files/rocscience_downloads/
        RS2_Groundwater-Analysis-Verification.zip
          └ Groundwater Analysis Verification/groundwater #006_04.fez   (a zip)
              ├ groundwater #006_04.fea    mesh: nodes + cst_element (tri3) table
              └ groundwater #006_04.slw    the seepage model: restraints, tractions

Override the archive path with ``$XSLOPE_GW_VENDOR_ZIP`` or pass it as the first
argument.  With the archive absent the script prints a note and exits 0, so it
never breaks a checkout that does not have it.

What is read, and how
---------------------
``.fea`` — ``nodes:`` gives ``<id> x: <x> y: <y>``; ``elements:`` gives
``<id> type: cst_element nodes: [a,b,c] material: <name> materialID: <n>``.  Only
``cst_element`` (the linear triangle) appears in this case, which is XSLOPE's
``tri3``.  Node ids are remapped to a 0-based contiguous index.

``.slw`` — the same node numbering.  ``restraints:`` carries the nodal SPECIFIED
TOTAL HEADS as ``n: <id>  tx: <head>``: ``tx: 10`` are the submerged upstream
face (the reservoir) and ``tx: 0`` are the base drain.  ``tractions:`` carries
the fourteen infiltration cards as ``e: <elem> side: <k> ... qn1: <q> qn2: <q>``,
where side *k* spans local nodes *k*−1 and *k* (0-based), cyclically: side 1 is
the edge between local nodes 0 and 1, side 2 between 1 and 2, side 3 between 2
and 0.  That mapping is not guessed — it is the only one of the three rotations
under which all fourteen cards land on the dam's exposed surface, and
``check_traction_edges`` re-derives it on every run rather than trusting it.

How the drain is posed
----------------------
The vendor writes the drain as a specified head of 0.  XSLOPE is given it as an
EXIT FACE instead, which is what ``gw006d.xlsx`` does and why: with no exit face
anywhere the model has no free surface to track, XSLOPE selects the confined
solver, and the unsaturated law this problem exists to exercise is dropped.  The
exit face is placed on exactly the nodes the vendor holds at head 0, so the two
statements cover the same boundary; only the condition on it differs, and it
differs in the direction that keeps the physics.

Run:  PYTHONPATH=. python3 benchmarks/rocscience/gw006_vendor_mesh.py [ZIP]
"""

import io
import os
import re
import sys
import zipfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np                                                  # noqa: E402

DEFAULT_ZIP = os.environ.get(
    "XSLOPE_GW_VENDOR_ZIP",
    os.path.expanduser("~/python_projects/vendor_files/rocscience_downloads/"
                       "RS2_Groundwater-Analysis-Verification.zip"))
MEMBER = "Groundwater Analysis Verification/groundwater #006_04.fez"
CORPUS = os.path.join(os.path.dirname(__file__), '..', '..', 'docs',
                      'verification', 'files', 'rocscience_gw', 'gw006d.xlsx')

#: The crest centerline — the manual's line 1-1 — and the elevations Fig 6.18
#: prints markers at.
LINE_X = 26.0
LINE_Y = (0.05, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0)
#: The element size the GW6 test tags mesh gw006d at.
CORPUS_SIZE = 1.0

#: 1-based traction side k -> the pair of 0-based local node indices it spans.
SIDE_NODES = {1: (0, 1), 2: (1, 2), 3: (2, 0)}

#: The dam's exposed surface, as the polyline every traction edge must lie on:
#: upstream toe, crest shoulders, downstream toe (Fig 6.1).
SURFACE = ((0.0, 0.0), (24.0, 12.0), (28.0, 12.0), (52.0, 0.0))


def read_fez(zip_path=DEFAULT_ZIP, member=MEMBER):
    """The ``.fea`` and ``.slw`` text out of the nested archive."""
    with zipfile.ZipFile(zip_path) as outer:
        fez = zipfile.ZipFile(io.BytesIO(outer.read(member)))
        text = {}
        for name in fez.namelist():
            if name.endswith((".fea", ".slw")):
                text[name[-3:]] = fez.read(name).decode("utf-8", "replace")
    return text["fea"], text["slw"]


def _node_table(lines, header):
    """``{id: (x, y)}`` from the block that follows ``header``."""
    start = next(i for i, l in enumerate(lines) if l.startswith(header))
    out = {}
    for line in lines[start + 1:]:
        m = re.match(r"\s*(\d+)\s+x:\s*([-\d.eE+]+)\s+y:\s*([-\d.eE+]+)", line)
        if m:
            out[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
        elif out and not line.strip():
            break
    return out


def parse_mesh(fea):
    """The vendor mesh as XSLOPE takes it: nodes, tri3 elements, types, materials.

    Element node ids are remapped to a 0-based contiguous index in ascending
    vendor-id order, and the element rows are padded to nine columns, which is the
    shape ``build_mesh_from_polygons`` produces for a mixed-order mesh.
    """
    lines = fea.splitlines()
    nid_xy = _node_table(lines, "nodes:")
    ids = sorted(nid_xy)
    id2idx = {nid: k for k, nid in enumerate(ids)}
    nodes = np.array([nid_xy[nid] for nid in ids], dtype=float)

    start = next(i for i, l in enumerate(lines) if l.startswith("elements:"))
    row = re.compile(r"\s*(\d+)\s+type:\s*(\S+)\s+nodes:\s*\[([\d,\s]+)\]"
                     r"\s+material:\s*\S+\s+materialID:\s*(\d+)")
    elems, mats, eid_row = [], [], {}
    for line in lines[start + 1:]:
        m = row.match(line)
        if not m:
            if elems and not line.strip():
                break
            continue
        if m.group(2) != "cst_element":
            raise ValueError(f"unexpected element type {m.group(2)!r}; this "
                             "converter handles the linear triangle only")
        tri = [id2idx[int(v)] for v in m.group(3).split(",")]
        eid_row[int(m.group(1))] = len(elems)
        elems.append(tri + [0] * 6)
        mats.append(int(m.group(4)))
    return {"nodes": nodes,
            "elements": np.array(elems, dtype=int),
            "element_types": np.full(len(elems), 3, dtype=int),
            "element_materials": np.array(mats, dtype=int),
            "id2idx": id2idx, "eid_row": eid_row}


def parse_bcs(slw, mesh):
    """The vendor's nodal heads and edge tractions, on the mesh's 0-based indexing.

    Returns ``(heads, tractions)`` where ``heads`` is ``{node_index: total head}``
    and ``tractions`` is a list of ``(element_row, (node_a, node_b), qn)``.
    """
    lines = slw.splitlines()
    id2idx, eid_row = mesh["id2idx"], mesh["eid_row"]

    start = next(i for i, l in enumerate(lines) if l.strip() == "restraints:")
    heads = {}
    for line in lines[start + 1:]:
        m = re.match(r"n:\s*(\d+)\s+tx:\s*([-\d.eE+]+)\s*$", line.strip())
        if m:
            heads[id2idx[int(m.group(1))]] = float(m.group(2))
        elif heads and not line.strip():
            break

    start = next(i for i, l in enumerate(lines) if l.strip() == "tractions:")
    tractions = []
    for line in lines[start + 1:]:
        m = re.match(r"e:\s*(\d+)\s+side:\s*(\d+).*?qn1:\s*([-\d.eE+]+)"
                     r".*?qn2:\s*([-\d.eE+]+)", line.strip())
        if m:
            if m.group(3) != m.group(4):
                raise ValueError(f"traction on element {m.group(1)} is not "
                                 "uniform along its edge; this converter takes "
                                 "a uniform qn only")
            r = eid_row[int(m.group(1))]
            a, b = SIDE_NODES[int(m.group(2))]
            tractions.append((r, (int(mesh["elements"][r][a]),
                                  int(mesh["elements"][r][b])),
                              float(m.group(3))))
        elif tractions and not line.strip():
            break
    return heads, tractions


def check_traction_edges(mesh, slw, tol=1e-6):
    """Re-derive the traction side convention instead of trusting it.

    A traction card names an element and a side; which pair of local nodes that
    side spans is a file-format convention, and getting it wrong puts the rain on
    an interior edge where it is silently wrong rather than obviously wrong.  All
    three rotations are tried and the one that puts EVERY card on the dam's
    exposed surface is the answer.  Raises if that is not exactly one of them.
    """
    from shapely.geometry import LineString, Point

    surface = LineString(SURFACE)
    lines = slw.splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip() == "tractions:")
    cards = []
    for line in lines[start + 1:]:
        m = re.match(r"e:\s*(\d+)\s+side:\s*(\d+)", line.strip())
        if m:
            cards.append((mesh["eid_row"][int(m.group(1))], int(m.group(2))))
        elif cards and not line.strip():
            break

    good = []
    for shift in range(3):
        mapping = {k: ((k - 1 + shift) % 3, (k + shift) % 3) for k in (1, 2, 3)}
        on = all(
            all(surface.distance(Point(mesh["nodes"][mesh["elements"][r][i]]))
                <= tol for i in mapping[side])
            for r, side in cards)
        if on:
            good.append(mapping)
    if len(good) != 1:
        raise ValueError(f"{len(good)} of the 3 side conventions put every "
                         "traction card on the dam surface; the file format "
                         "cannot be read unambiguously")
    if good[0] != SIDE_NODES:
        raise ValueError(f"traction side convention is {good[0]}, not the "
                         f"{SIDE_NODES} this module assumes")
    return len(cards)


def vendor_flux_nodal(mesh, tractions):
    """The vendor's tractions as a consistent nodal load vector.

    ``q·L/2`` at each end of a linear edge — the same rule
    ``xslope.seep.assemble_flux_nodal`` applies, computed here straight from the
    vendor's own cards so the two can be compared rather than assumed equal.
    """
    nodes = mesh["nodes"]
    f = np.zeros(len(nodes))
    for _r, (a, b), qn in tractions:
        length = float(np.hypot(*(nodes[b] - nodes[a])))
        f[a] += qn * length / 2.0
        f[b] += qn * length / 2.0
    return f


def build_vendor_seep_data(mesh, heads, tractions, drain_head=0.0):
    """A ``seep_data`` on the vendor mesh carrying the vendor's own conditions.

    Built through ``build_seep_data`` so the material properties, the assembly and
    every solver-facing field come from XSLOPE's normal path; only the boundary
    arrays are then replaced by the vendor's, node for node.  The nodes the vendor
    holds at ``drain_head`` become an exit face (see the module docstring).
    """
    from xslope.fileio import load_slope_data
    from xslope.seep import build_seep_data

    sd = load_slope_data(CORPUS)
    m = {k: mesh[k] for k in ("nodes", "elements", "element_types",
                              "element_materials")}
    seep_data = build_seep_data(m, sd, check_inputs=False)

    bc_type = np.zeros(len(mesh["nodes"]), dtype=int)
    bc_values = np.zeros(len(mesh["nodes"]))
    for idx, h in heads.items():
        if h == drain_head:
            bc_type[idx] = 2                      # exit face (the toe drain)
        else:
            bc_type[idx] = 1
            bc_values[idx] = h
    # What XSLOPE's own polyline path made of the corpus file's boundary
    # statement on this mesh, kept so the caller can compare it with the
    # vendor's cards: that comparison IS the transcription check.
    ours = {"bc_type": np.array(seep_data["bc_type"]),
            "flux_nodal": np.array(seep_data["flux_nodal"], dtype=float)}
    seep_data["bc_type"] = bc_type
    seep_data["bc_values"] = bc_values
    seep_data["flux_nodal"] = vendor_flux_nodal(mesh, tractions)
    return sd, seep_data, ours


def solve_corpus_fresh(size=None):
    """The same case on XSLOPE's own mesh, at the discretization the GW6 tags run
    at — the column the vendor-mesh run is compared against."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)
    from xslope.seep import build_seep_data, run_seepage_analysis

    sd = load_slope_data(CORPUS)
    mesh = build_mesh_from_polygons(get_material_polygons(sd),
                                    size or CORPUS_SIZE, "tri3",
                                    size_regions=extract_size_regions(sd))
    seep_data = build_seep_data(mesh, sd)
    solution = run_seepage_analysis(seep_data, tol=1e-5, max_iter=2000)
    head = np.asarray(solution["head"], dtype=float)
    return solution["flowrate"], line_profile(np.asarray(seep_data["nodes"]),
                                              head)


def line_profile(nodes, head, x=LINE_X, elevations=LINE_Y):
    """Pressure head up the vertical line at ``x``, by linear interpolation over
    the mesh triangulation."""
    import matplotlib.tri as mtri

    triang = mtri.Triangulation(nodes[:, 0], nodes[:, 1])
    interp = mtri.LinearTriInterpolator(triang, head)
    out = []
    for y in elevations:
        h = np.ma.filled(interp(x, y), np.nan)
        out.append(float(h) - y)
    return out


def main(argv):
    zip_path = argv[1] if len(argv) > 1 else DEFAULT_ZIP
    if not os.path.exists(zip_path):
        print(f"vendor archive not present at {zip_path} — nothing to do "
              "(this check is local-only; set $XSLOPE_GW_VENDOR_ZIP to point at it)")
        return 0

    from xslope.seep import run_seepage_analysis

    fea, slw = read_fez(zip_path)
    mesh = parse_mesh(fea)
    n_cards = check_traction_edges(mesh, slw)
    print(f"side convention verified against the dam surface on all "
          f"{n_cards} traction cards")
    heads, tractions = parse_bcs(slw, mesh)
    nodes = mesh["nodes"]
    print(f"vendor mesh   {len(nodes)} nodes · {len(mesh['elements'])} tri3 "
          f"elements · materials {sorted(set(mesh['element_materials']))}")
    fixed = sorted({h for h in heads.values()})
    print(f"vendor heads  {len(heads)} nodal cards at total head {fixed}")
    for h in fixed:
        xs = [nodes[i][0] for i, v in heads.items() if v == h]
        ys = [nodes[i][1] for i, v in heads.items() if v == h]
        print(f"   head {h:g}: {len(xs)} nodes, x {min(xs):g}–{max(xs):g}, "
              f"y {min(ys):g}–{max(ys):g}")
    print(f"vendor rain   {len(tractions)} edge cards at qn "
          f"{sorted({q for _r, _e, q in tractions})}")
    total = sum(q * float(np.hypot(*(nodes[b] - nodes[a])))
                for _r, (a, b), q in tractions)
    span = sorted({round(nodes[i][0], 6) for _r, e in
                   [(r, e) for r, e, _q in tractions] for i in e})
    print(f"   Σ q·L = {total:.6e} over x {span[0]:g}–{span[-1]:g}")

    sd, seep_data, ours = build_vendor_seep_data(mesh, heads, tractions)

    # The transcription check: what gw006d.xlsx's boundary POLYLINES produce on
    # this same mesh, against the vendor's own node and edge cards.
    theirs_type = np.asarray(seep_data["bc_type"])
    same_head = int(np.sum((ours["bc_type"] == 1) & (theirs_type == 1)))
    same_face = int(np.sum((ours["bc_type"] == 2) & (theirs_type == 2)))
    print(f"\ntranscription: gw006d's polylines on this mesh give "
          f"{int(np.sum(ours['bc_type'] == 1))} head and "
          f"{int(np.sum(ours['bc_type'] == 2))} exit-face nodes; the vendor's "
          f"cards give {int(np.sum(theirs_type == 1))} and "
          f"{int(np.sum(theirs_type == 2))} — {same_head} and {same_face} shared")
    print(f"   assembled rain: ours {ours['flux_nodal'].sum():.6e} vs vendor "
          f"cards {np.asarray(seep_data['flux_nodal']).sum():.6e} · worst nodal "
          f"difference {np.abs(ours['flux_nodal'] - np.asarray(seep_data['flux_nodal'])).max():.3e}")

    solution = run_seepage_analysis(seep_data, tol=1e-5, max_iter=2000)
    head = np.asarray(solution["head"], dtype=float)
    print(f"\nXSLOPE on the vendor mesh: Q = {solution['flowrate']:.6e} "
          f"m³/s per m · converged {solution.get('converged')}")

    # The same model, meshed by XSLOPE at the corpus tag's discretization, so the
    # only thing separating the two columns is the mesh.
    fresh_q, fresh_psi = solve_corpus_fresh()
    print(f"XSLOPE on its own {CORPUS_SIZE:g} m mesh:  Q = {fresh_q:.6e} "
          f"(vendor mesh is {100.0 * (solution['flowrate'] - fresh_q) / fresh_q:+.2f}%)")
    print("   line 1-1 (x = %g) pressure head, vendor mesh vs fresh mesh:"
          % LINE_X)
    print("      %-6s %10s %10s %8s" % ("y", "vendor", "fresh", "Δ"))
    for y, a, b in zip(LINE_Y, line_profile(nodes, head), fresh_psi):
        print("      %-6g %10.4f %10.4f %+8.4f" % (y, a, b, a - b))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
