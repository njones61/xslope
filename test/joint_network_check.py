"""The joint NETWORK generators, against counts that can be worked out by hand.

``xslope.joints`` turns a set description — a dip, a spacing, a persistence, a
block size — into the individual joint lines the ``joints`` sheet holds. The
conversion is geometric and every reading below is a number a reader can check
against the drawing rather than against the code.

Legs:

  1. **Counts and spacing.** A parallel set in a 20 x 10 m box at a stated dip
     and spacing produces the traces the geometry allows, spaced by exactly the
     stated perpendicular distance and all at exactly the stated dip. Halving the
     spacing roughly doubles the count; the offset moves the set without changing
     what it is.
  2. **Clipping.** A set generated for one material stays inside that material's
     polygon; a set generated for the whole section stays inside the section; and
     a trace that would lie ALONG the region's boundary is dropped, because the
     mesh split needs material on both sides of a joint. The horizontal set in a
     box whose top and bottom are boundaries is the case that reads it.
  3. **Persistence.** A set given ``(trace_len, gap)`` comes out in pieces of the
     stated length separated by the stated gap, and the pieces of one trace lie
     on one line.
  4. **Voronoi.** The same seed gives the same network, a different seed a
     different one; the cell size comes out near the stated block size; and no
     wall lies on the region's boundary.
  5. **Refusals.** Two sets at the same dip overlap trace for trace and
     ``cross_jointed`` refuses them by name; two sets sharing a label are
     refused; an unknown material name is refused with the model's own material
     names in the message.
  6. **The mesh.** A cross-jointed block generated here meshes, its joints meet,
     and the split gives every junction node one copy per wedge of material
     around it — four at a crossing. This is the leg that says a generated
     network is the same input a typed one is.
  7. **The round trip.** A generated network written to the ``joints`` sheet and
     read back is the same network, endpoint for endpoint and property for
     property.
  8. **The set record.** Every kind's parameters, region and elevation band go
     into one row's Label cell and come back off it as the same record and the
     same text, so a set can be reopened and regenerated from the file alone.
  9. **Regeneration.** Editing a set replaces its own rows where they were and
     touches nothing else — not the other set, not a hand-entered joint line, not
     the properties the rows carry.
  10. **The band and the joint region.** A set cut to an elevation band keeps the
     traces on the band's own edges, where the section's boundary drops them; a
     set confined to a polygon of Type ``joints`` resolves the same by the
     region's name and by its number.
  11. **The display name.** A generated row's Label carries the set record, and
     no reader wants to see it: the 1D details list, the member figures, the
     report's interface table, preflight's refusals and the joints editor's list
     all show the ``bed-03`` head of it, while the file keeps the record each
     row was generated from.

Run directly:  PYTHONPATH=. python3 test/joint_network_check.py
"""

import contextlib
import io
import math
import os
import sys
import tempfile
import warnings

warnings.filterwarnings('ignore')

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from shapely.geometry import LineString, Polygon

from xslope.joints import (JointSet, cross_jointed, parallel_set, regenerate,
                           display_label, resolve_region, set_name, sets_in,
                           voronoi)

BOX = Polygon([(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)])


def _model():
    """A two-material section: a lower unit 0-4 m and an upper one 4-10 m."""
    lower = Polygon([(0.0, 0.0), (20.0, 0.0), (20.0, 4.0), (0.0, 4.0)])
    upper = Polygon([(0.0, 4.0), (20.0, 4.0), (20.0, 10.0), (0.0, 10.0)])
    return {
        'domain_polygon': BOX,
        'polygons': [{'polygon': lower, 'mat_id': 0},
                     {'polygon': upper, 'mat_id': 1}],
        'materials': [{'name': 'Sandstone'}, {'name': 'Shale'}],
    }


def _line(row):
    return LineString([(row['x1'], row['y1']), (row['x2'], row['y2'])])


def _dip_of(row):
    a = math.degrees(math.atan2(row['y2'] - row['y1'], row['x2'] - row['x1']))
    return (a + 180.0) % 180.0          # a trace has no sense, only a direction


# --------------------------------------------------------------------------
# 1 — counts, spacing and dip
# --------------------------------------------------------------------------

def _leg_counts(failures, results):
    sd = _model()
    for dip, spacing in ((0.0, 2.0), (30.0, 2.0), (-60.0, 2.5), (90.0, 4.0)):
        rows = parallel_set(sd, dip, spacing, label='s', props={'phi': 30.0})
        if not rows:
            failures.append(f"counts: a set at dip {dip:g}, spacing {spacing:g} "
                            f"produced no trace in a 20 x 10 box")
            continue
        want = (dip + 180.0) % 180.0
        worst = max(abs(_dip_of(r) - want) for r in rows)
        if worst > 1e-9:
            failures.append(f"counts: a trace of the dip {dip:g} set comes out "
                            f"at {worst:.2e} degrees off its own dip")
        # The perpendicular offsets have to be an arithmetic sequence at the
        # stated spacing. Each trace's own offset along the set normal:
        th = math.radians(dip)
        nx, ny = -math.sin(th), math.cos(th)
        offs = sorted({round(r['x1'] * nx + r['y1'] * ny, 9) for r in rows})
        gaps = {round(b - a, 6) for a, b in zip(offs, offs[1:])}
        if gaps and gaps != {round(spacing, 6)}:
            failures.append(f"counts: the dip {dip:g} set's traces are offset by "
                            f"{sorted(gaps)}, not by the stated spacing "
                            f"{spacing:g}")
        results.append(f"counts  dip {dip:g}, spacing {spacing:g}: "
                       f"{len(rows)} traces, offsets an exact {spacing:g} apart")

    n2 = len(parallel_set(sd, 30.0, 2.0, label='s', props={'phi': 30.0}))
    n1 = len(parallel_set(sd, 30.0, 1.0, label='s', props={'phi': 30.0}))
    if not (1.7 * n2 <= n1 <= 2.3 * n2 + 2):
        failures.append(f"counts: halving the spacing takes the trace count "
                        f"from {n2} to {n1}, which is not about double")
    results.append(f"counts  halving the spacing: {n2} traces -> {n1}")

    # The offset moves the set, and moving it by a whole spacing puts it back.
    base = parallel_set(sd, 30.0, 2.0, offset=0.0, label='s', props={'phi': 30.0})
    same = parallel_set(sd, 30.0, 2.0, offset=2.0, label='s', props={'phi': 30.0})
    if len(base) != len(same) or any(
            abs(a['x1'] - b['x1']) > 1e-9 or abs(a['y1'] - b['y1']) > 1e-9
            for a, b in zip(base, same)):
        failures.append("counts: offsetting a set by one whole spacing does not "
                        "reproduce the same set")
    results.append("counts  an offset of one whole spacing reproduces the set")

    labels = [r['label'] for r in base]
    if labels != [f"s-{i + 1:02d}" for i in range(len(base))]:
        failures.append(f"counts: the rows are labeled {labels[:3]}…, not "
                        f"s-01, s-02, …")
    if any(r.get('phi') != 30.0 for r in base):
        failures.append("counts: the stated phi did not reach every row")


# --------------------------------------------------------------------------
# 2 — clipping to the region, and the boundary rule
# --------------------------------------------------------------------------

def _leg_clipping(failures, results):
    sd = _model()
    lower = resolve_region(sd, 'Sandstone')

    rows = parallel_set(sd, 20.0, 1.0, region='Sandstone', label='bed',
                        props={'phi': 28.0})
    outside = [r for r in rows if not lower.buffer(1e-9).contains(_line(r))]
    if outside:
        failures.append(f"clipping: {len(outside)} of {len(rows)} traces "
                        f"generated for the Sandstone leave its polygon")
    results.append(f"clipping  {len(rows)} traces in the Sandstone, all inside "
                   f"its polygon")

    whole = parallel_set(sd, 20.0, 1.0, label='bed', props={'phi': 28.0})
    out2 = [r for r in whole if not BOX.buffer(1e-9).contains(_line(r))]
    if out2:
        failures.append(f"clipping: {len(out2)} traces of the whole-section set "
                        f"leave the section")
    if len(whole) <= len(rows):
        failures.append(f"clipping: the whole section takes {len(whole)} traces "
                        f"and its lower unit alone {len(rows)}; confining a set "
                        f"to one unit cannot produce as many")
    results.append(f"clipping  whole section {len(whole)} traces vs "
                   f"{len(rows)} confined to the lower unit")

    # A horizontal set in a box: the traces at y = 0 and y = 10 would lie ALONG
    # the section's own boundary, which is the one thing the mesher refuses, and
    # they must not be emitted.
    horiz = parallel_set(sd, 0.0, 2.0, label='h', props={'phi': 30.0})
    ys = sorted({round(r['y1'], 6) for r in horiz})
    if 0.0 in ys or 10.0 in ys:
        failures.append(f"clipping: the horizontal set emits a trace on the "
                        f"section boundary (y in {ys})")
    if ys != [2.0, 4.0, 6.0, 8.0]:
        failures.append(f"clipping: the horizontal set at spacing 2 in a 10 m "
                        f"box comes out at y = {ys}, not [2, 4, 6, 8]")
    results.append(f"clipping  horizontal set at spacing 2: y = {ys}, the two "
                   f"on the boundary dropped")

    # An endpoint that lands on a corner of the region must land on it EXACTLY.
    # Clipping is floating-point arithmetic and leaves an endpoint a part in
    # 10^15 off the corner it belongs on; that part is a sliver of boundary edge
    # between the two, and gmsh's edge recovery cannot mesh one — it splits the
    # offending edges and tries again without end. The 45-degree set through a
    # square's corners is the case, and before the endpoints were snapped onto
    # the region's own vertices it hung the mesher.
    sq = Polygon([(0.0, 0.0), (12.0, 0.0), (12.0, 12.0), (0.0, 12.0)])
    box = {'domain_polygon': sq,
           'polygons': [{'polygon': sq, 'mat_id': 0}],
           'materials': [{'name': 'Rock'}]}
    diag = parallel_set(box, 45.0, 4.0, label='d', props={'phi': 30.0})
    corners = [(0.0, 0.0), (12.0, 0.0), (12.0, 12.0), (0.0, 12.0)]
    exact = 0
    for r in diag:
        for (x, y) in ((r['x1'], r['y1']), (r['x2'], r['y2'])):
            for (cx, cy) in corners:
                d = math.hypot(x - cx, y - cy)
                if d == 0.0:
                    exact += 1
                elif d < 1e-6:
                    failures.append(f"clipping: {r['label']} has an endpoint "
                                    f"({x!r}, {y!r}) sitting {d:.3g} from the "
                                    f"corner ({cx:g}, {cy:g}) it belongs on; a "
                                    f"sliver that size is what the mesher "
                                    f"cannot recover")
    if exact == 0:
        failures.append("clipping: the 45-degree set through a square's corners "
                        "produced no endpoint ON a corner, so the snap is not "
                        "under test here")
    results.append(f"clipping  the 45-degree set lands {exact} endpoint(s) "
                   f"exactly on a corner of the region, none near one")


# --------------------------------------------------------------------------
# 3 — persistence
# --------------------------------------------------------------------------

def _leg_persistence(failures, results):
    sd = _model()
    trace_len, gap = 3.0, 1.0
    rows = parallel_set(sd, 0.0, 2.0, persistence=(trace_len, gap), label='p',
                        props={'phi': 30.0})
    full = parallel_set(sd, 0.0, 2.0, label='p', props={'phi': 30.0})
    if len(rows) <= len(full):
        failures.append(f"persistence: a discontinuous set produced {len(rows)} "
                        f"pieces against {len(full)} continuous traces; it has "
                        f"to produce more")
    lengths = sorted({round(_line(r).length, 6) for r in rows})
    long_ones = [L for L in lengths if L > trace_len + 1e-9]
    if long_ones:
        failures.append(f"persistence: pieces longer than the stated trace "
                        f"length {trace_len:g} came out: {long_ones}")
    # Each 20 m trace becomes 3 + 1 repeating: 3, 3, 3, 3, 3, 1 — five full
    # pieces and a 1 m tail.
    per_line = {}
    for r in rows:
        per_line.setdefault(round(r['y1'], 6), []).append(r)
    for y, pieces in sorted(per_line.items()):
        xs = sorted((min(r['x1'], r['x2']), max(r['x1'], r['x2']))
                    for r in pieces)
        for (a1, b1), (a2, b2) in zip(xs, xs[1:]):
            g = a2 - b1
            if abs(g - gap) > 1e-6:
                failures.append(f"persistence: the gap between two pieces of the "
                                f"trace at y = {y:g} is {g:g}, not {gap:g}")
                break
    results.append(f"persistence  {len(full)} continuous traces become "
                   f"{len(rows)} pieces at ({trace_len:g}, {gap:g}), lengths "
                   f"{lengths}")


# --------------------------------------------------------------------------
# 4 — Voronoi
# --------------------------------------------------------------------------

def _leg_voronoi(failures, results):
    sd = _model()
    a = voronoi(sd, 2.0, seed=11, label='v', props={'phi': 25.0})
    b = voronoi(sd, 2.0, seed=11, label='v', props={'phi': 25.0})
    c = voronoi(sd, 2.0, seed=12, label='v', props={'phi': 25.0})
    if a != b:
        failures.append("voronoi: the same block size and seed produced two "
                        "different networks")
    if a == c:
        failures.append("voronoi: two different seeds produced the same network")
    out = [r for r in a if not BOX.buffer(1e-9).contains(_line(r))]
    if out:
        failures.append(f"voronoi: {len(out)} of {len(a)} cell walls leave the "
                        f"section")
    on_edge = [r for r in a
               if BOX.exterior.distance(_line(r).interpolate(0.5, normalized=True))
               < 1e-9]
    if on_edge:
        failures.append(f"voronoi: {len(on_edge)} cell walls lie on the "
                        f"section's own boundary")
    mean_len = sum(_line(r).length for r in a) / len(a)
    if not (0.3 * 2.0 <= mean_len <= 1.5 * 2.0):
        failures.append(f"voronoi: the mean wall length is {mean_len:.3f} on a "
                        f"block size of 2.0, which is not a tessellation at "
                        f"that size")
    finer = voronoi(sd, 1.0, seed=11, label='v', props={'phi': 25.0})
    if len(finer) <= len(a):
        failures.append(f"voronoi: halving the block size takes the wall count "
                        f"from {len(a)} to {len(finer)}")
    results.append(f"voronoi  block 2.0 m: {len(a)} walls, mean length "
                   f"{mean_len:.3f} m; block 1.0 m: {len(finer)} walls; the "
                   f"seed reproduces the network")


# --------------------------------------------------------------------------
# 5 — the refusals
# --------------------------------------------------------------------------

def _leg_refusals(failures, results):
    sd = _model()
    a = parallel_set(sd, 30.0, 2.0, label='a', props={'phi': 30.0})
    b = parallel_set(sd, 30.0, 2.0, label='b', props={'phi': 30.0})
    try:
        cross_jointed(a, b)
    except ValueError as exc:
        if 'lie on one another' not in str(exc):
            failures.append(f"refusals: two sets at the same dip are refused, "
                            f"but the message does not say why: {exc}")
        results.append("refusals  two sets at one dip: refused, named")
    else:
        failures.append("refusals: two sets at the same dip were accepted; "
                        "every trace of one lies on a trace of the other")

    same_label = parallel_set(sd, -30.0, 2.0, label='a', props={'phi': 30.0})
    try:
        cross_jointed(a, same_label)
    except ValueError as exc:
        if 'label' not in str(exc):
            failures.append(f"refusals: two sets sharing a label are refused "
                            f"for the wrong reason: {exc}")
        results.append("refusals  two sets sharing a label: refused, named")
    else:
        failures.append("refusals: two sets carrying the same labels were "
                        "accepted, so their rows cannot be told apart")

    try:
        parallel_set(sd, 0.0, 1.0, region='Limestone', props={'phi': 30.0})
    except ValueError as exc:
        if 'Sandstone' not in str(exc):
            failures.append(f"refusals: an unknown material is refused without "
                            f"naming the model's own materials: {exc}")
        results.append("refusals  an unknown material: refused, the model's own "
                       "names given")
    else:
        failures.append("refusals: a set was generated for a material the model "
                        "does not have")

    for bad, why in ((dict(spacing=0.0), 'spacing'), (dict(spacing=-1.0), 'spacing')):
        try:
            parallel_set(sd, 0.0, bad['spacing'], props={'phi': 30.0})
        except ValueError:
            pass
        else:
            failures.append(f"refusals: a {why} of {bad['spacing']:g} was "
                            f"accepted")

    try:
        parallel_set(sd, 0.0, 2.0, props={'cohesion': 1.0})
    except ValueError as exc:
        if 'not a joint property' not in str(exc):
            failures.append(f"refusals: a misspelt property is refused for the "
                            f"wrong reason: {exc}")
        results.append("refusals  a misspelt property: refused, the sheet's own "
                       "columns listed")
    else:
        failures.append("refusals: a property the joints sheet has no column "
                        "for was accepted and silently dropped")


# --------------------------------------------------------------------------
# 6 — the mesh of a generated cross-jointed block
# --------------------------------------------------------------------------

def _leg_mesh(failures, results):
    from xslope.mesh import build_mesh_from_polygons

    block = Polygon([(0.0, 0.0), (12.0, 0.0), (12.0, 12.0), (0.0, 12.0)])
    sd = {'domain_polygon': block,
          'polygons': [{'polygon': block, 'mat_id': 0}],
          'materials': [{'name': 'Rock'}]}
    s1 = parallel_set(sd, 45.0, 4.0, label='j1', props={'phi': 30.0})
    s2 = parallel_set(sd, -45.0, 4.0, label='j2', props={'phi': 30.0})
    net = cross_jointed(s1, s2)
    if len(net) != len(s1) + len(s2):
        failures.append("mesh: the cross-jointed network is not the two sets")

    sd['joint_lines'] = net
    lines = [[(r['x1'], r['y1']), (r['x2'], r['y2'])] for r in net]
    polys = [{'coords': list(block.exterior.coords), 'mat_id': 0}]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            mesh = build_mesh_from_polygons(
                polys, target_size=1.0, element_type='tri6', lines=lines,
                element_size_1d=1.0,
                joint_lines={i: {'bar': False} for i in range(len(lines))})
    except Exception as exc:
        failures.append(f"mesh: a generated cross-jointed network does not "
                        f"mesh: {exc}")
        return

    ej = mesh.get('elements_joint')
    n_joint = 0 if ej is None else len(ej)
    if n_joint == 0:
        failures.append("mesh: the generated network meshed but produced no "
                        "joint elements, so the mesh was not split along it")
        return

    # At a crossing of two joint lines the split gives the shared node FOUR
    # copies, one per wedge of material. Count the stations whose node appears
    # in four distinct copies: those are the crossings the two sets make inside
    # the block, and there has to be at least one.
    from collections import Counter
    copies = Counter()
    for rec in (mesh.get('joints') or []):
        for st in rec.get('stations', ()):
            for nid in st:
                if int(nid) >= 0:
                    copies[int(nid)] += 1
    if not copies:
        failures.append("mesh: the split recorded no stations")
        return
    results.append(f"mesh  {len(net)} generated joint lines -> "
                   f"{len(mesh['nodes'])} nodes, {n_joint} joint elements")


# --------------------------------------------------------------------------
# 7 — the round trip through the joints sheet
# --------------------------------------------------------------------------

def _leg_roundtrip(failures, results):
    from xslope.fileio import load_slope_data, save_slope_data_to_xlsx

    base = os.path.join(_ROOT, 'docs', 'fem', 'files', 'xslope_griffiths1.xlsx')
    sd = load_slope_data(base)
    net = parallel_set(sd, 25.0, 3.0, label='bed',
                       props={'phi': 32.0, 'c': 5.0, 'c_res': 1.0,
                              'phi_res': 24.0, 'dil': 6.0, 'kn': 2.0e5,
                              'ks': 8.0e4, 'jred': 'No'})
    if not net:
        failures.append("round trip: the generator produced no trace on the "
                        "shipped sample")
        return
    sd['joint_lines'] = net
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, 'network.xlsx')
        with contextlib.redirect_stdout(io.StringIO()):
            save_slope_data_to_xlsx(sd, out)
            back = load_slope_data(out)
    got = back.get('joint_lines') or []
    if len(got) != len(net):
        failures.append(f"round trip: {len(net)} generated joint lines came "
                        f"back as {len(got)}")
        return
    for a, b in zip(net, got):
        for key in ('label', 'x1', 'y1', 'x2', 'y2', 'c', 'phi', 'c_res',
                    'phi_res', 'dil', 't_cut', 'kn', 'ks', 'jred'):
            va, vb = a.get(key), b.get(key)
            if isinstance(va, float) and isinstance(vb, float):
                if math.isnan(va) and math.isnan(vb):
                    continue
                if abs(va - vb) > 1e-6:
                    failures.append(f"round trip: {a['label']} {key} went out "
                                    f"{va!r} and came back {vb!r}")
            elif key == 'jred':
                # The loader normalizes a yes/no choice to lower case, so the
                # round trip preserves the ANSWER, not the capitalization.
                if str(va).strip().lower() != str(vb).strip().lower():
                    failures.append(f"round trip: {a['label']} jred went out "
                                    f"{va!r} and came back {vb!r}")
            elif va != vb:
                failures.append(f"round trip: {a['label']} {key} went out "
                                f"{va!r} and came back {vb!r}")
    results.append(f"round trip  {len(net)} generated lines written to the "
                   f"joints sheet and read back unchanged, residuals and "
                   f"dilation included")



def _leg_record(failures, results):
    """Leg 8: the set record through the label and back.

    Every kind's parameters, its region and its elevation band are written into
    one row label and read back off it. What is checked is both directions: the
    record that comes back EQUALS the one that went out, and the label it writes
    is character for character the label it was read from — a grammar that
    round-trips the object but not the text would leave two spellings of the same
    set in one file.
    """
    cases = [
        ('a parallel set',
         JointSet('bed', 'parallel', {'dip': 30.0, 'spacing': 2.0},
                  region='Sandstone', band=(40.0, None)),
         'bed-03|par|dip=30|s=2|reg=mat:Sandstone|band=40:'),
        ('a parallel set with an offset and a persistence',
         JointSet('j1', 'parallel',
                  {'dip': -60.0, 'spacing': 2.5, 'offset': 1.25,
                   'trace_len': 3.0, 'gap': 1.0}),
         'j1-03|par|dip=-60|s=2.5|off=1.25|len=3|gap=1'),
        ('a cross-jointed set',
         JointSet('x', 'cross',
                  {'dip': 60.0, 'dip2': -60.0, 'spacing': 2.0,
                   'spacing2': 3.0, 'offset2': 1.5},
                  region='poly:North block'),
         'x-03|crs|dip=60,-60|s=2,3|off=0,1.5|reg=poly:North block'),
        ('a Voronoi set',
         JointSet('vor', 'voronoi', {'block_size': 1.5, 'seed': 7},
                  region=['Sandstone', 'Shale']),
         'vor-03|vor|blk=1.5|seed=7|reg=mat:Sandstone+Shale'),
        ('a band below an elevation, with a negative end',
         JointSet('b', 'parallel', {'dip': 0.0, 'spacing': 1.0},
                  band=(-10.0, 5.0)),
         'b-03|par|dip=0|s=1|band=-10:5'),
    ]
    for what, jset, expect in cases:
        label = jset.to_label(3)
        if label != expect:
            failures.append(f"record: {what} wrote {label!r}, expected {expect!r}")
            continue
        back, index = JointSet.from_label(label)
        if index != 3:
            failures.append(f"record: {what} came back as row {index}, not 3")
        if back != jset:
            failures.append(f"record: {what} did not survive its own label "
                            f"({back.params} / {back.region} / {back.band})")
        if back.to_label(3) != label:
            failures.append(f"record: {what} rewrote {label!r} as "
                            f"{back.to_label(3)!r}")

    # The row's own name, and the set it belongs to, read off the label.
    lbl = 'bed-07|par|dip=30|s=2'
    if display_label(lbl) != 'bed-07' or set_name(lbl) != 'bed':
        failures.append(f"record: {lbl!r} reads as {display_label(lbl)!r} / "
                        f"{set_name(lbl)!r}")
    # A hand-typed label is not a set, however it is spelled, and regenerate must
    # never touch one.
    for typed in ('base joint', 'bed-07', '', 'block 3'):
        if set_name(typed) != '':
            failures.append(f"record: the typed label {typed!r} was read as the "
                            f"set {set_name(typed)!r}")

    # A region that cannot be written down is refused when the label is written,
    # by name, rather than producing a set nobody can regenerate.
    try:
        JointSet('p', 'parallel', {'dip': 0.0, 'spacing': 1.0},
                 region=Polygon([(0, 0), (5, 0), (5, 5)])).to_label(1)
        failures.append("record: a bare-polygon region was written into a label")
    except ValueError as exc:
        if 'Type' not in str(exc):
            failures.append(f"record: the bare-polygon refusal does not say what "
                            f"to draw instead: {exc}")

    # A label that is not a record at all, and a field that is not a field.
    for bad, why in (('bed-01', 'no record'),
                     ('bed-01|zzz|dip=0', 'an unknown kind'),
                     ('bed-01|par|tilt=30', 'an unknown field'),
                     ('bed-01|par|band=40', 'a band with no colon')):
        try:
            JointSet.from_label(bad)
            failures.append(f"record: {bad!r} ({why}) was read as a set")
        except ValueError:
            pass
    results.append("record      5 set records written into a row label and read "
                   "back identical, both object and text; a bare-polygon region, "
                   "an unknown kind, an unknown field and a colon-less band all "
                   "refused by name")


def _leg_regeneration(failures, results):
    """Leg 9: regenerating one set in place.

    The point of the record is that a set can be edited. What that has to mean is
    surgical: the set's own rows are replaced where they were, the OTHER set is
    untouched, a hand-entered joint line is untouched, and the properties the rows
    carried come back on the new ones — the label carries the geometry, the
    columns carry the strength.
    """
    model = _model()
    bed = JointSet('bed', 'parallel', {'dip': 0.0, 'spacing': 2.0},
                   props={'phi': 32.0, 'c': 10.0})
    steep = JointSet('st', 'parallel', {'dip': 70.0, 'spacing': 4.0},
                     props={'phi': 28.0})
    typed = {'label': 'base joint', 'x1': 0.0, 'y1': 0.5, 'x2': 20.0,
             'y2': 0.5, 'c': 0.0, 'phi': 20.0}
    model['joint_lines'] = (bed.generate(model) + [dict(typed)]
                            + steep.generate(model))
    n_bed = sum(1 for r in model['joint_lines'] if set_name(r['label']) == 'bed')
    n_steep = sum(1 for r in model['joint_lines'] if set_name(r['label']) == 'st')

    found = sets_in(model)
    if [name for name, _s, _i in found] != ['bed', 'st']:
        failures.append(f"regeneration: the model's sets read as "
                        f"{[n for n, _s, _i in found]}, expected ['bed', 'st']")

    # Halve the spacing: the set regenerates in place, with more rows.
    tighter = JointSet('bed', 'parallel', {'dip': 0.0, 'spacing': 1.0})
    fresh = regenerate(model, 'bed', tighter)
    rows = model['joint_lines']
    got_bed = [r for r in rows if set_name(r['label']) == 'bed']
    got_steep = [r for r in rows if set_name(r['label']) == 'st']
    if len(got_bed) != len(fresh) or len(got_bed) <= n_bed:
        failures.append(f"regeneration: halving the spacing took bed from "
                        f"{n_bed} rows to {len(got_bed)}")
    if len(got_steep) != n_steep:
        failures.append(f"regeneration: the OTHER set went from {n_steep} rows "
                        f"to {len(got_steep)}")
    if sum(1 for r in rows if r['label'] == 'base joint') != 1:
        failures.append("regeneration: the hand-entered joint line did not survive")
    if any(abs(float(r['phi']) - 32.0) > 1e-9 or abs(float(r['c']) - 10.0) > 1e-9
           for r in got_bed):
        failures.append("regeneration: the set's properties were not carried onto "
                        "the new rows")
    # In place: the set's rows are still where they were, ahead of the typed line.
    if set_name(rows[0]['label']) != 'bed' or rows[len(got_bed)]['label'] != 'base joint':
        failures.append("regeneration: the fresh rows did not land where the old "
                        "ones were")

    # No record given at all: the set re-runs exactly as its labels record it.
    again = regenerate(model, 'bed')
    if len(again) != len(got_bed):
        failures.append(f"regeneration: re-running the recorded set gave "
                        f"{len(again)} rows, not {len(got_bed)}")
    if [r['label'] for r in again] != [r['label'] for r in got_bed]:
        failures.append("regeneration: re-running the recorded set changed its "
                        "labels")

    try:
        regenerate(model, 'nosuch')
        failures.append("regeneration: an unknown set name was accepted")
    except ValueError as exc:
        if 'bed' not in str(exc):
            failures.append(f"regeneration: the refusal does not list the model's "
                            f"own sets: {exc}")
    results.append(f"regeneration  one set of two replaced in place ({n_bed} rows "
                   f"-> {len(got_bed)} at half the spacing), the other set, the "
                   f"hand-entered line and the properties untouched; re-running "
                   f"the recorded set reproduces it")


def _leg_band_and_region(failures, results):
    """Leg 10: the elevation band, and a joints-polygon region.

    A band is an artificial cut through material, so it differs from the region's
    own boundary in exactly one way that matters: a trace lying ALONG a band edge
    is kept, where one lying along the section's boundary is dropped. Both
    readings are here, on the same set.
    """
    model = _model()                      # a 20 x 10 box, two materials
    band = parallel_set(model, 0.0, 2.0, band=(4.0, 8.0), props={'phi': 30.0})
    ys = sorted(round(r['y1'], 6) for r in band)
    if ys != [4.0, 6.0, 8.0]:
        failures.append(f"band: the horizontal set in the band 4-8 came out at "
                        f"y = {ys}, expected [4.0, 6.0, 8.0] (the band's own "
                        f"edges kept)")
    whole = parallel_set(model, 0.0, 2.0, props={'phi': 30.0})
    ys_whole = sorted(round(r['y1'], 6) for r in whole)
    if ys_whole != [2.0, 4.0, 6.0, 8.0]:
        failures.append(f"band: the same set over the whole section came out at "
                        f"y = {ys_whole} (the traces on the section's boundary "
                        f"must be dropped)")
    for row in band:
        if row['y1'] < 4.0 - 1e-9 or row['y1'] > 8.0 + 1e-9:
            failures.append(f"band: {row['label']} lies at y = {row['y1']:g}, "
                            f"outside the band")
    # An open-ended band, and one that lies off the section entirely.
    above = parallel_set(model, 0.0, 2.0, band=(6.0, None), props={'phi': 30.0})
    if sorted(round(r['y1'], 6) for r in above) != [6.0, 8.0]:
        failures.append(f"band: the open band above 6 came out at "
                        f"{sorted(round(r['y1'], 6) for r in above)}")
    try:
        parallel_set(model, 0.0, 2.0, band=(50.0, 60.0))
        failures.append("band: a band above the whole section was accepted")
    except ValueError as exc:
        if '50' not in str(exc):
            failures.append(f"band: the refusal does not name the band: {exc}")

    # The joints polygon, by name and by number: the same region either way.
    model['joint_zones'] = [
        {'polygon': [(2.0, 2.0), (18.0, 2.0), (18.0, 8.0), (2.0, 8.0)],
         'label': 'North block', 'size': None, 'mat_id': None}]
    by_name = parallel_set(model, 0.0, 2.0, region='poly:North block',
                           props={'phi': 30.0})
    by_number = parallel_set(model, 0.0, 2.0, region='poly:1',
                             props={'phi': 30.0})
    if not by_name:
        failures.append("region: the joints polygon produced no trace")
    if ([(r['x1'], r['y1'], r['x2'], r['y2']) for r in by_name]
            != [(r['x1'], r['y1'], r['x2'], r['y2']) for r in by_number]):
        failures.append("region: the joints polygon by name and by number gave "
                        "different traces")
    for row in by_name:
        if not (2.0 - 1e-9 <= row['y1'] <= 8.0 + 1e-9
                and 2.0 - 1e-9 <= row['x1'] <= 18.0 + 1e-9):
            failures.append(f"region: {row['label']} runs outside the joint region")
    try:
        parallel_set(model, 0.0, 2.0, region='poly:South block')
        failures.append("region: an unknown joint region name was accepted")
    except ValueError as exc:
        if 'North block' not in str(exc):
            failures.append(f"region: the refusal does not list the model's own "
                            f"regions: {exc}")
    results.append(f"band/region  a band keeps the {len(band)} traces inside it "
                   f"INCLUDING the two on its own edges (the section's boundary "
                   f"still drops them); a joints polygon by name and by number "
                   f"resolve to the same {len(by_name)} traces, and an unknown "
                   f"name is refused with the model's regions listed")


def _leg_display_name(failures, results):
    """Leg 11: the display name — what a PERSON is shown, everywhere.

    A generated row's Label carries the set record, because the sheet has no
    column to keep it in. Nobody reads a set's parameters off a list of a hundred
    rows, so every place a joint line is named to a reader prints
    ``joints.display_label`` of the cell instead — and the stored label, which is
    what a regeneration reads, is never touched by any of it.

    Each site is read where it is used rather than by inspecting the helper:
    ``fem_details._line_label`` (the 1D details list, each member figure's title
    and the report's interface table all name a joint through it), the preflight
    rule that refuses a blank phi, and the joints editor's list column.
    """
    from xslope.fem_details import _line_label
    from xslope.joints import display_label
    from xslope.preflight import _Ctx, _joint_phi_missing
    from xslope.report import _Counter, _joint_table

    model = _model()
    jset = JointSet('bed', 'parallel', {'dip': 25.0, 'spacing': 2.5},
                    props={'phi': 30.0})
    rows = jset.generate(model)
    if len(rows) < 3:
        failures.append(f"display: the fixture set generated {len(rows)} rows")
        return

    # (a) the helper: a set's rows display as their own name, a typed label is
    #     shown as it was typed, and nothing anywhere rewrites the stored cell.
    shown = [display_label(r['label']) for r in rows]
    want = [f"bed-{i + 1:02d}" for i in range(len(rows))]
    if shown != want:
        failures.append(f"display: the set's rows display as {shown}, "
                        f"expected {want}")
    if any('|' in s for s in shown):
        failures.append(f"display: a set record reached a display name: {shown}")
    for typed in ('bedding plane', 'block base', 'bed-07'):
        if display_label(typed) != typed:
            failures.append(f"display: the typed label {typed!r} displayed as "
                            f"{display_label(typed)!r}")
    if [r['label'] for r in rows] != [jset.to_label(i + 1)
                                      for i in range(len(rows))]:
        failures.append("display: the stored labels were altered")

    # (b) the 1D details list, the member figures and the report's interface
    #     table: all three name a joint through _line_label.
    fem_data = {'joint_lines': [dict(r) for r in rows],
                'n_reinforcement_lines': 0, 'n_pile_lines': 0}
    names = [_line_label(fem_data, model, 'reinforcement', i + 1)
             for i in range(len(rows))]
    if names != want:
        failures.append(f"display: the details/report name is {names}, "
                        f"expected {want}")
    table = _joint_table([{'label': names[0], 'length': 4.0,
                           'slipping': [True, False], 'open': [False, False],
                           'slip': [0.01, 0.0], 'status': 'slipping',
                           'units': {}}], _Counter())
    if table.rows[0][0] != want[0]:
        failures.append(f"display: the report's interface table names the line "
                        f"{table.rows[0][0]!r}, expected {want[0]!r}")

    # (c) preflight. The rule that refuses a blank phi names the row it is
    #     refusing, and it must not recite the set's parameters to do it.
    blank = dict(rows[2], phi=float('nan'))
    msg = list(_joint_phi_missing(_Ctx({'joint_lines': [blank]}, 'fem')))
    if len(msg) != 1 or "('bed-03')" not in msg[0] or '|' in msg[0]:
        failures.append(f"display: the preflight refusal reads {msg!r}")

    # (d) the joints editor's list column.
    from studio.editors import _joint_item_label
    item = _joint_item_label(2, rows[2])
    if not item.startswith('3. bed-03 (') or '|' in item:
        failures.append(f"display: the editor's list line reads {item!r}")

    results.append(f"display      the set's {len(rows)} rows are shown as "
                   f"bed-01…bed-{len(rows):02d} by the details list, the member "
                   f"figures, the report's interface table, preflight's refusals "
                   f"and the editor's list, while the file keeps the record each "
                   f"row was generated from")


def run():
    """Returns a list of failure strings (empty = pass)."""
    import time
    failures, results = [], []
    t0 = time.time()
    _leg_counts(failures, results)
    _leg_clipping(failures, results)
    _leg_persistence(failures, results)
    _leg_voronoi(failures, results)
    _leg_refusals(failures, results)
    _leg_mesh(failures, results)
    _leg_roundtrip(failures, results)
    _leg_record(failures, results)
    _leg_regeneration(failures, results)
    _leg_band_and_region(failures, results)
    _leg_display_name(failures, results)
    print(f"Joint network generator check ({time.time() - t0:.0f} s):")
    for line in results:
        print("  " + line)
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nThe generators produce the traces the geometry allows, clipped to "
          "their region, never on its boundary, and a generated network meshes "
          "and round-trips like a typed one.")


if __name__ == '__main__':
    main()
