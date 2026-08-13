"""Render the 4-panel FEM/SSRM figures for the RS2 SSR corpus page
(docs/verification/rs2.md).

Four panels per problem, arranged 2x2, giving the reader the inputs, the mesh,
and the two SSRM field solutions at a glance:

  upper-left  — plot_inputs(mode="fem"): the slope geometry with the FEM-relevant
                overlays (material zones, water table / piezo line, reinforcement,
                loads). This is the "what went in" panel.
  lower-left  — plot_fem_data: the mesh, coloured by material zone, with the
                displacement/boundary conditions drawn.
  upper-right — plot_fem_results('shear_strain'): the viscoplastic max shear
                strain AT FAILURE. Colour fill only (no mesh lines) so the
                failure mechanism reads at panel size.
  lower-right — plot_fem_results('displace_vector'): the viscoplastic displacement
                vectors AT FAILURE, over the mesh boundary outline — the
                kinematics of the same failure mechanism.

The two right panels render the AT-FAILURE (unconverged) field solve_ssrm
captures a margin beyond critical (result['failure_solution']) — the developed
rotational mechanism Griffiths & Lane plot, not the sub-critical last-converged
settlement — so the strain contours and the displacement arrows are two views of
one mechanism. Titles lead with the locked FS ("… at Failure  FS = X").

Every solve also exports the converged + at-failure field sidecars next to the
case xlsx via export_fem_solution, together with the mesh they were solved on
(export_mesh_to_json), so all future re-renders are solve-free and read the same
discretization the fields came off rather than a rebuilt approximation of it.

The cases are parsed straight out of the ``fem_ssrm`` test tags in rs2.md rather
than kept in a second list here. That is deliberate: the figure is then rendered
on the SAME mesh and F-bracket the regression lock uses, so a retagged mesh size
cannot silently leave a stale figure behind.

Rows the page reports without locking carry no tag, so they are registered in
EXTRA_CASES here instead; a row the page documents as reaching NO equilibrium is
registered with ``figure='inputs'`` and gets the inputs panel alone — no solve, and
no invented strain field. ``--audit`` lists every registered row with no rendered
PNG, so a locked row can never again go quietly figure-less. A row whose mechanism
the page already figures under another benchmark — the same model re-run under a
depth filter, an SSR polygon or a later drawdown frame — is named in
``EXPECTED_NO_FIGURE`` with its reason instead; the audit reports those separately
and fails on any exemption that no longer describes a real row (see that dict).
``--audit`` exits non-zero while any row is unexplained.

Each SSRM solve costs about a minute (more on a fine mesh), so a full run is slow.
Pass benchmark ids to render a subset.

Run from the repo root:
    python benchmarks/rocscience/make_rs2_figures.py            # all
    python benchmarks/rocscience/make_rs2_figures.py RS2-30 RS2-31a
    python benchmarks/rocscience/make_rs2_figures.py --audit    # coverage only
"""

import io
import os
import re
import sys
import time
import warnings
import contextlib

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import math

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from xslope.fileio import load_slope_data
from xslope.fem import build_fem_data, solve_ssrm, export_fem_solution
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                         extract_constraint_line_geometry, extract_point_constraints,
                         extract_size_regions, export_mesh_to_json)
from xslope.style import resolve_style, material_style
# The composite is ONE figure whose four panels are drawn into axes the composite
# owns (fixed positions), so all four plot rectangles come out pixel-identical.
# That means calling the LOWER-LEVEL ax-drawing helpers rather than the whole-figure
# wrappers (plot_inputs / plot_fem_data / plot_fem_results), which each own a figure
# and run their own tight_layout / legend / colorbar and so can't share one figure.
from xslope.plot import (
    plot_base_geometry, plot_piezo_line, plot_dloads, plot_tcrack_surface, plot_ssr_zones,
    plot_reinforcement_lines as _plot_input_reinf_lines, plot_piles, plot_line_loads,
    adaptive_colorbar_ticks, adaptive_edge_linewidth,
)
from xslope.plot_fem import (
    plot_shear_strain_contours, plot_displacement_vectors,
    plot_reinforcement_lines as _plot_mesh_reinf_lines, _plot_boundary_conditions,
)

ROOT = os.path.join(os.path.dirname(__file__), '..', '..')
RS2_MD = os.path.join(ROOT, 'docs', 'verification', 'rs2.md')
OUT = os.path.join(ROOT, 'docs', 'verification', 'images')

TAG_RE = re.compile(r'<!--\s*test:\s*(.*?)\s*-->')


# Rows the page documents but does NOT regression-lock, so they carry no fem_ssrm
# test tag for parse_tags to find. Registered here — in the same tag shape — so a
# reported-not-locked row still gets a figure instead of silently having none.
#
# The multi-tiered geotextile wall family (RS2-48–55, Leshchinsky & Han 2004) is the
# case in point: no row in it is locked. The baseline's SSR row is not attempted (RS2
# splits that mesh at the sheets and joins it with slip interfaces, so the page carries
# no figure for it), and the seven parametric variants are reported. All eight share
# what was the baseline's model settings — 1.0 m tri6 mesh, the
# vendor's isotropic at-rest field stress (k0 = 1) and static tensile caps
# (tension_srf off) — but the variants run the auto bracket rather than the baseline's
# narrow one, because the family spans 0.76–1.10.
#
# ``figure='inputs'`` marks a row the page documents as reaching NO equilibrium on the
# corpus mesh. There is no failure mechanism to draw for those, so they get the
# inputs panel alone; nothing is solved and no strain field is invented.
_WALL = dict(element_type='tri6', target_size='1.0', tolerance='0.02',
             f_min='0.5', f_max='3.0', max_iter='16000',
             tension_srf='false', k0='1')

EXTRA_CASES = [
    # RS2 Part IV VP65: solves and brackets, reported rather than locked (its factor
    # and RS2's describe different mechanisms), so it carries no tag. Same settings as
    # its VP64/VP66 siblings.
    {'file': 'files/rocscience/vp065.xlsx', 'benchmark': 'RS2-P4-VP65',
     'element_type': 'tri6', 'target_size': '6.0', 'tolerance': '0.02',
     'f_min': '1.4', 'f_max': '2.8', 'max_iter': '16000',
     'tension_srf': 'true', 'k0': '1'},
    {**_WALL, 'file': 'files/rocscience/vp088.xlsx', 'benchmark': 'RS2-49'},
    {**_WALL, 'file': 'files/rocscience/vp089.xlsx', 'benchmark': 'RS2-50',
     'figure': 'inputs'},
    {**_WALL, 'file': 'files/rocscience/vp090.xlsx', 'benchmark': 'RS2-51-wall',
     'figure': 'inputs'},
    {**_WALL, 'file': 'files/rocscience/vp091.xlsx', 'benchmark': 'RS2-52'},
    {**_WALL, 'file': 'files/rocscience/vp092.xlsx', 'benchmark': 'RS2-53'},
    {**_WALL, 'file': 'files/rocscience/vp093.xlsx', 'benchmark': 'RS2-54',
     'figure': 'inputs'},
    {**_WALL, 'file': 'files/rocscience/vp094.xlsx', 'benchmark': 'RS2-55'},
]


# Rows that deliberately get NO figure of their own, each with the reason.
#
# Most locked rows on this page are a model in their own right and earn a composite.
# A minority are a SECOND run of a model the page already figures — the same mesh
# under a depth filter, an SSR polygon, a later drawdown frame, a thicker layer — and
# for those the question is whether the mechanism the figure would draw differs from
# the one already on the row. Where it does, the variant is rendered (RS2-4-zone,
# RS2-40-deep, RS2-66a-deep, RS2-P4-VP68-zone all show a different band from their
# unconstrained twin). Where it does not, a second composite would repeat the first
# picture, and the row is named here instead.
#
# This is a claim about the page, not a way to empty the audit, so it is enforced
# both ways: an entry naming a row that is not registered, or a row that HAS a PNG,
# is reported as a DEAD exemption and fails the audit exactly like a missing figure.
EXPECTED_NO_FIGURE = {
    # RS2-64 — C3 and C5 are the second and third short-term Original slopes. All
    # three lock unconstrained on a simple convex profile and fail by the same deep
    # rotation; C1 is figured.
    'RS2-64c': 'same unconstrained deep rotation as C1 (RS2-64a.png); only the profile differs',
    'RS2-64e': 'same unconstrained deep rotation as C1 (RS2-64a.png); only the profile differs',

    # RS2-66 filter-off — the face skin is surface-parallel and depth-independent:
    # the page measures 1.056 at every h1, so the mechanism cannot vary with the
    # soft-band thickness. Both ends of the family are figured.
    'RS2-66b': 'depth-independent face skin, identical at every h1; figured at both ends (RS2-66a.png, RS2-66e.png)',
    'RS2-66c': 'depth-independent face skin, identical at every h1; figured at both ends (RS2-66a.png, RS2-66e.png)',
    'RS2-66d': 'depth-independent face skin, identical at every h1; figured at both ends (RS2-66a.png, RS2-66e.png)',

    # RS2-66 filtered — one deep basal squeeze through the soft band, figured at the
    # thinnest band where it separates furthest from the skin.
    'RS2-66b-deep': 'same deep basal squeeze as RS2-66a-deep.png, through a thicker soft band',
    'RS2-66c-deep': 'same deep basal squeeze as RS2-66a-deep.png, through a thicker soft band',
    'RS2-66d-deep': 'same deep basal squeeze as RS2-66a-deep.png, through a thicker soft band',
    'RS2-66e-deep': 'at h1 = 10 m the filter changes nothing (1.056 both ways) — RS2-66e.png IS this run',

    # RS2-67 — six drawdown stages of one dam, two mechanisms between them: the
    # unconstrained downstream face, and the upstream face when RS2's Search Area
    # confines reduction to it. One of each is figured.
    'RS2-67a': 'same downstream-face mechanism as RS2-67b.png, with no pore-pressure field',
    'RS2-67c': 'same unconstrained downstream-face mechanism as RS2-67b.png, at the 90 h field',
    'RS2-67e': 'same unconstrained downstream-face mechanism as RS2-67b.png, at the drained limit',
    'RS2-67f': 'same Search-Area-confined upstream mechanism as RS2-67d.png, at the drained limit',

    # RS2 Part IV VP102 — the drawdown mechanism is one downstream-face wedge at every
    # frame of a monotone sequence. One frame of each case is figured: the phi_b = 0
    # frame the page reports as its worst, and the phi_b = 37 frame that sets the dot.
    'RS2-P4-VP102-t-60-c2': 'same downstream-face wedge as RS2-P4-VP102-t-300-c2.png, at an earlier frame',
    'RS2-P4-VP102-t-1500-c2': 'same downstream-face wedge as RS2-P4-VP102-t-300-c2.png, at a later frame',
    'RS2-P4-VP102-t-60-c3': 'same downstream-face wedge as RS2-P4-VP102-t-1500-c3.png, at an earlier frame',
    'RS2-P4-VP102-t-300-c3': 'same downstream-face wedge as RS2-P4-VP102-t-1500-c3.png, at an earlier frame',
}


# Sidecar stem overrides. ``make_figure`` writes the exported field next to the case
# xlsx under the file's own stem, which is right when a file carries ONE run. Where a
# file carries two — an unconstrained lock and a constrained/filtered variant — the
# second run would overwrite the first one's sidecars, and a later --from-sidecar
# re-render would draw the wrong mechanism under the right caption. Naming the stem
# per benchmark keeps each run's field beside its own figure.
SIDECAR_STEM = {
    'RS2-4-zone': 'vp005_zone',
    'RS2-40-deep': 'vp077b_deep',
    'RS2-66a-deep': 'rs2_66a_deep',
    'RS2-P4-VP68-zone': 'vp068_zone',
    'RS2-P4-VP102-t-300-c2': 'vp102t_300_c2',
    'RS2-P4-VP102-t-1500-c3': 'vp102t_1500_c3',
}


def _sidecar_stem(tag, path):
    """The stem ``export_fem_solution`` / ``import_fem_solution`` use for this row —
    the case xlsx without its extension, unless the benchmark is one of the
    second-run-on-a-shared-file rows named in ``SIDECAR_STEM``."""
    bench = tag.get('benchmark', '')
    override = SIDECAR_STEM.get(bench)
    if override:
        return os.path.join(os.path.dirname(path), override)
    return os.path.splitext(path)[0]


def parse_tags(path=RS2_MD):
    """Every fem_ssrm tag on the page, as a list of key->value dicts. The tag is
    the single source of truth for mesh size, element type and F-bracket."""
    cases = []
    with open(path) as fh:
        for line in fh:
            m = TAG_RE.search(line)
            if not m:
                continue
            kv = {}
            for part in m.group(1).split(','):
                if '=' not in part:
                    continue
                k, v = part.split('=', 1)
                kv[k.strip()] = v.strip()
            if kv.get('type') != 'fem_ssrm':
                continue
            cases.append(kv)
    return cases


def _declare_dry_beyond_piezo(sd, mesh):
    """A piezo line that stops short of the MESH, made explicit for the FEM build.

    A piezometric line assigns pore pressure over its own x-extent and nothing
    beyond it. XSLOPE refuses to sample one from outside that extent rather than
    quietly reading zero (xslope.fem._check_piezo_extent) — silently deleting pore
    pressure below a water table over-predicts strength, so the model has to say
    which it means. VP65 is the case in point: its pool is upstream only, its line
    ends at x = 117.778 where pool elevation meets the downstream face, and RS2's
    own import of the same problem solves hydrostatic-to-el-20 there and exactly
    zero beyond. Its LEM circles sit inside the line and solve untouched; only the
    full-width FEM mesh reaches past it.

    So for the figure the line is carried on BELOW the mesh past its own end — the
    explicit spelling of "dry here", giving node-for-node the same zero field the
    truncated line implied. In memory and for the FEM build only: the corpus file
    keeps the vendor's line, and the inputs panel still draws it as authored.
    """
    line = sd.get('piezo_line') or []
    if len(line) < 2:
        return sd
    if not any(str(m.get('u', '')).strip().lower() == 'piezo'
               for m in sd.get('materials', [])):
        return sd
    pts = sorted(((float(x), float(y)) for x, y in line), key=lambda p: p[0])
    xs = mesh['nodes'][:, 0]
    x_lo, x_hi = float(np.min(xs)), float(np.max(xs))
    if pts[0][0] <= x_lo and pts[-1][0] >= x_hi:
        return sd                       # the line already covers every node
    # An elevation strictly below the mesh: no node can sit under it, so u = 0.
    y_dry = float(np.min(mesh['nodes'][:, 1])) - 1.0
    span = max(x_hi - x_lo, 1.0)
    eps = 1e-6 * span                   # the drop needs a non-zero x to interpolate over
    ext = list(pts)
    if pts[0][0] > x_lo:
        ext = [(x_lo, y_dry), (pts[0][0] - eps, y_dry)] + ext
    if pts[-1][0] < x_hi:
        ext = ext + [(pts[-1][0] + eps, y_dry), (x_hi, y_dry)]
    out = dict(sd)
    out['piezo_line'] = ext
    return out


def _uses_seep(sd):
    """True where a material reads its pore pressures off stored nodal seepage
    values. Such a model carries its mesh as part of its definition: the solve has
    to reuse that discretization, and the inputs panel has to show it."""
    return any(str(m.get('u', '')).strip().lower() == 'seep'
               for m in sd.get('materials', []))


def _build(tag):
    """Mesh + fem_data exactly as run_tests.run_fem_test does. Returns
    ``(sd, fem_data, path, mesh)`` — ``path`` is the resolved case xlsx, so the
    sidecars can be written next to it, and ``mesh`` is the discretization the
    fields were built on, which has to be exported beside them (see make_figure)."""
    # tag paths are relative to docs/verification/
    path = os.path.normpath(os.path.join(ROOT, 'docs', 'verification', tag['file']))
    sd = load_slope_data(path)

    # seepage-coupled models must reuse the stored mesh (nodal seep_u)
    if _uses_seep(sd) and sd.get('mesh') is not None:
        mesh = sd['mesh']
    else:
        target = tag.get('target_size')
        if target is None:
            xs = [x for x, _ in sd['ground_surface'].coords]
            target = (max(xs) - min(xs)) / 100
        else:
            target = float(target)
        lines, _nr, _np = extract_constraint_line_geometry(sd)
        polys = get_material_polygons(sd, reinf_lines=lines)
        # Feature-aware mesh refinement, matching run_tests._refine_kwargs: a tag
        # carrying refine_factor=N (and optional refine_features=a;b) must build the
        # SAME refined mesh the lock was recorded on, or the figure would show the
        # unrefined coarse mesh and a different FS.
        refine_kw = {}
        if tag.get('refine_factor') is not None:
            refine_kw['refine_factor'] = float(tag['refine_factor'])
            feats = tag.get('refine_features')
            if feats:
                refine_kw['refine_features'] = [s.strip() for s in str(feats).split(';')
                                                if s.strip()]
        # The model's own mesh refinement polygons (polygon sheet Type='refine')
        # travel with the file the same way its material zones do, so they belong
        # here for the same reason refine_factor does: a file that refines a band
        # is solved on the refined mesh by the suite, and a figure built without
        # the overlay would draw a coarser mesh and report a different FS.
        mesh = build_mesh_from_polygons(
            polys, target_size=target,
            element_type=tag.get('element_type', 'tri6'), lines=lines,
            point_constraints=extract_point_constraints(sd),
            size_regions=extract_size_regions(sd), **refine_kw)
    # `sd` (unmodified) is what the inputs panel draws; the FEM build gets the
    # dry-beyond-the-line spelling where a piezo line stops short of the mesh.
    return sd, build_fem_data(_declare_dry_beyond_piezo(sd, mesh), mesh), path, mesh


def _tag_ssr_zone(tag):
    """The tag's ``ssr_zone`` polygon as a vertex list, or None. Flat
    x1;y1;x2;y2;… (semicolon-separated, since tag key=value pairs are comma-split)."""
    if not tag.get('ssr_zone'):
        return None
    _z = [float(v) for v in str(tag['ssr_zone']).split(';') if v.strip() != '']
    return list(zip(_z[0::2], _z[1::2]))


def _clip_zone_to_mesh(polygon, fem_data, tol=0.15, pad=0.05):
    """The part of a search-area polygon that lies ON the model, as a vertex list.

    Vendors write some of these rectangles far outside the section — RS2's #067
    upstream area is 97 m tall on a 29 m dam — because only the intersection with the
    mesh can be reduced and the overhang costs them nothing. It costs the FIGURE a
    great deal: the inputs panel takes its extent from everything it draws, and the
    composite imposes that extent on all four panels, so one overhanging rectangle
    shrinks the dam to a third of the panel in every one of them. Clipping draws the
    same reduced region and frames the model. Display only — the solver always gets the
    polygon exactly as the tag writes it.

    Only a LARGE overhang is clipped: a polygon staying within ``tol`` of the mesh
    box on every side is left alone, so the many zones that merely graze or slightly
    overshoot a boundary (RS2-64's corridors, RS2-4's exclusion ring) are drawn to the
    vertex as authored, and only a rectangle that is actually distorting the frame is
    touched. Returns the polygon unchanged if the clip degenerates — nothing to draw
    would be worse than a loose frame."""
    from shapely.geometry import Polygon, box
    nodes = np.asarray(fem_data['nodes'])
    x0, x1 = float(np.min(nodes[:, 0])), float(np.max(nodes[:, 0]))
    y0, y1 = float(np.min(nodes[:, 1])), float(np.max(nodes[:, 1]))
    w, h = max(x1 - x0, 1e-9), max(y1 - y0, 1e-9)
    poly = Polygon(polygon)
    if not poly.is_valid:
        poly = poly.buffer(0)
    if poly.is_empty:
        return polygon
    if poly.within(box(x0 - w * tol, y0 - h * tol, x1 + w * tol, y1 + h * tol)):
        return polygon
    clipped = poly.intersection(box(x0 - w * pad, y0 - h * pad,
                                    x1 + w * pad, y1 + h * pad))
    if clipped.is_empty or clipped.geom_type != 'Polygon':
        return polygon
    return [(float(x), float(y)) for x, y in clipped.exterior.coords]


def _record_ssr_zone(sd, ssr_zone, fem_data=None):
    """A tag-carried ssr_zone is a constraint the reader must be able to see, exactly
    like a file-carried one. It is handed to solve_ssrm as a kwarg, not through the
    file, so slope_data knows nothing about it and the inputs panel would draw an
    unconstrained model beside a constrained mechanism. Record it for drawing only,
    never for solving — every caller applies this AFTER its solve (or instead of one,
    on the sidecar path), so the kwarg stays the single source of truth for the solver
    and no run ever sees both a kwarg polygon and a file polygon at once. The tag's
    semantics are solve_ssrm's: reduce inside, i.e. a search area.

    ``fem_data``, when given, clips the DRAWN polygon to the model — see
    ``_clip_zone_to_mesh`` for why an overhanging rectangle has to be clipped."""
    if ssr_zone and not sd.get('ssr_zones'):
        drawn = (_clip_zone_to_mesh(ssr_zone, fem_data) if fem_data is not None
                 else list(ssr_zone))
        sd = dict(sd)
        sd['ssr_zones'] = [{'kind': 'reduce', 'polygon': list(drawn),
                            'label': 'SSR reduce'}]
    return sd


def build_and_solve(tag):
    """Build the mesh, run the SSRM bracket, and return the pieces the figure and
    the sidecars need: (sd, fem_data, field, failure, FS, path, mesh).

    ``field`` is the last-CONVERGED field (the F just below critical); it is what
    the export writes as the converged sidecar. ``failure`` is the AT-FAILURE
    (unconverged) mechanism solve_ssrm captures a margin beyond critical
    (result['failure_solution']) — the developed rotational mechanism the strain
    and displacement-vector panels render, and the second sidecar pair. ``failure``
    is None if the capture was skipped/failed, in which case the panels fall back
    to the converged field. ``path`` is the case xlsx (for the sidecar stem), and
    ``mesh`` is the discretization the fields were solved on — exported beside them,
    since a reload that has to guess at the mesh reads a section that was never solved.
    """
    sd, fem_data, path, mesh = _build(tag)

    # SSR-exclusion material names (semicolon-separated within the tag value,
    # since tag key=value pairs are comma-split) — held at full strength.
    ssr_exclude = None
    if tag.get('ssr_exclude'):
        ssr_exclude = [s.strip() for s in str(tag['ssr_exclude']).split(';') if s.strip()]

    # SSR search-area polygon (RS2's "SSR Search Area"): flat x1;y1;x2;y2;... list.
    ssr_zone = _tag_ssr_zone(tag)

    # Elastic-materials (RS2's "Plasticity Specifications: None"): purely-elastic
    # can't-fail zones, semicolon-separated within the tag value. Must be passed
    # or the figure would solve unconstrained and show a different FS/mechanism.
    elastic_materials = None
    if tag.get('elastic_materials'):
        elastic_materials = [s.strip() for s in str(tag['elastic_materials']).split(';')
                             if s.strip()]

    # Solver settings that select a different initial state or a different
    # mechanism, and so change the answer the figure draws. Both must be passed for
    # the same reason elastic_materials must: k0 pins the in-situ stress the vendor
    # model starts from (RS2's isotropic Kx = Kz = 1), and min_slip_depth selects the
    # deep mechanism on a row whose unfiltered minimum is a surficial skin. Left out,
    # the figure would re-solve a different problem than its own lock.
    extra = {}
    if 'k0' in tag:
        extra['k0'] = float(tag['k0'])
    if 'min_slip_depth' in tag:
        extra['min_slip_depth'] = float(tag['min_slip_depth'])
    # Matric-suction strength (Fredlund extended MC), mirroring run_tests: a per-material
    # "Name:deg" list, semicolon-separated within the tag value. It is a STRENGTH term —
    # the phi_b = 37 twins on VP102 read a full 0.38 above their phi_b = 0 partners on the
    # same file and the same field — so a figure that dropped it would draw a different
    # mechanism at a different factor from its own lock.
    if tag.get('suction_phi_b'):
        sp = {}
        for tok in str(tag['suction_phi_b']).split(';'):
            tok = tok.strip()
            if not tok:
                continue
            name, sep, deg = tok.rpartition(':')
            if not sep or not name.strip() or deg.strip() == '':
                raise ValueError(f"suction_phi_b entry {tok!r} must be 'Name:degrees'")
            sp[name.strip()] = float(deg)
        extra['suction_phi_b'] = sp or None
    if tag.get('suction_cap'):
        extra['suction_cap'] = float(tag['suction_cap'])

    with contextlib.redirect_stdout(io.StringIO()):
        sol = solve_ssrm(fem_data,
                         F_min=float(tag.get('f_min', 0.5)),
                         F_max=float(tag.get('f_max', 3.0)),
                         tolerance=float(tag.get('tolerance', 0.02)),
                         max_iterations=int(tag.get('max_iter', 4000)),
                         ssr_exclude=ssr_exclude,
                         ssr_zone=ssr_zone,
                         elastic_materials=elastic_materials,
                         **extra,
                         # tension_srf mirrors run_tests: divide each tension cap by
                         # the trial F (RS2's tensilestrength_SRF=1). Without this the
                         # figure re-solve would disagree with the tag's lock on any
                         # benchmark that uses it (RS2-62c: 0.744 static vs 0.769 SRF).
                         # Absent tag -> the solver's own default (on); an explicit
                         # tension_srf=false in a tag is honored, so the figure and the
                         # lock always solve the same envelope.
                         tension_srf=(str(tag['tension_srf']).lower() in ('true','1','yes')
                                      if 'tension_srf' in tag else True),
                         # capture_failure_state is default-on; keep it explicit so
                         # the at-failure mechanism (the right-panel field) is always
                         # captured and exported.
                         capture_failure_state=True,
                         debug_level=0)
    if not sol.get('converged'):
        raise RuntimeError(f'SSRM did not converge: {sol.get("error")}')

    # Record the tag-carried zone for DRAWING only — after the solve, so the kwarg
    # stays the solver's single source of truth (see _record_ssr_zone).
    sd = _record_ssr_zone(sd, ssr_zone, fem_data)

    field = sol.get('last_solution')
    if field is None:
        raise RuntimeError('SSRM returned no last_solution to plot')
    # The at-failure (unconverged) mechanism for the strain/vector panels; None if
    # the capture was skipped/failed (the panels then fall back to ``field``).
    failure = sol.get('failure_solution')
    return sd, fem_data, field, failure, sol['FS'], path, mesh


# ── Composite geometry (inches). These are the figure's structural chrome —
# margins, gutters, the data-axes width, the colorbar-slot pitch — NOT per-figure
# fudge factors: every panel's plot rectangle is DERIVED from them plus the slope's
# own aspect, so all four come out identical by construction. Legend-band heights
# are MEASURED (two-pass), never guessed.
_WD_IN = 6.0        # data-axes width; height = _WD_IN × (padded-domain aspect)
_ML, _MR = 0.75, 0.30      # left / right outer margins
_CG = 0.70                 # gutter between the two columns
_MT, _MB = 0.55, 0.35      # top / bottom outer margins
_ROW_GAP = 0.55            # between rows — holds the lower row's title
_CAX_W = 0.16              # colorbar strip width
_CAX_GAP = 0.14            # gap from a data axes to its first colorbar
_CAX_PITCH = 0.82          # slot pitch (bar + its tick/rotated-label band)
_PAD_FRAC = 0.04           # uniform content cushion, both axes, in data units (4%)
_TITLE_FS = 12             # one title size for every panel
_TITLE_PAD = 8             # one title pad for every panel
_LEG_FS = 8                # one legend size for both left panels
_LEG_CLEAR = 4.0           # points: the legend hangs this far below the panel's
                           # x-tick labels — a FIXED gap (converted to each panel's
                           # own axes fraction at draw time), so the legend drop never
                           # scales with panel height. A constant -0.11 axes-fraction
                           # anchor dropped ~0.11×panel-height, which on tall panels
                           # pushed the inputs legend past _ROW_GAP into the lower
                           # row's mesh-panel title (the collision this closes).
# One rcParams size for ALL text in the figure (root cause of Norm's mismatched
# fonts was four figures rescaled to a common width, scaling their fonts apart).
_RC = {'font.size': 10, 'axes.titlesize': _TITLE_FS, 'axes.labelsize': 10,
       'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': _LEG_FS}


def _at_failure_field(field, failure, fs):
    """The field the strain + vector panels render, with the markers the ax-level
    drawers read for their "… at Failure  FS = X" titles — mirroring what
    plot_fem_results builds internally (deform_field / contour_field). The at-failure
    (unconverged) mechanism when captured, else the converged field."""
    base = failure if failure is not None else field
    d = {**base}
    if failure is not None:
        d['_at_failure'] = True
    if fs is not None:
        d['_ssrm_fs'] = fs
    return d


def _draw_inputs_panel(ax, sd, style):
    """Draw the FEM-inputs panel into ``ax`` by replaying the exact drawer sequence
    plot_inputs(mode='fem') uses — geometry / material zones / piezo line / dloads /
    tension crack / reinforcement / piles / line loads — but into a caller-owned
    axes (no figure ownership, no limits: the composite imposes uniform limits).
    Returns (handles, labels) for the legend the composite draws below the panel."""
    from matplotlib.collections import LineCollection
    from xslope.style import feature_style

    # Background mesh, drawn only where the mesh is part of the MODEL: a seep-coupled
    # case reads its pore pressures off stored nodal values, so the discretization is
    # an input the reader has to see. The test is what the model does with the mesh,
    # not whether a {stem}_mesh.json happens to lie beside the xlsx — since make_figure
    # writes exactly that file, keying on the file's presence made the panel depend on
    # whether the builder had run before: the first row on a shared file drew a clean
    # section and the next row on the SAME file drew a meshed one (RS2-4 / RS2-4-zone),
    # and a second pass over the corpus would mesh every panel that a first pass left
    # clean. Width finalized after layout so a dense grid stays a hairline.
    mesh = sd.get('mesh') if _uses_seep(sd) else None
    bg_lc, bg_segs = None, None
    if mesh is not None:
        m_nodes, m_elems, m_et = mesh['nodes'], mesh['elements'], mesh['element_types']
        segs = []
        for elem, et in zip(m_elems, m_et):
            if et in (3, 6):
                edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
            elif et in (4, 8, 9):
                edges = [(elem[0], elem[1]), (elem[1], elem[2]),
                         (elem[2], elem[3]), (elem[3], elem[0])]
            else:
                continue
            segs.extend(m_nodes[[a, b]] for a, b in edges)
        if segs:
            mfs = feature_style(style, 'mesh')
            bg_lc = LineCollection(segs, colors=mfs.get('color', 'gray'),
                                   alpha=mfs.get('alpha', 0.25),
                                   linewidths=mfs.get('linewidth', 0.5),
                                   linestyles=mfs.get('linestyle', '-'))
            ax.add_collection(bg_lc)
            bg_segs = segs

    plot_base_geometry(ax, sd, labels=True, style=style)
    # File-carried SSR zone overlays (v20 polygon-sheet sentinel rows). plot_inputs
    # draws these on mode='fem' and this panel is that plot, so omitting them hid a
    # real constraint: the vendor exclusion/search polygons a corpus file carries
    # decide which elements the SSRM may reduce, and a reader comparing the inputs
    # panel with the strain panel has to be able to see them.
    plot_ssr_zones(ax, sd, style=style)
    plot_piezo_line(ax, sd, style=style)
    plot_dloads(ax, sd, style=style)
    plot_tcrack_surface(ax, sd, style=style)
    _plot_input_reinf_lines(ax, sd, style=style)
    plot_piles(ax, sd, style=style)
    plot_line_loads(ax, sd, style=style)

    # Seismic coefficient badge (FEM: sign of k gives the direction) — same as
    # plot_inputs, so pseudo-static cases carry the k = … annotation.
    k_seismic = sd.get('k_seismic', 0.0)
    if k_seismic:
        arrow = '→' if k_seismic > 0 else '←'
        ax.text(0.98, 0.97, f'k = {abs(k_seismic):g} {arrow}'.strip(),
                transform=ax.transAxes, fontsize=10, fontweight='bold',
                ha='right', va='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                          edgecolor='orange', linewidth=1.0, alpha=0.9))

    # The legend is read off the DRAWN artists and nothing else, de-duplicated by
    # label — the discipline plot_inputs and _legend_below hold in the package. Every
    # layer above names what it draws, the load sets included: plot_dloads labels the
    # surface line of EACH block with its set's key, so a model with several typed
    # blocks put several identical entries in this list, and an entry appended here
    # from sd['dloads'] put one more on top of them. vp009 printed "Distributed Load"
    # three times, rs2_29clay twice. One drawn thing, one key.
    handles, labels = ax.get_legend_handles_labels()
    seen, dh, dl = set(), [], []
    for h, l in zip(handles, labels):
        if l in seen or (isinstance(l, str) and l.startswith('_')):
            continue
        seen.add(l)
        dh.append(h)
        dl.append(l)
    return dh, dl, (bg_lc, bg_segs)


def _draw_mesh_panel(ax, fem_data, style, alpha=0.6):
    """Draw the mesh / material-zone panel into ``ax`` by replaying plot_fem_data's
    fill + boundary + reinforcement + boundary-condition drawing (the one panel with
    no single ax-level helper — the material fills are inline in plot_fem_data), into
    a caller-owned axes. Returns (legend_handles, title). No limits/aspect set: the
    composite imposes uniform limits."""
    from matplotlib.collections import PatchCollection, LineCollection
    from matplotlib.patches import Polygon
    import matplotlib.patches as patches

    nodes = fem_data['nodes']
    elements = fem_data['elements']
    element_materials = fem_data['element_materials']
    element_types = fem_data.get('element_types')
    if element_types is None:
        element_types = np.full(len(elements), 3)
    materials = np.unique(element_materials)
    mat_to_color = {mat: material_style(style, int(mat) - 1)['color'] for mat in materials}

    # Batch fill polygons per material (quadratic elements subdivided into sub-tris /
    # sub-quads with no interior edges, so the fill reads as a clean zone).
    mat_fill = {mat: [] for mat in materials}
    edge_segments = []
    for idx, en in enumerate(elements):
        et = element_types[idx]
        mat = element_materials[idx]
        if et == 3:
            mat_fill[mat].append(nodes[en[:3]])
        elif et == 6:
            n0, n1, n2 = nodes[en[0]], nodes[en[1]], nodes[en[2]]
            n3, n4, n5 = nodes[en[3]], nodes[en[4]], nodes[en[5]]
            mat_fill[mat].extend([np.array([n0, n3, n5]), np.array([n3, n1, n4]),
                                  np.array([n5, n4, n2]), np.array([n3, n4, n5])])
            edge_segments.extend([[n0, n1], [n1, n2], [n2, n0]])
        elif et == 4:
            mat_fill[mat].append(nodes[en[:4]])
        elif et in (8, 9):
            c = nodes[en[:8]].mean(axis=0) if et == 8 else nodes[en[8]]
            n0, n1, n2, n3 = nodes[en[0]], nodes[en[1]], nodes[en[2]], nodes[en[3]]
            n4, n5, n6, n7 = nodes[en[4]], nodes[en[5]], nodes[en[6]], nodes[en[7]]
            mat_fill[mat].extend([np.array([n0, n4, c, n7]), np.array([n4, n1, n5, c]),
                                  np.array([c, n5, n2, n6]), np.array([n7, c, n6, n3])])
            edge_segments.extend([[n0, n1], [n1, n2], [n2, n3], [n3, n0]])
    for mat in materials:
        polys = mat_fill[mat]
        if not polys:
            continue
        has_edge = any(element_types[i] in (3, 4)
                       for i, m in enumerate(element_materials) if m == mat)
        ec = 'k' if has_edge else 'none'
        lw = 0.5 if has_edge else 0.0
        pc = PatchCollection([Polygon(p) for p in polys], facecolor=mat_to_color[mat],
                             edgecolor=ec, linewidth=lw, alpha=alpha, gid='MESH_FILL')
        ax.add_collection(pc)
    if edge_segments:
        ax.add_collection(LineCollection(edge_segments, colors='k', linewidths=0.5, gid='MESH'))

    material_names = fem_data.get('material_names', [])
    legend_handles = []
    for mat in materials:
        label = material_names[mat - 1] if material_names and mat <= len(material_names) else f'Material {mat}'
        legend_handles.append(patches.Patch(facecolor=mat_to_color[mat], alpha=alpha,
                                            edgecolor='none', label=label))

    # 1D elements (reinforcement truss / pile beam)
    elements_1d = fem_data.get('elements_1d', np.array([]).reshape(0, 3))
    pile_mask = fem_data.get('pile_elem_mask', np.zeros(len(elements_1d), dtype=bool))
    n_reinf = n_pile = 0
    if len(elements_1d) > 0:
        reinf_segs, pile_segs = [], []
        for i in range(len(elements_1d)):
            seg = [nodes[elements_1d[i][0]], nodes[elements_1d[i][1]]]
            (pile_segs if pile_mask[i] else reinf_segs).append(seg)
            n_pile += int(pile_mask[i]); n_reinf += int(not pile_mask[i])
        if reinf_segs:
            ax.add_collection(LineCollection(reinf_segs, colors='red', linewidths=2.5,
                                             zorder=5, gid='REINFORCEMENT'))
            legend_handles.append(plt.Line2D([0], [0], color='red', lw=2.5,
                                            label=f'Reinforcement ({n_reinf} elements)'))
        if pile_segs:
            ax.add_collection(LineCollection(pile_segs, colors='green', linewidths=3.5,
                                             zorder=5, gid='PILES'))
            legend_handles.append(plt.Line2D([0], [0], color='green', lw=3.5,
                                            label=f'Pile ({n_pile} elements)'))

    _plot_boundary_conditions(ax, nodes, fem_data['bc_type'], fem_data['bc_values'],
                              legend_handles, 0.03, fem_data.get('roller_x_nodes', set()))

    num_tri = int(np.sum((element_types == 3) | (element_types == 6)))
    num_quad = int(np.sum((element_types == 4) | (element_types == 8) | (element_types == 9)))
    parts = []
    if num_tri:
        parts.append(f'{num_tri} triangles')
    if num_quad:
        parts.append(f'{num_quad} quads')
    if n_reinf:
        parts.append(f'{n_reinf} reinforcement')
    if n_pile:
        parts.append(f'{n_pile} pile')
    title = f"FEM Mesh with Material Zones ({', '.join(parts)})"
    return legend_handles, title


def _fit_ncol(ax, fig, handles, labels, anchor, fs):
    """Fewest legend rows whose balanced columns fit within THIS panel's axes width
    (not the whole figure — _fit_legend_ncol budgets against the figure, too wide for
    a per-panel legend). Same try-1-row-then-2 search, measured against the drawn
    axes box. Falls back to two rows without a renderer."""
    n = len(labels)
    if n <= 1:
        return max(1, n)
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        budget = ax.get_window_extent().width
        for rows in range(1, n + 1):
            ncol = math.ceil(n / rows)
            leg = ax.legend(handles=handles, labels=labels, loc='upper center',
                            bbox_to_anchor=anchor, ncol=ncol, fontsize=fs,
                            frameon=False, handlelength=1.6, columnspacing=1.2,
                            handletextpad=0.5)
            fits = leg.get_window_extent(renderer).width <= budget
            leg.remove()
            if fits:
                return ncol
        return 1
    except Exception:
        return max(1, math.ceil(n / 2))


def _inputs_content_datalim(sd, style):
    """dataLim (x0, x1, y0, y1) of the FEM-inputs panel drawn on a throwaway axes. The
    inputs panel is the ONLY one whose content can exceed the mesh-node bbox: the three
    field panels are mesh-based by construction, but the inputs panel also draws
    distributed-load arrows above the ground surface and load lines that can reach past
    the toe. Read so the shared composite domain can contain ALL drawn content."""
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            _draw_inputs_panel(ax, sd, style)
        dl = ax.dataLim
        return float(dl.x0), float(dl.x1), float(dl.y0), float(dl.y1)
    finally:
        plt.close(fig)


def _composite_domain(fem_data, sd, style):
    """The padded data domain imposed on all four panels (returns X0, X1, Y0, Y1).

    Base = mesh-node bbox + a uniform _PAD_FRAC cushion — unchanged for the vast
    majority of cases. But the composite imposes ONE shared domain on all four panels
    for pixel-uniformity, and the inputs panel can draw content OUTSIDE the mesh bbox
    (distributed-load arrows above the surface; a load line past the toe). When that
    content would fall outside the base domain, the domain is re-padded around the
    UNION of the mesh bbox and the inputs content, so every panel keeps a strictly
    positive content cushion (the _verify_uniform gate) and the four plot rectangles
    stay identical. Byte-identical to the mesh-only domain whenever the inputs content
    already fits inside the base pad (76 of the 84 corpus cases)."""
    nodes = fem_data['nodes']
    x0, x1 = float(nodes[:, 0].min()), float(nodes[:, 0].max())
    y0, y1 = float(nodes[:, 1].min()), float(nodes[:, 1].max())
    pad = _PAD_FRAC * max(x1 - x0, y1 - y0)
    if pad <= 0:
        pad = 1.0
    X0, X1, Y0, Y1 = x0 - pad, x1 + pad, y0 - pad, y1 + pad
    ix0, ix1, iy0, iy1 = _inputs_content_datalim(sd, style)
    if ix0 < X0 or ix1 > X1 or iy0 < Y0 or iy1 > Y1:
        ux0, ux1 = min(x0, ix0), max(x1, ix1)
        uy0, uy1 = min(y0, iy0), max(y1, iy1)
        upad = _PAD_FRAC * max(ux1 - ux0, uy1 - uy0)
        if upad <= 0:
            upad = 1.0
        X0, X1, Y0, Y1 = ux0 - upad, ux1 + upad, uy0 - upad, uy1 + upad
    return X0, X1, Y0, Y1


def _build_composite(bench, sd, fem_data, afield, style, leg0_in, leg1_in, dpi, domain):
    """Build the whole 2×2 composite as ONE figure. Every panel is drawn into an
    axes the composite places at a FIXED rectangle, so with the same imposed limits
    and equal aspect the four plot rectangles are pixel-identical by construction.
    ``domain`` (X0, X1, Y0, Y1 from _composite_domain) is the shared padded data
    domain imposed on all four panels. Returns (fig, axes4, (ax_ul, ax_ll),
    field_specs, cbars) where axes4 = the four data axes in UL, UR, LL, LR order and
    field_specs are the strain-panel colorbar (mappable, label) pairs."""
    X0, X1, Y0, Y1 = domain
    DW, DH = X1 - X0, Y1 - Y0
    aspect = DH / DW if DW > 0 else 0.4

    # Colorbar slot count: field bar + any reinforcement-force bars on the strain
    # panel. Reserve the SAME horizontal budget in the right column either way.
    # (Determined by a throwaway probe below via the strain drawer's returned specs;
    # here we size to the worst case using elements_1d presence.)
    Hd = _WD_IN * aspect
    fig = plt.figure(figsize=(1, 1), dpi=dpi)   # size set after we know n_slots

    def rect(xi, yi, wi, hi, W, H):
        return [xi / W, yi / H, wi / W, hi / H]

    # First draw the strain panel on a scratch axes to learn n_slots (field + reinf).
    # Cheap (no solve); discarded. Then lay out the real figure.
    probe_ax = fig.add_axes([0, 0, 1, 1])
    sm_probe, reinf_probe = plot_shear_strain_contours(
        probe_ax, fem_data, afield, show_mesh=False, show_reinforcement=True,
        single_panel=True)
    n_slots = 1 + len(reinf_probe)
    fig.clear()

    Cb = _CAX_GAP + n_slots * _CAX_PITCH
    W = _ML + _WD_IN + _CG + _WD_IN + Cb + _MR
    H = _MB + leg1_in + Hd + _ROW_GAP + leg0_in + Hd + _MT
    fig.set_size_inches(W, H)

    x_c0 = _ML
    x_c1 = _ML + _WD_IN + _CG
    y_r1 = _MB + leg1_in
    y_r0 = _MB + leg1_in + Hd + _ROW_GAP + leg0_in

    ax_ul = fig.add_axes(rect(x_c0, y_r0, _WD_IN, Hd, W, H))
    ax_ur = fig.add_axes(rect(x_c1, y_r0, _WD_IN, Hd, W, H))
    ax_ll = fig.add_axes(rect(x_c0, y_r1, _WD_IN, Hd, W, H))
    ax_lr = fig.add_axes(rect(x_c1, y_r1, _WD_IN, Hd, W, H))

    in_h, in_l, _bg = _draw_inputs_panel(ax_ul, sd, style)
    mesh_handles, mesh_title = _draw_mesh_panel(ax_ll, fem_data, style)
    strain_mappable, reinf_specs = plot_shear_strain_contours(
        ax_ur, fem_data, afield, show_mesh=False, show_reinforcement=True,
        single_panel=True)
    plot_displacement_vectors(
        ax_lr, fem_data, afield, show_mesh=False, show_reinforcement=True,
        plot_boundary=True, displacement_tolerance=0.5, single_panel=True)

    # Impose the SAME padded limits + equal aspect on every panel. The fixed rect
    # already carries the padded-domain aspect, so 'box' is a no-op here (no shrink,
    # no re-centering) — the four data rectangles stay identical.
    for ax in (ax_ul, ax_ur, ax_ll, ax_lr):
        ax.set_xlim(X0, X1)
        ax.set_ylim(Y0, Y1)
        ax.set_aspect('equal', adjustable='box')
        ax.set_axisbelow(True)

    # One title size + pad for all four (the drawers set their own fontsize/pad —
    # re-apply uniformly so title-to-frame spacing matches across the figure).
    ax_ul.set_title('FEM Inputs')
    ax_ll.set_title(mesh_title)
    for ax in (ax_ul, ax_ur, ax_ll, ax_lr):
        ax.set_title(ax.get_title(), fontsize=_TITLE_FS, pad=_TITLE_PAD)

    # Colorbars live in fixed slots to the RIGHT of the strain panel — separate axes
    # that never steal width from any host, so all four hosts stay equal (the
    # invisible-slot intent of 1043fc3, realized by fixed geometry instead of
    # make_axes_locatable, which shrinks hosts under equal aspect).
    field_specs = ([(strain_mappable, 'VP Max Shear Strain')]
                   + [(m, l) for (m, l) in reinf_specs])
    cbars = []
    for j, (mp, lbl) in enumerate(field_specs):
        if mp is None:
            continue
        cx = x_c1 + _WD_IN + _CAX_GAP + j * _CAX_PITCH
        cax = fig.add_axes(rect(cx, y_r0, _CAX_W, Hd, W, H))
        cb = fig.colorbar(mp, cax=cax)
        cb.set_label(lbl, rotation=270, labelpad=14)
        cbars.append(cb)

    # Legends below the two left panels — anchored to their axes (fixed rects), so
    # they drop into the reserved leg-band without perturbing any plot rectangle. The
    # anchor's vertical offset is a FIXED distance (the panel's rendered x-tick-label
    # band plus _LEG_CLEAR), converted to THIS panel's axes fraction — so the drop is
    # constant in pixels and does NOT scale with panel height. _measure_legend_band
    # reserves the full occupied span (this drop + the legend height), so the legend
    # can never reach the lower row's title. (Root cause of the old collision: a fixed
    # -0.11 axes-fraction anchor dropped ~0.11×panel-height, unbounded on tall panels.)
    fig.canvas.draw()
    _rr = fig.canvas.get_renderer()

    def _leg_anchor(ax):
        ext = ax.get_window_extent()
        labels = [t for t in ax.get_xticklabels() if t.get_text()]
        band_px = (ext.y0 - min(t.get_window_extent(_rr).y0 for t in labels)
                   if labels else 0.0)
        drop_px = band_px + _LEG_CLEAR / 72.0 * fig.dpi
        return (0.5, -drop_px / ext.height)

    if in_h:
        anchor = _leg_anchor(ax_ul)
        ncol = _fit_ncol(ax_ul, fig, in_h, in_l, anchor, _LEG_FS)
        ax_ul.legend(in_h, in_l, loc='upper center', bbox_to_anchor=anchor, ncol=ncol,
                     frameon=False, fontsize=_LEG_FS, handlelength=1.6,
                     columnspacing=1.2, handletextpad=0.5)
    if mesh_handles:
        anchor = _leg_anchor(ax_ll)
        ncol = _fit_ncol(ax_ll, fig, mesh_handles, [h.get_label() for h in mesh_handles],
                         anchor, _LEG_FS)
        ax_ll.legend(handles=mesh_handles, loc='upper center', bbox_to_anchor=anchor,
                     ncol=ncol, frameon=False, fontsize=_LEG_FS, handlelength=1.6,
                     columnspacing=1.2, handletextpad=0.5)

    fig.canvas.draw()
    for cb in cbars:
        adaptive_colorbar_ticks(fig, cb)
    return fig, (ax_ul, ax_ur, ax_ll, ax_lr), (ax_ul, ax_ll), field_specs, cbars


def _measure_legend_band(ax, dpi):
    """Occupied vertical span (inches, + small pad) from ``ax``'s bottom edge down to
    the BOTTOM of the legend hung below it — i.e. the anchor drop PLUS the legend's own
    height, not the height alone. Sizing the row's leg-band to this full span keeps the
    panel below uncrowded (the legend can't reach the lower row's title even when the
    anchor drops it clear of tall x-tick labels) while leaving the right panel free of
    dead space (Norm: no extra whitespace on the right pair). Band-independent — the
    drop and the legend height are both fixed across the two layout passes — so pass 1's
    measurement is exact for pass 2."""
    leg = ax.get_legend()
    if leg is None:
        return 0.20
    try:
        span_px = ax.get_window_extent().y0 - leg.get_window_extent().y0
        return span_px / dpi + 0.10
    except Exception:
        return 0.60


def render_figure(bench, sd, fem_data, field, failure=None, fs=None,
                  out_dir=OUT, dpi=150):
    """Render the 2×2 composite as ONE figure and save <out_dir>/<bench>.png.

    Pure plotting: no solve happens here, so a cached (sd, fem_data, field, failure)
    can be re-rendered cheaply. Two passes: the first draws to MEASURE the two left
    legends' heights, the second lays the figure out with rows sized to them (so the
    inter-panel spacing is content-derived, not a fixed grid stretch). The four plot
    rectangles are pixel-identical by construction; render_figure asserts it before
    saving (see _verify_uniform).

    ``failure`` (result['failure_solution']) drives the strain + vector panels — the
    AT-FAILURE mechanism, with ``fs`` in the "… at Failure  FS = X" titles; None
    falls back to ``field``.
    """
    os.makedirs(out_dir, exist_ok=True)
    style = resolve_style(None)
    afield = _at_failure_field(field, failure, fs)
    # One shared padded domain for all four panels — mesh bbox, widened to contain any
    # above-surface / past-toe inputs chrome (see _composite_domain). Computed once so
    # both layout passes and the uniformity gate use the identical domain.
    domain = _composite_domain(fem_data, sd, style)

    with plt.rc_context(_RC):
        # Pass 1: provisional bands → measure the real legend heights.
        fig, _axes, (ax_ul, ax_ll), _fs, _cb = _build_composite(
            bench, sd, fem_data, afield, style, 1.0, 1.0, dpi, domain)
        leg0 = _measure_legend_band(ax_ul, dpi)
        leg1 = _measure_legend_band(ax_ll, dpi)
        plt.close(fig)

        # Pass 2: final layout with measured leg-bands.
        fig, axes, _left, field_specs, cbars = _build_composite(
            bench, sd, fem_data, afield, style, leg0, leg1, dpi, domain)

        rects, cushion = _verify_uniform(fig, axes, domain)
        print(f'  [{bench}] axes px (x,y,w,h) UL/UR/LL/LR: {rects}  '
              f'min content cushion (data units): {cushion}', flush=True)

        out = os.path.join(out_dir, f'{bench}.png')
        fig.savefig(out, dpi=dpi)
        plt.close(fig)
    return out


def _build_inputs_only(bench, sd, style, leg_in, dpi, domain):
    """Build the single-panel INPUTS figure: the composite's upper-left panel, in the
    composite's own chrome (same margins, same _WD_IN data width, same legend band),
    with no other panel beside it. Used for a row the page documents as reaching no
    equilibrium — there is no failure mechanism to draw, so nothing is drawn in its
    place. Returns (fig, ax)."""
    X0, X1, Y0, Y1 = domain
    DW, DH = X1 - X0, Y1 - Y0
    aspect = DH / DW if DW > 0 else 0.4
    Hd = _WD_IN * aspect

    W = _ML + _WD_IN + _MR
    H = _MB + leg_in + Hd + _MT
    fig = plt.figure(figsize=(W, H), dpi=dpi)
    ax = fig.add_axes([_ML / W, (_MB + leg_in) / H, _WD_IN / W, Hd / H])

    in_h, in_l, _bg = _draw_inputs_panel(ax, sd, style)
    ax.set_xlim(X0, X1)
    ax.set_ylim(Y0, Y1)
    ax.set_aspect('equal', adjustable='box')
    ax.set_axisbelow(True)
    ax.set_title('FEM Inputs', fontsize=_TITLE_FS, pad=_TITLE_PAD)

    fig.canvas.draw()
    _rr = fig.canvas.get_renderer()
    if in_h:
        ext = ax.get_window_extent()
        labels = [t for t in ax.get_xticklabels() if t.get_text()]
        band_px = (ext.y0 - min(t.get_window_extent(_rr).y0 for t in labels)
                   if labels else 0.0)
        anchor = (0.5, -(band_px + _LEG_CLEAR / 72.0 * fig.dpi) / ext.height)
        ncol = _fit_ncol(ax, fig, in_h, in_l, anchor, _LEG_FS)
        ax.legend(in_h, in_l, loc='upper center', bbox_to_anchor=anchor, ncol=ncol,
                  frameon=False, fontsize=_LEG_FS, handlelength=1.6,
                  columnspacing=1.2, handletextpad=0.5)
    return fig, ax


def render_inputs_figure(bench, sd, fem_data, out_dir=OUT, dpi=150):
    """Render the inputs-only figure and save <out_dir>/<bench>.png. No solve: the
    mesh is built only to fix the SAME padded domain the family's composites use, so
    an unlocked row's frame matches its locked siblings exactly. Two passes, like
    render_figure: measure the legend band, then lay out against it."""
    os.makedirs(out_dir, exist_ok=True)
    style = resolve_style(None)
    domain = _composite_domain(fem_data, sd, style)

    with plt.rc_context(_RC):
        fig, ax = _build_inputs_only(bench, sd, style, 1.0, dpi, domain)
        leg = _measure_legend_band(ax, dpi)
        plt.close(fig)

        fig, ax = _build_inputs_only(bench, sd, style, leg, dpi, domain)
        rects, cushion = _verify_uniform(fig, (ax,), domain)
        print(f'  [{bench}] axes px (x,y,w,h): {rects}  '
              f'min content cushion (data units): {cushion}', flush=True)

        out = os.path.join(out_dir, f'{bench}.png')
        fig.savefig(out, dpi=dpi)
        plt.close(fig)
    return out


def _verify_uniform(fig, axes, domain, wtol=2.0, htol=2.0):
    """Acceptance, MEASURED not eyeballed: the four data axes must be equal in pixel
    width and height (within tol), and every panel's content must sit strictly inside
    its frame (cushion present). Raises AssertionError with the numbers on failure;
    returns the per-axes (x, y, w, h) pixel rectangles for the report."""
    fig.canvas.draw()
    rects = [ax.get_window_extent() for ax in axes]
    ws = [r.width for r in rects]
    hs = [r.height for r in rects]
    dw, dh = max(ws) - min(ws), max(hs) - min(hs)
    assert dw <= wtol, f'panel widths differ by {dw:.2f}px (>{wtol}); widths={[round(w,1) for w in ws]}'
    assert dh <= htol, f'panel heights differ by {dh:.2f}px (>{htol}); heights={[round(h,1) for h in hs]}'
    # Cushion: content dataLim strictly inside the imposed viewLim on all sides.
    X0, X1, Y0, Y1 = domain
    min_cushion = None
    for ax in axes:
        dl = ax.dataLim
        left = dl.x0 - X0
        right = X1 - dl.x1
        bot = dl.y0 - Y0
        top = Y1 - dl.y1
        m = min(left, right, bot, top)
        min_cushion = m if min_cushion is None else min(min_cushion, m)
    assert min_cushion is not None and min_cushion > 0, \
        f'content touches a frame edge (min cushion {min_cushion}); expected > 0'
    return [(round(r.x0, 1), round(r.y0, 1), round(r.width, 1), round(r.height, 1))
            for r in rects], round(min_cushion, 4)


def make_figure(tag, dpi=150):
    """Solve the SSRM bracket, export the sidecars, render the composite.

    ``figure='inputs'`` short-circuits the solve entirely and renders the
    inputs-only figure (returns fs=None) — the mode for a row with no equilibrium
    to draw."""
    bench = tag.get('benchmark', os.path.basename(tag['file']).split('.')[0])
    if tag.get('figure') == 'inputs':
        sd, fem_data, _path, _mesh = _build(tag)
        return render_inputs_figure(bench, sd, fem_data, dpi=dpi), None
    sd, fem_data, field, failure, fs, path, mesh = build_and_solve(tag)

    # Sidecars next to the case xlsx (Norm directive): the converged field plus,
    # when captured, the at-failure mechanism, so every future re-render is
    # solve-free. Stem is the xlsx path without extension → {stem}_fem_*.csv,
    # except on a second run of a shared file (see SIDECAR_STEM).
    stem = _sidecar_stem(tag, path)
    meta = {'benchmark': bench, 'analysis': 'ssrm', 'FS': float(fs),
            'expected_fs': tag.get('expected_fs'), 'file': tag.get('file')}
    with contextlib.redirect_stdout(io.StringIO()):
        export_fem_solution(fem_data, field, stem, meta=meta,
                            failure_solution=failure)
        # The mesh the fields were solved on, beside them. Without it a reload
        # discovers whatever {stem}_mesh.json happens to sit there — or, for the
        # rows that ship none, rebuilds from the tag and hopes the discretization
        # comes out node-for-node the same, so everything downstream reads a
        # section that was never solved.
        export_mesh_to_json(mesh, f'{stem}_mesh.json')

    out = render_figure(bench, sd, fem_data, field, failure=failure, fs=fs, dpi=dpi)
    return out, fs


def make_figure_from_sidecar(tag, dpi=150):
    """Re-render the composite SOLVE-FREE from the committed sidecars. Rebuilds the
    mesh (no solve) and imports the converged + at-failure fields written by
    export_fem_solution, so a layout change can be re-rendered without touching the
    solver. Returns (out_path, fs) where fs comes from the run-meta sidecar."""
    from xslope.fem import import_fem_solution, import_fem_meta
    bench = tag.get('benchmark', os.path.basename(tag['file']).split('.')[0])
    if tag.get('figure') == 'inputs':
        sd, fem_data, _path, _mesh = _build(tag)
        return render_inputs_figure(bench, sd, fem_data, dpi=dpi), None
    sd, fem_data, path, _mesh = _build(tag)
    # The inputs panel must show the same constraint the solved field was produced
    # under, on this path as on the solving one — otherwise a --from-sidecar re-render
    # of a constrained row would draw an unconstrained model beside its confined
    # mechanism.
    sd = _record_ssr_zone(sd, _tag_ssr_zone(tag), fem_data)
    stem = _sidecar_stem(tag, path)
    with contextlib.redirect_stdout(io.StringIO()):
        solution = import_fem_solution(fem_data, stem)
    failure = solution.get('failure_solution')
    meta = import_fem_meta(stem) or {}
    fs = meta.get('FS', tag.get('expected_fs'))
    fs = float(fs) if fs is not None else None
    # The at-failure title needs the FS marker; import carries it via meta, not the
    # reconstructed field — thread it through render_figure's fs arg.
    out = render_figure(bench, sd, fem_data, solution, failure=failure, fs=fs, dpi=dpi)
    return out, fs


def audit(out_dir=OUT, verbose=True):
    """Coverage check: which registered rows have no rendered PNG. Every fem_ssrm
    tag on the page is a figure this script is supposed to produce, so a locked row
    with no <benchmark>.png is a silent coverage hole (the whole RS2-48–55 wall
    family was one). Reported-not-locked rows in EXTRA_CASES are checked too, listed
    separately.

    A row named in ``EXPECTED_NO_FIGURE`` is not a hole — it is a row whose mechanism
    the page already figures elsewhere, adjudicated once and given a reason. Those are
    listed with their reasons rather than counted as missing, and the exemption list
    is itself checked: an entry for a row that is not registered, or for a row that
    now HAS a figure, is a DEAD exemption and counts as a failure, so the list cannot
    quietly outlive the text it describes.

    Returns ``(missing_locked, missing_reported, dead)``. The audit is clean when all
    three are empty; ``__main__`` exits non-zero otherwise."""
    tags = parse_tags()
    cases = tags + EXTRA_CASES
    registered = {}
    missing_locked, missing_reported, exempt = [], [], []
    for tag in cases:
        bench = tag.get('benchmark', os.path.basename(tag['file']).split('.')[0])
        registered[bench] = os.path.exists(os.path.join(out_dir, f'{bench}.png'))
        if registered[bench]:
            continue
        if bench in EXPECTED_NO_FIGURE:
            exempt.append((bench, EXPECTED_NO_FIGURE[bench]))
            continue
        (missing_locked if tag.get('expected_fs') else missing_reported).append(
            (bench, tag.get('file', '?')))

    dead = []
    for bench in EXPECTED_NO_FIGURE:
        if bench not in registered:
            dead.append((bench, 'not a registered row'))
        elif registered[bench]:
            dead.append((bench, 'row now HAS a figure'))

    if verbose:
        print(f'{len(cases)} registered rows ({len(tags)} tagged, '
              f'{len(EXTRA_CASES)} reported-only); figures in '
              f'{os.path.normpath(out_dir)}')
        for label, rows in (('locked, NO figure', missing_locked),
                            ('reported, no figure', missing_reported)):
            print(f'  {label}: {len(rows)}')
            for bench, f in rows:
                print(f'    {bench:22s} {f}')
        print(f'  expected no figure (named): {len(exempt)}')
        for bench, why in exempt:
            print(f'    {bench:22s} {why}')
        print(f'  DEAD exemptions: {len(dead)}')
        for bench, why in dead:
            print(f'    {bench:22s} {why}')
        total = len(missing_locked) + len(missing_reported) + len(dead)
        print(f'  unexplained rows: {total}')
    return missing_locked, missing_reported, dead


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    args = sys.argv[1:]
    if '--audit' in args:
        ml, mr, dead = audit()
        sys.exit(1 if (ml or mr or dead) else 0)
    # --from-sidecar re-renders SOLVE-FREE from the committed sidecars (no solver).
    from_sidecar = '--from-sidecar' in args
    only = set(a for a in args if not a.startswith('--'))
    cases = parse_tags() + EXTRA_CASES
    print(f'{len(cases)} registered rows (fem_ssrm tags + EXTRA_CASES)'
          f"{'  (solve-free re-render from sidecars)' if from_sidecar else ''}")
    for tag in cases:
        bench = tag.get('benchmark', '?')
        if only and bench not in only:
            continue
        t0 = time.time()
        try:
            out, fs = (make_figure_from_sidecar(tag) if from_sidecar
                       else make_figure(tag))
            exp = tag.get('expected_fs')
            fsx = ('inputs-only' if tag.get('figure') == 'inputs'
                   else f'{fs:.3f}' if fs is not None else 'n/a')
            lock = (f'  lock={float(exp):.3f} d={fs-float(exp):+.3f}'
                    if exp and fs is not None else '')
            print(f'ok   {bench:10s} FS={fsx}{lock}  ({time.time()-t0:.0f}s)  '
                  f'{os.path.basename(out)}', flush=True)
        except Exception as e:
            import traceback
            print(f'FAIL {bench:10s} {type(e).__name__}: {e}', flush=True)
            traceback.print_exc()
