"""Build the earth-dam sample workbooks under docs/lem/files from one definition
of the dam.

  xslope_earth_dam_up         docs/lem/samples.md, problem 4 (upstream slope,
                              steady piezometric line)
  xslope_earth_dam_down       docs/lem/samples.md, problem 4 (downstream slope,
                              the book's piezometric line)
  xslope_earth_dam_down_seep  docs/lem/samples.md, problem 4 (downstream slope,
                              pore pressures from the seepage solution)
  xslope_gsat_seep            docs/lem/samples.md, problem 6 (gamma_sat with the
                              water table from a steady seepage solution)
  xslope_gsat_rapid           docs/lem/samples.md, problem 6 (gamma_sat under
                              rapid drawdown, staged piezometric lines; the same
                              model as xslope_rapid below)
  xslope_earth_dam_rapid      the rapid-drawdown model with two seepage stages
                              (seep-bc and seep-bc (2)); a fixture of the seepage
                              and report checks and of the Studio screenshots

and, under docs/inputs/slope (the driver scripts' download set, and fixtures of
the regression suite and of Studio's screenshots and checks):

  xslope_dam                  the steady upstream model (problem 4's
                              xslope_earth_dam_up, under the name the driver
                              scripts use)
  xslope_rapid                the staged-piezometric-line rapid-drawdown model

The dam is Example 2 of Section 6.13.6 of Duncan, Wright & Brandon, Soil
Strength and Slope Stability, 2nd ed. (2014), pp. 121-123: a shell and a clay
core on a clay layer over sand, crest at El. 317, reservoir at El. 302,
foundation surface at El. 227, clay to El. 197 and sand to El. 182, the base of
the model. The book analyzes the downstream slope with pore pressures from a
finite-element seepage solution and from a piezometric line taken as that
solution's zero-pressure line (its Figure 6.42), each for the overall critical
circle and for circles tangent to El. 197 (its Table 6.8).

The book states the elevations, the 3:1 faces, the 1H:2V core faces, the
strengths, one unit weight per zone and the permeabilities (ft/min, entered here
as printed with the time unit set to minutes). The rest is scaled from Figures
6.41 and 6.42 at their own scale bar, to the nearest foot: the downstream toe at
x = 581 with the upstream toe at 0, so a 41 ft crest; the core 237 to 345 at its
base and 278 to 304 at its top, El. 309; and the ends of the section at x = -158
and 638. The zero-pressure line is traced from Figure 6.42 to the nearest half
foot (_PIEZO_FULL).

The two downstream models are the book's analysis and carry its one unit weight
per zone. The others are the samples page's own analyses of the same dam (the
upstream slope and rapid drawdown, which the book does not analyze here), and
each zone there carries a saturated unit weight, its moist unit weight plus
5 pcf, for the part of it below the water table.

The rapid-drawdown strength parameters d and psi (the Duncan, Wright & Wong
total-stress envelope of a low-permeability material) are written only into the
rapid-drawdown runs, with the values of the course problem the model comes
from (CE 544, unit 2, rapid drawdown homework): core d = 300 psf, psi = 20 deg;
clay d = 100 psf, psi = 18 deg; shell and sand blank, as free-draining
materials. A blank d or psi loads as 0, and a slice whose material has d = 0 and
psi = 0 is treated as drained by the rapid-drawdown solver, which is how the
shell and sand have always been analyzed. The steady runs leave d and psi
blank on every material: nothing in a steady analysis reads them. The
drawn-down piezometric line (El. 250) is the course problem's, carried through
the core faces at the same elevations as before; downstream of the core it is
the full-pool line.

The three ``u = seep`` workbooks share one section, so they share one seepage
mesh: ``--meshes`` writes it beside each of them (``{base}_mesh.json``,
linear triangles at MESH_SIZE). The solutions (``{base}_seep.csv``,
``{base}_seep2.csv``) are solved on that mesh by tools/make_seep_sidecars.py.
A plain run never touches the meshes.

From 2026-10-18 the section, the zero-pressure line and the permeabilities are
the book's, as above. They replace an earlier reading of the section (downstream
toe at x = 590, a 50 ft crest, the core top at El. 307, ends at x = -150 and
740), a piezometric line drawn through that core, and permeabilities entered as
ft/s values under a time unit of days. The downstream models carry the book's one
unit weight; the others keep the saturated unit weights added on 2026-10-17.
Each model is written through :func:`xslope.fileio.save_slope_data_to_xlsx`,
which copies the current template and never writes a formula cell.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/build_lem_samples.py
                         PYTHONPATH=. python3 benchmarks/build_lem_samples.py --meshes
"""
import os
import sys

from xslope.fileio import save_slope_data_to_xlsx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'docs', 'lem', 'files')
#: The downloadable driver-script inputs (docs/inputs/slope), a second output
#: directory so verify_rebuild can check the two sets against their own folders.
OUT_INPUTS = os.path.join(ROOT, 'docs', 'inputs', 'slope')

_EMPTY_BC = {'specified_heads': [], 'specified_fluxes': [], 'exit_face': []}

# =============================================================================
# The dam
# =============================================================================
_X_LEFT, _X_RIGHT = -158.0, 638.0

#: Profile lines in material order: shell, core, clay, sand.
_PROFILE = [
    ([(0.0, 227.0), (270.0, 317.0), (311.0, 317.0), (581.0, 227.0)], 0),
    ([(237.0, 227.0), (278.0, 309.0), (304.0, 309.0), (345.0, 227.0)], 1),
    ([(_X_LEFT, 227.0), (_X_RIGHT, 227.0)], 2),
    ([(_X_LEFT, 197.0), (_X_RIGHT, 197.0)], 3),
]
_MAX_DEPTH = 182.0

#: name, gamma (pcf), c' (psf), phi' (deg), gamma_sat (pcf), k (ft/min)
_MATERIALS = [
    ('Shell', 125.0, 0.0, 34.0, 130.0, 1e-2),
    ('Core', 122.0, 100.0, 26.0, 127.0, 1e-6),
    ('Clay', 123.0, 0.0, 24.0, 128.0, 1e-5),
    ('Sand', 127.0, 0.0, 32.0, 132.0, 1e-3),
]

#: Rapid-drawdown total-stress envelope (d psf, psi deg) of the low-permeability
#: materials; every other material is free-draining (d and psi blank).
_RAPID_D_PSI = {'Core': (300.0, 20.0), 'Clay': (100.0, 18.0)}

#: Reservoir at El. 302: the book's zero-pressure line (Figure 6.42). It leaves
#: the pool at the core's upstream face, leaves the core near its downstream toe,
#: and meets the downstream face 1.5 ft above the toe.
_PIEZO_FULL = [(_X_LEFT, 302.0), (274.5, 302.0), (284.0, 294.5), (306.0, 279.5),
               (312.0, 274.0), (332.0, 243.5), (338.0, 235.5), (344.0, 232.5),
               (470.0, 232.0), (576.5, 228.5), (581.0, 227.0), (_X_RIGHT, 227.0)]
#: Reservoir at El. 250 after drawdown: the course problem's line, through the
#: core's upstream face at El. 250 and its downstream face at El. 233. Downstream
#: of the core it is the full-pool line, which the drawdown upstream does not move.
_PIEZO_DRAWN = [(_X_LEFT, 250.0), (248.5, 250.0), (311.0, 242.0),
                (342.0, 233.0)] + _PIEZO_FULL[7:]

#: Seepage boundaries: the reservoir face as a specified head, the downstream
#: face and ground as a potential exit face.
_EXIT_FACE = [(311.0, 317.0), (581.0, 227.0), (_X_RIGHT, 227.0)]
_SEEP_BC_FULL = {
    'specified_heads': [{'head': 302.0, 'kind': 'head',
                         'coords': [(_X_LEFT, 227.0), (0.0, 227.0), (225.0, 302.0)]}],
    'specified_fluxes': [], 'exit_face': list(_EXIT_FACE)}
_SEEP_BC_DRAWN = {
    'specified_heads': [{'head': 250.0, 'kind': 'head',
                         'coords': [(_X_LEFT, 227.0), (0.0, 227.0), (69.0, 250.0)]}],
    'specified_fluxes': [], 'exit_face': list(_EXIT_FACE)}

#: Starting circles (Depth option): one tangent to the top of the sand, one to
#: the base. Upstream and downstream differ only in the centre's x.
_YO = 500.0
_DEPTHS = (197.0, 182.0)
_XO_UP, _XO_DOWN = 100.0, 490.0

#: Element size (ft) of the seepage mesh the three u = seep models share.
MESH_SIZE = 3.5
_SEEP_MODELS = ('xslope_earth_dam_down_seep', 'xslope_gsat_seep',
                'xslope_earth_dam_rapid')


def _material(name, gamma, c, phi, gamma_sat, k, *, u, rapid, seep_props,
              E=700000.0, nu=0.3, t_cut=0.0):
    d, psi = _RAPID_D_PSI.get(name, (None, None)) if rapid else (None, None)
    return {
        'name': name, 'gamma': gamma, 'gamma_sat': gamma_sat, 'option': 'mc',
        'c': c, 'phi': phi, 'cp': 0.0, 'r_elev': 0.0, 'd': d, 'psi': psi,
        't_cut': t_cut, 'phi_b': None, 's_cap': None, 'Ss': None, 'Sy': None,
        'pow_a': 0.0, 'pow_b': 0.0, 'pow_c': 0.0, 'pow_d': 0.0,
        'u': u, 'ru': 0.0,
        'sigma_gamma': 0.0, 'sigma_c': 0.0, 'sigma_phi': 0.0,
        'sigma_cp': 0.0, 'sigma_d': 0.0, 'sigma_psi': 0.0,
        'k1': k if seep_props else 0.0, 'k2': k if seep_props else 0.0,
        'alpha': 0.0, 'unsat': 'lf',
        'kr0': 0.0001 if seep_props else 0.0, 'h0': -1.0 if seep_props else 0.0,
        'vg_a': 0.0, 'vg_n': 0.0, 'vg_l': 0.5,
        'E': E, 'nu': nu,
        'hb_sci': 0.0, 'hb_gsi': 0.0, 'hb_mi': 0.0, 'hb_d': 0.0,
    }


def _dam(*, u, rapid, seep_props, xo=_XO_UP, book_weights=False, **mat_extra):
    """The dam with the water and analysis options of one run.

    u             the pore-pressure option of every material ('piezo' or 'seep')
    rapid         write the rapid-drawdown d and psi (only a rapid run reads them)
    seep_props    write the permeabilities and the linear-front unsaturated law
                  (the runs that solve their own seepage), with minutes as the
                  time unit of the book's ft/min
    book_weights  the book's one unit weight per zone, with no saturated unit
                  weight (the downstream models, which are the book's analysis)
    """
    sd = {
        'unit_system': 'imperial',
        'gamma_water': 62.4,
        'tcrack_depth': 0.0,
        'tcrack_water': 0.0,
        'k_seismic': 0.0,
        'water_loads': 'auto',
        'materials': [_material(n, g, c, f, None if book_weights else gs, k, u=u,
                                rapid=rapid, seep_props=seep_props, **mat_extra)
                      for n, g, c, f, gs, k in _MATERIALS],
        'profile_lines': [{'coords': list(c), 'mat_id': m, 'size': None}
                          for c, m in _PROFILE],
        'max_depth': _MAX_DEPTH,
        'piezo_line': [],
        'piezo_line2': [],
        'piezo_phreatic': False,
        'piezo_phreatic2': False,
        'circles': [{'Xo': xo, 'Yo': _YO, 'Depth': dep, 'R': _YO - dep}
                    for dep in _DEPTHS],
        'non_circ': [],
        'seepage_bc': {'specified_heads': [dict(h, coords=list(h['coords']))
                                           for h in _SEEP_BC_FULL['specified_heads']],
                       'specified_fluxes': [],
                       'exit_face': list(_EXIT_FACE)},
        'seepage_bc2': dict(_EMPTY_BC),
        'has_seepage_bc2': False,
    }
    if seep_props:
        sd['time_unit'] = 'min'
    return sd


def _write(sd, name, out=None):
    dst = os.path.join(out or OUT, name)
    save_slope_data_to_xlsx(sd, dst)
    return dst


# =============================================================================
# docs/lem/samples.md, problem 4: the upstream slope, and the book's downstream
# analysis with the piezometric line and with the seepage solution
# =============================================================================
def build_earth_dam_up():
    sd = _dam(u='piezo', rapid=False, seep_props=False)
    sd['piezo_line'] = list(_PIEZO_FULL)
    return _write(sd, 'xslope_earth_dam_up.xlsx')


def build_earth_dam_down():
    sd = _dam(u='piezo', rapid=False, seep_props=False, xo=_XO_DOWN,
              book_weights=True)
    sd['piezo_line'] = list(_PIEZO_FULL)
    return _write(sd, 'xslope_earth_dam_down.xlsx')


def build_earth_dam_down_seep():
    sd = _dam(u='seep', rapid=False, seep_props=True, xo=_XO_DOWN,
              book_weights=True)
    return _write(sd, 'xslope_earth_dam_down_seep.xlsx')


# =============================================================================
# docs/lem/samples.md, problem 6: gamma_sat with a seepage water table
# =============================================================================
def build_gsat_seep():
    sd = _dam(u='seep', rapid=False, seep_props=True)
    return _write(sd, 'xslope_gsat_seep.xlsx')


# =============================================================================
# docs/lem/samples.md, problem 6: gamma_sat under rapid drawdown
# =============================================================================
def build_gsat_rapid():
    """Staged piezometric lines (pool at El. 302, then El. 250). E and nu are
    entered as 0 and t_cut is blank, as the shipped file has them: this file is
    solved by limit equilibrium only, which reads none of the three."""
    sd = _dam(u='piezo', rapid=True, seep_props=False,
              E=0.0, nu=0.0, t_cut=None)
    sd['piezo_line'] = list(_PIEZO_FULL)
    sd['piezo_line2'] = list(_PIEZO_DRAWN)
    return _write(sd, 'xslope_gsat_rapid.xlsx')


# =============================================================================
# Rapid drawdown from two seepage stages
# =============================================================================
def build_earth_dam_rapid():
    sd = _dam(u='seep', rapid=True, seep_props=True)
    sd['seepage_bc2'] = {
        'specified_heads': [dict(h, coords=list(h['coords']))
                            for h in _SEEP_BC_DRAWN['specified_heads']],
        'specified_fluxes': [], 'exit_face': list(_EXIT_FACE)}
    sd['has_seepage_bc2'] = True
    return _write(sd, 'xslope_earth_dam_rapid.xlsx')


# =============================================================================
# docs/inputs/slope: the dam as the driver scripts and the test fixtures load it
# =============================================================================
def build_input_dam():
    """The steady upstream model (problem 4's xslope_earth_dam_up) under the
    download set's name. As shipped, the sand row's E and nu are entered as 0
    where every other row carries 700,000 and 0.3; kept as found and reported,
    since only limit equilibrium is run on this file."""
    sd = _dam(u='piezo', rapid=False, seep_props=False)
    sd['piezo_line'] = list(_PIEZO_FULL)
    sand = sd['materials'][3]
    sand['E'], sand['nu'] = 0.0, 0.0
    return _write(sd, 'xslope_dam.xlsx', OUT_INPUTS)


def build_input_rapid():
    """The staged-piezometric-line rapid-drawdown model, the same model as
    problem 6's rapid run (xslope_gsat_rapid). E and nu entered as 0 and t_cut
    blank, as shipped."""
    sd = _dam(u='piezo', rapid=True, seep_props=False,
              E=0.0, nu=0.0, t_cut=None)
    sd['piezo_line'] = list(_PIEZO_FULL)
    sd['piezo_line2'] = list(_PIEZO_DRAWN)
    return _write(sd, 'xslope_rapid.xlsx', OUT_INPUTS)


BUILDERS = [
    build_earth_dam_up,
    build_earth_dam_down,
    build_earth_dam_down_seep,
    build_gsat_seep,
    build_gsat_rapid,
    build_earth_dam_rapid,
]

INPUT_BUILDERS = [
    build_input_dam,
    build_input_rapid,
]


def build_meshes():
    """Mesh the section once and write the mesh beside each u = seep model."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import (build_mesh_from_polygons, export_mesh_to_json,
                             get_material_polygons)
    sd = load_slope_data(os.path.join(OUT, _SEEP_MODELS[0] + '.xlsx'))
    mesh = build_mesh_from_polygons(get_material_polygons(sd), MESH_SIZE, 'tri3')
    for stem in _SEEP_MODELS:
        dst = os.path.join(OUT, stem + '_mesh.json')
        export_mesh_to_json(mesh, dst)
        print('meshed', dst, len(mesh['nodes']), 'nodes')


def main(argv=()):
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(OUT_INPUTS, exist_ok=True)
    for fn in BUILDERS + INPUT_BUILDERS:
        print('built', fn())
    if '--meshes' in argv:
        build_meshes()


if __name__ == '__main__':
    main(sys.argv[1:])
