"""Which piezometric lines the input and solution plots draw (``xslope.plot``).

A piezometric line is an input sheet a model may carry without the analysis
consuming it, and a line on the section that nothing reads is a statement about
the model that the run does not honor. The plots therefore draw a line only where
it is READ, and there are exactly two ways it is:

  * **pore pressure** -- some material's ``u`` option is ``piezo``, which takes the
    pressure on a slice base from the line;
  * **the weight split** -- a material carrying ``gamma_sat`` weighs its slice
    saturated below the water table and moist above it, and that water table is
    the piezometric line wherever no seepage solution supplies one;
  * **the water load** -- the line is what states the pool standing on the ground.
    Which sheet states that pool is ``xslope.water.water_line_for_stage``'s
    precedence: the seepage head boundaries wherever a seepage analysis is
    defined, otherwise the stage's piezometric line. So a model whose materials
    all read a seepage solution still shows its piezometric line when no boundary
    set was defined and the line is what puts water on the slope -- and a model
    that does define head boundaries does not, because they are the source and the
    line is inert.

This module pins the rule from both ends -- the predicate
(:func:`xslope.plot.piezo_line_used`) and what actually reaches the figure, read
off the rendered legend, which is where a reader learns a line is there:

  * a piezo-option model draws its line, and both lines of a rapid-drawdown deck;
  * the same model with its materials moved off ``piezo``, one unit weight each
    and no pool on the ground draws nothing;
  * a seep-option model with no boundary set still draws the line the pool comes
    from, and the same model WITH its boundary set does not;
  * a model whose material carries gamma_sat draws the line its weight split is
    measured from, and stops when that unit weight is removed;
  * a model that defines no piezometric line at all is unaffected either way;
  * ``only_if_used=False`` still draws whatever the model defines, which is what
    an editor previewing the sheet being typed needs.

File-light: five shipped workbooks, one 20-slice solve, no search and no seepage
run.

Run directly:  PYTHONPATH=. python3 test/piezo_visibility_check.py
"""

import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from xslope.fileio import load_slope_data
from xslope.plot import plot_inputs, plot_solution, piezo_line_used
from xslope.slice import generate_slices
from xslope.solve import solve_selected

ROOT = Path(__file__).resolve().parent.parent / 'docs'

#: A layered slope whose materials take pore pressure from the line, and whose
#: line runs INSIDE the slope -- it loads nothing, so the u option is the only
#: thing that reads it.
PIEZO = ROOT / 'lem' / 'files' / 'xslope_method_slices_problem.xlsx'
#: A rapid-drawdown deck stated with two piezometric lines: full pool and
#: drawn-down pool, both standing on the upstream face.
RAPID = ROOT / 'tutorials' / 'files' / 'xslope_johnson_rapid_start.xlsx'
#: The same reservoir solved by seepage instead: materials on u = seep, both
#: boundary sets defined, no piezometric line anywhere.
SEEPED = ROOT / 'tutorials' / 'files' / 'xslope_johnson_rapid.xlsx'
#: An embankment against a reservoir: a piezometric line that DOES stand above
#: the ground, with head boundaries defined beside it.
POOL = ROOT / 'inputs' / 'slope' / 'xslope_dam.xlsx'
#: A slope whose one material carries a saturated unit weight: u = none, no pool,
#: and the line is what splits each slice's weight -- moist above, saturated below.
GSAT = ROOT / 'lem' / 'files' / 'xslope_gsat_sidecar.xlsx'

LABEL = {1: 'Piezometric Line', 2: 'Piezometric Line 2'}


def _model(path, u=None, drop_bc=False, drop_gamma_sat=False):
    """A shipped model, optionally with every material moved to one ``u`` option,
    its saturated unit weights removed, and its seepage boundary sets deleted."""
    sd = load_slope_data(str(path))
    for m in sd['materials']:
        if u is not None:
            m['u'] = u
        if drop_gamma_sat:
            m['gamma_sat'] = None
    if drop_bc:
        sd['seepage_bc'] = {}
        sd['seepage_bc2'] = {}
    return sd


def _legend_labels(fig):
    """Every label the figure offers a reader, legend entries included."""
    out = set()
    for ax in fig.axes:
        out.update(ax.get_legend_handles_labels()[1])
        leg = ax.get_legend()
        if leg is not None:
            out.update(t.get_text() for t in leg.get_texts())
    return out


def _drawn_on_inputs(sd, **kw):
    fig = plot_inputs(sd, **kw)
    try:
        return {s for s in _legend_labels(fig) if s.startswith('Piezometric')}
    finally:
        plt.close(fig)


def _drawn_on_solution(sd):
    circle = sd['circles'][0] if sd.get('circular') else None
    ok, res = generate_slices(sd, circle=circle,
                              non_circ=None if circle else sd.get('non_circ'),
                              num_slices=20)
    if not ok:
        raise RuntimeError(f"could not slice the model: {res}")
    slice_df, failure_surface = res
    results = solve_selected('bishop', slice_df, rapid=False)
    fig = plot_solution(sd, slice_df, failure_surface, results)
    try:
        return {s for s in _legend_labels(fig) if s.startswith('Piezometric')}
    finally:
        plt.close(fig)


def _expect(failures, what, got, want):
    if got != want:
        failures.append(f"{what}: drew {sorted(got) or 'nothing'}, "
                        f"expected {sorted(want) or 'nothing'}")


# ---------------------------------------------------------------------------
# The predicate
# ---------------------------------------------------------------------------

def check_used_pore_pressure(failures):
    """A material on u = piezo reads the line, whatever the water loads do."""
    sd = _model(PIEZO)
    if not piezo_line_used(sd, 1):
        failures.append("u = piezo does not mark Line 1 as read")
    sd = _model(RAPID)
    for stage in (1, 2):
        if not piezo_line_used(sd, stage):
            failures.append(f"u = piezo does not mark Line {stage} as read "
                            f"on a rapid-drawdown deck")


def check_used_water_load(failures):
    """With no material reading the line, the pool it states still does."""
    sd = _model(POOL, u='seep', drop_bc=True)
    if not piezo_line_used(sd, 1):
        failures.append("a piezometric line that is the model's only water-load "
                        "source is not marked as read")
    # The same model with its head boundaries back: they are the source, by the
    # precedence in water_line_for_stage, and the line is inert.
    sd = _model(POOL, u='seep')
    if piezo_line_used(sd, 1):
        failures.append("a piezometric line is marked as read while the seepage "
                        "head boundaries are the water-load source")


def check_used_weight_split(failures):
    """A material carrying gamma_sat weighs its slice against the line."""
    sd = _model(GSAT)
    if not piezo_line_used(sd, 1):
        failures.append("the water table a gamma_sat weight split is measured "
                        "from is not marked as read")
    # Strip the saturated unit weight and nothing reads the line any more.
    for m in sd['materials']:
        m['gamma_sat'] = None
    if piezo_line_used(sd, 1):
        failures.append("the line is marked as read on a model with no saturated "
                        "unit weight, no piezo option and no water load")


def check_unused_line(failures):
    """u = none, one unit weight, no pool: nothing in the run reads the line."""
    for u in ('none', 'ru'):
        if piezo_line_used(_model(PIEZO, u=u, drop_gamma_sat=True), 1):
            failures.append(f"u = {u} with one unit weight and no water load "
                            f"marks the line as read")
    # The u option alone does not settle it: put the saturated unit weights back
    # and the same model reads the line again, for the weight split.
    if not piezo_line_used(_model(PIEZO, u='none'), 1):
        failures.append("moving the materials off u = piezo hid a line the "
                        "gamma_sat weight split still measures from")


def check_no_line_at_all(failures):
    """A model that defines no piezometric line answers False without deriving."""
    sd = _model(SEEPED)
    for stage in (1, 2):
        if piezo_line_used(sd, stage):
            failures.append(f"Line {stage} is marked as read on a model that "
                            f"defines no piezometric line")


# ---------------------------------------------------------------------------
# What reaches the figure
# ---------------------------------------------------------------------------

def check_inputs_plot(failures):
    """The four cases as the input plot draws them."""
    _expect(failures, "inputs, u = piezo",
            _drawn_on_inputs(_model(PIEZO)), {LABEL[1]})
    _expect(failures, "inputs, rapid-drawdown pair on u = piezo",
            _drawn_on_inputs(_model(RAPID)), {LABEL[1], LABEL[2]})
    _expect(failures, "inputs, seepage solution with both boundary sets",
            _drawn_on_inputs(_model(SEEPED)), set())
    _expect(failures, "inputs, u = seep with no boundary set (line is the load "
                      "source)",
            _drawn_on_inputs(_model(POOL, u='seep', drop_bc=True)), {LABEL[1]})
    _expect(failures, "inputs, u = seep with head boundaries (they are the load "
                      "source)",
            _drawn_on_inputs(_model(POOL, u='seep')), set())
    _expect(failures, "inputs, u = none with one unit weight and no water load",
            _drawn_on_inputs(_model(PIEZO, u='none', drop_gamma_sat=True)), set())
    _expect(failures, "inputs, gamma_sat weight split",
            _drawn_on_inputs(_model(GSAT)), {LABEL[1]})


def check_solution_plot(failures):
    """The solution plot answers the same way as the input plot."""
    _expect(failures, "solution, u = piezo",
            _drawn_on_solution(_model(PIEZO)), {LABEL[1]})
    _expect(failures, "solution, u = none with one unit weight and no water load",
            _drawn_on_solution(_model(PIEZO, u='none', drop_gamma_sat=True)), set())


def check_editor_preview_opt_out(failures):
    """``only_if_used=False`` draws the sheet as typed, for an editor preview."""
    sd = _model(PIEZO, u='none', drop_gamma_sat=True)
    fig = plt.figure()
    try:
        from xslope.plot import plot_piezo_line
        ax = fig.add_subplot(111)
        plot_piezo_line(ax, sd, only_if_used=False)
        got = {s for s in ax.get_legend_handles_labels()[1]
               if s.startswith('Piezometric')}
    finally:
        plt.close(fig)
    _expect(failures, "only_if_used=False", got, {LABEL[1]})


CHECKS = [check_used_pore_pressure, check_used_water_load,
          check_used_weight_split, check_unused_line,
          check_no_line_at_all, check_inputs_plot, check_solution_plot,
          check_editor_preview_opt_out]


def run():
    """Returns a list of failure descriptions (empty on success)."""
    failures = []
    for chk in CHECKS:
        try:
            chk(failures)
        except Exception as e:                                 # noqa: BLE001
            failures.append(f"{chk.__name__} raised {type(e).__name__}: {e}")
    return failures


if __name__ == '__main__':
    import sys
    bad = run()
    for f in bad:
        print(f"FAIL  {f}")
    print("piezo_visibility_check: "
          + (f"{len(bad)} failure(s)" if bad else f"{len(CHECKS)} checks passed"))
    sys.exit(1 if bad else 0)
