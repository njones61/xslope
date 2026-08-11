"""The search window a sweep searches inside (``xslope.sensitivity``).

A model may declare a search window on its circles sheet -- entry and exit
ranges, a center box, a tangent-depth band, a minimum slip depth. It is how an
engineer says *this* is the mechanism I am studying, and Studio's Run LEM path
has always honoured it.

A parametric study is exactly where ignoring it does the most damage. A sweep
exists to answer "how much does the factor of safety move when this input
moves", and it answers by re-searching at every value. Search unconstrained and
a step is free to settle in a different surface FAMILY from its neighbour -- a
deep foundation circle where the last point found a shallow fill one -- and the
step change in FS between them is then read as sensitivity to the parameter when
it is really a change in what was measured. The same reasoning holds for every
study built on the sweep: a tornado compares bars against each other, a design
sweep interpolates a crossing along the curve, and neither means anything if the
points are not on one mechanism.

So the window is read from the model by default, everywhere a sweep searches,
exactly as the interface reads it. What this module pins:

  * a model that declares no window is searched unconstrained, as before -- the
    default must not invent a constraint nobody entered;
  * a declared window reaches the search at EVERY point, base case included;
  * a half-declared range (one end only) is not a window and is not forwarded;
  * an explicit ``search_opts`` from the caller beats the file, and
    ``use_file_window=False`` turns the file off entirely;
  * the studies built on the sweep -- ``design``, ``back_analysis``,
    ``tornado``, ``scaled_sensitivity`` -- all inherit it;
  * a NON-circular sweep takes the one limit its search understands
    (``min_slip_depth``) and drops the circle-shaped ones rather than raising;
  * and the PARITY that makes all of the above one behaviour rather than three:
    the same model reaches the same search with the same constraints whether it
    is swept, run through a reliability engine, or run from Studio's Run LEM
    button -- because all of them read the window through the single helper in
    ``xslope.search``, on both the circular and the non-circular branch.

File-light and fast: it loads the shipped ACADS sample for a real model, injects
the window as a dict, and never runs a search -- the search functions are
replaced by recorders, because what is under test is which constraints reach
them.

Run directly:  PYTHONPATH=. python3 test/sweep_window_check.py
"""

import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')

from xslope import search as SEARCH
from xslope import sensitivity as S
from xslope.fileio import load_slope_data

SAMPLE = Path(__file__).resolve().parent.parent / 'docs' / 'lem' / 'files' / \
    'xslope_acads_simple.xlsx'

#: A full window, one of every kind the circles sheet can declare.
FULL_WINDOW = {
    'entry_x_min': 10.0, 'entry_x_max': 30.0,
    'exit_x_min': 0.0, 'exit_x_max': 8.0,
    'center_box_x_min': 5.0, 'center_box_y_min': 20.0,
    'center_box_x_max': 40.0, 'center_box_y_max': 60.0,
    'max_tangent_depth': -5.0,
    'min_slip_depth': 1.5,
}

_WINDOW_KEYS = ('entry_range', 'exit_range', 'center_box', 'tangent_depth',
                'min_slip_depth')


def _model(window=None, circular=True):
    sd = load_slope_data(str(SAMPLE))
    if window is not None:
        sd['search_window'] = dict(window)
    sd['circular'] = circular
    return sd


class _Recorder:
    """Stands in for a search function: records the kwargs, answers a fixed
    minimum. A sweep only reads FS and the circle out of the result, so a canned
    one is enough -- and it makes every point cost nothing."""

    def __init__(self, circular=True):
        self.calls = []
        self.circular = circular

    def __call__(self, slope_data, method_name, **kw):
        self.calls.append(kw)
        best = {'FS': 1.234, 'Xo': 20.0, 'Yo': 40.0, 'Depth': 5.0}
        if self.circular:
            return [best], True, [], []
        return [best], True, []

    def windows(self):
        """The window constraints seen at each call, one dict per call."""
        return [{k: v for k, v in kw.items() if k in _WINDOW_KEYS}
                for kw in self.calls]


def _sweep(sd, failures, circular=True, **kw):
    """Run a 3-point searched sweep with the search replaced by a recorder."""
    rec = _Recorder(circular=circular)
    name = 'circular_search' if circular else 'noncircular_search'
    orig = getattr(SEARCH, name)
    setattr(SEARCH, name, rec)
    try:
        kw.setdefault('param', 'mat:Soil:c')
        kw.setdefault('values', [10.0, 20.0])
        kw.setdefault('methods', ('bishop',))
        kw.setdefault('num_slices', 20)
        ok, res = S.sensitivity(sd, search=True, **kw)
    finally:
        setattr(SEARCH, name, orig)
    if not ok:
        failures.append(f"the sweep was refused: {res}")
        return rec, None
    return rec, res


# ---------------------------------------------------------------------------

def check_no_window_stays_unconstrained(failures):
    """A model that declares nothing is searched exactly as it was before."""
    rec, res = _sweep(_model(), failures)
    if res is None:
        return
    if not rec.calls:
        failures.append("no window: the sweep never reached the search")
        return
    for i, w in enumerate(rec.windows()):
        if w:
            failures.append(f"no window declared, but point {i} was constrained "
                            f"by {sorted(w)}")


def check_declared_window_reaches_every_point(failures):
    """Every constraint the file declares, at every point of the sweep."""
    rec, res = _sweep(_model(FULL_WINDOW), failures)
    if res is None:
        return
    if len(rec.calls) != 3:
        failures.append(f"a base case plus 2 values should be 3 searches, "
                        f"got {len(rec.calls)}")
    for i, w in enumerate(rec.windows()):
        if w.get('entry_range') != (10.0, 30.0):
            failures.append(f"point {i}: entry_range is {w.get('entry_range')!r}, "
                            f"not the file's (10.0, 30.0)")
        if w.get('exit_range') != (0.0, 8.0):
            failures.append(f"point {i}: exit_range is {w.get('exit_range')!r}, "
                            f"not the file's (0.0, 8.0)")
        if w.get('center_box') != (5.0, 20.0, 40.0, 60.0):
            failures.append(f"point {i}: center_box is {w.get('center_box')!r}, "
                            f"not the file's (5.0, 20.0, 40.0, 60.0)")
        if w.get('min_slip_depth') != 1.5:
            failures.append(f"point {i}: min_slip_depth is "
                            f"{w.get('min_slip_depth')!r}, not the file's 1.5")
        td = w.get('tangent_depth')
        if not (isinstance(td, tuple) and td[0] == -5.0 and td[1] > -5.0):
            failures.append(f"point {i}: tangent_depth is {td!r}; the sheet's "
                            f"single elevation should become a band from it up "
                            f"to the top of the model")


def check_half_declared_range_is_not_a_window(failures):
    """One end of a range is not a window; inventing the other end would
    constrain the search in a direction nobody asked for."""
    rec, res = _sweep(_model({'entry_x_min': 10.0, 'exit_x_max': 8.0}), failures)
    if res is None:
        return
    for i, w in enumerate(rec.windows()):
        if 'entry_range' in w or 'exit_range' in w:
            failures.append(f"point {i}: a half-declared range was forwarded as "
                            f"{sorted(w)}")


def check_caller_beats_the_file(failures):
    """An explicit search_opts wins; the file fills only what it leaves open."""
    rec, res = _sweep(_model(FULL_WINDOW), failures,
                      search_opts={'entry_range': (100.0, 200.0)})
    if res is None:
        return
    for i, w in enumerate(rec.windows()):
        if w.get('entry_range') != (100.0, 200.0):
            failures.append(f"point {i}: the caller's entry_range was overwritten "
                            f"by the file ({w.get('entry_range')!r})")
        if w.get('exit_range') != (0.0, 8.0):
            failures.append(f"point {i}: the file's exit_range should still fill "
                            f"the gap the caller left, got {w.get('exit_range')!r}")


def check_use_file_window_false_turns_it_off(failures):
    """The opt-out searches unconstrained regardless of what the file says."""
    rec, res = _sweep(_model(FULL_WINDOW), failures, use_file_window=False)
    if res is None:
        return
    for i, w in enumerate(rec.windows()):
        if w:
            failures.append(f"point {i}: use_file_window=False still applied "
                            f"{sorted(w)}")


def check_studies_inherit_it(failures):
    """design, back_analysis, tornado and scaled_sensitivity all search inside
    the same window their sweep does."""
    for name, call in (
            ('design', lambda sd: S.design(sd, param='mat:Soil:c', low=10.0,
                                           high=20.0, steps=2, method='bishop',
                                           num_slices=20)),
            ('back_analysis', lambda sd: S.back_analysis(
                sd, param='mat:Soil:c', low=10.0, high=20.0, steps=2,
                method='bishop', num_slices=20)),
            ('tornado', lambda sd: S.tornado(sd, ['mat:Soil:c'], method='bishop',
                                             num_slices=20)),
            ('scaled_sensitivity', lambda sd: S.scaled_sensitivity(
                sd, ['mat:Soil:c'], method='bishop', num_slices=20)),
    ):
        rec = _Recorder()
        orig = SEARCH.circular_search
        SEARCH.circular_search = rec
        try:
            ok, res = call(_model(FULL_WINDOW))
        finally:
            SEARCH.circular_search = orig
        if not ok:
            failures.append(f"{name}: refused the run: {res}")
            continue
        if not rec.calls:
            failures.append(f"{name}: never reached the search")
            continue
        for i, w in enumerate(rec.windows()):
            if w.get('entry_range') != (10.0, 30.0):
                failures.append(f"{name} point {i}: the file's entry_range did "
                                f"not reach the search ({w.get('entry_range')!r})")


def check_noncircular_takes_only_what_it_understands(failures):
    """A non-circular search has no entry/exit/center limits. The one bound it
    does understand is forwarded; the rest are dropped, not raised."""
    sd = _model(FULL_WINDOW, circular=False)
    if not sd.get('non_circ'):
        # the sample carries no non-circular surface; give it a trivial one so
        # the sweep reaches the non-circular branch at all
        sd['non_circ'] = [{'X': 0.0, 'Y': 0.0, 'Movement': 'free'},
                          {'X': 10.0, 'Y': 0.0, 'Movement': 'free'}]
    rec, res = _sweep(sd, failures, circular=False)
    if res is None:
        return
    if not rec.calls:
        failures.append("non-circular: the sweep never reached the search")
        return
    for i, w in enumerate(rec.windows()):
        if w.get('min_slip_depth') != 1.5:
            failures.append(f"non-circular point {i}: min_slip_depth is "
                            f"{w.get('min_slip_depth')!r}, not the file's 1.5")
        extra = sorted(k for k in w if k != 'min_slip_depth')
        if extra:
            failures.append(f"non-circular point {i}: forwarded {extra}, which "
                            f"noncircular_search does not accept")


# ---------------------------------------------------------------------------
# One window, every consumer
# ---------------------------------------------------------------------------

class _Stop(Exception):
    """Raised by the parity recorder once it has the kwargs it came for."""


class _FirstCall(_Recorder):
    """Records the first search call and stops the run there.

    Parity is about what reaches the search, not about what the caller does with
    the answer, so there is no reason to let a reliability campaign run its
    remaining 2N searches (or a Monte Carlo its ten thousand evaluations) to find
    out."""

    def __call__(self, slope_data, method_name, **kw):
        self.calls.append(kw)
        raise _Stop()


def _first_window(fn, circular=True):
    """The window constraints the first search of ``fn()`` is handed."""
    rec = _FirstCall(circular=circular)
    name = 'circular_search' if circular else 'noncircular_search'
    orig = getattr(SEARCH, name)
    setattr(SEARCH, name, rec)
    try:
        fn()
    except _Stop:
        pass
    finally:
        setattr(SEARCH, name, orig)
    if not rec.calls:
        return None
    return rec.windows()[0]


def _studio_window(sd, circular=True):
    """The window Studio's Run LEM path forwards to its search.

    ``analysis_search_kwargs`` is the whole of that path's window handling —
    Studio's Run LEM button and the report's own solve both run through
    :func:`xslope.search.run_lem_analysis`, which reads the window there and
    nowhere else — so asking it with no tolerances set answers exactly what a
    real run would."""
    from xslope.search import analysis_search_kwargs
    kw = analysis_search_kwargs(sd, circular=circular, announce=False)
    return {k: v for k, v in kw.items() if k in _WINDOW_KEYS}


def check_every_consumer_reads_one_window(failures):
    """A sweep, both reliability engines and Studio's Run LEM path all search the
    same model inside the same window.

    The consistency this is here to hold: the window is a property of the MODEL,
    so which entry point happens to be driving the search cannot change which
    mechanism is measured. A divergent private copy of the reading in any one of
    them shows up here as a disagreement.
    """
    from xslope import reliability as R

    sd = _model(FULL_WINDOW)
    # Monte Carlo needs at least one declared sigma before it will search at all;
    # what it searches WITH is what is under test here, not the sampling.
    sd_mc = _model(FULL_WINDOW)
    sd_mc['materials'] = [dict(m) for m in sd_mc['materials']]
    sd_mc['materials'][0]['sigma_c'] = 1.0

    seen = {}

    rec, res = _sweep(_model(FULL_WINDOW), failures)
    seen['a sweep'] = rec.windows()[0] if rec.calls else None
    seen['a Taylor reliability run'] = _first_window(
        lambda: R.reliability_taylor(sd, 'bishop', check_inputs=False))
    seen['a Monte Carlo reliability run'] = _first_window(
        lambda: R.reliability_mc(sd_mc, 'bishop', n_samples=4,
                                 check_inputs=False))
    seen["Studio's Run LEM"] = _studio_window(_model(FULL_WINDOW))

    missing = sorted(k for k, v in seen.items() if v is None)
    if missing:
        failures.append(f"never reached a search: {', '.join(missing)}")
        return
    ref_name, ref = 'a sweep', seen['a sweep']
    if not ref:
        failures.append("a sweep applied no window at all, so there is nothing "
                        "for the other consumers to match")
        return
    for name, w in seen.items():
        if w != ref:
            failures.append(
                f"{name} searched inside {sorted(w.items())}, but {ref_name} "
                f"searched inside {sorted(ref.items())} — the same model must "
                f"reach every search with the same window")


def check_noncircular_parity(failures):
    """The non-circular branch is consistent too: a sweep and Studio both hand a
    non-circular search the one limit it understands, and nothing else.

    This is the branch the two used to disagree on -- Studio applied no window at
    all to a non-circular search while a sweep forwarded ``min_slip_depth`` -- so
    it is pinned on both sides rather than inferred from the circular case.
    """
    sd = _model(FULL_WINDOW, circular=False)
    if not sd.get('non_circ'):
        sd['non_circ'] = [{'X': 0.0, 'Y': 0.0, 'Movement': 'free'},
                          {'X': 10.0, 'Y': 0.0, 'Movement': 'free'}]
    rec, res = _sweep(sd, failures, circular=False)
    if res is None or not rec.calls:
        failures.append("non-circular parity: the sweep never reached the search")
        return
    swept = rec.windows()[0]
    studio = _studio_window(_model(FULL_WINDOW, circular=False), circular=False)
    if swept != studio:
        failures.append(f"non-circular: Studio hands the search "
                        f"{sorted(studio.items())} but a sweep hands it "
                        f"{sorted(swept.items())}")
    if swept != {'min_slip_depth': 1.5}:
        failures.append(f"non-circular: the branch should carry the file's "
                        f"min_slip_depth and nothing else, got "
                        f"{sorted(swept.items())}")


CHECKS = [
    check_no_window_stays_unconstrained,
    check_declared_window_reaches_every_point,
    check_half_declared_range_is_not_a_window,
    check_caller_beats_the_file,
    check_use_file_window_false_turns_it_off,
    check_studies_inherit_it,
    check_noncircular_takes_only_what_it_understands,
    check_every_consumer_reads_one_window,
    check_noncircular_parity,
]


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
    print("sweep_window_check: "
          + (f"{len(bad)} failure(s)" if bad else f"{len(CHECKS)} checks passed"))
    sys.exit(1 if bad else 0)
