"""What happens to a model when its reinforcement is EDITED.

Three ways an edit to the reinforce sheet used to be absorbed instead of applied,
each of them silent:

* A ``Dir`` or ``Appl`` outside the vocabulary. The slicer tested these with
  ``== "tangent"`` and ``== "active"`` and took the other branch for everything
  else, so ``Dir = 0.0`` was applied AXIALLY and a nonsense ``Appl`` as PASSIVE
  capacity -- a different force at a different inclination, reported as if it
  were the model that had been asked for.

* Clearing the lines. ``reinforce_lines`` is DERIVED from
  ``reinforcement_lines``; the slicer fell back to the derived list whenever the
  source was empty, so deleting every reinforcement row left the analysis
  reinforced by the point lists built before the deletion.

* Saving the file. ``Dir`` and ``Appl`` are filled in the sheet by a VLOOKUP on
  ``Type``, and the two columns are one shared-formula group each -- the formula
  itself lives on the master cell at H3/I3 and every row below holds a
  back-reference. Writing a literal into the master deleted the only copy, and
  Dir rows 9-22 read back as a bare ``'='``.

Each leg asserts the fix AND the behavior it must not have disturbed: a blank
still takes the documented default, a model that carries no source key at all
still runs the legacy point-list path, and a row whose Dir/Appl are exactly what
its Type derives is saved as the Type alone with the formula left in place.
"""

import os
import sys
import tempfile
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

REINF_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          'docs', 'lem', 'files', 'xslope_reinforce.xlsx')
TEMPLATE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'docs', 'inputs', 'input_template.xlsx')

#: The Dir column's shared-formula group spans H3:H22 (Appl's spans I3:I22), so
#: rows 9-22 are the ones no reinforcement row of the test file reaches: they are
#: only ever formulas, and are where an orphaned group shows up.
UNTOUCHED_ROWS = (9, 15, 22)


def _model():
    from xslope.fileio import load_slope_data
    return load_slope_data(REINF_FILE)


def _reinforcement_force(slope_data):
    """Total tangent-active reinforcement pull the slices carry, in one number.

    ``check_inputs=False`` is deliberate: the preflight refuses a bad vocabulary
    before the engine sees it, and this check is about what the ENGINE does when
    it is handed one anyway (an importer, a sweep, a model assembled in memory).
    """
    from xslope.slice import generate_slices
    _, (df, _fs) = generate_slices(slope_data, circle=slope_data['circles'][0],
                                   num_slices=40, check_inputs=False, debug=False)
    return float(df['p'].sum())


def _vocabulary_leg():
    """A Dir or Appl that is not a word from the vocabulary must stop the run."""
    failures = []
    for key, value, expect in (('dir', 0.0, 'Dir'), ('dir', 'sideways', 'Dir'),
                               ('appl', 0.0, 'Appl'), ('appl', 'later', 'Appl')):
        sd = _model()
        for r in sd['reinforcement_lines']:
            r[key] = value
        try:
            _reinforcement_force(sd)
        except ValueError as exc:
            text = str(exc)
            if expect not in text or repr(value) not in text:
                failures.append(f"{key}={value!r}: refused, but the message names "
                                f"neither the field nor the value: {text[:90]!r}")
            if "'Line 1'" not in text:
                failures.append(f"{key}={value!r}: refused without naming the row: "
                                f"{text[:90]!r}")
        else:
            failures.append(f"{key}={value!r} was accepted and applied silently")

    # A blank still means "take the documented default" -- the vocabulary guard
    # must not have turned an empty cell into an error.
    sd = _model()
    reference = _reinforcement_force(_model())
    for r in sd['reinforcement_lines']:
        r['dir'] = None
        r['appl'] = ''
    try:
        blank = _reinforcement_force(sd)
    except ValueError as exc:
        failures.append(f"a blank Dir/Appl was refused: {str(exc)[:90]!r}")
    else:
        if abs(blank - reference) > 1e-9:
            failures.append(f"a blank Dir/Appl no longer defaults to tangent/active: "
                            f"{blank} vs {reference}")
    return failures


def _cleared_lines_leg():
    """An explicit empty reinforcement list must mean NO reinforcement."""
    failures = []
    reference = _reinforcement_force(_model())
    if reference <= 0.0:
        return [f"the fixture carries no reinforcement force to remove ({reference})"]

    # The source is emptied and the DERIVED list is deliberately left in place:
    # that is exactly the state an edit produces, and the state the slicer used
    # to fall back into.
    sd = _model()
    sd['reinforcement_lines'] = []
    if not sd.get('reinforce_lines'):
        failures.append("the fixture's derived reinforce_lines were already empty, "
                        "so this leg proves nothing")
    cleared = _reinforcement_force(sd)
    if cleared != 0.0:
        failures.append(f"clearing reinforcement_lines left {cleared} of "
                        f"reinforcement force in the slices (was {reference})")

    # ensure_reinforce_pullout is where a model on its way into an engine is
    # brought up to date; it must clear the derived list too, since the canvas
    # and the report read that one.
    from xslope.fileio import ensure_reinforce_pullout
    sd = _model()
    sd['reinforcement_lines'] = []
    ensure_reinforce_pullout(sd)
    if sd.get('reinforce_lines'):
        failures.append(f"{len(sd['reinforce_lines'])} derived reinforce_lines "
                        f"survived an emptied source list")

    # The legacy point-list path is for models that carry no source key at all,
    # and must still work.
    sd = _model()
    del sd['reinforcement_lines']
    legacy = _reinforcement_force(sd)
    if legacy <= 0.0:
        failures.append(f"a hand-built model carrying only reinforce_lines lost its "
                        f"reinforcement ({legacy})")
    return failures


def _formula_cells(path):
    """{cell: formula or value} for the reinforce sheet's Dir and Appl columns."""
    import openpyxl
    wb = openpyxl.load_workbook(path)
    ws = wb['reinforce']
    out = {}
    for col in ('G', 'H', 'I'):
        for row in range(3, 23):
            out[f'{col}{row}'] = ws[f'{col}{row}'].value
    wb.close()
    return out


def _round_trip(slope_data):
    from xslope.fileio import load_slope_data, save_slope_data_to_xlsx
    tmp = tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False).name
    try:
        save_slope_data_to_xlsx(slope_data, tmp, template=TEMPLATE)
        return load_slope_data(tmp), _formula_cells(tmp)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def _saved_file_leg():
    """A saved file with reinforcement rows must reload, and keep its formulas."""
    failures = []
    fields = ('x1', 'y1', 'x2', 'y2', 't_max', 'lp1', 'lp2', 'type', 'dir', 'appl')

    def _compare(before, after, tag):
        if len(before) != len(after):
            failures.append(f"{tag}: {len(before)} line(s) saved, {len(after)} reloaded")
            return
        for i, (b, a) in enumerate(zip(before, after)):
            for f in fields:
                if isinstance(b[f], float):
                    same = abs(b[f] - a[f]) <= 1e-9
                else:
                    same = b[f] == a[f]
                if not same:
                    failures.append(f"{tag}: line {i + 1} {f} {b[f]!r} -> {a[f]!r}")

    def _formulas_intact(cells, tag, skip_rows=()):
        for cell, value in cells.items():
            row = int(cell[1:])
            if cell[0] == 'G' or row in skip_rows:
                continue
            if not isinstance(value, str) or not value.startswith('='):
                failures.append(f"{tag}: {cell} is no longer a formula ({value!r})")
            elif 'VLOOKUP' not in value:
                failures.append(f"{tag}: {cell} lost its formula body ({value!r}) -- "
                                f"the shared group was orphaned")

    # 1. Dir/Appl exactly as the Type derives them: written as the Type alone,
    #    with every derived cell in both columns left to the formula.
    sd = _model()
    reloaded, cells = _round_trip(sd)
    _compare(sd['reinforcement_lines'], reloaded['reinforcement_lines'], 'preset')
    _formulas_intact(cells, 'preset')

    # 2. Dir/Appl OVERRIDDEN. The two override cells become literals -- that is
    #    what an override is -- and every row the override does not occupy keeps
    #    a working formula rather than a back-reference to a deleted master.
    sd = _model()
    n = len(sd['reinforcement_lines'])
    for r in sd['reinforcement_lines']:
        r['dir'] = 'axial'
        r['appl'] = 'passive'
    reloaded, cells = _round_trip(sd)
    _compare(sd['reinforcement_lines'], reloaded['reinforcement_lines'], 'override')
    _formulas_intact(cells, 'override', skip_rows=tuple(range(3, 3 + n)))
    for row in UNTOUCHED_ROWS:
        for col in ('H', 'I'):
            value = cells[f'{col}{row}']
            if value == '=':
                failures.append(f"override: {col}{row} reads back as a bare '=' -- "
                                f"the shared-formula master was deleted")
    for row in range(3, 3 + n):
        if cells[f'H{row}'] != 'Axial' or cells[f'I{row}'] != 'Passive':
            failures.append(f"override: row {row} did not record the override "
                            f"({cells[f'H{row}']!r}, {cells[f'I{row}']!r})")
    return failures


def run():
    warnings.simplefilter('ignore')
    if not os.path.exists(REINF_FILE):
        return [f"missing fixture {REINF_FILE}"]
    failures = []
    failures += _vocabulary_leg()
    failures += _cleared_lines_leg()
    failures += _saved_file_leg()
    return failures


if __name__ == '__main__':
    problems = run()
    for p in problems:
        print('FAIL:', p)
    print('OK' if not problems else f'{len(problems)} failure(s)')
    sys.exit(1 if problems else 0)
