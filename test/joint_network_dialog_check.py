"""The Studio's "Build network…" dialog, and what the button does with it.

The dialog is where a joint network is described — a kind, its parameters, the
region it exists in, one property set — and the button beside the joints editor's
list is what opens it. Between them they have to hold three promises:

  1. **What the preview shows is what OK writes.** The dialog generates the set
     as the fields change, counts it, and hands back exactly those rows — named
     for the set and nothing more.
  2. **The next build starts where the last one left off.** A set is not stored,
     so nothing reopens one; within a session the dialog offers the parameters it
     was last used with, under a name the model is not already using.
  3. **The button appends.** The rows land after the ones already there, and
     nothing that was there — the other set, a hand-entered joint line — moves or
     changes.

Offscreen Qt, no file and no solver work.

Run directly:  PYTHONPATH=. python3 test/joint_network_dialog_check.py
"""

import os
import sys
import warnings

warnings.filterwarnings('ignore')
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from shapely.geometry import Polygon                          # noqa: E402

from xslope.joints import JointSet, set_name                  # noqa: E402


def _model():
    """A 40 x 20 section of one material, with a joint region over its left half
    and two joint lines already on the sheet: one hand-entered, one a set."""
    domain = Polygon([(0.0, 0.0), (40.0, 0.0), (40.0, 20.0), (0.0, 20.0)])
    sd = {
        'domain_polygon': domain,
        'polygons': [{'polygon': domain, 'mat_id': 0}],
        'materials': [{'name': 'Sandstone'}],
        'joint_zones': [{'polygon': [(2.0, 2.0), (18.0, 2.0), (18.0, 18.0),
                                     (2.0, 18.0)],
                         'label': 'North block', 'size': None, 'mat_id': None}],
        'joint_lines': [],
    }
    typed = {'label': 'base joint', 'x1': 0.0, 'y1': 1.0, 'x2': 40.0, 'y2': 1.0,
             'c': 0.0, 'phi': 20.0}
    steep = JointSet('st', 'parallel', {'dip': 80.0, 'spacing': 6.0},
                     props={'phi': 28.0})
    sd['joint_lines'] = [dict(typed)] + steep.generate(sd)
    return sd


def _leg_dialog(failures, results):
    """Leg 1: the dialog describes a set, counts it, and hands back its rows."""
    from studio.network_dialog import BuildNetworkDialog

    sd = _model()
    dlg = BuildNetworkDialog(sd, sd['joint_lines'], None)
    if dlg.result_rows():
        failures.append("dialog: an empty form generated rows")
    if dlg._buttons.button(dlg._buttons.StandardButton.Ok).isEnabled():
        failures.append("dialog: OK is offered on a form with no parameters")
    if set_name(dlg._name.text() + '-01') == '':
        failures.append(f"dialog: the name it opens on, {dlg._name.text()!r}, is "
                        f"not a usable set name")
    if dlg._name.text() in {set_name(r['label']) for r in sd['joint_lines']}:
        failures.append(f"dialog: it opens on {dlg._name.text()!r}, which the "
                        f"model is already using")

    edits = dlg._param_edits['parallel']
    edits['dip'].setText('30')
    edits['spacing'].setText('2')
    dlg._prop_edits['phi'].setText('32')
    dlg._prop_edits['dil'].setText('5')
    rows = dlg.result_rows()
    want = JointSet('set1', 'parallel', {'dip': 30.0, 'spacing': 2.0},
                    props={'phi': 32.0, 'dil': 5.0}).generate(sd)
    if len(rows) != len(want) or not rows:
        failures.append(f"dialog: {len(rows)} rows for a set the generator makes "
                        f"{len(want)} of")
    elif [r['label'] for r in rows] != [r['label'] for r in want]:
        failures.append("dialog: the rows it hands back are not the set's own")
    if not str(len(rows)) in dlg._status.text():
        failures.append(f"dialog: the count is not reported ({dlg._status.text()!r})")
    if any(float(r['phi']) != 32.0 or float(r['dil']) != 5.0 for r in rows):
        failures.append("dialog: the property set did not reach every row")

    # The region and the band narrow it, and both reach the label.
    dlg._region.setCurrentIndex(dlg._region.count() - 1)      # the joint region
    narrowed = dlg.result_rows()
    if not narrowed or len(narrowed) >= len(rows):
        failures.append(f"dialog: confining the set to the joint region took it "
                        f"from {len(rows)} rows to {len(narrowed)}")
    if set_name(narrowed[0]['label']) != dlg._name.text():
        failures.append(f"dialog: the rows are not named for the set "
                        f"({narrowed[0]['label']!r})")
    dlg._band_lo.setText('8')
    banded = dlg.result_rows()
    if not banded or len(banded) >= len(narrowed):
        failures.append(f"dialog: the elevation band took the set from "
                        f"{len(narrowed)} rows to {len(banded)}")
    if any(ch in banded[0]['label'] for ch in '|='):
        failures.append(f"dialog: a row carries more than the set's name "
                        f"({banded[0]['label']!r})")

    # A name another set is already using is refused, by name, before OK.
    # (Regeneration is debounced; reading it back flushes rather than waits.)
    dlg._name.setText('st')
    dlg._flush()
    if dlg._buttons.button(dlg._buttons.StandardButton.Ok).isEnabled():
        failures.append("dialog: a set name the model already uses was accepted")
    if 'st' not in dlg._status.text():
        failures.append(f"dialog: the clash does not name the set "
                        f"({dlg._status.text()!r})")
    dlg._name.setText('bed')

    # A refusal from the generator is shown, not raised.
    dlg._kind.setCurrentIndex(1)
    for key, value in (('dip', '45'), ('dip2', '45'), ('spacing', '2'),
                       ('spacing2', '2')):
        dlg._param_edits['cross'][key].setText(value)
    if dlg.result_rows():
        failures.append("dialog: two sets at one dip generated a network")
    if 'lie on one another' not in dlg._status.text():
        failures.append(f"dialog: the generator's refusal is not shown "
                        f"({dlg._status.text()!r})")
    results.append(f"dialog      an empty form offers no OK; a dip and a spacing "
                   f"give {len(rows)} rows, the joint region {len(narrowed)} and "
                   f"the band above 8 {len(banded)}, every row named for the set "
                   f"and nothing more; a name in use and a pair at one dip are "
                   f"refused in the dialog rather than on OK")


def _leg_session_memory(failures, results):
    """Leg 2: the next build starts where the last one left off.

    Nothing reopens a set — the model does not keep one — so what the dialog owes
    somebody building a second network is the numbers they just typed, under a
    name that cannot be written over the first set.
    """
    from studio import network_dialog
    from studio.network_dialog import BuildNetworkDialog

    sd = _model()
    first = BuildNetworkDialog(sd, sd['joint_lines'])
    for key, value in (('dip', '25'), ('spacing', '2.5'), ('offset', '0.5')):
        first._param_edits['parallel'][key].setText(value)
    first._region.setCurrentIndex(first._region.count() - 1)   # the joint region
    first._band_lo.setText('6')
    for key, value in (('c', '0'), ('phi', '34'), ('dil', '5')):
        first._prop_edits[key].setText(value)
    first._prop_edits['jred'].setCurrentIndex(
        first._prop_edits['jred'].findText('No'))
    described = first.result_set()
    rows = first.result_rows()
    if described is None or not rows:
        failures.append("session: the first dialog describes no set")
        return
    first.accept()                       # what OK does, without a window

    sd['joint_lines'] = sd['joint_lines'] + rows
    second = BuildNetworkDialog(sd, sd['joint_lines'])
    back = second.result_set()
    if back is None:
        failures.append("session: the second dialog opened on nothing")
        return
    if (back.kind != described.kind or back.params != described.params
            or back.region != described.region or back.band != described.band):
        failures.append(f"session: it opened on {back.params} / {back.region} / "
                        f"{back.band}, not on what was last built")
    if back.props != described.props:
        failures.append(f"session: the properties came back as {back.props}")
    if back.name == described.name:
        failures.append(f"session: it opened on the name {back.name!r}, which "
                        f"the model is already using")
    if back.name in {set_name(r['label']) for r in sd['joint_lines']}:
        failures.append(f"session: {back.name!r} is a name already in the model")
    if not second._buttons.button(second._buttons.StandardButton.Ok).isEnabled():
        failures.append("session: the second set cannot be written")
    if network_dialog._LAST_SET is None:
        failures.append("session: nothing was remembered")
    network_dialog._LAST_SET = None      # a check leaves no state behind
    results.append(f"session     a second dialog opens on the first's dip, "
                   f"spacing, offset, region, band and properties, under "
                   f"{back.name!r} rather than {described.name!r}")


def _leg_button(failures, results):
    """Leg 3: what the editor's button does with what the dialog returns."""
    from PySide6.QtWidgets import QDialog

    from studio import network_dialog
    from studio.editors import JointsEditor, _build_joint_network

    sd = _model()
    editor = JointsEditor().build(sd, None)
    before = editor.result_rows()
    n_typed = sum(1 for r in before if not set_name(r.get('label')))
    n_steep = sum(1 for r in before if set_name(r.get('label')) == 'st')

    class _Stub(network_dialog.BuildNetworkDialog):
        """The dialog, driven instead of shown: the fields the user would fill
        in are set before exec(), which accepts without a window."""
        fields = {}

        def exec(self):
            for key, value in self.fields.items():
                if not key.startswith('_'):
                    self._param_edits[self._kind_at()][key].setText(value)
            self._name.setText(self.fields.get('_name', self._name.text()))
            self._prop_edits['phi'].setText(self.fields.get('_phi', '30'))
            self._refresh()
            return QDialog.Accepted

    network_dialog.BuildNetworkDialog = _Stub
    try:
        # A set: appended, nothing that was there touched.
        _Stub.fields = {'dip': '0', 'spacing': '4', '_name': 'bed', '_phi': '32'}
        _build_joint_network(sd, editor)
        after = editor.result_rows()
        n_bed = sum(1 for r in after if set_name(r.get('label')) == 'bed')
        if n_bed == 0:
            failures.append("button: the new set was not written")
        if len(after) != len(before) + n_bed:
            failures.append(f"button: {len(before)} rows + {n_bed} generated came "
                            f"to {len(after)}")
        if [r['label'] for r in after[:len(before)]] != [r['label'] for r in before]:
            failures.append("button: a new set did not land after the rows that "
                            "were there")
        if any(float(r['phi']) != 32.0
               for r in after if set_name(r.get('label')) == 'bed'):
            failures.append("button: the set's properties did not reach its rows")

        # A SECOND set beside it: appended in turn, with the first still whole.
        _Stub.fields = {'dip': '70', 'spacing': '5', '_name': 'cross', '_phi': '30'}
        _build_joint_network(sd, editor)
        both = editor.result_rows()
        n_cross = sum(1 for r in both if set_name(r.get('label')) == 'cross')
        if n_cross == 0:
            failures.append("button: the second set was not written")
        if sum(1 for r in both if set_name(r.get('label')) == 'bed') != n_bed:
            failures.append("button: building a second set changed the first")
        if sum(1 for r in both if not set_name(r.get('label'))) != n_typed:
            failures.append("button: the hand-entered joint line did not survive")
        if sum(1 for r in both if set_name(r.get('label')) == 'st') != n_steep:
            failures.append("button: the set already on the sheet did not survive")
        if len(both) != len(after) + n_cross:
            failures.append(f"button: {len(after)} rows + {n_cross} generated came "
                            f"to {len(both)}")
        if [r['label'] for r in both[:len(after)]] != [r['label'] for r in after]:
            failures.append("button: the second set moved the rows before it")
    finally:
        network_dialog.BuildNetworkDialog = _Stub.__bases__[0]
        network_dialog._LAST_SET = None
    results.append(f"button      a set appends {n_bed} rows after the "
                   f"{len(before)} already there, and a second appends "
                   f"{n_cross} more, leaving the first set, the other set and "
                   f"the hand-entered line exactly where they were")


def run():
    """Returns a list of failure strings (empty = pass)."""
    import time
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    failures, results = [], []
    t0 = time.time()
    _leg_dialog(failures, results)
    _leg_session_memory(failures, results)
    _leg_button(failures, results)
    app.processEvents()
    print(f"Joint network dialog check ({time.time() - t0:.0f} s):")
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
    print("\nThe dialog describes a set, previews exactly what it writes, "
          "opens the next one on the last one's parameters, and the button "
          "appends without disturbing anything already on the sheet.")


if __name__ == '__main__':
    main()
