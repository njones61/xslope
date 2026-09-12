"""The Studio's "Build network…" dialog, and what the button does with it.

The dialog is where a joint network is described — a kind, its parameters, the
region it exists in, one property set — and the button beside the joints editor's
list is what opens it. Between them they have to hold three promises:

  1. **What the preview shows is what OK writes.** The dialog generates the set
     as the fields change, counts it, and hands back exactly those rows.
  2. **A set reopens as itself.** Selecting one of a set's rows and pressing the
     button fills the dialog in from the row's own label, identically enough that
     pressing OK without touching anything reproduces the set.
  3. **Editing a set replaces it in place.** The regenerated rows land where the
     old ones were; the other set, a hand-entered joint line and the properties
     the rows carry are untouched.

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
    if set_name(dlg._name.text() + '-01|par') == '':
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
    if 'poly:North block' not in narrowed[0]['label']:
        failures.append(f"dialog: the region is not in the rows' label "
                        f"({narrowed[0]['label']!r})")
    dlg._band_lo.setText('8')
    banded = dlg.result_rows()
    if not banded or len(banded) >= len(narrowed):
        failures.append(f"dialog: the elevation band took the set from "
                        f"{len(narrowed)} rows to {len(banded)}")
    if 'band=8:' not in banded[0]['label']:
        failures.append(f"dialog: the band is not in the rows' label "
                        f"({banded[0]['label']!r})")

    # A name another set is already using is refused, by name, before OK.
    dlg._name.setText('st')
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
                   f"the band above 8 {len(banded)}, each with the region and the "
                   f"band in its label; a name in use and a pair at one dip are "
                   f"refused in the dialog rather than on OK")


def _leg_reopen(failures, results):
    """Leg 2: a set reopens as itself."""
    from studio.network_dialog import BuildNetworkDialog

    sd = _model()
    original = JointSet('bed', 'parallel',
                        {'dip': 25.0, 'spacing': 2.5, 'offset': 0.5},
                        region='poly:North block', band=(6.0, None),
                        props={'c': 0.0, 'phi': 34.0, 'dil': 5.0, 'jred': 'No'})
    rows = original.generate(sd)
    reopened = BuildNetworkDialog(sd, sd['joint_lines'] + rows,
                                  JointSet.from_rows(rows))
    back = reopened.result_set()
    if back is None:
        failures.append("reopen: the dialog does not describe a set at all")
        return
    if back != original:
        failures.append(f"reopen: the set came back as {back.to_label(1)!r}, not "
                        f"{original.to_label(1)!r}")
    if reopened.editing() != 'bed':
        failures.append(f"reopen: the dialog is editing {reopened.editing()!r}")
    if [r['label'] for r in reopened.result_rows()] != [r['label'] for r in rows]:
        failures.append("reopen: pressing OK unchanged would not reproduce the set")
    if back.props.get('jred') != 'No' or float(back.props.get('phi')) != 34.0:
        failures.append(f"reopen: the properties came back as {back.props}")
    # Its own name is not a clash with itself.
    if not reopened._buttons.button(reopened._buttons.StandardButton.Ok).isEnabled():
        failures.append("reopen: the set's own name was read as a clash")
    results.append("reopen      a set with a region, a band, an offset and a "
                   "property set reopens as the same record, down to the row "
                   "labels OK would write")


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
        # A new set: appended, nothing else touched.
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

        # Reopened on that set, at half the spacing: replaced in place.
        editor.replace_rows(after, select=len(before))       # a row of 'bed'
        _Stub.fields = {'dip': '0', 'spacing': '2', '_name': 'bed', '_phi': '32'}
        _build_joint_network(sd, editor)
        edited = editor.result_rows()
        n_bed2 = sum(1 for r in edited if set_name(r.get('label')) == 'bed')
        if n_bed2 <= n_bed:
            failures.append(f"button: halving the spacing took the set from "
                            f"{n_bed} rows to {n_bed2}")
        if sum(1 for r in edited if not set_name(r.get('label'))) != n_typed:
            failures.append("button: the hand-entered joint line did not survive "
                            "a regeneration")
        if sum(1 for r in edited if set_name(r.get('label')) == 'st') != n_steep:
            failures.append("button: the OTHER set did not survive a regeneration")
        if len(edited) != len(before) + n_bed2:
            failures.append(f"button: regenerating left {len(edited)} rows, not "
                            f"{len(before) + n_bed2} — the old set was not replaced")
        if [r['label'] for r in edited[:len(before)]] != [r['label'] for r in before]:
            failures.append("button: regenerating moved the rows around it")
        if any(float(r['phi']) != 32.0
               for r in edited if set_name(r.get('label')) == 'bed'):
            failures.append("button: the regenerated rows lost their properties")
    finally:
        network_dialog.BuildNetworkDialog = _Stub.__bases__[0]
    results.append(f"button      a new set appends {n_bed} rows after the "
                   f"{len(before)} already there; reopening it at half the "
                   f"spacing replaces exactly those rows in place with "
                   f"{n_bed2}, leaving the hand-entered line, the other set and "
                   f"the properties alone")


def run():
    """Returns a list of failure strings (empty = pass)."""
    import time
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    failures, results = [], []
    t0 = time.time()
    _leg_dialog(failures, results)
    _leg_reopen(failures, results)
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
          "reopens one as itself, and the button replaces an edited set in "
          "place.")


if __name__ == '__main__':
    main()
