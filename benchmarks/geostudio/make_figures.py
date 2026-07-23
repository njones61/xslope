"""Render the side-by-side (inputs | solution) figures for the GeoStudio-only
verification corpus problems (docs/verification/geostudio.md, gs2_*.xlsx).

Reuses the rocscience corpus figure machinery with SRC pointed at
docs/verification/files/geostudio/. Run from the repo root:

    PYTHONPATH=. python3 benchmarks/geostudio/make_figures.py [gs2_NN ...]
"""
import importlib.util
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.join(_HERE, '..', '..')
sys.path.insert(0, _REPO)

_spec = importlib.util.spec_from_file_location(
    'rs_make_figures', os.path.join(_HERE, '..', 'rocscience', 'make_figures.py'))
_mf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mf)
_mf.SRC = os.path.join(_REPO, 'docs', 'verification', 'files', 'geostudio')

# (file stem, surface kind, method) — same vocabulary as the rocscience maker.
CASES = [
    ('gs2_18', 'circle', 'mprice'),    # B&C Case 2: SLOPE/W's critical circle stored
    ('gs2_26', 'noncirc', 'spencer'),  # Baker planar: the critical plane (M-P declines it)
    ('gs2_33', 'noncirc', 'mprice'),   # Priest rigid block: specified plane + tcrack
    ('gs2_45', 'csearch', 'spencer'),  # EC7 cutting (factored parameters)
    ('gs2_46', 'csearch', 'mprice'),   # EC7 dam (factored parameters + FE seepage)
]

if __name__ == '__main__':
    only = set(sys.argv[1:])
    os.makedirs(_mf.OUT, exist_ok=True)
    for stem, kind, method in CASES:
        if only and stem not in only:
            continue
        try:
            out = _mf.make_figure(stem, kind, method)
            print(f'{stem}: {out}')
        except Exception as e:
            print(f'{stem}: FAILED — {type(e).__name__}: {e}')
