"""Figure generator for the VP102 / RS2 Part IV #102 transient rapid-drawdown port.

Plots the factor-of-safety-vs-time drawdown curve — XSLOPE's own transient result
against the published Slide2 Spencer (Table 102.3) and RS2 SSR (Tables 102.3/102.4)
columns — for the section in rocscience.md#vp102 (and cross-referenced from
rs2.md#p4-vp102). The XSLOPE points are the locked regression values; the published
columns are transcribed from the Huang & Jia (2008) / RS2 Part IV manual tables.

Run:  python benchmarks/rocscience/make_vp102t_figure.py
Out:  docs/verification/images/vp102t_curve.png
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IMG = os.path.join(os.path.dirname(__file__), '..', '..',
                   'docs', 'verification', 'images')

# Published Slide2 Spencer, Table 102.3 (phi_b = 0), incl. the initial steady state.
_T_PUB = [0, 60, 70, 75, 80, 85, 90, 100, 300, 600, 1000, 1500]
_SLIDE_SPENCER = [1.745, 1.804, 1.820, 1.828, 1.835, 1.842, 1.851, 1.867,
                  2.092, 2.242, 2.329, 2.373]
_RS2_SSR_C2 = [1.70, 1.77, 1.78, 1.79, 1.80, 1.81, 1.82, 1.83, 2.06, 2.19, 2.26, 2.29]

# XSLOPE locked regression values (initial = vp102b Spencer; transient = vp102t_*).
_XS_T = [0, 60, 100, 300, 600, 1500]
_XS_SPENCER = [1.719, 1.724, 1.765, 1.967, 2.140, 2.299]
# XSLOPE SSRM (Case 2, phi_b = 0): dry end-member context is 2.370; transient below.
_XS_SSR_T = [60, 300, 1500]
_XS_SSR_C2 = [1.670, 1.888, 2.195]


def vp102t_curve():
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(_T_PUB, _SLIDE_SPENCER, '-', color='#1f6feb', lw=1.6,
            marker='o', ms=4, label='Slide2 Spencer (Table 102.3)')
    ax.plot(_T_PUB, _RS2_SSR_C2, '--', color='#8250df', lw=1.4,
            marker='^', ms=4, label='RS2 SSR, $\\varphi^b=0$ (Table 102.3)')
    ax.plot(_XS_T, _XS_SPENCER, '-', color='#cf222e', lw=1.8,
            marker='s', ms=6, label='XSLOPE Spencer (LEM)')
    ax.plot(_XS_SSR_T, _XS_SSR_C2, ':', color='#bc4c00', lw=1.6,
            marker='D', ms=6, label='XSLOPE SSRM, $\\varphi^b=0$')

    ax.set_xlabel('Time after drawdown (hours)')
    ax.set_ylabel('Factor of safety')
    ax.set_title('VP102 rapid drawdown: factor of safety vs time\n'
                 'homogeneous earth dam, reservoir el. 24 $\\to$ 7 at $t=0$ '
                 '(Huang & Jia 2008)', fontsize=10)
    ax.set_xlim(-30, 1560)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(loc='lower right', fontsize=8, framealpha=0.95)
    ax.annotate('initial steady state\n(governing minimum)',
                xy=(0, 1.72), xytext=(120, 1.62), fontsize=7.5,
                arrowprops=dict(arrowstyle='->', color='0.4', lw=0.8), color='0.3')
    fig.tight_layout()
    out = os.path.join(IMG, 'vp102t_curve.png')
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


if __name__ == '__main__':
    print(vp102t_curve())
