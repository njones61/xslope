"""PARKED / BLOCKED builder for the SEEP/W transient verification example
"Infiltration into multi-layered system" (SEEPW-T06).

    *** This model is NOT built into the corpus. It is blocked by the transient
    solver / lock-runner envelope and writes NO .xlsx by default. ***

The example ponds water on top of a fourteen-layer laboratory soil profile and
watches the wetting front descend (the "Infiltration Experiment" analysis), then
lets it drain (the "Drainage Experiment" child).  The census flagged the DRAINAGE
leg as out of envelope up front -- it is hysteretic, and XSLOPE carries a single
retention curve per material (no wetting/drying hysteresis).  The plan was to port
the INFILTRATION leg only.  At build time two further, harder walls appeared in
the infiltration leg itself, and neither is a builder change:

  1. NON-STEADY, PER-LAYER INITIAL CONDITION.  The vendor sets the initial state
     through a per-material ``InitPWP`` attribute, layer by layer, and it is not a
     hydrostatic or steady field -- it is an imposed measurement profile:

         top  y 1.02-1.10  Layer 1   InitPWP  -8 kPa
              y 0.96-1.02  Layer 2            -8
              y 0.92-0.96  Layer 3            -5
              y 0.86-0.92  Layer 4           -50   <- a suction SPIKE
              y 0.80-0.86  Layer 5           -30
              y 0.70-0.80  Layer 6           -30
              y 0.64-0.70  Layer 7           -30
              y 0.50-0.64  Layer 8           -20
              ...          Layers 9-14       -20
         base y 0.00-0.10  Layer 14          -20

     There is no steady solve that returns this field (a -50 kPa spike wedged
     between -5 and -30 kPa layers is nobody's equilibrium).  XSLOPE's transient
     solver CAN take an arbitrary field through its ``h_init`` argument -- but the
     corpus lock runner (``run_tseep_head_test`` in run_tests.py) computes the
     initial condition as a t = 0 STEADY solve of the boundary configuration; it
     has no way to inject a per-node initial head from a tag.  So even though the
     forward solve is expressible in code, it cannot be locked as a corpus tag.

  2. UNIT-GRADIENT (FREE-DRAINAGE) BASE BOUNDARY.  The base is a "Unit Gradient"
     boundary (dH/dn set so the gradient is gravity-only, q = kr*Ksat out).
     XSLOPE's seepage BC set is specified head, specified flux and potential
     seepage (exit) face -- there is no unit-gradient boundary.  An exit face
     clamps the base to psi = 0 rather than letting it drain under a unit
     gradient, so it does not reproduce this outlet.

Either wall alone would block a faithful lock; together they put the infiltration
leg out of the current envelope.  Unblocking wall (1) is the same non-submerged /
per-node initial-state capability the lock machinery would need for any imposed-IC
transient (a solver+run_tests change, gated on Norm, cf. SEEPW-T07); wall (2) is a
new boundary-condition type in seep.py.  Until both land, T06 stays PLANNED.  The
van Genuchten storage + kr path this example exercises is already covered against
clean oracles by SEEPW-T02 (infiltration into dry soil) and SEEPW-T05 (leach
column), so nothing is lost by parking it.

Model (from the vendor .gsz, read-only oracle -- never committed):
    1-D column 1.1 m tall, 14 SatUnsat layers (each a 66-point tabulated VWC +
    conductivity function; a per-layer van Genuchten fit would be needed for a
    port).  Top (y = 1.1): ponded, specified pressure head +0.1 m.  Base (y = 0):
    unit-gradient free drainage.  Infiltration duration 3600 s, 120 steps; the
    published reference is a field / HYDRUS-1D VWC profile (Zettl 2011, Huang
    2011), with SEEP/W's own node.csv as the per-step oracle.

Run:  PYTHONPATH=. python3 benchmarks/geostudio/build_gs2_mlayer.py
"""

import sys

# NOT registered as a corpus builder: T06 is blocked by the two envelope walls in
# the module docstring (non-steady per-layer InitPWP that the steady-IC lock
# runner cannot reproduce, and an unsupported unit-gradient base boundary). No
# corpus .xlsx is emitted.
BUILDERS = []  # intentionally empty

if __name__ == '__main__':
    print('SEEPW-T06 is PARKED/BLOCKED. The infiltration leg needs (1) a per-node '
          'imposed initial condition the corpus lock runner cannot express (its IC '
          'is a steady solve) and (2) a unit-gradient free-drainage base boundary '
          'the solver does not have. No corpus file written. See the module '
          'docstring; the vG storage+kr path is already covered by SEEPW-T02/T05.')
    sys.exit(0)
