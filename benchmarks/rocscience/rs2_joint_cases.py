"""The RS2 joint-corpus rows that are registered without a test tag.

A row the page (docs/verification/rs2_joints.md) reports without locking carries
no ``fem_ssrm`` / ``fem_tilt`` tag, so the settings it is run at are stated here
instead. Two readers take them from this one table, so they cannot drift apart:

* ``make_rs2_joint_figures`` solves and draws each row at these settings;
* ``build_joint_problems`` writes each row's K0 into its workbook
  (``benchmarks.tag_k0.declared_k0``), so the file a user opens holds the
  initial stress the reported factor was computed under.

A row that later locks gains a tag, and the tag wins over the entry here in both
readers. Data only: no imports, so the builder can read it without loading the
figure machinery.
"""

#: The settings every reported-only joint row shares.
_JOINT = dict(element_type='tri6', tolerance='0.02', f_min='0.5', f_max='3.0',
              max_iter='250000', tension_srf='false', k0='1')

EXTRA_CASES = [
    # Rows that bracket a factor but cannot lock, because at least one trial of
    # the bracket that defines it -- or of the refinement step that confirms it --
    # reaches the sweep budget without a verdict. They carry no tag, so they are
    # registered here at the settings the page states, and their figures read the
    # mechanism the same way a locked row's does.
    {**_JOINT, 'file': 'files/rocscience/joints/rj003.xlsx',
     'target_size': '12.0', 'benchmark': 'RJ-3'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj005.xlsx',
     'target_size': '12.0', 'benchmark': 'RJ-5'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj006.xlsx',
     'target_size': '12.0', 'benchmark': 'RJ-6'},
    # Problem 20's Voronoi mass. Corpus size is the mean BLOCK width measured on
    # the vendor's own traces, this row's stand-in for a joint spacing. The rest
    # of the termination family — the four problem-1 cases and problems 9 and 11
    # to 14 — have locked and carry tags of their own, so they are not listed
    # here: `registered` would ignore a duplicate entry anyway, and `--audit`
    # names one that is left behind.
    {**_JOINT, 'file': 'files/rocscience/joints/rj020.xlsx',
     'target_size': '2.895', 'benchmark': 'RJ-20'},
]

#: Rows measured by a SWEEP rather than by a bracket. Problem 16 is scored in the
#: tilt angle at which a block grid topples, which XSLOPE reaches by pushing the
#: model with a seismic coefficient at full strength, so there is no bracket for
#: the ordinary producer to draw and no factor of safety to title the panels
#: with. Each entry names the two coefficients the sweep closed on — the last one
#: the stack stands at and the first one it goes at — and the figure draws the
#: second with both angles in its title. See
#: ``make_rs2_joint_figures.make_sweep_figure``.
#:
#: A row the page LOCKS carries a ``type=fem_tilt`` tag with the same keys, and
#: the tag wins over the entry here (see ``make_rs2_joint_figures.sweep_cases``),
#: exactly as a ``fem_ssrm`` tag wins over an ``EXTRA_CASES`` entry.
SWEEP_CASES = [
    {'file': 'files/rocscience/joints/rj016.xlsx', 'benchmark': 'RJ-16',
     'target_size': '0.09', 'max_iter': '250000', 'element_type': 'tri6',
     'tension_srf': 'false', 'k0': '1',
     'k_stand': '0.179846', 'k_fail': '0.181396'},
]
