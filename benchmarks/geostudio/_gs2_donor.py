"""Donor material handling shared by the GeoStudio-only corpus builders.

Every ``build_gs2_*.py`` starts from ``docs/lem/files/xslope_acads_simple.xlsx``
purely because it is a valid metric single-material input file: the builder then
overwrites the geometry, the material table, the water condition and the trial
surfaces, keeping nothing of the ACADS problem.  Nothing, that is, except what a
``dict(sd['materials'][0])`` copies silently.

That donor is also an RS2 SSRM benchmark (RS2-1), so its one material carries the
elastic constants and tensile cap transcribed from the vendor .fez for the ACADS
soil -- E = 50,000 kPa, nu = 0.4, t_cut = 3 kPa.  Copying the material dict
wholesale pulled those three into every GS-2 problem, where they describe nothing:
they were the ACADS soil's numbers, not a Cannon Dam shell's or a heap-leach ore
pile's.  They are also invisible to every GS-2 answer -- the rows on
docs/verification/geostudio.md are limit-equilibrium factor-of-safety rows and
transient seepage head rows, and neither solver reads E, nu or t_cut -- so the
contamination was silent, except that an inherited 3 kPa cap sits below most of
these materials' Mohr-Coulomb apex and made ``load_slope_data`` print a spurious
"the cutoff never binds" warning.

:func:`donor_material` therefore hands back the donor material with the FEM-only
fields reset to what the GS-2 corpus carries: the template default E/nu and no
tensile cutoff.  Should a GS-2 row ever be solved with the FEM, that row's file
needs real constants -- the vendor's if the vendor specified them, otherwise
``benchmarks/rocscience/elastic_props.assign_elastic_props`` -- and must stop
using this default.
"""

# Template-default elastic constants, which is what every checked-in gs2_*.xlsx
# carries.  Inert on every GS-2 row (no FEM row exists); see the module docstring
# before reusing them for one that is.
E_DEFAULT = 100000.0
NU_DEFAULT = 0.3

# The donor's OTHER inheritance: it is also the design/reliability sample, so its one
# material carries the ACADS soil's standard deviations (s(g) = 1.2, s(c) = 1.8,
# s(f) = 2.744).  Those describe no GS-2 problem, and every GS-2 row is a
# deterministic factor of safety or a head, so a copied sigma changes no answer -- but
# it does put a bogus one-click +/-sigma range in front of anyone who opens the file,
# and would silently corrupt a probabilistic run on it.  Zeroed with the rest.
_SIGMA_KEYS = ('sigma_gamma', 'sigma_c', 'sigma_phi',
               'sigma_cp', 'sigma_d', 'sigma_psi')


def donor_material(slope_data, **overrides):
    """A copy of the donor's material 0 with the donor's FEM-only fields and its
    reliability standard deviations dropped.

    ``**overrides`` are applied on top, so a builder can write
    ``donor_material(sd, name='Clay', c=5.0, phi=30.0)`` in one call.
    """
    m = dict(slope_data['materials'][0])
    m.update(E=E_DEFAULT, nu=NU_DEFAULT, t_cut=None)
    m.update({k: 0.0 for k in _SIGMA_KEYS})
    m.update(overrides)
    return m
