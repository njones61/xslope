# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np


#: What every place that shows a rapid drawdown result says when the surface
#: crosses no material with d and ψ — the Log, the Analysis Report and Studio's
#: LEM Solution tab — so the three say it in the same words. ``surface`` is the
#: noun for the surface analyzed: "critical circle", "circle", "surface", ...
NO_LOW_K_SENTENCE = ("The {surface} crosses no material with d and ψ, so stages 1 "
                     "and 2 used drained strengths (c′, φ′) and stage 3 was not "
                     "required.")


def crosses_no_low_k(results):
    """True when a rapid drawdown result's surface crossed no material with d and
    ψ, so every stage was a drained analysis.

    Read from the per-slice record :func:`rapid_drawdown` leaves on its result
    (``rapid_low_k``); False for any other result, and for a rapid drawdown
    result made before that record existed.
    """
    if not isinstance(results, dict) or "stage1_FS" not in results:
        return False
    flags = results.get("rapid_low_k")
    return flags is not None and not any(flags)


def no_low_k_sentence(results, surface):
    """:data:`NO_LOW_K_SENTENCE` for ``surface`` where it applies, else ``""``."""
    return NO_LOW_K_SENTENCE.format(surface=surface) if crosses_no_low_k(results) else ""


def validate_rapid_drawdown(slope_data):
    """
    Validates that slope_data has the required inputs for rapid drawdown analysis.
    Prints warnings for missing optional inputs and raises ValueError for missing required inputs.
    """
    materials = slope_data['materials']

    # Hard error: at least one material must have non-zero d or psi
    has_d_psi = any(m.get('d', 0) != 0 or m.get('psi', 0) != 0 for m in materials)
    if not has_d_psi:
        raise ValueError("Rapid drawdown requires at least one material with non-zero d or psi values. Check your input template.")

    # Warning: no second set of distributed loads
    from .water import water_loads_mode
    if water_loads_mode(slope_data) != "auto" and not slope_data.get('dloads2'):
        print("[WARNING] Rapid drawdown: no second set of distributed loads (dloads2) found.")

    # Warning: piezo method selected but no second piezo line
    has_piezo = any(m.get('u') == 'piezo' for m in materials)
    if has_piezo and not slope_data.get('piezo_line2'):
        print("[WARNING] Rapid drawdown: piezo method selected but no second piezo line found.")

    # Warning: seep method selected but no second seep solution
    has_seep = any(m.get('u') == 'seep' for m in materials)
    if has_seep and 'seep_u2' not in slope_data:
        print("[WARNING] Rapid drawdown: seep method selected but no second seep solution found.")


#: The two ways stage 2 of the two-stage seismic analysis can take its undrained
#: strength, and the words the Log and the report use for each.
SEISMIC2_STRENGTHS = {
    "interpolated": "interpolated between the K_c = 1 and drained envelopes",
    "kc1": "from the K_c = 1 envelope",
}


def seismic2_strength(seismic2):
    """The stage 2 strength choice a ``seismic2`` argument names, or None when it
    asks for no two-stage seismic analysis.

    ``seismic2`` is False (no two-stage analysis), True or ``"interpolated"`` (the
    default: Duncan, Wright & Wong's interpolation, as in rapid drawdown), or
    ``"kc1"`` (the K_c = 1 envelope at the stage 1 stresses, with no
    interpolation, the way Duncan, Wright & Brandon (2014) computed their published
    two-stage factors of safety (four soils on one slope at k = 0.15)). It may
    also be a dict, ``{"strength": <one of those>, "factor": <0 < f <= 1>}``, which
    adds the stage 2 strength factor (:func:`seismic2_factor`); a dict without
    ``"strength"`` means ``"interpolated"``. Any other value raises ValueError."""
    if isinstance(seismic2, dict):
        seismic2_factor(seismic2)
        unknown = set(seismic2) - {"strength", "factor"}
        if unknown:
            raise ValueError(f"seismic2 takes the keys 'strength' and 'factor', not "
                             f"{sorted(unknown)}.")
        strength = seismic2.get("strength", "interpolated")
        if strength is False or strength is None:
            raise ValueError("A seismic2 dict names a stage 2 strength: 'interpolated' "
                             "or 'kc1'.")
        return seismic2_strength(strength)
    if seismic2 is False or seismic2 is None:
        return None
    if seismic2 is True:
        return "interpolated"
    if isinstance(seismic2, str) and seismic2.lower() in SEISMIC2_STRENGTHS:
        return seismic2.lower()
    raise ValueError(f"seismic2 must be False, True, 'interpolated', 'kc1', or a dict "
                     f"with 'strength' and 'factor', not {seismic2!r}.")


def seismic2_factor(seismic2):
    """The stage 2 strength factor a ``seismic2`` argument carries: the number every
    stage 2 undrained strength (d/psi and c/p materials) is multiplied by, 1 unless
    ``seismic2`` is a dict with a ``"factor"``. The screening methods' strength
    reduction (Duncan, Wright & Brandon 2014, Table 10.1; 0.8 to 0.85 of the static
    undrained strength for soils that don't liquefy) goes here: in a two-stage run
    it can't be applied to the inputs, which also set stage 1. A factor outside
    0 < f <= 1 raises ValueError."""
    if not isinstance(seismic2, dict) or seismic2.get("factor") is None:
        return 1.0
    try:
        f = float(seismic2["factor"])
    except (TypeError, ValueError):
        raise ValueError(f"The stage 2 strength factor must be a number, not "
                         f"{seismic2['factor']!r}.") from None
    if not (0.0 < f <= 1.0):
        raise ValueError(f"The stage 2 strength factor must be greater than 0 and at "
                         f"most 1, not {f:g}.")
    return f


def seismic2_tag(seismic2):
    """The label a run's announcement puts after the method's name: ``""`` with no
    two-stage analysis, else " (two-stage seismic)" with the K_c = 1 choice and
    the stage 2 strength factor named when they apply."""
    strength = seismic2_strength(seismic2)
    if not strength:
        return ""
    words = ", K_c = 1 envelope" if strength == "kc1" else ""
    return f" (two-stage seismic{words}{seismic2_factor_words(seismic2_factor(seismic2))})"


def seismic2_factor_words(factor):
    """How the Log, the report, and run labels name a stage 2 strength factor:
    ``""`` at 1, ``", times 0.80"`` otherwise."""
    return "" if factor == 1.0 else f", times {factor:.2f}"


def seismic_stages_line(results, where):
    """The Log's line for a two-stage seismic result: each stage's factor of
    safety on ``where`` ("the critical surface", "the surface"), and how stage 2
    took its undrained strength."""
    how = SEISMIC2_STRENGTHS[results.get("seismic_strength", "interpolated")]
    how += seismic2_factor_words(results.get("seismic_factor", 1.0))
    return (f"Two-stage seismic stages on {where}: Stage 1 (static) FS = "
            f"{results['seismic_stage1_FS']:.4f}, Stage 2 (seismic, undrained strength "
            f"{how}) FS = {results['seismic_stage2_FS']:.4f}")


def _undrained_strengths(df, sigma_fc, tau_fc, label, debug_level=0, kc1=False):
    """Give every low-permeability slice (nonzero d or psi) its undrained strength
    from the consolidation stresses sigma'_fc and tau_fc of the first stage, by
    Duncan, Wright & Wong's interpolation between the K_c = 1 envelope (d, psi)
    and the K_c = K_f envelope (c', phi'), with the negative-stress rule; set
    c = tau_ff and phi = 0 on those slices of ``df`` in place. The other slices
    keep c', phi'. Shared by the rapid drawdown and the two-stage seismic
    analyses, whose second stages differ only in the loads they apply.

    With ``kc1`` the strength is the K_c = 1 envelope alone, tau_ff = d +
    sigma'_fc tan psi, with no interpolation and so no fallback (the two-stage
    seismic analysis's "kc1" choice).

    Returns (error message or None, low_k_flags, fallback_flags), the flags plain
    bools in slice order."""
    low_k_flags = [False] * len(df)
    fallback_flags = [False] * len(df)
    for i in range(len(df)):
        # Check if this slice has low-K material (d or psi are not zero)
        d_val = df.iloc[i]['d']
        psi_val = df.iloc[i]['psi']

        if d_val > 0 or psi_val > 0:
            low_k_flags[i] = True
            # Low-K material - calculate undrained strength
            if debug_level >= 2:
                print(f"Processing low-K material for slice {i+1}")

            # Get consolidation stresses for this slice
            sigma_fc_i = sigma_fc.iloc[i]
            tau_fc_i = tau_fc.iloc[i]
            phi_deg = df.iloc[i]['phi1']  # Use original phi for calculations
            c_val = df.iloc[i]['c1']      # Use original c for calculations

            phi_rad = np.radians(phi_deg)
            cos_phi = np.cos(phi_rad)
            sin_phi = np.sin(phi_rad)

            # tau_ff for the two envelopes: d-psi (Kc=1) and c'-phi' (Kc=Kf).
            tau_ff_k1 = d_val + sigma_fc_i * np.tan(np.radians(psi_val))  # d-psi curve (Kc=1)
            tau_ff_kf = c_val + sigma_fc_i * np.tan(phi_rad)             # c-phi curve (Kc=Kf)

            if kc1:
                if not np.isfinite(tau_ff_k1):
                    return (f"{label}: undrained strength is not finite for slice {i+1} "
                            f"(sigma'_fc={sigma_fc_i:.4g})."), low_k_flags, fallback_flags
                df.iloc[i, df.columns.get_loc('c')] = float(max(0.0, tau_ff_k1))
                df.iloc[i, df.columns.get_loc('phi')] = 0.0
                if debug_level >= 2:
                    print(f"  Slice {i+1}: K_c = 1 envelope -> tau_ff={max(0.0, tau_ff_k1):.3f}")
                continue

            # Fall back to the lower of the two curves (the doc's negative-stress rule)
            # whenever the K1/Kf interpolation is ill-conditioned: cos(phi) ~ 0, the Kf
            # denominator factor ~ 0, or a non-positive minor principal stress sigma'_3c
            # on either envelope (eqs 7, 8). Previously these cases hit `continue` and the
            # slice silently kept its drained strength in the undrained Stage-2 solve.
            #
            # sigma'_3c (eq 7) is also the DENOMINATOR of K1 (eq 4), so the test must
            # exclude zero, not just negatives -- the source says "negative (or zero)".
            # At exactly zero K1 is +inf and tau_ff is NaN, which max(0.0, NaN) would
            # silently return as 0.0.
            kf_first = sigma_fc_i - c_val * cos_phi   # Kf denominator factor (eq 6)
            use_fallback = abs(cos_phi) < 1e-12 or abs(kf_first) < 1e-12
            if not use_fallback:
                sigma3_k1 = sigma_fc_i + tau_fc_i * (sin_phi - 1) / cos_phi          # eq (7)
                sigma3_kf = kf_first * (1 - sin_phi) / (cos_phi ** 2)                # eq (8)
                use_fallback = sigma3_k1 <= 0 or sigma3_kf <= 0

            if use_fallback:
                fallback_flags[i] = True
                tau_ff = min(tau_ff_k1, tau_ff_kf)
            else:
                K1 = (sigma_fc_i + tau_fc_i * (sin_phi + 1) / cos_phi) / \
                     (sigma_fc_i + tau_fc_i * (sin_phi - 1) / cos_phi)               # eq (4)
                Kf = ((sigma_fc_i + c_val * cos_phi) * (1 + sin_phi)) / \
                     (kf_first * (1 - sin_phi))                                       # eq (6)
                if abs(Kf - 1) < 1e-12:
                    tau_ff = tau_ff_k1
                else:
                    tau_ff = ((Kf - K1) * tau_ff_k1 + (K1 - 1) * tau_ff_kf) / (Kf - 1)  # eq (5)

            # The Kc=1 and Kc=Kf envelopes bound the physically possible states, so
            # eq (5) interpolates and never extrapolates. The Stage-1 FS >= 1 guard of
            # the callers already assures K1 <= Kf; a non-finite tau_ff would mean that
            # reasoning failed, and must not be laundered into 0.0 by the clamp below.
            if not np.isfinite(tau_ff):
                return (f"{label}: undrained strength is not finite for slice {i+1} "
                        f"(sigma'_fc={sigma_fc_i:.4g}, tau_fc={tau_fc_i:.4g}). The K1/Kf "
                        "interpolation is degenerate."), low_k_flags, fallback_flags

            tau_ff = max(0.0, tau_ff)   # undrained shear strength cannot be negative

            if debug_level >= 2:
                print(f"  Slice {i+1}: tau_ff_k1={tau_ff_k1:.3f}, tau_ff_kf={tau_ff_kf:.3f}, "
                      f"fallback={use_fallback} -> tau_ff={tau_ff:.3f}")

            # Set undrained strength parameters
            df.iloc[i, df.columns.get_loc('c')] = float(tau_ff)
            df.iloc[i, df.columns.get_loc('phi')] = 0.0
        else:
            # High-K material - keep original c and phi
            if debug_level >= 2:
                print(f"Slice {i+1}: High-K material, keeping original c and phi")
    return None, low_k_flags, fallback_flags


def seismic_two_stage(df, method_name, debug_level=1, strength="interpolated"):
    """The two-stage pseudostatic seismic analysis of Duncan, Wright & Brandon
    (2014), Sec. 10.4: the rapid drawdown procedure with the earthquake in place of
    the drawdown, for slopes that have reached consolidated equilibrium.

    Stage 1 is the static slope before the earthquake: the slice table as given,
    with its pore pressures and loads, drained strengths (c', phi'), and NO
    seismic force, which gives the consolidation stresses sigma'_fc and tau_fc on
    each slice base. Stage 2 gives every slice in a material with d and psi its
    undrained strength from those stresses (the interpolation of
    :func:`_undrained_strengths`, the same as rapid drawdown's stage 2), keeps
    c', phi' on the free-draining slices, keeps the same pore pressures and loads,
    and adds the seismic force kW: its factor of safety is the answer. There is
    no third stage, because drainage during an earthquake is unlikely (Sec. 10.4).

    Duncan, Wright & Brandon (2014, Sec. 10.4.2) take the stage 2 undrained
    strength "using the same procedures as for rapid drawdown", which is the
    default here. They computed their published two-stage factors of safety
    (four soils on one slope at k = 0.15) from the K_c = 1 envelope alone,
    tau_ff = d + sigma'_fc tan psi at the stage 1 stresses; ``strength="kc1"``
    does that.

    Parameters
    ----------
    df : pandas.DataFrame
        The slice table, with its ``kw`` column (the seismic force on each slice,
        from the model's k), ``d``, ``psi``, and the drained ``c1``, ``phi1``.
    method_name : str
        The limit equilibrium method ('spencer', 'bishop', ...).
    debug_level : int
        0: no output, 1: each stage's factor of safety, >1: per-slice detail.
    strength : str, True, or dict
        How stage 2 takes the undrained strength: ``"interpolated"`` (or True, the
        default) or ``"kc1"``; see :func:`seismic2_strength`. A dict
        ``{"strength": ..., "factor": f}`` also multiplies every undrained
        strength in stage 2 by f (:func:`seismic2_factor`), the screening methods'
        strength reduction: the low-permeability slices' tau_ff and the c/p slices'
        strength (``cp_flag``). The free-draining slices keep their drained
        strength, since the reduction is on undrained strength (Duncan, Wright &
        Brandon 2014, p. 182).

    Returns
    -------
    (bool, dict or str)
        (True, the stage 2 solver result with the stage record added) or
        (False, an error message). The record:

        - ``seismic_stage1_FS``: the static factor of safety (stage 1).
        - ``seismic_stage2_FS``: the pseudostatic factor of safety (stage 2), the
          same as ``FS``.
        - ``seismic_low_k``: per slice, True where the base material has d or psi
          and takes an undrained strength in stage 2.
        - ``seismic_fallback``: per slice, True where a consolidation stress came
          out zero or negative and stage 2 used the lower of the two envelopes
          (always False with ``"kc1"``).
        - ``seismic_strength``: ``"interpolated"`` or ``"kc1"``.
        - ``seismic_factor``: the stage 2 strength factor (1 unless one was given).

        The keys are distinct from rapid drawdown's, so nothing that reads a
        rapid drawdown result mistakes this for one.
    """
    for _flag, _opt, _name in (('pow_flag', 'pow', 'power-curve'),
                               ('hb_flag', 'hb', 'Hoek-Brown')):
        if _flag in df.columns and bool(df[_flag].any()):
            return False, (f"One or more slices use the {_name} strength option "
                           f"(option='{_opt}'), which is not supported in the "
                           "two-stage seismic analysis.")

    factor = seismic2_factor(strength)
    strength = seismic2_strength(strength)
    if strength is None:
        return False, "seismic_two_stage needs a stage 2 strength: 'interpolated' or 'kc1'."

    from . import solve
    method_func = getattr(solve, method_name)

    caller_df = df
    df = df.copy()
    kw = df['kw'].copy()
    kvw = df['kvw'].copy() if 'kvw' in df.columns else None

    if debug_level >= 1:
        print("=== TWO-STAGE SEISMIC ANALYSIS ===")
        print("Stage 1: static, before the earthquake...")
    df['kw'] = 0.0
    if kvw is not None:
        df['kvw'] = 0.0
    success, result_stage1 = method_func(df)
    if not success:
        return False, f"Stage 1 failed: {result_stage1}"
    stage1_FS = result_stage1['FS']
    if debug_level >= 1:
        print(f"Stage 1 FS = {stage1_FS:.4f}")

    # As in rapid drawdown, stage 1 only supplies the consolidation stresses, and
    # they exist only for a slope that stands: at FS < 1 the mobilized shear lies
    # above the failure envelope (K1 > Kf) and the interpolation would extrapolate.
    if stage1_FS < 1.0:
        return False, (
            f"The two-stage seismic analysis requires a slope that is stable before the "
            f"earthquake, but the Stage 1 (static) factor of safety is {stage1_FS:.4f} < 1. "
            f"The consolidation stresses are undefined because the stress state lies "
            f"above the failure envelope (K1 > Kf).")
    if 'n_eff' not in df.columns:
        return False, "Stage 1 did not compute n_eff values"

    sigma_fc = df['n_eff'] / df['dl']                                              # eq (2)
    tau_fc = (1.0 / stage1_FS) * (df['c'] + sigma_fc * np.tan(np.radians(df['phi'])))  # eq (3)

    if debug_level >= 1:
        print("Stage 2: the earthquake, with undrained strengths...")
    err, low_k_flags, fallback_flags = _undrained_strengths(
        df, sigma_fc, tau_fc, "Two-stage seismic", debug_level, kc1=(strength == "kc1"))
    if err:
        return False, err
    if factor != 1.0:
        # Every undrained strength in stage 2 (471a): the low-permeability slices,
        # which now carry c = tau_ff, phi = 0, and the c/p slices, whose undrained
        # strength is the same in both stages. Neither has a friction angle here.
        low = np.asarray(low_k_flags, dtype=bool)
        if 'cp_flag' in df.columns:
            low = low | df['cp_flag'].to_numpy(dtype=bool)
        df.loc[low, 'c'] = df.loc[low, 'c'].to_numpy() * factor
    df['kw'] = kw.values
    if kvw is not None:
        df['kvw'] = kvw.values
    success, result = method_func(df)
    if not success:
        return False, f"Stage 2 failed: {result}"
    stage2_FS = result['FS']
    if debug_level >= 1:
        print(f"Stage 2 FS = {stage2_FS:.4f}")
        print("=== END TWO-STAGE SEISMIC ANALYSIS ===")

    for col in ("n_eff", "u", "c", "phi"):
        if col in df.columns and col in caller_df.columns:
            caller_df[col] = df[col].values

    result['seismic_stage1_FS'] = stage1_FS
    result['seismic_stage2_FS'] = stage2_FS
    result['seismic_low_k'] = low_k_flags
    result['seismic_fallback'] = fallback_flags
    result['seismic_strength'] = strength
    result['seismic_factor'] = factor
    return True, result


def rapid_drawdown(df, method_name, debug_level=1):
    """
    Performs rapid drawdown analysis using a three-stage approach.
    
    Parameters:
        df: pandas.DataFrame
            Slice data with all required columns including rapid drawdown specific data:
            - c, phi: current strength parameters
            - c1, phi1: original strength parameters (for stage 3)
            - d, psi: rapid drawdown parameters for low-K materials
            - u: pore pressure (stage 1)
            - u2: pore pressure for lowered pool (stage 2)
            - dload, d_x, d_y: distributed loads (stage 1)
            - dload2, d_x2, d_y2: distributed loads for lowered pool (stage 2)
        method_name: str
            The method name to use ('oms', 'bishop', 'spencer', etc.)
        debug_level: int
            0: no output, 1: print FS at each stage, >1: detailed debug info
    
    Returns:
        Tuple(bool, dict): (True, result_dict) or (False, error_message)

        result_dict is the governing stage's solver result, with the stage
        record added:
            - stage1_FS, stage2_FS, stage3_FS: the stage factors of safety.
              When Stage 3 is not needed, stage3_FS repeats stage2_FS.
            - stage3_run: True when Stage 3 was solved.
            - rapid_low_k: per slice, True where the base material has a
              nonzero d or psi and so takes an undrained strength in Stage 2.
            - stage2_fallback: per slice, True where a consolidation stress came
              out zero or negative and Stage 2 used the lower of the two
              envelopes instead of interpolating between them.
            - stage3_drained: per slice, True where the drained strength was
              lower and replaced the undrained one in Stage 3.
        The per-slice lists are plain bools in slice-table order. A surface that
        crosses no material with d or psi is a drained analysis in every stage:
        rapid_low_k is all False, stage3_run is False, and the factor of safety
        is Stage 2's (see :func:`crosses_no_low_k`).
    """
    # Curved-envelope materials are incompatible with the staged strength
    # overrides of the drawdown procedure (both mutate per-slice c/phi); refuse
    # clearly rather than silently letting one clobber the other.
    for _flag, _opt, _name in (('pow_flag', 'pow', 'power-curve'),
                               ('hb_flag', 'hb', 'Hoek-Brown')):
        if _flag in df.columns and bool(df[_flag].any()):
            return False, (f"One or more slices use the {_name} strength option "
                           f"(option='{_opt}'), which is not supported in rapid "
                           "drawdown analysis.")

    # Import solve module and get the method function
    from . import solve
    method_func = getattr(solve, method_name)

    # Work on a copy: the analysis overwrites strength and load columns (Stage 2
    # swaps in the drawdown pore pressures / loads and undrained strengths), so do
    # not mutate the caller's slice DataFrame. We DO write the winning stage's
    # per-slice results (n_eff, u, c, phi) back to it at the end, so the caller —
    # and the plots/search cache that read it — see the rapid-drawdown stresses
    # rather than the stale Stage-0 values.
    caller_df = df
    df = df.copy()

    # A surface that crosses no low-permeability material (no slice with a
    # nonzero d or psi) is still a surface: every material it cuts drains during
    # the drawdown, so its answer is the drained analysis under the water after
    # drawdown. It runs through the same three stages — Stage 2 finds no slice to
    # give an undrained strength and keeps c', phi' throughout, and Stage 3 has
    # nothing to substitute, so it is not required — and the governing factor is
    # Stage 2's. Refusing it, as this function once did, made a search discard
    # such circles as unsolved trials and miss a critical surface in a weak
    # free-draining shell. Whether the MODEL carries any d/psi material at all
    # is a separate question, refused up front by validate_rapid_drawdown and by
    # preflight (rapid.no_low_k_material).
    if debug_level >= 1:
        print("=== RAPID DRAWDOWN ANALYSIS ===")
    
    # Stage 1: Pre-drawdown conditions
    if debug_level >= 1:
        print("Stage 1: Pre-drawdown conditions...")
    
    # Use original conditions (c, phi, u, dload, d_x, d_y)
    success, result_stage1 = method_func(df)
    if not success:
        return False, f"Stage 1 failed: {result_stage1}"
    
    stage1_FS = result_stage1['FS']
    if debug_level >= 1:
        print(f"Stage 1 FS = {stage1_FS:.4f}")

    # The three-stage procedure presumes the slope is stable BEFORE drawdown: stage 1
    # exists only to supply consolidation stresses. If FS < 1 the mobilized shear
    # tau_fc = (1/FS)(c' + sigma'_fc tan phi') exceeds the failure envelope, so the
    # consolidation stress state lies above it and K1 > Kf. Since K1 rises
    # monotonically with tau_fc (hence with 1/FS) and equals Kf exactly at FS = 1,
    # `stage1_FS < 1` is precisely the condition `K1 > Kf`, uniformly over slices.
    #
    # Eq (5) then EXTRAPOLATES beyond the Kc=Kf envelope -- which the source defines as
    # the physical extreme -- driving tau_ff below it and eventually negative, where the
    # max(0, tau_ff) clamp below turns it into a silent zero-strength slice. A search
    # calling this on trial surfaces would let such a surface win on a fictitious FS~0.
    # The negative-stress fallback does NOT catch this: sigma'_3c stays positive.
    if stage1_FS < 1.0:
        return False, (
            f"Rapid drawdown requires a slope that is stable before drawdown, but the "
            f"Stage 1 (full pool) factor of safety is {stage1_FS:.4f} < 1. The "
            f"consolidation stresses are undefined because the stress state lies above "
            f"the failure envelope (K1 > Kf)."
        )

    # Calculate consolidation stresses for each slice
    # N_eff should be available from the method function
    if 'n_eff' not in df.columns:
        return False, "Stage 1 did not compute n_eff values"
    
    # Calculate sigma_fc and tau_fc for each slice
    sigma_fc = df['n_eff'] / df['dl']  # Equation (2)
    tau_fc = (1.0 / stage1_FS) * (df['c'] + sigma_fc * np.tan(np.radians(df['phi'])))  # Equation (3)
    
    if debug_level >= 2:
        print("Stage 1 consolidation stresses:")
        for i in range(len(df)):
            print(f"  Slice {i+1}: sigma_fc = {sigma_fc.iloc[i]:.2f}, tau_fc = {tau_fc.iloc[i]:.2f}")
    
    # Stage 2: Post-drawdown conditions with undrained strengths
    if debug_level >= 1:
        if ((df['d'] > 0) | (df['psi'] > 0)).any():
            print("Stage 2: Post-drawdown conditions with undrained strengths...")
        else:
            print("Stage 2: Post-drawdown conditions; no slice is in a material "
                  "with d and ψ, so every slice keeps c′, φ′...")
   
    # Update pore pressures and distributed loads for stage 2
    df['u'] = df['u2']
    df['dload'] = df['dload2']
    df['d_x'] = df['d_x2']
    df['d_y'] = df['d_y2']
    if 'beta2' in df.columns:
        # The stage-2 load carries its OWN inclination: the two dloads sheets declare
        # their directions independently, so a stage-1 vertical surcharge beside a
        # stage-2 normal pool (or the reverse) must not inherit the other's beta.
        df['beta'] = df['beta2']
    
    # Per-slice record of what each stage did, in slice order, for the report:
    # which slices are low-permeability (they take an undrained strength in
    # Stage 2), which of those fell back to the lower envelope under the
    # negative-stress rule, and which took the drained strength in Stage 3.
    # Plain Python bools, so a result dict stays picklable and JSON-safe.
    drained_flags = [False] * len(df)
    err, low_k_flags, fallback_flags = _undrained_strengths(
        df, sigma_fc, tau_fc, "Rapid drawdown", debug_level)
    if err:
        return False, err

    # Calculate Stage 2 FS
    success, result_stage2 = method_func(df)
    if not success:
        return False, f"Stage 2 failed: {result_stage2}"
    
    stage2_FS = result_stage2['FS']
    if debug_level >= 1:
        print(f"Stage 2 FS = {stage2_FS:.4f}")

    stage2_state = df.copy()   # per-slice n_eff/u/c/phi for the Stage-2 result

    # Stage 3: Check drained strengths
    if debug_level >= 1:
        print("Stage 3: Checking drained strengths...")
    
    # Check if any low-K slices need drained strength
    need_stage3 = False
    
    for i in range(len(df)):
        d_val = df.iloc[i]['d']
        psi_val = df.iloc[i]['psi']
        
        if d_val > 0 or psi_val > 0:
            # This is a low-K material slice
            if 'n_eff' not in df.columns:
                return False, "Stage 2 did not compute n_eff values"
            
            # Calculate drained strength using equations (9) and (10)
            sigma_prime = df.iloc[i]['n_eff'] / df.iloc[i]['dl']  # Equation (9)
            tau_drained = df.iloc[i]['c1'] + sigma_prime * np.tan(np.radians(df.iloc[i]['phi1']))  # Equation (10)
            
            # Compare with undrained strength (current c value)
            tau_undrained = df.iloc[i]['c']
            
            if debug_level >= 2:
                print(f"Slice {i+1}: tau_drained = {tau_drained:.4f}, tau_undrained = {tau_undrained:.4f}")
            
            if tau_drained < tau_undrained:
                # Use drained strength
                df.iloc[i, df.columns.get_loc('c')] = float(df.iloc[i]['c1'])
                df.iloc[i, df.columns.get_loc('phi')] = float(df.iloc[i]['phi1'])
                need_stage3 = True
                drained_flags[i] = True
                
                if debug_level >= 2:
                    print(f"  Using drained strength for slice {i+1}")
    
    if need_stage3:
        if debug_level >= 1:
            print("Stage 3: Recalculating FS with drained strengths...")
        
        success, result_stage3 = method_func(df)
        if not success:
            return False, f"Stage 3 failed: {result_stage3}"
        
        stage3_FS = result_stage3['FS']
        if debug_level >= 1:
            print(f"Stage 3 FS = {stage3_FS:.4f}")
    else:
        stage3_FS = stage2_FS
        result_stage3 = result_stage2
        if debug_level >= 1:
            print("Stage 3: No drained strength adjustments needed")

    stage3_state = df.copy()   # per-slice state after Stage 3 (== Stage 2 if skipped)

    # Final FS is the lower of Stage 2 and Stage 3
    if stage2_FS < stage3_FS:
        final_FS = stage2_FS
        result = result_stage2
        winning_state = stage2_state
    else:
        final_FS = stage3_FS
        result = result_stage3
        winning_state = stage3_state

    # Hand the winning stage's per-slice results back to the caller's DataFrame so
    # the cached/plotted slices carry real rapid-drawdown stresses (not Stage-0).
    for col in ("n_eff", "u", "c", "phi"):
        if col in winning_state.columns and col in caller_df.columns:
            caller_df[col] = winning_state[col].values
    
    if debug_level >= 1:
        print(f"Final rapid drawdown FS = {final_FS:.4f}")
        print("=== END RAPID DRAWDOWN ANALYSIS ===")
    
    # Append stage FS to result
    result['stage1_FS'] = stage1_FS
    result['stage2_FS'] = stage2_FS
    # stage3_FS repeats stage2_FS when Stage 3 was not needed (other code reads
    # it that way); stage3_run is what says whether Stage 3 was actually solved.
    result['stage3_FS'] = stage3_FS
    result['stage3_run'] = bool(need_stage3)
    # Per-slice flags, one entry per slice in slice-table order.
    result['rapid_low_k'] = low_k_flags
    result['stage2_fallback'] = fallback_flags
    result['stage3_drained'] = drained_flags

    return True, result


# ---------------------------------------------------------------------------
# The reliability family moved to :mod:`xslope.reliability`. It is re-exported
# here so every historical import path keeps working unchanged, e.g.
#   from xslope.advanced import reliability, reliability_mc, reliability_fem
# Prefer importing these from xslope.reliability in new code.
# ---------------------------------------------------------------------------
from .reliability import (  # noqa: F401  (re-exported for backward compatibility)
    reliability,
    reliability_taylor,
    reliability_mc,
    reliability_rs,
    reliability_fem,
    MC_DEFAULT_SEED,
    MC_DEFAULT_SAMPLES,
    _MC_FLOOR,
    _strength_param_mapping,
    _reliability_param_info,
    _perturbed_slope_data,
    _finalize_reliability,
    _mc_param_info,
    _mc_sampled_slope_data,
)
