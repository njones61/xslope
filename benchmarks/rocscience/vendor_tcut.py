"""RS2 vendor tensile-strength caps (``t_cut``) for the RS2 benchmark corpus.

Every RS2 verification model assigns each Mohr-Coulomb material a TENSILE STRENGTH
(the .fez ``Tensile Strength`` field, read by ``xslope.rs2.read_fez`` as ``t_cut``)
and solves its SSR with that cap in place.  The corpus builders originally dropped
it, which is not "no cap": an uncapped Mohr-Coulomb material carries an implicit
tensile strength of c/tan(phi) — 78 kPa where the vendor allowed 28.5, and UNBOUNDED
where phi = 0 — unreduced by the SRF.  That fictitious tension holds crest cuts shut
and inflates SSRM FS (RS2-62 was the case that exposed it; see build_rs2.rs2_62c).

The caps below are transcribed VERBATIM from the vendor .fez that each corpus file
reproduces.  Almost all of them are T = c (RS2's default when a Slide2 model is
imported), but NOT all — the exceptions are called out inline, and they are the
load-bearing ones: RS2 #061 and #067 allow NO tension at all on a cohesive soil,
#045_01/#046 do the same on their c = 200/300 phi = 0 foundations, and #066 caps
every material at a flat 50 kPa regardless of c.

A material the vendor leaves elastic ("Plasticity Specifications: None") or models
with a non-Mohr-Coulomb envelope (power curve, generalized Hoek-Brown) carries no
tensile field and gets no cap here — ``None`` entries record that the material WAS
checked against the vendor model, not that it was overlooked.

Whether the cap is reduced along with c and tan(phi) during the SSR is a separate
vendor switch (``tensilestrength_SRF``); it rides on the doc-page test tag as
``tension_srf``, not here.  It is 1 on every Slide2-Import model and on RS2 #059-#068,
and 0 on the rest of the main set.

Caps are per-FILE, never inherited: most builders start from a loaded base file
(``xslope_acads_simple.xlsx``) and copy its material dict, so
``apply_vendor_t_cut`` CLEARS every material's t_cut before applying this table.
That is why the table must list every corpus file that carries a cap — a file
missing from it is written uncapped.

Sources: RS2_Slope-Stability-Verification-RS2-and-Slide2-Import (264 .fez).
"""

import os
import sys

# {output file name: {material name: cap in the file's own stress unit, or None}}
VENDOR_T_CUT = {
    # RS2-1
    #   LEM sample file; the RS2-1 SSRM tag runs on it. t_cut is FEM-only,
    #   so the LEM locks on this file are untouched.
    'xslope_acads_simple.xlsx': {
        'Soil': 3.0,
    },
    # RS2-2
    'vp003.xlsx': {
        'Soil #1': 0.0,
        'Soil #2': 5.3,
        'Soil #3': 7.2,
    },
    # RS2-3
    'vp004.xlsx': {
        'Soil #1': 0.0,
        'Soil #2': 5.3,
        'Soil #3': 7.2,
    },
    # RS2-4
    'vp005.xlsx': {
        'Rockfill': 0.0,
        'Transitions': 0.0,
        'Filter': 0.0,
        'Core': 85.0,
    },
    # RS2-5
    #   LEM sample file (RS2-5 SSRM tag). FEM-only field.
    'xslope_acads_weak_layer.xlsx': {
        'Soil 1': 28.5,
        'Weak Layer': 0.0,
    },
    # RS2-6
    'vp009.xlsx': {
        'Soil 1': 28.5,
        'Weak Layer': 0.0,
    },
    # RS2-7
    'vp010.xlsx': {
        'Soil': 11.0,
    },
    # RS2-10
    #   LEM sample file (RS2-10 SSRM tag). FEM-only field.
    'xslope_arai_tagyo.xlsx': {
        'Soil': 41.65,
    },
    # RS2-11
    'vp015.xlsx': {
        'Upper Layer': 29.4,
        'Middle Layer': 9.8,
        'Lower Layer': 294.0,
    },
    # RS2-12
    'vp016.xlsx': {
        'A&T ex3 soil': 41.65,
    },
    # RS2-13
    'vp017.xlsx': {
        'Y&U soil': 9.8,
    },
    # RS2-14
    'vp018.xlsx': {
        'Spencer 1969 soil': 10.8,
    },
    # RS2-15
    'vp019.xlsx': {
        'Upper Layer': 49.0,
        'Layer 2': 0.0,
        'Layer 3': 7.84,
        'Bottom Layer': 0.0,
    },
    # RS2-16
    'vp020.xlsx': {
        'Layer 1': 9.8,
        'Layer 2': 58.8,
        'Layer 3': 19.8,
        'Layer 4 (seam)': 9.8,
    },
    # RS2-17
    'vp021a.xlsx': {
        'F&K soil': 600.0,
    },
    # RS2-17b
    'vp021b.xlsx': {
        'F&K soil': 600.0,
    },
    # RS2-18
    'vp022a.xlsx': {
        'Upper soil': 600.0,
        'Weak soil': 0.0,
    },
    # RS2-18b
    'vp022b.xlsx': {
        'Upper soil': 600.0,
        'Weak soil': 0.0,
    },
    # RS2-19
    'vp024.xlsx': {
        'Upper Layer': 30.0,
        'Middle Layer': 20.0,
        'Bottom Layer': 150.0,
    },
    # RS2-20
    'vp025.xlsx': {
        'weightless clay': 49.0,
    },
    # RS2-21
    'vp026.xlsx': {
        'soil': 20.0,
    },
    # RS2-22
    'vp027_fem.xlsx': {
        'soil': 500.0,
    },
    # RS2-23 — '#023' rock1 carries T: 100 psf (= c at the datum); the two elastic
    #   outer zones mirror rock2, Plasticity:None, which carries no tensile field.
    'vp029_split.xlsx': {
        'Bay Mud': 100.0,
        'Bay Mud (elastic outer 1)': None,
        'Bay Mud (elastic outer 2)': None,
    },
    # RS2 #024_01 / VP32a — caps transcribed from the vendor model. The RS2 SSR
    #   rows for this problem are not attempted (the vendor mesh is split across a
    #   slip interface), so the file serves the Slide2 LEM row; t_cut is FEM-only.
    'vp032a.xlsx': {
        'Upper embankment': 0.0,
        'Lower embankment': 0.0,
        'Clay 1': 43.0,
        'Clay 2': 31.0,
        'Clay 3': 30.0,
        'Clay 4': 32.0,
        'Clay 5': 32.0,
    },
    # RS2 #024_02 / VP32c — same, for the H = 8.75 case.
    'vp032c.xlsx': {
        'Upper embankment': 0.0,
        'Lower embankment': 0.0,
        'Clay 1': 43.0,
        'Clay 2': 31.0,
        'Clay 3': 30.0,
        'Clay 4': 32.0,
        'Clay 5': 32.0,
    },
    # RS2-25 — the Syncrude tailings dyke, vendor model 'slope stability #025.fez',
    #   which states T: 0 on all five materials.
    #   The file-number join missed this one for the same reason VENDOR_E_NU records
    #   above: vp033 is named for its Slide2 problem number (VP33) while its vendor
    #   model is RS2 #025, so the join looked at #033.fez — a different problem — and
    #   correctly declined rather than transcribe the wrong constants.
    'vp033.xlsx': {
        'Tailing sand (TS)': 0.0,
        'Glacio-fluvial sand (Pf4)': 0.0,
        'Sandy till (Pgs)': 0.0,
        'Clayey till (Pgc)': 0.0,
        'Disturbed clay-shale (Kca)': 0.0,
    },
    # RS2-26
    #   PARTIAL: this row is built from the manual tables, not the .fez. Only the
    #   two fills exist in the vendor model; the two c = 0 sands have no
    #   counterpart (and a c = 0 material has zero tensile apex either way).
    'vp034.xlsx': {
        'Phase II fill': 2901.6,
        'Sand drain': None,
        'Phase I fill': 2230.0,
        'Foundation sand': None,
    },
    # RS2-27
    'vp036.xlsx': {
        'Li & Lumb soil': 18.0,
    },
    # RS2-28 — '#028_01/02/03' rock1 carries T: 10 below its own c/tan(phi) apex of
    #   12.8 kPa, so the cap binds; rock2 is "Plasticity: Non" and carries none.
    'rs2_28a.xlsx': {
        'Cut soil': 10.0,
        'Elastic outer': None,
    },
    'rs2_28b.xlsx': {
        'Cut soil': 10.0,
        'Elastic outer': None,
    },
    'rs2_28c.xlsx': {
        'Cut soil': 10.0,
        'Elastic outer': None,
    },
    # RS2-29
    'vp039c.xlsx': {
        'Fill': 0.0,
        'Soft Clay': 20.0,
    },
    # RS2-29 clay case — '#029_clay' caps both materials at T: 20 (= c, phi = 0, so
    #   an uncapped material would carry UNBOUNDED tension). The vendor also drops the
    #   cap to Tr: 0 on tensile failure; XSLOPE's cap is constant, so the file carries
    #   the peak value and the section brackets the brittle drop.
    'rs2_29clay.xlsx': {
        'Clay Fill': 20.0,
        'Soft Clay': 20.0,
    },
    # RS2-31a
    'vp044b.xlsx': {
        'clay (MC)': 11.64,
    },
    # RS2-31b
    'vp044c.xlsx': {
        'clay (LLA)': 0.39,
    },
    # RS2-32
    'vp045a.xlsx': {
        'clay (MC)': 11.64,
    },
    # RS2-33
    'vp056.xlsx': {
        'Sandy clay': 300.0,
    },
    # RS2-34
    'vp061b.xlsx': {
        'London clay (MC)': 6.0,
    },
    # RS2-36a
    'vp071a.xlsx': {
        'Material 1': 200.0,
    },
    # RS2-36b
    'vp071b.xlsx': {
        'Material 1': 200.0,
    },
    # RS2-38
    #   from #038.fez (tensilestrength_SRF = 0). The #038-1/#038-2 sub-models
    #   in the same section DO set the SRF switch, but this row is not built
    #   from them — hence no tension_srf on the RS2-38 tag.
    'vp074.xlsx': {
        'Sand': 0.0,
        'Saturated clay': 2500.0,
    },
    # RS2-41
    'vp079.xlsx': {
        'Embankment': 0.0,
        'Foundation': 450.0,
    },
    # RS2-43
    'vp081.xlsx': {
        'Embankment': 0.0,
        'Foundation': 500.0,
    },
    # RS2-42
    'vp075.xlsx': {
        'Fill': 0.0,
        'Clay crust': 41.0,
        'Marine clay': 34.5,
        'Lacustrine clay': 31.2,
    },
    # RS2-44
    'vp082.xlsx': {
        'Embankment': 600.0,
        'Foundation': 0.0,
    },
    # RS2-45a
    #   T != c: #045_01 caps the c = 200, phi = 0 foundation at ZERO — the vendor
    #   allows it no tension at all (its uncapped MC apex is infinite).
    'vp083a.xlsx': {
        'Embankment': 0.0,
        'Foundation': 0.0,
    },
    # RS2-45b
    #   #045_02 (the c = 300 case) DOES cap at T = c, unlike #045_01 above.
    'vp083b.xlsx': {
        'Embankment': 0.0,
        'Foundation': 300.0,
    },
    # RS2-46a
    #   T != c: #046 caps the c = 300, phi = 0 foundation at ZERO, in all four
    #   cz cases (verified identical across #046_01..04).
    'vp084a.xlsx': {
        'Embankment': 0.0,
        'Foundation': 0.0,
    },
    # RS2-46b
    'vp084b.xlsx': {
        'Embankment': 0.0,
        'Foundation': 0.0,
    },
    # RS2-46c
    'vp084c.xlsx': {
        'Embankment': 0.0,
        'Foundation': 0.0,
    },
    # RS2-46d
    'vp084d.xlsx': {
        'Embankment': 0.0,
        'Foundation': 0.0,
    },
    # RS2-40 — vendor models 'slope stability #040_01/02.fez' state T: 0 on both
    #   materials (C: 0 phi: 38 for the Shell, C: 0 phi: 20 for the Core).
    #   vp077b is authored from the Slide2 manual's figures rather than from a .fez,
    #   so the vendor-first pass never reached it.
    'vp077b.xlsx': {
        'Shell': 0.0,
        'Core': 0.0,
    },
    # RS2-47
    'vp078.xlsx': {
        'Material 1': 1000.0,
    },
    # RS2-47b
    'vp078b.xlsx': {
        'Material 1': 1000.0,
    },
    # RS2-47c
    'vp078c.xlsx': {
        'Material 1': 1000.0,
    },
    # RS2 #048 / VP87 — caps transcribed from the vendor model. The RS2 SSR row
    #   for this wall is not attempted (the vendor mesh is split at the sheets and
    #   joined by slip interfaces), so the file serves the Slide2 LEM row.
    'vp087.xlsx': {
        'Reinforced and retained fill': 0.0,
        'Foundation soil': 10.0,
        'Blocks': 2.5,
    },
    # RS2-52 (the SSRM file; vendor #052 gives the weak foundation c = 0, so its
    # T = 0 too — the only live cap on this variant is the blocks').
    'vp091_fem.xlsx': {
        'Reinforced and retained fill': 0.0,
        'Foundation soil': 0.0,
        'Blocks': 2.5,
    },
    # RS2-56a
    'rs2_56a.xlsx': {
        'soil': 5.0,
    },
    # RS2-56b
    'rs2_56b.xlsx': {
        'soil': 20.0,
    },
    # RS2-57a
    'rs2_57a.xlsx': {
        'soil': 5.0,
    },
    # RS2-57b
    'rs2_57b.xlsx': {
        'soil': 20.0,
    },
    # RS2-58a
    'rs2_58a.xlsx': {
        'soil': 5.0,
    },
    # RS2-56/57/58, the eleven non-bracket cases. Vendor T = c on all seventeen
    # Pruska models (a real cap: for case 5 of #56 the uncapped apex is
    # 20/tan30 = 34.6 kPa against the vendor's 20), tensilestrength_SRF: 0.
    'rs2_56c1.xlsx': {'soil': 20.0},
    'rs2_56c3.xlsx': {'soil': 20.0},
    'rs2_56c4.xlsx': {'soil': 5.0},
    'rs2_57c2.xlsx': {'soil': 20.0},
    'rs2_57c3.xlsx': {'soil': 5.0},
    'rs2_57c4.xlsx': {'soil': 20.0},
    'rs2_57c5.xlsx': {'soil': 5.0},
    'rs2_58c2.xlsx': {'soil': 20.0},
    'rs2_58c3.xlsx': {'soil': 5.0},
    'rs2_58c4.xlsx': {'soil': 20.0},
    'rs2_58c5.xlsx': {'soil': 5.0},
    # RS2-58b
    'rs2_58b.xlsx': {
        'soil': 20.0,
    },
    # RS2-59
    'rs2_59.xlsx': {
        'YellowClay/Debris': 50.0,
        'Waste': 1.0,
        'GreyClay': 250.0,
    },
    # RS2-61-case2
    #   T != c and LOAD-BEARING: #061 allows NO tension on a c = 5, phi = 30 soil
    #   (uncapped MC apex 8.66).
    'rs2_61a.xlsx': {
        'soil': 0.0,
    },
    # RS2-63
    'rs2_63.xlsx': {
        'soil': 10.0,
    },
    # RS2-64a
    #   RS2 #064 cases 1-12 all use T = c; the blocked cases (f/h/i/j/l) and the
    #   material-partition split files carry the same cap on their rock1/soil.
    'rs2_64a.xlsx': {
        'soil': 40.9,
    },
    # RS2-64c
    'rs2_64c.xlsx': {
        'soil': 33.6,
    },
    # RS2-64e
    'rs2_64e.xlsx': {
        'soil': 33.6,
    },
    # RS2-64b
    'rs2_64b.xlsx': {
        'soil': 27.8,
    },
    # RS2-64d
    'rs2_64d.xlsx': {
        'soil': 28.4,
    },
    # RS2-64g
    'rs2_64g.xlsx': {
        'soil': 12.0,
    },
    # RS2-64k
    'rs2_64k.xlsx': {
        'soil': 7.0,
    },
    # RS2-65
    'rs2_65.xlsx': {
        'Rockfill Lyulyaka': 20.0,
        'Fill': 22.5,
        'Rockfill G.Sakar': 20.0,
        'Counterfill': 22.5,
        'Tailings': 0.0,
        'Alluvial Clay': 0.0,
        'Marly Clay': 0.0,
        'Marl': 30.0,
    },
    # RS2-66a
    #   T != c: #066 caps EVERY material at a flat 50 kPa regardless of c
    #   (identical across all five h cases). Non-binding on the c = 0 fill,
    #   binding on both phi = 0 clays, whose uncapped apex is infinite.
    'rs2_66a.xlsx': {
        'embankment': 50.0,
        'soft ground': 50.0,
        'bearing stratum': 50.0,
    },
    # RS2-66b
    'rs2_66b.xlsx': {
        'embankment': 50.0,
        'soft ground': 50.0,
        'bearing stratum': 50.0,
    },
    # RS2-66c
    'rs2_66c.xlsx': {
        'embankment': 50.0,
        'soft ground': 50.0,
        'bearing stratum': 50.0,
    },
    # RS2-66d
    'rs2_66d.xlsx': {
        'embankment': 50.0,
        'soft ground': 50.0,
        'bearing stratum': 50.0,
    },
    # RS2-66e
    'rs2_66e.xlsx': {
        'embankment': 50.0,
        'soft ground': 50.0,
        'bearing stratum': 50.0,
    },
    # RS2-67a
    #   T != c and LOAD-BEARING: #067 allows NO tension on the c = 13.8,
    #   phi = 37 dam material (uncapped MC apex 18.31).
    'rs2_67a.xlsx': {
        'rock1': 0.0,
    },
    # RS2-67c
    'rs2_67c.xlsx': {
        'rock1': 0.0,
    },
    # RS2-67d
    'rs2_67d.xlsx': {
        'rock1': 0.0,
    },
    # RS2-67b
    'rs2_67b.xlsx': {
        'rock1': 0.0,
    },
    # RS2-67e
    'rs2_67e.xlsx': {
        'rock1': 0.0,
    },
    # RS2-67f
    'rs2_67f.xlsx': {
        'rock1': 0.0,
    },
    # RS2-9
    #   Cubzac-les-Ponts. rock4 is the elastic face skin ("Plasticity Specifications:
    #   None") and carries no cap; tensilestrength_SRF = 0, so the caps below are
    #   static through the SSR (the tag runs tension_srf=false).
    'rs2_9.xlsx': {
        'Embankment': 0.0,
        'Embankment (elastic face skin)': None,
        'Upper Clay': 10.0,
        'Lower Clay': 10.0,
    },
    # RS2-P4-VP2
    #   #002 carries a SECOND c=32/phi=10 material with T = 0 (rock2) filling a
    #   near-surface strip over the T = 32 body (rock1) — how RS2 imports Slide2's
    #   tension crack. BOTH zones transcribe; the builder lays the strip as its own
    #   material. Same pattern in vp057/vp060/vp064.
    'vp002.xlsx': {
        'ACADS 1(b) soil (crack zone)': 0.0,
        'ACADS 1(b) soil': 32.0,
    },
    # RS2-P4-VP6
    'vp006.xlsx': {
        'Rockfill': 0.0,
        'Transitions': 0.0,
        'Filter': 0.0,
        'Core': 85.0,
    },
    # RS2-P4-VP57
    #   rock3 is the T = 0 tension-crack twin of rock1 (see vp002).
    'vp057.xlsx': {
        'Sandy clay (crack zone)': 0.0,
        'Sandy clay': 300.0,
        'Highly plastic clay': 0.0,
    },
    # RS2-P4-VP60
    #   rock3 is the T = 0 tension-crack twin of rock1; 'Firm soil' is XSLOPE's
    #   stand-in for the vendor's elastic rock2 and takes no cap. The twin's zone is
    #   MEASURED and recorded rather than carried (build_problems.vp060, rs2.md):
    #   transcribed it costs 7.6% of the SSRM factor and moves the row away from
    #   RS2's own SSR, so the crack stays stated the LEM way and only rock1's cap
    #   appears here.
    'vp060.xlsx': {
        'Sandy clay': 800.0,
        'Firm soil': None,
    },
    # RS2-P4-VP64
    #   rock5 is the T = 0 tension-crack twin of rock1 (see vp002). As on vp060 the
    #   twin's crest zone is recorded rather than carried (build_problems.vp064,
    #   rs2.md): it reads 2.331 against the locked 2.369 and off RS2's own 2.37.
    'vp064.xlsx': {
        'Embankment': 1000.0,
        'Sand': 0.0,
        'Foundation Clay': 3000.0,
        'Rock': 0.0,
    },
    # RS2-P4-VP67
    'vp067.xlsx': {
        'Embankment': 1780.0,
        'Foundation': 1600.0,
    },
    # RS2-P4-VP67c
    'vp067c.xlsx': {
        'Embankment': 1780.0,
        'Foundation': 1600.0,
        'Foundation lower': 1600.0,
    },
    # RS2-P4-VP68
    'vp068.xlsx': {
        'Soil 1': 600.0,
        'Soil 2': 400.0,
        'Soil 3': 500.0,
    },
    # RS2-P4-VP69 — '#069' gives every material T: 0. Both fills are c = 0, so the
    #   cap coincides with the Mohr-Coulomb apex and is inert; transcribed anyway.
    'vp069.xlsx': {
        'Embankment': 0.0,
        'Foundation': 0.0,
    },
    # RS2-P4-VP70
    'vp070a.xlsx': {
        'Material 1': 100.0,
    },
    # RS2-P4-VP102
    'vp102a.xlsx': {
        'Material 1': 13.8,
    },
    # RS2-P4-VP102-t-60-c2, RS2-P4-VP102-t-60-c3
    'vp102t_60.xlsx': {
        'Material 1': 13.8,
    },
    # RS2-P4-VP102-t-300-c2, RS2-P4-VP102-t-300-c3
    'vp102t_300.xlsx': {
        'Material 1': 13.8,
    },
    # RS2-P4-VP102-t-1500-c2, RS2-P4-VP102-t-1500-c3
    'vp102t_1500.xlsx': {
        'Material 1': 13.8,
    },
    # (no test tag — blocked row)
    'rs2_64f.xlsx': {
        'soil': 28.4,
    },
    # (no test tag — blocked row)
    'rs2_64h.xlsx': {
        'soil': 3.0,
    },
    # (no test tag — blocked row)
    'rs2_64i.xlsx': {
        'soil': 7.0,
    },
    # (no test tag — blocked row)
    'rs2_64j.xlsx': {
        'soil': 4.0,
    },
    # (no test tag — blocked row)
    'rs2_64l.xlsx': {
        'soil': 2.0,
    },
    # (no test tag — blocked row)
    #   the elastic outer zones ('rock2*') are the vendor's Plasticity:None
    #   material and take no cap.
    'rs2_64h_split.xlsx': {
        'rock1': 3.0,
    },
    # (no test tag — blocked row)
    'rs2_64l_split.xlsx': {
        'rock1': 2.0,
    },
    # RS2-62a
    #   #062 (the row that exposed the whole defect): T = 20/0/10 kPa = c, applied
    #   with tensilestrength_SRF = 1. Restoring these turned a phantom F >= 1.3
    #   equilibrium into a crisp failure boundary at 0.769. See rs2_62c().
    'rs2_62a.xlsx': {
        'Soil 1': 20.0,
        'Soil 2 (soft band)': 0.0,
        'Soil 3': 10.0,
    },
    # RS2-62b
    'rs2_62b.xlsx': {
        'Soil 1': 20.0,
        'Soil 2 (soft band)': 0.0,
        'Soil 3': 10.0,
    },
    # RS2-62c
    'rs2_62c.xlsx': {
        'Soil 1': 20.0,
        'Soil 2 (soft band)': 0.0,
        'Soil 3': 10.0,
    },
}


def apply_vendor_t_cut(materials, path):
    """Set each material's ``t_cut`` to the vendor cap for ``path``'s output file.

    Clears t_cut on EVERY material first — a builder that starts from a loaded base
    file inherits that file's cap, which belongs to a different problem — then applies
    this file's row of VENDOR_T_CUT.  Raises if the table names a material the model
    does not have (a rename in the builder must not silently drop a cap).
    """
    names = [str(m.get('name', '')).strip() for m in materials]
    for m in materials:
        m['t_cut'] = None
    caps = VENDOR_T_CUT.get(os.path.basename(str(path)))
    if not caps:
        return
    unknown = [n for n in caps if n not in names]
    if unknown:
        raise KeyError(f'vendor_tcut: {os.path.basename(str(path))} has no material(s) '
                       f'{unknown}; model materials are {names}')
    dupes = {n for n in names if names.count(n) > 1 and n in caps}
    if dupes:
        raise KeyError(f'vendor_tcut: duplicate material name(s) {sorted(dupes)} in '
                       f'{os.path.basename(str(path))} — caps cannot be assigned by name')
    for m, n in zip(materials, names):
        if caps.get(n) is not None:
            m['t_cut'] = float(caps[n])


# ---------------------------------------------------------------------------------------
# Vendor ELASTIC constants (E, nu)
# ---------------------------------------------------------------------------------------
# Same doctrine as the caps above: a value the vendor model specifies is an INPUT, and
# inputs transcribe verbatim.  Every RS2 verification model assigns each material a
# linear-elastic pair, read off the .fez ``material types:`` block::
#
#     material 1: rock1
#      Elastic Properties: LinearElastic
#       nu: 0.4 E: 50000
#
# and solves its published SSR with those constants.  The corpus builders dropped them
# and an auto-classifier (elastic_props.py) filled the gap by soil type.  That was the
# right fix for a file with no vendor model behind it, but on a vendor-backed row it
# SUBSTITUTES for a specified input — and while the SSRM factor of safety is invariant
# to E on a plain slope, it is NOT invariant to nu once reinforcement is present
# (measured: ~4% on a reinforced model), quite apart from every displacement the FEM
# reports.
#
# So: vendor value where the vendor specifies one, classifier only where it does not.
# The pairs below come from the SAME file -> .fez -> material correspondence that
# produced VENDOR_T_CUT — the caps that mapping yields were re-derived through this
# path and matched the reviewed table exactly, so the elastic pairs ride on an
# already-reviewed mapping.  A material whose (c, phi) has no counterpart in the vendor
# model falls back to the model's uniform pair when it has one (many RS2 models give
# every zone the same constants), and is otherwise left to the classifier.
#
# NOT applied over a value the BUILDER sets explicitly: those are published constants
# transcribed from the source paper (the Pruska rs2_56/57/58 Poisson-ratio study,
# hammah_hb1's Hoek-Brown rock, rs2_59's Case-1 pair, rs2_60/61/64's E = 14000) and are
# verbatim transcriptions in their own right.  Where such a value differs from the .fez,
# the builder's docstring cites its source; see apply_vendor_e_nu.
#
# {file name: {material name: (nu, E) in the file's own units}}
VENDOR_E_NU = {
    # RS2-1
    'xslope_acads_simple.xlsx': {
        'Soil': (0.4, 50000.0),
    },
    # RS2-2
    'vp003.xlsx': {
        'Soil #1': (0.4, 50000.0),
        'Soil #2': (0.4, 50000.0),
        'Soil #3': (0.4, 50000.0),
    },
    # RS2-3
    'vp004.xlsx': {
        'Soil #1': (0.4, 50000.0),
        'Soil #2': (0.4, 50000.0),
        'Soil #3': (0.4, 50000.0),
    },
    # RS2-4
    'vp005.xlsx': {
        'Rockfill': (0.4, 50000.0),
        'Transitions': (0.4, 50000.0),
        'Filter': (0.4, 50000.0),
        'Core': (0.4, 50000.0),
    },
    # RS2-5
    'xslope_acads_weak_layer.xlsx': {
        'Soil 1': (0.4, 50000.0),
        'Weak Layer': (0.4, 50000.0),
    },
    # RS2-6
    'vp009.xlsx': {
        'Soil 1': (0.4, 50000.0),
        'Weak Layer': (0.4, 50000.0),
    },
    # RS2-7
    'vp010.xlsx': {
        'Soil': (0.4, 50000.0),
    },
    # RS2-10
    'xslope_arai_tagyo.xlsx': {
        'Soil': (0.4, 50000.0),
    },
    # RS2-11
    'vp015.xlsx': {
        'Upper Layer': (0.4, 50000.0),
        'Middle Layer': (0.4, 50000.0),
        'Lower Layer': (0.4, 50000.0),
    },
    # RS2-12
    'vp016.xlsx': {
        'A&T ex3 soil': (0.4, 50000.0),
    },
    # RS2-13
    'vp017.xlsx': {
        'Y&U soil': (0.4, 50000.0),
    },
    # RS2-14
    'vp018.xlsx': {
        'Spencer 1969 soil': (0.4, 50000.0),
    },
    # RS2-15
    'vp019.xlsx': {
        'Upper Layer': (0.4, 50000.0),
        'Layer 2': (0.4, 50000.0),
        'Layer 3': (0.4, 50000.0),
        'Bottom Layer': (0.4, 50000.0),
    },
    # RS2-16
    'vp020.xlsx': {
        'Layer 1': (0.4, 50000.0),
        'Layer 2': (0.4, 50000.0),
        'Layer 3': (0.4, 50000.0),
        'Layer 4 (seam)': (0.4, 50000.0),
    },
    # RS2-17
    'vp021a.xlsx': {
        'F&K soil': (0.4, 1000000.0),
    },
    # RS2-17b
    'vp021b.xlsx': {
        'F&K soil': (0.4, 1000000.0),
    },
    # RS2-18
    'vp022a.xlsx': {
        'Upper soil': (0.4, 1000000.0),
        'Weak soil': (0.4, 1000000.0),
    },
    # RS2-18b
    'vp022b.xlsx': {
        'Upper soil': (0.4, 1000000.0),
        'Weak soil': (0.4, 1000000.0),
    },
    # RS2-19
    'vp024.xlsx': {
        'Upper Layer': (0.4, 50000.0),
        'Middle Layer': (0.4, 50000.0),
        'Bottom Layer': (0.4, 50000.0),
    },
    # RS2-20
    'vp025.xlsx': {
        'weightless clay': (0.4, 50000.0),
    },
    # RS2-21
    'vp026.xlsx': {
        'soil': (0.4, 50000.0),
    },
    # RS2-22
    'vp027_fem.xlsx': {
        'soil': (0.4, 1000000.0),
    },
    # RS2-23 — '#023' gives rock1 and rock2 the same pair, in psf.
    'vp029_split.xlsx': {
        'Bay Mud': (0.4, 1000000.0),
        'Bay Mud (elastic outer 1)': (0.4, 1000000.0),
        'Bay Mud (elastic outer 2)': (0.4, 1000000.0),
    },
    # RS2 #024_01 / VP32a
    'vp032a.xlsx': {
        'Upper embankment': (0.4, 50000.0),
        'Lower embankment': (0.4, 50000.0),
        'Clay 1': (0.4, 50000.0),
        'Clay 2': (0.4, 50000.0),
        'Clay 3': (0.4, 50000.0),
        'Clay 4': (0.4, 50000.0),
        'Clay 5': (0.4, 50000.0),
    },
    # RS2 #024_02 / VP32c
    'vp032c.xlsx': {
        'Upper embankment': (0.4, 50000.0),
        'Lower embankment': (0.4, 50000.0),
        'Clay 1': (0.4, 50000.0),
        'Clay 2': (0.4, 50000.0),
        'Clay 3': (0.4, 50000.0),
        'Clay 4': (0.4, 50000.0),
        'Clay 5': (0.4, 50000.0),
    },
    # RS2-25 — the Syncrude tailings dyke. Its vendor model is the RS2-NATIVE
    # 'slope stability #025.fez' (five zones, gamma 20/17/17/17/17, phi
    # 34/34/34/7.5/7.5 — the dyke, matching this file zone for zone), NOT
    # '#033.fez', which is a single-material imperial slope (c = 300 psf,
    # phi = 30, gamma = 120 pcf) and a different problem altogether. The corpus
    # file is named for its Slide2 problem number (VP33) while the vendor model
    # is numbered by its RS2 problem number (#25), so the file-number join that
    # built this table looked at #033.fez, found materials that do not match, and
    # correctly declined to assign rather than transcribe the wrong model's
    # constants — which is why this row was missing until it was added by hand.
    # Rows whose two numbering spaces disagree have to be matched by content.
    'vp033.xlsx': {
        'Tailing sand (TS)': (0.4, 50000.0),
        'Glacio-fluvial sand (Pf4)': (0.4, 50000.0),
        'Sandy till (Pgs)': (0.4, 50000.0),
        'Clayey till (Pgc)': (0.4, 50000.0),
        'Disturbed clay-shale (Kca)': (0.4, 50000.0),
    },
    # RS2-26
    'vp034.xlsx': {
        'Phase II fill': (0.4, 1000000.0),
        'Sand drain': (0.4, 1000000.0),
        'Phase I fill': (0.4, 1000000.0),
        'Foundation sand': (0.4, 1000000.0),
    },
    # RS2-27
    'vp036.xlsx': {
        'Li & Lumb soil': (0.4, 50000.0),
    },
    # RS2-28 — '#028_01/02/03' give rock1 and rock2 the SAME elastic pair, so the
    #   corridor and the elastic outer zone carry no stiffness contrast.
    'rs2_28a.xlsx': {
        'Cut soil': (0.4, 50000.0),
        'Elastic outer': (0.4, 50000.0),
    },
    'rs2_28b.xlsx': {
        'Cut soil': (0.4, 50000.0),
        'Elastic outer': (0.4, 50000.0),
    },
    'rs2_28c.xlsx': {
        'Cut soil': (0.4, 50000.0),
        'Elastic outer': (0.4, 50000.0),
    },
    # RS2-29 clay case
    'rs2_29clay.xlsx': {
        'Clay Fill': (0.4, 50000.0),
        'Soft Clay': (0.4, 50000.0),
    },
    # RS2-29
    'vp039c.xlsx': {
        'Fill': (0.4, 50000.0),
        'Soft Clay': (0.4, 50000.0),
    },
    # RS2-31a
    'vp044b.xlsx': {
        'clay (MC)': (0.4, 50000.0),
    },
    # RS2-31b
    'vp044c.xlsx': {
        'clay (LLA)': (0.4, 50000.0),
    },
    # RS2-31c
    'vp044a.xlsx': {
        'clay (power)': (0.4, 50000.0),
    },
    # RS2-31d — RS2's own Generalized Hoek-Brown rendering of the same slope
    # ('slope stability #031-powecurve.fea', the model behind RS2's published 1.11)
    'vp044d.xlsx': {
        'clay (GHB)': (0.4, 50000.0),
    },
    # RS2-32
    'vp045a.xlsx': {
        'clay (MC)': (0.4, 50000.0),
    },
    # RS2-32b
    'vp045b.xlsx': {
        'clay (power)': (0.4, 50000.0),
    },
    # RS2-33
    'vp056.xlsx': {
        'Sandy clay': (0.4, 1000000.0),
    },
    # RS2-34
    'vp061b.xlsx': {
        'London clay (MC)': (0.4, 50000.0),
    },
    # RS2-34b
    'vp061a.xlsx': {
        'London clay (power)': (0.4, 50000.0),
    },
    # RS2-36a
    'vp071a.xlsx': {
        'Material 1': (0.4, 1000000.0),
    },
    # RS2-36b
    'vp071b.xlsx': {
        'Material 1': (0.4, 1000000.0),
    },
    # RS2-38
    'vp074.xlsx': {
        'Sand': (0.4, 1000000.0),
        'Saturated clay': (0.4, 1000000.0),
    },
    # RS2-41
    'vp079.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-43
    'vp081.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-42
    'vp075.xlsx': {
        'Fill': (0.4, 50000.0),
        'Clay crust': (0.4, 50000.0),
        'Marine clay': (0.4, 50000.0),
        'Lacustrine clay': (0.4, 50000.0),
    },
    # RS2-44
    'vp082.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-45a
    'vp083a.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-45b
    'vp083b.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-46a
    'vp084a.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-46b
    'vp084b.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-46c
    'vp084c.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-46d
    'vp084d.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-47
    'vp078.xlsx': {
        'Material 1': (0.4, 1000000.0),
    },
    # RS2-47b
    'vp078b.xlsx': {
        'Material 1': (0.4, 1000000.0),
    },
    # RS2-47c
    'vp078c.xlsx': {
        'Material 1': (0.4, 1000000.0),
    },
    # RS2 #048 / VP87
    'vp087.xlsx': {
        'Reinforced and retained fill': (0.4, 50000.0),
        'Foundation soil': (0.4, 50000.0),
        'Blocks': (0.4, 50000.0),
    },
    # RS2-52 — #052 gives all three materials the same pair, as #048 does.
    'vp091_fem.xlsx': {
        'Reinforced and retained fill': (0.4, 50000.0),
        'Foundation soil': (0.4, 50000.0),
        'Blocks': (0.4, 50000.0),
    },
    # RS2-56a
    'rs2_56a.xlsx': {
        'soil': (0.35, 10000.0),
    },
    # RS2-56b
    'rs2_56b.xlsx': {
        'soil': (0.3, 5000.0),
    },
    # RS2-57a
    'rs2_57a.xlsx': {
        'soil': (0.35, 10000.0),
    },
    # RS2-57b
    'rs2_57b.xlsx': {
        'soil': (0.3, 5000.0),
    },
    # RS2-58a
    'rs2_58a.xlsx': {
        'soil': (0.35, 10000.0),
    },
    # RS2-58b
    'rs2_58b.xlsx': {
        'soil': (0.3, 5000.0),
    },
    # RS2-56/57/58, the eleven non-bracket cases (see _PRUSKA_ELASTIC in build_rs2).
    'rs2_56c1.xlsx': {'soil': (0.3, 5000.0)},
    'rs2_56c3.xlsx': {'soil': (0.3, 5000.0)},
    'rs2_56c4.xlsx': {'soil': (0.35, 10000.0)},
    'rs2_57c2.xlsx': {'soil': (0.3, 5000.0)},
    'rs2_57c3.xlsx': {'soil': (0.35, 10000.0)},
    'rs2_57c4.xlsx': {'soil': (0.3, 5000.0)},
    'rs2_57c5.xlsx': {'soil': (0.35, 10000.0)},
    'rs2_58c2.xlsx': {'soil': (0.3, 5000.0)},
    'rs2_58c3.xlsx': {'soil': (0.35, 10000.0)},
    'rs2_58c4.xlsx': {'soil': (0.3, 5000.0)},
    'rs2_58c5.xlsx': {'soil': (0.35, 10000.0)},
    # RS2-61-case2
    'rs2_61a.xlsx': {
        'soil': (0.4, 50000.0),
    },
    # RS2-63
    'rs2_63.xlsx': {
        'soil': (0.3, 14000.0),
    },
    # RS2-64a
    'rs2_64a.xlsx': {
        'soil': (0.3, 50000.0),
    },
    # RS2-64c
    'rs2_64c.xlsx': {
        'soil': (0.3, 50000.0),
    },
    # RS2-64e
    'rs2_64e.xlsx': {
        'soil': (0.3, 50000.0),
    },
    # RS2-64b
    'rs2_64b.xlsx': {
        'soil': (0.4, 50000.0),
    },
    # RS2-64d
    'rs2_64d.xlsx': {
        'soil': (0.3, 50000.0),
    },
    # RS2-64g
    'rs2_64g.xlsx': {
        'soil': (0.4, 50000.0),
    },
    # RS2-64k
    'rs2_64k.xlsx': {
        'soil': (0.3, 50000.0),
    },
    # RS2-64f / h / i / j / l and the two split builds. Every #064 .fea gives BOTH
    # rock1 and the Plasticity:None rock2 the same LinearElastic pair, E = 50000
    # throughout, with nu = 0.4 on #064_02 / _07 / _12 (cases 2, 7 and 12) and 0.3
    # elsewhere — so the split files' elastic outer zones take their case's pair too.
    'rs2_64f.xlsx': {
        'soil': (0.3, 50000.0),
    },
    'rs2_64h.xlsx': {
        'soil': (0.3, 50000.0),
    },
    'rs2_64i.xlsx': {
        'soil': (0.3, 50000.0),
    },
    'rs2_64j.xlsx': {
        'soil': (0.3, 50000.0),
    },
    'rs2_64l.xlsx': {
        'soil': (0.4, 50000.0),
    },
    'rs2_64h_split.xlsx': {
        'rock1': (0.3, 50000.0),
        'rock2a': (0.3, 50000.0),
        'rock2b': (0.3, 50000.0),
        'rock2c': (0.3, 50000.0),
    },
    'rs2_64l_split.xlsx': {
        'rock1': (0.4, 50000.0),
        'rock2a': (0.4, 50000.0),
        'rock2b': (0.4, 50000.0),
    },
    # RS2-65
    'rs2_65.xlsx': {
        'Rockfill Lyulyaka': (0.3, 75000.0),
        'Fill': (0.31, 70000.0),
        'Rockfill G.Sakar': (0.3, 75000.0),
        'Counterfill': (0.31, 70000.0),
        'Tailings': (0.35, 16100.0),
        'Alluvial Clay': (0.34, 16300.0),
        'Marly Clay': (0.33, 38000.0),
        'Marl': (0.3, 75000.0),
    },
    # RS2-66a
    'rs2_66a.xlsx': {
        'embankment': (0.3, 20000.0),
        'soft ground': (0.3, 20000.0),
        'bearing stratum': (0.3, 20000.0),
    },
    # RS2-66b
    'rs2_66b.xlsx': {
        'embankment': (0.3, 20000.0),
        'soft ground': (0.3, 20000.0),
        'bearing stratum': (0.3, 20000.0),
    },
    # RS2-66c
    'rs2_66c.xlsx': {
        'embankment': (0.3, 20000.0),
        'soft ground': (0.3, 20000.0),
        'bearing stratum': (0.3, 20000.0),
    },
    # RS2-66d
    'rs2_66d.xlsx': {
        'embankment': (0.3, 20000.0),
        'soft ground': (0.3, 20000.0),
        'bearing stratum': (0.3, 20000.0),
    },
    # RS2-66e
    'rs2_66e.xlsx': {
        'embankment': (0.3, 20000.0),
        'soft ground': (0.3, 20000.0),
        'bearing stratum': (0.3, 20000.0),
    },
    # RS2-67a
    'rs2_67a.xlsx': {
        'rock1': (0.3, 100000.0),
    },
    # RS2-67c
    'rs2_67c.xlsx': {
        'rock1': (0.3, 100000.0),
    },
    # RS2-67d
    'rs2_67d.xlsx': {
        'rock1': (0.3, 100000.0),
    },
    # RS2-67b
    'rs2_67b.xlsx': {
        'rock1': (0.3, 100000.0),
    },
    # RS2-67e
    'rs2_67e.xlsx': {
        'rock1': (0.3, 100000.0),
    },
    # RS2-67f
    'rs2_67f.xlsx': {
        'rock1': (0.3, 100000.0),
    },
    # RS2-9
    'rs2_9.xlsx': {
        'Embankment': (0.4, 50000.0),
        'Embankment (elastic face skin)': (0.4, 50000.0),
        'Upper Clay': (0.4, 50000.0),
        'Lower Clay': (0.4, 50000.0),
    },
    # RS2-P4-VP2
    'vp002.xlsx': {
        'ACADS 1(b) soil (crack zone)': (0.4, 50000.0),
        'ACADS 1(b) soil': (0.4, 50000.0),
    },
    # RS2-P4-VP6
    'vp006.xlsx': {
        'Rockfill': (0.4, 50000.0),
        'Transitions': (0.4, 50000.0),
        'Filter': (0.4, 50000.0),
        'Core': (0.4, 50000.0),
    },
    # RS2-P4-VP57
    'vp057.xlsx': {
        'Sandy clay (crack zone)': (0.4, 1000000.0),
        'Sandy clay': (0.4, 1000000.0),
        'Highly plastic clay': (0.4, 1000000.0),
    },
    # RS2-P4-VP60
    'vp060.xlsx': {
        'Sandy clay': (0.4, 1000000.0),
        'Firm soil': (0.4, 1000000.0),
    },
    # RS2-P4-VP64
    'vp064.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Sand': (0.4, 1000000.0),
        'Foundation Clay': (0.4, 1000000.0),
        'Rock': (0.4, 1000000.0),
    },
    # RS2-P4-VP67
    'vp067.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-P4-VP67c
    'vp067c.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
        'Foundation lower': (0.4, 1000000.0),
    },
    # RS2-P4-VP68
    'vp068.xlsx': {
        'Soil 1': (0.4, 1000000.0),
        'Soil 2': (0.4, 1000000.0),
        'Soil 3': (0.4, 1000000.0),
    },
    # RS2-P4-VP69 — '#069' gives all eight of its materials the same pair, in psf.
    'vp069.xlsx': {
        'Embankment': (0.4, 1000000.0),
        'Foundation': (0.4, 1000000.0),
    },
    # RS2-P4-VP70
    'vp070a.xlsx': {
        'Material 1': (0.4, 1000000.0),
    },
    # RS2-P4-VP102
    'vp102a.xlsx': {
        'Material 1': (0.4, 50000.0),
    },
    # RS2-P4-VP102-t-60-c2 / RS2-P4-VP102-t-60-c3
    'vp102t_60.xlsx': {
        'Material 1': (0.4, 50000.0),
    },
    # RS2-P4-VP102-t-300-c2 / RS2-P4-VP102-t-300-c3
    'vp102t_300.xlsx': {
        'Material 1': (0.4, 50000.0),
    },
    # RS2-P4-VP102-t-1500-c2 / RS2-P4-VP102-t-1500-c3
    'vp102t_1500.xlsx': {
        'Material 1': (0.4, 50000.0),
    },
}


def apply_vendor_e_nu(materials, path):
    """Set each material's ``(nu, E)`` from VENDOR_E_NU for ``path``'s output file.

    **The vendor table wins.** Where a row exists for this output file and material,
    its pair is written unconditionally — over an unset value, over the historical
    unit-blind default, and over a literal the builder set for itself. The published
    number the row is scored against was produced by the vendor's own model with the
    vendor's own constants, so nothing in a builder can outrank them; a builder that
    disagrees is a transcription to correct, not a preference to honour. When the two
    disagree the override is announced on stderr rather than made silently, so the
    disagreement shows up in the build log instead of being discovered by an audit.

    Unlike the tensile caps this does NOT clear first: a material with no row here
    keeps whatever the builder set, and falls through to the soil-type classifier
    only if the builder said nothing.

    Returns the set of material names this call assigned. Pass it to
    ``elastic_props.assign_elastic_props(pinned=...)`` so the classifier cannot
    reclassify a vendor pair that happens to equal its "nobody set this" sentinel —
    the collision that silently rewrote all six RS2-67 files from the vendor's
    100000/0.3 to a soil-type 32000/0.4.

    Raises if the table names a material the model does not have, so a rename in a
    builder cannot silently drop the vendor constants.
    """
    from elastic_props import finite
    props = VENDOR_E_NU.get(os.path.basename(str(path)))
    if not props:
        return set()
    names = [str(m.get('name', '')).strip() for m in materials]
    unknown = [n for n in props if n not in names]
    if unknown:
        raise KeyError(f'vendor_tcut: {os.path.basename(str(path))} has no material(s) '
                       f'{unknown}; model materials are {names}')
    dupes = {n for n in names if names.count(n) > 1 and n in props}
    if dupes:
        raise KeyError(f'vendor_tcut: duplicate material name(s) {sorted(dupes)} in '
                       f'{os.path.basename(str(path))} — elastic constants cannot be '
                       f'assigned by name')
    assigned = set()
    for m, n in zip(materials, names):
        pair = props.get(n)
        if pair is None:
            continue
        nu, E = float(pair[0]), float(pair[1])
        E_old, nu_old = finite(m.get('E')), finite(m.get('nu'))
        if E_old > 0.0 and (round(E_old, 3), round(nu_old, 3)) != (round(E, 3), round(nu, 3)):
            print(f'vendor_tcut: {os.path.basename(str(path))} / {n}: builder set '
                  f'(E={E_old:g}, nu={nu_old:g}); the vendor model says '
                  f'(E={E:g}, nu={nu:g}) — using the vendor pair.', file=sys.stderr)
        m['nu'], m['E'] = nu, E
        assigned.add(n)
    return assigned
