"""Can the recommended experiment actually be performed? Checked late, and it matters.

This module exists because the project spent eight stages refining an
experimental recommendation without once asking whether the material it
requires can be made. It cannot — not as originally stated.

Two independent feasibility walls:

1. SYNTHESIS. Mirror-image proteins are made by total chemical synthesis, and
   the largest ever reported is a 358-residue D-Dpo4 polymerase (Nature 2016;
   the prior record was 312-residue DapA). The triple-helical domain of
   collagen alpha1(I) is ~1014 residues and the full pro-chain ~1464 — roughly
   3-4x beyond the record — and collagen additionally needs three chains and
   4-hydroxyproline, which chemical synthesis would have to install. **All-D
   collagen is not currently makeable.**

   Collagen-mimetic peptides (CMPs) are ~36 residues, which is routine
   solid-phase synthesis, and they traverse the same hierarchy: triple helix ->
   sticky-ended nanofibre -> hydrogel, with gels degraded by collagenase at
   rates comparable to natural collagen. So the experiment is feasible with
   D-CMP, not with D-collagen. Everything in this project should be read as
   applying to CMPs.

   What survives the substitution: the level-1 argument (both are Gly-X-Y
   triple helices, so heterochiral helices are excluded either way) and the
   level-2 question (lateral and end-to-end packing of homochiral helices).
   What changes: a CMP triple helix is ~10 nm rather than ~300 nm, so
   `lambda_c` counts many more helices for a fibril of given length.

2. MECHANICS. Native CMP hydrogels are soft — reported storage moduli run from
   tens to hundreds of Pa, i.e. 0.01-1 kPa. The cell-based discrimination
   protocol wants ~2 kPa, and it degrades badly below that: worst-case
   separation falls from 0.316 at 2 kPa to 0.045 at 1 kPa and to ~0 below
   0.5 kPa, because a cell on a very soft gel is already in the low-YAP state
   and further softening cannot move it. There is no headroom downward.

   Stiffening a CMP gel (higher concentration, crosslinking) is possible but
   carries its own risk: crosslinks that bridge broken fibrils would supply
   mechanical continuity independently of fibril integrity, which is precisely
   the signal being measured.

The conclusion is not that the project is unusable, but that its two
recommended experiments have very different standing:

  * The CELL-FREE rheology series is feasible today with D-CMP. Rheology reads
    modulus versus composition directly and does not care that the modulus is
    low. It separates the four assembly hypotheses by curve shape.
  * The CELL experiment is the problematic one and should not be attempted
    before the rheology has been done and a gel formulation reaching ~2 kPa
    without masking crosslinks has been demonstrated.
"""

from __future__ import annotations

#: Largest mirror-image protein made by total chemical synthesis (residues).
#: D-Dpo4 DNA polymerase; the previous record was 312-residue DapA.
D_PROTEIN_SYNTHESIS_RECORD = 358

#: Residue counts for candidate targets.
TARGETS = {
    "collagen alpha1(I) triple-helical domain": 1014,
    "full pro-alpha1(I) chain": 1464,
    "collagen-mimetic peptide (CMP)": 36,
}

#: Reported storage moduli for native CMP hydrogels, in kPa.
CMP_STIFFNESS_RANGE_KPA = (0.01, 1.0)


def synthesis_feasibility(record: int = D_PROTEIN_SYNTHESIS_RECORD) -> dict:
    """How far each target sits from the state of the art in D-protein synthesis."""
    rows = []
    for name, n in TARGETS.items():
        rows.append({
            "target": name,
            "residues": n,
            "ratio_to_record": n / record,
            "feasible": n <= record,
        })
    return {"record": record, "rows": rows,
            "collagen_feasible": TARGETS["collagen alpha1(I) triple-helical domain"] <= record,
            "cmp_feasible": TARGETS["collagen-mimetic peptide (CMP)"] <= record}


def cell_experiment_feasibility(p=None, max_gel_kPa: float = None) -> dict:
    """Does the cell protocol still discriminate within reachable gel stiffness?

    Compares the best separation achievable anywhere against the best achievable
    within the stiffness range a CMP hydrogel actually reaches.
    """
    from .parameters import get_preset
    from .analysis import two_point_protocol

    p = p if p is not None else get_preset("clutch_calibrated")
    cap = max_gel_kPa if max_gel_kPa is not None else CMP_STIFFNESS_RANGE_KPA[1]

    t = two_point_protocol(
        p=p, stiffnesses=(0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0))
    reachable = [r for r in t["rows"] if r["stiffness"] <= cap]
    best_reachable = max(reachable, key=lambda r: r["worst_separation"])

    return {
        "rows": t["rows"],
        "best_anywhere": t["best_separation"],
        "best_anywhere_stiffness": t["best_stiffness"],
        "best_reachable": best_reachable["worst_separation"],
        "best_reachable_stiffness": best_reachable["stiffness"],
        "penalty": (t["best_separation"] / best_reachable["worst_separation"]
                    if best_reachable["worst_separation"] > 0 else float("inf")),
        "max_gel_kPa": cap,
        "cell_experiment_viable_natively": best_reachable["worst_separation"] > 0.1,
    }
