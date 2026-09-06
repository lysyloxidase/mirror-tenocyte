#!/usr/bin/env python3
"""What published enantiomeric-mixture systems did — and why it settles little.

    ./.venv/bin/python scripts/run_evidence.py

No synthesis is available to this project, so the remaining move was to use
other people's data. This reports that survey together with its provenance and
fragility, both of which are bad enough that the survey should be treated as a
list of leads rather than as evidence.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from mirror_tenocyte.evidence import (
    SURVEY, MAJORITY_RULES_EVIDENCE, FULL_TEXT,
    prior, verified_prior, fragility, folds, unverified_magnitudes,
)

NICE = {"racemic_enhanced": "racemic enhancement (peak)",
        "self_sorting": "self-sorting (flat)",
        "frustrated_kinetic": "weakened on mixing (dip)",
        "majority_rules": "majority rules (dip at centre)"}


def rule(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def main():
    rule("1  The survey, with provenance")
    for e in SURVEY:
        flag = "" if e["verification"] == FULL_TEXT else "   <-- NOT VERIFIED"
        series = "composition series" if e["composition_series"] else "endpoints only"
        fold = f"{e['fold']:.0f}x" if e["fold"] else "—"
        print(f"  {e['system']}")
        print(f"      {e['class']}")
        print(f"      outcome: {NICE[e['outcome']]:32} {fold:>5}   [{series}]")
        print(f"      source:  {e['verification']}{flag}")
        print(f"      {e['detail']}")
        print(f"      ref: {e['ref']}")
        print()

    rule("2  Counts — deliberately not percentages")
    pr = prior()
    for k, v in sorted(pr["counts"].items(), key=lambda x: -x[1]):
        print(f"  {NICE[k]:32} {v} of {pr['n_systems']}")
    print(f"\n  full text read for {pr['n_full_text']} of {pr['n_systems']} entries.")
    print("  Percentages are withheld on purpose: with n = 6 they would imply a")
    print("  precision this survey does not have.")

    print("\n  reported racemic/enantiopure modulus ratios:")
    print("  (outcome DIRECTION and effect MAGNITUDE are sourced separately —")
    print("   an abstract can confirm one while saying nothing about the other)")
    for s, f, v, ok in folds():
        mark = "" if ok else "   <-- MAGNITUDE UNSUPPORTED by any readable source"
        print(f"      {f:5.1f}x   {s}{mark}")
    um = unverified_magnitudes()
    if um:
        print(f"\n  {len(um)} of {len(folds())} quoted magnitudes rest on search")
        print("  summaries alone. The phenylalanine abstract, reached via an")
        print("  institutional repository, confirms only that mixing gives")
        print("  'mechanically more robust' structures — it contains no numbers.")

    rule("3  Fragility — the leader flips on a single reclassification")
    fr = fragility()
    print(f"  leading outcome: {NICE[fr['leader']]}")
    print(f"  entries whose reclassification would change that leader: "
          f"{fr['n_entries_whose_change_flips_leader']} of {fr['n_entries']}")
    for s, a, b in fr["flips"]:
        print(f"      {s[:44]:46} {a} -> {b}")

    vp = verified_prior()
    print(f"\n  Restricted to entries actually read in full text: {vp['counts']}")
    if not vp["decisive"]:
        print("  => a tie. The verified subset of this survey says NOTHING about")
        print("     which outcome to expect.")

    rule("4  A retracted argument")
    print("  An earlier version of this script argued: solution phase favours")
    print("  homochiral assembly, solid-state packing favours heterochiral (the")
    print("  Abeta reconciliation, plus Wallach's rule for racemic crystals);")
    print("  collagen's level-2 lateral packing is solid-state; therefore")
    print("  `racemic_enhanced` belongs at the top of the prior.")
    print()
    print("  THAT ARGUMENT IS WITHDRAWN. Wallach's rule concerns molecular")
    print("  crystals of small molecules. A collagen fibril is a hydrated,")
    print("  quasi-hexagonal array of ~300 nm rods holding a large water")
    print("  fraction — not a molecular crystal. Whether the rule transfers")
    print("  across that boundary was never checked.")
    print()
    print("  It is the same class of error this project had just diagnosed when")
    print("  it refused to transfer the MAX1 result to collagen: carrying a")
    print("  result across a structural boundary without justification. The")
    print("  rigour was also asymmetric — the level-1 demotion rests on a number")
    print("  and a Boltzmann calculation (7.87 kcal/mol, underflowing after ~30")
    print("  residues); the re-promotion rested on an analogy.")
    print()
    print("  >> `racemic_enhanced` is UNRESOLVED for level 2. Neither excluded")
    print("     nor favoured. The same is true of the other three.")

    rule("5  What this survey is actually for")
    print("  Not a prior. A list of leads, and one concrete design consequence:")
    print()
    print(f"  `majority_rules` has the best quantitative backing of the four")
    print(f"  (exact 1D Ising; MMP {MAJORITY_RULES_EVIDENCE['mismatch_penalty_kJ_mol']} "
          f"and HRP {MAJORITY_RULES_EVIDENCE['helix_reversal_penalty_kJ_mol']} kJ/mol)")
    print("  and contributes ZERO mechanical data points — every study of it")
    print("  reads out circular dichroism. That is independent of anything above")
    print("  and is why PROTOCOL.md adds CD to the same samples: a rheology-only")
    print("  experiment would be blind to the one mode with a real theory behind it.")


if __name__ == "__main__":
    main()
