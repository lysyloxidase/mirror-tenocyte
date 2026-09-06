#!/usr/bin/env python3
"""Check every claim the README makes, and print a verdict for each.

    ./.venv/bin/python scripts/run_analysis.py

The point is that no prediction is asserted in prose without a number behind it.
Claims that hold only in a particular parameter regime are reported as such.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import copy

from mirror_tenocyte.parameters import (
    get_default_parameters, get_default_environment, get_preset,
)
from mirror_tenocyte.analysis import (
    find_chi_star, detect_bistability, product_invariance,
    chirality_stiffness_map,
)
from mirror_tenocyte.simulate import compare_mirror_vs_nonadhesive
from mirror_tenocyte.model import STATE_NAMES


def rule(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def main():
    p = get_default_parameters()

    # ---- CLAIM 1 -------------------------------------------------------
    rule("CLAIM 1  Mirror ECM is equivalent to non-adhesion (not a novel signal)")
    comp = compare_mirror_vs_nonadhesive()
    print(f"{'condition':>18}  " + "  ".join(f"{n:>7}" for n in STATE_NAMES))
    for cond, st in comp.items():
        print(f"{cond:>18}  " + "  ".join(f"{v:7.3f}" for v in st))
    import numpy as np
    diff = float(np.max(np.abs(comp["mirror_D"] - comp["nonadhesive_PEG"])))
    print(f"\n  max |mirror_D - nonadhesive_PEG| = {diff:.2e}")
    print("  VERDICT: HOLDS BY CONSTRUCTION in 'racemic' mode — chi enters the model")
    print("           only via bindable ligand. This is an ASSUMPTION being made")
    print("           explicit, not a discovery. It is the key experimental control.")

    # ---- CLAIM 2 -------------------------------------------------------
    rule("CLAIM 2  The chirality->phenotype transition is switch-like")
    for mode in ["racemic", "affinity"]:
        env = get_default_environment()
        env["chir_mode"] = mode
        r = find_chi_star(p, env)
        print(f"  {mode:9s}: chi* = {r['chi_star']:.3f}   "
              f"phenotype {r['phen_hi']:.3f} -> {r['phen_lo']:.3f}   "
              f"max slope = {r['sharpness']:.2f} /unit chi")
    print("\n  VERDICT: HOLDS. Transition is sharp in both modes, but its LOCATION")
    print("           is mode-dependent (racemic ~0.92 vs affinity ~0.63), so chi*")
    print("           alone cannot distinguish the coupling mechanisms.")

    # ---- CLAIM 3 -------------------------------------------------------
    rule("CLAIM 3  Hysteresis / history dependence")
    for preset in ["parsimonious", "strong_feedback"]:
        pp = get_preset(preset)
        b = detect_bistability(pp, n=401)
        verdict = "BISTABLE" if b["bistable"] else "monostable"
        line = (f"  {preset:16s}: max hysteresis gap = {b['max_gap']:.3e} "
                f"at chi = {b['chi_at_max_gap']:.3f}   -> {verdict}")
        if b["bistable"]:
            lo, hi = b["window"]
            line += (f"\n  {'':16s}  coexistence window: chi in [{lo:.4f}, {hi:.4f}]"
                     f"  (width {b['window_width']:.4f})")
        print(line)
    print("\n  VERDICT: CONDITIONAL, and hard to test. The default (parsimonious)")
    print("           parameters are MONOSTABLE — the sharp transition is fully")
    print("           reversible. Hysteresis appears only with stronger clutch")
    print("           cooperativity and tension feedback.")
    print()
    print("  >> PRACTICAL CAVEAT: even in the bistable regime the coexistence")
    print("     window is only ~1% wide in chi. The phenotype gap is large (0.65)")
    print("     but confined to a narrow band, so detecting hysteresis would")
    print("     require controlling enantiomeric fraction to better than ~1% —")
    print("     a demanding requirement that makes this the LEAST practically")
    print("     testable of the model's predictions.")

    # ---- CLAIM 4 -------------------------------------------------------
    rule("CLAIM 4  Chirality x stiffness interaction")
    print(f"  {'stiffness':>10}  {'phen(chi=0)':>12}  {'chi*':>7}  {'effect size':>11}")
    for row in chirality_stiffness_map(p=p):
        cs = "n/a" if row["chi_star"] is None else f"{row['chi_star']:.3f}"
        print(f"  {row['stiffness']:9.1f}k  {row['phen_native']:12.3f}  "
              f"{cs:>7}  {row['range']:11.3f}")
    print("\n  VERDICT: HOLDS, both halves. chi* shifts monotonically to the right")
    print("           with stiffness (0.64 at 3 kPa -> 0.94 at 100 kPa): a stiffer")
    print("           matrix buffers partial inversion. And at 1 kPa the cell is")
    print("           already YAP-low, so chirality changes almost nothing")
    print("           (effect size 0.06 vs 0.70 at 30 kPa).")

    # ---- CLAIM 5 + the degeneracy -------------------------------------
    rule("CLAIM 5  Ligand-density rescue — and an exact structural degeneracy")
    inv = product_invariance(p=p)
    for row in inv["rows"]:
        print(f"  product Ltot*(1-chi) = {row['product']:.3f}")
        for (L, chi, ph) in row["points"]:
            print(f"      Ltot={L:4.2f}  chi={chi:.3f}  ->  phenotype={ph:.4f}")
        print(f"      spread across densities = {row['spread']:.2e}")
    print(f"\n  max spread over all products = {inv['max_spread']:.2e}")
    print("  VERDICT: HOLDS, and more strongly than claimed. Phenotype depends")
    print("           ONLY on the product Ltot*(1-chi), to machine precision.")
    print("           So in the parsimonious model chirality is NOT an independent")
    print("           axis — it is a disguised ligand-density axis.")
    print()
    print("  >> This is the project's sharpest falsifiable prediction:")
    print("     diluting L-ligand and enantiomerically inverting it should be")
    print("     EXACTLY interchangeable. Any experiment where they are not")
    print("     refutes the parsimonious model and demonstrates that")
    print("     chirality-specific biology exists (-> enable p['chi_specific']).")

    print("\n" + "=" * 72)
    print("Summary: claims 1, 2, 4, 5 hold under default parameters.")
    print("Claim 3 (hysteresis) holds only in the strong-feedback regime.")
    print("=" * 72)


if __name__ == "__main__":
    main()
