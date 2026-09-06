#!/usr/bin/env python3
"""Derive the assembly competence curve instead of assuming it, and see what breaks.

    ./.venv/bin/python scripts/run_kinetics.py

The assembly layer rested on an ASSUMED symmetric quadratic, 1 - s*4*chi*(1-chi),
chosen only because it was the simplest form with the right endpoints.
Predictions 6-8 all depend on it. This script replaces it with a curve derived
from growth-poisoning kinetics and re-tests them.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np

from mirror_tenocyte.parameters import get_preset, get_default_environment
from mirror_tenocyte.assembly import assembly_competence
from mirror_tenocyte.analysis import (
    detect_rebound, decoupled_scan, design_sensitivity,
    competence_halving_point, kinetics_sensitivity,
)

MODES = ["frustrated", "frustrated_kinetic"]
NAMES = {"frustrated": "assumed quadratic",
         "frustrated_kinetic": "derived from capping kinetics"}


def rule(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def env(mode, **kw):
    e = get_default_environment()
    e["assembly_mode"] = mode
    e.update(kw)
    return e


def main():
    p = get_preset("clutch_calibrated")

    # ---- 1. the two competence curves -----------------------------------
    rule("1  Assumed vs derived competence curve")
    print("  Capping kinetics: an L fibril grows at rate ~(1-chi) and is capped")
    print("  at rate ~chi*p_cap, so mean length before termination is")
    print("  (1-chi)/(chi*p_cap). Load bearing needs length > lambda_c.")
    print()
    print(f"  {'chi':>7} {'assumed':>10} {'derived':>10}")
    for c in [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50]:
        a = assembly_competence({"assembly_mode": "frustrated", "chi_struct": c}, p)
        b = assembly_competence({"assembly_mode": "frustrated_kinetic",
                                 "chi_struct": c}, p)
        print(f"  {c:7.2f} {a:10.4f} {b:10.4f}")
    h = competence_halving_point(p)
    print(f"\n  competence halves at chi = {h:.4f}  ({100 * h:.1f}% D)  [derived]")
    print("  competence halves at chi = 0.1667  (16.7% D)  [assumed quadratic]")
    print("\n  Capping is FAR more efficient than the quadratic guess implied.")
    print("  This is the classic chiral-poisoning behaviour of crystal growth:")
    print("  trace amounts of the opposing enantiomer arrest growth entirely.")

    # ---- 2. what happens to the predictions -----------------------------
    rule("2  Re-testing predictions 6-8 under the derived curve")
    res = {}
    for mode in MODES:
        res[mode] = (detect_rebound(p=p, env=env(mode), n=401),
                     decoupled_scan(p=p, env=env(mode), n=41),
                     design_sensitivity(p=p, env=env(mode)))
    a, b = res["frustrated"], res["frustrated_kinetic"]
    print(f"  {'prediction':>28} {'assumed':>10} {'derived':>10}")
    print(f"  {'6  rebound size':>28} {a[0]['rise']:10.4f} {b[0]['rise']:10.4f}"
          "   <-- COLLAPSES")
    print(f"  {'7  decoupled dip':>28} {a[1]['dip']:10.4f} {b[1]['dip']:10.4f}"
          "   <-- STRENGTHENS")
    print(f"  {'8  best gel stiffness':>28} {a[2]['best_stiffness']:10.0f} "
          f"{b[2]['best_stiffness']:10.0f}   <-- unchanged")
    print("\n  Prediction 6 nearly vanishes. With efficient capping, competence")
    print("  collapses at low chi and only recovers as chi -> 1, by which point")
    print("  the ligand (~1-chi) is gone. There is no window where mechanics")
    print("  have recovered but adhesion still exists — so no rebound.")

    # ---- 3. is that conclusion robust? -----------------------------------
    rule("3  Sensitivity to the capping kinetics (which are NOT measured here)")
    s = kinetics_sensitivity(p=p)
    print(f"  {'lambda_c':>9} {'p_cap':>6} {'6 rebound':>10} {'7 dip':>8} "
          f"{'halve@chi':>10}")
    for r in s["rows"]:
        hc = "n/a" if r["halving_chi"] is None else f"{r['halving_chi']:.4f}"
        print(f"  {r['lambda_c']:9.0f} {r['cap_efficiency']:6.1f} "
              f"{r['rebound']:10.4f} {r['dip']:8.4f} {hc:>10}")
    print(f"\n  rebound spans {s['rebound_range'][0]:.4f} to {s['rebound_range'][1]:.4f}"
          "  <-- SWINGS WILDLY")
    print(f"  dip     spans {s['dip_range'][0]:.4f} to {s['dip_range'][1]:.4f}"
          "  <-- mostly saturated high")

    rule("4  CONCLUSION — which experiment to actually run")
    print("  Prediction 6 (the coupled dose-response 'rebound') is NOT a reliable")
    print("  test. Its magnitude swings from ~0 to ~0.53 across plausible values")
    print("  of two parameters nobody has measured. It was an artefact of the")
    print("  assumed quadratic as much as a consequence of the physics.")
    print()
    print("  Prediction 7 (clamped adhesion, constant L-RGD on a D/L gel) is")
    print("  robust: it stays near-maximal across most of that same range.")
    print("  It should be the primary cell experiment.")
    print()
    print("  >> BUT DO THIS FIRST, AND IT NEEDS NO CELLS AT ALL:")
    print("     rheology on the D/L gel series alone. That single measurement")
    print("       * separates self-sorting (flat modulus) from frustration (dip)")
    print("       * measures the halving composition, which pins lambda_c/p_cap")
    print("       * therefore predicts whether a rebound should exist")
    print("     It is cheaper than any cell experiment and constrains the model")
    print("     more than either of them. Everything downstream is conditional")
    print("     on it.")
    print()
    print("  Note also where the action is: under the derived curve the")
    print("  transition sits at a few percent D, not near 50:50. A composition")
    print("  series sampled at 0, 25, 50, 75, 100% would MISS IT ENTIRELY.")

    _figure(p, os.path.join(ROOT, "results"))


def _figure(p, results_dir):
    os.makedirs(results_dir, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chi = np.linspace(0.0, 1.0, 601)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))

    ax = axes[0]
    for mode, c in zip(MODES, ["tab:blue", "tab:green"]):
        y = [assembly_competence({"assembly_mode": mode, "chi_struct": float(x)}, p)
             for x in chi]
        ax.plot(chi, y, lw=2.2, color=c, label=NAMES[mode])
    ax.set_xlabel("structural chirality  chi")
    ax.set_ylabel("assembly competence")
    ax.set_title("(a) the curve that was assumed")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    for mode, c in zip(MODES, ["tab:blue", "tab:green"]):
        r = detect_rebound(p=p, env=env(mode), n=401)
        ax.plot(r["chi"], r["phen"], lw=2.2, color=c,
                label=f"{NAMES[mode]}  (rebound {r['rise']:+.3f})")
    ax.set_xlabel("environmental chirality  chi")
    ax.set_ylabel("phenotype index")
    ax.set_title("(b) prediction 6 collapses")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[2]
    s = kinetics_sensitivity(p=p)
    for cap, mk in [(1.0, "-o"), (0.3, "--s")]:
        rows = [r for r in s["rows"] if r["cap_efficiency"] == cap]
        ax.semilogx([r["lambda_c"] for r in rows], [r["rebound"] for r in rows],
                    mk, ms=4, lw=2, color="tab:red",
                    label=f"6 rebound (p_cap={cap})")
        ax.semilogx([r["lambda_c"] for r in rows], [r["dip"] for r in rows],
                    mk, ms=4, lw=2, color="tab:green",
                    label=f"7 dip (p_cap={cap})")
    ax.set_xlabel("lambda_c  (monomers needed to bear load)")
    ax.set_ylabel("effect size")
    ax.set_title("(c) 6 is fragile, 7 is robust")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)

    fig.suptitle("Deriving the assembly curve kills prediction 6 and promotes 7",
                 y=1.03, fontsize=12)
    fig.tight_layout()
    out = os.path.join(results_dir, "kinetics.png")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Figure written to {out}")


if __name__ == "__main__":
    main()
