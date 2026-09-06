#!/usr/bin/env python3
"""The assembly (fibrillogenesis) layer: does chirality become a real axis?

    ./.venv/bin/python scripts/run_assembly.py

Compares the three assembly hypotheses, tests whether each breaks the
Ltot*(1-chi) degeneracy, finds the non-monotonic rebound that distinguishes
frustrated assembly, and reports the gel stiffness that maximises the effect.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import copy
import numpy as np

from mirror_tenocyte.parameters import get_default_parameters, get_default_environment
from mirror_tenocyte.assembly import effective_stiffness
from mirror_tenocyte.analysis import (
    product_invariance, detect_rebound, decoupled_scan, design_sensitivity,
)

MODES = ["none", "self_sorting", "frustrated"]


def rule(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def env_with(mode, **kw):
    e = get_default_environment()
    e["assembly_mode"] = mode
    e.update(kw)
    return e


def main():
    p = get_default_parameters()
    results_dir = os.path.join(ROOT, "results")
    os.makedirs(results_dir, exist_ok=True)

    # ---- 1. does the degeneracy break? ---------------------------------
    rule("1  Does the assembly layer break the Ltot*(1-chi) degeneracy?")
    for mode in MODES:
        inv = product_invariance(p=p, env=env_with(mode))
        verdict = ("DEGENERATE — chirality is disguised ligand density"
                   if inv["invariant"] else "DEGENERACY BROKEN — chirality is a real axis")
        print(f"  {mode:13s}: max spread = {inv['max_spread']:.3e}   {verdict}")
    print("\n  Self-sorting reduces EXACTLY to the parsimonious model: if D and L")
    print("  segregate into separate perfect fibrils, every fibril is mechanically")
    print("  normal and only the bindable fraction falls. So the degeneracy is the")
    print("  fingerprint of self-sorting, and its violation is the fingerprint of")
    print("  assembly frustration.")

    # ---- 2. mechanical competence -------------------------------------
    rule("2  Effective stiffness vs structural chirality (frustrated mode)")
    print(f"  {'chi_struct':>11}  {'E_eff (kPa)':>12}")
    for cs in [0.0, 0.25, 0.5, 0.75, 1.0]:
        e = env_with("frustrated", chi_struct=cs)
        print(f"  {cs:11.2f}  {effective_stiffness(e, p):12.3f}")
    print("\n  Symmetric about 0.5 and fully recovered at chi = 1: a PURE mirror")
    print("  matrix is mechanically normal (enantiomers are isoenergetic and")
    print("  elastic moduli are parity-invariant). Only MIXTURES are frustrated.")

    # ---- 3. the coupled dose-response signature ------------------------
    rule("3  Coupled dose-response: the rebound signature")
    for mode in MODES:
        r = detect_rebound(p=p, env=env_with(mode))
        if r["rebound"]:
            print(f"  {mode:13s}: REBOUND  chi {r['chi_min']:.3f} -> {r['chi_max']:.3f}, "
                  f"phenotype {r['phen_min']:.4f} -> {r['phen_max']:.4f} "
                  f"(rise {r['rise']:+.4f})")
        else:
            print(f"  {mode:13s}: monotonic decline, no rebound")
    print("\n  This is the discriminating experiment, and it is the EASY one: just")
    print("  make D/L collagen gels across a composition series and read out")
    print("  nuclear YAP / collagen. A monotonic sigmoid supports recognition-loss")
    print("  only; a dip near racemic followed by partial recovery around")
    print("  chi ~ 0.75-0.8 supports assembly frustration.")

    # ---- 4. the decoupled experiment -----------------------------------
    rule("4  Decoupled experiment: constant L-RGD, varying structural chirality")
    for mode in MODES:
        d = decoupled_scan(p=p, env=env_with(mode))
        shape = "U-shaped" if d["dip"] > 1e-3 else "flat"
        print(f"  {mode:13s}: L={d['phen_pure_L']:.4f}  racemic={d['phen_racemic']:.4f}  "
              f"D={d['phen_pure_D']:.4f}   dip={d['dip']:.4f}  [{shape}, "
              f"symmetric={d['symmetric']}]")
    print("\n  Clamping adhesion removes the confound entirely: under recognition")
    print("  loss this curve must be FLAT, under frustration it is U-shaped and")
    print("  symmetric, with pure-D indistinguishable from pure-L.")

    # ---- 5. experimental design ----------------------------------------
    rule("5  Which gel stiffness makes the effect detectable?")
    ds = design_sensitivity(p=p)
    print(f"  {'E_base':>8} {'E at racemic':>13} {'phen(pure)':>11} "
          f"{'phen(racemic)':>14} {'dip':>8}")
    for r in ds["rows"]:
        print(f"  {r['stiffness']:8.0f} {r['stiffness']*0.1:13.2f} "
              f"{r['phen_pure']:11.4f} {r['phen_racemic']:14.4f} {r['dip']:8.4f}")
    print(f"\n  Deepest contrast at E_base = {ds['best_stiffness']:.0f} kPa "
          f"(dip = {ds['best_dip']:.3f}).")
    print("  PRACTICAL POINT: at tendon-like 30 kPa the dip is only ~0.08, because")
    print("  even a frustrated gel stays above the mechanosensing threshold. Around")
    print("  8 kPa the racemic midpoint falls below threshold and the effect is")
    print("  ~8x larger. Choosing the wrong gel stiffness could hide the phenomenon.")

    _figures(p, results_dir)

    print("\n" + "=" * 72)
    print("Bottom line: chirality becomes an independent variable only if mixed")
    print("D/L matrices fail to assemble. That hypothesis is directly testable,")
    print("and it predicts a rebound the recognition-loss model cannot produce.")
    print("=" * 72)


def _figures(p, results_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chi = np.linspace(0.0, 1.0, 201)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))

    # (a) assembly competence / effective stiffness
    ax = axes[0]
    for mode in ["self_sorting", "frustrated"]:
        ys = []
        for c in chi:
            e = env_with(mode, chi_struct=float(c))
            ys.append(effective_stiffness(e, p))
        ax.plot(chi, ys, lw=2.2, label=mode)
    ax.set_xlabel("structural chirality  chi_struct")
    ax.set_ylabel("effective stiffness (kPa)")
    ax.set_title("(a) mechanical competence of the matrix")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    # (b) coupled dose-response — the rebound
    ax = axes[1]
    for mode in MODES:
        r = detect_rebound(p=p, env=env_with(mode))
        lbl = mode + (f"  (rebound {r['rise']:+.3f})" if r["rebound"] else "  (monotonic)")
        ax.plot(r["chi"], r["phen"], lw=2.2, label=lbl)
    ax.set_xlabel("environmental chirality  chi")
    ax.set_ylabel("phenotype index")
    ax.set_title("(b) coupled dose-response")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)

    # (c) decoupled scan — adhesion clamped
    ax = axes[2]
    for mode in MODES:
        d = decoupled_scan(p=p, env=env_with(mode), n=81)
        ax.plot(d["chi_struct"], d["phen"], lw=2.2,
                label=f"{mode}  (dip {d['dip']:.3f})")
    ax.set_xlabel("structural chirality  chi_struct   (adhesion clamped)")
    ax.set_ylabel("phenotype index")
    ax.set_title("(c) decoupled: constant L-RGD")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)

    fig.suptitle("Does mixed-chirality matrix fail to assemble? Three testable signatures",
                 y=1.03, fontsize=12)
    fig.tight_layout()
    out = os.path.join(results_dir, "assembly.png")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Figure written to {out}")


if __name__ == "__main__":
    main()
