#!/usr/bin/env python3
"""Four competing assembly hypotheses, and the protocol that separates them.

    ./.venv/bin/python scripts/run_modes.py

The model started with two hypotheses (self-sorting, frustration). Literature
checks added two more: racemic coassembly that is STIFFER than either pure form
(inverting the predicted sign), and majority-rules chiral copolymerisation, in
which the minority enantiomer is ABSORBED rather than excluded.

The fourth one breaks the previous design advice. Frustration and majority-rules
are both collapsed at racemic composition and differ there by only ~0.008, so no
single racemic-vs-pure contrast can separate them at any stiffness. They differ
in the SHAPE of the composition curve, so the protocol needs two compositions.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np

from mirror_tenocyte.parameters import get_preset, get_default_environment
from mirror_tenocyte.assembly import effective_stiffness
from mirror_tenocyte.analysis import (
    decoupled_scan, discriminate_assembly_modes, ASSEMBLY_MODES,
)

LABEL = {
    "self_sorting": "self-sorting (separate perfect fibrils)",
    "frustrated_kinetic": "frustration (capping arrests growth)",
    "majority_rules": "majority rules (minority absorbed)",
    "racemic_enhanced": "racemic coassembly (heterochiral, stiffer)",
}
EVIDENCE = {
    "self_sorting": "documented in peptide amphiphiles; Abeta forms homochiral sheets",
    "frustrated_kinetic": "inferred from chiral crystal-growth poisoning",
    "majority_rules": "MEASURED: MMP 0.94 vs HRP 7.8 kJ/mol; exact 1D Ising theory",
    "racemic_enhanced": "MEASURED: MAX1/DMAX1 800 Pa vs 200 Pa (4x)",
}


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

    rule("1  The hypothesis space — four-sided")
    for m in ASSEMBLY_MODES:
        print(f"  {m:20s} {LABEL[m]}")
        print(f"  {'':20s}   evidence: {EVIDENCE[m]}")
    print("\n  All four are empirically realised in DIFFERENT systems. Collagen is")
    print("  a PPII triple helix, structurally unlike a beta-hairpin, so neither")
    print("  the MAX1 result nor the peptide-amphiphile result transfers to it")
    print("  automatically. Which mode collagen follows is genuinely unknown.")

    rule("1b  Does the MAX1 result actually transfer to collagen? No.")
    from mirror_tenocyte.assembly import heterochiral_helix_probability
    print("  Collagen assembles hierarchically: chain -> triple helix ->")
    print("  microfibril -> fibril. Chirality mismatch can act at more than one")
    print("  level, and the levels are NOT equivalent.")
    print()
    print("  LEVEL 1 (chain -> triple helix). Three chains interdigitate with Gly")
    print("  at the core. MD puts a single D-Ala substitution at ~7.87 kcal/mol —")
    print("  comparable to the Gly->Ala mutations causing osteogenesis imperfecta")
    print("  — kinking the helix and breaking its register.")
    print(f"\n  {'mismatched residues':>22} {'P(heterochiral helix)':>22}")
    for n in [1, 3, 10, 30]:
        print(f"  {n:22} {heterochiral_helix_probability(p, n):22.3e}")
    print("\n  A triple-helical domain has ~1000 residues, so this underflows long")
    print("  before a full chain. Heterochiral triple helices DO NOT FORM:")
    print("  level 1 sorts perfectly, with no free parameter.")
    print()
    print("  >> WHY THAT MATTERS: the MAX1/DMAX1 stiffening works through a")
    print("     'rippled' beta-sheet of ALTERNATING L and D strands — a level-1")
    print("     mixing mechanism. Collagen forbids level-1 mixing, so the")
    print("     strongest empirical support for racemic enhancement DOES NOT")
    print("     TRANSFER. Enhancement via level-2 helix co-packing is not")
    print("     excluded, but it is now unsupported rather than well-evidenced.")
    print()
    print("  LEVEL 2 (triple helix -> fibril) is where the open question lives.")
    print("  Note the unit change: the capping 'monomer' is then a whole triple")
    print("  helix (~300 nm rod), so lambda_c counts helices, not residues.")

    rule("2  Effective stiffness vs composition (base 30 kPa)")
    cols = [0.0, 0.25, 0.5, 0.75, 1.0]
    print(f"  {'mode':>20} " + " ".join(f"{c:>8.2f}" for c in cols))
    for m in ASSEMBLY_MODES:
        vals = [effective_stiffness(env(m, chi_struct=c), p) for c in cols]
        print(f"  {m:>20} " + " ".join(f"{v:8.2f}" for v in vals))
    print("\n  Note all four agree at the endpoints — a pure mirror matrix is")
    print("  mechanically identical to a pure natural one in every hypothesis.")
    print("  They differ ONLY at intermediate composition. Comparing pure-L with")
    print("  pure-D therefore cannot distinguish them at all.")

    rule("3  Signs at racemic — and why sign alone is not enough")
    print(f"  {'mode':>20} {'pure L':>8} {'racemic':>8} {'pure D':>8} "
          f"{'signed':>9}  shape")
    for m in ASSEMBLY_MODES:
        d = decoupled_scan(p=p, env=env(m), n=41)
        signed = d["phen_racemic"] - d["phen_pure_L"]
        shape = ("FLAT" if abs(signed) < 1e-3
                 else ("DIP" if signed < 0 else "PEAK"))
        print(f"  {m:>20} {d['phen_pure_L']:8.4f} {d['phen_racemic']:8.4f} "
              f"{d['phen_pure_D']:8.4f} {signed:9.4f}  {shape}")
    print("\n  Sign is more robust than magnitude — but note that frustration and")
    print("  majority-rules give the SAME sign and nearly the same value here.")
    print("  A racemic-vs-pure contrast cannot tell them apart.")

    rule("4  No single composition separates all four")
    disc = discriminate_assembly_modes(p=p)
    print(f"  {'E (kPa)':>8} " + " ".join(f"{m[:14]:>15}" for m in ASSEMBLY_MODES)
          + "   min gap")
    for r in disc["rows"]:
        print(f"  {r['stiffness']:8.1f} "
              + " ".join(f"{r['effects'][m]:15.4f}" for m in ASSEMBLY_MODES)
              + f"   {r['min_separation']:.4f}")
    print(f"\n  Best single-point separation is only {disc['best_separation']:.4f} "
          f"(at {disc['best_stiffness']:.0f} kPa),")
    print("  limited by the frustration vs majority-rules pair.")
    print("  An earlier version of this script recommended 2 kPa on the strength")
    print("  of a THREE-mode analysis, where the gap was 0.290. Adding the fourth")
    print("  hypothesis cuts that to 0.008 at the same stiffness. The stiffness")
    print("  was never the problem — the single contrast was.")

    rule("5  Two compositions do separate all four")
    from mirror_tenocyte.analysis import two_point_protocol
    tp = two_point_protocol(p=p)
    lo, hi = tp["compositions"]
    print(f"  Measure at chi = {lo} and chi = {hi}, adhesion clamped.")
    print(f"  {'E (kPa)':>8} {'worst-case sep':>15}   limiting pair")
    for r in tp["rows"]:
        a, b = r["limiting_pair"]
        print(f"  {r['stiffness']:8.1f} {r['worst_separation']:15.4f}   {a} vs {b}")
    print(f"\n  BEST: {tp['best_stiffness']:.0f} kPa, worst-case separation "
          f"{tp['best_separation']:.4f}")
    print()
    print(f"  {'mode':>20} {'at chi=' + str(lo):>12} {'at chi=' + str(hi):>12}   pattern")
    patterns = {"self_sorting": "flat, flat", "frustrated_kinetic": "DOWN, down",
                "majority_rules": "flat, DOWN", "racemic_enhanced": "up,   up"}
    for m in ASSEMBLY_MODES:
        a, b = tp["best_signatures"][m]
        print(f"  {m:>20} {a:12.4f} {b:12.4f}   {patterns[m]}")
    print()
    print("  >> RECOMMENDATION: clamped-adhesion on a soft (~2 kPa) gel, at TWO")
    print("     compositions (10% and 50% D). Each hypothesis has a distinct")
    print("     two-point pattern; error bars must beat "
          f"{tp['best_separation']:.2f}.")
    print()
    print("  And still cheaper and prior to all of it: RHEOLOGY on the gel series")
    print("  with no cells. Sample BOTH low D fraction and near-racemic — the")
    print("  four modes put their signal in different places.")

    rule("6  Does any of this depend on my made-up phenotype index?")
    from mirror_tenocyte.analysis import (
        observable_discrimination, observable_ranking_two_point,
    )
    od = observable_discrimination(p=p, stiffness=2.0)
    print("  Every number above is in `phenotype_index`, a composite whose")
    print("  weights I chose by hand. Sign-consistency across all six states:")
    print(f"    weight-independent: {od['weight_independent']}")
    print()
    print("  That sounds like robustness, but the Jacobian of this model has zero")
    print("  negative off-diagonal couplings — it is a MONOTONE chain, so sign")
    print("  agreement is forced by the architecture, not evidence about biology.")

    rule("7  What to measure — and a trap")
    two = observable_ranking_two_point(p=p)
    print(f"  Under the TWO-POINT protocol (chi = {two['compositions'][0]} and "
          f"{two['compositions'][1]} at {two['stiffness']:.0f} kPa):")
    print()
    print(f"  {'worst-case sep':>15}  {'observable':>9}  assay")
    for r in two["rows"]:
        print(f"  {r['worst_separation']:15.4f}  {r['observable']:>9}  {r['assay']}")
    print()
    print(f"  Best: {two['best']}.  Worst: {two['worst']} — avoid it.")
    print()
    print("  >> THE TRAP: under the INADEQUATE single-point protocol the ranking")
    print(f"     INVERTS — it would tell you to measure {od['best']} and that")
    print(f"     {od['worst']} is useless, the opposite of the truth. The limiting")
    print("     pair there is frustration vs majority-rules, which differ most")
    print("     UPSTREAM. Fix the protocol before choosing the assay.")

    _figure(p, disc, os.path.join(ROOT, "results"))


def _figure(p, disc, results_dir):
    os.makedirs(results_dir, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chi = np.linspace(0.0, 1.0, 401)
    colors = {"self_sorting": "tab:blue",
              "frustrated_kinetic": "tab:green",
              "majority_rules": "tab:orange",
              "racemic_enhanced": "tab:red"}

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))

    ax = axes[0]
    for m in ASSEMBLY_MODES:
        y = [effective_stiffness(env(m, chi_struct=float(c)), p) for c in chi]
        ax.semilogy(chi, y, lw=2.2, color=colors[m], label=LABEL[m].split(" (")[0])
    ax.set_xlabel("structural chirality  chi")
    ax.set_ylabel("effective stiffness (kPa, log)")
    ax.set_title("(a) three hypotheses, same endpoints")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)

    ax = axes[1]
    for m in ASSEMBLY_MODES:
        e = env(m)
        e["Estiff"] = 2.0
        d = decoupled_scan(p=p, env=e, n=61)
        signed = d["phen_racemic"] - d["phen_pure_L"]
        ax.plot(d["chi_struct"], d["phen"], lw=2.2, color=colors[m],
                label=f"{LABEL[m].split(' (')[0]}  ({signed:+.3f})")
    ax.set_xlabel("structural chirality  (adhesion clamped)")
    ax.set_ylabel("phenotype index")
    ax.set_title("(b) clamped adhesion at 2 kPa: three signs")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)

    ax = axes[2]
    for m in ASSEMBLY_MODES:
        ax.semilogx([r["stiffness"] for r in disc["rows"]],
                    [r["effects"][m] for r in disc["rows"]],
                    "-o", ms=4, lw=2, color=colors[m],
                    label=LABEL[m].split(" (")[0])
    ax.axvline(disc["best_stiffness"], color="k", ls=":", lw=2,
               label=f"best separation: {disc['best_stiffness']:.0f} kPa")
    ax.axvline(8, color="gray", ls="--", lw=1.5, label="old advice: 8 kPa")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("base gel stiffness (kPa)")
    ax.set_ylabel("signed effect at racemic")
    ax.set_title("(c) 8 kPa would hide the enhancement")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)

    fig.suptitle("A literature check added a third hypothesis that inverts the "
                 "prediction — and moves the optimal gel from 8 to 2 kPa",
                 y=1.03, fontsize=12)
    fig.tight_layout()
    out = os.path.join(results_dir, "modes.png")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Figure written to {out}")


if __name__ == "__main__":
    main()
