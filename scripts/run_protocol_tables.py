#!/usr/bin/env python3
"""Regenerate every number quoted in PROTOCOL.md.

    ./.venv/bin/python scripts/run_protocol_tables.py

The protocol is the project's actual deliverable, so its tables must be
reproducible rather than transcribed. This prints them in the same order they
appear there.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import itertools

import numpy as np

from mirror_tenocyte.parameters import get_preset
from mirror_tenocyte.assembly import assembly_competence

#: CMP triple helix ~10 nm vs collagen ~300 nm, so ~30x more helices are needed
#: to span the same physical fibril length. lambda_c ~ 20 (collagen) -> ~600.
LAMBDA_CMP = 600.0

MODES = ["self_sorting", "frustrated_kinetic", "majority_rules",
         "racemic_enhanced"]
NICE = {"self_sorting": "self-sorting", "frustrated_kinetic": "capping",
        "majority_rules": "majority rules", "racemic_enhanced": "enhancement"}


def competence(mode, chi, p, lam=None):
    """Normalised modulus G'(chi)/G'(0), as the experiment reports it.

    Each hypothesis is divided by its OWN value at chi = 0, because that is what
    normalising to the pure-L gel does experimentally. It matters for
    majority-rules, whose net helicity at chi = 0 is 0.9916 rather than exactly
    1 (a finite-temperature Ising chain always carries a few reversals).
    """
    q = dict(p)
    if lam is not None:
        q["lambda_c"] = lam
    ref = assembly_competence({"assembly_mode": mode, "chi_struct": 0.0}, q)
    val = assembly_competence({"assembly_mode": mode, "chi_struct": chi}, q)
    return val / ref


def main():
    p = get_preset("clutch_calibrated")

    print("=" * 68)
    print("PROTOCOL section 3 — why the low arm must be log-spaced")
    print("=" * 68)
    print("  Predicted G'(chi)/G'(0) under CAPPING at three critical lengths")
    print(f"  {'chi':>8} " + " ".join(f"{'lc=' + str(l):>12}"
                                      for l in [20, 100, 600]))
    for c in [0.001, 0.005, 0.01, 0.05]:
        vals = [competence("frustrated_kinetic", c, p, l) for l in [20, 100, 600]]
        print(f"  {100*c:7.1f}% " + " ".join(f"{v:12.3f}" for v in vals))
    print("\n  lambda_c = 600 is the CMP-appropriate value: a CMP triple helix is")
    print("  ~10 nm against ~300 nm for collagen, so ~30x more helices span the")
    print("  same physical length. The collapse then sits BELOW 1% D.")

    print("\n" + "=" * 68)
    print("PROTOCOL section 5 — classification table")
    print("=" * 68)
    print(f"  {'hypothesis':>16} {'at chi=1%':>12} {'at chi=50%':>12}   pattern")
    pattern = {"self_sorting": "flat throughout",
               "frustrated_kinetic": "collapses, stays dead",
               "majority_rules": "flat, then dips at racemic",
               "racemic_enhanced": "rises, peaks at racemic"}
    sig = {}
    for m in MODES:
        lam = LAMBDA_CMP if m == "frustrated_kinetic" else None
        a = competence(m, 0.01, p, lam)
        b = competence(m, 0.50, p, lam)
        sig[m] = np.array([a, b])
        print(f"  {NICE[m]:>16} {a:12.2f} {b:12.2f}   {pattern[m]}")

    worst = min((float(np.linalg.norm(sig[a] - sig[b])), NICE[a], NICE[b])
                for a, b in itertools.combinations(MODES, 2))
    print(f"\n  worst-case pairwise separation: {worst[0]:.3f} "
          f"({worst[1]} vs {worst[2]})")
    print("  Instrument precision is not the limit; gel-to-gel reproducibility is.")

    print("\n" + "=" * 68)
    print("PROTOCOL section 6 — symmetry control")
    print("=" * 68)
    print("  Every hypothesis predicts G'(chi) = G'(1-chi).")
    print(f"  {'hypothesis':>16} {'G(10%)':>10} {'G(90%)':>10} {'|diff|':>10}")
    ok = True
    for m in MODES:
        lam = LAMBDA_CMP if m == "frustrated_kinetic" else None
        a = competence(m, 0.10, p, lam)
        b = competence(m, 0.90, p, lam)
        ok &= abs(a - b) < 1e-9
        print(f"  {NICE[m]:>16} {a:10.4f} {b:10.4f} {abs(a-b):10.2e}")
    print(f"\n  symmetry holds for all hypotheses: {ok}")
    print("  So any measured asymmetry is a problem with the samples, or a")
    print("  finding in its own right — never an expected model behaviour.")

    print("\n" + "=" * 68)
    print("Full recommended composition series")
    print("=" * 68)
    series = [0, 0.001, 0.003, 0.01, 0.03, 0.10, 0.20, 0.30, 0.40, 0.50,
              0.90, 0.99, 1.00]
    print(f"  {'chi':>8} " + " ".join(f"{NICE[m][:13]:>14}" for m in MODES))
    for c in series:
        vals = []
        for m in MODES:
            lam = LAMBDA_CMP if m == "frustrated_kinetic" else None
            vals.append(competence(m, c, p, lam))
        print(f"  {100*c:7.1f}% " + " ".join(f"{v:14.3f}" for v in vals))
    print(f"\n  {len(series)} compositions x 3 replicates = {3*len(series)} gels")

    print("\n" + "=" * 68)
    print("PROTOCOL section 3 — material budget and pilot staging")
    print("=" * 68)
    n_rep, vol_uL, wt_pct = 3, 200, 1.0
    mg_per_gel = vol_uL / 1000 * 10 * wt_pct
    d_mg = sum(series) * n_rep * mg_per_gel
    l_mg = sum(1 - c for c in series) * n_rep * mg_per_gel
    print(f"  {vol_uL} uL per gel at {wt_pct} wt% -> {mg_per_gel:.1f} mg peptide per gel")
    print(f"  {'':22} {'gels only':>12} {'3x overage':>12}")
    print(f"  {'D-CMP':>22} {d_mg:11.0f} mg {3*d_mg:11.0f} mg")
    print(f"  {'L-CMP':>22} {l_mg:11.0f} mg {3*l_mg:11.0f} mg")
    low = sum(c for c in series if c <= 0.10) * n_rep * mg_per_gel
    print(f"\n  low arm (0-10% D, {sum(1 for c in series if c <= 0.10)} compositions,"
          f" {n_rep*sum(1 for c in series if c <= 0.10)} gels) uses only "
          f"{low:.2f} mg of D-CMP")
    print("  => run it first: it separates capping from the other three for a")
    print("     fraction of a milligram of the expensive enantiomer.")


if __name__ == "__main__":
    main()
