#!/usr/bin/env python3
"""Generate the core figures and a summary table for the mirror-tenocyte model.

Run from the repository root (with the venv active):

    ./.venv/bin/python scripts/run_dose_response.py

Outputs three PNGs into results/ and prints a chirality dose-response table plus
the mirror-vs-non-adhesive comparison.
"""

import os
import sys

# make the package importable when run as a plain script
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np

from mirror_tenocyte.simulate import (
    dose_response,
    stiffness_grid,
    hysteresis,
    compare_mirror_vs_nonadhesive,
    plot_dose_response,
    plot_stiffness_grid,
    plot_hysteresis,
    plot_preset_comparison,
)
from mirror_tenocyte.model import STATE_NAMES, classify_phenotype


def main():
    results_dir = os.path.join(ROOT, "results")
    os.makedirs(results_dir, exist_ok=True)

    chi = np.linspace(0.0, 1.0, 41)

    dr = dose_response(chi)
    plot_dose_response(dr, os.path.join(results_dir, "dose_response.png"))

    grid = stiffness_grid(chi)
    plot_stiffness_grid(grid, os.path.join(results_dir, "stiffness_grid.png"))

    hys = hysteresis(chi)
    plot_hysteresis(hys, os.path.join(results_dir, "hysteresis.png"))

    plot_preset_comparison(os.path.join(results_dir, "preset_comparison.png"))

    # ---- text summary -------------------------------------------------
    print("\nChirality dose-response (default 30 kPa, racemic coupling)")
    print(f"{'chi':>5}  {'Cint':>6} {'FAKp':>6} {'Tension':>7} {'YAPn':>6} "
          f"{'Col':>6} {'phenotype':>9}")
    for i in range(0, len(dr["chi"]), 4):
        print(f"{dr['chi'][i]:5.2f}  {dr['Cint'][i]:6.3f} {dr['FAKp'][i]:6.3f} "
              f"{dr['Tension'][i]:7.3f} {dr['YAPn'][i]:6.3f} {dr['Col'][i]:6.3f} "
              f"{dr['phenotype'][i]:9.3f}")

    # locate the steepest drop (candidate critical chirality chi*)
    dphen = np.diff(dr["phenotype"])
    i_star = int(np.argmin(dphen))
    chi_star = 0.5 * (dr["chi"][i_star] + dr["chi"][i_star + 1])
    print(f"\nSteepest phenotype transition near chi* ~ {chi_star:.2f} "
          f"(delta phenotype = {dphen[i_star]:.3f} per {chi[1]-chi[0]:.3f} chi step)")

    comp = compare_mirror_vs_nonadhesive()
    print("\nMirror-D vs. non-adhesive (PEG) control — steady states:")
    header = "  ".join(f"{n:>7}" for n in STATE_NAMES)
    print(f"{'condition':>18}  {header}")
    for cond, st in comp.items():
        vals = "  ".join(f"{v:7.3f}" for v in st)
        print(f"{cond:>18}  {vals}")
    print("\nPhenotype calls:")
    for cond, st in comp.items():
        print(f"  {cond:>18}: {classify_phenotype(st)}")

    print(f"\nFigures written to {results_dir}/")


if __name__ == "__main__":
    main()
