#!/usr/bin/env python3
"""Time-course trajectories of every node at a few chirality levels.

    ./.venv/bin/python scripts/run_timecourse.py

Writes results/timecourse.png showing how the network relaxes to its steady
state under increasingly mirror-image environments.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import copy
import numpy as np

from mirror_tenocyte.parameters import get_default_parameters, get_default_environment
from mirror_tenocyte.model import run_timecourse, STATE_NAMES


def main():
    results_dir = os.path.join(ROOT, "results")
    os.makedirs(results_dir, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = get_default_parameters()
    chi_levels = [0.0, 0.4, 0.6, 0.8, 1.0]

    fig, axes = plt.subplots(1, len(chi_levels), figsize=(16, 3.4), sharey=True)
    for ax, chi in zip(axes, chi_levels):
        env = copy.deepcopy(get_default_environment())
        env["chi"] = chi
        t, Y = run_timecourse(p, env, t_max=30.0)
        for i, name in enumerate(STATE_NAMES):
            ax.plot(t, Y[i], lw=1.8, label=name)
        ax.set_title(f"chi = {chi:g}")
        ax.set_xlabel("time (a.u.)")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("activity [0, 1]")
    axes[-1].legend(fontsize=7, loc="center right")
    fig.suptitle("Node trajectories vs. environmental chirality", y=1.02)
    fig.tight_layout()
    out = os.path.join(results_dir, "timecourse.png")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
