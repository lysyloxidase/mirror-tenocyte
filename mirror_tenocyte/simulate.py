"""Experiments on the model: chirality dose-response, stiffness interaction,
hysteresis, and time courses. Pure-numeric functions here; plotting helpers at
the bottom use a non-interactive backend so scripts run headless.
"""

from __future__ import annotations

import copy

import numpy as np

from .parameters import get_default_parameters, get_default_environment
from .model import (
    run_to_steady_state,
    run_timecourse,
    phenotype_index,
    STATE_NAMES,
)


def dose_response(chi_values=None, p=None, env=None):
    """Steady state as a function of environmental chirality chi.

    Returns dict with 'chi', each state variable's steady-state array, and
    'phenotype'. Each chi is solved from the same naive initial condition.
    """
    if chi_values is None:
        chi_values = np.linspace(0.0, 1.0, 41)
    p = p or get_default_parameters()
    env = copy.deepcopy(env or get_default_environment())

    states = []
    for chi in chi_values:
        env["chi"] = float(chi)
        states.append(run_to_steady_state(p, env))
    states = np.array(states)  # (n_chi, n_state)

    out = {"chi": np.asarray(chi_values, dtype=float)}
    for i, name in enumerate(STATE_NAMES):
        out[name] = states[:, i]
    out["phenotype"] = np.array([phenotype_index(s) for s in states])
    return out


def stiffness_grid(chi_values=None, stiffness_values=(1.0, 10.0, 30.0, 100.0),
                   p=None, env=None):
    """Chirality dose-response at several substrate stiffnesses (kPa)."""
    if chi_values is None:
        chi_values = np.linspace(0.0, 1.0, 41)
    p = p or get_default_parameters()
    base_env = copy.deepcopy(env or get_default_environment())

    curves = {}
    for e in stiffness_values:
        env_e = copy.deepcopy(base_env)
        env_e["Estiff"] = float(e)
        curves[e] = dose_response(chi_values, p=p, env=env_e)
    return curves


def hysteresis(chi_values=None, p=None, env=None):
    """Forward (L->D) then backward (D->L) sweep, each seeded from the previous
    steady state. Divergence between branches indicates bistability/hysteresis.
    """
    if chi_values is None:
        chi_values = np.linspace(0.0, 1.0, 41)
    p = p or get_default_parameters()
    env = copy.deepcopy(env or get_default_environment())

    up = np.asarray(chi_values, dtype=float)
    down = up[::-1]

    def sweep(seq):
        phen, y = [], None
        for chi in seq:
            env["chi"] = float(chi)
            y = run_to_steady_state(p, env, y0=y)
            phen.append(phenotype_index(y))
        return np.array(phen)

    up_phen = sweep(up)
    down_phen = sweep(down)
    return {"chi_up": up, "phen_up": up_phen,
            "chi_down": down, "phen_down": down_phen}


def compare_mirror_vs_nonadhesive(p=None, env=None):
    """Contrast a fully mirror matrix (chi=1) with a genuinely non-adhesive
    substrate (Ltot=0). In 'racemic' mode they coincide *by construction* —
    which is the model's core prediction and the key experimental control.
    """
    p = p or get_default_parameters()
    env = copy.deepcopy(env or get_default_environment())

    mirror = copy.deepcopy(env); mirror["chi"] = 1.0
    peg = copy.deepcopy(env); peg["chi"] = 0.0; peg["Ltot"] = 0.0
    native = copy.deepcopy(env); native["chi"] = 0.0

    return {
        "native_L": run_to_steady_state(p, native),
        "mirror_D": run_to_steady_state(p, mirror),
        "nonadhesive_PEG": run_to_steady_state(p, peg),
    }


# --------------------------------------------------------------------------
# plotting helpers (matplotlib, Agg backend)
# --------------------------------------------------------------------------

def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_dose_response(dr, path):
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for name in ["Cint", "FAKp", "Tension", "YAPn", "Col"]:
        ax.plot(dr["chi"], dr[name], label=name, lw=2)
    ax.plot(dr["chi"], dr["phenotype"], "k--", lw=2.5, label="phenotype index")
    ax.set_xlabel("environmental chirality  chi  (0 = all-L  ->  1 = all-D mirror)")
    ax.set_ylabel("steady-state activity  [0, 1]")
    ax.set_title("Tenocyte mechanotransduction vs. ECM chirality")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_stiffness_grid(curves, path):
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for e, dr in curves.items():
        ax.plot(dr["chi"], dr["phenotype"], lw=2, label=f"{e:g} kPa")
    ax.set_xlabel("environmental chirality  chi")
    ax.set_ylabel("phenotype index  [0, 1]")
    ax.set_title("Chirality x stiffness interaction")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(title="substrate stiffness", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_preset_comparison(path):
    """Side-by-side hysteresis sweeps for the monostable and bistable regimes.

    Makes the conditional nature of the hysteresis prediction visible: the
    parsimonious defaults give a sharp but fully reversible transition, while
    stronger feedback opens a history-dependent gap.
    """
    from .parameters import get_preset
    from .analysis import detect_bistability

    plt = _mpl()
    presets = ["parsimonious", "strong_feedback"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, name in zip(axes, presets):
        b = detect_bistability(get_preset(name))
        ax.plot(b["chi"], b["forward"], "-o", ms=3, lw=2, label="forward  L -> D")
        ax.plot(b["chi"], b["reverse"], "-s", ms=3, lw=2, label="reverse  D -> L")
        tag = ("bistable" if b["bistable"] else "monostable")
        ax.set_title(f"{name}  ({tag}, gap = {b['max_gap']:.3f})")
        ax.set_xlabel("environmental chirality  chi")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("phenotype index  [0, 1]")
    fig.suptitle("Is the chirality transition reversible? Depends on feedback strength",
                 y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_hysteresis(h, path):
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(h["chi_up"], h["phen_up"], "-o", ms=3, lw=2, label="forward  L -> D")
    ax.plot(h["chi_down"], h["phen_down"], "-s", ms=3, lw=2, label="reverse  D -> L")
    ax.set_xlabel("environmental chirality  chi")
    ax.set_ylabel("phenotype index  [0, 1]")
    ax.set_title("Hysteresis sweep (gap => bistable, history-dependent phenotype)")
    ax.set_ylim(-0.02, 1.02)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
