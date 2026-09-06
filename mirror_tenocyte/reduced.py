"""The model's actual content, stated plainly: a surface over two scalars.

Why this module exists
----------------------
A critical audit of the project found that the six-ODE cascade does far less
work than its size suggests. Inspect `model.derivatives`: the environment
enters ONLY through two quantities,

    L_eff = chirality.effective_ligand(env, p)      bindable ligand density
    E_eff = assembly.effective_stiffness(env, p)    stiffness after assembly

Nothing else about `env` is read. Therefore every steady state — and every
prediction in this project — is a function of exactly those two numbers:

    phenotype = F(L_eff, E_eff)

That is a structural fact, not an approximation. It has three consequences,
and they are the reason to keep this module rather than hide it:

1. Honesty. Three assembly hypotheses that produce the same E_eff produce
   *bit-identical* phenotypes. The cascade cannot distinguish mechanism, only
   magnitude. Claims about "integrin/FAK/YAP signalling" distinguishing
   chirality hypotheses would be overclaiming; the discriminating content lives
   entirely in the assembly layer.

2. Speed. Evaluating F by interpolation is ~1000x faster than integrating the
   ODEs, which makes dense parameter sweeps cheap.

3. Falsifiability. If a future version of the model adds a channel that is NOT
   mediated by these two scalars — a chirality-sensing receptor, a
   topology-sensitive term, a non-monotone feedback — then the reduction will
   BREAK, and `verify_reduction` will say so. The reduction failing is the
   signal that the cascade has started to earn its keep.

The cascade is not useless: it supplies the *shape* of F (saturation, the
switch, the relative weighting of adhesion and tension) and it names
intermediate observables that an experimentalist can actually stain for. But it
adds no degrees of freedom beyond F itself.
"""

from __future__ import annotations

import numpy as np

from .parameters import get_default_parameters, get_default_environment
from .chirality import effective_ligand
from .assembly import effective_stiffness
from .model import run_to_steady_state, phenotype_index, STATE_NAMES


def env_coordinates(env: dict, p: dict) -> tuple:
    """The only two numbers the ODE model reads out of `env`."""
    return float(effective_ligand(env, p)), float(effective_stiffness(env, p))


def _solve_at(L_eff: float, E_eff: float, p: dict):
    """Steady state at given effective coordinates, via a canonical env."""
    env = get_default_environment()
    env["chir_mode"] = "racemic"
    env["chi"] = 0.0            # so effective_ligand == Ltot
    env["Ltot"] = float(L_eff)
    env["assembly_mode"] = "none"   # so effective_stiffness == Estiff
    env["Estiff"] = float(E_eff)
    return run_to_steady_state(p, env)


def verify_reduction(p=None, cases=None, tol=1e-9) -> dict:
    """Check that different environments with equal (L_eff, E_eff) agree exactly.

    This is the test that keeps the reduction honest. If someone adds a
    chirality channel that bypasses the two scalars, this fails — which is the
    desired behaviour, not a regression.
    """
    p = p or get_default_parameters()

    if cases is None:
        cases = []
        # same E_eff reached through three different assembly hypotheses
        for chi_s, base in [(0.05, 30.0), (0.2, 30.0), (0.5, 8.0)]:
            for mode in ["frustrated_kinetic", "racemic_enhanced", "self_sorting"]:
                e = get_default_environment()
                e["assembly_mode"] = mode
                e["Estiff"] = base
                e["chi"] = 0.0
                e["chi_struct"] = chi_s
                cases.append(e)
        # same L_eff reached through different chi/Ltot combinations
        for prod in [0.3, 0.085]:
            for Ltot in [0.5, 1.0, 2.0]:
                e = get_default_environment()
                e["Ltot"] = Ltot
                e["chi"] = 1.0 - prod / Ltot
                cases.append(e)

    worst = 0.0
    checked = 0
    for env in cases:
        L, E = env_coordinates(env, p)
        full = run_to_steady_state(p, env)
        red = _solve_at(L, E, p)
        worst = max(worst, float(np.max(np.abs(full - red))))
        checked += 1
    return {"n_cases": checked, "max_deviation": worst,
            "reduction_holds": bool(worst < tol)}


class ReducedModel:
    """Interpolated F(L_eff, E_eff) — the whole model, as a lookup surface."""

    def __init__(self, p=None, n_L=36, n_E=36,
                 L_range=(1e-3, 4.0), E_range=(1e-2, 300.0)):
        self.p = p or get_default_parameters()
        self.L_grid = np.concatenate([[0.0],
                                      np.logspace(np.log10(L_range[0]),
                                                  np.log10(L_range[1]), n_L)])
        self.E_grid = np.concatenate([[0.0],
                                      np.logspace(np.log10(E_range[0]),
                                                  np.log10(E_range[1]), n_E)])
        self.phen = np.zeros((self.L_grid.size, self.E_grid.size))
        self.states = np.zeros((self.L_grid.size, self.E_grid.size,
                                len(STATE_NAMES)))
        for i, L in enumerate(self.L_grid):
            for j, E in enumerate(self.E_grid):
                s = _solve_at(L, E, self.p)
                self.states[i, j] = s
                self.phen[i, j] = phenotype_index(s)

    def _bilinear(self, table, L, E):
        i = np.clip(np.searchsorted(self.L_grid, L) - 1, 0, self.L_grid.size - 2)
        j = np.clip(np.searchsorted(self.E_grid, E) - 1, 0, self.E_grid.size - 2)
        L0, L1 = self.L_grid[i], self.L_grid[i + 1]
        E0, E1 = self.E_grid[j], self.E_grid[j + 1]
        tl = 0.0 if L1 == L0 else (L - L0) / (L1 - L0)
        te = 0.0 if E1 == E0 else (E - E0) / (E1 - E0)
        tl = min(max(tl, 0.0), 1.0)
        te = min(max(te, 0.0), 1.0)
        return ((1 - tl) * (1 - te) * table[i, j]
                + tl * (1 - te) * table[i + 1, j]
                + (1 - tl) * te * table[i, j + 1]
                + tl * te * table[i + 1, j + 1])

    def phenotype(self, env: dict) -> float:
        L, E = env_coordinates(env, self.p)
        return float(self._bilinear(self.phen, L, E))

    def state(self, env: dict):
        L, E = env_coordinates(env, self.p)
        return np.array([self._bilinear(self.states[:, :, k], L, E)
                         for k in range(len(STATE_NAMES))])
