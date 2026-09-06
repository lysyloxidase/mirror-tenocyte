"""The coupled ODE model of tenocyte mechanotransduction under chiral ECM.

State vector (all normalised to [0, 1]):
    y = [Cint, FAKp, Rho, Tension, YAPn, Col]

The right-hand side is intentionally simple: each node relaxes toward a
set-point determined by its upstream driver, with two positive-feedback loops
that give the network its interesting (potentially switch-like) behaviour:

  (1) Tension -> integrin engagement   (talin unfolding / clutch reinforcement)
  (2) Nuclear YAP -> its own retention  (cytoskeletal-gene feedback)

Stiffness enters twice, as in real mechanosensing: it is needed both to load
the clutch (KE_adh) and to let actomyosin build tension against a resistant
substrate (KE_ten). On a very soft matrix, tension collapses even with plenty
of ligand — matching the classic soft-substrate loss of nuclear YAP.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .chirality import effective_ligand, chirality_specific_gain
from .assembly import effective_stiffness

STATE_NAMES = ["Cint", "FAKp", "Rho", "Tension", "YAPn", "Col"]


def _hill(x: float, k: float, n: float) -> float:
    """Hill activation; guards against tiny negatives from the integrator."""
    x = max(x, 0.0)
    xn = x ** n
    return xn / (k ** n + xn)


def derivatives(t: float, y, p: dict, env: dict):
    """ODE right-hand side. Signature matches scipy.integrate.solve_ivp."""
    Cint, FAKp, Rho, Tension, YAPn, Col = y

    # Stiffness is what the matrix actually delivers *after* assembly quality,
    # so a frustrated D/L mixture is mechanically softer than either pure form.
    stiff = effective_stiffness(env, p)
    # Hill exponents default to 1, reproducing the original E/(E+K) form. The
    # 'clutch_calibrated' preset supplies values fitted to the explicit
    # motor-clutch simulation (see clutch.py), which the Hill form describes
    # with R^2 > 0.99 — so the coarse-grained shape is justified, and only its
    # constants needed correcting.
    fE_adh = _hill(stiff, p["KE_adh"], p.get("hE_adh", 1.0))
    fE_ten = _hill(stiff, p["KE_ten"], p.get("hE_ten", 1.0))

    # --- integrin engagement (chirality-gated ligand + tension feedback) ---
    l_eff = effective_ligand(env, p)
    reinforce = 1.0 + p["a_reinforce"] * max(Tension, 0.0)
    drive = l_eff * reinforce * fE_adh
    cint_set = _hill(drive, p["K_adh"], p["n_adh"])
    dCint = p["kon"] * cint_set * (1.0 - Cint) - p["koff"] * Cint

    # --- FAK (Y397) ---
    dFAK = p["kf"] * Cint * (1.0 - FAKp) - p["kdf"] * FAKp

    # --- RhoA/ROCK ---
    dRho = p["kr"] * FAKp * (1.0 - Rho) - p["kdr"] * Rho

    # --- cytoskeletal tension (needs a resistant substrate) ---
    dTen = p["kt"] * Rho * fE_ten * (1.0 - Tension) - p["kdt"] * Tension

    # --- YAP/TAZ nuclear translocation (switch + self-reinforcement) ---
    signal = (p["w_ten"] * Tension + p["w_fak"] * FAKp)
    signal *= (1.0 + p["a_yap"] * max(YAPn, 0.0))
    signal += chirality_specific_gain(env, p)   # optional non-adhesion channel
    yap_set = _hill(signal, p["K_yap"], p["n_yap"])
    dYAP = p["ky"] * yap_set * (1.0 - YAPn) - p["kdy"] * YAPn

    # --- collagen-I production proxy ---
    dCol = p["kc"] * YAPn - p["kdc"] * Col

    return [dCint, dFAK, dRho, dTen, dYAP, dCol]


def run_to_steady_state(p: dict, env: dict, y0=None, t_max: float = 200.0):
    """Integrate to (approximate) steady state; return the final state vector.

    y0 lets you seed the initial condition — useful for hysteresis sweeps where
    you continue from the previous environment's resting state.
    """
    if y0 is None:
        y0 = np.full(len(STATE_NAMES), 0.05)
    sol = solve_ivp(
        derivatives,
        (0.0, t_max),
        y0,
        args=(p, env),
        method="LSODA",
        rtol=1e-8,
        atol=1e-10,
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    final = np.clip(sol.y[:, -1], 0.0, 1.0)
    return final


def run_timecourse(p: dict, env: dict, y0=None, t_max: float = 40.0, n: int = 400):
    """Integrate and return (t, Y) with Y shape (len(STATE_NAMES), n)."""
    if y0 is None:
        y0 = np.full(len(STATE_NAMES), 0.05)
    t_eval = np.linspace(0.0, t_max, n)
    sol = solve_ivp(
        derivatives,
        (0.0, t_max),
        y0,
        args=(p, env),
        method="LSODA",
        t_eval=t_eval,
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    return sol.t, np.clip(sol.y, 0.0, 1.0)


def phenotype_index(state) -> float:
    """Collapse a steady state into a single 'tenogenic activation' score in [0,1].

    Weighted toward the transcriptional/matrix output (nuclear YAP and collagen)
    but including tension, since a spread/contractile tenocyte maintaining matrix
    is the phenotype of interest. High score ~ activated matrix-producing
    tenocyte; low score ~ rounded / quiescent / dedifferentiating cell.
    """
    state = np.asarray(state, dtype=float)
    _, _, _, tension, yapn, col = state
    return float(np.clip(0.5 * col + 0.3 * yapn + 0.2 * tension, 0.0, 1.0))


def classify_phenotype(state) -> str:
    """Human-readable label for a steady state (illustrative thresholds)."""
    idx = phenotype_index(state)
    if idx >= 0.6:
        return "activated tenogenic (nuclear YAP, high collagen)"
    if idx >= 0.3:
        return "intermediate / mechanically primed"
    return "quiescent / dedifferentiated (cytoplasmic YAP, low collagen)"
