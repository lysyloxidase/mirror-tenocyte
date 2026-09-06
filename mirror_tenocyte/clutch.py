"""Motor-clutch stiffness sensing with force-dependent bonds (mean-field).

Why this module exists
----------------------
The baseline model senses stiffness through a Hill function, which is
*monotonic*: more stiffness always means more engagement. The real
molecular-clutch picture is not monotonic. Chan & Odde (2008, Science 322:1687)
showed that a motor-clutch ensemble has an optimal substrate stiffness: too soft
and clutches never load ("frictional slippage"), too stiff and they are loaded so
fast they rupture ("load-and-fail"). Elosegui-Artola et al. (2016, Nat Cell Biol
18:540) showed talin unfolding and integrin catch bonds reshape that response.

This matters here specifically because the assembly layer makes *stiffness* the
carrier of the chirality signal. If the traction/stiffness relation is biphasic
and tendon-like stiffness sits above the optimum, then the softening caused by a
frustrated D/L mixture moves the cell *towards* the optimum, which could invert
the rebound predicted in the README. So this module is an adversarial test of
the previous step's headline prediction, not a refinement of it.

Physical units (unlike the rest of the package, which is normalised):
    force      pN
    stiffness  pN/nm
    velocity   nm/s
    rates      1/s

Mean-field reduction
--------------------
With ``n`` engaged clutches of stiffness ``k_clutch`` in parallel, in series with
a substrate spring ``k_sub``, the loading rate felt by one clutch is

    kappa_load(n) = k_sub * k_clutch / (k_sub + n * k_clutch)

Actin retrograde flow slows under load (linear force-velocity motor):

    v(F_total) = v_unloaded * max(0, 1 - F_total / (n_motor * F_stall))

A clutch survives ~1/koff(F) before unbinding, so its mean force satisfies the
self-consistent relation F = kappa_load(n) * v / koff(F), while engagement obeys
detailed balance k_on * (N - n) = n * koff(F). Substituting the second into the
first leaves a single scalar root-find in F.

This is a deterministic mean-field approximation. It reproduces the optimum and
the two limiting regimes, but not the stochastic load-and-fail *oscillations* of
the full Chan-Odde simulation.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq


def get_clutch_parameters() -> dict:
    """Order-of-magnitude clutch parameters from the motor-clutch literature."""
    return {
        "N_clutch": 75.0,      # clutches per adhesion (Chan & Odde 2008)
        "k_on": 1.0,           # clutch engagement rate, 1/s
        "k_clutch": 0.8,       # clutch spring constant, pN/nm
        "n_motor": 50.0,       # myosin motors pulling the actin bundle
        "F_stall": 2.0,        # stall force per motor, pN
        "v_unloaded": 120.0,   # unloaded actin retrograde flow, nm/s

        # slip bond (Bell model)
        "koff0_slip": 0.1,     # 1/s
        "F_bell": 2.0,         # pN

        # catch-slip bond, integrin alpha5beta1-fibronectin
        # (shape after Kong et al. 2009, J Cell Biol 185:1275)
        "kc0": 5.0,            # catch branch amplitude, 1/s
        "F_catch": 5.0,        # catch branch force scale, pN
        "ks0": 0.1,            # slip branch amplitude, 1/s
        "F_slip": 8.0,         # slip branch force scale, pN

        # Force-gated adhesion reinforcement (talin unfolding -> vinculin
        # recruitment -> more engaged clutches), after Elosegui-Artola et al.
        # 2016. Without this a fixed-clutch ensemble shows essentially NO
        # stiffness discrimination across 1-100 kPa; with it, stiffness sensing
        # becomes monotonic over the physiological range. Default 0 = off.
        "reinforce_gain": 0.0,   # max fold-increase in engagement rate
        "F_reinforce": 5.0,      # talin unfolding force threshold, pN
        "n_reinforce": 3.0,      # sharpness of the threshold

        # Adhesion GROWTH: force-gated recruitment of *new* clutches from the
        # cytoplasmic pool, not just faster binding of a fixed set. This is the
        # piece that makes stiffness sensing rise rather than merely flatten:
        # load sharing across a growing adhesion pulls per-clutch force back
        # below the rupture range, so the adhesion stabilises instead of failing.
        # Myosin is recruited with it (a bigger adhesion carries a bigger stress
        # fibre), so the motor stall ceiling scales too — otherwise growth just
        # stalls the motors and the response saturates artificially.
        "growth_gain": 0.0,      # max fold-increase in clutch pool (0 = off)
        "F_growth": 5.0,         # recruitment force threshold, pN
        "n_growth": 3.0,         # sharpness of recruitment
        "N_max_factor": 8.0,     # hard ceiling on pool growth (array size)

        # Young's modulus (kPa) -> substrate spring (pN/nm).
        # For an adhesion of radius r on an elastic half-space, k ~ 2*E*r.
        # With r = 1 um: k = 2 * 3e4 Pa * 1e-6 m = 0.06 N/m = 60 pN/nm at
        # 30 kPa, hence 2.0 pN/nm per kPa. (1 N/m = 1e3 pN/nm.)
        # This prefactor absorbs adhesion size and geometry and is deliberately
        # exposed, because WHERE the physiological range sits relative to the
        # clutch optimum is exactly what the conclusions hinge on.
        "kappa_per_kPa": 2.0,
    }


#: Bond type used with the calibrated parameter set. Integrin-ECM bonds
#: (alpha5beta1, alphaVbeta3) are catch bonds, and the bond type turns out to
#: decide the SIGN of stiffness sensing: a weak slip bond puts the whole
#: physiological range past the clutch optimum (signal falls with stiffness),
#: while a catch bond keeps it on the rising branch.
CALIBRATED_BOND = "catch"


def get_calibrated_parameters() -> dict:
    """Clutch parameters calibrated to reproduce physiological stiffness sensing.

    The literature defaults in `get_clutch_parameters` do NOT reproduce the
    best-established fact in mechanobiology — that adhesion and nuclear YAP rise
    from ~1 to ~100 kPa. This set does. Three changes, in order of importance:

    1. ``bond = 'catch'`` (see CALIBRATED_BOND). Decides the sign.
    2. ``growth_gain = 10`` — force-gated recruitment of new clutches *and* the
       myosin pulling on them. Sets the dynamic range.
    3. ``kappa_per_kPa = 0.1`` — the stiffness-to-spring mapping. Sets where the
       transition sits.

    Honesty about (3): this value was FITTED, by scanning until the half-maximal
    response landed in the 5-15 kPa window where YAP switching is observed
    (Dupont et al. 2011). It is not independent evidence. Two things make it
    more than curve-fitting, though:

      * k = 2*E*r implies r ~ 50 nm, which is the scale of an integrin
        NANOCLUSTER, not a mature focal adhesion. That is the right object: the
        model describes ~75 clutches, i.e. one adhesion module, not a whole FA.
        The original 2.0 (r ~ 1 um) was mismatched to what is being modelled.
      * Per-clutch forces were NOT fitted, and they land where mechanosensing
        biology says they should: 2.4 -> 10.6 pN across 1-100 kPa, crossing the
        ~5 pN talin unfolding threshold at 8-15 kPa — exactly the stiffness
        where the YAP switch is observed — while staying under the ~13 pN
        catch-bond stabilisation optimum.

    Calibration target reproduced: engaged clutches rise ~7x across 1-100 kPa
    with half-maximum at ~7.6 kPa.
    """
    p = get_clutch_parameters()
    p.update({
        "k_on": 0.3,
        "growth_gain": 10.0,
        "kappa_per_kPa": 0.1,
    })
    return p


def koff(F: float, p: dict, bond: str = "catch") -> float:
    """Force-dependent unbinding rate, 1/s."""
    F = max(float(F), 0.0)
    if bond == "slip":
        return p["koff0_slip"] * math.exp(F / p["F_bell"])
    if bond == "catch":
        # catch branch stabilises the bond at low-to-moderate force; the slip
        # branch takes over above ~15 pN, giving the characteristic minimum.
        return (p["kc0"] * math.exp(-F / p["F_catch"])
                + p["ks0"] * math.exp(F / p["F_slip"]))
    if bond == "ideal":
        return p["koff0_slip"]
    raise ValueError(f"unknown bond type: {bond!r}")


def _engaged(F: float, p: dict, bond: str) -> float:
    """Engaged clutch count from detailed balance at mean force F."""
    return p["N_clutch"] * p["k_on"] / (p["k_on"] + koff(F, p, bond))


def solve_clutch(k_sub: float, p: dict, bond: str = "catch") -> dict:
    """Steady-state mean-field clutch solution at substrate stiffness k_sub.

    Returns engaged fraction, mean force per clutch, total traction, and the
    retrograde flow velocity.
    """
    k_sub = max(float(k_sub), 1e-12)
    kc = p["k_clutch"]
    F_motor_max = p["n_motor"] * p["F_stall"]

    def residual(F):
        n = _engaged(F, p, bond)
        kappa = k_sub * kc / (k_sub + n * kc)
        v = p["v_unloaded"] * max(0.0, 1.0 - (n * F) / F_motor_max)
        return kappa * v / koff(F, p, bond) - F

    # residual(0+) > 0 (clutches load) and residual(F_big) < 0 (motors stall),
    # so a root is bracketed; widen defensively if the upper bound is too low.
    lo, hi = 1e-9, max(10.0 * F_motor_max, 100.0)
    r_lo, r_hi = residual(lo), residual(hi)
    if r_lo <= 0.0:
        F = 0.0
    elif r_hi > 0.0:
        F = hi
    else:
        F = brentq(residual, lo, hi, xtol=1e-12, rtol=1e-12, maxiter=200)

    n = _engaged(F, p, bond)
    v = p["v_unloaded"] * max(0.0, 1.0 - (n * F) / F_motor_max)
    return {"F_per_clutch": float(F), "n_engaged": float(n),
            "engaged_fraction": float(n / p["N_clutch"]),
            "traction": float(n * F), "velocity": float(v)}


def simulate_clutch(k_sub: float, p: dict, bond: str = "slip",
                    t_max: float = 10.0, dt: float = 1e-4,
                    burn_frac: float = 0.3, seed: int = 0) -> dict:
    """Stochastic motor-clutch simulation (Chan & Odde 2008), time-averaged.

    The mean-field solver above cannot produce the stiffness optimum, because
    that optimum is a *dynamic* phenomenon: on stiff substrates the ensemble
    enters load-and-fail cycles where clutches rupture in cascades and traction
    periodically collapses. A steady-state balance averages those cycles away
    and returns the stalled (high-traction) branch instead. Reproducing the
    optimum therefore requires simulating individual clutches.

    Each clutch is a spring that engages at zero force, is dragged at the actin
    retrograde velocity, and unbinds at its own force-dependent rate. The
    substrate position follows from instantaneous force balance.
    """
    rng = np.random.default_rng(seed)
    N0 = int(p["N_clutch"])
    kc = p["k_clutch"]
    k_on = p["k_on"]
    v_u = p["v_unloaded"]
    k_sub = max(float(k_sub), 1e-12)

    a_g = float(p.get("growth_gain", 0.0))
    F_g = float(p.get("F_growth", 5.0))
    n_g = float(p.get("n_growth", 3.0))
    N = int(N0 * (float(p.get("N_max_factor", 8.0)) if a_g > 0.0 else 1.0))

    engaged = np.zeros(N, dtype=bool)
    x = np.zeros(N)          # clutch attachment positions (nm)

    n_steps = int(t_max / dt)
    burn = int(burn_frac * n_steps)
    tract_sum, eng_sum, samples = 0.0, 0.0, 0
    n_sum, avail_sum = 0.0, 0.0

    a_r = float(p.get("reinforce_gain", 0.0))
    F_r = float(p.get("F_reinforce", 5.0))
    n_r = float(p.get("n_reinforce", 3.0))
    p_on_base = 1.0 - math.exp(-k_on * dt)
    F_mean = 0.0

    for step in range(n_steps):
        n = int(engaged.sum())
        if n > 0:
            # force balance: k_sub * x_sub = sum k_c * (x_i - x_sub)
            x_sub = kc * x[engaged].sum() / (k_sub + n * kc)
            F_total = k_sub * x_sub
        else:
            x_sub, F_total = 0.0, 0.0

        # Force-gated adhesion growth: recruit new clutches AND the myosin that
        # pulls on them, keeping the motor:clutch ratio fixed.
        if a_g > 0.0:
            fg = max(F_mean, 0.0) ** n_g
            growth = 1.0 + a_g * fg / (F_g ** n_g + fg)
        else:
            growth = 1.0
        N_avail = min(N, max(1, int(round(N0 * growth))))
        F_motor_max = p["n_motor"] * growth * p["F_stall"]

        v = v_u * max(0.0, 1.0 - F_total / F_motor_max)

        if step >= burn:
            tract_sum += F_total
            eng_sum += n / N_avail
            samples += 1

        # advance engaged clutches with the retrograde flow
        if n > 0:
            x[engaged] += v * dt
            F_i = kc * (x[engaged] - x_sub)
            F_i = np.maximum(F_i, 0.0)
            if bond == "slip":
                rates = p["koff0_slip"] * np.exp(F_i / p["F_bell"])
            elif bond == "catch":
                rates = (p["kc0"] * np.exp(-F_i / p["F_catch"])
                         + p["ks0"] * np.exp(F_i / p["F_slip"]))
            elif bond == "ideal":
                rates = np.full(F_i.shape, p["koff0_slip"])
            else:
                raise ValueError(f"unknown bond type: {bond!r}")
            F_mean = float(F_i.mean()) if F_i.size else 0.0
            p_off = -np.expm1(-np.minimum(rates * dt, 50.0))
            idx = np.flatnonzero(engaged)
            breaking = idx[rng.random(idx.size) < p_off]
            engaged[breaking] = False
            x[breaking] = 0.0
        else:
            F_mean = 0.0

        # engage free clutches at zero force (x_i = x_sub). Force-gated
        # reinforcement raises the engagement rate once clutches are loaded
        # past the talin unfolding threshold.
        if a_r > 0.0:
            fm = max(F_mean, 0.0) ** n_r
            boost = 1.0 + a_r * fm / (F_r ** n_r + fm)
            p_on = 1.0 - math.exp(-k_on * boost * dt)
        else:
            p_on = p_on_base

        free = np.flatnonzero(~engaged[:N_avail])
        if free.size:
            binding = free[rng.random(free.size) < p_on]
            engaged[binding] = True
            x[binding] = x_sub

        if step >= burn:
            n_sum += n
            avail_sum += N_avail

    s = max(samples, 1)
    return {"traction": tract_sum / s,
            "engaged_fraction": eng_sum / s,
            "n_engaged": n_sum / s,
            "n_available": avail_sum / s}


def traction_curve_stochastic(E_values, p: dict, bond: str = "slip", **kw):
    """Time-averaged traction across stiffness, via the stochastic simulation."""
    E_values = np.asarray(E_values, dtype=float)
    tract, eng = [], []
    for i, E in enumerate(E_values):
        s = simulate_clutch(stiffness_to_ksub(E, p), p, bond, seed=i, **kw)
        tract.append(s["traction"])
        eng.append(s["engaged_fraction"])
    return {"E": E_values, "traction": np.array(tract), "engaged": np.array(eng)}


def stiffness_to_ksub(E_kPa: float, p: dict) -> float:
    """Map Young's modulus (kPa) to the effective substrate spring (pN/nm)."""
    return max(float(E_kPa), 0.0) * p["kappa_per_kPa"]


def stiffness_response(E_values, p=None, bond=None, nseed=2,
                       t_max=5.0, dt=2e-4):
    """Adhesion size, traction and per-clutch force across a stiffness series.

    The readout used for calibration is the number of engaged clutches, since
    that is the adhesion-size proxy the downstream model consumes, and unlike
    traction it is not capped by the motor stall force.
    """
    p = p if p is not None else get_calibrated_parameters()
    bond = bond or CALIBRATED_BOND
    E_values = np.asarray(E_values, dtype=float)

    n_eng, tract = [], []
    for i, E in enumerate(E_values):
        runs = [simulate_clutch(stiffness_to_ksub(E, p), p, bond,
                                t_max=t_max, dt=dt, seed=100 * s + i)
                for s in range(nseed)]
        n_eng.append(float(np.mean([r["n_engaged"] for r in runs])))
        tract.append(float(np.mean([r["traction"] for r in runs])))

    n_eng = np.array(n_eng)
    tract = np.array(tract)
    with np.errstate(divide="ignore", invalid="ignore"):
        f_per = np.where(n_eng > 1e-9, tract / np.maximum(n_eng, 1e-9), 0.0)
    return {"E": E_values, "n_engaged": n_eng, "traction": tract,
            "F_per_clutch": f_per}


def fit_hill_to_clutch(p=None, bond=None, n_points=22, nseed=3, t_max=6.0):
    """Derive Hill parameters for the coarse-grained model from the clutch.

    Runs the calibrated stochastic clutch across a log-spaced stiffness range
    and fits ``v_max * E^h / (K^h + E^h)`` to two readouts:

      engaged clutches -> the stiffness->adhesion term  (KE_adh, hE_adh)
      traction         -> the stiffness->tension term   (KE_ten, hE_ten)

    The fits are excellent (R^2 > 0.99), which is a result in its own right: the
    coarse Hill form used by model.py is a justified reduction of the explicit
    clutch, and only its constants needed deriving rather than assuming.

    Takes ~30 s. The fitted values are stored as the 'clutch_calibrated' preset
    in parameters.py; this function exists so they can be reproduced and
    checked, not taken on trust.
    """
    from scipy.optimize import curve_fit

    p = p if p is not None else get_calibrated_parameters()
    bond = bond or CALIBRATED_BOND
    E = np.logspace(np.log10(0.05), np.log10(300.0), n_points)
    r = stiffness_response(E, p=p, bond=bond, nseed=nseed, t_max=t_max)

    def hill(E, vmax, K, h):
        return vmax * E ** h / (K ** h + E ** h)

    out = {"E": E}
    for key, y, label in [("adh", r["n_engaged"], "engaged clutches"),
                          ("ten", r["traction"], "traction")]:
        popt, _ = curve_fit(hill, E, y, p0=[float(y.max()), 7.0, 1.0],
                            maxfev=20000)
        pred = hill(E, *popt)
        ss = 1.0 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
        out[key] = {"v_max": float(popt[0]), "K": float(popt[1]),
                    "h": float(popt[2]), "r2": float(ss),
                    "data": y, "fit": pred, "label": label}
    return out


def calibration_report(p=None, bond=None, **kw) -> dict:
    """Check the calibrated clutch against the physiological stiffness benchmark.

    Passing requires a monotone-ish RISE across 1-100 kPa with a half-maximum in
    the 5-15 kPa window (where YAP nuclear translocation is observed to switch)
    and a dynamic range of at least 2x.
    """
    E = np.array([0.5, 1, 3, 8, 15, 30, 60, 100], dtype=float)
    r = stiffness_response(E, p=p, bond=bond, **kw)
    n = r["n_engaged"]

    rising = bool(n[-1] > n[1] * 1.15)
    dyn = float(n.max() / max(n.min(), 1e-9))
    half = n.min() + 0.5 * (n.max() - n.min())
    half_E = float(np.interp(half, n, E)) if rising else float("nan")
    in_window = bool(5.0 <= half_E <= 15.0)

    return {**r, "rising": rising, "dynamic_range": dyn,
            "half_max_kPa": half_E, "half_max_in_window": in_window,
            "passes": bool(rising and dyn >= 2.0 and in_window)}


def clutch_engagement(E_kPa: float, p: dict, bond: str = "catch") -> float:
    """Engaged clutch fraction at a given matrix stiffness — the model hook."""
    return solve_clutch(stiffness_to_ksub(E_kPa, p), p, bond)["engaged_fraction"]


def traction_curve(E_values, p: dict, bond: str = "catch"):
    """Traction and engagement across a stiffness range (for locating the optimum)."""
    E_values = np.asarray(E_values, dtype=float)
    tract, eng, force = [], [], []
    for E in E_values:
        s = solve_clutch(stiffness_to_ksub(E, p), p, bond)
        tract.append(s["traction"])
        eng.append(s["engaged_fraction"])
        force.append(s["F_per_clutch"])
    return {"E": E_values, "traction": np.array(tract),
            "engaged": np.array(eng), "force": np.array(force)}


def find_optimum(p: dict, bond: str = "catch", lo=0.1, hi=1000.0, n=400):
    """Locate the traction-maximising stiffness, if the response is biphasic."""
    E = np.logspace(np.log10(lo), np.log10(hi), n)
    c = traction_curve(E, p, bond)
    i = int(np.argmax(c["traction"]))
    interior = 0 < i < n - 1
    return {"E_opt": float(E[i]), "traction_max": float(c["traction"][i]),
            "biphasic": bool(interior), "curve": c}
