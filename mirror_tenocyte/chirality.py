"""How environmental chirality couples into the model.

The central biophysical fact: integrin-ligand recognition is *stereospecific*.
Integrins read chiral epitopes (e.g. the GFOGER motif of collagen via the
alpha1/alpha2/alpha11 I-domain, or RGD in fibronectin/vitronectin) through a
lock-and-key interface. A mirror-image (all-D) ligand presents the enantiomeric
surface, which the natural, chiral integrin cannot productively engage.

Consequences encoded here:

  * `effective_ligand` collapses the chirality knob `chi` into the *bindable*
    ligand density that the integrin module actually sees. In 'racemic' mode a
    fully mirror matrix (chi = 1) yields zero bindable ligand — i.e. it is,
    to first order, indistinguishable from a non-adhesive substrate. That
    equivalence is itself a testable prediction and a required experimental
    control (D-ECM vs. passivated/PEG surface).

  * `chirality_specific_gain` is a *separate, optional* channel for the
    hypothesis that a D-environment does something beyond loss of recognition
    (e.g. altered supramolecular fibril assembly changing local mechanics/
    topology, engagement of chirality-sensing receptors, or differential
    protein adsorption). It defaults to zero so the baseline model stays
    parsimonious; turn it up to explore "chirality is not merely adhesion loss".
"""

from __future__ import annotations

import math


def effective_ligand(env: dict, p: dict) -> float:
    """Return the productively-bindable ligand density seen by integrins.

    Depends on total ligand density `Ltot` and the chirality composition `chi`
    through the selected `chir_mode`.
    """
    chi = float(env["chi"])
    chi = min(max(chi, 0.0), 1.0)
    ltot = float(env["Ltot"])
    mode = env.get("chir_mode", "racemic")

    if mode == "racemic":
        # Fraction chi of ligands are D and non-bindable.
        return ltot * (1.0 - chi)
    if mode == "affinity":
        # Every ligand's productive affinity decays smoothly with inversion.
        return ltot * math.exp(-float(p["beta_chi"]) * chi)
    raise ValueError(f"unknown chir_mode: {mode!r}")


def chirality_specific_gain(env: dict, p: dict) -> float:
    """Optional additive signal representing D-ECM effects *beyond* adhesion loss.

    Returns a non-negative term (scaled by p['chi_specific']) that is fed into
    the YAP-activating signal. With the default chi_specific = 0 this is inert.
    Use it to encode a hypothesis where reversed chirality actively perturbs
    signaling independently of integrin occupancy.
    """
    gain = float(p.get("chi_specific", 0.0))
    if gain == 0.0:
        return 0.0
    chi = min(max(float(env["chi"]), 0.0), 1.0)
    return gain * chi
