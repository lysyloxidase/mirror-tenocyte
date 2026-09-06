"""Quantitative analyses that turn the README's claims into checkable numbers.

Each function here answers one question about the model's behaviour, so that
predictions are *measured* from the simulation rather than asserted in prose.
"""

from __future__ import annotations

import copy

import numpy as np

from .parameters import get_default_parameters, get_default_environment
from .model import run_to_steady_state, phenotype_index


def _phen_curve(p, env, chi_values):
    """Phenotype index across chi, each solved from a naive initial condition."""
    env = copy.deepcopy(env)
    out = []
    for c in chi_values:
        env["chi"] = float(c)
        out.append(phenotype_index(run_to_steady_state(p, env)))
    return np.array(out)


def find_chi_star(p=None, env=None, n=201):
    """Locate the critical chirality chi* (half-maximal phenotype crossing).

    Returns a dict with chi*, the phenotype range, and a sharpness measure
    (max |d phenotype / d chi|). A large range with a large sharpness means a
    switch-like transition; a small range means chirality barely matters here.
    """
    p = p or get_default_parameters()
    env = env or get_default_environment()
    chi = np.linspace(0.0, 1.0, n)
    phen = _phen_curve(p, env, chi)

    lo, hi = float(phen.min()), float(phen.max())
    rng = hi - lo
    if rng < 1e-3:
        return {"chi_star": None, "phen_lo": lo, "phen_hi": hi,
                "range": rng, "sharpness": 0.0, "chi": chi, "phen": phen}

    half = lo + 0.5 * rng
    idx = int(np.argmin(np.abs(phen - half)))
    sharp = float(np.max(np.abs(np.diff(phen))) / (chi[1] - chi[0]))
    return {"chi_star": float(chi[idx]), "phen_lo": lo, "phen_hi": hi,
            "range": rng, "sharpness": sharp, "chi": chi, "phen": phen}


def detect_bistability(p=None, env=None, n=81):
    """Forward/reverse continuation; report the largest hysteresis gap.

    A gap above ~1e-3 means two stable states coexist over a range of chi, so
    the phenotype depends on the environment's *history*, not just its value.
    """
    p = p or get_default_parameters()
    env = copy.deepcopy(env or get_default_environment())
    chi = np.linspace(0.0, 1.0, n)

    def sweep(seq):
        y, out = None, []
        for c in seq:
            env["chi"] = float(c)
            y = run_to_steady_state(p, env, y0=y)
            out.append(phenotype_index(y))
        return np.array(out)

    fwd = sweep(chi)
    rev = sweep(chi[::-1])[::-1]
    gap = np.abs(fwd - rev)
    i = int(np.argmax(gap))

    # Width of the chi-interval over which two states coexist. This is the
    # quantity that decides whether hysteresis is *detectable*: a large
    # phenotype gap confined to a very narrow chi window demands enantiomeric
    # control finer than that window.
    mask = gap > 1e-3
    if mask.any():
        window = (float(chi[mask].min()), float(chi[mask].max()))
        width = window[1] - window[0]
    else:
        window, width = (None, None), 0.0

    return {"chi": chi, "forward": fwd, "reverse": rev,
            "max_gap": float(gap.max()), "chi_at_max_gap": float(chi[i]),
            "bistable": bool(gap.max() > 1e-3),
            "window": window, "window_width": float(width)}


def product_invariance(p=None, env=None, products=(0.30, 0.12, 0.085, 0.06),
                       densities=(0.5, 1.0, 2.0, 4.0)):
    """Test whether phenotype depends only on the product Ltot * (1 - chi).

    In 'racemic' mode the chirality knob enters the model *only* through the
    bindable-ligand density, so enantiomeric excess and ligand density should be
    exactly interchangeable. This measures the residual spread; a spread at
    machine precision confirms the degeneracy.

    This is the model's sharpest falsifiable prediction: an experiment in which
    diluting L-ligand and enantiomerically inverting it are NOT interchangeable
    would refute the parsimonious model and establish chirality-specific biology.
    """
    p = p or get_default_parameters()
    base = copy.deepcopy(env or get_default_environment())

    rows = []
    for prod in products:
        vals = []
        for L in densities:
            chi = 1.0 - prod / L
            if not (0.0 <= chi <= 1.0):
                continue
            e = copy.deepcopy(base)
            e["Ltot"] = float(L)
            e["chi"] = float(chi)
            vals.append((L, chi, phenotype_index(run_to_steady_state(p, e))))
        if vals:
            phen = [v[2] for v in vals]
            rows.append({"product": prod, "points": vals,
                         "spread": float(max(phen) - min(phen))})
    max_spread = max((r["spread"] for r in rows), default=0.0)
    return {"rows": rows, "max_spread": max_spread,
            "invariant": bool(max_spread < 1e-6)}


def detect_rebound(p=None, env=None, n=201, tol=1e-6):
    """Find a non-monotonic 'rebound' in the phenotype-vs-chirality curve.

    A monotonic decline is the signature of pure recognition loss (chirality as
    disguised ligand density). A dip followed by partial recovery is the
    signature of assembly frustration: near-racemic matrices fail to build
    fibrils, but as composition approaches the pure mirror the matrix becomes
    mechanically competent again — even though its epitopes stay unreadable.

    Returns the location and size of the largest rising stretch.
    """
    p = p or get_default_parameters()
    env = env or get_default_environment()
    chi = np.linspace(0.0, 1.0, n)
    phen = _phen_curve(p, env, chi)

    d = np.diff(phen)
    rising = np.where(d > tol)[0]
    if rising.size == 0:
        return {"chi": chi, "phen": phen, "rebound": False,
                "rise": 0.0, "chi_min": None, "chi_max": None}

    i0, i1 = int(rising.min()), int(rising.max()) + 1
    return {"chi": chi, "phen": phen, "rebound": True,
            "rise": float(phen[i1] - phen[i0]),
            "chi_min": float(chi[i0]), "chi_max": float(chi[i1]),
            "phen_min": float(phen[i0]), "phen_max": float(phen[i1])}


def decoupled_scan(p=None, env=None, n=41, chi_adhesive=0.0):
    """Phenotype vs *structural* chirality with adhesive chirality held fixed.

    Models the cleanest discriminating experiment: a D/L collagen gel
    functionalised with a constant density of natural L-RGD, so integrin
    engagement is clamped while the mechanical environment is inverted. Under
    pure recognition loss this curve is flat; under assembly frustration it is
    U-shaped and symmetric about chi = 0.5.
    """
    p = p or get_default_parameters()
    env = copy.deepcopy(env or get_default_environment())
    env["chi"] = float(chi_adhesive)

    chi_s = np.linspace(0.0, 1.0, n)
    phen = []
    for c in chi_s:
        env["chi_struct"] = float(c)
        phen.append(phenotype_index(run_to_steady_state(p, env)))
    phen = np.array(phen)

    mid = phen[len(phen) // 2]
    return {"chi_struct": chi_s, "phen": phen,
            "phen_pure_L": float(phen[0]), "phen_racemic": float(mid),
            "phen_pure_D": float(phen[-1]),
            "dip": float(phen[0] - mid),
            "symmetric": bool(abs(phen[0] - phen[-1]) < 1e-6)}


def design_sensitivity(stiffnesses=(3, 5, 8, 10, 15, 20, 30, 50, 100),
                       p=None, env=None):
    """Find the base stiffness at which the frustration dip is easiest to see.

    The contrast is largest when the frustrated midpoint falls *below* the
    cell's mechanosensing threshold while the pure compositions stay above it.
    Picking the wrong gel stiffness can shrink the effect by an order of
    magnitude, so this is a practical experimental-design output.
    """
    p = p or get_default_parameters()
    base = copy.deepcopy(env or get_default_environment())
    base["assembly_mode"] = "frustrated"
    base["chi"] = 0.0

    rows = []
    for E in stiffnesses:
        e = copy.deepcopy(base)
        e["Estiff"] = float(E)
        e["chi_struct"] = 0.0
        p0 = phenotype_index(run_to_steady_state(p, e))
        e["chi_struct"] = 0.5
        p5 = phenotype_index(run_to_steady_state(p, e))
        rows.append({"stiffness": float(E), "phen_pure": p0,
                     "phen_racemic": p5, "dip": p0 - p5})
    best = max(rows, key=lambda r: r["dip"])
    return {"rows": rows, "best_stiffness": best["stiffness"],
            "best_dip": best["dip"]}


def competence_halving_point(p=None, mode="frustrated_kinetic", n=4001):
    """D fraction at which mechanical competence falls to half.

    This is a *cell-free* observable: it can be read off gel rheology across a
    D/L composition series, with no cells involved. It pins down the capping
    kinetics (lambda_c, cap_efficiency) that everything else depends on.
    """
    from .assembly import assembly_competence

    p = p or get_default_parameters()
    chi = np.linspace(0.0, 0.5, n)
    comp = np.array([assembly_competence(
        {"assembly_mode": mode, "chi_struct": float(c)}, p) for c in chi])
    if comp.min() > 0.5 or comp.max() < 0.5:
        return None
    return float(np.interp(0.5, comp[::-1], chi[::-1]))


def kinetics_sensitivity(p=None, lambda_values=(3, 5, 10, 20, 50, 100),
                         cap_values=(1.0, 0.3)):
    """How the predictions depend on the unknown capping kinetics.

    lambda_c (fibril length needed to bear load) and cap_efficiency are not
    measured here, so any prediction that swings wildly across their plausible
    range should not be the basis of an experiment. This quantifies which
    predictions are robust and which are not.
    """
    from .parameters import get_default_environment

    p = p or get_default_parameters()
    env = get_default_environment()
    env["assembly_mode"] = "frustrated_kinetic"

    rows = []
    for lam in lambda_values:
        for cap in cap_values:
            q = dict(p)
            q["lambda_c"] = float(lam)
            q["cap_efficiency"] = float(cap)
            reb = detect_rebound(p=q, env=env, n=401)
            dec = decoupled_scan(p=q, env=env, n=41)
            rows.append({"lambda_c": float(lam), "cap_efficiency": float(cap),
                         "rebound": reb["rise"], "dip": dec["dip"],
                         "halving_chi": competence_halving_point(q)})

    reb_vals = [r["rebound"] for r in rows]
    dip_vals = [r["dip"] for r in rows]
    return {"rows": rows,
            "rebound_range": (min(reb_vals), max(reb_vals)),
            "dip_range": (min(dip_vals), max(dip_vals))}


ASSEMBLY_MODES = ("self_sorting", "frustrated_kinetic", "majority_rules",
                  "racemic_enhanced")


def discriminate_assembly_modes(p=None, stiffnesses=(0.5, 1, 2, 3, 5, 8, 15, 30, 60),
                                modes=ASSEMBLY_MODES, n=31):
    """Find the gel stiffness that best separates the competing assembly modes.

    The three hypotheses predict *opposite signs* for the clamped-adhesion
    experiment: flat (self-sorting), down (frustration), up (racemic
    enhancement). Sign is a far more robust readout than magnitude — but only
    if the base stiffness leaves headroom in both directions.

    That is the catch: detecting a DIP wants a stiff starting point (so
    softening hurts), while detecting a PEAK wants a soft one (so stiffening
    helps, rather than pushing an already-saturated cell). Optimising for one
    hypothesis can make another invisible. This maximises the minimum pairwise
    separation instead.
    """
    from .parameters import get_default_environment

    p = p or get_default_parameters()
    rows = []
    for E in stiffnesses:
        effects = {}
        for m in modes:
            env = get_default_environment()
            env["assembly_mode"] = m
            env["Estiff"] = float(E)
            d = decoupled_scan(p=p, env=env, n=n)
            # signed: negative = dip, positive = peak
            effects[m] = d["phen_racemic"] - d["phen_pure_L"]
        vals = list(effects.values())
        sep = min(abs(vals[i] - vals[j])
                  for i in range(len(vals)) for j in range(i + 1, len(vals)))
        rows.append({"stiffness": float(E), "effects": effects,
                     "min_separation": float(sep)})
    best = max(rows, key=lambda r: r["min_separation"])
    return {"rows": rows, "best_stiffness": best["stiffness"],
            "best_separation": best["min_separation"]}


#: What each model state corresponds to at the bench. The model's own
#: `phenotype_index` is a composite with weights I chose; experimentalists
#: measure these individually, so predictions should be stated in them.
OBSERVABLE_ASSAYS = {
    "Cint": "focal adhesion size/number (paxillin or vinculin immunostaining)",
    "FAKp": "phospho-FAK Y397 (western or immunofluorescence)",
    "Rho": "RhoA activity (GTP pulldown or FRET biosensor)",
    "Tension": "traction force microscopy; stress-fibre density",
    "YAPn": "YAP nuclear:cytoplasmic ratio  <-- the standard readout",
    "Col": "COL1A1 expression / procollagen staining",
}


def observable_discrimination(p=None, stiffness=2.0, chi_racemic=0.5,
                              modes=ASSEMBLY_MODES):
    """Discriminating power of each individual observable, not the composite.

    Every headline number in this project is expressed in `phenotype_index`, a
    composite whose weights were chosen by hand. This re-expresses the
    three-way mode discrimination in each measurable state variable.

    The important output is the sign-consistency flag. If every observable
    gives the same sign for a given mode, then ANY non-negative weighting of
    them gives that sign too — so the qualitative conclusion is independent of
    the composite's weights by construction, not by lucky choice. That is a
    stronger guarantee than sampling random weights.
    """
    from .parameters import get_default_environment
    from .model import run_to_steady_state, STATE_NAMES, phenotype_index

    p = p or get_default_parameters()

    def steady(mode, chi_struct):
        env = get_default_environment()
        env["assembly_mode"] = mode
        env["Estiff"] = float(stiffness)
        env["chi"] = 0.0                 # adhesion clamped
        env["chi_struct"] = float(chi_struct)
        return run_to_steady_state(p, env)

    base = {m: steady(m, 0.0) for m in modes}
    rac = {m: steady(m, chi_racemic) for m in modes}

    rows = []
    for i, name in enumerate(STATE_NAMES):
        effects = {m: float(rac[m][i] - base[m][i]) for m in modes}
        vals = list(effects.values())
        gap = min(abs(vals[a] - vals[b])
                  for a in range(len(vals)) for b in range(a + 1, len(vals)))
        rows.append({"observable": name, "effects": effects,
                     "min_separation": float(gap),
                     "assay": OBSERVABLE_ASSAYS.get(name, "")})

    comp = {m: float(phenotype_index(rac[m]) - phenotype_index(base[m]))
            for m in modes}

    # sign consistency: does every observable agree on the sign per mode?
    def sign(x, tol=1e-6):
        return 0 if abs(x) < tol else (1 if x > 0 else -1)

    consistent = {}
    for m in modes:
        signs = {sign(r["effects"][m]) for r in rows}
        consistent[m] = len(signs) == 1

    rows.sort(key=lambda r: -r["min_separation"])
    return {"rows": rows, "composite": comp,
            "sign_consistent": consistent,
            "weight_independent": all(consistent.values()),
            "best": rows[0]["observable"], "worst": rows[-1]["observable"]}


def two_point_protocol(p=None, compositions=(0.10, 0.50),
                       stiffnesses=(1, 2, 3, 5, 8, 15), modes=ASSEMBLY_MODES):
    """Classify the assembly hypotheses from TWO compositions, not one.

    An earlier version of this analysis compared racemic against pure and
    searched for the stiffness maximising the spread. That was the wrong
    question. The four hypotheses differ in the *shape* of competence versus
    composition, and a single contrast throws the shape away:

      self_sorting        flat everywhere
      frustrated_kinetic  collapses at a few percent D — signal at the EXTREMES
      majority_rules      flat until the majority is lost — signal at the CENTRE
      racemic_enhanced    rises toward racemic

    At racemic, frustration and majority-rules are both collapsed and nearly
    indistinguishable (separation ~0.008). At low D fraction, majority-rules is
    nearly indistinguishable from self-sorting. Neither point works alone; the
    pair does, because each mode has a distinct (low-chi, racemic) signature.

    Returns per-mode signatures and the worst-case pairwise distance, which is
    the quantity a real experiment has to beat with its error bars.
    """
    from .parameters import get_default_environment
    from .model import run_to_steady_state, phenotype_index
    import itertools

    p = p or get_default_parameters()

    def phen(mode, E, chi_struct):
        env = get_default_environment()
        env["assembly_mode"] = mode
        env["Estiff"] = float(E)
        env["chi"] = 0.0                 # adhesion clamped
        env["chi_struct"] = float(chi_struct)
        return phenotype_index(run_to_steady_state(p, env))

    rows = []
    for E in stiffnesses:
        base = {m: phen(m, E, 0.0) for m in modes}
        sig = {m: tuple(phen(m, E, c) - base[m] for c in compositions)
               for m in modes}
        worst, pair = None, None
        for a, b in itertools.combinations(modes, 2):
            d = float(np.linalg.norm(np.array(sig[a]) - np.array(sig[b])))
            if worst is None or d < worst:
                worst, pair = d, (a, b)
        rows.append({"stiffness": float(E), "signatures": sig,
                     "worst_separation": worst, "limiting_pair": pair})

    best = max(rows, key=lambda r: r["worst_separation"])
    return {"rows": rows, "compositions": tuple(compositions),
            "best_stiffness": best["stiffness"],
            "best_separation": best["worst_separation"],
            "best_signatures": best["signatures"],
            "limiting_pair": best["limiting_pair"]}


def observable_ranking_two_point(p=None, stiffness=2.0,
                                 compositions=(0.10, 0.50), modes=ASSEMBLY_MODES):
    """Rank readouts under the two-point protocol, which is the real design.

    Readout choice and protocol choice are coupled, and getting the protocol
    wrong inverts the answer. Under a single racemic-vs-pure contrast with four
    hypotheses in play, the limiting pair is frustration vs majority-rules, and
    those two differ most in the UPSTREAM nodes — the saturating steps below
    compress the difference — so pFAK ranks first and nuclear YAP last. Under
    the two-point protocol the ranking reverts to adhesion and transcriptional
    readouts. The inversion is an artefact of the inadequate design, not a fact
    about the biology.
    """
    from .parameters import get_default_environment
    from .model import run_to_steady_state, STATE_NAMES
    import itertools

    p = p or get_default_parameters()

    def state(mode, chi_struct):
        env = get_default_environment()
        env["assembly_mode"] = mode
        env["Estiff"] = float(stiffness)
        env["chi"] = 0.0
        env["chi_struct"] = float(chi_struct)
        return run_to_steady_state(p, env)

    base = {m: state(m, 0.0) for m in modes}
    pts = {m: {c: state(m, c) for c in compositions} for m in modes}

    rows = []
    for i, name in enumerate(STATE_NAMES):
        sig = {m: np.array([pts[m][c][i] - base[m][i] for c in compositions])
               for m in modes}
        worst = min(float(np.linalg.norm(sig[a] - sig[b]))
                    for a, b in itertools.combinations(modes, 2))
        rows.append({"observable": name, "signatures": sig,
                     "worst_separation": worst,
                     "assay": OBSERVABLE_ASSAYS.get(name, "")})
    rows.sort(key=lambda r: -r["worst_separation"])
    return {"rows": rows, "best": rows[0]["observable"],
            "worst": rows[-1]["observable"],
            "compositions": tuple(compositions), "stiffness": float(stiffness)}


def chirality_stiffness_map(stiffnesses=(1.0, 3.0, 10.0, 30.0, 100.0),
                            p=None, env=None):
    """How chi* and the size of the chirality effect vary with stiffness."""
    p = p or get_default_parameters()
    base = copy.deepcopy(env or get_default_environment())
    out = []
    for E in stiffnesses:
        e = copy.deepcopy(base)
        e["Estiff"] = float(E)
        r = find_chi_star(p, e)
        out.append({"stiffness": float(E), "chi_star": r["chi_star"],
                    "phen_native": float(r["phen"][0]), "range": r["range"]})
    return out
