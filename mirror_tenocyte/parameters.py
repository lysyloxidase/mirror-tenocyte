"""Default parameters and environment for the mirror-tenocyte model.

All state variables are normalised activities/occupancies in [0, 1]:

    Cint     fraction of engaged / clustered integrins (molecular clutch)
    FAKp     active (Y397-phosphorylated) FAK fraction
    Rho      active RhoA/ROCK signal (contractility driver)
    Tension  cytoskeletal / stress-fibre tension level
    YAPn     nuclear YAP/TAZ fraction
    Col      collagen-I production proxy (Col1a1 program)

Rate constants are in units of 1/time (arbitrary time units; steady states are
the object of study, so absolute timescale is not calibrated). Values are
order-of-magnitude, chosen so that the natural all-L environment (chi = 0) on a
tendon-like stiffness (~ tens of kPa) sits in the "high adhesion / nuclear YAP /
high collagen" tenogenic regime, and are *not* fitted to data.

Biological motivation (qualitative):
  * Molecular clutch / stiffness sensing:
      Chan & Odde 2008; Elosegui-Artola et al. 2016 Nat Cell Biol.
  * FAK Y397 as the adhesion-proximal switch:
      Mitra, Hanks & Schlaepfer 2005 Nat Rev Mol Cell Biol.
  * RhoA/ROCK -> actomyosin tension -> YAP/TAZ:
      Dupont et al. 2011 Nature; Totaro, Panciera & Piccolo 2018 Nat Cell Biol.
  * YAP/TAZ as mechanosensitive transcriptional co-activators (CTGF/CCN2, CYR61)
    and tenocyte matrix output / tenogenic maintenance:
      Jones et al. 2019; Kaneko et al.; Piccolo, Dupont & Cordenonsi 2014.
"""

from __future__ import annotations


def get_default_parameters() -> dict:
    """Return a fresh dict of default kinetic + shape parameters."""
    return {
        # ---- integrin engagement / molecular clutch --------------------
        "kon": 4.0,        # engagement rate toward clutch set-point
        "koff": 1.0,       # disengagement rate
        "K_adh": 0.55,     # half-saturation ligand drive for engagement
        "n_adh": 3.0,      # cooperativity of integrin clustering (Hill)
        "a_reinforce": 2.0,  # tension -> adhesion positive feedback (talin/clutch)
        "KE_adh": 8.0,     # stiffness (kPa) for half-max clutch engagement

        # ---- FAK -------------------------------------------------------
        "kf": 3.0,
        "kdf": 1.0,

        # ---- RhoA/ROCK -------------------------------------------------
        "kr": 3.0,
        "kdr": 1.0,

        # ---- cytoskeletal tension -------------------------------------
        "kt": 3.0,
        "kdt": 1.0,
        "KE_ten": 6.0,     # stiffness (kPa) needed to "resist" and build tension
        "hE_adh": 1.0,     # Hill exponent for the stiffness->adhesion term
        "hE_ten": 1.0,     # Hill exponent for the stiffness->tension term

        # ---- YAP/TAZ nuclear translocation ----------------------------
        "ky": 3.0,
        "kdy": 1.0,
        "K_yap": 0.5,      # half-max signal for nuclear YAP
        "n_yap": 3.0,      # switch sharpness (bistability when combined w/ a_yap)
        "w_ten": 0.7,      # weight of tension in the YAP-activating signal
        "w_fak": 0.3,      # weight of FAK in the YAP-activating signal
        "a_yap": 0.6,      # YAP self-reinforcement (nuclear retention feedback)

        # ---- collagen / tenogenic output ------------------------------
        "kc": 1.0,         # keep Col ~ tracking YAPn (avoid saturating the readout)
        "kdc": 1.0,

        # ---- chirality mapping (see chirality.py) ---------------------
        "beta_chi": 4.0,   # affinity-decay steepness for 'affinity' mode
        "chi_specific": 0.0,  # optional D-ECM-specific signaling gain (hypothesis knob)

        # ---- supramolecular assembly (see assembly.py) ----------------
        "assembly_severity": 0.9,  # how completely a 50:50 mixture fails to assemble
        "assembly_exponent": 1.0,  # network stiffness ~ competence ** exponent

        # 'frustrated_kinetic': growth-poisoning (capping) kinetics.
        # Signal sits at the EXTREMES — competence halves at a few percent D.
        "cap_efficiency": 1.0,     # P(a mismatched encounter terminates growth)
        "lambda_c": 20.0,          # triple helices needed for a fibril to bear load
        "percolation_sharpness": 2.0,

        # 'majority_rules': 1D Ising chiral copolymerisation. The minority
        # enantiomer is ABSORBED into the majority screw sense rather than
        # excluded, so the signal sits at the CENTRE of the composition range.
        # Values measured for C3-symmetric disks (JACS 2005, 127, 5490),
        # NOT for collagen — the form is principled, the numbers are borrowed.
        "mismatch_penalty": 0.94,          # kJ/mol
        "helix_reversal_penalty": 7.8,     # kJ/mol
        "majority_severity": 0.7,          # competence lost at zero net helicity

        # 'racemic_enhanced': fold stiffness gain of the racemic coassembly.
        # 4.0 is the measured MAX1/DMAX1 value (800 Pa vs 200 Pa) — but that is
        # a beta-hairpin, and its mechanism does not transfer to collagen; see
        # assembly.heterochiral_helix_probability.
        "racemic_enhancement": 4.0,

        # Level-1 hierarchy: cost of a D residue inside an L triple helix.
        # 7.87 kcal/mol from MD (J. Phys. Chem. B 2009, 113, 8983).
        "ddG_heterochiral": 7.87,
        "temperature_K": 310.0,
    }


#: Named parameter regimes. The default ('parsimonious') is deliberately the
#: weakest-assumption setting and turns out to be MONOSTABLE — it shows a sharp
#: but reversible, history-independent transition. 'strong_feedback' raises the
#: clutch cooperativity and the tension->adhesion reinforcement into the regime
#: where the network becomes bistable and shows hysteresis.
#:
#: Whether real tenocytes sit in the monostable or the bistable regime is an
#: EMPIRICAL question, not something this model can settle. Presets exist so
#: that both hypotheses can be simulated and told apart experimentally, rather
#: than one being silently baked into the defaults.
PRESETS = {
    "parsimonious": {},  # defaults as-is
    "strong_feedback": {
        "a_reinforce": 5.0,
        "n_adh": 4.0,
        "koff": 0.7,
        "K_adh": 0.6,
    },
    # Stiffness terms DERIVED from the calibrated stochastic motor-clutch
    # (clutch.py) rather than assumed. The clutch response is fitted by Hill
    # functions with R^2 = 0.992 (engaged clutches -> adhesion) and 0.999
    # (traction -> tension), which vindicates the coarse-grained functional
    # form and supplies its constants.
    #
    # The adhesion constant was nearly right (8 assumed vs 5.6 derived), but
    # the TENSION constant was ~3x too low (6 assumed vs 17.0 derived): the
    # baseline model builds cytoskeletal tension far too easily on soft
    # matrix. Everything downstream of tension is affected.
    "clutch_calibrated": {
        "KE_adh": 5.62,
        "hE_adh": 0.869,
        "KE_ten": 16.97,
        "hE_ten": 1.208,
    },
}


def get_preset(name: str = "parsimonious") -> dict:
    """Return default parameters with a named preset's overrides applied."""
    if name not in PRESETS:
        raise ValueError(f"unknown preset {name!r}; choose from {sorted(PRESETS)}")
    p = get_default_parameters()
    p.update(PRESETS[name])
    return p


def get_default_environment() -> dict:
    """Return a fresh dict describing the pericellular environment.

    chir_mode selects how chirality couples to integrin recognition:
      'racemic'  : fraction chi of ligands are D and non-bindable, so the
                   effective bindable ligand density scales as (1 - chi).
                   Under this (parsimonious) view chi = 1 is identical to a
                   non-adhesive substrate.
      'affinity' : every ligand's affinity degrades smoothly with chi as
                   exp(-beta_chi * chi); models a graded loss of productive
                   binding rather than an all-or-none population swap.

    Assembly coupling (see assembly.py) is selected by `assembly_mode`:
      'none'         : stiffness is whatever Estiff says, independent of chi
                       (the baseline used for the original analysis).
      'self_sorting' : D and L segregate into separate perfect fibrils, so
                       mechanics stay normal at every chi. Reduces exactly to
                       the parsimonious model.
      'frustrated'   : mismatched monomers poison fibril growth, so effective
                       stiffness dips at intermediate chi and recovers at both
                       pure compositions.

    `chi_struct` sets the chirality of the *structural* matrix separately from
    the adhesive ligand; None means it follows `chi`.
    """
    return {
        "chi": 0.0,          # environmental chirality: 0 = all-L, 1 = all-D (mirror)
        "Ltot": 1.0,         # total ECM ligand density (normalised)
        "Estiff": 30.0,      # substrate/matrix stiffness in kPa (tendon-like)
        "chir_mode": "racemic",
        "assembly_mode": "none",
        "chi_struct": None,  # None -> follows chi (one matrix, one composition)
    }
