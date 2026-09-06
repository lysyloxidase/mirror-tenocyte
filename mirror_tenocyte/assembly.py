"""Supramolecular assembly of a mixed-chirality matrix — the mechanical channel.

This module exists to answer the central weakness found in the baseline model:
if chirality enters only through integrin recognition, then it is not an
independent variable at all — phenotype collapses onto the single product
``Ltot * (1 - chi)`` (see analysis.product_invariance). Chirality only becomes a
genuine axis if D/L composition changes something *other* than bindable ligand
density. Fibrillogenesis is that something.

The physics being encoded
-------------------------
Collagen I is a right-handed superhelix of three left-handed, PPII-like chains.
Its exact mirror image is a left-handed superhelix of right-handed chains. Two
facts follow, and they pull in opposite directions:

1. A *pure* mirror matrix assembles perfectly well. Enantiomeric structures are
   isoenergetic, and elastic moduli are parity-invariant, so an all-D fibril has
   the same stiffness as an all-L fibril. (Empirical precedent: chemically
   synthesised D-proteins fold into mirror structures with identical
   thermodynamics and activity — Milton et al. 1992, Science 256:1445, D-HIV-1
   protease.) A pure mirror tendon would be mechanically normal and merely twist
   the other way.

2. A *mixed* matrix may not assemble at all. A D monomer cannot pack into a
   growing L triple helix; mismatched subunits cap growing fibrils. Assembly is
   therefore expected to fail worst at 50:50 and to recover at both pure
   compositions.

That makes mechanical competence **non-monotonic and symmetric about chi = 0.5**,
whereas bindable ligand density is monotonic in chi. The two channels are no
longer expressible as one variable, and the degeneracy breaks.

FOUR competing assembly hypotheses are provided, because they are
experimentally distinguishable and predict different SHAPES of competence versus
composition — not merely different magnitudes:

``self_sorting``
    D and L segregate into separate, internally perfect fibrils (conglomerate
    behaviour). Every fibril is mechanically normal, so competence stays 1 at
    all chi; only the *bindable* fraction falls. This reduces **exactly** to the
    parsimonious model — so the Ltot*(1-chi) degeneracy is the fingerprint of
    self-sorting.

``frustrated``
    Mismatched monomers poison growth (racemic-compound / capping behaviour).
    Competence dips at intermediate chi, dragging effective stiffness down with
    it, and the degeneracy is violated. Uses an ASSUMED symmetric quadratic.

``frustrated_kinetic``
    The same physical hypothesis, but with the competence curve DERIVED from
    growth-poisoning kinetics rather than assumed. See below — it turns out to
    be dramatically more sensitive to small enantiomeric contamination than the
    quadratic guess, which changes what experiment to run.

``racemic_enhanced``
    The opposite outcome: L and D coassemble into a heterochiral fibril that is
    STIFFER than either pure form. Racemic MAX1/DMAX1 beta-hairpin hydrogels reach G' ~ 800 Pa
    versus ~200 Pa for either enantiopure gel — a 4-fold enhancement — via a
    "rippled" beta-sheet with alternating L and D strands, exactly as Pauling
    and Corey predicted in 1953. The reported Job plot has a clean symmetric
    maximum at 0.5 mole fraction, which is the shape used here.
        Schneider/Pochan et al., "Enhanced mechanical rigidity of hydrogels
        formed from enantiomeric peptide assemblies" (PMC3202337);
        "Molecular, Local, and Network-Level Basis for the Enhanced Stiffness
        of Hydrogel Networks Formed from Coassembled Racemic Peptides",
        ACS Cent. Sci. 2017, 3, 586.

    IMPORTANT: this mode was added after checking the literature, and it
    inverts the sign of the model's headline prediction. A later version then
    argued it should LEAD the prior for level 2, on the grounds that lateral
    fibril packing is "solid-state" and Wallach's rule favours heterochiral
    packing there. That argument was WITHDRAWN — see evidence.py. A collagen
    fibril is a hydrated quasi-hexagonal array of ~300 nm rods, not a molecular
    crystal, and the transfer was never checked.

    All four outcomes are empirically realised in different systems, none of
    them a triple helix. `racemic_enhanced` is UNRESOLVED for level 2, exactly
    like the other three.

Deriving the competence curve
-----------------------------
Treat fibril growth as monomer addition that terminates when a mismatched
monomer caps the growing end. For an L fibril in a matrix with D fraction chi,
correct additions occur at rate proportional to (1 - chi) and capping events at
rate proportional to chi * p_cap, so the mean length before termination is

    lambda_L(chi) = (1 - chi) / (chi * p_cap)

and by mirror symmetry lambda_D(chi) = chi / ((1 - chi) * p_cap). Both
populations grow at once, each poisoned by the other's monomers.

A fibril only bears load once it is long enough to span between crosslinks, so
competence follows a percolation-like threshold g(lambda) with critical length
lambda_c. Weighting each population by its mass fraction:

    A(chi) = (1 - chi) * g(lambda_L) + chi * g(lambda_D)

This still satisfies A(0) = A(1) = 1 and dips at chi = 0.5, but the *shape* is
quite different from the quadratic guess: capping is highly efficient, so
competence collapses at low contamination. With lambda_c = 100 monomers,
competence halves at ~1% D, versus ~17% for the assumed quadratic. That is the
classic chiral-poisoning behaviour seen in crystal growth, where trace amounts
of the opposing enantiomer arrest growth entirely (Addadi & Lahav).
"""

from __future__ import annotations


def heterochiral_helix_probability(p: dict, n_residues: float = 30.0) -> float:
    """Boltzmann probability of incorporating a D chain into an L triple helix.

    This settles the FIRST level of the collagen assembly hierarchy, and it
    settles it without any free parameter.

    Collagen assembles in stages: chain -> triple helix -> microfibril ->
    fibril. Chirality mismatch can in principle act at any of them, and the
    levels behave very differently:

      Level 1 (chain -> triple helix). Three chains must interdigitate with Gly
        at the helix core. Molecular dynamics puts the cost of a single D-Ala
        substitution at ~7.87 kcal/mol — comparable to the Gly->Ala mutations
        that cause osteogenesis imperfecta — producing a kink that breaks the
        helical register (J. Phys. Chem. B 2009, 113, 8983; PMID 19518060).
        At 37 C that is exp(-7.87/RT) ~ 3e-6 per residue, and a triple-helical
        domain has ~1000 residues. The probability underflows after a few dozen.
        => heterochiral triple helices DO NOT FORM. Chains sort perfectly.

        This is the one claim in the project with EXPERIMENTAL support, not just
        a calculation. Host-guest peptides carrying a single D-Asp residue
        (Shah et al., Biopolymers 1999, 49, 297; PMID 10079768) fail to form a
        triple helix at all in phosphate-buffered saline — one wrong stereocentre
        is enough. Two qualifications come with it, and neither rescues an all-D
        chain:
          * Solvent-dependent. In 67% aqueous ethylene glycol, a strongly
            stabilising solvent, helices do form — with melting temperatures
            more than 30 C lower. The exclusion is a physiological-buffer
            statement, which is the condition the proposed experiment uses.
          * Sequence-context-dependent. Mixed D-/L-Asp peptides formed
            heterotrimers for the Gly-Asp-Ala guest triplet but not for
            Gly-Asp-Hyp. So a chain carrying ONE D residue can sometimes join a
            mixed trimer.
        Both concern single substitutions. An all-D 36-mer carries ~24
        stereocentres, so the argument against heterochiral helices is stronger
        than the MD number alone implies, not weaker.

      Level 2 (triple helix -> fibril). Homochiral helices, now already sorted,
        pack laterally via charge pairs and hydrophobic contacts in a staggered
        array. A D helix is the mirror image of an L helix: charges are achiral
        and could still pair, but the stagger geometry is chiral. THIS is where
        the open question actually lives, and it is what `assembly_mode` models.

    The practical consequence is a warning about transferring literature: the
    MAX1/DMAX1 racemic stiffening (see `racemic_enhanced`) works through a
    "rippled" beta-sheet of ALTERNATING L and D strands — that is a level-1
    mixing mechanism. Collagen forbids level-1 mixing, so the strongest
    empirical support for racemic enhancement does not carry over. Enhancement
    via level-2 helix co-packing is not excluded, but it is unsupported.

    Also note the unit change this forces: at level 2 the capping "monomer" of
    `frustrated_kinetic` is a whole triple helix (~300 nm rod), not a residue,
    so `lambda_c` must be read in triple-helix counts.
    """
    import math
    ddG = float(p.get("ddG_heterochiral", 7.87))   # kcal/mol per D residue
    T = float(p.get("temperature_K", 310.0))
    RT = 1.987e-3 * T
    return math.exp(-min(n_residues, 200.0) * ddG / RT)


def assembly_competence(env: dict, p: dict) -> float:
    """Fraction of the structural matrix forming mechanically competent fibrils.

    Returns a value in (0, 1]. Depends on the *structural* chirality
    ``chi_struct``, which defaults to following the adhesive chirality ``chi``
    but can be set independently — decoupling them is the basis of the model's
    cleanest proposed experiment (see analysis.decoupled_scan).
    """
    mode = env.get("assembly_mode", "none")
    if mode == "none":
        return 1.0

    chi_s = structural_chirality(env)
    if mode == "self_sorting":
        # Separate D-fibrils and L-fibrils, each internally perfect. Mirror
        # fibrils are mechanically equivalent to natural ones, so the composite
        # network is mechanically normal at every composition.
        return 1.0
    if mode == "frustrated":
        # Symmetric about chi = 0.5, unity at both pure compositions. The
        # 4*chi*(1-chi) factor is the simplest such form; `assembly_severity`
        # sets how completely a racemic mixture fails to assemble.
        sev = float(p.get("assembly_severity", 0.9))
        comp = 1.0 - sev * 4.0 * chi_s * (1.0 - chi_s)
        return max(comp, 1e-3)   # a failed gel is soft, not massless
    if mode == "frustrated_kinetic":
        return max(_competence_from_capping(chi_s, p), 1e-3)
    if mode == "majority_rules":
        return max(_competence_from_majority_rules(chi_s, p), 1e-3)
    if mode == "racemic_enhanced":
        # Heterochiral coassembly that is STIFFER than either pure form.
        # Symmetric parabola peaking at chi = 0.5, matching the Job-plot shape
        # reported for MAX1/DMAX1 (see module docstring).
        gain = float(p.get("racemic_enhancement", 4.0))
        return 1.0 + (gain - 1.0) * 4.0 * chi_s * (1.0 - chi_s)
    raise ValueError(f"unknown assembly_mode: {mode!r}")


def net_helicity(chi: float, p: dict) -> float:
    """Net screw-sense bias of a mixed-chirality supramolecular polymer.

    Chiral copolymerisation has a mature quantitative theory that this project
    initially missed: a **one-dimensional two-component Ising model**, exact in
    the long-chain limit (van Gestel/Palmans/Meijer; see
    Nat. Commun. 2011, 2, 509 for the general treatment). Two energies matter:

      HRP  helix reversal penalty — cost of a domain wall between opposite
           screw senses
      MMP  mismatch penalty — cost of one monomer sitting in a helix of its
           non-preferred sense

    Measured for a C3-symmetric disk system: MMP = 0.94 kJ/mol and
    HRP = 7.8 kJ/mol, i.e. reversal costs ~8x more than mismatch
    (J. Am. Chem. Soc. 2005, 127, 5490; 2009, 131, 18784).

    That inequality is the whole point. When reversal is dearer than mismatch,
    the polymer does NOT expel or cap the minority enantiomer and does NOT
    break into domains — it **absorbs the minority into the majority's screw
    sense**. This is the "majority-rules" effect, and it is a fourth outcome
    distinct from self-sorting, capping and racemic coassembly.

    With field H = MMP*ee and coupling J = HRP/2, the exact 1D Ising
    magnetisation gives the net helicity

        m(ee) = sinh(bH) / sqrt(sinh^2(bH) + exp(-4bJ))

    which amplifies by exp(2bJ) ~ 21 at 37 C, so helicity saturates by
    ee ~ 0.27. The consequence for this project is a signature localised at the
    CENTRE of the composition range rather than at its extremes.

    Caveat, and it is a large one: these constants come from small-molecule
    disks in apolar solvent, not from collagen. Whether a collagen fibril's
    supertwist behaves as a cooperative 1D chain with comparable energies is
    unknown. The functional form is principled; the numbers are borrowed.
    """
    import math
    R = 8.314e-3          # kJ/(mol K)
    T = float(p.get("temperature_K", 310.0))
    RT = R * T
    mmp = float(p.get("mismatch_penalty", 0.94))     # kJ/mol
    hrp = float(p.get("helix_reversal_penalty", 7.8))  # kJ/mol

    ee = abs(1.0 - 2.0 * min(max(chi, 0.0), 1.0))    # enantiomeric excess
    bH = mmp * ee / RT
    bJ = (hrp / 2.0) / RT
    s = math.sinh(bH)
    return s / math.sqrt(s * s + math.exp(-4.0 * bJ))


def _competence_from_majority_rules(chi: float, p: dict) -> float:
    """Competence when the minority enantiomer is absorbed, not excluded.

    Fibrils still form and still grow through the minority monomer, so there is
    no length truncation. What degrades near racemic composition is screw-sense
    coherence: domain walls accumulate as the net helicity falls, and they are
    structural defects. Competence therefore tracks |m| rather than collapsing.
    """
    sev = float(p.get("majority_severity", 0.7))
    return 1.0 - sev * (1.0 - abs(net_helicity(chi, p)))


def _competence_from_capping(chi: float, p: dict) -> float:
    """Mechanical competence derived from growth-poisoning kinetics.

    Both enantiomeric populations grow simultaneously, each capped by the
    other's monomers. Mass-weighted sum of their load-bearing fractions.
    """
    eps = 1e-12
    p_cap = float(p.get("cap_efficiency", 1.0))
    lam_c = float(p.get("lambda_c", 20.0))
    m = float(p.get("percolation_sharpness", 2.0))

    lam_L = (1.0 - chi) / (chi * p_cap + eps)
    lam_D = chi / ((1.0 - chi) * p_cap + eps)

    def g(lam):
        lam = max(lam, 0.0)
        return lam ** m / (lam_c ** m + lam ** m)

    return (1.0 - chi) * g(lam_L) + chi * g(lam_D)


def structural_chirality(env: dict) -> float:
    """Chirality of the bulk/structural matrix.

    Falls back to the adhesive chirality ``chi`` when not set, which is the
    physically coupled case (one matrix, one composition). Setting it separately
    represents a designed experiment: e.g. a D/L collagen gel functionalised
    with a *constant* density of natural L-RGD, so adhesion is held fixed while
    the mechanical environment is inverted.
    """
    chi_s = env.get("chi_struct", None)
    if chi_s is None:
        chi_s = env.get("chi", 0.0)
    return min(max(float(chi_s), 0.0), 1.0)


def effective_stiffness(env: dict, p: dict) -> float:
    """Substrate stiffness actually felt by the cell, after assembly quality.

    Network stiffness is taken to scale with competence raised to
    ``assembly_exponent`` (default 1). Fibre networks generally stiffen
    superlinearly with connectivity, so an exponent > 1 is defensible; it is
    left at 1 to avoid smuggling in extra nonlinearity.
    """
    base = float(env["Estiff"])
    comp = assembly_competence(env, p)
    m = float(p.get("assembly_exponent", 1.0))
    return base * (comp ** m)
