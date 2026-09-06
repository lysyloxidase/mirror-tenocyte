"""Published outcomes for enantiomeric (D/L) assembly — a weak, fragile survey.

Nobody has measured what mixed-chirality collagen or collagen-mimetic peptide
does. No synthesis is available to this project either, so this module collects
what other people have measured in *other* systems.

READ THE LIMITS BEFORE THE CONTENT
----------------------------------
This table is held as data and covered by tests, which makes it look more
rigorous than it is. The tests check *structural* consistency — that every
outcome maps to a modelled hypothesis, that every entry carries a provenance
level, that the fragility flag is set. **They cannot check whether the entries
are true.** Specifically:

* **Only 2 of 6 entries were verified from full text.** The rest come from
  search-result summaries or indirect citation in a review. Each entry carries a
  `verification` field; consult it before relying on the entry. The
  phenylalanine numbers in particular (53.1 / 5.8 / 1.8 GPa) come from a summary
  of a paywalled paper. Its abstract, reached through an institutional
  repository, confirms only the DIRECTION of the effect and contains no
  numbers at all — so the magnitude is unsupported by anything readable here.
  Outcome direction and effect magnitude are tracked separately for this
  reason (`verification` vs `fold_verified`).
* **n = 6, and the leading outcome flips if any single entry is
  reclassified** (see `fragility`). Reporting this as percentages would imply a
  precision that does not exist, so `prior()` deliberately returns counts and a
  fragility flag rather than fractions.
* **Zero entries are triple helices.** All are beta-sheets, amphipathic
  peptides, or single amino acids. Two of the three "enhanced" entries are not
  aqueous peptide gels at all — phenylalanine is a solidified amino acid with a
  modulus in GPa, three orders of magnitude from a hydrogel.
* **Publication bias is severe and directional.** "Racemic is 9x stiffer"
  publishes; "we mixed them and nothing happened" does not. The self-sorting
  count is almost certainly too low.

What the survey does NOT license
--------------------------------
An earlier version of this module argued that because solution phase favours
homochiral assembly while solid-state packing favours heterochiral (the Abeta
reconciliation, and Wallach's rule for racemic crystals), collagen's level-2
lateral packing must sit in the heterochiral-favouring regime — and used that to
promote `racemic_enhanced` to the top of the prior.

**That argument is withdrawn.** Wallach's rule concerns molecular crystals of
small molecules; the Abeta observation concerns beta-sheet peptides in solution
versus the solid state. A collagen fibril is a hydrated, quasi-hexagonal array
of ~300 nm rods holding a large water fraction — it is not a molecular crystal,
and whether the rule transfers across that boundary was never checked.

This matters because it is the *same class of error* the project had just
diagnosed when it rejected transferring the MAX1 result to collagen: carrying a
result across a structural boundary without justification. The rigour is also
asymmetric. The level-1 demotion rests on a number and a Boltzmann calculation
(7.87 kcal/mol, underflowing after ~30 residues). The re-promotion rested on an
analogy. **`racemic_enhanced` is therefore UNRESOLVED for level 2** — neither
excluded nor favoured.
"""

from __future__ import annotations

#: Outcome classes, matching analysis.ASSEMBLY_MODES.
ENHANCED = "racemic_enhanced"
SORTING = "self_sorting"
WEAKENED = "frustrated_kinetic"      # any mechanical degradation on mixing
MAJORITY = "majority_rules"

#: How well each entry is sourced. Anything below FULL_TEXT should be treated
#: as a lead to check, not as a data point.
FULL_TEXT = "full text read"
ABSTRACT = "search summary / abstract only"
INDIRECT = "cited indirectly in a review; primary not seen"

SURVEY = [
    {
        "system": "MAX1 / DMAX1 β-hairpin hydrogel",
        "class": "β-sheet fibril, hydrogel",
        "outcome": ENHANCED,
        "fold": 4.0,
        "verification": FULL_TEXT,
        "detail": "G' 800 Pa racemic vs ~200 Pa enantiopure; heterochiral "
                  "'rippled' sheet with nested Val packing",
        "composition_series": True,     # 3:1, 1:1, 1:3 Job plot, symmetric max
        "ref": "Nagy et al. JACS 2011; Schneider/Pochan ACS Cent. Sci. 2017, 3, 586",
    },
    {
        "system": "Phenylalanine D/L racemate",
        "class": "single amino acid, crystalline assembly (GPa, not a hydrogel)",
        "outcome": ENHANCED,
        "fold": 9.0,                    # 53.1 GPa racemic vs 5.8 (D) / 1.8 (L)
        "verification": ABSTRACT,
        "fold_verified": False,
        "detail": "The abstract (institutional repository) confirms only the "
                  "DIRECTION — mixing gives 'nanostructures that are "
                  "mechanically more robust'. The numbers 53.1 / 5.8 / 1.8 GPa "
                  "come from a search-result summary and could not be "
                  "corroborated from any text available here. Also three orders "
                  "of magnitude stiffer than any hydrogel in this survey",
        "composition_series": False,
        "ref": "ACS Nano 2020, 14, 1694",
    },
    {
        "system": "KFE8 amphipathic peptide",
        "class": "β-sheet fibril",
        "outcome": ENHANCED,
        "fold": None,
        "verification": INDIRECT,
        "detail": "reported to coassemble into rippled β-sheets with alternating "
                  "enantiomers; primary source not consulted",
        "composition_series": False,
        "ref": "cited in Chem. Rev. 2021, 121, 13869",
    },
    {
        "system": "Aβ-derived peptides (solution)",
        "class": "β-sheet oligomer, aqueous solution",
        "outcome": SORTING,
        "fold": 1.0,
        "verification": FULL_TEXT,
        "detail": "homochiral pleated sheets dominate; for one peptide 'no "
                  "detectable heterochiral assembly'. Heterochiral tetramers, "
                  "where seen, are two homochiral dimers",
        "composition_series": False,
        "ref": "PMC9258340",
    },
    {
        "system": "β-amyloid-inspired peptide amphiphiles",
        "class": "peptide amphiphile nanofibre",
        "outcome": SORTING,
        "fold": 1.0,
        "verification": ABSTRACT,
        "detail": "chirality-mediated self-sorting into chemically distinct "
                  "nanostructures; enantioselective enzymatic degradation",
        "composition_series": False,
        "ref": "PMID 32970093",
    },
    {
        "system": "Tripeptide Dff added to opposite-enantiomer hFF",
        "class": "tripeptide hydrogel",
        "outcome": WEAKENED,
        "fold": None,
        "verification": ABSTRACT,
        "detail": "coassembled gel reached LOWER elastic modulus and markedly "
                  "slower gelation; same-stereoconfiguration mixing was neutral",
        "composition_series": False,
        "ref": "ChemRxiv 2025-qfj9j / Faraday Discuss. 2025",
    },
]

#: The one entry here that bears directly on collagen, and the only claim in the
#: project with experimental rather than computational support. It concerns
#: LEVEL 1 (does a mismatched chain enter a triple helix), not level 2.
LEVEL1_EXPERIMENTAL = {
    "system": "Host-guest collagen peptides with a single D-Asp residue",
    "framework": "acetyl-(Gly-Pro-Hyp)3-Gly-Xaa-Yaa-(Gly-Pro-Hyp)4-Gly-Gly-amide",
    "finding": "D-Asp PREVENTS triple-helix formation in phosphate-buffered "
               "saline. Helices form only in 67% aqueous ethylene glycol, with "
               "melting temperatures >30 C lower than the L-Asp peptides.",
    "qualification": "Sequence-context-dependent: mixed D-/L-Asp peptides formed "
                     "heterotrimers for Gly-Asp-Ala but NOT for Gly-Asp-Hyp. "
                     "Both observations concern SINGLE substitutions; an all-D "
                     "36-mer carries ~24 stereocentres, so this strengthens the "
                     "case against heterochiral helices rather than weakening it.",
    "verification": ABSTRACT,
    "ref": "Shah et al., Biopolymers 1999, 49, 297 (PMID 10079768)",
}

#: Chiral supramolecular polymers (BTA, C3-disks) show majority-rules cleanly,
#: but the readout is CD, not mechanics, so they are listed separately: they
#: establish that the mode is real without contributing a modulus data point.
MAJORITY_RULES_EVIDENCE = {
    "system": "C3-symmetric disks / benzene-1,3,5-tricarboxamides",
    "readout": "circular dichroism (net helicity), not modulus",
    "mismatch_penalty_kJ_mol": 0.94,
    "helix_reversal_penalty_kJ_mol": 7.8,
    "verification": ABSTRACT,
    "ref": "JACS 2005, 127, 5490; JACS 2009, 131, 18784; Nat. Commun. 2011, 2, 509",
}


def prior() -> dict:
    """Outcome counts. Deliberately does NOT return fractions.

    With n = 6 and a leader that flips on any single reclassification (see
    `fragility`), percentages would imply a precision the survey does not have.
    Use `verified_prior` to see what survives if unverified entries are dropped.
    """
    counts = {}
    for e in SURVEY:
        counts[e["outcome"]] = counts.get(e["outcome"], 0) + 1
    return {
        "counts": counts,
        "n_systems": len(SURVEY),
        "n_full_text": sum(1 for e in SURVEY if e["verification"] == FULL_TEXT),
        "n_with_composition_series": sum(1 for e in SURVEY
                                         if e["composition_series"]),
        "n_triple_helical": 0,
        "leader": max(counts, key=counts.get),
        "fragile": fragility()["leader_flips_on_one_change"],
    }


def verified_prior() -> dict:
    """Counts using only entries read in full text — which is almost nothing."""
    counts = {}
    for e in SURVEY:
        if e["verification"] == FULL_TEXT:
            counts[e["outcome"]] = counts.get(e["outcome"], 0) + 1
    return {"counts": counts,
            "n_systems": sum(v for v in counts.values()),
            "decisive": len(counts) == 1}


def fragility() -> dict:
    """How many single-entry reclassifications change the leading outcome."""
    base = prior_counts = {}
    for e in SURVEY:
        prior_counts[e["outcome"]] = prior_counts.get(e["outcome"], 0) + 1
    leader = max(prior_counts, key=prior_counts.get)

    flips = []
    for e in SURVEY:
        for alt in (ENHANCED, SORTING, WEAKENED, MAJORITY):
            if alt == e["outcome"]:
                continue
            c = dict(prior_counts)
            c[e["outcome"]] -= 1
            c[alt] = c.get(alt, 0) + 1
            if max(c, key=c.get) != leader:
                flips.append((e["system"], e["outcome"], alt))
                break
    return {"leader": leader, "n_entries": len(SURVEY),
            "n_entries_whose_change_flips_leader": len(flips),
            "flips": flips,
            "leader_flips_on_one_change": len(flips) > 0}


def folds() -> list:
    """Reported racemic/enantiopure modulus ratios, with provenance.

    `fold_verified` is separate from `verification`: an abstract can confirm
    the direction of an effect while saying nothing about its magnitude. The
    phenylalanine entry is exactly that case.
    """
    return [(e["system"], e["fold"], e["verification"],
             e.get("fold_verified", e["verification"] == FULL_TEXT))
            for e in SURVEY if e["fold"] is not None]


def unverified_magnitudes() -> list:
    """Fold values quoted in this survey that no readable source supports."""
    return [(s, f) for s, f, _, ok in folds() if not ok]
