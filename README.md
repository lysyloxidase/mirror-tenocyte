# mirror-tenocyte

An *in silico* model of how **mixed-chirality (D/L) collagen matrices** would
affect tenocyte mechanotransduction — integrin engagement → FAK → RhoA/ROCK
tension → YAP/TAZ → collagen output.

> **What this project actually turned out to be.** It began as "what happens to
> a tenocyte in a mirror-image environment". That question answered itself
> almost immediately and unexcitingly: integrin recognition is stereospecific,
> so a *pure* mirror matrix presents no bindable epitope and the cell behaves as
> if on a non-adhesive surface. Nothing further follows.
>
> Everything interesting here concerns **mixtures**. That is a chiral
> biomaterials question with a cell-biological readout, and the name is now a
> historical accident. Read it as such.

> Scope note: purely computational, concerning tenocyte cell biology and gel
> mechanics. Unrelated to the biosecurity discussion around synthesising
> mirror-image *organisms*. Nothing here involves building chiral matter.

> Not validated against tenocyte data — **zero** tenocyte-specific measurements
> are used. Not suitable for clinical or therapeutic inference. It is a tool
> for designing experiments.

**New here? Read [SUMMARY.md](SUMMARY.md)** — what the project produced, what it
retracted, and whether it is publishable, in one page.

---

## Bottom line

1. **A pure mirror matrix ≡ a non-adhesive substrate.** Structural, not a
   discovery. The required control is D-ECM vs. passivated (PEG) surface.
2. **Whether chirality is an independent variable at all depends on one
   unmeasured thing:** how D and L collagen behave when mixed. **Four** outcomes
   are possible, all realised in *other* systems. They differ in the *shape* of
   the composition curve, not just its sign — so a single contrast cannot
   separate them.
3. **The decisive first experiment needs no cells** — rheology on a D/L gel
   composition series, sampled at both low D fraction *and* near racemic. It is
   also the only one currently feasible: **all-D collagen cannot be
   synthesised** (3× beyond the record), so the work must use collagen-mimetic
   peptides, whose gels are too soft for the cell protocol.
4. **The cell model contributes less than its size suggests** (see below). The
   discriminating content lives entirely in the assembly layer.

---

## What the model actually is

The six-ODE cascade reads the environment through **exactly two scalars**:

```
L_eff = bindable ligand density      (chirality.effective_ligand)
E_eff = stiffness after assembly     (assembly.effective_stiffness)

phenotype = F(L_eff, E_eff)
```

This is a structural fact, verified to `0.00e+00`: three different assembly
hypotheses that produce the same `E_eff` produce **bit-identical** states. The
Jacobian has **zero negative off-diagonal couplings** — it is a monotone chain.

Consequences, stated plainly:

- The cascade **cannot distinguish mechanism, only magnitude**. Any claim that
  "integrin/FAK/YAP signalling distinguishes the chirality hypotheses" would be
  overclaiming.
- `reduced.ReducedModel` evaluates `F` by interpolation, **~600× faster** than
  integrating the ODEs.
- If a future channel is added that the two scalars do *not* mediate — a
  chirality-sensing receptor, a topology term, a non-monotone feedback — the
  reduction **breaks**, and `test_model_reduces_to_two_scalars` fails. That
  failure would be good news: it would mean the cascade had started to earn
  its keep.

What the cascade *does* contribute: the shape of `F` (saturation, the switch,
the relative weighting of adhesion vs. tension), and named intermediate
observables an experimentalist can stain for.

---

## The open question: what do D and L collagen do when mixed?

Collagen assembles hierarchically: **chain → triple helix → microfibril →
fibril**. Chirality mismatch can act at more than one level, and the levels are
not equivalent.

**Level 1 (chain → triple helix) — settled, no free parameter.** A single D-Ala
substitution costs **7.87 kcal/mol** ([MD, J. Phys. Chem. B 2009](https://pubs.acs.org/doi/10.1021/jp808690m)),
comparable to the Gly→Ala mutations causing osteogenesis imperfecta. At 37 °C
the probability of a heterochiral helix underflows after ~30 mismatched
residues; a triple-helical domain has ~1000. **Heterochiral triple helices do
not form.**

**Level 2 (triple helix → fibril) — completely unconstrained by any data I
could find.** Already-sorted homochiral helices pack laterally via charge pairs
(achiral, could still pair) in a staggered array (chiral, geometry may not
match). Four outcomes, all documented in *other* systems:

| mode | mechanism | evidence | where the signal sits |
|---|---|---|---|
| `self_sorting` | D and L form separate perfect fibrils | peptide amphiphiles; Aβ forms homochiral sheets | nowhere — flat |
| `frustrated_kinetic` | mismatched monomers cap growing fibrils | inferred from chiral crystal-growth poisoning | **extremes** — halves at ~5% D |
| `majority_rules` | minority enantiomer **absorbed** into the majority screw sense | measured: MMP 0.94 vs HRP 7.8 kJ/mol; exact 1D Ising theory | **centre** — window χ∈[0.36, 0.64] |
| `racemic_enhanced` | heterochiral coassembly, stiffer than either pure form | MAX1/DMAX1: 800 vs 200 Pa — but see below | centre, opposite sign |

**`majority_rules` is the best-grounded of the four.** Chiral copolymerisation
has a mature quantitative theory — a one-dimensional two-component Ising model,
exact in the long-chain limit. Two measured energies govern it: the **mismatch
penalty** (0.94 kJ/mol, a monomer in its non-preferred screw sense) and the
**helix reversal penalty** (7.8 kJ/mol, a domain wall). Because reversal costs
~8× more than mismatch, the polymer does **not** expel or cap the minority
enantiomer — it absorbs it. Net helicity follows
`m = sinh(βH)/√(sinh²(βH) + e^(−4βJ))` with `H = MMP·ee`, `J = HRP/2`,
amplifying by ~21× at 37 °C so helicity saturates by ee ≈ 0.27.

This directly contradicts the capping picture: **10% D barely perturbs
competence (0.991) under majority-rules, versus 0.152 under capping.** The
constants are borrowed from C3-symmetric disks in apolar solvent, not collagen —
the functional form is principled, the numbers are not collagen's.

**`racemic_enhanced` — demoted for level 1, UNRESOLVED for level 2.** The MAX1
stiffening works through a "rippled" β-sheet of *alternating L and D strands* —
a **level-1** mixing mechanism, which collagen forbids. That argument stands and
level 1 remains settled.

An earlier version of this README then *re-promoted* it for level 2, arguing
that lateral packing is "solid-state" and so falls in the regime where
heterochiral assembly is favoured. **That argument is withdrawn** — see the
survey below. Level 2 is unresolved, as are the other three.

---

## What published D/L systems do — a list of leads, not a prior

Since no synthesis is available, the remaining move was other people's data.
`scripts/run_evidence.py` surveys published enantiomeric-mixture systems
(`mirror_tenocyte/evidence.py` holds the table with per-entry provenance).

**The survey is too weak to function as a prior, and is reported as counts
rather than percentages for that reason:**

| outcome | systems | note |
|---|---|---|
| racemic enhancement (peak) | 3 of 6 | MAX1 4× (verified); phenylalanine 9× (**magnitude unsupported**) |
| self-sorting (flat) | 2 of 6 | one verified |
| weakened on mixing (dip) | 1 of 6 | unverified |
| majority rules | **0** mechanical studies | measured by CD only |

Four reasons not to lean on it:

- **Only 2 of 6 entries were read in full text.** The rest come from
  search-result summaries or indirect citation. Outcome *direction* and effect
  *magnitude* are tracked separately, because an abstract can confirm one while
  saying nothing about the other: the phenylalanine abstract (reached via an
  institutional repository) confirms only that mixing gives "mechanically more
  robust" structures and **contains no numbers at all**, so the quoted
  53.1 / 5.8 / 1.8 GPa rest on a search summary alone. **2 of the 4 quoted
  magnitudes** are unsupported this way.
- **The leading outcome flips if any single entry is reclassified** — 3 of the
  6 entries have that power. Restricted to the two verified entries the survey
  is a **1–1 tie** and says nothing.
- **Zero entries are triple helices.** Two of the three "enhanced" entries are
  not aqueous gels at all; phenylalanine is a solidified amino acid with a
  modulus in GPa, three orders of magnitude from a hydrogel.
- **Publication bias is directional.** "Racemic is 9× stiffer" publishes; "we
  mixed them and nothing happened" does not.

> **A retracted argument, kept visible.** This section previously argued that
> solution phase favours homochiral assembly while solid-state packing favours
> heterochiral (the Aβ reconciliation plus Wallach's rule), that collagen's
> level-2 lateral packing is solid-state, and therefore that enhancement led the
> prior. **Withdrawn.** Wallach's rule concerns molecular crystals of small
> molecules; a collagen fibril is a hydrated quasi-hexagonal array of ~300 nm
> rods, not a molecular crystal, and the transfer was never checked. It is the
> same error the project had just diagnosed in refusing to carry MAX1 across to
> collagen — and the rigour was asymmetric: the level-1 demotion rests on a
> Boltzmann calculation, the re-promotion rested on an analogy.

**What the survey is actually good for** is one design consequence, independent
of everything above: `majority_rules` has the best quantitative backing of the
four hypotheses (exact 1D Ising, two measured energies) and contributes **zero**
mechanical data points, because every study of it reads out circular dichroism.
A rheology-only experiment would be blind to it. Hence the CD arm in
[PROTOCOL.md](PROTOCOL.md) — same samples, cheap, orthogonal axis.

---

## Feasibility — checked late, and it changes the recommendation

Eight stages refined an experimental recommendation without once asking whether
the material can be made. It cannot, as originally stated.

**All-D collagen is not synthesisable.** Mirror-image proteins are made by total
chemical synthesis, and the record is a **358-residue** D-Dpo4 polymerase
([Nature 2016](https://www.nature.com/articles/celldisc201737); prior record
312-residue DapA). Collagen α1(I)'s triple-helical domain is ~1014 residues:

| target | residues | × record | makeable |
|---|---|---|---|
| collagen α1(I) helical domain | 1014 | 2.8 | **no** |
| full pro-α1(I) chain | 1464 | 4.1 | **no** |
| **collagen-mimetic peptide (CMP)** | **36** | 0.10 | **yes, routine** |

Collagen additionally needs three chains and 4-hydroxyproline installed
chemically. **Read this whole project as applying to CMPs, not collagen.** CMPs
traverse the same hierarchy (triple helix → sticky-ended nanofibre → hydrogel)
and their gels are degraded by collagenase at comparable rates.

*What survives the substitution:* the level-1 argument (both are Gly-X-Y triple
helices) and the level-2 question. *What changes:* a CMP triple helix is ~10 nm
rather than ~300 nm, so `λ_c` counts many more helices per unit fibril length.

**And the cell experiment does not survive it.** Native CMP hydrogels reach only
~10–1000 Pa. The two-point protocol degrades sharply below 2 kPa:

| gel stiffness | worst-case separation |
|---|---|
| 2.0 kPa | **0.316** |
| 1.5 kPa | 0.182 |
| 1.0 kPa (top of CMP range) | **0.045** |
| ≤0.5 kPa | ~0 |

A **7× penalty** for staying within reachable stiffness. The reason is
structural: on a very soft gel the cell is already in the low-YAP state, so
further softening cannot move it — there is no headroom downward. Stiffening by
crosslinking is possible but risks masking the signal, since crosslinks bridging
broken fibrils supply mechanical continuity independently of fibril integrity —
exactly what is being measured.

---

## The experiments

> **Standing differs sharply between the two.** The cell-free one is feasible
> today; the cell one is not, until a stiffened gel is demonstrated.

### First, and the only one currently feasible: rheology on a D/L CMP gel series

> **→ [PROTOCOL.md](PROTOCOL.md) is the bench-ready version of this**: peptide
> sequences, the composition series and why its low arm must be log-spaced,
> replicate counts, the symmetry control, and a table mapping curve shape to
> hypothesis.

One measurement that (a) separates the four modes by the *shape* of modulus
vs. composition, (b) locates the halving composition, which pins `λ_c` and
`p_cap`, and (c) determines the sign every cell experiment should expect.
Cheaper than any cell work and constrains the model more.

> **Sampling warning.** The four hypotheses put their signal in *different
> places*: capping at a few percent D, majority-rules and enhancement near
> 50:50. A series sampled only at 0/25/50/75/100% would miss capping entirely;
> a series sampled only near the extremes would miss the other two. **Sample
> both regions.**

### Then — *blocked on gel formulation*: clamped-adhesion, two compositions, ~2 kPa

Functionalise the D/L gel with a **constant** density of natural L-RGD, so
integrin engagement is fixed while mechanics vary. Measure at **χ = 0.10 and
χ = 0.50**. Each hypothesis has a distinct two-point signature:

| mode | at χ=0.10 | at χ=0.50 | pattern |
|---|---|---|---|
| `self_sorting` | +0.000 | +0.000 | flat, flat |
| `frustrated_kinetic` | **−0.323** | −0.324 | **down, down** |
| `majority_rules` | −0.001 | **−0.316** | **flat, down** |
| `racemic_enhanced` | +0.206 | +0.290 | up, up |

Worst-case pairwise separation **0.316** — the number a real experiment's error
bars must beat.

**Why two points and not a better stiffness.** A single racemic-vs-pure contrast
cannot do it: at χ=0.5 frustration and majority-rules are both collapsed and
differ by only **0.008**. The best achievable single-point separation over all
stiffnesses is 0.115, ~3× worse. The modes differ in the *shape* of the
composition curve, and one contrast throws the shape away.

All four modes agree at **both** endpoints — comparing pure-L with pure-D
discriminates nothing at all.

### What to measure

Under the **two-point protocol** (worst-case pairwise separation per readout):

| separation | observable | assay |
|---|---|---|
| **0.397** | `Cint` | focal adhesion size/number (paxillin, vinculin) |
| 0.370 | `YAPn` | **YAP nuclear:cytoplasmic ratio** — standard readout |
| 0.370 | `Col` | COL1A1 expression / procollagen |
| 0.154 | `FAKp` | phospho-FAK Y397 |
| 0.103 | `Tension` | traction force microscopy |
| 0.055 | `Rho` | RhoA activity ← **avoid** |

Signal is compressed by the saturating steps mid-cascade — a property of the
network, not the assays.

> **Readout choice and protocol choice are coupled, and the coupling is a trap.**
> Under the *inadequate* single-point protocol the ranking **inverts**: pFAK
> becomes best (0.117) and nuclear YAP worst (0.006), because the limiting pair
> is then frustration vs majority-rules, which differ most upstream. Choosing
> the readout before fixing the protocol would lead to exactly the wrong assay.
> A test (`test_readout_ranking_depends_on_protocol`) pins both halves of this.

---

## Quick start

Python 3.9+, three dependencies, no install step — everything runs from source.

```bash
git clone <this repo> && cd mirror-tenocyte
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt

./.venv/bin/python tests/test_model.py        # 112 checks, ~2 min
./.venv/bin/python scripts/run_modes.py       # the current headline result
./.venv/bin/python scripts/run_kinetics.py    # assembly curve derivation
./.venv/bin/python scripts/run_protocol_tables.py   # regenerates every number in PROTOCOL.md
```

Most scripts run in seconds. Two are slower because they run stochastic
simulations: `run_clutch.py` (~30 s) and `run_calibration.py --refit` (~30 s).
All scripts are read-only apart from writing figures to `results/`.

```python
from mirror_tenocyte.parameters import get_preset
from mirror_tenocyte.analysis import discriminate_assembly_modes
discriminate_assembly_modes(p=get_preset("clutch_calibrated"))["best_stiffness"]  # 2.0
```

Figures land in `results/`.

| file | role |
|------|------|
| `mirror_tenocyte/model.py` | the six ODEs, integrators, phenotype index |
| `mirror_tenocyte/reduced.py` | **what the model actually is**: `F(L_eff, E_eff)` |
| `mirror_tenocyte/assembly.py` | D/L composition → mechanical competence; the hierarchy |
| `mirror_tenocyte/chirality.py` | `chi` → bindable ligand |
| `mirror_tenocyte/clutch.py` | motor-clutch in pN/nm/s; **diagnostic only**, not imported by `model.py` |
| `mirror_tenocyte/analysis.py` | mode discrimination, observables, sensitivity |
| `mirror_tenocyte/feasibility.py` | **can the experiment be done at all?** synthesis + stiffness walls |
| `mirror_tenocyte/evidence.py` | published D/L outcomes with per-entry provenance — leads, not a prior |
| `scripts/run_evidence.py` | the survey, its fragility, and a retracted argument |
| **`PROTOCOL.md`** | **the deliverable: a bench protocol for the decisive experiment** |
| `scripts/run_protocol_tables.py` | regenerates every number in `PROTOCOL.md` |
| `scripts/run_modes.py` | the four hypotheses and the two-point protocol |
| `scripts/run_kinetics.py` | deriving the assembly curve |
| `scripts/run_calibration.py` | deriving the stiffness constants (`--refit`, ~30 s) |
| `tests/test_model.py` | including tests that enforce the honest framing |

---

## What is assumed, and how much

**34 model parameters. 6 have any external anchor. 0 tenocyte measurements.**

Of those 6, two are partly circular: `KE_adh` and `KE_ten` came from fitting a
clutch whose stiffness mapping was itself tuned to a YAP-switching benchmark.
Varying that knob 4× moves `KE_adh` from 16.3 to 2.7 — it is essentially the
tuned quantity. What survives is **scale-free**: the *ratio* `KE_ten/KE_adh ≈ 3`
holds (2.8–3.6) across that range. That ratio is the real content of the
calibration — the baseline model had 0.75, wrong by ~4×.

Genuinely assumed: the assembly competence functional form, `λ_c`, `p_cap`, the
phenotype-index weights, and all downstream kinetics.

**A caveat on the weight-independence result.** All six observables agree on
sign, so any non-negative weighting agrees too. That sounds like robustness, but
the Jacobian check shows it is **forced by the monotone cascade** — a property
of the architecture, not evidence about biology. Stated correctly: the
conclusion is weight-independent, and that fact carries less information than
it first appears to.

## Explored and not load-bearing

Kept in the code as honest negative results, not used by any current prediction:

- **Bistability / hysteresis** (`get_preset("strong_feedback")`). Exists, but
  only with hand-tuned feedback, and the coexistence window is ~1% wide in χ —
  requiring enantiomeric control better than 1% to detect. Not usable.
- **`chi_specific`** — an off-by-default channel for chirality-specific biology
  beyond adhesion loss. Conceptually the alternative to point 1 above; never
  used in a quantitative prediction.
- **The `affinity` coupling mode** — an alternative to `racemic`; changes where
  the transition sits but no conclusion.

## Honest next steps

Both are subtractive or external — the computational branch is close to
exhausted:

- **Get level-2 data.** Rheology on a D/L **CMP** series — searched three times,
  no published data found. Everything downstream is conditional on it, and it is
  now the only clearly feasible experiment in the project.
- **Solve the stiffness problem** before any cell work: a CMP formulation
  reaching ~2 kPa whose stiffness comes from fibril integrity rather than from
  crosslinks that would mask the signal.
- **Stop modelling.** Three consecutive rounds of added rigour *changed*
  conclusions; the next only reinforced them, and the audit then found the
  cascade to be decorative. Further grounding of the cell layer has low
  expected value.

A frank retrospective — including which predictions were retracted and why — is
in [HISTORY.md](HISTORY.md).
