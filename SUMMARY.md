# mirror-tenocyte — project summary

**What it produced:** one cheap, well-specified experiment that nobody has run,
plus an honest record of four things that were claimed and then withdrawn.

**What it did not produce:** a mirror tenocyte, a validated model, or a
publishable finding.

Scale: 2 400 lines of package code, 1 600 of analysis scripts, 112 automated
checks, 10 figures. **37 model parameters, 8 with any external anchor, 0
tenocyte measurements** — and of those 8, two are partly circular (fitted to a
benchmark, then "derived" from it), one is transferred from a system shown not
to apply, and two are borrowed from small-molecule disks in apolar solvent.

---

## 1. The question dissolved, and that is the first result

The brief was: simulate how inverting the chirality of a tenocyte's environment
changes mechanotransduction and collagen output.

That answered itself almost immediately. Integrin recognition is stereospecific,
so a **pure** mirror matrix presents no bindable epitope and the cell behaves as
if on a non-adhesive surface. There is no distinctive "mirror tenocyte
phenotype" — there is a cell on nothing.

Worse, the model made that trivially explicit: with chirality entering only
through recognition, the phenotype depends on the single product
`Ltot·(1−χ)` to machine precision (spread **2.7e-14**). Chirality was not an
independent variable at all; it was a disguised ligand-density axis.

Everything afterwards concerns **mixtures**, which is a chiral-biomaterials
question with a cell-biological readout. The project name is a historical
accident.

---

## 2. What survived every round of criticism

Without exception, these are the claims that arrived with a number and a
calculation. Nothing that arrived as an analogy survived.

| result | number | status |
|---|---|---|
| **Heterochiral triple helices do not form** — 7.87 kcal/mol per D residue (MD), and **experimentally**, a single D-Asp prevents triple-helix formation in PBS (Shah 1999) | 10⁻¹⁶⁷ | settled; the only claim here with experimental support. Solvent- and sequence-context-dependent, but both caveats concern *single* substitutions |
| **All-D collagen cannot be synthesised** — record mirror-image protein is 358 residues; collagen α1(I) helical domain is 1014 | 2.8× beyond | forces the work onto collagen-mimetic peptides |
| **The six-ODE cascade reduces exactly to `F(L_eff, E_eff)`** — three assembly hypotheses at equal effective stiffness give bit-identical states; Jacobian has zero negative off-diagonal couplings | 0.00e+00 | proved, with a test that *should* fail if the model ever earns its complexity |
| **Two D-hydroxyprolines are sold and only one is the enantiomer** — *trans* (CAS 139262-20-7) vs *cis* (214852-45-6) | — | silent-failure trap in the protocol |
| **The low-composition pilot is nearly free** — 0–10 % D uses 0.86 mg of the expensive enantiomer | 0.86 mg | staging recommendation |

---

## 3. What was retracted

The retractions are the most informative part of the record, and are preserved
in [HISTORY.md](HISTORY.md).

| claim | why it was withdrawn |
|---|---|
| ⑥ "dip and rebound" in the coupled dose-response | Largely an artefact of an **assumed** quadratic assembly curve. Deriving the curve from capping kinetics cut the rebound from 0.344 to 0.003, and across plausible parameters it swings 0.00–0.53. Unusable. |
| ⑧ "run the cell experiment at 8 kPa" | Optimal for detecting frustration **alone**. A fourth hypothesis, found later, made it 8× worse. |
| Bistability / hysteresis as a prediction | Real only in a hand-tuned regime, with a coexistence window ~1 % wide in χ. Kept as a labelled negative result. |
| "Racemic enhancement leads the prior" | The promoting argument (Wallach's rule, solid-state packing) was an **analogy**, while the demoting argument was a Boltzmann calculation — and it was the same cross-boundary transfer error the project had diagnosed one turn earlier. Now **unresolved**. |

A pattern worth stating plainly: **every time something looked exciting it was
overrated in the same turn and corrected in the next.** The corrections came
from looking outward — has anyone measured this, can it be built, can the
reagent be bought — not from refining the model further.

---

## 4. The open question, and the experiment that closes it

Collagen assembles hierarchically: chain → triple helix → microfibril → fibril.
**Level 1 is settled** (above). **Level 2 — do already-sorted D and L triple
helices pack into the same fibril? — is unconstrained by any data I could
find**, after three searches.

Four outcomes, each documented in *other* systems, none in a triple helix:

| mode | predicted modulus curve |
|---|---|
| self-sorting | flat |
| capping | collapses at a few % D, stays dead |
| majority rules | flat, then dips only near racemic |
| racemic enhancement | rises, peaks at racemic |

They differ in **shape**, so a single pure-vs-racemic contrast cannot separate
them — all four agree at *both* endpoints.

**[PROTOCOL.md](PROTOCOL.md) is the deliverable.** Two 36-mer peptides
(L and D `(PKG)₄(POG)₄(DOG)₄`), 13 compositions × 3 replicates, rheology plus
CD on the same samples. Weeks of a peptide chemist's time. Key design points,
each of which is easy to get wrong:

- **Sample log-spaced below 1 % D.** Capping may collapse the gel there; a
  conventional 0/25/50/75/100 % series would miss it entirely.
- **Add CD.** `majority_rules` has the best quantitative backing of the four
  (exact 1D Ising, two measured energies) and **zero** mechanical data points —
  every study of it reads out CD. Rheology alone is blind to it.
- **Fix the protocol before choosing the assay.** Under an inadequate
  single-point design the readout ranking *inverts* — it would say measure pFAK
  and that nuclear YAP is useless, the opposite of the truth.

---

## 5. Is it publishable?

**Not as a research paper.** There is no result. The predictions that survive
are restatements of "cells sense stiffness"; the model has zero validation and
29 of 37 parameters assumed outright; the cell layer is provably a monotone
reparametrisation of one scalar, so a reviewer would rightly ask why the ODEs
exist; and the strongest argument in the project (level 1) is someone else's
molecular dynamics, correctly applied (though it does now carry independent
experimental support, found late).

**Legitimately postable as a preprint.** This is the honest answer, and it costs
nothing. A preprint would claim only what is true: *here is a decomposition of
the problem, here is why level 1 is settled and level 2 is not, here are four
distinguishable hypotheses, and here is the cheapest experiment that separates
them — including three ways to design it wrongly.* That is a design contribution,
not a finding, and a preprint is the venue where that distinction is acceptable.
It also stakes the design as citable prior art if anyone later runs it.

**One transferable methodological point** may be worth more than the biology: a
signalling cascade built as a monotone chain **cannot distinguish upstream
mechanism, only magnitude**, and this is checkable in a few lines
(`reduced.verify_reduction`, Jacobian sign pattern). Anyone with a
multi-node mechanotransduction ODE model can run that check on their own work.
It is a paragraph, not a paper — but it generalises.

**What would make it a real paper:** the rheology curve. With it, the model
becomes the interpretive framework for a first measurement of how chirality
controls collagen-mimetic matrix mechanics, and the four hypotheses collapse to
one settled answer. Without it, there is nothing to report — and that is a
limitation no amount of further modelling can remove.

---

## 6. Why it stopped

Three exits, all closed: **no laboratory** (no peptide synthesis available),
**no literature access** (4 of 6 survey entries behind paywalls, only 2 read in
full), **and the model is exhausted** — its cell layer is proved decorative and
its remaining parameters cannot be constrained by more computation.

The evidence layer reflects this honestly: 2 of 6 entries verified, the leading
outcome flips on any single reclassification, and restricted to verified entries
it is a **1–1 tie** that says nothing.

---

| file | what it is |
|---|---|
| [PROTOCOL.md](PROTOCOL.md) | the deliverable — bench-ready, every number regenerable |
| [README.md](README.md) | current state of the model and what is assumed |
| [HISTORY.md](HISTORY.md) | chronological record, including all retractions |
| `mirror_tenocyte/` | model, assembly hypotheses, clutch, reduction, evidence, feasibility |
| `scripts/` | one per claim; each regenerates its own numbers |
| `tests/test_model.py` | 112 checks, several of which enforce honest framing |
