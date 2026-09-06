# Development history

This is the chronological record of how the model was built, including
**predictions that were later retracted or reversed**. It is kept deliberately:
the retractions are the most informative part, and a reader who wants to know
how much to trust a claim should be able to see how it survived.

For the current position, see [README.md](README.md). Nothing here is
authoritative where it conflicts with the README.

Summary of what happened to each stage:

| stage | outcome |
|---|---|
| 1 — baseline model | found the `Ltot·(1−χ)` degeneracy: chirality was a disguised ligand axis |
| 2 — bistability | hysteresis only in a hand-tuned regime, 1% wide window — **abandoned** |
| 3 — motor clutch | challenged prediction ⑥; clutch failed its own calibration check |
| 4 — clutch calibration | corrected the tension constant; predictions got ~3× stronger |
| 5 — assembly kinetics | **retracted ⑥** — the rebound was largely an artefact of an assumed curve |
| 6 — literature check | found a third hypothesis that **inverts ⑦**; moved optimal gel 8 → 2 kPa |
| 7 — hierarchy | level-1 thermodynamics excludes the level-1 mixing mechanism for collagen |
| 8 — observables | conclusions are weight-independent — but that turned out to be structural |
| audit | the six-ODE cascade reduces exactly to a function of two scalars |

---

# mirror-tenocyte

An *in silico*, hypothesis-generating model of how **environmental chirality
inversion** ("mirror" ECM / biomolecules) would propagate through tenocyte
**mechanotransduction** — integrin engagement → FAK → RhoA/ROCK tension →
YAP/TAZ → collagen output.

The goal is **not** a quantitatively validated simulator. It is a compact,
transparent dynamical model whose purpose is to (1) force the underlying
assumptions to be explicit, and (2) generate falsifiable predictions that can
later be tested at the bench.

> Scope note: this project is purely computational and concerns tenocyte cell
> biology. It is unrelated to — and should not be conflated with — the
> biosecurity discussion around synthesizing mirror-image *organisms*
> ("mirror life"). Nothing here involves building chiral matter.

---

## The core biological idea (and the one subtlety that reshapes the model)

Integrin–ligand recognition is **stereospecific**. Integrins read chiral
epitopes (collagen's GFOGER motif via the α1/α2/α11 I-domain; RGD in
fibronectin/vitronectin) through a lock-and-key interface. A mirror-image
(all-D) ligand presents the *enantiomeric* surface, which the natural, chiral
integrin cannot productively engage.

**Consequence:** you cannot "half-invert" a molecule. A graded chirality knob is
only meaningful as an **enantiomeric composition** of a population — a racemic
fraction `chi` of D vs. L ligands. And to first order, a *fully* mirror matrix
does not create a novel signal; it **removes the recognizable ligand**, so the
cell should behave as if on a **non-adhesive substrate**. That reframing is
baked into the model and is itself the first prediction (see below).

Two coupling modes are provided (`chir_mode` in `parameters.py`):

- **`racemic`** (default): fraction `chi` of ligands are D and non-bindable, so
  bindable ligand density scales as `(1 - chi)`. Here `chi = 1` ≡ non-adhesive.
- **`affinity`**: every ligand's productive affinity decays smoothly as
  `exp(-beta_chi · chi)` — a graded loss of binding rather than a population swap.

A separate, **off-by-default** channel (`chi_specific`) lets you encode the
competing hypothesis that a D-environment does something *beyond* adhesion loss
(altered fibril assembly / local mechanics, chirality-sensing receptors,
differential protein adsorption — all of which have empirical precedent in the
chiral-biomaterials literature).

---

## Model structure

Five coupled nodes, each relaxing toward an upstream-driven set-point, with two
positive-feedback loops (tension→adhesion clutch reinforcement; nuclear
YAP→its own retention) that let the network behave in a switch-like way.

```
 chirality chi ─┐
                ▼
   ECM ligand (bindable) ──► integrin clutch (Cint) ──► FAK ──► RhoA/ROCK
                ▲                     ▲                            │
   stiffness ───┘                     │ tension feedback           ▼
   (loads clutch)                     └──────────────── cytoskeletal Tension
                                                              │
   stiffness (resists) ───────────────────────────────────►  │
                                                              ▼
                                                   YAP/TAZ nuclear (YAPn)
                                                              │
                                                              ▼
                                                   collagen-I output (Col)
```

State vector (all normalised to `[0,1]`): `[Cint, FAKp, Rho, Tension, YAPn, Col]`.
Stiffness enters **twice** — to load the clutch and to let actomyosin build
tension against a resistant substrate — so a soft matrix collapses tension even
with abundant ligand (the classic soft-substrate loss of nuclear YAP).

Files:

| file | role |
|------|------|
| `mirror_tenocyte/parameters.py` | default kinetics + environment, with literature notes |
| `mirror_tenocyte/chirality.py`  | `chi` → bindable ligand; optional D-specific channel |
| `mirror_tenocyte/assembly.py`   | fibrillogenesis: D/L composition → mechanical competence |
| `mirror_tenocyte/clutch.py`     | motor-clutch in pN/nm/s (mean-field + stochastic); diagnostic only |
| `mirror_tenocyte/model.py`      | the ODE right-hand side, integrators, phenotype index |
| `mirror_tenocyte/simulate.py`   | dose-response, stiffness grid, hysteresis, plots |
| `mirror_tenocyte/analysis.py`   | χ* localisation, bistability detection, invariance test |
| `scripts/run_dose_response.py`  | core figures + summary table |
| `scripts/run_timecourse.py`     | node trajectories at several `chi` |
| `scripts/run_analysis.py`       | **checks every claim below and prints a verdict** |
| `scripts/run_assembly.py`       | **the assembly layer: three discriminating signatures** |
| `scripts/run_clutch.py`         | **adversarial test of ⑥ (~30 s runtime)** |
| `scripts/run_calibration.py`    | **derives the stiffness constants, re-tests ⑥–⑧** |
| `scripts/run_kinetics.py`       | **derives the assembly curve; retracts ⑥, promotes ⑦** |
| `scripts/run_modes.py`          | **three hypotheses, the two-level hierarchy; revises ⑧ to 2 kPa** |
| `tests/test_model.py`           | sanity tests (bounds, monotonicity, presets, assembly) |

---

## Quick start

```bash
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt
./.venv/bin/python scripts/run_dose_response.py
./.venv/bin/python scripts/run_timecourse.py
./.venv/bin/python scripts/run_analysis.py    # claim-by-claim verdicts
./.venv/bin/python scripts/run_assembly.py    # assembly layer, discriminating tests
./.venv/bin/python tests/test_model.py        # sanity tests
```

Figures land in `results/`. Minimal programmatic use:

```python
from mirror_tenocyte.simulate import dose_response
dr = dose_response()          # sweep chi in [0,1] at default 30 kPa
print(dr["phenotype"])        # tenogenic activation index vs chi

from mirror_tenocyte.parameters import get_preset
from mirror_tenocyte.analysis import detect_bistability
detect_bistability(get_preset("strong_feedback"))["max_gap"]   # 0.646
```

---

## Predictions — with the numbers behind them

Everything below is *measured* from the simulation by
`scripts/run_analysis.py`, not asserted. Run it to reproduce these figures.
Crucially, the claims are separated into what follows from the model's
**structure** (assumptions made explicit) and what is a genuine **result**.

### ⓪ The sharpest prediction: chirality is a disguised ligand-density axis

In `racemic` mode the phenotype depends **only on the product `Ltot · (1 − χ)`**,
to machine precision (max spread across ligand densities: **2.7e-14**):

| `Ltot·(1−χ)` | `Ltot`=0.5 | 1.0 | 2.0 | 4.0 | spread |
|---|---|---|---|---|---|
| 0.300 | 0.7011 | 0.7011 | 0.7011 | 0.7011 | 0 |
| 0.085 | 0.4185 | 0.4185 | 0.4185 | 0.4185 | 2.7e-14 |
| 0.060 | 0.0193 | 0.0193 | 0.0193 | 0.0193 | 3.8e-16 |

So diluting L-ligand and enantiomerically inverting it should be **exactly
interchangeable**. → **This is the experiment to run.** Titrate an L/D ratio
series against a matched L-dilution series. If the two curves superimpose, the
parsimonious model stands and "mirror ECM" is simply adhesion loss. If they
diverge, the parsimonious model is *refuted* and chirality-specific biology is
demonstrated — which is the interesting outcome, and what `chi_specific` exists
to model.

### ① Mirror ECM ≡ non-adhesion — an assumption, not a discovery

A fully D matrix and a passivated (PEG) surface give an identical steady state
(max difference **0.0e+00**). This **holds by construction**: χ enters the model
only through bindable ligand. It is stated as a prediction because it is the
required **experimental control** — D-ECM vs. PEG. Enabling `chi_specific`
breaks the equivalence (difference 0.707), which is how the competing
hypothesis is represented.

### ② The transition is switch-like — but χ* can't identify the mechanism

| coupling mode | χ* | phenotype | max slope |
|---|---|---|---|
| `racemic`  | 0.915 | 0.704 → 0.000 | 18.6 /unit χ |
| `affinity` | 0.625 | 0.704 → 0.000 | 18.2 /unit χ |

Sharp in both, but χ* is **mode-dependent**, so measuring a threshold alone
does not distinguish the two coupling mechanisms — the density-titration test
(⓪) is needed for that.

### ③ Hysteresis is *conditional* — the default model does not predict it

| preset | hysteresis gap | coexistence window | regime |
|---|---|---|---|
| `parsimonious` (default) | 2.5e-11 | — | **monostable** — fully reversible |
| `strong_feedback` | 0.646 at χ=0.91 | χ ∈ [0.9125, 0.9225], width **0.010** | **bistable** |

The parsimonious defaults are monostable: the sharp drop is completely
reversible. Bistability only emerges with stronger clutch cooperativity and
tension feedback. Which regime real tenocytes occupy is an **empirical
question this model cannot settle**. See `results/preset_comparison.png`.

> **Practical caveat — this is the weakest prediction.** Even in the bistable
> regime, the coexistence window is only **~1% wide in χ**. The phenotype gap is
> large (0.65) but confined to a narrow band, so detecting hysteresis would
> require controlling enantiomeric fraction to better than ~1%. Combined with
> the fact that the regime itself was found by parameter search, claim ③ should
> be treated as a motivating possibility, **not** as something worth designing a
> first experiment around. Start with ⓪.

### ④ Chirality × stiffness — a genuine interaction

| stiffness | phenotype at χ=0 | χ* | effect size |
|---|---|---|---|
| 1 kPa | 0.061 | 0.120 | **0.061** |
| 3 kPa | 0.623 | 0.635 | 0.623 |
| 10 kPa | 0.684 | 0.865 | 0.684 |
| 30 kPa | 0.704 | 0.915 | 0.704 |
| 100 kPa | 0.711 | 0.935 | 0.711 |

Both halves hold: χ* shifts monotonically right with stiffness (a stiff matrix
**buffers** partial inversion), and on soft matrix the cell is already YAP-low
so chirality barely matters (effect size 0.06 vs 0.70). → Predicts that a
chirality effect will be **hardest to detect on soft substrates** — a practical
warning for experimental design.

### ⑤ Ligand-density rescue

Raising total ligand rescues an intermediate χ: at χ=0.95, phenotype goes
0.008 → 0.596 → 0.694 → 0.703 as `Ltot` goes 1 → 2 → 4 → 8. This is the same
degeneracy as ⓪ viewed along the other axis.

---

## Resolving ⓪ — when does chirality become a *real* variable?

Result ⓪ is a problem for the whole project: if χ only removes recognisable
ligand, there is no chirality biology to study. Chirality becomes an independent
axis only if D/L composition changes something *other* than bindable ligand
density. `assembly.py` supplies that something — **fibrillogenesis**.

**The physics.** Collagen I is a right-handed superhelix of three left-handed
PPII-like chains; its mirror image is the left-handed superhelix of right-handed
chains. Two consequences pull in opposite directions:

1. A **pure** mirror matrix assembles perfectly. Enantiomeric structures are
   isoenergetic and elastic moduli are parity-invariant, so an all-D fibril is
   exactly as stiff as an all-L one. (Precedent: synthetic D-proteins fold to
   mirror structures with identical thermodynamics — Milton et al. 1992,
   *Science* 256:1445, D-HIV-1 protease.) A pure mirror tendon would be
   mechanically normal and merely twist the other way.
2. A **mixed** matrix may not assemble at all: a D monomer cannot pack into a
   growing L triple helix, so mismatched subunits cap fibrils.

So mechanical competence is **non-monotonic and symmetric about χ=0.5**, while
bindable ligand is monotonic. Two channels, not one — the degeneracy breaks.

Two competing assembly hypotheses, `assembly_mode`:

| mode | physical picture | degeneracy | verdict |
|---|---|---|---|
| `self_sorting` | D and L segregate into separate, internally perfect fibrils | spread 2.7e-14 → **intact** | reduces *exactly* to the parsimonious model |
| `frustrated` | mismatched monomers poison fibril growth | spread **0.478** → **broken** | chirality is a genuine axis |

> **The conceptual payoff:** the `Ltot·(1−χ)` degeneracy is the *fingerprint of
> self-sorting*, and its violation is the *fingerprint of frustrated assembly*.
> Result ⓪ stops being an embarrassment and becomes a readout.

### Three signatures that discriminate the hypotheses

See `results/assembly.png`. Effective stiffness under frustration runs
30 → 9.75 → **3.0** → 9.75 → 30 kPa across χ — a symmetric V with full recovery
at pure mirror.

**⑥ The rebound — ⚠️ RETRACTED as a primary experiment. See Stage 5.**
In the *coupled* case — D/L collagen gels across a composition series — the
phenotype under frustration dips near χ≈0.53 and rebounds before collapsing as
adhesion vanishes. Recognition-loss models cannot produce a rebound, so this
looked like the easy discriminating experiment. Its *direction* survived an
explicit motor-clutch challenge (Stage 3–4), but its *existence* did not survive
deriving the assembly curve (Stage 5): the rebound swings from ~0 to ~0.53
across capping parameters nobody has measured, and nearly vanishes under the
derived curve. It was substantially an artefact of the assumed quadratic.
Do not design around it. Its *size* remains a useful readout of capping
efficiency once the other parameters are pinned down.

**⑦ The clamped-adhesion experiment (the clean one).** Functionalise a D/L
collagen gel with a *constant* density of natural L-RGD, so integrin engagement
is fixed while mechanics are inverted. Under recognition loss this curve must be
**flat** (dip 0.0000); under frustration it is **U-shaped and symmetric**
(dip 0.081), with pure-D indistinguishable from pure-L to 6 decimals. This
isolates the mechanical channel with no adhesion confound.

**⑧ Pick the right gel stiffness or you will miss it.** The contrast depends
sharply on base stiffness, because the effect requires the frustrated midpoint
to fall *below* the mechanosensing threshold while the pure compositions stay
above it:

| base stiffness | stiffness at racemic | dip |
|---|---|---|
| 3 kPa | 0.30 kPa | 0.622 |
| **8 kPa** | **0.80 kPa** | **0.662** ← best |
| 15 kPa | 1.50 kPa | 0.263 |
| 30 kPa (tendon-like) | 3.00 kPa | **0.081** |
| 100 kPa | 10.0 kPa | 0.027 |

At tendon-like 30 kPa the effect is ~8× smaller than at 8 kPa, because even a
frustrated gel stays above threshold. **Running this at physiological tendon
stiffness would likely hide the phenomenon** — a concrete design warning that
would be easy to get wrong.

---

## Adversarial test of ⑥ — an explicit motor-clutch (`clutch.py`)

Predictions ⑥–⑧ all inherit one assumption from the baseline Hill function:
**that softening the matrix reduces the mechanotransduction signal.** If the
real stiffness response is biphasic and tenocytes sit *above* their clutch
optimum, softening moves them *towards* it, the signal rises, and the dip
becomes a peak. `scripts/run_clutch.py` tests that assumption rather than
trusting it.

A mean-field motor-clutch plus a stochastic Chan–Odde simulation
(Chan & Odde 2008, *Science* 322:1687; catch-slip bond after Kong et al. 2009;
force-gated talin reinforcement after Elosegui-Artola et al. 2016), in real
units (pN, nm, s).

**Stage 1 — the uncalibrated module said ⑥ inverts.** With literature defaults
(slip bond, `F_bell`=2 pN) the stochastic simulation reproduces the published
clutch optimum (peak at k_sub ≈ 0.46 pN/nm), and the physiological range lands
*entirely above it*: engaged fraction **falls** with stiffness (0.72×). Taken at
face value that inverts ⑥ into a peak at racemic.

Also worth recording: the mean-field approximation **does** find an optimum when
clutch binding is slow (biphasic at `k_on`=0.3, optimum 1.67 kPa) — an earlier
claim here that it was *structurally* incapable was wrong. Its real failure is
quantitative: it misplaces the optimum ~6× and overestimates the stiff branch.

**Stage 2 — but the module failed its own benchmark, so its verdict was void.**
Cell spreading, adhesion size and nuclear YAP all *increase* from ~1 to ~100 kPa.
The uncalibrated module predicts the opposite. A verdict from a model that
contradicts the best-established fact in the field is not evidence.

**Stage 3 — calibration, and the verdict reverses.** Scanning bond type,
adhesion growth and the stiffness mapping isolated what matters:

| change | effect on the sign |
|---|---|
| adhesion growth (clutch + myosin recruitment) | scales everything, sign **unchanged** |
| slip → **catch bond** | **sign flips to rising** |
| stiffness→spring mapping | sets *where* the transition sits, not its direction |

Integrin–ECM bonds (α5β1, αVβ3) *are* catch bonds, so the weak slip bond was a
wrong choice rather than a neutral simplification. `get_calibrated_parameters()`
now **passes the benchmark**: engaged clutches rise **7×** across 1–100 kPa with
half-maximum at **7.6 kPa**, inside the observed 5–15 kPa YAP-switch window.

> **An independent check that was not fitted:** per-clutch force runs
> 2.4 → 10.6 pN and crosses the ~5 pN talin-unfolding threshold at 8–15 kPa —
> exactly the stiffness where the YAP switch is observed — while staying below
> the ~13 pN catch-bond stabilisation optimum. The forces land where
> mechanosensing biology says they should, and none of that was tuned.

### Revised verdict: ⑥ survives

The calibrated signal **rises** 5.4× across 1–100 kPa, so softening the matrix
*lowers* it. A frustrated D/L matrix therefore gives a **dip** at racemic
composition, not a peak. The earlier inversion was an artefact of a weak slip
bond. See `results/clutch.png`.

The calibration also **independently reinforces ⑧**: the response saturates
above ~30 kPa (62 → 66 → 68 engaged clutches), so at tendon-like stiffness the
cell sits near the top of its sensing range and softening must be large before
anything moves. Two unrelated arguments now point to running the experiment on
a **softer (~8 kPa) gel**.

> **Still worth measuring first.** The sign is now supported by a calibrated
> model rather than assumed, but it rests on a *fitted* stiffness mapping
> (`kappa_per_kPa` = 0.1, chosen so the half-max landed in the YAP window; it
> corresponds to an integrin-nanocluster radius of ~50 nm, which is the right
> scale for a ~75-clutch module). A stiffness series on *these* tenocytes
> (FA size or YAP N:C ratio) would confirm the working point directly, and costs
> far less than a chirality series read backwards.

### Stage 4 — folding the clutch back in: constants *derived*, not assumed

Predictions ⑥–⑧ were all computed with **assumed** Hill constants. The final
step derives them from the calibrated clutch instead (`scripts/run_calibration.py`,
`--refit` re-derives in ~30 s and reproduces the stored values to the digit).

**The clutch response is fitted by a Hill function to better than 1% variance**
— R²=0.992 (engaged clutches → adhesion), R²=0.999 (traction → tension). That is
a result in itself: the coarse-grained form in `model.py` was a *justified*
reduction all along. Only its constants were wrong.

| term | assumed | derived from clutch | verdict |
|---|---|---|---|
| adhesion | K=8.0, h=1 | K=5.62, h=0.869 | ~right |
| tension | K=6.0, h=1 | **K=16.97**, h=1.208 | **3× too low** |

The baseline builds cytoskeletal tension far too easily on soft matrix. Since
the assembly layer acts *by softening the matrix*, everything downstream of
tension is affected — so this is not a cosmetic correction.

**Do the predictions survive?** All three, and two get stronger:

| prediction | assumed | derived | change |
|---|---|---|---|
| ⑥ rebound size | 0.122 | **0.344** | 2.8× stronger |
| ⑦ decoupled dip | 0.081 | **0.227** | 2.8× stronger |
| ⑧ best gel stiffness | 8 kPa | **8 kPa** | unchanged |

Grounding the stiffness response made the effects *larger and easier to detect*,
and the usable window widened: the dip at 15 kPa goes 0.263 → 0.474, and even at
tendon-like 30 kPa it goes 0.081 → 0.227. See `results/calibration.png`.

Use the derived constants with `get_preset("clutch_calibrated")`. They are not
the default, because the defaults are what every earlier number in this README
was computed with and changing them silently would invalidate the record.

The clutch module itself is still **not imported by `model.py`** — a test
enforces this. Its contribution is the four constants, not a runtime dependency.

---

## Stage 5 — deriving the assembly curve, and retracting ⑥

The assembly layer rested on an **assumed** quadratic `1 − s·4χ(1−χ)`, picked
only for having the right endpoints. Predictions ⑥–⑧ all depended on it.
`scripts/run_kinetics.py` replaces it with a curve derived from growth-poisoning
kinetics (`assembly_mode="frustrated_kinetic"`).

**The derivation.** An L fibril adds monomers at rate ∝(1−χ) and is capped by a
D monomer at rate ∝χ·`p_cap`, so mean length before termination is
`(1−χ)/(χ·p_cap)`; by mirror symmetry `λ_D = χ/((1−χ)·p_cap)`. A fibril bears
load only above a critical length `λ_c`, so, weighting each population by mass:

```
A(χ) = (1−χ)·g(λ_L) + χ·g(λ_D),      g(λ) = λ^m / (λ_c^m + λ^m)
```

Endpoints and symmetry are preserved, but the **shape is radically different**:

| | competence halves at |
|---|---|
| assumed quadratic | χ = 0.167 (**17%** D) |
| derived, λ_c=20 | χ = 0.046 (4.6%) |
| derived, λ_c=100 | χ = 0.010 (**1%**) |

Capping is far more efficient than the quadratic implied — the classic
chiral-poisoning behaviour of crystal growth, where trace opposite enantiomer
arrests growth entirely.

**Consequences.** ⑥ collapses (rebound 0.344 → 0.003) while ⑦ strengthens
(dip 0.227 → **0.689**, near-total phenotype collapse). With efficient capping,
competence recovers only as χ→1, by which point the ligand (~1−χ) is gone — so
there is no window where mechanics have returned but adhesion remains.

**Sensitivity to the unmeasured kinetics** (`λ_c`, `p_cap`) — this is the part
that matters:

| | range across plausible parameters | verdict |
|---|---|---|
| ⑥ rebound | 0.000 – 0.528 | **fragile, unusable** |
| ⑦ decoupled dip | 0.023 – 0.689, ≥0.67 for λ_c≥20 | **robust** |

See `results/kinetics.png`.

### The experiment to run first needs no cells

> **Rheology on the D/L gel series alone.** One measurement that
> (a) separates self-sorting (flat modulus) from frustration (dip),
> (b) locates the halving composition, which pins `λ_c` and `p_cap`, and
> (c) therefore predicts whether a rebound should exist at all.
> It is cheaper than any cell experiment and constrains the model more than
> either of them. Everything downstream is conditional on it.

> **Sampling warning.** Under the derived curve the transition sits at a *few
> percent* D, not near 50:50. A composition series sampled at 0/25/50/75/100%
> — the obvious choice — would **miss it entirely**.

Then ⑦ (clamped adhesion) as the primary cell experiment, at ~8 kPa (⑧).

---

## Stage 6 — a literature check adds a third hypothesis that inverts ⑦

The model had two assembly hypotheses. A search of the peptide self-assembly
literature found a third, **better documented than either**, which reverses the
sign of the headline prediction: L and D can **coassemble into a fibril stiffer
than either pure form**.

Racemic MAX1/DMAX1 β-hairpin hydrogel reaches **G′ ≈ 800 Pa vs ≈200 Pa** for
either enantiopure gel — a 4× enhancement — via a heterochiral "rippled" β-sheet
with alternating L and D strands, exactly as [Pauling and Corey predicted in
1953](https://pmc.ncbi.nlm.nih.gov/articles/PMC3202337/). The reported Job plot
has a clean symmetric maximum at 0.5 mole fraction. Modelled as
`assembly_mode="racemic_enhanced"`.

| hypothesis | evidence | prediction ⑦ |
|---|---|---|
| `self_sorting` | documented in peptide amphiphiles; Aβ forms homochiral sheets | **flat** (0.000) |
| `frustrated_kinetic` | inferred from chiral crystal-growth poisoning | **dip** (−0.324) |
| `racemic_enhanced` | **measured**: MAX1/DMAX1, 4× stiffer | **peak** (+0.290) |

All three are empirically realised in *different* systems. Collagen is a PPII
triple helix, structurally unlike a β-hairpin, so **neither result transfers to
it automatically** — which mode collagen follows is genuinely unknown.

**Why this is good news for the experiment.** The three hypotheses differ in the
*sign* of ⑦, and sign is a far more robust readout than magnitude — it does not
depend on the capping parameters that made ⑥ unusable.

**Why it changes the design.** All three modes agree at *both* endpoints — a
pure mirror matrix is mechanically identical to a pure natural one under every
hypothesis. So comparing pure-L with pure-D **cannot discriminate at all**; the
information is entirely at intermediate composition.

### ⑧ revised: 8 kPa → **2 kPa**

Detecting a *dip* wants a stiff starting point; detecting a *peak* wants a soft
one, or the cell is already saturated and stiffening does nothing. Optimising
for one hypothesis hides another:

| base stiffness | self-sorting | frustration | enhancement | min gap |
|---|---|---|---|---|
| 1 kPa | 0.000 | −0.047 | +0.476 | 0.047 |
| **2 kPa** | 0.000 | **−0.324** | **+0.290** | **0.290** |
| 8 kPa (old advice) | 0.000 | −0.614 | +0.076 | 0.076 |
| 30 kPa | 0.000 | −0.689 | +0.021 | 0.021 |

The earlier 8 kPa recommendation was optimal for detecting **frustration alone**,
because enhancement had not been modelled. There the peak is only 0.076 and
would likely be missed. See `results/modes.png`.

> **Revised recommendation:** run the clamped-adhesion experiment on a **soft
> (~2 kPa) gel**. This costs sensitivity to frustration but is the only setting
> that keeps all three hypotheses distinguishable.

> And still prior to all of it: **rheology on the gel series with no cells**.
> It reads the stiffness-vs-composition curve directly and settles which mode
> operates before a single cell is plated.

---

## Stage 7 — collagen is hierarchical, and that settles half the question

Stage 6 imported an alarming result from β-hairpin peptides. Stage 7 asks
whether it actually transfers, and the answer is **no** — for a specific
structural reason with a quantitative basis.

Collagen assembles in stages: **chain → triple helix → microfibril → fibril**.
Chirality mismatch can act at more than one level, and the levels are not
equivalent.

**Level 1 (chain → triple helix) — settled, no free parameter.** Three chains
interdigitate with Gly at the helix core. MD puts a single D-Ala substitution at
**7.87 kcal/mol** — comparable to the Gly→Ala mutations that cause osteogenesis
imperfecta — producing a kink that breaks the helical register
([J. Phys. Chem. B 2009, 113, 8983](https://pubs.acs.org/doi/10.1021/jp808690m)).
At 37 °C:

| mismatched residues | P(heterochiral helix) |
|---|---|
| 1 | 2.8×10⁻⁶ |
| 3 | 2.3×10⁻¹⁷ |
| 30 | 3.4×10⁻¹⁶⁷ |

A triple-helical domain has ~1000 residues, so this underflows long before a
full chain. **Heterochiral triple helices do not form; chains sort perfectly.**

**Why that defuses Stage 6.** The MAX1/DMAX1 stiffening works through a
"rippled" β-sheet of *alternating L and D strands* — a **level-1 mixing**
mechanism. Collagen forbids level-1 mixing, so the strongest empirical support
for `racemic_enhanced` **does not carry over**. Enhancement via level-2 helix
co-packing is not excluded, but it is now *unsupported* rather than
*well-evidenced*. The hypothesis stays in the code, demoted.

**Level 2 (triple helix → fibril) is where the open question actually lives.**
Already-sorted homochiral helices pack laterally via charge pairs (achiral, so
they could still pair) in a staggered array (chiral, so geometry may not match).

> **Unit correction this forces:** at level 2 the capping "monomer" of
> `frustrated_kinetic` is a whole **triple helix (~300 nm rod)**, not a residue.
> `λ_c` must be read as a count of triple helices. The numbers are unchanged;
> their interpretation is not.

---

## Stage 8 — the phenotype index turns out not to matter

Every headline number so far is expressed in `phenotype_index`, a composite
(`0.5·Col + 0.3·YAPn + 0.2·Tension`) whose weights **I chose by hand**. No
experimentalist measures that. Stage 8 re-expresses the three-way mode
discrimination in each individually measurable state.

| min gap | observable | what it is at the bench |
|---|---|---|
| **0.311** | `Cint` | focal adhesion size/number (paxillin, vinculin) ← best |
| 0.300 | `Col` | COL1A1 expression / procollagen staining |
| 0.300 | `YAPn` | **YAP nuclear:cytoplasmic ratio** — the standard readout |
| 0.118 | `Tension` | traction force microscopy; stress-fibre density |
| 0.117 | `FAKp` | phospho-FAK Y397 |
| 0.041 | `Rho` | RhoA activity (GTP pulldown, FRET) ← avoid |
| 0.290 | *composite* | the hand-made index, beaten by three single readouts |

**The weights don't matter — and that's a guarantee, not luck.** All six
observables agree on the *sign* for every mode (flat / negative / positive).
If every observable agrees on a sign, then **any** non-negative weighting of
them agrees too. So the qualitative conclusion is weight-independent *by
construction* — a stronger statement than randomly sampling weights would give.
A test enforces this sign-consistency.

The composite was adding nothing and can be dropped in favour of a direct
measurement.

### What to measure

> Use **focal-adhesion size, nuclear YAP ratio, or COL1A1** — all ≈0.30.
> **Do not rely on RhoA activity (0.041) or pFAK (0.117):** the signal is
> compressed mid-cascade, and those assays would need ~7× better precision to
> resolve the same difference.

That mid-cascade compression is itself worth noting — it is a property of the
network, not of the assays. Signals entering a saturating step lose dynamic
range, and Rho sits in the most saturated part of this one.

> Still assumed rather than derived: all downstream kinetics. `λ_c` and `p_cap`
> are exposed parameters, not measured values. Level-2 behaviour is
> unconstrained by any data I could find. And the clutch is calibrated to a
> generic mechanobiology benchmark, not to tenocyte data.

---

## Limitations (read before over-interpreting)

- **The headline result is a tautology made visible, not a discovery.** Because
  χ enters only through bindable ligand, the exact `Ltot·(1−χ)` degeneracy (⓪)
  is a *consequence of the model's structure*. Its value is that it makes the
  assumption falsifiable and tells you precisely which experiment breaks it —
  not that it reveals anything about tenocytes on its own.
- **Parameters were tuned, and that matters.** The `strong_feedback` preset was
  found by searching for a regime that produces hysteresis. It therefore
  demonstrates that the network *can* be bistable; it is **not** evidence that
  tenocytes *are*. Reported as a conditional prediction for exactly this reason.
- Phenomenological, **not fitted**; parameters are illustrative order-of-magnitude
  values. Absolute timescales are uncalibrated — steady states are the object of
  study.
- Deterministic and single-cell; no stochasticity, no spatial ECM/topology, no
  explicit talin/vinculin catch-bond mechanics (the clutch is coarse-grained).
- Tenogenic identity is reduced to a scalar; Scx/Mkx/Tnmd programs and
  collagen-III are not resolved.
- **The assembly layer is a hypothesis with a chosen functional form.** The
  `1 − s·4χ(1−χ)` competence curve is the simplest symmetric form with the right
  endpoints; it is not derived from nucleation-growth kinetics. Its *qualitative*
  consequences (symmetry, recovery at pure mirror, a rebound) are robust to that
  choice; the specific dip depths and the χ≈0.78 rebound location are not.
- Whether real D/L collagen self-sorts or frustrates is **unknown to me** — the
  model deliberately encodes both and shows how to tell them apart, rather than
  assuming one.
- **The clutch is calibrated to a generic benchmark, not to tenocytes.** It now
  reproduces both the published motor-clutch optimum and physiological stiffness
  sensing, but `kappa_per_kPa` was fitted to put the half-max in the YAP-switch
  window. Its support for ⑥ is therefore a consistency argument, not
  independent evidence. It is deliberately not imported by `model.py`.
- Not validated against data and **not suitable for clinical or therapeutic
  inference** — it is a tool for designing experiments, nothing more.

## Sensible next steps

- ~~Add a supramolecular-assembly term so mixed D/L collagen alters *mechanics*~~
  — **done**, see `assembly.py` and the section above.
- ~~Replace the coarse clutch with an explicit molecular-clutch / catch-bond layer~~
  — **done** (`clutch.py`), and ~~it needs calibrating before its verdict counts~~
  — **also done**: it now passes the stiffness-sensing benchmark and upholds ⑥.
- ~~consider wiring the clutch into `model.py` to replace the Hill function~~
  — **resolved differently and better**: the clutch response *is* a Hill
  function (R²>0.99), so it was folded in as derived constants rather than a
  runtime dependency. See Stage 4.
- Calibrate against **tenocyte** data rather than a generic mechanobiology
  benchmark. This is now the binding constraint on the whole project.
- The `kappa_per_kPa` mapping is fitted; measuring adhesion size directly on the
  substrates of interest would replace the fit with a constraint.
- ~~Derive the assembly competence curve from nucleation-growth kinetics~~
  — **done** (Stage 5). It retracted ⑥ and promoted ⑦.
- Measure `lambda_c` and `cap_efficiency` by rheology on a D/L gel series. This
  is now the single highest-value experiment and needs no cells.
- Derive assembly competence from nucleation-growth kinetics with explicit
  D-capping, instead of the assumed symmetric form.
- Calibrate YAP nuclear/cytoplasmic dynamics to published stiffness sweeps.
- Bifurcation analysis over (`chi`, stiffness) to map the switch/hysteresis
  region formally.
