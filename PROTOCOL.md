# Experimental protocol — chirality and collagen-mimetic peptide assembly

**One measurement, no cells.** Rheology on a D/L collagen-mimetic peptide (CMP)
hydrogel composition series. It decides between four mutually exclusive
hypotheses about what happens when mirror-image triple helices are mixed, and it
is the gate on everything else in this project.

Estimated effort: routine solid-phase synthesis of two 36-mers, **13
compositions × 3 replicates = 39 gels**, standard oscillatory rheology. Weeks,
not months.

---

## 1. The question, and the decision it makes

Homochiral triple helices form readily; heterochiral ones do not. This is the
one claim here with **experimental** backing rather than only a calculation:
host-guest collagen peptides carrying a *single* D-Asp residue fail to form a
triple helix at all in phosphate-buffered saline
([Shah et al., Biopolymers 1999, 49, 297](https://pubmed.ncbi.nlm.nih.gov/10079768/)),
consistent with the ~7.87 kcal/mol per-residue penalty from molecular dynamics.
So in a D/L mixture, **chains sort into pure-D and pure-L triple helices**.

> Two qualifications, neither of which rescues an all-D chain. The exclusion is
> **solvent-dependent** — in 67 % aqueous ethylene glycol helices do form, with
> melting points >30 °C lower — so it is a statement about physiological buffer,
> which is the condition used here. And it is **sequence-context-dependent**:
> mixed D-/L-Asp peptides formed heterotrimers for a Gly-Asp-Ala guest triplet
> but not for Gly-Asp-Hyp. Both concern *single* substitutions; an all-D 36-mer
> carries ~24 stereocentres, so the argument is stronger for whole chains, not
> weaker.

What is *not* known is the next level up: **do D and L triple helices pack into
the same fibril, and what does that do to the modulus?** Four outcomes are each
documented in other supramolecular systems. Nobody has measured which one a
collagen-like helix does.

They differ in the *shape* of modulus-vs-composition, so one curve decides it.

---

## 2. Peptides

Sequence (O = 4-hydroxyproline), 36 residues, from the canonical
multi-hierarchical self-assembling CMP
([O'Leary et al., Nat. Chem. 2011, 3, 821](https://www.nature.com/articles/nchem.1123)):

```
L-CMP :  (Pro-Lys-Gly)4 (Pro-Hyp-Gly)4 (Asp-Hyp-Gly)4      "KOD"
D-CMP :  exact enantiomer — invert every stereocentre
```

For D-CMP use D-Pro, D-Lys, D-Asp and (2R,4S)-4-hydroxyproline. **Glycine is
achiral**, so 12 of the 36 positions are unchanged — the synthesis is a standard
Fmoc SPPS run with mirror-image building blocks, not a special technique. All
required Fmoc building blocks are commercially available.

> ⚠ **Order the right hydroxyproline — two are sold and only one is the
> enantiomer.** Collagen's residue is (2S,4R), i.e. *trans*-4-hydroxy-**L**-proline.
> Its mirror image is (2R,4S) = ***trans*-4-hydroxy-D-proline, Fmoc-D-*trans*-Hyp-OH,
> CAS 139262-20-7**. The also-catalogued **Fmoc-*cis*-D-4-hydroxyproline
> (CAS 214852-45-6) is a diastereomer, not an enantiomer** — using it yields a
> peptide that is *not* the mirror image of L-CMP. That failure is silent: the
> gel would still form, the symmetry control in §6 would break, and the
> asymmetry could easily be misread as a real result. Confirm by CD that the two
> spectra are equal in magnitude and opposite in sign before mixing anything.

Purify both to the same HPLC standard and confirm by CD: the two spectra must be
**mirror images** (equal magnitude, opposite sign). Any asymmetry in magnitude
means unequal purity or incomplete folding, and will propagate into every
subsequent modulus ratio.

> **Do not confuse two different "D fractions."** The experimental variable χ is
> the mole fraction of *whole D-CMP chains* deliberately mixed in. Residue-level
> racemisation during synthesis is a different thing: in physiological buffer a
> chain carrying a stray D residue does not enter a triple helix (Shah 1999), so
> it drops out of the assembly rather than acting as a mixed-chirality
> component. It lowers effective concentration; it does not set χ. Note the
> caveat above, though — in one sequence context single-D chains did join mixed
> trimers, so this is a strong tendency rather than an absolute.

---

## 3. Composition series — the part that is easy to get wrong

Mix L-CMP and D-CMP at mole fraction χ of D, holding **total peptide
concentration constant** (0.5 or 1.0 wt % in phosphate buffer, the conditions
under which this peptide is known to gel).

Sample **log-spaced at the low end and linearly near racemic**:

| arm | χ (mol % D) | why |
|---|---|---|
| low, log-spaced | 0, **0.1, 0.3, 1, 3, 10** | capping may collapse the gel below 1% D |
| mid/high, linear | 20, 30, **40**, 50 | majority-rules and enhancement act here; 30–45% is also where CD separates them best (§4) |
| symmetry check | 90, 99, 100 | internal control (see §6) |

**Why the low arm must be logarithmic.** If capping operates, the collapse
position measures the critical fibril length λ_c, and for CMP geometry λ_c is
large — a CMP triple helix is ~10 nm against ~300 nm for collagen, so ~30× more
helices are needed to span the same physical length. Predicted normalised
modulus under capping:

| χ | λ_c=20 | λ_c=100 | λ_c=600 (CMP-appropriate) |
|---|---|---|---|
| 0.1% | 0.999 | 0.989 | **0.734** |
| 0.5% | 0.985 | 0.794 | **0.099** |
| 1% | 0.951 | 0.490 | **0.026** |
| 5% | 0.451 | 0.033 | 0.001 |

A series starting at 5% would see a dead gel at every point and learn only *that*
it died, not *where* — losing the one number that pins λ_c.

### Material required, and how to stage it

At 1 wt % and ~200 µL per rheometer loading (20 mm plate, ~0.5 mm gap), each gel
takes ~2 mg of peptide. For the full 39 gels:

| | gels only | with 3× overage (HPLC, CD, repeats) |
|---|---|---|
| **D-CMP** | 27 mg | ~80 mg |
| **L-CMP** | 51 mg | ~155 mg |

Both are routine custom-synthesis scales.

> **Run the low arm first as a pilot — it costs almost nothing.** All six
> low-arm compositions (0–10 % D, 18 gels) together consume **0.86 mg of
> D-CMP**, because at low χ the mixture is nearly all L. That pilot alone
> separates capping from the other three hypotheses, since capping predicts
> G′/G′(0) ≈ 0.03 at 1 % D while every other hypothesis predicts ≈ 1.0.
> The expensive D peptide is spent almost entirely on the 50–100 % points,
> which are the symmetry control — commit to those only after the pilot.

---

## 4. Measurement

- Oscillatory shear rheology, plate–plate. Report **G′ in the linear
  viscoelastic plateau**, after a strain sweep to confirm you are inside it.
- Gel in situ, same thermal history for every sample; this peptide assembles on
  the hours timescale, so fix an incubation time and hold it constant.
- **n = 3 independent gels per composition** (independent preparations, not
  three measurements of one gel).
- Report everything **normalised: G′(χ)/G′(0)**. The model predicts the ratio,
  not the absolute value, and normalisation removes batch-to-batch offsets.

**Required precision is undemanding.** The hypotheses differ by factors of
3–1000× at the diagnostic compositions. The limiting factor is *gel-to-gel
reproducibility*, not instrument resolution — if replicate G′(0) values scatter
within ±20%, that is comfortably enough.

### Also run CD on the same samples — it catches what rheology misses

Take a circular dichroism spectrum of every composition. It costs almost
nothing (same samples, an instrument most peptide labs already have) and adds an
**orthogonal** axis, because it measures net screw sense rather than network
mechanics.

This matters for one specific reason: **`majority rules` has the best
quantitative backing of the four hypotheses and contributes zero mechanical data
points** — every published study of it reads out CD, not modulus. A
rheology-only experiment could see a modest dip near racemic and misassign it.

Under self-sorting and capping the fibrils stay homochiral, so net CD is just
the population difference — **linear** in composition. Under majority rules the
minority helix is absorbed into the majority screw sense, so net CD follows the
Ising magnetisation — **strongly nonlinear and saturating**:

| χ | linear (sorting / capping) | majority rules | ratio |
|---|---|---|---|
| 0.20 | 0.600 | 0.977 | 1.6× |
| 0.30 | 0.400 | 0.949 | 2.4× |
| **0.40** | **0.200** | **0.833** | **4.2×** |
| 0.45 | 0.100 | 0.601 | 6.0× |
| 0.50 | 0 | 0 | — |

Both go to zero at exactly racemic, so **sample χ = 0.3–0.45**, where the
divergence is largest. A CD magnitude that stays high while the population
balance approaches 50:50 is chiral amplification, and it is unambiguous.

Racemic enhancement should additionally show a *changed spectral shape* (a
heterochiral fibril is a different structure), not merely a changed magnitude.

---

## 5. Classification — read the answer off two compositions

| hypothesis | G′(χ)/G′(0) at **χ = 1%** | at **χ = 50%** | pattern |
|---|---|---|---|
| **self-sorting** — separate perfect D and L fibrils | 1.00 | 1.00 | flat throughout |
| **capping** — mismatched helices terminate growth | **0.03** | 0.00 | collapses immediately, stays dead |
| **majority rules** — minority absorbed into majority screw sense | 1.00 | **0.30** | flat, then dips only near racemic |
| **racemic enhancement** — heterochiral co-packing, stiffer | 1.12 | **4.00** | rises, maximal at racemic |

Adding the CD axis makes the classification two-dimensional and removes the one
genuine ambiguity — `capping` and `majority rules` both produce a dip in modulus,
but only `majority rules` shows nonlinear CD:

| hypothesis | modulus at racemic | CD at χ=0.40 |
|---|---|---|
| self-sorting | 1.00 | 0.20 (linear) |
| capping | 0.00 | 0.20 (linear) |
| majority rules | 0.30 | **0.83 (amplified)** |
| racemic enhancement | 4.00 | changed *shape*, not just magnitude |

Worst-case pairwise separation across these two points is **0.70** (self-sorting
vs majority rules); no two hypotheses are close. The remaining compositions refine the curve: the low arm
locates the capping collapse (→ λ_c), the mid arm locates the majority-rules
window (→ mismatch and helix-reversal penalties).

---

## 6. Controls and pitfalls

**Symmetry is the internal control.** Every hypothesis predicts
G′(χ) = G′(1−χ), because a pure mirror gel is mechanically identical to a pure
natural one — enantiomers are isoenergetic and elastic moduli are
parity-invariant. So **G′(100% D) must equal G′(0% D)** and G′(90%) must equal
G′(10%). If they do not, something is wrong before any interpretation begins:
unequal peptide purity, unequal concentration, or a genuine asymmetry that would
itself be the finding.

**Hold concentration constant, not volume.** The variable is composition; a
concentration gradient masquerading as a chirality effect is the easiest way to
get a spurious curve.

**Comparing pure-L with pure-D tells you nothing.** All four hypotheses agree at
both endpoints. All the information is at intermediate composition. (It is also
the obvious first experiment, which is why it is worth stating.)

**Do not start with cells.** A companion cell experiment is *not* currently
feasible: native CMP gels reach only ~10–1000 Pa, and the cell-based
discrimination degrades ~7× below 2 kPa because a cell on a very soft gel is
already in the low-YAP state and further softening cannot move it. Stiffening by
crosslinking risks masking the very signal being measured, since crosslinks
bridging broken fibrils supply mechanical continuity independently of fibril
integrity.

---

## 7. What each outcome licenses

| result | interpretation | next step |
|---|---|---|
| **flat** | self-sorting; chirality acts only through integrin recognition | the cell-biology question collapses to "non-adhesive substrate" — stop |
| **collapse below a few % D** | capping; fibril growth is poisoned by trace opposite enantiomer | measure the halving point → λ_c; chirality is a genuine mechanical axis |
| **flat then dip near racemic** | majority rules; the minority helix is absorbed | fit the 1D Ising form → mismatch and helix-reversal penalties for a collagen-like helix (**new numbers, nobody has these**) |
| **rise, peak at racemic** | heterochiral co-packing stiffens the fibril | mechanistically surprising given level-1 exclusion — worth structural follow-up (cryo-EM, CD of the fibril) |

Any of the last three makes the cell experiment worth designing properly, with a
gel formulation that reaches ~2 kPa without masking crosslinks.

---

*Every number here is regenerated by `scripts/run_protocol_tables.py`; run it
to reproduce the tables. The
functional forms are principled, but the constants for majority-rules and
enhancement are borrowed from other supramolecular systems (C3-symmetric disks;
MAX1 β-hairpin) — no collagen-like system has been measured. That is exactly
what this experiment fixes.*
