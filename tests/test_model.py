#!/usr/bin/env python3
"""Sanity tests for the mirror-tenocyte model.

Dependency-light: run directly, no pytest required.

    ./.venv/bin/python tests/test_model.py
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import copy
import numpy as np

from mirror_tenocyte.parameters import (
    get_default_parameters, get_default_environment, get_preset,
)
from mirror_tenocyte.chirality import effective_ligand
from mirror_tenocyte.model import run_to_steady_state, phenotype_index
from mirror_tenocyte.analysis import detect_bistability, product_invariance

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  — {detail}" if detail else ""))
    if not cond:
        FAILURES.append(name)


def test_states_bounded():
    print("\nstate variables stay within [0, 1]")
    p = get_default_parameters()
    for chi in [0.0, 0.25, 0.5, 0.75, 0.9, 1.0]:
        for E in [0.5, 5.0, 30.0, 200.0]:
            env = get_default_environment()
            env["chi"], env["Estiff"] = chi, E
            s = run_to_steady_state(p, env)
            ok = np.all(s >= -1e-9) and np.all(s <= 1.0 + 1e-9)
            if not ok:
                check(f"bounded at chi={chi}, E={E}", False, f"state={s}")
                return
    check("bounded across chi x stiffness grid", True)


def test_monotonic_in_chirality():
    print("\nphenotype is non-increasing as the matrix is inverted")
    p = get_default_parameters()
    env = get_default_environment()
    chi = np.linspace(0, 1, 61)
    phen = []
    for c in chi:
        env["chi"] = float(c)
        phen.append(phenotype_index(run_to_steady_state(p, env)))
    d = np.diff(np.array(phen))
    check("monotonically non-increasing", bool(np.all(d <= 1e-6)),
          f"largest increase = {d.max():.2e}")


def test_mirror_equals_nonadhesive():
    print("\nfully mirror matrix == non-adhesive substrate (racemic mode)")
    p = get_default_parameters()
    m = get_default_environment(); m["chi"] = 1.0
    peg = get_default_environment(); peg["Ltot"] = 0.0
    a, b = run_to_steady_state(p, m), run_to_steady_state(p, peg)
    check("steady states coincide", float(np.max(np.abs(a - b))) < 1e-9,
          f"max diff = {np.max(np.abs(a-b)):.2e}")


def test_effective_ligand_modes():
    print("\nchirality -> bindable ligand mapping")
    p = get_default_parameters()
    env = get_default_environment()
    env["chi"] = 0.0
    check("racemic: chi=0 gives full ligand",
          abs(effective_ligand(env, p) - env["Ltot"]) < 1e-12)
    env["chi"] = 1.0
    check("racemic: chi=1 gives zero ligand",
          abs(effective_ligand(env, p)) < 1e-12)
    env["chir_mode"] = "affinity"
    check("affinity: chi=1 leaves a small residual, not exactly zero",
          0.0 < effective_ligand(env, p) < 0.05)
    env["chi"] = 1.5
    env["chir_mode"] = "racemic"
    check("chi is clamped above 1", effective_ligand(env, p) >= 0.0)


def test_soft_substrate_kills_yap():
    print("\nsoft matrix suppresses nuclear YAP even with full L-ligand")
    p = get_default_parameters()
    soft = get_default_environment(); soft["Estiff"] = 0.5
    stiff = get_default_environment(); stiff["Estiff"] = 30.0
    ys = run_to_steady_state(p, soft)[4]
    yh = run_to_steady_state(p, stiff)[4]
    check("YAP_nuc(soft) << YAP_nuc(stiff)", ys < 0.2 * yh,
          f"soft={ys:.3f}, stiff={yh:.3f}")


def test_product_degeneracy():
    print("\nphenotype depends only on Ltot*(1-chi) in racemic mode")
    inv = product_invariance()
    check("exact degeneracy holds", inv["invariant"],
          f"max spread = {inv['max_spread']:.2e}")


def test_presets_differ_in_bistability():
    print("\npresets occupy different dynamical regimes")
    mono = detect_bistability(get_preset("parsimonious"))
    bi = detect_bistability(get_preset("strong_feedback"))
    check("parsimonious preset is monostable", not mono["bistable"],
          f"gap = {mono['max_gap']:.2e}")
    check("strong_feedback preset is bistable", bi["bistable"],
          f"gap = {bi['max_gap']:.3f} at chi = {bi['chi_at_max_gap']:.2f}")


def test_chi_specific_channel():
    print("\noptional chirality-specific channel breaks the equivalence")
    p = get_default_parameters()
    p["chi_specific"] = 0.8
    m = get_default_environment(); m["chi"] = 1.0
    peg = get_default_environment(); peg["Ltot"] = 0.0
    a, b = run_to_steady_state(p, m), run_to_steady_state(p, peg)
    check("mirror != non-adhesive once chi_specific is on",
          float(np.max(np.abs(a - b))) > 1e-3,
          f"max diff = {np.max(np.abs(a-b)):.3f}")


def _env(mode, **kw):
    e = get_default_environment()
    e["assembly_mode"] = mode
    e.update(kw)
    return e


def test_assembly_symmetry():
    print("\nassembly competence is symmetric and unity at pure compositions")
    from mirror_tenocyte.assembly import assembly_competence, effective_stiffness
    p = get_default_parameters()
    check("competence = 1 at chi_struct = 0",
          abs(assembly_competence(_env("frustrated", chi_struct=0.0), p) - 1.0) < 1e-12)
    check("competence = 1 at chi_struct = 1 (pure mirror assembles fine)",
          abs(assembly_competence(_env("frustrated", chi_struct=1.0), p) - 1.0) < 1e-12)
    a, b = (assembly_competence(_env("frustrated", chi_struct=c), p) for c in (0.3, 0.7))
    check("symmetric about chi = 0.5", abs(a - b) < 1e-12, f"{a:.6f} vs {b:.6f}")
    lo = assembly_competence(_env("frustrated", chi_struct=0.5), p)
    check("minimum at racemic composition", lo < a, f"racemic={lo:.4f} < 0.3/0.7={a:.4f}")
    check("pure mirror is mechanically identical to pure natural",
          abs(effective_stiffness(_env("frustrated", chi_struct=1.0), p)
              - effective_stiffness(_env("frustrated", chi_struct=0.0), p)) < 1e-12)


def test_assembly_mode_default_is_inert():
    print("\ndefault assembly mode does not perturb the baseline model")
    p = get_default_parameters()
    a = run_to_steady_state(p, get_default_environment())
    b = run_to_steady_state(p, _env("none"))
    check("assembly_mode='none' matches the untouched default",
          float(np.max(np.abs(a - b))) < 1e-12)


def test_self_sorting_reduces_to_parsimonious():
    print("\nself-sorting is indistinguishable from pure recognition loss")
    p = get_default_parameters()
    for chi in [0.0, 0.3, 0.7, 1.0]:
        a = run_to_steady_state(p, _env("none", chi=chi))
        b = run_to_steady_state(p, _env("self_sorting", chi=chi))
        if float(np.max(np.abs(a - b))) > 1e-9:
            check(f"self_sorting == none at chi={chi}", False)
            return
    check("self_sorting == none across chi", True)


def test_degeneracy_breaks_only_under_frustration():
    print("\nthe Ltot*(1-chi) degeneracy breaks only under assembly frustration")
    p = get_default_parameters()
    for mode, expect_invariant in [("none", True), ("self_sorting", True),
                                   ("frustrated", False)]:
        inv = product_invariance(p=p, env=_env(mode))
        check(f"{mode}: invariant == {expect_invariant}",
              inv["invariant"] == expect_invariant,
              f"max spread = {inv['max_spread']:.2e}")


def test_rebound_signature():
    print("\nnon-monotonic rebound appears only under frustration")
    from mirror_tenocyte.analysis import detect_rebound
    p = get_default_parameters()
    for mode, expect in [("none", False), ("self_sorting", False),
                         ("frustrated", True)]:
        r = detect_rebound(p=p, env=_env(mode))
        check(f"{mode}: rebound == {expect}", r["rebound"] == expect,
              f"rise = {r['rise']:+.4f}")


def test_decoupled_scan_shape():
    print("\ndecoupled scan is flat without frustration, U-shaped with it")
    from mirror_tenocyte.analysis import decoupled_scan
    p = get_default_parameters()
    flat = decoupled_scan(p=p, env=_env("self_sorting"), n=21)
    check("self_sorting: flat when adhesion is clamped", flat["dip"] < 1e-9,
          f"dip = {flat['dip']:.2e}")
    u = decoupled_scan(p=p, env=_env("frustrated"), n=21)
    check("frustrated: U-shaped", u["dip"] > 1e-3, f"dip = {u['dip']:.4f}")
    check("frustrated: pure-D matches pure-L", u["symmetric"],
          f"L={u['phen_pure_L']:.6f}, D={u['phen_pure_D']:.6f}")


def test_clutch_physics():
    print("\nmotor-clutch module reproduces known behaviour")
    from mirror_tenocyte.clutch import (
        get_clutch_parameters, koff, solve_clutch, simulate_clutch, find_optimum,
    )
    p = get_clutch_parameters()

    F = np.linspace(0, 40, 400)
    kc = np.array([koff(f, p, "catch") for f in F])
    i = int(np.argmin(kc))
    check("catch bond has an interior koff minimum", 0 < i < len(F) - 1,
          f"minimum at F = {F[i]:.1f} pN")
    check("slip bond unbinding increases with force",
          koff(10, p, "slip") > koff(0, p, "slip"))

    s_soft = solve_clutch(0.01, p, "slip")
    s_stiff = solve_clutch(100.0, p, "slip")
    check("force per clutch rises with substrate stiffness",
          s_stiff["F_per_clutch"] > s_soft["F_per_clutch"],
          f"{s_soft['F_per_clutch']:.2f} -> {s_stiff['F_per_clutch']:.2f} pN")
    check("traction never exceeds motor stall force",
          s_stiff["traction"] <= p["n_motor"] * p["F_stall"] + 1e-6)

    q = get_clutch_parameters()
    q["k_on"] = 0.3
    tr = [simulate_clutch(k, q, "slip", t_max=4.0, dt=2e-4, seed=3)["traction"]
          for k in [0.03, 0.5, 50.0]]
    check("stochastic clutch is biphasic (peak at intermediate stiffness)",
          tr[1] > tr[0] and tr[1] > tr[2],
          f"traction {tr[0]:.1f} / {tr[1]:.1f} / {tr[2]:.1f} pN")

    check("mean-field biphasicity is parameter-dependent, not structural",
          find_optimum(q, "slip")["biphasic"],
          "biphasic at k_on=0.3")


def test_clutch_calibration():
    print("\ncalibrated clutch reproduces physiological stiffness sensing")
    from mirror_tenocyte.clutch import (
        get_calibrated_parameters, calibration_report, stiffness_response,
    )
    cal = calibration_report(nseed=1, t_max=4.0)
    check("signal rises across 1-100 kPa", cal["rising"],
          f"{cal['n_engaged'][1]:.1f} -> {cal['n_engaged'][-1]:.1f}")
    check("dynamic range >= 2x", cal["dynamic_range"] >= 2.0,
          f"{cal['dynamic_range']:.1f}x")
    check("half-max within the observed YAP switch window (5-15 kPa)",
          cal["half_max_in_window"], f"{cal['half_max_kPa']:.1f} kPa")

    # Not fitted: per-clutch force should straddle the ~5 pN talin threshold.
    f = cal["F_per_clutch"]
    check("per-clutch force crosses the talin unfolding threshold",
          f[1] < 5.0 < f[-1], f"{f[1]:.2f} -> {f[-1]:.2f} pN")
    check("per-clutch force stays physically sane (< 30 pN)", f.max() < 30.0)


def test_bond_type_sets_the_sign():
    print("\nbond type, not adhesion growth, decides the sign of stiffness sensing")
    from mirror_tenocyte.clutch import get_calibrated_parameters, stiffness_response
    p = get_calibrated_parameters()
    E = np.array([1.0, 100.0])
    slip = stiffness_response(E, p=p, bond="slip", nseed=1, t_max=4.0)["n_engaged"]
    catch = stiffness_response(E, p=p, bond="catch", nseed=1, t_max=4.0)["n_engaged"]
    check("catch bond gives rising stiffness sensing", catch[1] > catch[0],
          f"{catch[0]:.1f} -> {catch[1]:.1f}")
    check("slip bond does not", slip[1] <= slip[0] * 1.15,
          f"{slip[0]:.1f} -> {slip[1]:.1f}")


def test_hill_exponent_defaults_are_inert():
    print("\ngeneralised Hill exponents default to the original E/(E+K) form")
    p = get_default_parameters()
    check("hE_adh and hE_ten default to 1", p["hE_adh"] == 1.0 and p["hE_ten"] == 1.0)
    s = run_to_steady_state(p, get_default_environment())
    check("default steady state unchanged by the generalisation",
          np.allclose(s, [0.7953, 0.7047, 0.6789, 0.6292, 0.7223, 0.7223], atol=2e-3),
          f"{np.round(s, 4)}")


def test_clutch_calibrated_preset():
    print("\nclutch-derived preset corrects the tension constant")
    from mirror_tenocyte.parameters import PRESETS
    cal = PRESETS["clutch_calibrated"]
    check("adhesion constant stays the same order as assumed",
          3.0 < cal["KE_adh"] < 12.0, f"KE_adh = {cal['KE_adh']}")
    check("tension constant is ~3x the assumed 6.0",
          cal["KE_ten"] > 2.5 * 6.0, f"KE_ten = {cal['KE_ten']}")


def test_predictions_survive_calibration():
    print("\npredictions 6-8 survive the clutch-derived constants")
    from mirror_tenocyte.analysis import (
        detect_rebound, decoupled_scan, design_sensitivity,
    )
    p = get_preset("clutch_calibrated")
    e = _env("frustrated")

    reb = detect_rebound(p=p, env=e)
    check("6: rebound still present", reb["rebound"], f"rise = {reb['rise']:+.4f}")
    check("6: rebound is larger than with assumed constants",
          reb["rise"] > 0.1216, f"{reb['rise']:.4f} vs 0.1216")

    dec = decoupled_scan(p=p, env=e, n=21)
    check("7: decoupled dip still U-shaped and symmetric",
          dec["dip"] > 1e-3 and dec["symmetric"], f"dip = {dec['dip']:.4f}")

    des = design_sensitivity(p=p)
    check("8: optimal gel stiffness unchanged at 8 kPa",
          abs(des["best_stiffness"] - 8.0) < 1e-9,
          f"{des['best_stiffness']:.0f} kPa")


def test_kinetic_competence_endpoints_and_symmetry():
    print("\nkinetic competence keeps the required endpoints and symmetry")
    from mirror_tenocyte.assembly import assembly_competence
    p = get_default_parameters()

    def A(c):
        return assembly_competence(
            {"assembly_mode": "frustrated_kinetic", "chi_struct": c}, p)

    check("A(0) = 1 (pure natural assembles)", abs(A(0.0) - 1.0) < 1e-9)
    check("A(1) = 1 (pure mirror assembles just as well)", abs(A(1.0) - 1.0) < 1e-9)
    check("symmetric about chi = 0.5", abs(A(0.2) - A(0.8)) < 1e-9,
          f"{A(0.2):.6f} vs {A(0.8):.6f}")
    check("minimum at racemic", A(0.5) < A(0.2))


def test_capping_is_sharper_than_quadratic():
    print("\ncapping kinetics collapse competence faster than the quadratic guess")
    from mirror_tenocyte.analysis import competence_halving_point
    p = get_default_parameters()
    h = competence_halving_point(p)
    check("derived curve halves well below the quadratic's 16.7%",
          h is not None and h < 0.10, f"halves at chi = {h:.4f}")


def test_rebound_is_fragile_dip_is_robust():
    print("\nprediction 6 is parameter-fragile; prediction 7 is robust")
    from mirror_tenocyte.analysis import kinetics_sensitivity
    s = kinetics_sensitivity(p=get_preset("clutch_calibrated"))
    lo, hi = s["rebound_range"]
    check("6: rebound swings across the plausible kinetic range",
          hi - lo > 0.3, f"spans {lo:.4f} to {hi:.4f}")
    robust = [r["dip"] for r in s["rows"] if r["lambda_c"] >= 20]
    check("7: dip stays large for lambda_c >= 20",
          min(robust) > 0.5, f"min dip = {min(robust):.4f}")


def test_racemic_enhanced_mode():
    print("\nracemic coassembly mode matches the measured MAX1 behaviour")
    from mirror_tenocyte.assembly import assembly_competence, effective_stiffness
    p = get_default_parameters()

    def A(c):
        return assembly_competence(
            {"assembly_mode": "racemic_enhanced", "chi_struct": c}, p)

    check("A(0) = A(1) = 1 (pure forms are the reference)",
          abs(A(0.0) - 1.0) < 1e-12 and abs(A(1.0) - 1.0) < 1e-12)
    check("racemic is 4x stiffer, as measured for MAX1/DMAX1",
          abs(A(0.5) - 4.0) < 1e-9, f"A(0.5) = {A(0.5):.3f}")
    check("symmetric about 0.5 (matches the reported Job plot)",
          abs(A(0.25) - A(0.75)) < 1e-12)
    check("maximum at racemic, not at the pure forms", A(0.5) > A(0.25) > A(0.0))


def test_three_modes_give_three_signs():
    print("\nthe three assembly hypotheses predict three different signs")
    from mirror_tenocyte.analysis import decoupled_scan, ASSEMBLY_MODES
    p = get_preset("clutch_calibrated")
    signs = {}
    for m in ASSEMBLY_MODES:
        e = _env(m)
        e["Estiff"] = 2.0
        d = decoupled_scan(p=p, env=e, n=31)
        signs[m] = d["phen_racemic"] - d["phen_pure_L"]
    check("self-sorting is flat", abs(signs["self_sorting"]) < 1e-6,
          f"{signs['self_sorting']:+.4f}")
    check("frustration gives a DIP", signs["frustrated_kinetic"] < -0.05,
          f"{signs['frustrated_kinetic']:+.4f}")
    check("racemic coassembly gives a PEAK", signs["racemic_enhanced"] > 0.05,
          f"{signs['racemic_enhanced']:+.4f}")


def test_single_point_protocol_is_inadequate():
    print("\nsingle racemic-vs-pure contrast cannot separate four hypotheses")
    from mirror_tenocyte.analysis import discriminate_assembly_modes
    d = discriminate_assembly_modes(p=get_preset("clutch_calibrated"))
    check("best achievable single-point separation is poor",
          d["best_separation"] < 0.15,
          f"{d['best_separation']:.4f} at {d['best_stiffness']:.0f} kPa")
    at2 = next(r for r in d["rows"] if r["stiffness"] == 2.0)
    e = at2["effects"]
    check("at 2 kPa frustration and majority-rules are nearly identical",
          abs(e["frustrated_kinetic"] - e["majority_rules"]) < 0.02,
          f"{e['frustrated_kinetic']:.4f} vs {e['majority_rules']:.4f}")
    print("      This is why the protocol needs two compositions, not a better")
    print("      stiffness: at racemic both dip modes are equally collapsed.")


def test_endpoints_cannot_discriminate():
    print("\ncomparing pure-L with pure-D cannot distinguish any hypothesis")
    from mirror_tenocyte.assembly import effective_stiffness
    from mirror_tenocyte.analysis import ASSEMBLY_MODES
    p = get_default_parameters()
    for m in ASSEMBLY_MODES:
        a = effective_stiffness(_env(m, chi_struct=0.0), p)
        b = effective_stiffness(_env(m, chi_struct=1.0), p)
        if abs(a - b) > 1e-9:
            check(f"{m}: endpoints agree", False, f"{a:.4f} vs {b:.4f}")
            return
    check("all modes agree at both pure compositions", True,
          "so the discriminating information is at intermediate chi only")


def test_heterochiral_triple_helix_is_forbidden():
    print("\nlevel 1: heterochiral triple helices cannot form")
    from mirror_tenocyte.assembly import heterochiral_helix_probability
    p = get_default_parameters()
    p1 = heterochiral_helix_probability(p, 1)
    p30 = heterochiral_helix_probability(p, 30)
    check("a single D residue is already strongly disfavoured", p1 < 1e-4,
          f"P = {p1:.3e}")
    check("30 mismatched residues is effectively impossible", p30 < 1e-100,
          f"P = {p30:.3e}")
    check("probability falls monotonically with mismatch count", p30 < p1)
    check("uses the MD-derived penalty, not a guess",
          abs(p["ddG_heterochiral"] - 7.87) < 1e-9,
          "7.87 kcal/mol, J. Phys. Chem. B 2009")


def test_conclusions_do_not_depend_on_phenotype_weights():
    print("\nmode discrimination is independent of the composite index weights")
    from mirror_tenocyte.analysis import observable_discrimination
    od = observable_discrimination(p=get_preset("clutch_calibrated"), stiffness=2.0)
    check("all observables agree on the sign for every mode",
          od["weight_independent"], f"{od['sign_consistent']}")
    check("frustration is negative in every observable",
          all(r["effects"]["frustrated_kinetic"] < 0 for r in od["rows"]))
    check("racemic enhancement is positive in every observable",
          all(r["effects"]["racemic_enhanced"] > 0 for r in od["rows"]))
    check("self-sorting is flat in every observable",
          all(abs(r["effects"]["self_sorting"]) < 1e-6 for r in od["rows"]))


def test_readout_ranking_depends_on_protocol():
    print("\nreadout choice and protocol choice are coupled")
    from mirror_tenocyte.analysis import (
        observable_discrimination, observable_ranking_two_point,
    )
    p = get_preset("clutch_calibrated")

    single = observable_discrimination(p=p, stiffness=2.0)
    two = observable_ranking_two_point(p=p)

    check("under the two-point protocol, adhesion/transcriptional readouts win",
          two["best"] in ("Cint", "YAPn", "Col"),
          f"best = {two['best']} ({two['rows'][0]['worst_separation']:.3f})")
    check("mid-cascade RhoA remains the worst readout under the real protocol",
          two["worst"] == "Rho",
          f"worst = {two['worst']} ({two['rows'][-1]['worst_separation']:.3f})")
    check("best readout is several-fold better than the worst",
          two["rows"][0]["worst_separation"] > 3 * two["rows"][-1]["worst_separation"])
    check("an inadequate single-point protocol INVERTS the ranking",
          single["best"] == "FAKp",
          f"single-point would say measure {single['best']} — an artefact")


def test_model_reduces_to_two_scalars():
    print("\nthe six-ODE cascade reduces exactly to F(L_eff, E_eff)")
    from mirror_tenocyte.reduced import verify_reduction
    v = verify_reduction(p=get_preset("clutch_calibrated"))
    check("environments with equal (L_eff, E_eff) give identical states",
          v["reduction_holds"],
          f"{v['n_cases']} cases, max deviation {v['max_deviation']:.2e}")
    check("deviation is exactly zero, not merely small",
          v["max_deviation"] == 0.0)
    print("      NOTE: if this test ever FAILS, that is good news — it means a")
    print("      channel was added that the two scalars do not mediate, and the")
    print("      signalling cascade has started to carry real information.")


def test_assembly_mode_is_invisible_downstream():
    print("\nthe cascade cannot distinguish assembly mechanism, only magnitude")
    from mirror_tenocyte.assembly import assembly_competence
    p = get_preset("clutch_calibrated")
    # pick a base stiffness for each mode that lands on the same E_eff
    target = 6.0
    states = []
    for mode, chi_s in [("frustrated_kinetic", 0.1),
                        ("racemic_enhanced", 0.5),
                        ("self_sorting", 0.5)]:
        e = _env(mode, chi_struct=chi_s)
        A = assembly_competence(e, p)
        e["Estiff"] = target / A
        e["chi"] = 0.0
        states.append(run_to_steady_state(p, e))
    spread = float(np.max(np.abs(np.array(states) - states[0])))
    check("three different mechanisms at equal stiffness are indistinguishable",
          spread == 0.0, f"max spread = {spread:.2e}")


def test_majority_rules_mode():
    print("\nmajority-rules mode: minority absorbed, signal at the centre")
    from mirror_tenocyte.assembly import assembly_competence, net_helicity
    p = get_default_parameters()

    def A(c):
        return assembly_competence(
            {"assembly_mode": "majority_rules", "chi_struct": c}, p)

    check("A(0) = A(1) ~ 1 (pure forms assemble)",
          A(0.0) > 0.99 and abs(A(0.0) - A(1.0)) < 1e-12)
    check("symmetric about racemic", abs(A(0.2) - A(0.8)) < 1e-12)
    check("minimum at racemic", A(0.5) < A(0.2))
    check("chiral amplification: 10% D barely perturbs competence",
          A(0.1) > 0.95, f"A(0.1) = {A(0.1):.4f}")
    check("this is the OPPOSITE of capping, which collapses by 10% D",
          A(0.1) > 5 * assembly_competence(
              {"assembly_mode": "frustrated_kinetic", "chi_struct": 0.1}, p))
    check("net helicity saturates well before pure composition",
          net_helicity(0.35, p) > 0.9, f"m(chi=0.35) = {net_helicity(0.35, p):.3f}")


def test_two_point_protocol_beats_single_point():
    print("\ntwo compositions are needed to separate four hypotheses")
    from mirror_tenocyte.analysis import (
        two_point_protocol, discriminate_assembly_modes,
    )
    p = get_preset("clutch_calibrated")
    single = discriminate_assembly_modes(p=p)
    two = two_point_protocol(p=p)
    check("single-point discrimination is poor once four modes are in play",
          single["best_separation"] < 0.15,
          f"best single-point gap = {single['best_separation']:.4f}")
    check("two-point protocol is substantially better",
          two["best_separation"] > 2 * single["best_separation"],
          f"{two['best_separation']:.4f} vs {single['best_separation']:.4f}")
    check("optimal stiffness for the two-point protocol is soft",
          two["best_stiffness"] <= 3.0, f"{two['best_stiffness']:.0f} kPa")


def test_four_modes_have_distinct_patterns():
    print("\neach hypothesis has a distinct (low-chi, racemic) signature")
    from mirror_tenocyte.analysis import two_point_protocol, ASSEMBLY_MODES
    t = two_point_protocol(p=get_preset("clutch_calibrated"))
    sig = t["best_signatures"]
    lo, mid = sig["self_sorting"]
    check("self-sorting: flat, flat", abs(lo) < 1e-6 and abs(mid) < 1e-6)
    lo, mid = sig["frustrated_kinetic"]
    check("frustration: down at BOTH compositions", lo < -0.1 and mid < -0.1,
          f"({lo:+.3f}, {mid:+.3f})")
    lo, mid = sig["majority_rules"]
    check("majority rules: flat at low chi, down at racemic",
          abs(lo) < 0.05 and mid < -0.1, f"({lo:+.3f}, {mid:+.3f})")
    lo, mid = sig["racemic_enhanced"]
    check("enhancement: up at both", lo > 0.05 and mid > 0.05,
          f"({lo:+.3f}, {mid:+.3f})")


def test_experiment_feasibility():
    print("\nis the recommended experiment physically possible?")
    from mirror_tenocyte.feasibility import (
        synthesis_feasibility, cell_experiment_feasibility,
    )
    s = synthesis_feasibility()
    check("all-D collagen is beyond the D-protein synthesis record",
          not s["collagen_feasible"],
          f"1014 residues vs a {s['record']}-residue record")
    check("a D collagen-mimetic peptide is routine",
          s["cmp_feasible"], "36 residues")

    c = cell_experiment_feasibility()
    check("cell protocol does NOT work on native CMP gel stiffness",
          not c["cell_experiment_viable_natively"],
          f"best reachable separation {c['best_reachable']:.4f} at "
          f"{c['best_reachable_stiffness']:.1f} kPa")
    check("stiffness feasibility costs several-fold in separation",
          c["penalty"] > 3, f"{c['penalty']:.1f}x")
    print("      => the cell-free rheology series is the feasible experiment;")
    print("         the cell experiment needs a stiffened gel first.")


def test_protocol_numbers_are_reproducible():
    print("\nPROTOCOL.md classification table matches the model")
    from mirror_tenocyte.assembly import assembly_competence
    p = get_preset("clutch_calibrated")
    LAM_CMP = 600.0

    def norm(mode, chi, lam=None):
        q = dict(p)
        if lam is not None:
            q["lambda_c"] = lam
        ref = assembly_competence({"assembly_mode": mode, "chi_struct": 0.0}, q)
        val = assembly_competence({"assembly_mode": mode, "chi_struct": chi}, q)
        return val / ref

    # values quoted in PROTOCOL.md section 5
    expected = {
        ("self_sorting", 0.01): 1.00, ("self_sorting", 0.50): 1.00,
        ("frustrated_kinetic", 0.01): 0.03, ("frustrated_kinetic", 0.50): 0.00,
        ("majority_rules", 0.01): 1.00, ("majority_rules", 0.50): 0.30,
        ("racemic_enhanced", 0.01): 1.12, ("racemic_enhanced", 0.50): 4.00,
    }
    ok = True
    for (mode, chi), want in expected.items():
        lam = LAM_CMP if mode == "frustrated_kinetic" else None
        got = norm(mode, chi, lam)
        if abs(got - want) > 0.005:
            check(f"PROTOCOL {mode} at chi={chi}", False, f"{got:.3f} vs {want}")
            ok = False
    check("all eight quoted classification values reproduce", ok)

    check("capping collapse sits below 1% D at CMP geometry",
          norm("frustrated_kinetic", 0.01, LAM_CMP) < 0.05,
          f"G'(1%)/G'(0) = {norm('frustrated_kinetic', 0.01, LAM_CMP):.3f}")
    for mode in ("self_sorting", "frustrated_kinetic", "majority_rules",
                 "racemic_enhanced"):
        lam = LAM_CMP if mode == "frustrated_kinetic" else None
        if abs(norm(mode, 0.10, lam) - norm(mode, 0.90, lam)) > 1e-9:
            check(f"symmetry control holds for {mode}", False)
            return
    check("symmetry control G'(chi) = G'(1-chi) holds for every hypothesis", True,
          "so measured asymmetry is never an expected model behaviour")


def test_evidence_survey_is_well_formed():
    print("\npublished-evidence survey and the prior it yields")
    from mirror_tenocyte.evidence import SURVEY, prior, folds
    from mirror_tenocyte.analysis import ASSEMBLY_MODES

    for e in SURVEY:
        if e["outcome"] not in ASSEMBLY_MODES:
            check(f"outcome class valid for {e['system']}", False, e["outcome"])
            return
    check("every survey outcome maps to a modelled hypothesis", True,
          f"{len(SURVEY)} systems")

    pr = prior()
    check("prior reports counts, not fractions (n=6 is too small)",
          "fractions" not in pr)
    check("NO surveyed system is triple-helical — the key caveat",
          pr["n_triple_helical"] == 0)
    check("most entries are NOT verified from full text",
          pr["n_full_text"] < len(SURVEY) / 2,
          f"{pr['n_full_text']} of {pr['n_systems']}")
    check("almost none report a composition series",
          pr["n_with_composition_series"] <= 1,
          f"{pr['n_with_composition_series']} of {pr['n_systems']}")

    from mirror_tenocyte.evidence import fragility, verified_prior
    fr = fragility()
    check("the leading outcome is fragile to a single reclassification",
          fr["leader_flips_on_one_change"],
          f"{fr['n_entries_whose_change_flips_leader']} of {fr['n_entries']} entries can flip it")
    vp = verified_prior()
    check("restricted to verified entries the survey is NOT decisive",
          not vp["decisive"], f"{vp['counts']}")
    check("every entry carries a provenance level",
          all("verification" in e for e in SURVEY))

    from mirror_tenocyte.evidence import unverified_magnitudes, folds
    um = unverified_magnitudes()
    check("quoted effect magnitudes are flagged separately from direction",
          len(um) > 0, f"{len(um)} of {len(folds())} magnitudes unsupported")
    check("the 9x phenylalanine figure is marked unsupported",
          any("henylalanine" in s for s, _ in um),
          "its abstract confirms direction only, and contains no numbers")


def test_clutch_is_not_wired_into_model():
    print("\nclutch module stays out of the main model by default")
    import mirror_tenocyte.model as m
    src = open(m.__file__).read()
    check("model.py does not import the clutch module",
          "from .clutch" not in src and "import clutch" not in src)


def main():
    print("mirror-tenocyte test suite")
    test_states_bounded()
    test_monotonic_in_chirality()
    test_mirror_equals_nonadhesive()
    test_effective_ligand_modes()
    test_soft_substrate_kills_yap()
    test_product_degeneracy()
    test_presets_differ_in_bistability()
    test_chi_specific_channel()
    test_assembly_symmetry()
    test_assembly_mode_default_is_inert()
    test_self_sorting_reduces_to_parsimonious()
    test_degeneracy_breaks_only_under_frustration()
    test_rebound_signature()
    test_decoupled_scan_shape()
    test_clutch_physics()
    test_clutch_calibration()
    test_bond_type_sets_the_sign()
    test_hill_exponent_defaults_are_inert()
    test_clutch_calibrated_preset()
    test_predictions_survive_calibration()
    test_kinetic_competence_endpoints_and_symmetry()
    test_capping_is_sharper_than_quadratic()
    test_rebound_is_fragile_dip_is_robust()
    test_racemic_enhanced_mode()
    test_three_modes_give_three_signs()
    test_single_point_protocol_is_inadequate()
    test_endpoints_cannot_discriminate()
    test_heterochiral_triple_helix_is_forbidden()
    test_conclusions_do_not_depend_on_phenotype_weights()
    test_readout_ranking_depends_on_protocol()
    test_model_reduces_to_two_scalars()
    test_assembly_mode_is_invisible_downstream()
    test_majority_rules_mode()
    test_two_point_protocol_beats_single_point()
    test_four_modes_have_distinct_patterns()
    test_experiment_feasibility()
    test_protocol_numbers_are_reproducible()
    test_evidence_survey_is_well_formed()
    test_clutch_is_not_wired_into_model()

    print("\n" + "-" * 60)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): " + ", ".join(FAILURES))
        return 1
    print("All tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
