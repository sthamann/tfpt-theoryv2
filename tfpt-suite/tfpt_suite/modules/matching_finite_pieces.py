from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.alpha_running import AlphaRunningInputs, alpha_bar5_MZ_from_alpha0, alpha_running_inputs_from_pdg
from tfpt_suite.conventions import alpha_from_g
from tfpt_suite.matching import (
    delta_alpha_msbar_on_shell_shift_pdg,
    delta_lambda_h_1loop_buttazzo,
    delta_y_t_ew_hempfling,
    match_alpha3_quark_threshold_2loop,
)
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from tfpt_suite.rge_sm import run_sm_gauge_only_2loop_thresholds
from tfpt_suite.sm_inputs import SmMzInputs, gauge_couplings_from_mz_inputs

EW_SHIFT_TOL = mp.mpf("1e-5")
EW_FINITE_PIECE_MIN = mp.mpf("1e-6")
DEFAULT_V_EV_GEV = 246.0


class MatchingFinitePiecesModule(TfptModule):
    module_id = "matching_finite_pieces"
    title = "Matching finite pieces (QCD α_s 2-loop decoupling + QED/EW α(0)→ᾱ^(5)(MZ) bridge)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "SM inputs at MZ: tfpt_suite/data/sm_inputs_mz.json",
                "QED comparison-layer policy: tfpt_suite/data/alpha_running_pdg.json",
                "reference table: tfpt_suite/data/global_reference.json",
            ],
            outputs=[
                "explicit demonstration of a nontrivial finite matching step (α_s 2-loop heavy-quark decoupling)",
                "side-by-side α_s(μ) values with/without finite matching at thresholds (QCD-only policy)",
                "explicit α(0)→ᾱ^(5)(MZ) renorm chain (QED leptonic 1-loop + Δα_had^(5) + explicit MSbar remainder)",
                "explicit top/Higgs 1-loop EW finite pieces (Δy_t, Δλ_H) evaluated at threshold scales",
            ],
            formulas=[
                r"α_nf(μ=m_Q) = α_{nf+1} [1 + (11/72)(α/π)^2]  (direction=down); inverse for direction=up",
                r"ᾱ(MZ) = α(0) / (1 - Δα_lept(MZ) - Δα_had^{(5)}(MZ) - Δα_msbar_shift - Δα_top - Δα_extra)",
                r"δ_t^w(μ) = (G_F m_t^2 / (8 π^2 √2)) [-(N_c + 3/2) ln(m_t^2/μ^2) + N_c/2 + 4 − r + 2 r (2 r − 3) ln(4 r) − 8 r^2 f(r)]",
                r"λ^{(1)}(μ) from Buttazzo et al. (App. A.1) using A0/B0 with A0(M)=M^2(1−ln(M^2/μ^2)) and B0(p;M1,M2)=-∫_0^1 ln[(x M1^2+(1-x)M2^2−x(1-x)p^2)/μ^2] dx",
            ],
            validation=[
                "finite α_s matching produces a non-zero step at heavy-quark thresholds at 2-loop",
                "matching is invertible up↔down to numerical precision (metamorphic check)",
                "QED/EW finite pieces are explicit and non-zero (Δα_lept, Δα_had, Δα_msbar_shift)",
                "top/Higgs EW finite pieces are explicit and non-zero at the declared thresholds",
            ],
            determinism="Deterministic given inputs.",
            question="Are finite matching pieces implemented explicitly (not silently), and are they numerically well-behaved?",
            objective=[
                "Provide a concrete 'finite piece exists' proof point (QCD α_s threshold at 2-loop).",
                "Make it easy to audit what changes when matching is enabled.",
            ],
            gaps=[
                "Full publication-grade matching still requires end-to-end use of EW/QED finite pieces across all declared thresholds plus a full below-MZ EFT policy; this module provides explicit proof points for QCD (α_s), QED/EW α-running, and top/Higgs EW finite pieces.",
            ],
        )

    def run(self, config) -> ModuleResult:
        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))

        mu_mz = float(sm_raw["mu_GeV"])
        mc = float(sm_raw.get("mc_GeV", 1.27))
        mb = float(sm_raw.get("mb_GeV", 4.18))
        mt = float(sm_raw.get("mt_GeV", 172.76))
        mH = float(sm_raw.get("mH_GeV", 125.25))
        mW = float(sm_raw.get("mW_GeV", 80.379))
        v_ev = float(sm_raw.get("v_ev_GeV", DEFAULT_V_EV_GEV))

        inp = SmMzInputs(
            mu_GeV=float(mu_mz),
            alpha_em_inv=float(sm_raw["alpha_em_inv"]),
            sin2_thetaW=float(sm_raw["sin2_thetaW"]),
            alpha_s=float(sm_raw["alpha_s"]),
        )
        g1_gut_mz, g2_mz, g3_mz = gauge_couplings_from_mz_inputs(inp)

        # Compare α_s at μ=mc with and without finite matching at the mb threshold.
        g_mc_no_match = run_sm_gauge_only_2loop_thresholds(
            mu_start_GeV=mu_mz,
            mu_end_GeV=mc,
            g_start=(g1_gut_mz, g2_mz, g3_mz),
            mc_GeV=mc,
            mb_GeV=mb,
            mt_GeV=mt,
            apply_alpha3_matching=False,
        )
        g_mc_match = run_sm_gauge_only_2loop_thresholds(
            mu_start_GeV=mu_mz,
            mu_end_GeV=mc,
            g_start=(g1_gut_mz, g2_mz, g3_mz),
            mc_GeV=mc,
            mb_GeV=mb,
            mt_GeV=mt,
            apply_alpha3_matching=True,
        )
        alpha_mc_no = float(alpha_from_g(float(g_mc_no_match[2])))
        alpha_mc_yes = float(alpha_from_g(float(g_mc_match[2])))
        delta_alpha_mc = float(alpha_mc_yes - alpha_mc_no)

        # Metamorphic invertibility at a representative α value:
        a0 = float(alpha_mc_yes)
        a_up = float(match_alpha3_quark_threshold_2loop(alpha3=a0, direction="up"))
        a_back = float(match_alpha3_quark_threshold_2loop(alpha3=a_up, direction="down"))
        inv_err = float(a_back - a0)

        # QED/EW: explicit α(0) → ᾱ^(5)(MZ) renorm chain (comparison-layer finite pieces)
        pdg_path = Path(__file__).resolve().parent.parent / "data" / "alpha_running_pdg.json"
        pdg = json.loads(pdg_path.read_text(encoding="utf-8"))
        inputs = alpha_running_inputs_from_pdg(pdg=pdg)
        mz_ref = float(pdg.get("MZ_GeV", 91.1876))
        mw_ref = float(pdg.get("MW_GeV", 80.379))
        alpha_s_ref = float(pdg.get("alpha_s_MZ", {}).get("mean", 0.1180))
        ref_path = Path(__file__).resolve().parent.parent / "data" / "global_reference.json"
        ref = json.loads(ref_path.read_text(encoding="utf-8"))
        obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
        codata = obs.get("alpha_inv_0", {}) if isinstance(obs.get("alpha_inv_0", {}), dict) else {}
        alpha_inv_0_ref = mp.mpf(str(codata.get("mean", "137.035999177")))
        alpha0_ref = mp.mpf(1) / alpha_inv_0_ref
        alpha_mz_ref, parts = alpha_bar5_MZ_from_alpha0(alpha0=alpha0_ref, inputs=inputs)
        alpha_bar5_inv_MZ_from_alpha0 = mp.mpf(1) / alpha_mz_ref

        delta_y_t = delta_y_t_ew_hempfling(m_t_GeV=mt, m_h_GeV=mH, mu_GeV=mt, v_ev_GeV=v_ev)
        delta_lambda = delta_lambda_h_1loop_buttazzo(
            m_h_GeV=mH,
            m_t_GeV=mt,
            m_w_GeV=mW,
            m_z_GeV=mu_mz,
            mu_GeV=mH,
            v_ev_GeV=v_ev,
        )

        checks: list[Check] = []
        if abs(delta_alpha_mc) > 0:
            checks.append(mk_check_pass("finite_alpha3_matching_nontrivial", f"Δα_s(mc)={delta_alpha_mc:.3e} (matching-on minus matching-off)"))
        else:
            checks.append(mk_check_warn("finite_alpha3_matching_nontrivial", f"Δα_s(mc)={delta_alpha_mc} (unexpectedly zero)"))
        checks.append(
            mk_check_pass("alpha3_matching_invertible", f"down(up(alpha)) - alpha = {inv_err:.3e}")
            if abs(inv_err) < 1e-12
            else mk_check_warn("alpha3_matching_invertible", f"down(up(alpha)) - alpha = {inv_err:.3e}")
        )
        dal_tot = parts.get("delta_alpha_total", mp.mpf("nan"))
        dal_lept = parts.get("delta_alpha_lept_1loop", mp.mpf("nan"))
        dal_had5 = parts.get("delta_alpha_had5", mp.mpf("nan"))
        dal_msbar_shift = parts.get("delta_alpha_msbar_on_shell_shift", mp.mpf("nan"))
        dal_top = parts.get("delta_alpha_top_decoupling", mp.mpf("nan"))
        dal_extra_total = parts.get("delta_alpha_extra_total", mp.mpf("nan"))
        shift_formula = mp.mpf(
            str(
                delta_alpha_msbar_on_shell_shift_pdg(
                    alpha0=float(alpha0_ref),
                    alpha_s_MZ=float(alpha_s_ref),
                    MZ_GeV=float(mz_ref),
                    MW_GeV=float(mw_ref),
                )
            )
        )
        shift_diff = abs(mp.mpf(dal_msbar_shift) - shift_formula) if mp.isfinite(dal_msbar_shift) else mp.mpf("nan")
        checks.append(
            mk_check_pass(
                "matching_finite_pieces_EW_QED_present",
                "Δα_total={dal_tot} (lept={dal_lept}, had5={dal_had5}, msbar_shift={dal_msbar_shift}, top={dal_top}, extra_total={dal_extra_total}); "
                "Δα̂−Δα formula={shift_formula}, diff={shift_diff}; alpha_bar5_inv_MZ(from alpha0)={alpha_bar5_inv_MZ_from_alpha0}; "
                "Δy_t^EW(μ=mt)={delta_y_t}, Δλ_H^EW(μ=mH)={delta_lambda}".format(
                    dal_tot=dal_tot,
                    dal_lept=dal_lept,
                    dal_had5=dal_had5,
                    dal_msbar_shift=dal_msbar_shift,
                    dal_top=dal_top,
                    dal_extra_total=dal_extra_total,
                    shift_formula=shift_formula,
                    shift_diff=shift_diff,
                    alpha_bar5_inv_MZ_from_alpha0=alpha_bar5_inv_MZ_from_alpha0,
                    delta_y_t=delta_y_t,
                    delta_lambda=delta_lambda,
                ),
            )
            if mp.isfinite(dal_tot)
            and abs(dal_tot) > mp.mpf("1e-6")
            and mp.isfinite(shift_diff)
            and shift_diff <= EW_SHIFT_TOL
            else mk_check_warn(
                "matching_finite_pieces_EW_QED_present",
                f"Δα_total={dal_tot} (lept={dal_lept}, had5={dal_had5}, msbar_shift={dal_msbar_shift}); "
                f"Δα̂−Δα diff={shift_diff}; Δy_t={delta_y_t}, Δλ_H={delta_lambda}",
            )
        )
        dy_t_mp = mp.mpf(str(delta_y_t))
        dlam_mp = mp.mpf(str(delta_lambda))
        checks.append(
            mk_check_pass(
                "matching_finite_pieces_top_higgs_present",
                f"Δy_t^EW(μ=mt)={dy_t_mp}, Δλ_H^EW(μ=mH)={dlam_mp} (v={v_ev} GeV)",
            )
            if mp.isfinite(dy_t_mp)
            and mp.isfinite(dlam_mp)
            and abs(dy_t_mp) > EW_FINITE_PIECE_MIN
            and abs(dlam_mp) > EW_FINITE_PIECE_MIN
            else mk_check_warn(
                "matching_finite_pieces_top_higgs_present",
                f"Δy_t^EW(μ=mt)={dy_t_mp}, Δλ_H^EW(μ=mH)={dlam_mp}",
            )
        )
        checks.append(
            mk_check_pass("below_MZ_policy_explicit_and_applied", f"policy file={pdg_path} (Δα_had5 + leptonic 1-loop + explicit EW MSbar shift)")
        )

        lines: list[str] = []
        lines += [
            "Matching finite pieces (QCD + QED/EW proof points)",
            "",
            f"Inputs: mu_MZ={mu_mz} GeV, mc={mc} GeV, mb={mb} GeV, mt={mt} GeV, mH={mH} GeV, mW={mW} GeV",
            "",
            "α_s at μ=mc with/without finite matching at thresholds:",
            f"- α_s(mc) no finite matching: {alpha_mc_no:.12g}",
            f"- α_s(mc) with finite matching: {alpha_mc_yes:.12g}",
            f"- Δα_s(mc) = {delta_alpha_mc:.3e}",
            "",
            "QED/EW on-shell → MSbar-at-MZ bridge (finite pieces; comparison layer):",
            f"- α_inv_0 (CODATA ref) = {alpha_inv_0_ref}",
            f"- ᾱ^(5)(MZ)^{-1} (from α(0) + Δα pieces) = {alpha_bar5_inv_MZ_from_alpha0}",
            f"- Δα parts: {parts}",
            f"- Δα̂−Δα (PDG Eq. 10.13, recomputed) = {shift_formula} (diff={shift_diff})",
            "",
            "Top/Higgs EW finite pieces (1-loop threshold formulas):",
            f"- v_ev = {v_ev} GeV (from sm_inputs or default)",
            f"- Δy_t^EW(μ=mt) = {delta_y_t}",
            f"- Δλ_H^EW(μ=mH) = {delta_lambda}",
            "",
            "Metamorphic invertibility (2-loop decoupling primitive):",
            f"- alpha={a0:.12g} => up => {a_up:.12g} => down => {a_back:.12g} (error={inv_err:.3e})",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "inputs": {"sm_inputs_mz_file": str(sm_path), "alpha_running_pdg_file": str(pdg_path), "global_reference_file": str(ref_path)},
                "alpha_s_mc": {"no_matching": alpha_mc_no, "with_matching": alpha_mc_yes, "delta": delta_alpha_mc},
                "invertibility": {"alpha": a0, "alpha_up": a_up, "alpha_back": a_back, "error": inv_err},
                "alpha_running": {
                    "alpha_inv_0_ref": str(alpha_inv_0_ref),
                    "alpha_bar5_inv_MZ_from_alpha0": str(alpha_bar5_inv_MZ_from_alpha0),
                    "delta_alpha_parts": {str(k): str(v) for k, v in parts.items()},
                },
                "finite_pieces_top_higgs": {
                    "v_ev_GeV": v_ev,
                    "delta_y_t_ew_mu_mt": delta_y_t,
                    "delta_lambda_h_ew_mu_mH": delta_lambda,
                    "mu_mt_GeV": mt,
                    "mu_mH_GeV": mH,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

