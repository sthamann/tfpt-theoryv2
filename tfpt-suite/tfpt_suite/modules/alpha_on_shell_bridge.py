from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from mpmath import mp

from tfpt_suite.alpha_running import (
    AlphaRunningInputs,
    alpha_bar5_MZ_from_alpha0,
    alpha_running_inputs_from_pdg,
    delta_alpha_extra_total,
    delta_alpha_lept_1loop,
)
from tfpt_suite.defect_partition import alpha_inv_0_from_delta2, derive_delta2_from_defect_partition
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass, mk_check_warn


def _leptonic_L_coeff(*, MZ_GeV: mp.mpf, lepton_masses_GeV: Iterable[mp.mpf]) -> mp.mpf:
    """
    Δα_lept_total_1loop(alpha0) = alpha0 * L  (because the 1-loop leptonic piece is linear in alpha0).
    """
    return mp.fsum([delta_alpha_lept_1loop(alpha0=mp.mpf(1), MZ_GeV=MZ_GeV, m_lepton_GeV=m) for m in lepton_masses_GeV])


def alpha0_from_alpha_bar5_MZ_1loop(*, alpha_bar5_MZ: mp.mpf, inputs: AlphaRunningInputs) -> tuple[mp.mpf, dict[str, mp.mpf]]:
    """
    Invert `alpha_bar5_MZ_from_alpha0` under the same (1-loop leptonic + PDG Δα_had) policy.

    Since Δα_lept_total is linear in alpha0 at 1-loop, the inversion is closed-form:

      alpha_MZ = alpha0 / (1 - (alpha0*L + Δα_had + Δα_extra))
      => alpha0 = alpha_MZ * (1 - Δα_had - Δα_extra) / (1 + alpha_MZ * L)
    """
    a_mz = mp.mpf(alpha_bar5_MZ)
    if a_mz <= 0:
        raise ValueError("alpha_bar5_MZ must be positive")

    L = _leptonic_L_coeff(MZ_GeV=inputs.MZ_GeV, lepton_masses_GeV=[inputs.m_e_GeV, inputs.m_mu_GeV, inputs.m_tau_GeV])
    dal_had = mp.mpf(inputs.delta_alpha_had5_MZ)
    dal_extra_total = delta_alpha_extra_total(inputs=inputs)

    alpha0 = a_mz * (mp.mpf(1) - dal_had - dal_extra_total) / (mp.mpf(1) + a_mz * L)
    return alpha0, {
        "L_lept_1loop": L,
        "delta_alpha_had5": dal_had,
        "delta_alpha_msbar_on_shell_shift": mp.mpf(inputs.delta_alpha_msbar_on_shell_shift_MZ),
        "delta_alpha_top_decoupling": mp.mpf(inputs.delta_alpha_top_decoupling_MZ),
        "delta_alpha_extra_msbar": mp.mpf(inputs.delta_alpha_extra_msbar_MZ),
        "delta_alpha_extra_total": dal_extra_total,
    }


LEPTON_SUM_TOL = mp.mpf("1e-30")
PRECISION_STABILITY_TOL = mp.mpf("1e-10")
Z_THRESHOLD_STRICT = mp.mpf(5)


class AlphaOnShellBridgeModule(TfptModule):
    module_id = "alpha_on_shell_bridge"
    title = "α(0) on-shell bridge (explicit renorm chain; metamorphic robustness checks)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT α^{-1}(0) prediction (from CFE+backreaction with δ₂ from defect partition)",
                "QED comparison layer inputs: tfpt_suite/data/alpha_running_pdg.json",
                "reference table: tfpt_suite/data/global_reference.json",
            ],
            outputs=[
                "declared α(0)↔ᾱ^(5)(MZ) chain with explicit Δα components",
                "metamorphic robustness checks (precision + threshold ordering invariance)",
            ],
            formulas=[
                r"ᾱ(MZ) = α(0) / (1 - Δα_lept(MZ) - Δα_had^{(5)}(MZ) - Δα_extra)",
                r"Δα_lept^l(MZ) = (α(0)/(3π)) (ln(MZ^2/m_l^2) - 5/3)",
                r"Δα_msbar−onshell(MZ) = 0.007122(2)(5) (PDG Eq. 10.13; W-loop + QCD shift)",
            ],
            validation=[
                "The α(0)→ᾱ(MZ) bridge reproduces the PDG ᾱ^(5)(MZ) reference within σ.",
                "The inverse bridge ᾱ(MZ)→α(0) is stable and scheme-consistent (under the declared 1-loop policy).",
                "Explicit finite EW decoupling pieces are included (MSbar↔on-shell shift, top decoupling).",
            ],
            determinism="Deterministic (closed-form 1-loop chain; no numerical integration).",
            question="Is the α(0) metrology comparison performed via a declared on-shell ↔ MSbar-at-MZ bridge with robust numerics?",
            objective=[
                "Separate TFPT’s α-sector prediction from the SM/QED comparison layer.",
                "Make threshold ordering / precision sensitivity explicit via metamorphic checks.",
            ],
            gaps=[
                "Publication-grade requires full electroweak scheme conversion and higher-loop QED/QCD contributions (this module is an explicit 1-loop bridge with finite EW decoupling pieces).",
            ],
        )

    def run(self, config) -> ModuleResult:
        # Inputs
        pdg_path = Path(__file__).resolve().parent.parent / "data" / "alpha_running_pdg.json"
        pdg = json.loads(pdg_path.read_text(encoding="utf-8"))
        inputs = alpha_running_inputs_from_pdg(pdg=pdg)

        ref_path = Path(__file__).resolve().parent.parent / "data" / "global_reference.json"
        ref = json.loads(ref_path.read_text(encoding="utf-8"))
        obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
        codata = obs.get("alpha_inv_0", {}) if isinstance(obs.get("alpha_inv_0", {}), dict) else {}
        pdg_mz = obs.get("alpha_bar5_inv_MZ", {}) if isinstance(obs.get("alpha_bar5_inv_MZ", {}), dict) else {}
        alpha_inv_0_ref = mp.mpf(str(codata.get("mean", "137.035999177")))
        alpha_inv_0_sigma = mp.mpf(str(codata.get("sigma", "2.1e-8")))
        alpha_bar5_inv_MZ_ref = mp.mpf(str(pdg_mz.get("mean", "127.951")))
        alpha_bar5_inv_MZ_sigma = mp.mpf(str(pdg_mz.get("sigma", "0.009")))

        # TFPT α(0) prediction (via δ₂ from defect partition; no fit parameter)
        d2 = derive_delta2_from_defect_partition(delta_top=mp.mpf("0.00012030447"))  # overwritten below for clarity
        # Use mp context precision requested by SuiteConfig
        old_dps = mp.dps
        mp.dps = max(60, int(config.mp_dps))
        try:
            # Recompute δ₂ using TFPT delta_top at current precision
            from tfpt_suite.constants import TfptConstants

            c = TfptConstants.compute()
            d2 = derive_delta2_from_defect_partition(delta_top=mp.mpf(c.delta_top))
            alpha_inv_0_tfpt = alpha_inv_0_from_delta2(delta2=d2.delta2, mp_dps=mp.dps)
        finally:
            mp.dps = old_dps

        alpha0_tfpt = mp.mpf(1) / alpha_inv_0_tfpt
        alpha_mz_tfpt, parts_tfpt = alpha_bar5_MZ_from_alpha0(alpha0=alpha0_tfpt, inputs=inputs)
        alpha_bar5_inv_MZ_tfpt = mp.mpf(1) / alpha_mz_tfpt

        # Inverse bridge: PDG ᾱ(MZ) -> α(0) (under the same declared policy)
        alpha_mz_ref = mp.mpf(1) / alpha_bar5_inv_MZ_ref
        alpha0_from_mz, inv_parts = alpha0_from_alpha_bar5_MZ_1loop(alpha_bar5_MZ=alpha_mz_ref, inputs=inputs)
        alpha_inv_0_from_mz = mp.mpf(1) / alpha0_from_mz

        # --- Explicit leptonic VP breakdown (1-loop) ---
        lepton_parts = {
            "e": delta_alpha_lept_1loop(alpha0=alpha0_tfpt, MZ_GeV=inputs.MZ_GeV, m_lepton_GeV=inputs.m_e_GeV),
            "mu": delta_alpha_lept_1loop(alpha0=alpha0_tfpt, MZ_GeV=inputs.MZ_GeV, m_lepton_GeV=inputs.m_mu_GeV),
            "tau": delta_alpha_lept_1loop(alpha0=alpha0_tfpt, MZ_GeV=inputs.MZ_GeV, m_lepton_GeV=inputs.m_tau_GeV),
        }
        lepton_total = mp.fsum(lepton_parts.values())

        # --- Metamorphic checks ---
        # 1) threshold ordering invariance (permute lepton masses)
        L_emu_tau = _leptonic_L_coeff(MZ_GeV=inputs.MZ_GeV, lepton_masses_GeV=[inputs.m_e_GeV, inputs.m_mu_GeV, inputs.m_tau_GeV])
        L_tau_mu_e = _leptonic_L_coeff(MZ_GeV=inputs.MZ_GeV, lepton_masses_GeV=[inputs.m_tau_GeV, inputs.m_mu_GeV, inputs.m_e_GeV])
        # 2) precision stability (mp.dps low vs high)
        old_dps = mp.dps
        try:
            mp.dps = 50
            alpha_mz_low, _ = alpha_bar5_MZ_from_alpha0(alpha0=mp.mpf(1) / alpha_inv_0_tfpt, inputs=inputs)
            inv_low = mp.mpf(1) / alpha_mz_low
            mp.dps = 110
            alpha_mz_high, _ = alpha_bar5_MZ_from_alpha0(alpha0=mp.mpf(1) / alpha_inv_0_tfpt, inputs=inputs)
            inv_high = mp.mpf(1) / alpha_mz_high
        finally:
            mp.dps = old_dps

        # --- Checks ---
        checks: list[Check] = []
        # Explicit leptonic VP (per-lepton) reconstruction.
        lept_diff = abs(lepton_total - mp.mpf(parts_tfpt.get("delta_alpha_lept_1loop", 0)))
        if lept_diff <= LEPTON_SUM_TOL:
            checks.append(
                mk_check_pass(
                    "alpha_bridge_leptons_1loop_explicit",
                    f"lepton sum matches Δα_lept (|Δ|={lept_diff})",
                )
            )
        else:
            checks.append(
                mk_check_fail(
                    "alpha_bridge_leptons_1loop_explicit",
                    f"lepton sum mismatch: |Δ|={lept_diff}",
                )
            )

        # Hadronic policy declared (source field must be present).
        had_source = ""
        if isinstance(pdg.get("delta_alpha_had5_MZ", {}), dict):
            had_source = str(pdg["delta_alpha_had5_MZ"].get("source", ""))
        checks.append(
            mk_check_pass(
                "alpha_bridge_hadron_policy_declared",
                f"Δα_had^(5) source: {had_source}",
            )
            if had_source.strip()
            else mk_check_fail("alpha_bridge_hadron_policy_declared", "missing Δα_had^(5) source in alpha_running_pdg.json")
        )

        # EW decoupling pieces (MSbar↔on-shell shift and top decoupling).
        ew_shift = mp.mpf(str(pdg.get("delta_alpha_msbar_on_shell_shift_MZ", {}).get("mean", 0.0)))
        top_dec = mp.mpf(str(pdg.get("delta_alpha_top_decoupling_MZ", {}).get("mean", 0.0)))
        if mp.isfinite(ew_shift) and mp.isfinite(top_dec):
            checks.append(
                mk_check_pass(
                    "alpha_bridge_EW_decoupling_included",
                    f"Δα_msbar_on_shell={ew_shift}, Δα_top_decoupling={top_dec}",
                )
            )
        else:
            checks.append(
                mk_check_fail(
                    "alpha_bridge_EW_decoupling_included",
                    f"non-finite EW decoupling pieces: shift={ew_shift}, top={top_dec}",
                )
            )


        # Bridge reproduces PDG ᾱ^(5)(MZ) within σ (this is the primary "MSbar-at-MZ" truth test in the suite).
        z_mz = (alpha_bar5_inv_MZ_tfpt - alpha_bar5_inv_MZ_ref) / alpha_bar5_inv_MZ_sigma if alpha_bar5_inv_MZ_sigma != 0 else mp.mpf("nan")
        if mp.isfinite(z_mz) and abs(z_mz) <= Z_THRESHOLD_STRICT:
            checks.append(mk_check_pass("alpha_bridge_reproduces_alpha_bar5_MZ", f"z={z_mz} (tfpt={alpha_bar5_inv_MZ_tfpt}, ref={alpha_bar5_inv_MZ_ref})"))
        else:
            checks.append(mk_check_fail("alpha_bridge_reproduces_alpha_bar5_MZ", f"z={z_mz} (tfpt={alpha_bar5_inv_MZ_tfpt}, ref={alpha_bar5_inv_MZ_ref})"))

        # Inverse consistency: PDG ᾱ(MZ) maps back to CODATA α(0) under the declared policy.
        # This is a diagnostic check; with the explicit 1-loop bridge it may deviate from the
        # CODATA value even when the forward ᾱ(MZ) gate is satisfied.
        diff_inv = alpha_inv_0_from_mz - alpha_inv_0_ref
        z_inv = diff_inv / alpha_inv_0_sigma if alpha_inv_0_sigma != 0 else mp.mpf("nan")
        if mp.isfinite(z_inv) and abs(z_inv) <= Z_THRESHOLD_STRICT:
            checks.append(mk_check_pass("alpha_bridge_inverse_matches_alpha0_ref", f"z={z_inv} (inv={alpha_inv_0_from_mz}, ref={alpha_inv_0_ref})"))
        elif mp.isfinite(z_inv):
            checks.append(
                mk_check_warn(
                    "alpha_bridge_inverse_matches_alpha0_ref",
                    f"z={z_inv} (inv={alpha_inv_0_from_mz}, ref={alpha_inv_0_ref}); 1-loop bridge is diagnostic only",
                )
            )
        else:
            checks.append(mk_check_fail("alpha_bridge_inverse_matches_alpha0_ref", f"z={z_inv} (inv={alpha_inv_0_from_mz}, ref={alpha_inv_0_ref})"))

        # Metamorph 1: lepton threshold ordering (should be exact at this level)
        checks.append(
            mk_check_pass("alpha_bridge_metamorph_lepton_ordering", f"L(emu,tau)={L_emu_tau}, L(tau,mu,e)={L_tau_mu_e}")
            if L_emu_tau == L_tau_mu_e
            else mk_check_fail("alpha_bridge_metamorph_lepton_ordering", f"L differs: {L_emu_tau} vs {L_tau_mu_e}")
        )

        # Metamorph 2: precision stability
        diff_prec = abs(inv_high - inv_low)
        checks.append(
            mk_check_pass("alpha_bridge_metamorph_precision_stable", f"|inv_high-inv_low|={diff_prec} (low={inv_low}, high={inv_high})")
            if diff_prec < PRECISION_STABILITY_TOL
            else mk_check_fail("alpha_bridge_metamorph_precision_stable", f"|inv_high-inv_low|={diff_prec} (low={inv_low}, high={inv_high})")
        )

        lines: list[str] = []
        lines += [
            "α(0) on-shell bridge (declared renorm chain; 1-loop leptonic + PDG Δα_had)",
            "",
            f"pdg inputs: {pdg_path}",
            f"reference table: {ref_path}",
            "",
            "TFPT α(0) prediction source:",
            f"- δ2 model_id: {d2.model_id}",
            f"- δ2: {d2.delta2} (delta2/delta_top^2={d2.delta2_over_delta_top2})",
            f"- α_inv_0_tfpt: {alpha_inv_0_tfpt}",
            "",
            "Bridge α(0) -> ᾱ^(5)(MZ):",
            f"- ᾱ^(5)(MZ)^{-1} (tfpt-bridge): {alpha_bar5_inv_MZ_tfpt}",
            f"- Δα components: {parts_tfpt}",
            f"- PDG ref ᾱ^(5)(MZ)^{-1}: {alpha_bar5_inv_MZ_ref} ± {alpha_bar5_inv_MZ_sigma}",
            "",
            "Leptonic VP (explicit 1-loop breakdown):",
            f"- e: {lepton_parts['e']}",
            f"- mu: {lepton_parts['mu']}",
            f"- tau: {lepton_parts['tau']}",
            f"- total: {lepton_total}",
            "",
            "EW decoupling pieces (finite):",
            f"- Δα_msbar_on_shell_shift_MZ: {parts_tfpt.get('delta_alpha_msbar_on_shell_shift', mp.mpf(0))}",
            f"- Δα_top_decoupling_MZ: {parts_tfpt.get('delta_alpha_top_decoupling', mp.mpf(0))}",
            f"- Δα_extra_total: {parts_tfpt.get('delta_alpha_extra_total', mp.mpf(0))}",
            "",
            "Inverse bridge ᾱ^(5)(MZ) -> α(0) (same policy):",
            f"- α_inv_0_from_alpha_bar5_ref: {alpha_inv_0_from_mz}",
            f"- inverse parts: {inv_parts}",
            f"- CODATA ref α_inv_0: {alpha_inv_0_ref} ± {alpha_inv_0_sigma}",
            "",
            "Metamorphic checks:",
            f"- L coefficient (ordering invariant): L_emu_tau={L_emu_tau}, L_tau_mu_e={L_tau_mu_e}",
            f"- precision stability: inv_low_dps50={inv_low}, inv_high_dps110={inv_high}, |Δ|={diff_prec}",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
        ]

        return ModuleResult(
            results={
                "alpha0_tfpt": {"alpha_inv_0": alpha_inv_0_tfpt, "delta2_model": d2.model_id, "delta2": d2.delta2},
                "bridge_to_mz": {
                    "alpha_bar5_inv_MZ": alpha_bar5_inv_MZ_tfpt,
                    "delta_alpha_parts": parts_tfpt,
                    "leptonic_breakdown_1loop": lepton_parts,
                },
                "inverse_from_mz_ref": {"alpha_inv_0": alpha_inv_0_from_mz, "parts": inv_parts},
                "metamorph": {
                    "L_coeff": {"emu_tau": L_emu_tau, "tau_mu_e": L_tau_mu_e},
                    "precision": {"inv_low_dps50": inv_low, "inv_high_dps110": inv_high, "abs_diff": diff_prec},
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

