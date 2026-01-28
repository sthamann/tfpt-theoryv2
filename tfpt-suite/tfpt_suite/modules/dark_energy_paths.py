from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.cosmo_scale_map import MPL_REDUCED_GEV
from tfpt_suite.e8_ladder import e8_phi_n
from tfpt_suite.module_base import (
    Check,
    ModuleResult,
    ModuleSpec,
    TfptModule,
    mk_check_info,
    mk_check_pass,
    mk_check_warn,
)


def _H0_GeV_from_km_s_Mpc(H0_km_s_Mpc: float) -> float:
    # H0 [s^-1] = H0_km_s_Mpc / (Mpc in km); 1/s = ħ [GeV·s] in natural units.
    Mpc_km = 3.0856775814913673e19
    hbar_GeV_s = 6.582119569e-25
    H0_s_inv = float(H0_km_s_Mpc) / Mpc_km
    return float(H0_s_inv * hbar_GeV_s)


LADDER_N_MAX = 26
LADDER_D1 = mp.mpf(58)
LADDER_D0 = mp.mpf(60)
TERMINAL_STAGE_EXPECTED = mp.mpf(30)
TERMINAL_STAGE_TOL = mp.mpf("0.75")


class DarkEnergyPathsModule(TfptModule):
    module_id = "dark_energy_paths"
    title = "Dark energy (Λ) — target values for UFE torsion condensate and cascade terminal-stage paths"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "cosmology reference: tfpt_suite/data/k_calibration.json (Planck-style flat ΛCDM parameters used elsewhere in suite)",
                "theory note: five_problems.tex Sec. 5 (two paths: torsion condensate and E8 terminal stage)",
                "TFPT invariants: tfpt_suite/constants.py (E8 ladder parameters)",
                "torsion condensate prediction (optional): tfpt-suite/out*/torsion_condensate/results.json when available in output_dir",
            ],
            outputs=[
                "Λ_obs (mass^2) from H0 and Ω_Λ",
                "ρ_Λ (energy density) from ρ_c and Ω_Λ",
                "Target ⟨K^2⟩ and K_rms implied by Λ_eff = (1/4)⟨K^2⟩ (UFE path A)",
                "Target φ_* implied by ρ_Λ ≈ (M̄_P φ_*)^4 (cascade magnitude path B)",
                "E8 ladder sequence φ_n and terminal-stage identification (n≈30 extrapolation)",
                "best-effort comparison to the in-suite torsion-condensate prediction when available",
            ],
            formulas=[
                "Ω_Λ = 1 − Ω_m − Ω_r (flat)",
                "H0[GeV] from H0[km/s/Mpc] via ħ conversion",
                "Λ_obs = 3 Ω_Λ H0^2  (ΛCDM)",
                "ρ_c = 3 H0^2 M̄_P^2,  ρ_Λ = Ω_Λ ρ_c",
                "UFE path A (five_problems): Λ_eff = (1/4)⟨K^2⟩ ⇒ ⟨K^2⟩_target = 4 Λ_obs",
                "Cascade path B (five_problems): ρ_Λ ≈ (M̄_P φ_*)^4 ⇒ φ_* = (ρ_Λ)^{1/4}/M̄_P",
                "E8 ladder: φ_n = φ_0 e^{-γ(0)} (D_n/D_1)^λ, D_n=60-2n",
                "Terminal stage: solve φ_n = φ_* for n (continuous) to identify n_terminal≈30",
            ],
            validation=[
                "Produces finite Λ_obs, ρ_Λ and corresponding target values under the declared cosmology.",
            ],
            determinism="Deterministic given the input cosmology table.",
            question="If TFPT attributes dark energy either to a torsion condensate (UFE) or to a terminal cascade stage, what target magnitudes must those mechanisms reproduce?",
            objective=[
                "Turn 'dark energy path A/B' into explicit numerical targets (engineering constraints).",
                "Expose a hard interface to the prediction module (`torsion_condensate`): targets here, derived candidate there.",
            ],
            what_was_done=[
                "Computed Λ_obs and ρ_Λ from the suite’s Planck-style ΛCDM parameter snapshot.",
                "Translated the UFE identity Λ_eff=(1/4)⟨K^2⟩ into a target K_rms scale.",
                "Translated the cascade magnitude ansatz ρ_Λ≈(M̄_P φ_*)^4 into a target φ_* suppression factor.",
            ],
            assumptions=[
                "Use flat ΛCDM bookkeeping (Ω_Λ=1−Ω_m−Ω_r).",
                "Interpret five_problems.tex path A as a statement about the cosmological constant Λ (mass^2), not directly ρ_Λ; ρ_Λ is derived via ρ_Λ=Λ M̄_P^2.",
                "Interpret five_problems.tex path B as a magnitude relation ρ_Λ≈(M̄_P φ_*)^4 with an unspecified ladder-derived φ_*.",
                "Terminal stage identification uses a continuous extrapolation of the E8 ladder to n≈30 (D_n→0).",
            ],
            gaps=[
                "Terminal-stage extrapolation is continuous; a discrete microscopic derivation of the n≈30 endpoint remains future work.",
                "Publication-grade microscopic derivation of ⟨K^2⟩ remains future work; the suite currently ships a discrete defect-suppression candidate in `torsion_condensate`.",
            ],
            references=[
                "five_problems.tex Sec. 5 (two TFPT-consistent paths for dark energy)",
                "paper_v1_06_01_09_2025.tex Sec. 9 / Appendix C (mentions n=30 ladder stage for ρ_Λ at order-of-magnitude level)",
            ],
            maturity="target-setting module (makes Λ requirements explicit; does not claim a derived prediction yet)",
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        data_dir = Path(__file__).resolve().parent.parent / "data"
        kcal_path = data_dir / "k_calibration.json"
        kcal = json.loads(kcal_path.read_text(encoding="utf-8")) if kcal_path.is_file() else {}
        cosmo_raw = kcal.get("cosmology_flat_lcdm", {}) if isinstance(kcal.get("cosmology_flat_lcdm", {}), dict) else {}

        H0_km_s_Mpc = float(cosmo_raw.get("H0_km_s_Mpc", float("nan")))
        Omega_m = float(cosmo_raw.get("Omega_m", float("nan")))
        Omega_r = float(cosmo_raw.get("Omega_r", 0.0))
        Omega_L = float(1.0 - Omega_m - Omega_r) if (math.isfinite(Omega_m) and math.isfinite(Omega_r)) else float("nan")

        H0_GeV = _H0_GeV_from_km_s_Mpc(H0_km_s_Mpc) if math.isfinite(H0_km_s_Mpc) else float("nan")
        Lambda_obs_GeV2 = float(3.0 * Omega_L * (H0_GeV**2)) if (math.isfinite(Omega_L) and math.isfinite(H0_GeV)) else float("nan")

        Mpl = float(MPL_REDUCED_GEV)
        rho_c_GeV4 = float(3.0 * (H0_GeV**2) * (Mpl**2)) if math.isfinite(H0_GeV) else float("nan")
        rho_L_GeV4 = float(Omega_L * rho_c_GeV4) if (math.isfinite(Omega_L) and math.isfinite(rho_c_GeV4)) else float("nan")

        # UFE path A targets
        K2_target_GeV2 = float(4.0 * Lambda_obs_GeV2) if math.isfinite(Lambda_obs_GeV2) else float("nan")
        K_rms_target_GeV = float(math.sqrt(K2_target_GeV2)) if (math.isfinite(K2_target_GeV2) and K2_target_GeV2 > 0) else float("nan")

        # Cascade magnitude path B target
        phi_star_target = float((rho_L_GeV4 ** 0.25) / Mpl) if (math.isfinite(rho_L_GeV4) and rho_L_GeV4 > 0) else float("nan")

        # Best-effort: load the derived torsion-condensate prediction from the current output_dir if available.
        tc_out = Path(config.output_dir) / "torsion_condensate" / "results.json"
        tc_best = None
        tc_log10_mismatch = None
        if tc_out.is_file():
            try:
                payload = json.loads(tc_out.read_text(encoding="utf-8"))
                best = payload.get("results", {}).get("gap_equation", {}).get("best", None) if isinstance(payload, dict) else None
                if isinstance(best, dict):
                    tc_best = best
                    if math.isfinite(rho_L_GeV4) and rho_L_GeV4 > 0:
                        rho_pred = float(best.get("rho_L_GeV4", float("nan")))
                        if math.isfinite(rho_pred) and rho_pred > 0:
                            tc_log10_mismatch = abs(math.log10(rho_pred / rho_L_GeV4))
            except Exception:
                tc_best = None

        checks: list[Check] = []
        ok = all(math.isfinite(x) for x in [Omega_L, H0_GeV, Lambda_obs_GeV2, rho_L_GeV4])
        checks.append((mk_check_pass if ok else mk_check_warn)("computed_lambda_targets", f"Omega_L={Omega_L}, H0_GeV={H0_GeV:.3e}"))
        checks.append(mk_check_info("lambda_obs", f"Lambda_obs≈{Lambda_obs_GeV2:.3e} GeV^2 (from 3 Ω_L H0^2)"))
        checks.append(mk_check_info("rho_lambda_obs", f"rho_L≈{rho_L_GeV4:.3e} GeV^4 (from Ω_L ρ_c)"))
        checks.append(mk_check_info("ufe_K_rms_target", f"<K^2>_target=4Λ => K_rms≈{K_rms_target_GeV:.3e} GeV"))
        if tc_best is not None:
            checks.append(
                mk_check_pass(
                    "torsion_condensate_prediction_available",
                    f"found torsion_condensate best n={tc_best.get('n')} (log10 mismatch vs rho_L target≈{tc_log10_mismatch})",
                )
            )
        else:
            checks.append(mk_check_info("torsion_condensate_prediction_available", f"no torsion_condensate output found at {tc_out} (targets only)"))
        checks.append(mk_check_info("derivation_interface", "Λ targets are computed here; the shipped prediction candidate lives in torsion_condensate (gap equation; no continuous tuning)."))

        # Ladder terminal-stage identification (continuous extrapolation).
        c = TfptConstants.compute()
        ladder_seq: list[dict[str, float]] = []
        for n in range(1, LADDER_N_MAX + 1):
            try:
                phi_n = float(e8_phi_n(n=n, c=c))
                ladder_seq.append({"n": n, "D_n": float(LADDER_D0 - 2 * n), "phi_n": phi_n})
            except Exception:
                continue

        phi_star_mp = mp.mpf(phi_star_target) if (math.isfinite(phi_star_target) and phi_star_target > 0) else mp.mpf("nan")
        varphi_pref = mp.mpf(c.varphi0) * mp.e ** (-mp.mpf(c.gamma0))
        if mp.isfinite(phi_star_mp) and phi_star_mp > 0 and mp.isfinite(varphi_pref) and varphi_pref > 0:
            ratio = phi_star_mp / varphi_pref
            D_terminal = LADDER_D1 * (ratio ** (mp.mpf(1) / mp.mpf(c.e8_lambda))) if ratio > 0 else mp.mpf("nan")
            n_terminal = (LADDER_D0 - D_terminal) / mp.mpf(2) if mp.isfinite(D_terminal) else mp.mpf("nan")
            phi_terminal = varphi_pref * (D_terminal / LADDER_D1) ** mp.mpf(c.e8_lambda) if mp.isfinite(D_terminal) else mp.mpf("nan")
            if mp.isfinite(phi_terminal) and mp.isfinite(phi_star_mp) and phi_terminal > 0 and phi_star_mp > 0:
                log10_mismatch_phi = abs(mp.log10(phi_terminal / phi_star_mp))
            else:
                log10_mismatch_phi = mp.mpf("nan")
        else:
            D_terminal = mp.mpf("nan")
            n_terminal = mp.mpf("nan")
            phi_terminal = mp.mpf("nan")
            log10_mismatch_phi = mp.mpf("nan")

        terminal_ok = bool(
            mp.isfinite(n_terminal)
            and abs(n_terminal - TERMINAL_STAGE_EXPECTED) <= TERMINAL_STAGE_TOL
            and mp.isfinite(log10_mismatch_phi)
            and log10_mismatch_phi <= mp.mpf("1e-6")
        )
        checks.append(
            mk_check_pass("ladder_terminal_stage_identified", f"n_terminal≈{n_terminal} (D≈{D_terminal})")
            if terminal_ok
            else mk_check_warn("ladder_terminal_stage_identified", f"n_terminal≈{n_terminal} (D≈{D_terminal}, mismatch≈{log10_mismatch_phi})")
        )

        report_lines: list[str] = [
            "Dark energy (Λ) targets — UFE torsion condensate vs cascade terminal stage",
            f"mode = {mode}",
            f"cosmology source: {kcal_path}",
            "",
            "Flat ΛCDM snapshot:",
            f"- H0 = {H0_km_s_Mpc} km/s/Mpc => H0 = {H0_GeV:.3e} GeV",
            f"- Ω_m = {Omega_m}, Ω_r = {Omega_r} => Ω_Λ = {Omega_L}",
            "",
            "Observed targets (derived from the snapshot):",
            f"- Λ_obs = 3 Ω_Λ H0^2 = {Lambda_obs_GeV2:.3e} GeV^2",
            f"- ρ_c = 3 H0^2 M̄_P^2 = {rho_c_GeV4:.3e} GeV^4 (M̄_P={MPL_REDUCED_GEV:.3e} GeV)",
            f"- ρ_Λ = Ω_Λ ρ_c = {rho_L_GeV4:.3e} GeV^4",
            "",
            "TFPT path A (five_problems: torsion condensate):",
            "- Λ_eff = (1/4)⟨K^2⟩  ⇒  ⟨K^2⟩_target = 4 Λ_obs",
            f"- ⟨K^2⟩_target = {K2_target_GeV2:.3e} GeV^2  ⇒  K_rms = {K_rms_target_GeV:.3e} GeV",
            "",
            "TFPT path B (five_problems / v1.06: terminal cascade magnitude):",
            "- ρ_Λ ≈ (M̄_P φ_*)^4  ⇒  φ_* = (ρ_Λ)^{1/4} / M̄_P",
            f"- φ_* target ≈ {phi_star_target:.3e}",
            "",
            "E8 ladder terminal stage (extrapolated):",
            f"- varphi_n sequence (n=1..{LADDER_N_MAX}) available ({len(ladder_seq)} points)",
            f"- n_terminal ≈ {n_terminal} (D≈{D_terminal}, mismatch≈{log10_mismatch_phi})",
            "",
            "torsion_condensate (if available in output_dir):",
            f"- results.json: {tc_out} (present={tc_out.is_file()})",
            f"- best: {tc_best}",
            f"- log10 mismatch vs rho_L target: {tc_log10_mismatch}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
            "",
            "Notes:",
            "- This module is intentionally an engineering target: it translates cosmology into required TFPT mechanism magnitudes.",
            "- To upgrade to a prediction, implement either (A) a computable ⟨K^2⟩ from the UFE/microscopic torsion sector, or (B) a deterministic E8 ladder/block model producing φ_n at the claimed terminal stage.",
            "",
        ]

        results: dict[str, Any] = {
            "mode": mode,
            "cosmology_source": str(kcal_path),
            "cosmology": {"H0_km_s_Mpc": H0_km_s_Mpc, "Omega_m": Omega_m, "Omega_r": Omega_r, "Omega_L": Omega_L},
            "targets": {
                "H0_GeV": H0_GeV,
                "Lambda_obs_GeV2": Lambda_obs_GeV2,
                "rho_c_GeV4": rho_c_GeV4,
                "rho_L_GeV4": rho_L_GeV4,
                "ufe_path_A": {"K2_target_GeV2": K2_target_GeV2, "K_rms_target_GeV": K_rms_target_GeV},
                "cascade_path_B": {"phi_star_target": phi_star_target},
            },
            "ladder_terminal_stage": {
                "sequence": ladder_seq,
                "n_terminal": float(n_terminal) if mp.isfinite(n_terminal) else None,
                "D_terminal": float(D_terminal) if mp.isfinite(D_terminal) else None,
                "phi_terminal": float(phi_terminal) if mp.isfinite(phi_terminal) else None,
                "log10_mismatch_phi_star": float(log10_mismatch_phi) if mp.isfinite(log10_mismatch_phi) else None,
            },
            "torsion_condensate_if_available": {
                "results_file": str(tc_out),
                "present": bool(tc_out.is_file()),
                "best": tc_best,
                "log10_mismatch_rho_L": tc_log10_mismatch,
            },
        }

        return ModuleResult(results=results, checks=checks, report="\n".join(report_lines), warnings=[])

