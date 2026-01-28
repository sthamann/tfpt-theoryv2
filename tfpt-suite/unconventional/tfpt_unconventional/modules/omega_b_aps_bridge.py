from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.modules.aps_eta_gluing import _spectral_flow_dirac_1d


def _plot_omega_b_aps_bridge(
    *,
    out_dir: Path,
    K_candidate: mp.mpf,
    implied_K: mp.mpf,
    sigma_implied_K: mp.mpf,
    omega_b_pred: mp.mpf,
    omega_b_ref: mp.mpf,
    sigma_omega_b_ref: mp.mpf,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"omega_b_aps_bridge_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(11.2, 4.2))

        # Panel 1: coefficient K
        ax1.errorbar(
            [0],
            [float(implied_K)],
            yerr=[float(sigma_implied_K)],
            fmt="o",
            color="#1f77b4",
            label="implied K (Planck-derived)",
            capsize=4,
        )
        ax1.scatter([0], [float(K_candidate)], s=140, marker="*", color="#d62728", edgecolor="black", linewidth=0.7, label="K_candidate")
        ax1.set_xticks([0])
        ax1.set_xticklabels(["K"])
        ax1.set_ylabel("value")
        ax1.set_title("Ω_b coefficient K := Ω_b/β_rad")
        ax1.grid(True, ls=":", alpha=0.35)
        ax1.legend(loc="best")

        # Panel 2: Ω_b
        ax2.errorbar(
            [0],
            [float(omega_b_ref)],
            yerr=[float(sigma_omega_b_ref)],
            fmt="o",
            color="#1f77b4",
            label="Ω_b ref (Planck-derived)",
            capsize=4,
        )
        ax2.scatter([0], [float(omega_b_pred)], s=140, marker="*", color="#d62728", edgecolor="black", linewidth=0.7, label="Ω_b pred")
        ax2.set_xticks([0])
        ax2.set_xticklabels([r"$\Omega_b$"])
        ax2.set_ylabel("value")
        ax2.set_title("Ω_b prediction vs reference (diagnostic)")
        ax2.grid(True, ls=":", alpha=0.35)
        ax2.legend(loc="best")

        fig.tight_layout()
        path = out_dir / "omega_b_aps_bridge.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["omega_b_aps_bridge_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class _Assumption:
    id: str
    statement: str
    status: str


class OmegaBApsBridgeModule(TfptModule):
    module_id = "ux_omega_b_aps_bridge"
    title = "Unconventional: Ω_b coefficient bridge via APS seam term (ΔΓ)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT beta_rad = varphi0/(4π) from constants.py",
                "APS seam toy-model: spectral flow SF(U_Γ) for U_Γ(θ)=exp(i m θ) (module `aps_eta_gluing`)",
                "Planck 2018 reference (Ω_b h^2 and H0) for a diagnostic implied coefficient K := Ω_b/beta_rad",
            ],
            outputs=[
                "ΔΓ from APS spectral flow (m=1)",
                "candidate coefficient K_candidate := 2·ΔΓ − 1",
                "Ω_b prediction using K_candidate (diagnostic vs Planck-derived Ω_b)",
            ],
            formulas=[
                "APS seam term (toy-model): ΔΓ = 2π · SF(U_Γ), with minimal nontrivial class m=1 ⇒ SF=1 ⇒ ΔΓ=2π",
                "bridge (conditional): K_candidate = 2·ΔΓ − 1  (⇒ 4π − 1)",
                "Ω_b_pred = K_candidate · beta_rad",
            ],
            validation=[
                "show numerically that K_candidate equals (4π−1) once ΔΓ=2π is fixed by the spectral-flow computation",
                "report agreement vs Planck-derived Ω_b as diagnostic only",
            ],
            determinism="Deterministic; no randomness.",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        beta_rad = mp.mpf(c.beta_rad)

        # --- APS seam term (reuse the same toy-model as `aps_eta_gluing`) ---
        epsilon = 1e-3
        row = _spectral_flow_dirac_1d(spin="periodic", m=1, epsilon=epsilon)
        sf = mp.mpf(int(row.numeric_sf))
        delta_gamma = mp.mpf(2) * mp.pi * sf  # ΔΓ := 2π·SF

        # --- Bridge conjecture (conditional; to be replaced by an operator/anomaly-level derivation) ---
        K_candidate = mp.mpf(2) * delta_gamma - mp.mpf(1)
        omega_b_pred = K_candidate * beta_rad

        # --- Reference Ω_b from Planck (same reference file as the main Ω_b module) ---
        # `.../tfpt-suite/unconventional/tfpt_unconventional/modules/<this_file>`
        # -> parents[3] is `tfpt-suite/`
        ref_path = Path(__file__).resolve().parents[3] / "tfpt_suite" / "data" / "global_reference_minimal.json"
        ref = json.loads(ref_path.read_text(encoding="utf-8"))
        obs = ref.get("observables", {}) if isinstance(ref, dict) else {}

        omega_b_h2 = mp.mpf(str(obs["omega_b_h2_planck2018"]["mean"]))
        sigma_omega_b_h2 = mp.mpf(str(obs["omega_b_h2_planck2018"]["sigma"]))
        H0 = mp.mpf(str(obs["H0_planck2018_km_s_Mpc"]["mean"]))
        sigma_H0 = mp.mpf(str(obs["H0_planck2018_km_s_Mpc"]["sigma"]))

        h = H0 / mp.mpf(100)
        sigma_h = sigma_H0 / mp.mpf(100)
        omega_b_ref = omega_b_h2 / (h**2)
        sigma_omega_b_ref = mp.sqrt((sigma_omega_b_h2 / (h**2)) ** 2 + ((omega_b_h2 * mp.mpf(2) * sigma_h) / (h**3)) ** 2)

        implied_K = omega_b_ref / beta_rad
        sigma_implied_K = sigma_omega_b_ref / beta_rad if beta_rad != 0 else mp.mpf("inf")
        z_ref = abs(omega_b_pred - omega_b_ref) / sigma_omega_b_ref if sigma_omega_b_ref != 0 else mp.mpf("inf")

        assumptions = [
            _Assumption(
                id="aps_seam_term_toy_model",
                statement="Use the 1D Dirac-family toy model of the seam operator (Appendix seam) where ΔΓ := 2π·SF(U_Γ) for U_Γ(θ)=exp(i m θ).",
                status="implemented (see aps_eta_gluing)",
            ),
            _Assumption(
                id="minimal_nontrivial_class_m_eq_1",
                statement="Use the minimal nontrivial winding class m=1 (so SF=1 in the toy model).",
                status="checked in aps_eta_gluing (minimal m with SF>0 is 1)",
            ),
            _Assumption(
                id="omega_b_coeff_equals_2_delta_gamma_minus_1",
                statement="Bridge postulate: the baryon coefficient K := Ω_b/β_rad equals (2·ΔΓ − 1) in the same normalization.",
                status="conditional / to be replaced by an operator/anomaly-level derivation",
            ),
        ]

        checks = [
            Check(
                check_id="aps_spectral_flow_m1",
                passed=bool(int(row.numeric_sf) == 1 and int(row.winding_det_u) == 1),
                detail=f"m=1 periodic: numeric_sf={row.numeric_sf}, winding={row.winding_det_u}",
            ),
            Check(
                check_id="delta_gamma_equals_2pi",
                passed=bool(abs(delta_gamma - (mp.mpf(2) * mp.pi)) < mp.mpf("1e-40")),
                detail=f"ΔΓ=2π·SF with SF=1 ⇒ ΔΓ={delta_gamma}",
            ),
            Check(
                check_id="K_candidate_equals_4pi_minus_1",
                passed=bool(abs(K_candidate - (mp.mpf(4) * mp.pi - mp.mpf(1))) < mp.mpf("1e-40")),
                detail=f"K_candidate=2ΔΓ−1={K_candidate} (vs 4π−1)",
            ),
            Check(
                check_id="omega_b_pred_positive",
                passed=bool(omega_b_pred > 0),
                detail=f"Ω_b_pred={omega_b_pred}",
            ),
            Check(
                check_id="omega_b_pred_close_to_planck_ref_diagnostic",
                passed=bool(z_ref < mp.mpf(3)),
                detail=f"Ω_b_ref={omega_b_ref} ± {sigma_omega_b_ref}; z={z_ref}",
            ),
        ]

        lines: list[str] = []
        lines += [
            "Unconventional: Ω_b bridge via APS seam term",
            "",
            "Idea:",
            "- Replace a raw (4π−1) insertion by a quantity already computed in-suite: the APS seam term ΔΓ.",
            "- This remains conditional until an operator/anomaly-level derivation connects Ω_b to ΔΓ.",
            "",
            "APS seam toy-model (periodic spin):",
            f"- epsilon = {epsilon}",
            f"- m = 1 ⇒ SF(U_Γ) = {row.numeric_sf} (analytic={row.analytic_sf}), winding(det U_Γ) = {row.winding_det_u}",
            f"- ΔΓ := 2π·SF = {delta_gamma}",
            "",
            "Bridge conjecture (conditional):",
            f"- K_candidate := 2·ΔΓ − 1 = {K_candidate}",
            f"- beta_rad = varphi0/(4π) = {beta_rad}",
            f"- Ω_b_pred := K_candidate·beta_rad = {omega_b_pred}",
            "",
            "Reference (Planck 2018 base-ΛCDM; derived Ω_b):",
            f"- Ω_b h^2 = {omega_b_h2} ± {sigma_omega_b_h2}",
            f"- H0 = {H0} ± {sigma_H0} km/s/Mpc ⇒ h={h} ± {sigma_h}",
            f"- Ω_b_ref = Ω_b h^2 / h^2 = {omega_b_ref} ± {sigma_omega_b_ref}",
            f"- implied K = Ω_b_ref / beta_rad = {implied_K} ± {sigma_implied_K}",
            f"- |K_candidate - implied K| = {abs(K_candidate - implied_K)} (diagnostic)",
            "",
            "Assumptions (explicit):",
            *[f"- {a.id}: {a.statement}  [{a.status}]" for a in assumptions],
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module is intentionally conservative: it does not claim an operator-level Ω_b derivation.",
            "- It provides a concrete *docking point*: Ω_b coefficient ↔ ΔΓ (APS seam), which can be upgraded to an anomaly/inflow computation later.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"omega_b_aps_bridge_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_omega_b_aps_bridge(
                out_dir=out_dir,
                K_candidate=K_candidate,
                implied_K=implied_K,
                sigma_implied_K=sigma_implied_K,
                omega_b_pred=omega_b_pred,
                omega_b_ref=omega_b_ref,
                sigma_omega_b_ref=sigma_omega_b_ref,
            )
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "tfpt": {"beta_rad": beta_rad},
                "aps_seam": {
                    "epsilon": float(epsilon),
                    "spin": "periodic",
                    "m": 1,
                    "spectral_flow": int(row.numeric_sf),
                    "winding": int(row.winding_det_u),
                    "delta_gamma": delta_gamma,
                },
                "bridge": {
                    "K_candidate": K_candidate,
                    "omega_b_pred": omega_b_pred,
                    "assumptions": [a.__dict__ for a in assumptions],
                    "conditional": True,
                },
                "reference": {
                    "file": str(ref_path),
                    "omega_b_ref": omega_b_ref,
                    "sigma_omega_b_ref": sigma_omega_b_ref,
                    "implied_K": implied_K,
                    "sigma_implied_K": sigma_implied_K,
                    "z_score": z_ref,
                },
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

