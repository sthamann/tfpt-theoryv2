from __future__ import annotations

import json
import math
from pathlib import Path

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import (
    Check,
    ModuleResult,
    ModuleSpec,
    TfptModule,
    mk_check_info,
    mk_check_pass,
    mk_check_warn,
)
from tfpt_suite.reference_ledger import get_dataset
from theoryv3_suite.utils import LOG10_DEX_TOL_DEFAULT, Z_TOL_DEFAULT, ensure_ascii, safe_log10


ETA10_COEFF = 273.9
M_STAR_EV = 1.08e-3


def _severity_from_z(z: float) -> str:
    if not math.isfinite(z):
        return "WARN"
    az = abs(z)
    if az <= float(Z_TOL_DEFAULT):
        return "PASS"
    if az <= 3.0:
        return "WARN"
    return "FAIL"


def _kappa_eff_strong_washout(K: float) -> float:
    if not math.isfinite(K) or K <= 0:
        return float("nan")
    if K < 3.0:
        return float(1.0 / max(2.0, K))
    return float(0.3 / (K * (math.log(K) ** 0.6)))


def _plot_baryon_summary(
    *,
    out_dir: Path,
    omega_b_pred: float,
    omega_b_ref: float,
    H0_pred: float,
    H0_ref: float,
    eta_b_pred: float,
    eta_b_ref: float,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"baryon_consistency_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8.5, 5.6))

        axes[0].bar(["Omega_b pred", "Omega_b ref"], [omega_b_pred, omega_b_ref], color=["#2f855a", "#c05621"])
        axes[0].set_ylabel("Omega_b")
        axes[0].set_title("Omega_b identity and derived H0")
        axes[0].grid(True, axis="y", ls=":", alpha=0.4)

        axes[1].bar(
            ["H0 pred", "H0 ref", "eta_b pred", "eta_b ref"],
            [H0_pred, H0_ref, eta_b_pred * 1e10, eta_b_ref * 1e10],
            color=["#2f855a", "#c05621", "#2f855a", "#c05621"],
        )
        axes[1].set_ylabel("H0 (km/s/Mpc) and eta_b (x1e-10)")
        axes[1].grid(True, axis="y", ls=":", alpha=0.4)

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "baryon_consistency.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["baryon_consistency_png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class BaryonConsistencyAuditModule(TfptModule):
    module_id = "baryon_consistency_audit"
    title = "Baryon consistency audit (Omega_b, eta_b, derived H0)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants (beta_rad)",
                "Planck omega_b_h2 and H0 references",
                "Alternative H0 reference (SH0ES)",
                "seesaw thresholds and neutrino masses (eta_b proxy)",
            ],
            outputs=[
                "Omega_b = (4*pi - 1) * beta_rad",
                "eta_b proxy (leptogenesis-style)",
                "eta10_pred = 1e10 * eta_b_pred",
                "omega_b_h2_pred = eta10_pred / 273.9",
                "derived H0 from Omega_b and eta_b",
            ],
            formulas=[
                "Omega_b = (4*pi - 1) * beta_rad",
                "eta10 = 273.9 * omega_b_h2",
                "eta_b = 0.96e-2 * eps1 * kappa_eff",
                "eta10_pred = 1e10 * eta_b_pred",
                "omega_b_h2_pred = eta10_pred / 273.9",
                "H0 = 100 * sqrt((eta10_pred/273.9) / Omega_b_pred)",
            ],
            validation=[
                "Omega_b within 2 sigma of reference",
                "eta_b within 0.5 dex of Planck anchor",
                "H0 within 2 sigma of Planck anchor (derived check)",
            ],
            determinism="Deterministic given inputs.",
            question="Do baryon observables self-consistently close using beta_rad?",
            objective=[
                "Expose the Omega_b identity and baryogenesis proxy in one audit.",
                "Quantify the implied H0 from internal TFPT quantities.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        beta_rad = float(c.beta_rad)
        omega_b_pred = float((4 * mp.pi - 1) * mp.mpf(beta_rad))

        data_dir = Path(__file__).resolve().parents[3] / "tfpt_suite" / "data"
        thr_path = data_dir / "rge_thresholds_v25.json"
        tex_path = data_dir / "flavor_texture_v24.json"

        thr = json.loads(thr_path.read_text(encoding="utf-8"))
        tex = json.loads(tex_path.read_text(encoding="utf-8"))

        omega_ref = get_dataset("omega_b_h2_planck2018")
        H0_planck_ref = get_dataset("H0_planck2018")
        H0_alt_ref = get_dataset("H0_sh0es_2022")

        omega_b_h2 = float(omega_ref["value"])
        omega_b_h2_sigma = float(omega_ref["sigma"])
        H0_ref = float(H0_planck_ref["value"])
        H0_sigma = float(H0_planck_ref["sigma"])
        h_ref = H0_ref / 100.0
        omega_b_ref = omega_b_h2 / (h_ref**2)

        # Omega_b sigma propagation
        sigma_h = H0_sigma / 100.0
        dOmega_domega = 1.0 / (h_ref**2)
        dOmega_dh = -2.0 * omega_b_h2 / (h_ref**3)
        omega_b_sigma = math.sqrt((dOmega_domega * omega_b_h2_sigma) ** 2 + (dOmega_dh * sigma_h) ** 2)
        z_omega = (omega_b_pred - omega_b_ref) / omega_b_sigma if omega_b_sigma > 0 else float("nan")

        # Baryogenesis proxy (same as baryogenesis_mechanism)
        thresholds = thr.get("thresholds_GeV", {})
        M1_GeV = float(thresholds.get("MNR1", float("nan")))
        nm = (tex.get("neutrino_mechanism", {}) or {}).get("neutrino_masses_input_eV", [0.0, 0.0086, 0.05])
        m_eV = [float(x) for x in nm] if isinstance(nm, list) and len(nm) == 3 else [0.0, 0.0086, 0.05]
        m1_eV, m2_eV, m3_eV = m_eV

        v_GeV = 246.0
        delta_m_eV = max(0.0, m3_eV - m1_eV)
        delta_m_GeV = delta_m_eV * 1e-9
        eps1_max = (3.0 / (16.0 * math.pi)) * (M1_GeV * delta_m_GeV / (v_GeV**2)) if (math.isfinite(M1_GeV) and M1_GeV > 0) else float("nan")
        eps1 = beta_rad * eps1_max if math.isfinite(eps1_max) else float("nan")
        K = m3_eV / M_STAR_EV if M_STAR_EV > 0 else float("nan")
        kappa = _kappa_eff_strong_washout(K)
        eta_b_pred = 0.96e-2 * eps1 * kappa if (math.isfinite(eps1) and math.isfinite(kappa)) else float("nan")

        # Reference eta_b from omega_b_h2
        eta10_obs = ETA10_COEFF * omega_b_h2
        eta_b_ref = eta10_obs * 1e-10
        log10_mismatch = safe_log10(float(eta_b_pred / eta_b_ref)) if (eta_b_pred > 0 and eta_b_ref > 0) else float("nan")

        # Derived H0 from Omega_b + eta_b
        eta10_pred = eta_b_pred / 1e-10 if eta_b_pred > 0 else float("nan")
        omega_b_h2_pred = eta10_pred / ETA10_COEFF if eta10_pred > 0 else float("nan")
        h_pred = math.sqrt(omega_b_h2_pred / omega_b_pred) if (omega_b_pred > 0 and omega_b_h2_pred > 0) else float("nan")
        H0_pred = 100.0 * h_pred if math.isfinite(h_pred) else float("nan")
        z_H0_planck = (H0_pred - H0_ref) / H0_sigma if math.isfinite(H0_pred) and H0_sigma > 0 else float("nan")
        H0_alt_value = float(H0_alt_ref["value"])
        H0_alt_sigma = float(H0_alt_ref["sigma"])
        z_H0_alt = (H0_pred - H0_alt_value) / H0_alt_sigma if math.isfinite(H0_pred) and H0_alt_sigma > 0 else float("nan")

        checks: list[Check] = []
        checks.append(
            mk_check_pass("omega_b_within_2sigma", f"z={z_omega}")
            if math.isfinite(z_omega) and abs(z_omega) <= float(Z_TOL_DEFAULT)
            else mk_check_warn("omega_b_within_2sigma", f"z={z_omega}")
        )
        checks.append(
            mk_check_pass("eta_b_within_0p5_dex", f"log10 mismatch={log10_mismatch}")
            if math.isfinite(log10_mismatch) and abs(log10_mismatch) <= float(LOG10_DEX_TOL_DEFAULT)
            else mk_check_warn("eta_b_within_0p5_dex", f"log10 mismatch={log10_mismatch}")
        )
        if math.isfinite(z_H0_planck):
            sev = _severity_from_z(z_H0_planck)
            checks.append(
                Check(
                    check_id="H0_within_2sigma",
                    passed=sev != "FAIL",
                    detail=f"z_planck={z_H0_planck}",
                    severity=sev,
                )
            )
        if math.isfinite(z_H0_alt):
            checks.append(
                mk_check_info(
                    "H0_tension_indicator",
                    f"z_alt={z_H0_alt} (alt={H0_alt_value} ± {H0_alt_sigma}; planck={H0_ref} ± {H0_sigma})",
                )
            )

        lines = [
            "Baryon consistency audit",
            "",
            f"beta_rad = {beta_rad}",
            f"Omega_b_pred = {omega_b_pred}",
            f"Omega_b_ref = {omega_b_ref} (z={z_omega})",
            "",
            f"eta_b_pred = {eta_b_pred}",
            f"eta10_pred = {eta10_pred}",
            f"omega_b_h2_pred = {omega_b_h2_pred}",
            f"eta_b_ref = {eta_b_ref} (log10 mismatch={log10_mismatch})",
            "",
            f"H0_pred = {H0_pred} km/s/Mpc (z_planck={z_H0_planck}, z_alt={z_H0_alt})",
            f"H0_ref_planck = {H0_ref} ± {H0_sigma} km/s/Mpc ({H0_planck_ref.get('version')})",
            f"H0_ref_alt = {H0_alt_value} ± {H0_alt_sigma} km/s/Mpc ({H0_alt_ref.get('version')})",
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"baryon_consistency_png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_baryon_summary(
                out_dir=out_dir,
                omega_b_pred=omega_b_pred,
                omega_b_ref=omega_b_ref,
                H0_pred=H0_pred,
                H0_ref=H0_ref,
                eta_b_pred=eta_b_pred,
                eta_b_ref=eta_b_ref,
            )
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "references": {
                    "omega_b_h2_planck2018": omega_ref,
                    "H0_planck2018": H0_planck_ref,
                    "H0_sh0es_2022": H0_alt_ref,
                },
                "omega_b_pred": omega_b_pred,
                "omega_b_ref": omega_b_ref,
                "omega_b_sigma": omega_b_sigma,
                "omega_b_z": z_omega,
                "eta_b_pred": eta_b_pred,
                "eta10_pred": eta10_pred,
                "omega_b_h2_pred": omega_b_h2_pred,
                "h_pred": h_pred,
                "eta_b_ref": eta_b_ref,
                "eta_b_log10_mismatch": log10_mismatch,
                "H0_pred": H0_pred,
                "H0_ref_set": {
                    "planck": {"value": H0_ref, "sigma": H0_sigma},
                    "alt": {"value": H0_alt_value, "sigma": H0_alt_sigma},
                },
                "H0_z_planck": z_H0_planck,
                "H0_z_alt": z_H0_alt,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
