from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.cosmo_scale_map import MPL_REDUCED_GEV
from tfpt_suite.defect_partition import alpha_inv_0_from_delta2, derive_delta2_from_defect_partition
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from theoryv3_suite.utils import LOG10_DEX_TOL_DEFAULT, ensure_ascii, safe_log10


TOL_DEX = mp.mpf(LOG10_DEX_TOL_DEFAULT)
NORM_CANDIDATES = [
    ("n=1", mp.mpf(1)),
    ("n=1/sqrt(2)", mp.mpf(1) / mp.sqrt(2)),
    ("n=1/sqrt(3)", mp.mpf(1) / mp.sqrt(3)),
    ("n=1/2", mp.mpf(1) / 2),
    ("n=1/3^(1/4)", mp.mpf(1) / (mp.mpf(3) ** (mp.mpf(1) / 4))),
]


def _H0_GeV_from_km_s_Mpc(H0_km_s_Mpc: float) -> float:
    Mpc_km = 3.0856775814913673e19
    hbar_GeV_s = 6.582119569e-25
    H0_s_inv = float(H0_km_s_Mpc) / Mpc_km
    return float(H0_s_inv * hbar_GeV_s)


def _plot_rho_candidates(*, out_dir: Path, candidates: list[dict[str, Any]], rho_target: float) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"rho_lambda_candidates_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels = [c["label"] for c in candidates]
        values = [safe_log10(float(c["rho_L_GeV4"])) for c in candidates]
        target_log = safe_log10(float(rho_target))

        fig, ax = plt.subplots(figsize=(8.5, 4.0))
        ax.bar(labels, values, color="#805ad5")
        ax.axhline(float(target_log), color="black", lw=1.0, ls="--", alpha=0.8, label="target log10 rho_L")
        ax.set_ylabel("log10 rho_L [GeV^4]")
        ax.set_title("Dark energy candidates (exp(-alpha_inv/2))")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.legend(loc="best")

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "rho_lambda_candidates.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["rho_lambda_candidates_png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class DarkEnergyExponentialAuditModule(TfptModule):
    module_id = "dark_energy_exponential_audit"
    title = "Dark energy exponential audit (phi_star from exp(-alpha_inv/2))"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "alpha_inv_0 from g=5 delta2 model",
                "k_calibration cosmology snapshot (Omega_m, Omega_r, H0)",
            ],
            outputs=[
                "phi_star_base = exp(-alpha_inv_0/2)",
                "rho_L for discrete normalization candidates",
                "best candidate and log10 mismatch vs target",
            ],
            formulas=[
                "phi_star_base = exp(-alpha_inv_0 / 2)",
                "rho_L = (Mpl_bar * phi_star)^4",
                "Omega_L = 1 - Omega_m - Omega_r",
                "rho_L_target = 3 * Omega_L * H0^2 * Mpl_bar^2",
            ],
            validation=[
                "best candidate within 0.5 dex of target",
                "n=1/2 branch is preferred under current scan",
            ],
            determinism="Deterministic (finite scan, no continuous tuning).",
            question="Is the dark energy scale explained by exp(-alpha_inv/2) with a discrete normalization?",
            objective=[
                "Quantify the exponential suppression mechanism with explicit candidates.",
                "Track mismatch against the cosmology target ledger.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        d2 = derive_delta2_from_defect_partition(delta_top=mp.mpf(c.delta_top))
        alpha_inv_0 = alpha_inv_0_from_delta2(delta2=d2.delta2, mp_dps=int(getattr(config, "mp_dps", 80)))
        phi_star_base = mp.e ** (-(mp.mpf(1) / 2) * alpha_inv_0)

        data_dir = Path(__file__).resolve().parents[3] / "tfpt_suite" / "data"
        kcal_path = data_dir / "k_calibration.json"
        kcal = json.loads(kcal_path.read_text(encoding="utf-8"))
        cosmo = kcal.get("cosmology_flat_lcdm", {})
        H0_km_s_Mpc = float(cosmo.get("H0_km_s_Mpc", float("nan")))
        Omega_m = float(cosmo.get("Omega_m", float("nan")))
        Omega_r = float(cosmo.get("Omega_r", 0.0))
        Omega_L = float(1.0 - Omega_m - Omega_r) if (math.isfinite(Omega_m) and math.isfinite(Omega_r)) else float("nan")

        H0_GeV = _H0_GeV_from_km_s_Mpc(H0_km_s_Mpc) if math.isfinite(H0_km_s_Mpc) else float("nan")
        rho_L_target = float(3.0 * Omega_L * (H0_GeV**2) * (float(MPL_REDUCED_GEV) ** 2)) if math.isfinite(Omega_L) else float("nan")

        candidates: list[dict[str, Any]] = []
        best: dict[str, Any] | None = None
        for label, nrm in NORM_CANDIDATES:
            phi_star = phi_star_base * nrm
            rho_pred = (mp.mpf(MPL_REDUCED_GEV) * phi_star) ** 4
            log10_mismatch = abs(mp.log10(rho_pred / mp.mpf(rho_L_target))) if rho_L_target > 0 else mp.mpf("nan")
            row = {
                "label": str(label),
                "norm_factor": str(nrm),
                "phi_star": str(phi_star),
                "rho_L_GeV4": str(rho_pred),
                "log10_mismatch_rho_L": str(log10_mismatch),
            }
            candidates.append(row)
            if mp.isfinite(log10_mismatch):
                if best is None or mp.mpf(best["log10_mismatch_rho_L"]) > log10_mismatch:
                    best = dict(row)

        best_mis = mp.mpf(best["log10_mismatch_rho_L"]) if best is not None else mp.mpf("nan")
        ok = bool(best is not None and mp.isfinite(best_mis) and best_mis <= TOL_DEX)
        best_label = str(best["label"]) if isinstance(best, dict) and "label" in best else "n/a"
        if isinstance(best, dict):
            best["label_str"] = best_label

        checks: list[Check] = []
        checks.append(
            mk_check_pass("rho_L_within_tol", f"best={best_label}, log10_mismatch={best_mis}")
            if ok
            else mk_check_warn("rho_L_within_tol", f"best={best_label}, log10_mismatch={best_mis}")
        )
        checks.append(
            mk_check_pass("n_half_preferred", f"best={best_label}")
            if best_label == "n=1/2"
            else mk_check_warn("n_half_preferred", f"best={best_label} (expected n=1/2)")
        )

        lines = [
            "Dark energy exponential audit",
            "",
            f"alpha_inv_0 (g/4 delta2) = {alpha_inv_0}",
            f"phi_star_base = exp(-alpha_inv_0/2) = {phi_star_base}",
            "",
            "Target ledger (from k_calibration.json):",
            f"- H0 = {H0_km_s_Mpc} km/s/Mpc",
            f"- Omega_L = {Omega_L}",
            f"- rho_L_target = {rho_L_target} GeV^4",
            "",
            "Candidates:",
            *[f"- {c['label']}: log10 mismatch={c['log10_mismatch_rho_L']} (phi_star={c['phi_star']})" for c in candidates],
            "",
            f"Best: {best}",
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"rho_lambda_candidates_png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_rho_candidates(out_dir=out_dir, candidates=candidates, rho_target=rho_L_target)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "alpha_inv_0": str(alpha_inv_0),
                "phi_star_base": str(phi_star_base),
                "targets": {"H0_km_s_Mpc": H0_km_s_Mpc, "Omega_L": Omega_L, "rho_L_target": rho_L_target},
                "candidates": candidates,
                "best": best,
                "best_label_str": best_label,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
