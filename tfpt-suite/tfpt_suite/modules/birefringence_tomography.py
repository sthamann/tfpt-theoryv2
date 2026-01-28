from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.optimize import minimize

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


@dataclass(frozen=True)
class FitResult:
    model: str
    params: dict[str, float]
    logL: float
    bic: float


def _is_cmb_like_point(*, dataset: str, z: float) -> bool:
    """
    Classify a birefringence point as "CMB-like" (single offset at z≈1100) vs "low-z tomography".

    This is a *data-family* heuristic:
    - use the dataset string when possible
    - fall back to a redshift threshold for robustness
    """
    ds = str(dataset).strip().lower()
    if any(tag in ds for tag in ("cmb", "planck", "wmap", "act")):
        return True
    return bool(float(z) >= 100.0)


def _weighted_mean_and_sigma(*, x: np.ndarray, sigma: np.ndarray) -> tuple[float, float]:
    w = 1.0 / np.square(sigma)
    wsum = float(np.sum(w))
    if not np.isfinite(wsum) or wsum <= 0:
        return float("nan"), float("nan")
    mu = float(np.sum(w * x) / wsum)
    sig = float(np.sqrt(1.0 / wsum))
    return mu, sig


def _beta_step(z: np.ndarray, *, beta_inf: float, z_c: float, w: float) -> np.ndarray:
    """
    Step model for *cumulative* birefringence to a source at redshift z:
    - z→0 (very nearby): beta → 0
    - z→∞ (very distant): beta → beta_inf

    Implemented as a smooth tanh transition in z.
    """
    s = 0.5 * (1.0 + np.tanh((z - z_c) / w))
    return beta_inf * s


def _beta_drift(z: np.ndarray, *, beta_inf: float, p: float) -> np.ndarray:
    """
    Smooth drift proxy for cumulative birefringence:
      beta(z) = beta_inf * (1 - (1+z)^(-p))
    so beta(0)=0 and beta(∞)=beta_inf.
    """
    return beta_inf * (1.0 - np.power(1.0 + z, -p))


def _loglike_marginalized_calibration(
    *,
    beta_obs: np.ndarray,
    beta_pred: np.ndarray,
    sigma: np.ndarray,
    beta_cal_mean: float,
    beta_cal_sigma: float,
) -> float:
    """
    Gaussian likelihood with a *global calibration offset* beta_cal integrated out:

      beta_obs = beta_pred + beta_cal + noise
      beta_cal ~ N(beta_cal_mean, beta_cal_sigma^2)
      noise_i ~ N(0, sigma_i^2) independent

    The marginal distribution is multivariate normal with covariance:
      diag(sigma_i^2) + beta_cal_sigma^2 * 1 1^T

    This is computed deterministically using Sherman–Morrison identities.
    """
    y = (beta_obs - beta_pred) - beta_cal_mean
    var = sigma * sigma
    inv_var = 1.0 / var

    a = float(np.sum(inv_var))  # 1^T S^{-1} 1
    b = float(np.sum(y * inv_var))  # 1^T S^{-1} y
    c = float(np.sum(y * y * inv_var))  # y^T S^{-1} y

    n = int(beta_obs.shape[0])
    logdet = float(np.sum(np.log(var)))
    quad = c

    if beta_cal_sigma > 0:
        denom = 1.0 + (beta_cal_sigma * beta_cal_sigma) * a
        quad = c - ((beta_cal_sigma * beta_cal_sigma) * (b * b)) / denom
        logdet = logdet + float(np.log(denom))

    return float(-0.5 * (quad + logdet + n * np.log(2.0 * np.pi)))


def _fit_model(
    *,
    z: np.ndarray,
    beta_obs: np.ndarray,
    sigma: np.ndarray,
    beta_cal_mean: float,
    beta_cal_sigma: float,
    beta_model: Callable[..., np.ndarray],
    x0: np.ndarray,
    bounds: list[tuple[float, float]],
    param_names: list[str],
    fixed_kwargs: dict[str, float],
    model_name: str,
) -> FitResult:
    def neg_loglike_of(x: np.ndarray) -> float:
        kwargs = dict(fixed_kwargs)
        kwargs.update({name: float(val) for name, val in zip(param_names, x)})
        pred = beta_model(z, **kwargs)
        ll = _loglike_marginalized_calibration(
            beta_obs=beta_obs, beta_pred=pred, sigma=sigma, beta_cal_mean=beta_cal_mean, beta_cal_sigma=beta_cal_sigma
        )
        return -ll

    res = minimize(
        neg_loglike_of,
        x0=x0,
        bounds=bounds,
        method="L-BFGS-B",
        options={"maxiter": 2000},
    )
    x_best = res.x
    logL = -neg_loglike_of(x_best)
    k = len(param_names)
    n = len(z)
    bic = k * float(np.log(n)) - 2.0 * logL
    return FitResult(
        model=model_name,
        params={name: float(val) for name, val in zip(param_names, x_best)},
        logL=logL,
        bic=bic,
    )


class BirefringenceTomographyModule(TfptModule):
    module_id = "birefringence_tomography"
    title = "Birefringence tomography: step vs drift (data ingestion + calibration nuisance)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT beta_inf = varphi0/(4π) (radians, n=1)",
                "dataset/config: tfpt_suite/data/birefringence_tomography.json (mode=data|synthetic)",
                "two phenomenological models: step vs drift",
            ],
            outputs=["best-fit parameters (step/drift)", "BIC comparison", "approx Bayes factor"],
            formulas=[
                "step: beta(z)=beta_inf * 0.5*(1 + tanh((z-z_c)/w))  (so beta(0)=0, beta(∞)=beta_inf)",
                "drift: beta(z)=beta_inf * (1 - (1+z)^(-p))  (so beta(0)=0, beta(∞)=beta_inf)",
                "likelihood: beta_obs = beta_model + beta_cal + noise, with beta_cal as calibration nuisance",
                "beta_cal prior: N(beta_cal_mean, beta_cal_sigma^2) and analytic marginalization",
                "BIC = k ln(n) - 2 log L_max",
                "approx Bayes factor (step over drift): exp(0.5*(BIC_drift - BIC_step))",
            ],
            validation=[
                "on synthetic step-generated data, step should be preferred by BIC/logL for reasonable noise",
            ],
            determinism="Deterministic given config seed (synthetic noise) and solver bounds.",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        beta_inf = float(c.beta_rad)  # radians

        cfg_path = Path(__file__).resolve().parent.parent / "data" / "birefringence_tomography.json"
        cfg = json.loads(cfg_path.read_text(encoding="utf-8")) if cfg_path.exists() else {"mode": "synthetic"}

        mode = str(cfg.get("mode", "synthetic")).strip().lower()
        used_mode = mode
        cfg_vetted = bool(cfg.get("vetted", False))

        # dataset containers (all in radians internally)
        z: np.ndarray
        beta_obs: np.ndarray
        sigma: np.ndarray
        beta_cal_mean: float
        beta_cal_sigma: float
        datasets: list[str] = []
        truth_model_used: str | None = None

        if mode == "data":
            data = cfg.get("data", {})
            units = str(data.get("units", "degrees")).strip().lower()
            points = list(data.get("points", []))
            if not points:
                used_mode = "synthetic"
            else:
                z = np.array([float(p["z"]) for p in points], dtype=float)
                datasets = [str(p.get("dataset", f"point_{i}")) for i, p in enumerate(points)]
                beta_vals = np.array([float(p["beta"]) if "beta" in p else float(p["beta_deg"]) for p in points], dtype=float)
                sigma_vals = np.array([float(p["sigma"]) if "sigma" in p else float(p["sigma_deg"]) for p in points], dtype=float)

                if units.startswith("deg"):
                    beta_obs = np.deg2rad(beta_vals)
                    sigma = np.deg2rad(sigma_vals)
                    beta_cal_mean = float(np.deg2rad(float(data.get("calibration_prior_mean_deg", 0.0))))
                    beta_cal_sigma = float(np.deg2rad(float(data.get("calibration_prior_sigma_deg", 0.0))))
                else:
                    beta_obs = beta_vals
                    sigma = sigma_vals
                    beta_cal_mean = float(data.get("calibration_prior_mean", 0.0))
                    beta_cal_sigma = float(data.get("calibration_prior_sigma", 0.0))
        if used_mode == "synthetic":
            synth = cfg.get("synthetic", {})
            z_max = float(synth.get("z_max", 5.0))
            n_pts = int(synth.get("n_pts", 60))
            z = np.linspace(0.0, z_max, n_pts)
            datasets = ["synthetic"] * int(z.shape[0])

            truth_model = str(synth.get("truth_model", "step")).strip().lower()
            truth_model_used = truth_model
            truth_params = dict(synth.get("truth_params", {}))
            if truth_model == "drift":
                true_p = float(truth_params.get("p", 2.0))
                beta_true = _beta_drift(z, beta_inf=beta_inf, p=true_p)
            else:
                true_zc = float(truth_params.get("z_c", 0.7))
                true_w = float(truth_params.get("w", 0.15))
                beta_true = _beta_step(z, beta_inf=beta_inf, z_c=true_zc, w=true_w)

            noise_frac = float(synth.get("noise_sigma_fraction_of_beta_inf", 0.06))
            sigma = np.full_like(z, noise_frac * beta_inf, dtype=float)

            beta_cal_mean = 0.0
            beta_cal_sigma = float(np.deg2rad(float(synth.get("calibration_prior_sigma_deg", 0.0))))
            # draw one calibration offset (deterministic via suite seed)
            beta_cal = float(np.random.normal(loc=beta_cal_mean, scale=beta_cal_sigma)) if beta_cal_sigma > 0 else 0.0

            beta_obs = beta_true + beta_cal + np.random.normal(loc=0.0, scale=sigma, size=n_pts)

        z_max = float(np.max(z)) if len(z) else 0.0

        # --- Data-family split: CMB-only is an offset check; "tomography" is low-z only. ---
        is_cmb = np.array([_is_cmb_like_point(dataset=ds, z=float(zi)) for ds, zi in zip(datasets, z)], dtype=bool)
        lowz_mask = ~is_cmb
        cmb_mask = is_cmb

        z_low = z[lowz_mask]
        beta_low = beta_obs[lowz_mask]
        sigma_low = sigma[lowz_mask]
        datasets_low = [datasets[i] for i in range(len(datasets)) if bool(lowz_mask[i])]

        z_cmb = z[cmb_mask]
        beta_cmb = beta_obs[cmb_mask]
        sigma_cmb = sigma[cmb_mask]
        datasets_cmb = [datasets[i] for i in range(len(datasets)) if bool(cmb_mask[i])]

        # Fit both models (beta_inf fixed by TFPT; only shape params are fitted).
        # IMPORTANT: we do *not* use CMB-only points for step-vs-drift discrimination (both saturate at z≈1100).
        fit_scope = "low_z_only" if (used_mode == "data" and int(z_low.shape[0]) >= 3) else "all_points"
        z_fit = z_low if fit_scope == "low_z_only" else z
        beta_fit = beta_low if fit_scope == "low_z_only" else beta_obs
        sigma_fit = sigma_low if fit_scope == "low_z_only" else sigma

        fit_step = _fit_model(
            z=z_fit,
            beta_obs=beta_fit,
            sigma=sigma_fit,
            beta_cal_mean=beta_cal_mean,
            beta_cal_sigma=beta_cal_sigma,
            beta_model=_beta_step,
            x0=np.array([1.0, 0.2]),
            bounds=[(0.0, z_max), (0.02, 2.0)],
            param_names=["z_c", "w"],
            fixed_kwargs={"beta_inf": beta_inf},
            model_name="step",
        )

        fit_drift = _fit_model(
            z=z_fit,
            beta_obs=beta_fit,
            sigma=sigma_fit,
            beta_cal_mean=beta_cal_mean,
            beta_cal_sigma=beta_cal_sigma,
            beta_model=_beta_drift,
            x0=np.array([1.0]),
            bounds=[(0.05, 10.0)],
            param_names=["p"],
            fixed_kwargs={"beta_inf": beta_inf},
            model_name="drift",
        )

        delta_bic = fit_drift.bic - fit_step.bic  # >0 favors step
        bayes_factor_step_over_drift = float(np.exp(0.5 * delta_bic))
        delta_logL = fit_step.logL - fit_drift.logL

        # Optional "TFPT prior" smooth-step hypothesis: fixed (z_c, w) from config (no fitted shape parameters).
        tfpt_prior = cfg.get("tfpt_prior", {}) if isinstance(cfg.get("tfpt_prior", {}), dict) else {}
        zc_tfpt = float(tfpt_prior.get("z_c", 0.7))
        w_tfpt = float(tfpt_prior.get("w", 0.15))
        beta_tfpt = _beta_step(z_fit, beta_inf=beta_inf, z_c=zc_tfpt, w=w_tfpt)
        logL_tfpt = _loglike_marginalized_calibration(
            beta_obs=beta_fit,
            beta_pred=beta_tfpt,
            sigma=sigma_fit,
            beta_cal_mean=beta_cal_mean,
            beta_cal_sigma=beta_cal_sigma,
        )
        bic_tfpt = float(-2.0 * logL_tfpt)  # k=0 (shape fixed by hypothesis)

        # CMB offset check (no tomography): compare each CMB constraint to TFPT beta_inf.
        # Note: these published constraints are not independent; we therefore report them individually
        # plus a naive weighted mean (for orientation only).
        cmb_rows: list[dict[str, float | str]] = []
        if int(z_cmb.shape[0]) >= 1:
            for ds, bb, ss in zip(datasets_cmb, beta_cmb, sigma_cmb):
                zscore = float((float(bb) - float(beta_inf)) / float(ss)) if float(ss) > 0 else float("nan")
                cmb_rows.append(
                    {
                        "dataset": str(ds),
                        "z": float(1100.0),
                        "beta_obs_rad": float(bb),
                        "sigma_rad": float(ss),
                        "beta_inf_rad": float(beta_inf),
                        "z_score_vs_tfpt_beta_inf": zscore,
                    }
                )
        cmb_mean_rad, cmb_sigma_rad = _weighted_mean_and_sigma(x=beta_cmb, sigma=sigma_cmb) if int(z_cmb.shape[0]) >= 1 else (float("nan"), float("nan"))

        # Sanity checks: synthetic truth preference
        checks: list[Check] = []
        checks.append(
            Check(
                check_id="mode",
                passed=True,
                detail=f"mode={used_mode} (configured={mode})",
            )
        )
        checks.append(
            Check(
                check_id="real_data_ingested",
                passed=bool(used_mode == "data"),
                detail="requires mode='data' with non-empty tfpt_suite/data/birefringence_tomography.json['data']['points']; synthetic mode is not considered a verification test",
            )
        )
        checks.append(
            Check(
                check_id="data_families_split",
                passed=True,
                detail=f"n_lowz={int(z_low.shape[0])}, n_cmb={int(z_cmb.shape[0])} (CMB treated as offset check; step-vs-drift fit uses {fit_scope})",
            )
        )
        if used_mode == "synthetic":
            if truth_model_used == "drift":
                checks.append(
                    Check(
                        check_id="prefers_drift_on_drift_synthetic",
                        passed=bool(delta_bic < 0.0),
                        detail=f"truth=drift; ΔBIC(drift-step)={delta_bic:.3f} (negative favors drift); ΔlogL(step-drift)={delta_logL:.3f}",
                    )
                )
            else:
                checks.append(
                    Check(
                        check_id="prefers_step_on_step_synthetic",
                        passed=bool(delta_bic > 0.0),
                        detail=f"truth=step; ΔBIC(drift-step)={delta_bic:.3f} (positive favors step); ΔlogL(step-drift)={delta_logL:.3f}",
                    )
                )
        if used_mode == "data":
            checks.append(
                Check(
                    check_id="tomography_is_not_decisive_by_default",
                    passed=True,
                    detail=f"low-z step/drift Bayes factor is reported for orientation only (no publication-grade covariance; mixed-family points are not used for discrimination). bayes(step/drift)≈{bayes_factor_step_over_drift:.3g}",
                )
            )

        # Plot (optional)
        plot_path = None
        if config.plot:
            try:
                import matplotlib.pyplot as plt

                out_dir = self.output_dir(config)
                out_dir.mkdir(parents=True, exist_ok=True)
                plot_path = out_dir / "beta_tomography.png"

                z_fine = np.linspace(0.0, z_max, 400)
                step_curve = _beta_step(z_fine, beta_inf=beta_inf, z_c=fit_step.params["z_c"], w=fit_step.params["w"])
                drift_curve = _beta_drift(z_fine, beta_inf=beta_inf, p=fit_drift.params["p"])

                fig = plt.figure(figsize=(7, 4))
                ax = fig.add_subplot(1, 1, 1)
                if used_mode == "data" and int(z_low.shape[0]) >= 1 and int(z_cmb.shape[0]) >= 1:
                    ax.errorbar(z_low, beta_low, yerr=sigma_low, fmt="o", ms=3.5, lw=0.8, alpha=0.7, label="low-z (tomography)")
                    ax.errorbar(z_cmb, beta_cmb, yerr=sigma_cmb, fmt="s", ms=4.0, lw=0.8, alpha=0.8, label="CMB (offset checks)")
                else:
                    ax.errorbar(z, beta_obs, yerr=sigma, fmt="o", ms=3.5, lw=0.8, alpha=0.7, label=f"{used_mode} data")
                ax.plot(z_fine, step_curve, lw=2.0, label="best-fit step")
                ax.plot(z_fine, drift_curve, lw=2.0, label="best-fit drift")
                ax.set_xlabel("redshift z")
                ax.set_ylabel("beta(z) [rad]")
                ax.set_title("Birefringence tomography (step vs drift; TFPT beta_inf fixed; calibration marginalized)")
                ax.grid(True, ls=":", alpha=0.4)
                ax.legend(loc="best")
                fig.tight_layout()
                fig.savefig(plot_path, dpi=160)
                plt.close(fig)
            except Exception:
                plot_path = None

        report_lines = [
            "Birefringence tomography (data ingestion + calibration nuisance)",
            "",
            f"beta_inf (TFPT) = varphi0/(4π) = {beta_inf:.10g} rad",
            f"mode = {used_mode} (config file: {cfg_path})",
            f"n = {len(z)}",
            f"z range: [{float(np.min(z)) if len(z) else 0.0}, {z_max}]",
            f"calibration prior: mean={beta_cal_mean:.3e} rad, sigma={beta_cal_sigma:.3e} rad",
            "",
            "Data-family split:",
            f"- low-z tomography points: n={int(z_low.shape[0])}, z∈[{float(np.min(z_low)) if int(z_low.shape[0]) else 0.0:.3g}, {float(np.max(z_low)) if int(z_low.shape[0]) else 0.0:.3g}]",
            f"- CMB points (z≈1100): n={int(z_cmb.shape[0])} (treated as offset checks; no step-vs-drift discrimination)",
            "",
            "Fits (beta_inf fixed):",
            f"- fit_scope = {fit_scope}",
            f"- step: z_c={fit_step.params['z_c']:.4f}, w={fit_step.params['w']:.4f}, logL={fit_step.logL:.3f}, BIC={fit_step.bic:.2f}",
            f"- TFPT-prior step (fixed): z_c={zc_tfpt:.4f}, w={w_tfpt:.4f}, logL={logL_tfpt:.3f}, BIC={bic_tfpt:.2f}",
            f"- drift: p={fit_drift.params['p']:.4f}, logL={fit_drift.logL:.3f}, BIC={fit_drift.bic:.2f}",
            "",
            f"ΔBIC (drift - step) = {delta_bic:.3f}",
            f"ΔlogL (step - drift) = {delta_logL:.3f}",
            f"approx Bayes factor (step/drift) ≈ exp(ΔBIC/2) = {bayes_factor_step_over_drift:.3g}",
            "",
            "CMB offset checks (each is a single z≈1100 constraint; not independent; no tomography):",
        ]
        if cmb_rows:
            for row in cmb_rows:
                report_lines.append(
                    f"- {row['dataset']}: beta={float(row['beta_obs_rad']):.6g}±{float(row['sigma_rad']):.3g} rad; z-score vs TFPT beta_inf = {float(row['z_score_vs_tfpt_beta_inf']):.3g}"
                )
            report_lines += [
                f"- naive weighted mean (orientation only): beta={cmb_mean_rad:.6g}±{cmb_sigma_rad:.3g} rad; z-score={((cmb_mean_rad-beta_inf)/cmb_sigma_rad) if cmb_sigma_rad>0 else float('nan'):.3g}",
                "",
            ]
        else:
            report_lines.append("- (none)")
            report_lines.append("")
        report_lines += [
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
        ]
        if plot_path is not None:
            report_lines += [f"Plot: {plot_path}", ""]

        return ModuleResult(
            results={
                "tfpt": {"beta_inf_rad": beta_inf},
                "config_file": str(cfg_path),
                "mode": used_mode,
                "synthetic_truth": {"truth_model": truth_model_used} if used_mode == "synthetic" else None,
                "data": {
                    "z": [float(v) for v in z],
                    "datasets": datasets,
                    "beta_obs_rad": [float(v) for v in beta_obs],
                    "sigma_rad": [float(v) for v in sigma],
                    "calibration_prior_mean_rad": beta_cal_mean,
                    "calibration_prior_sigma_rad": beta_cal_sigma,
                    "families": {
                        "low_z": {"indices": [int(i) for i, ok in enumerate(lowz_mask.tolist()) if ok], "n": int(z_low.shape[0])},
                        "cmb": {"indices": [int(i) for i, ok in enumerate(cmb_mask.tolist()) if ok], "n": int(z_cmb.shape[0])},
                    },
                },
                "fit_step": {"params": fit_step.params, "logL": fit_step.logL, "bic": fit_step.bic},
                "fit_step_tfpt_prior": {"params": {"z_c": zc_tfpt, "w": w_tfpt}, "logL": logL_tfpt, "bic": bic_tfpt},
                "fit_drift": {"params": fit_drift.params, "logL": fit_drift.logL, "bic": fit_drift.bic},
                "comparison": {
                    "delta_bic_drift_minus_step": delta_bic,
                    "delta_bic_drift_minus_tfpt_prior_step": float(fit_drift.bic - bic_tfpt),
                    "delta_bic_step_minus_tfpt_prior_step": float(fit_step.bic - bic_tfpt),
                    "delta_logL_step_minus_drift": delta_logL,
                    "bayes_factor_step_over_drift": bayes_factor_step_over_drift,
                    "fit_scope": fit_scope,
                },
                "cmb_offset_checks": {
                    "rows": cmb_rows,
                    "naive_weighted_mean_rad": cmb_mean_rad,
                    "naive_weighted_sigma_rad": cmb_sigma_rad,
                },
                "plot": {"beta_tomography_png": str(plot_path) if plot_path is not None else None},
            },
            checks=checks,
            report="\n".join(report_lines),
            warnings=[],
        )

