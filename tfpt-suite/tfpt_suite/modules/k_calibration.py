from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from mpmath import mp
from scipy.integrate import quad

from tfpt_suite.cosmo_scale_map import (
    As_from_ln10_As,
    CosmoScaleInputs,
    MPL_REDUCED_GEV,
    N_reheat_from_rho_ratio,
    a0_over_a_reheat_from_entropy,
    a0_over_a_transition as cosmo_a0_over_a_transition,
    deltaN_from_rho_ratio,
    ell_from_k_hat,
    gev_to_mpc_inv,
    rho_end_GeV4_from_As_r,
    rho_reheat_GeV4,
    starobinsky_N_from_ns,
    starobinsky_r_from_N,
)
from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _plot_k_calibration_scaling(
    *,
    out_dir: Path,
    scaling_s: list[dict[str, float]],
    scaling_t: list[dict[str, float]],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"k_calibration_scaling_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        def xy(rows: list[dict[str, float]]) -> tuple[list[float], list[float]]:
            xs: list[float] = []
            ys: list[float] = []
            for r in rows:
                try:
                    xs.append(float(r["ell_target"]))
                    ys.append(float(r["a0_over_a_transition_needed"]))
                except Exception:
                    continue
            return xs, ys

        x_s, y_s = xy(scaling_s)
        x_t, y_t = xy(scaling_t)

        fig, ax = plt.subplots(figsize=(8.5, 4.5))
        if x_s and y_s:
            ax.plot(x_s, y_s, marker="o", lw=2.0, label="scalar")
        if x_t and y_t:
            ax.plot(x_t, y_t, marker="o", lw=2.0, label="tensor")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"target multipole $\ell$")
        ax.set_ylabel(r"required scaling $a_0/a_t$")
        ax.set_title("k→ℓ calibration: required scaling to place bounce feature")
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(loc="best")
        fig.tight_layout()

        path = out_dir / "k_calibration_scaling.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["k_calibration_scaling_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


FEASIBILITY_N_INFL_HALF_WIDTH = 10.0
FEASIBILITY_N_INFL_POINTS = 41
FEASIBILITY_T_POINTS_MAX = 80


def _plot_k_to_ell_feasibility(
    *,
    out_dir: Path,
    k_bounce_s_raw: float,
    k_bounce_t_raw: float,
    M_GeV: float,
    chi_star_Mpc: float,
    N_inflation_center: float,
    N_reheat: float,
    g_star_s_reheat: float,
    g_star_s_today: float,
    log10_Tmin: float,
    log10_Tmax: float,
    T_points: int,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"k_to_ell_feasibility_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)
        N_min = float(N_inflation_center - FEASIBILITY_N_INFL_HALF_WIDTH)
        N_max = float(N_inflation_center + FEASIBILITY_N_INFL_HALF_WIDTH)
        N_grid = np.linspace(N_min, N_max, FEASIBILITY_N_INFL_POINTS)
        logT_grid = np.linspace(log10_Tmin, log10_Tmax, T_points)
        T_grid = 10.0 ** logT_grid

        ell_s = np.full((logT_grid.size, N_grid.size), np.nan, dtype=float)
        ell_t = np.full_like(ell_s, np.nan)
        for i, logT in enumerate(logT_grid):
            T_reh = float(10.0 ** logT)
            for j, N_infl in enumerate(N_grid):
                a0_over_a_tr, _ = cosmo_a0_over_a_transition(
                    CosmoScaleInputs(
                        transition="inflation_start",
                        N_inflation_from_transition=float(N_infl),
                        N_reheat=float(N_reheat),
                        T_reheat_GeV=T_reh,
                        g_star_s_reheat=float(g_star_s_reheat),
                        g_star_s_today=float(g_star_s_today),
                    )
                )
                if a0_over_a_tr is None or not np.isfinite(a0_over_a_tr) or a0_over_a_tr <= 0:
                    continue
                ell_s[i, j] = ell_from_k_hat(
                    k_hat=float(k_bounce_s_raw),
                    M_GeV=float(M_GeV),
                    a0_over_a_tr=float(a0_over_a_tr),
                    chi_star_Mpc=float(chi_star_Mpc),
                )
                ell_t[i, j] = ell_from_k_hat(
                    k_hat=float(k_bounce_t_raw),
                    M_GeV=float(M_GeV),
                    a0_over_a_tr=float(a0_over_a_tr),
                    chi_star_Mpc=float(chi_star_Mpc),
                )

        fig, axes = plt.subplots(ncols=2, figsize=(12, 5), sharey=True)
        for ax, data, title in [
            (axes[0], ell_s, "Scalar ℓ_bounce"),
            (axes[1], ell_t, "Tensor ℓ_bounce"),
        ]:
            log_ell = np.log10(data)
            cf = ax.contourf(N_grid, logT_grid, log_ell, levels=24, cmap="viridis")
            ax.contour(N_grid, logT_grid, data, levels=[2.0, 2500.0], colors=["lime", "lime"], linewidths=1.5)
            ax.set_title(title)
            ax.set_xlabel(r"$N_{\mathrm{infl}}$")
            ax.grid(True, ls=":", alpha=0.3)
            fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04, label=r"$\log_{10}(\ell_\mathrm{bounce})$")
        axes[0].set_ylabel(r"$\log_{10}(T_{\mathrm{reh}}/\mathrm{GeV})$")
        fig.suptitle("k→ℓ feasibility map (CMB window contour: 2–2500)")
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        path = out_dir / "k_to_ell_feasibility.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        plot["k_to_ell_feasibility_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class FlatLcdm:
    H0_km_s_Mpc: float
    Omega_m: float
    Omega_r: float
    Omega_L: float
    z_star: float


def _chi_star_Mpc(cosmo: FlatLcdm) -> float:
    """
    Comoving distance to last scattering (flat ΛCDM):
      χ(z) = (c/H0) ∫_0^z dz' / E(z')
      E(z)=sqrt(Ω_r(1+z)^4 + Ω_m(1+z)^3 + Ω_Λ)
    """
    c_km_s = 299_792.458
    H0 = float(cosmo.H0_km_s_Mpc)
    pref = c_km_s / H0

    Om = float(cosmo.Omega_m)
    Or = float(cosmo.Omega_r)
    OL = float(cosmo.Omega_L)

    def Ez(z: float) -> float:
        zp1 = 1.0 + z
        return float(np.sqrt(Or * zp1**4 + Om * zp1**3 + OL))

    val, err = quad(lambda z: 1.0 / Ez(z), 0.0, float(cosmo.z_star), epsabs=0.0, epsrel=1e-9, limit=500)
    return float(pref * val)


class KCalibrationModule(TfptModule):
    module_id = "k_calibration"
    title = "k calibration: map bounce k-scale to CMB multipoles (assumption-explicit)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "bounce diagnostics: out/bounce_perturbations/results.json (preferred) or tfpt_suite/data/k_calibration.json (fallback)",
                "TFPT R^2 scale: M/Mpl (to set x=Mη units)",
                "flat ΛCDM distance model for χ_*",
            ],
            outputs=[
                "χ_*(z_*)",
                "naive ℓ_bounce estimate",
                "required a0/a_transition scaling to place bounce features at target ℓ",
                "threshold-derived reheating inputs when available (T_reh, N_reh)",
            ],
            formulas=[
                "x = M η, so dimensionless k = k_com/M in the mode equation in x-units",
                "ℓ ≈ k(Mpc^{-1}) · χ_*(Mpc)",
                "If absolute scale-factor normalization is unknown, infer the required a0/a_transition to map ℓ_bounce to a target ℓ",
            ],
            validation=[
                "produces χ_* and a finite ℓ_bounce estimate",
                "reports required scaling for target multipoles",
            ],
            determinism="Deterministic given inputs (no stochastic sampling).",
            question="How do the dimensionless bounce k-scales map to observable CMB multipoles once we make the scale-factor normalization (a0/a_transition) explicit?",
            objective=[
                "Quantify the missing scale-factor budget needed to place bounce features into a chosen ℓ window (e.g. CMB).",
                "Provide an assumption-explicit a0/a_transition estimate; optionally derive (N, N_reh) via the v1.06 reheating policy.",
            ],
            what_was_done=[
                "Load bounce diagnostics (preferred: live module output) and compute the naive ℓ mapping.",
                "Compute χ_* in a flat ΛCDM snapshot and translate k to ℓ via ℓ≈kχ_*.",
                "Prefer threshold-derived reheating inputs from `cosmo_threshold_history` when present; otherwise use the v1.06 reheating policy.",
                "Optionally derive reheating expansion using `assumptions.reheating_policy_v106` (n_s→N and ρ_end/ρ_reh→N_reh), then compute an entropy-based a0/a_transition estimate.",
            ],
            assumptions=[
                "Flat ΛCDM distance model for χ_* (no full Boltzmann transfer function).",
                "Bounce solver outputs k_bounce in x=Mη units; absolute normalization requires a0/a_transition policy.",
                "Reheating policy (when enabled) uses v1.06 ΔN/N_reh formulas and Planck (n_s, A_s) as external anchors.",
                "If `cosmo_threshold_history` output exists, its threshold-derived T_reh/N_reh override policy inputs.",
            ],
            gaps=[
                "Full TFPT closure still requires a first-principles reheating history; the threshold-derived policy provides a deterministic anchor but remains simplified.",
                "A publication-grade ℓ prediction needs a full transfer-function calculation and explicit observational target policy (CMB vs small-scale probes).",
            ],
            references=[
                "tfpt_suite/cosmo_scale_map.py (entropy mapping + v1.06 reheating helpers)",
                "tfpt_suite/data/k_calibration.json (explicit assumptions incl. reheating_policy_v106)",
                "tfpt_suite/modules/cosmo_threshold_history.py (threshold-derived reheating inputs)",
            ],
            maturity="assumption-explicit bridge (policy-derived; not yet threshold-derived)",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        data_path = Path(__file__).resolve().parent.parent / "data" / "k_calibration.json"
        cfg = json.loads(data_path.read_text(encoding="utf-8"))

        cosmo_raw = dict(cfg.get("cosmology_flat_lcdm", {}))
        Om = float(cosmo_raw["Omega_m"])
        Or = float(cosmo_raw.get("Omega_r", 0.0))
        OL = float(1.0 - Om - Or)
        cosmo = FlatLcdm(
            H0_km_s_Mpc=float(cosmo_raw["H0_km_s_Mpc"]),
            Omega_m=Om,
            Omega_r=Or,
            Omega_L=OL,
            z_star=float(cosmo_raw.get("z_star", 1090.0)),
        )
        chi_star = _chi_star_Mpc(cosmo)

        # Prefer live bounce outputs if present in the current output_dir.
        k_bounce_s_raw = None
        k_bounce_t_raw = None
        k_grid_hat: list[float] | None = None
        bounce_source = None
        bounce_path = config.output_dir / "bounce_perturbations" / "results.json"
        if bounce_path.exists():
            raw = json.loads(bounce_path.read_text(encoding="utf-8"))
            res = raw.get("results", {}) if isinstance(raw, dict) else {}
            diag = res.get("diagnostics", {}) if isinstance(res, dict) else {}
            if isinstance(diag, dict) and "k_bounce_s_est_raw" in diag and "k_bounce_t_est_raw" in diag:
                k_bounce_s_raw = float(diag["k_bounce_s_est_raw"])
                k_bounce_t_raw = float(diag["k_bounce_t_est_raw"])
                bounce_source = f"live output: {bounce_path}"
            if isinstance(res, dict) and "k_grid" in res and isinstance(res.get("k_grid"), list):
                try:
                    k_grid_hat = [float(x) for x in res.get("k_grid", [])]
                except Exception:
                    k_grid_hat = None

        if k_bounce_s_raw is None or k_bounce_t_raw is None:
            fb = dict(cfg.get("fallback_bounce_diagnostics", {}))
            k_bounce_s_raw = float(fb["k_bounce_s_est_raw"])
            k_bounce_t_raw = float(fb["k_bounce_t_est_raw"])
            bounce_source = f"fallback: {data_path}"

        if not k_grid_hat:
            # Bounce module default (k_hat ≡ k/k_bounce in rescaled units).
            k_grid_hat = [0.1, 10.0]

        # Absolute scale: convert M/Mpl to GeV using reduced Planck mass (Mpl_bar).
        # Mpl_bar ≈ 2.435e18 GeV (standard reduced Planck mass; action uses Mpl^2/2).
        Mpl_reduced_GeV = 2.435e18
        M_GeV = float(mp.mpf(c.M_over_Mpl) * mp.mpf(Mpl_reduced_GeV))

        gev_to_mpc = gev_to_mpc_inv()

        # Naive mapping assumes the simulation's absolute scale-factor normalization already corresponds
        # to a0=1 today (this is generally false; we therefore also report the required scaling).
        k_bounce_s_Mpc_inv_naive = float(k_bounce_s_raw * M_GeV * gev_to_mpc)
        k_bounce_t_Mpc_inv_naive = float(k_bounce_t_raw * M_GeV * gev_to_mpc)
        ell_bounce_s_naive = float(k_bounce_s_Mpc_inv_naive * chi_star)
        ell_bounce_t_naive = float(k_bounce_t_Mpc_inv_naive * chi_star)

        assumptions = dict(cfg.get("assumptions", {}))
        a0_over_a_transition = assumptions.get("a0_over_a_transition", None)
        if a0_over_a_transition is not None:
            a0_over_a_transition = float(a0_over_a_transition)
            if a0_over_a_transition <= 0:
                a0_over_a_transition = None

        ell_targets = list(assumptions.get("ell_targets", [2, 30, 700]))
        ell_targets_f = [float(x) for x in ell_targets]

        threshold_meta: dict[str, Any] | None = None
        threshold_path = config.output_dir / "cosmo_threshold_history" / "results.json"

        # Optional: deterministic reheating policy (v1.06) to reduce free parameters in k→ℓ mapping.
        # This derives:
        # - N_pivot from n_s via n_s≈1-2/N,
        # - r≈12/N^2,
        # - ρ_end proxy from (A_s,r),
        # - N_reheat from ρ scaling for a chosen w_reh and T_reh.
        policy_v106 = assumptions.get("reheating_policy_v106", None)
        v106_meta: dict[str, Any] | None = None

        # Cosmo mapping assumptions (explicit): transition label + N + reheating entropy mapping.
        transition = str(assumptions.get("transition", "inflation_start")).strip()
        N_inflation = float(assumptions.get("N_inflation_from_transition", 0.0) or 0.0)
        N_reheat = float(assumptions.get("N_reheat", 0.0) or 0.0)
        T_reh_GeV = float(assumptions.get("T_reheat_GeV", 0.0) or 0.0)
        g_star_s_reh = float(assumptions.get("g_star_s_reheat", 106.75) or 106.75)
        g_star_s_0 = float(assumptions.get("g_star_s_today", 3.91) or 3.91)

        if threshold_path.exists():
            try:
                raw = json.loads(threshold_path.read_text(encoding="utf-8"))
                res = raw.get("results", {}) if isinstance(raw, dict) else {}
                reheating = res.get("reheating", {}) if isinstance(res.get("reheating", {}), dict) else {}
                pivot = res.get("pivot", {}) if isinstance(res.get("pivot", {}), dict) else {}
                if "N_reheat" in reheating:
                    N_reheat = float(reheating.get("N_reheat", N_reheat) or N_reheat)
                if "T_reheat_GeV" in reheating:
                    T_reh_GeV = float(reheating.get("T_reheat_GeV", T_reh_GeV) or T_reh_GeV)
                if "g_star_s_reheat" in reheating:
                    g_star_s_reh = float(reheating.get("g_star_s_reheat", g_star_s_reh) or g_star_s_reh)
                N_pivot_thr = float(pivot.get("N_pivot", N_inflation) or N_inflation)
                if np.isfinite(N_pivot_thr):
                    N_inflation = float(N_pivot_thr)
                threshold_meta = {
                    "source": str(threshold_path),
                    "reheating": reheating,
                    "pivot": pivot,
                }
            except Exception as e:
                threshold_meta = {"source": str(threshold_path), "error": str(e)}

        if isinstance(policy_v106, dict) and bool(policy_v106.get("enabled", False)) and threshold_meta is None:
            try:
                ref_path = Path(__file__).resolve().parent.parent / "data" / "global_reference_minimal.json"
                ref = json.loads(ref_path.read_text(encoding="utf-8"))
                obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
                ns = (
                    float(obs.get("n_s_planck2018", {}).get("mean", float("nan")))
                    if isinstance(obs.get("n_s_planck2018", {}), dict)
                    else float("nan")
                )
                ln10_As = (
                    float(obs.get("ln10_As_planck2018", {}).get("mean", float("nan")))
                    if isinstance(obs.get("ln10_As_planck2018", {}), dict)
                    else float("nan")
                )
                As = As_from_ln10_As(ln10_As)
                N_pivot = starobinsky_N_from_ns(ns)
                r = starobinsky_r_from_N(N_pivot)

                w_reh = float(policy_v106.get("w_reh", 0.0))
                g_star = float(policy_v106.get("g_star_reheat", 120.0))
                c_end = float(policy_v106.get("c_end", 0.35))
                g_star_s_reh = float(policy_v106.get("g_star_s_reheat", g_star_s_reh))
                # Allow the policy to override the canonical T_reh (otherwise keep the existing assumptions value).
                if "T_reheat_GeV" in policy_v106:
                    T_reh_GeV = float(policy_v106.get("T_reheat_GeV", T_reh_GeV) or T_reh_GeV)

                rho_end = rho_end_GeV4_from_As_r(As=As, r=r, c_end=c_end, Mpl_reduced_GeV=MPL_REDUCED_GEV)
                rho_reh = rho_reheat_GeV4(T_reheat_GeV=T_reh_GeV, g_star=g_star)
                N_reheat = N_reheat_from_rho_ratio(w_reh=w_reh, rho_reh=rho_reh, rho_end=rho_end)
                dN = deltaN_from_rho_ratio(w_reh=w_reh, rho_reh=rho_reh, rho_end=rho_end)

                # Override N_inflation with the derived pivot value (explicit policy).
                N_inflation = float(N_pivot)

                v106_meta = {
                    "reference_file": str(ref_path),
                    "n_s": ns,
                    "ln10_As": ln10_As,
                    "A_s": As,
                    "N_pivot": N_pivot,
                    "r": r,
                    "w_reh": w_reh,
                    "g_star_reheat": g_star,
                    "c_end": c_end,
                    "rho_end_GeV4": rho_end,
                    "rho_reh_GeV4": rho_reh,
                    "deltaN": dN,
                }
            except Exception as e:
                v106_meta = {"error": str(e)}

        cosmo_inp = CosmoScaleInputs(
            transition=transition,
            N_inflation_from_transition=N_inflation,
            N_reheat=N_reheat,
            T_reheat_GeV=T_reh_GeV,
            g_star_s_reheat=g_star_s_reh,
            g_star_s_today=g_star_s_0,
        )
        a0_over_a_transition_est, a0_est_note = cosmo_a0_over_a_transition(cosmo_inp)

        a0_over_a_reh_est = (
            a0_over_a_reheat_from_entropy(
                T_reheat_GeV=T_reh_GeV,
                g_star_s_reheat=g_star_s_reh,
                g_star_s_today=g_star_s_0,
                T0_K=cosmo_inp.T0_K,
            )
            if T_reh_GeV > 0
            else None
        )
        a0_over_a_end_est = (
            float(a0_over_a_reh_est * float(np.exp(N_reheat))) if (a0_over_a_reh_est is not None and np.isfinite(a0_over_a_reh_est)) else None
        )

        def scaling_needed(ell_naive: float, ell_target: float) -> dict[str, float]:
            s = ell_naive / ell_target if ell_target > 0 else float("nan")
            N = float(np.log(s)) if s > 0 else float("nan")
            return {"ell_target": ell_target, "a0_over_a_transition_needed": float(s), "N_needed": N}

        scaling_s = [scaling_needed(ell_bounce_s_naive, L) for L in ell_targets_f]
        scaling_t = [scaling_needed(ell_bounce_t_naive, L) for L in ell_targets_f]

        # If user provides a0/a_transition, report calibrated ℓ.
        ell_bounce_s_cal = None
        ell_bounce_t_cal = None
        if a0_over_a_transition is not None:
            ell_bounce_s_cal = ell_bounce_s_naive / a0_over_a_transition
            ell_bounce_t_cal = ell_bounce_t_naive / a0_over_a_transition

        # Also report a calibrated ℓ under the expansion-budget estimate (if provided).
        ell_bounce_s_est = None
        ell_bounce_t_est = None
        if a0_over_a_transition_est is not None and a0_over_a_transition_est > 0:
            ell_bounce_s_est = ell_from_k_hat(
                k_hat=float(k_bounce_s_raw),
                M_GeV=float(M_GeV),
                a0_over_a_tr=float(a0_over_a_transition_est),
                chi_star_Mpc=float(chi_star),
            )
            ell_bounce_t_est = ell_from_k_hat(
                k_hat=float(k_bounce_t_raw),
                M_GeV=float(M_GeV),
                a0_over_a_tr=float(a0_over_a_transition_est),
                chi_star_Mpc=float(chi_star),
            )

        # Acceptance-style diagnostic: does the default k_hat grid map to ℓ≈[2,2500] under the estimate?
        ell_range_est = None
        if a0_over_a_transition_est is not None and a0_over_a_transition_est > 0 and k_grid_hat:
            kh_min = float(min(k_grid_hat))
            kh_max = float(max(k_grid_hat))
            k_hat_grid_source = "bounce_output"
            kh_cfg = assumptions.get("k_hat_grid_for_ell_range", None)
            if isinstance(kh_cfg, list) and len(kh_cfg) == 2:
                try:
                    a = float(kh_cfg[0])
                    b = float(kh_cfg[1])
                    if np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0 and a != b:
                        kh_min = float(min(a, b))
                        kh_max = float(max(a, b))
                        k_hat_grid_source = "config"
                except Exception:
                    pass
            ell_s_min = ell_from_k_hat(k_hat=kh_min * float(k_bounce_s_raw), M_GeV=M_GeV, a0_over_a_tr=a0_over_a_transition_est, chi_star_Mpc=chi_star)
            ell_s_max = ell_from_k_hat(k_hat=kh_max * float(k_bounce_s_raw), M_GeV=M_GeV, a0_over_a_tr=a0_over_a_transition_est, chi_star_Mpc=chi_star)
            ell_t_min = ell_from_k_hat(k_hat=kh_min * float(k_bounce_t_raw), M_GeV=M_GeV, a0_over_a_tr=a0_over_a_transition_est, chi_star_Mpc=chi_star)
            ell_t_max = ell_from_k_hat(k_hat=kh_max * float(k_bounce_t_raw), M_GeV=M_GeV, a0_over_a_tr=a0_over_a_transition_est, chi_star_Mpc=chi_star)
            ell_range_est = {
                "k_hat_grid": [kh_min, kh_max],
                "k_hat_grid_source": k_hat_grid_source,
                "ell_scalar_range": [float(ell_s_min), float(ell_s_max)],
                "ell_tensor_range": [float(ell_t_min), float(ell_t_max)],
            }

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="bounce_scale_loaded",
                passed=bool(k_bounce_s_raw is not None and k_bounce_t_raw is not None),
                detail=f"k_bounce_s_raw={k_bounce_s_raw}, k_bounce_t_raw={k_bounce_t_raw} ({bounce_source})",
            )
        )
        checks.append(
            Check(
                check_id="chi_star_computed",
                passed=bool(np.isfinite(chi_star) and chi_star > 0),
                detail=f"chi_star(z*={cosmo.z_star}) = {chi_star:.3f} Mpc (flat ΛCDM)",
            )
        )
        checks.append(
            Check(
                check_id="ell_bounce_finite",
                passed=bool(np.isfinite(ell_bounce_s_naive) and np.isfinite(ell_bounce_t_naive)),
                detail=f"ell_bounce_naive scalar={ell_bounce_s_naive:.3e}, tensor={ell_bounce_t_naive:.3e}",
            )
        )
        if a0_over_a_transition_est is not None:
            checks.append(
                Check(
                    check_id="expansion_budget_estimate_present",
                    passed=bool(np.isfinite(a0_over_a_transition_est) and a0_over_a_transition_est > 0),
                    detail=f"{a0_est_note}; a0/a_transition(est)={a0_over_a_transition_est:.3e} (N_infl={N_inflation}, N_reh={N_reheat}, T_reh={T_reh_GeV:.3e} GeV)",
                )
            )
        if ell_range_est is not None:
            ell_min_target = 2.0
            ell_max_target = 2500.0
            checks.append(
                Check(
                    check_id="ell_range_covers_cmb_scales_scalar_est",
                    passed=bool(float(ell_range_est["ell_scalar_range"][0]) <= ell_min_target and float(ell_range_est["ell_scalar_range"][1]) >= ell_max_target),
                    detail=f"scalar ell range ≈[{ell_range_est['ell_scalar_range'][0]:.3g}, {ell_range_est['ell_scalar_range'][1]:.3g}] from k_hat∈[{ell_range_est['k_hat_grid'][0]}, {ell_range_est['k_hat_grid'][1]}] (source={ell_range_est.get('k_hat_grid_source')}) vs target [{ell_min_target},{ell_max_target}]",
                )
            )
            checks.append(
                Check(
                    check_id="ell_range_covers_cmb_scales_tensor_est",
                    passed=bool(float(ell_range_est["ell_tensor_range"][0]) <= ell_min_target and float(ell_range_est["ell_tensor_range"][1]) >= ell_max_target),
                    detail=f"tensor ell range ≈[{ell_range_est['ell_tensor_range'][0]:.3g}, {ell_range_est['ell_tensor_range'][1]:.3g}] from k_hat∈[{ell_range_est['k_hat_grid'][0]}, {ell_range_est['k_hat_grid'][1]}] (source={ell_range_est.get('k_hat_grid_source')}) vs target [{ell_min_target},{ell_max_target}]",
                )
            )

        lines: list[str] = []
        lines += [
            "k calibration (assumption-explicit)",
            "",
            f"config: {data_path}",
            f"bounce diagnostics source: {bounce_source}",
            "",
            "Cosmology model (flat ΛCDM):",
            f"- H0 = {cosmo.H0_km_s_Mpc} km/s/Mpc",
            f"- Ω_m = {cosmo.Omega_m}",
            f"- Ω_r = {cosmo.Omega_r}",
            f"- Ω_Λ = {cosmo.Omega_L}",
            f"- z_* = {cosmo.z_star}",
            f"- χ_* = {chi_star:.3f} Mpc",
            "",
            "TFPT scale:",
            f"- M/Mpl = {c.M_over_Mpl}",
            f"- Mpl(reduced) = {Mpl_reduced_GeV:.3e} GeV",
            f"- M = {M_GeV:.3e} GeV",
            "",
            "Bounce scale in x=Mη units (dimensionless):",
            f"- k_bounce_s_raw ≈ {k_bounce_s_raw:.6g}",
            f"- k_bounce_t_raw ≈ {k_bounce_t_raw:.6g}",
            "",
            "Naive projection (assumes absolute normalization already corresponds to a0=1 today):",
            f"- k_bounce_s ≈ {k_bounce_s_Mpc_inv_naive:.3e} Mpc^-1  =>  ℓ_bounce_s ≈ {ell_bounce_s_naive:.3e}",
            f"- k_bounce_t ≈ {k_bounce_t_Mpc_inv_naive:.3e} Mpc^-1  =>  ℓ_bounce_t ≈ {ell_bounce_t_naive:.3e}",
            "",
            "Required overall scaling to place bounce features at target multipoles:",
            "  (interpretation: needed a0/a_transition so that ℓ_bounce = ℓ_bounce_naive/(a0/a_transition))",
            "",
            "Scalar:",
            *[f"- ℓ_target={d['ell_target']:.0f}: a0/a_t needed={d['a0_over_a_transition_needed']:.3e}  (N_needed=ln(a0/a_t)={d['N_needed']:.3f})" for d in scaling_s],
            "",
            "Tensor:",
            *[f"- ℓ_target={d['ell_target']:.0f}: a0/a_t needed={d['a0_over_a_transition_needed']:.3e}  (N_needed=ln(a0/a_t)={d['N_needed']:.3f})" for d in scaling_t],
            "",
        ]
        if a0_over_a_transition_est is not None:
            lines += [
                "Expansion-budget estimate (optional assumptions):",
                f"- transition = {transition}",
                f"- N_inflation_from_transition = {N_inflation:.3f}",
                f"- N_reheat = {N_reheat:.3f} (a_reh/a_end = exp(N_reheat))",
                f"- T_reheat = {T_reh_GeV:.3e} GeV (instantaneous reheating assumed)",
                f"- g*_s(reheat)={g_star_s_reh}, g*_s(today)={g_star_s_0}",
                f"- a0/a_reh ≈ {a0_over_a_reh_est:.3e}" if a0_over_a_reh_est is not None else "- a0/a_reh: n/a",
                f"- a0/a_end ≈ {a0_over_a_end_est:.3e}" if a0_over_a_end_est is not None else "- a0/a_end: n/a",
                f"- a0/a_transition(est) ≈ {a0_over_a_transition_est:.3e}  (N_est=ln(a0/a_t)={np.log(a0_over_a_transition_est):.3f})",
                f"- ℓ_bounce_s(est) ≈ {ell_bounce_s_est:.3e}" if ell_bounce_s_est is not None else "- ℓ_bounce_s(est): n/a",
                f"- ℓ_bounce_t(est) ≈ {ell_bounce_t_est:.3e}" if ell_bounce_t_est is not None else "- ℓ_bounce_t(est): n/a",
                *( [f"- ℓ range (est, scalar) ≈ [{ell_range_est['ell_scalar_range'][0]:.3g}, {ell_range_est['ell_scalar_range'][1]:.3g}]"] if ell_range_est is not None else [] ),
                "",
            ]
        if a0_over_a_transition is not None:
            lines += [
                f"Using provided a0/a_transition = {a0_over_a_transition:.6g}:",
                f"- ℓ_bounce_s(calibrated) = {ell_bounce_s_cal:.6g}",
                f"- ℓ_bounce_t(calibrated) = {ell_bounce_t_cal:.6g}",
                "",
            ]
        if threshold_meta is not None:
            lines += [
                "Threshold-derived reheating inputs (cosmo_threshold_history):",
                json.dumps(threshold_meta, indent=2, sort_keys=True),
                "",
            ]
        if v106_meta is not None:
            lines += [
                "Reheating policy v1.06 (derived inputs; overrides N and N_reheat when enabled):",
                json.dumps(v106_meta, indent=2, sort_keys=True),
                "",
            ]
        lines += [
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"k_calibration_scaling_png": None, "k_to_ell_feasibility_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_k_calibration_scaling(out_dir=out_dir, scaling_s=scaling_s, scaling_t=scaling_t)
            policy_path = Path(__file__).resolve().parent.parent / "data" / "cosmo_reheating_policy_v106.json"
            log10_Tmin = -3.0
            log10_Tmax = 15.0
            T_points = 60
            if policy_path.is_file():
                try:
                    policy = json.loads(policy_path.read_text(encoding="utf-8"))
                    scan = policy.get("scan", {}) if isinstance(policy.get("scan", {}), dict) else {}
                    log10_Tmin = float(scan.get("log10_Tmin_GeV", log10_Tmin))
                    log10_Tmax = float(scan.get("log10_Tmax_GeV", log10_Tmax))
                    n_points = int(scan.get("n_points", T_points))
                    if n_points > 0:
                        T_points = min(n_points, FEASIBILITY_T_POINTS_MAX)
                except Exception:
                    pass
            plot_feas, plot_feas_warnings = _plot_k_to_ell_feasibility(
                out_dir=out_dir,
                k_bounce_s_raw=k_bounce_s_raw,
                k_bounce_t_raw=k_bounce_t_raw,
                M_GeV=M_GeV,
                chi_star_Mpc=chi_star,
                N_inflation_center=N_inflation,
                N_reheat=N_reheat,
                g_star_s_reheat=g_star_s_reh,
                g_star_s_today=g_star_s_0,
                log10_Tmin=log10_Tmin,
                log10_Tmax=log10_Tmax,
                T_points=T_points,
            )
            plot.update(plot_feas)
            warnings.extend(plot_feas_warnings)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "config_file": str(data_path),
                "bounce_source": bounce_source,
                "bounce": {"k_bounce_s_est_raw": k_bounce_s_raw, "k_bounce_t_est_raw": k_bounce_t_raw},
                "cosmology": {
                    "H0_km_s_Mpc": cosmo.H0_km_s_Mpc,
                    "Omega_m": cosmo.Omega_m,
                    "Omega_r": cosmo.Omega_r,
                    "Omega_L": cosmo.Omega_L,
                    "z_star": cosmo.z_star,
                    "chi_star_Mpc": chi_star,
                },
                "constants": {
                    "M_over_Mpl": c.M_over_Mpl,
                    "Mpl_reduced_GeV": Mpl_reduced_GeV,
                    "M_GeV": M_GeV,
                    "GeV_to_Mpc_inv": gev_to_mpc,
                },
                "naive": {
                    "k_bounce_s_Mpc_inv": k_bounce_s_Mpc_inv_naive,
                    "k_bounce_t_Mpc_inv": k_bounce_t_Mpc_inv_naive,
                    "ell_bounce_s": ell_bounce_s_naive,
                    "ell_bounce_t": ell_bounce_t_naive,
                },
                "assumptions": {
                    "transition": transition,
                    "a0_over_a_transition": a0_over_a_transition,
                    "ell_targets": ell_targets_f,
                    "N_inflation_from_transition": N_inflation,
                    "N_reheat": N_reheat,
                    "T_reheat_GeV": T_reh_GeV,
                    "g_star_s_reheat": g_star_s_reh,
                    "g_star_s_today": g_star_s_0,
                    "threshold_history": threshold_meta,
                    "reheating_policy_v106": v106_meta,
                },
                "expansion_budget_estimate": {
                    "transition": transition,
                    "N_inflation_from_transition": N_inflation,
                    "N_reheat": N_reheat,
                    "T_reheat_GeV": T_reh_GeV,
                    "g_star_s_reheat": g_star_s_reh,
                    "g_star_s_today": g_star_s_0,
                    "note": a0_est_note,
                    "a0_over_a_reheat": a0_over_a_reh_est,
                    "a0_over_a_end": a0_over_a_end_est,
                    "a0_over_a_transition": a0_over_a_transition_est,
                    "ell_range_est": ell_range_est,
                    "ell_bounce_s": ell_bounce_s_est,
                    "ell_bounce_t": ell_bounce_t_est,
                },
                "required_scaling": {"scalar": scaling_s, "tensor": scaling_t},
                "calibrated": {
                    "ell_bounce_s": ell_bounce_s_cal,
                    "ell_bounce_t": ell_bounce_t_cal,
                },
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

