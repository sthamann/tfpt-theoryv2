from __future__ import annotations

import json
import math
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from mpmath import mp
from scipy.integrate import quad

from tfpt_suite.cosmo_scale_map import (
    As_from_ln10_As,
    CosmoScaleInputs,
    N_reheat_from_rho_ratio,
    a0_over_a_transition,
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


def _plot_cosmo_history_sampler(
    *,
    out_dir: Path,
    scenario_logs: dict[str, dict[str, list[float]]],
    scenario_counts: dict[str, dict[str, int]],
    cmb_window: tuple[float, float],
    ell_targets: list[float],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {
        "cosmo_history_sampler_log10ell_png": None,
        "cosmo_history_sampler_cmb_fraction_png": None,
    }
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        # --- Plot 1: log10(ell) histograms (scenario × sector) ---
        scenario_ids = list(scenario_logs.keys())
        scenario_ids.sort()

        # Choose a robust x-range based on quantiles across all series.
        all_x: list[float] = []
        for sid in scenario_ids:
            for sec in ("scalar", "tensor"):
                all_x.extend([float(x) for x in scenario_logs.get(sid, {}).get(sec, []) if np.isfinite(float(x))])
        if all_x:
            x_lo = float(np.quantile(np.array(all_x, dtype=float), 0.01))
            x_hi = float(np.quantile(np.array(all_x, dtype=float), 0.99))
            pad = 0.08 * max(1.0, abs(x_hi - x_lo))
            x_lo -= pad
            x_hi += pad
        else:
            x_lo, x_hi = -10.0, 10.0

        fig, axes = plt.subplots(nrows=len(scenario_ids), ncols=2, figsize=(11, 3.6 * max(1, len(scenario_ids))), sharex=True)
        if len(scenario_ids) == 1:
            axes = np.array([axes])  # type: ignore[assignment]

        log_cmb_lo = float(math.log10(float(cmb_window[0])))
        log_cmb_hi = float(math.log10(float(cmb_window[1])))

        for r, sid in enumerate(scenario_ids):
            for c, sec in enumerate(("scalar", "tensor")):
                ax = axes[r, c]
                xs = [float(x) for x in scenario_logs.get(sid, {}).get(sec, []) if np.isfinite(float(x))]
                ax.hist(xs, bins=70, color=("#1f77b4" if sec == "scalar" else "#ff7f0e"), alpha=0.85)

                ax.axvline(log_cmb_lo, color="black", lw=1.0, ls="--", alpha=0.8)
                ax.axvline(log_cmb_hi, color="black", lw=1.0, ls="--", alpha=0.8)

                for t in ell_targets:
                    if t > 0:
                        ax.axvline(float(math.log10(float(t))), color="gray", lw=0.8, ls=":", alpha=0.6)

                ax.set_xlim(x_lo, x_hi)
                ax.set_ylabel("count")
                ax.set_title(f"{sid}: {sec}")
                ax.grid(True, ls=":", alpha=0.35)

                cnt = scenario_counts.get(sid, {})
                if sec == "scalar":
                    k = "scalar_in_cmb"
                    n = "scalar_total"
                else:
                    k = "tensor_in_cmb"
                    n = "tensor_total"
                num = int(cnt.get(k, 0))
                den = int(cnt.get(n, max(1, len(xs))))
                frac = float(num) / float(den) if den > 0 else float("nan")
                ax.text(
                    0.02,
                    0.92,
                    f"in CMB window: {num}/{den} = {frac:.3g}",
                    transform=ax.transAxes,
                    fontsize=10,
                    va="top",
                    ha="left",
                    bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
                )

        axes[-1, 0].set_xlabel(r"$\log_{10}\,\ell_{\mathrm{bounce}}$")
        axes[-1, 1].set_xlabel(r"$\log_{10}\,\ell_{\mathrm{bounce}}$")
        fig.suptitle(r"Bounce multipole location under explicit priors (histograms of $\log_{10}\ell$)")
        fig.tight_layout()

        p1 = out_dir / "cosmo_history_sampler_log10ell.png"
        fig.savefig(p1, dpi=160)
        plt.close(fig)
        plot["cosmo_history_sampler_log10ell_png"] = str(p1)

        # --- Plot 2: CMB-window fractions (by scenario and sector) ---
        fig2, ax2 = plt.subplots(figsize=(8.5, 4.2))
        x = np.arange(len(scenario_ids), dtype=float)
        w = 0.36
        frac_s = []
        frac_t = []
        for sid in scenario_ids:
            cnt = scenario_counts.get(sid, {})
            s = float(cnt.get("scalar_in_cmb", 0)) / float(max(1, int(cnt.get("scalar_total", 1))))
            t = float(cnt.get("tensor_in_cmb", 0)) / float(max(1, int(cnt.get("tensor_total", 1))))
            frac_s.append(s)
            frac_t.append(t)

        ax2.bar(x - w / 2.0, frac_s, width=w, label="scalar", color="#1f77b4", alpha=0.9)
        ax2.bar(x + w / 2.0, frac_t, width=w, label="tensor", color="#ff7f0e", alpha=0.9)
        ax2.set_xticks(x)
        ax2.set_xticklabels(scenario_ids)
        ax2.set_ylim(0.0, 1.0)
        ax2.set_ylabel("fraction in CMB window")
        ax2.set_title(f"CMB-window hit fraction (ℓ∈[{cmb_window[0]:.0f},{cmb_window[1]:.0f}])")
        ax2.grid(True, axis="y", ls=":", alpha=0.35)
        ax2.legend(loc="best")
        fig2.tight_layout()

        p2 = out_dir / "cosmo_history_sampler_cmb_fraction.png"
        fig2.savefig(p2, dpi=160)
        plt.close(fig2)
        plot["cosmo_history_sampler_cmb_fraction_png"] = str(p2)

    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class _Scenario:
    id: str
    note: str
    N_inflation_range: tuple[float, float]
    N_reheat_range: tuple[float, float]
    log10_T_reheat_GeV_range: tuple[float, float]
    samples: int


@dataclass(frozen=True)
class _BestCandidate:
    target_ell: float
    sector: str  # scalar|tensor
    ell_bounce: float
    abs_log10_error: float
    params: dict[str, float]


def _tfpt_suite_dir() -> Path:
    # .../tfpt-suite/unconventional/tfpt_unconventional/modules/<this_file>
    return Path(__file__).resolve().parents[3]


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _chi_star_Mpc_flat_lcdm(*, H0_km_s_Mpc: float, Omega_m: float, Omega_r: float, z_star: float) -> float:
    """
    Comoving distance to last scattering (flat ΛCDM):
      χ(z) = (c/H0) ∫_0^z dz' / E(z')
      E(z)=sqrt(Ω_r(1+z)^4 + Ω_m(1+z)^3 + Ω_Λ)
    """
    c_km_s = 299_792.458
    H0 = float(H0_km_s_Mpc)
    pref = c_km_s / H0

    Om = float(Omega_m)
    Or = float(Omega_r)
    OL = float(1.0 - Om - Or)

    def Ez(z: float) -> float:
        zp1 = 1.0 + z
        return float(np.sqrt(Or * zp1**4 + Om * zp1**3 + OL))

    val, _err = quad(lambda z: 1.0 / Ez(z), 0.0, float(z_star), epsabs=0.0, epsrel=1e-9, limit=500)
    return float(pref * val)


def _quantiles(xs: list[float], qs: list[float]) -> dict[str, float]:
    if not xs:
        return {str(q): float("nan") for q in qs}
    arr = np.array(xs, dtype=float)
    out: dict[str, float] = {}
    for q in qs:
        out[str(q)] = float(np.quantile(arr, q))
    return out


class CosmoHistorySamplerModule(TfptModule):
    module_id = "ux_cosmo_history_sampler"
    title = "Unconventional: cosmology history sampler (k→ℓ feasibility + robustness)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "bounce diagnostics: prefer tfpt-suite/out/bounce_perturbations/results.json; else fallback in tfpt_suite/data/k_calibration.json",
                "cosmology assumptions: tfpt_suite/data/k_calibration.json (H0, Ω_m, Ω_r, z_*)",
                "explicit scale-factor mapping model: tfpt_suite/cosmo_scale_map.py",
            ],
            outputs=[
                "distribution of ℓ_bounce under explicit priors over (N_inflation, N_reheat, T_reheat)",
                "fraction of samples that land in CMB-like ℓ windows (diagnostic)",
                "best candidates (closest to target ℓ in log-space)",
            ],
            formulas=[
                "a0/a_transition = (a0/a_reh) * exp(N_reheat) * exp(N_inflation_from_transition)",
                "a0/a_reh from entropy conservation (instantaneous reheating)",
                "ℓ ≈ k(Mpc^{-1}) · χ_*;  k(Mpc^{-1}) = k_hat * M_GeV * (GeV→Mpc^{-1}) / (a0/a_transition)",
            ],
            validation=[
                "report whether plausible priors allow ℓ_bounce to sit in the CMB window without hidden tuning",
            ],
            determinism="Deterministic given config seed + precision.",
        )

    def run(self, config) -> ModuleResult:
        suite_dir = _tfpt_suite_dir()
        kc_path = suite_dir / "tfpt_suite" / "data" / "k_calibration.json"
        kc = _read_json(kc_path)
        pol_path = suite_dir / "tfpt_suite" / "data" / "cosmo_reheating_policy_v106.json"
        ref_min_path = suite_dir / "tfpt_suite" / "data" / "global_reference_minimal.json"

        cosmo_raw = dict(kc.get("cosmology_flat_lcdm", {}))
        H0 = float(cosmo_raw["H0_km_s_Mpc"])
        Om = float(cosmo_raw["Omega_m"])
        Or = float(cosmo_raw.get("Omega_r", 0.0))
        z_star = float(cosmo_raw.get("z_star", 1090.0))
        chi_star = _chi_star_Mpc_flat_lcdm(H0_km_s_Mpc=H0, Omega_m=Om, Omega_r=Or, z_star=z_star)

        # Bounce diagnostics: prefer the main-suite output folder if present.
        bounce_path = suite_dir / "out" / "bounce_perturbations" / "results.json"
        k_bounce_s_raw: float | None = None
        k_bounce_t_raw: float | None = None
        bounce_source: str
        if bounce_path.exists():
            raw = _read_json(bounce_path)
            res = raw.get("results", {}) if isinstance(raw, dict) else {}
            diag = res.get("diagnostics", {}) if isinstance(res, dict) else {}
            if isinstance(diag, dict) and "k_bounce_s_est_raw" in diag and "k_bounce_t_est_raw" in diag:
                k_bounce_s_raw = float(diag["k_bounce_s_est_raw"])
                k_bounce_t_raw = float(diag["k_bounce_t_est_raw"])
                bounce_source = f"live: {bounce_path}"
            else:
                bounce_source = f"live-but-missing-diagnostics: {bounce_path}"
        else:
            bounce_source = "fallback (bounce_perturbations output not found)"

        if k_bounce_s_raw is None or k_bounce_t_raw is None:
            fb = dict(kc.get("fallback_bounce_diagnostics", {}))
            k_bounce_s_raw = float(fb["k_bounce_s_est_raw"])
            k_bounce_t_raw = float(fb["k_bounce_t_est_raw"])

        # Absolute scale: M from TFPT.
        c = TfptConstants.compute()
        Mpl_reduced_GeV = 2.435e18
        M_GeV = float(mp.mpf(c.M_over_Mpl) * mp.mpf(Mpl_reduced_GeV))

        gev_to_mpc = float(gev_to_mpc_inv())
        ell_naive_s = float(k_bounce_s_raw * M_GeV * gev_to_mpc * chi_star)
        ell_naive_t = float(k_bounce_t_raw * M_GeV * gev_to_mpc * chi_star)

        # Targets from the k_calibration config (keeps this aligned with the core suite).
        assm = dict(kc.get("assumptions", {}))
        transition = str(assm.get("transition", "horizon_exit_of_pivot"))
        ell_targets = [float(x) for x in list(assm.get("ell_targets", [2, 30, 700]))]

        # Policy layer (v1.06 reheating window): make a threshold/policy-driven scenario explicit.
        pol = _read_json(pol_path) if pol_path.exists() else {}
        pol_ass = pol.get("assumptions", {}) if isinstance(pol.get("assumptions", {}), dict) else {}
        pol_scan = pol.get("scan", {}) if isinstance(pol.get("scan", {}), dict) else {}
        w_reh = float(pol_ass.get("w_reh", 0.0))
        g_star = float(pol_ass.get("g_star_reheat", 120.0))
        c_end = float(pol_ass.get("c_end", 0.35))
        log10_Tmin_pol = float(pol_scan.get("log10_Tmin_GeV", -3.0))
        log10_Tmax_pol = float(pol_scan.get("log10_Tmax_GeV", 15.0))
        T_floor_pol = float(pol_ass.get("T_reheat_GeV_floor", 6.0e-3))
        deltaN_floor = float(pol_ass.get("deltaN_floor", -6.0))

        ref_min = _read_json(ref_min_path) if ref_min_path.exists() else {}
        obs = ref_min.get("observables", {}) if isinstance(ref_min.get("observables", {}), dict) else {}
        ns = float(obs.get("n_s_planck2018", {}).get("mean", float("nan"))) if isinstance(obs.get("n_s_planck2018", {}), dict) else float("nan")
        ln10_As = float(obs.get("ln10_As_planck2018", {}).get("mean", float("nan"))) if isinstance(obs.get("ln10_As_planck2018", {}), dict) else float("nan")
        As = As_from_ln10_As(ln10_As)
        N_pivot = starobinsky_N_from_ns(ns)
        r = starobinsky_r_from_N(N_pivot)
        rho_end = rho_end_GeV4_from_As_r(As=As, r=r, c_end=c_end)

        # Explicit priors (three scenarios: “policy_v106”, “plausible” and “extended”), to make the conclusion robust.
        # You can reduce runtime for quick iterations by setting:
        #   TFPT_UX_COSMO_SAMPLES_PLAUSIBLE=200
        #   TFPT_UX_COSMO_SAMPLES_EXTENDED=300
        # (The suite remains deterministic given `seed`.)
        samples_policy = int(os.getenv("TFPT_UX_COSMO_SAMPLES_POLICY", "2500").strip() or "2500")
        samples_plausible = int(os.getenv("TFPT_UX_COSMO_SAMPLES_PLAUSIBLE", "4000").strip() or "4000")
        samples_extended = int(os.getenv("TFPT_UX_COSMO_SAMPLES_EXTENDED", "6000").strip() or "6000")
        scenarios = [
            _Scenario(
                id="policy_v106_threshold_driven",
                note="Threshold/policy-driven: N_pivot from n_s; N_reheat derived from T_reheat via ρ scaling (v1.06 ΔN model).",
                N_inflation_range=(float(N_pivot), float(N_pivot)),
                N_reheat_range=(0.0, 0.0),
                log10_T_reheat_GeV_range=(log10_Tmin_pol, log10_Tmax_pol),
                samples=samples_policy,
            ),
            _Scenario(
                id="plausible",
                note="Conservative ranges around standard inflation+reheating expectations.",
                N_inflation_range=(50.0, 70.0),
                N_reheat_range=(0.0, 15.0),
                log10_T_reheat_GeV_range=(-2.0, 15.0),  # 10 MeV .. 1e15 GeV
                samples=samples_plausible,
            ),
            _Scenario(
                id="extended",
                note="Extended ranges to test whether CMB visibility requires extreme expansion budgets.",
                N_inflation_range=(50.0, 140.0),
                N_reheat_range=(0.0, 40.0),
                log10_T_reheat_GeV_range=(-2.0, 15.0),
                samples=samples_extended,
            ),
        ]

        cmb_window = (2.0, 2500.0)

        scenario_rows: list[dict[str, object]] = []
        scenario_logs: dict[str, dict[str, list[float]]] = {}
        scenario_counts: dict[str, dict[str, int]] = {}
        all_checks: list[Check] = []

        for sc in scenarios:
            ell_s_list: list[float] = []
            ell_t_list: list[float] = []

            in_cmb_s = 0
            in_cmb_t = 0

            # Track best “closest to target” candidates (log-space).
            best: list[_BestCandidate] = []

            for _ in range(int(sc.samples)):
                logT = float(
                    sc.log10_T_reheat_GeV_range[0]
                    + (sc.log10_T_reheat_GeV_range[1] - sc.log10_T_reheat_GeV_range[0]) * random.random()
                )
                T_reh = float(10.0**logT)

                # BBN sanity floor (avoid meaningless tiny reheating).
                T_floor = float(max(1e-2, T_floor_pol))
                if T_reh < T_floor:
                    T_reh = T_floor

                if sc.id == "policy_v106_threshold_driven":
                    # Derive N_infl and N_reheat from the v1.06 policy model.
                    N_infl = float(N_pivot)
                    rho_reh = rho_reheat_GeV4(T_reheat_GeV=float(T_reh), g_star=float(g_star))
                    N_reh = float(N_reheat_from_rho_ratio(w_reh=float(w_reh), rho_reh=float(rho_reh), rho_end=float(rho_end)))
                    # Optional ΔN floor gate: skip points outside the v1.06 reheating window.
                    dN = float(deltaN_from_rho_ratio(w_reh=float(w_reh), rho_reh=float(rho_reh), rho_end=float(rho_end)))
                    if math.isfinite(dN) and dN < deltaN_floor:
                        continue
                    if not (math.isfinite(N_reh) and N_reh >= 0):
                        continue
                else:
                    N_infl = float(sc.N_inflation_range[0] + (sc.N_inflation_range[1] - sc.N_inflation_range[0]) * random.random())
                    N_reh = float(sc.N_reheat_range[0] + (sc.N_reheat_range[1] - sc.N_reheat_range[0]) * random.random())

                inp = CosmoScaleInputs(
                    transition=transition,
                    N_inflation_from_transition=N_infl,
                    N_reheat=N_reh,
                    T_reheat_GeV=T_reh,
                    g_star_s_reheat=float(assm.get("g_star_s_reheat", 106.75) or 106.75),
                    g_star_s_today=float(assm.get("g_star_s_today", 3.91) or 3.91),
                )
                a0_over_at, _note = a0_over_a_transition(inp)
                if a0_over_at is None or not math.isfinite(float(a0_over_at)) or float(a0_over_at) <= 0:
                    continue

                ell_s = float(
                    ell_from_k_hat(
                        k_hat=float(k_bounce_s_raw),
                        M_GeV=float(M_GeV),
                        a0_over_a_tr=float(a0_over_at),
                        chi_star_Mpc=float(chi_star),
                    )
                )
                ell_t = float(
                    ell_from_k_hat(
                        k_hat=float(k_bounce_t_raw),
                        M_GeV=float(M_GeV),
                        a0_over_a_tr=float(a0_over_at),
                        chi_star_Mpc=float(chi_star),
                    )
                )

                if not (math.isfinite(ell_s) and math.isfinite(ell_t) and ell_s > 0 and ell_t > 0):
                    continue

                ell_s_list.append(float(ell_s))
                ell_t_list.append(float(ell_t))

                if cmb_window[0] <= ell_s <= cmb_window[1]:
                    in_cmb_s += 1
                if cmb_window[0] <= ell_t <= cmb_window[1]:
                    in_cmb_t += 1

                for tgt in ell_targets:
                    for sector, ell_val in (("scalar", ell_s), ("tensor", ell_t)):
                        abs_log_err = abs(math.log10(float(ell_val)) - math.log10(float(tgt)))
                        best.append(
                            _BestCandidate(
                                target_ell=float(tgt),
                                sector=sector,
                                ell_bounce=float(ell_val),
                                abs_log10_error=float(abs_log_err),
                                params={
                                    "N_inflation_from_transition": float(N_infl),
                                    "N_reheat": float(N_reh),
                                    "T_reheat_GeV": float(T_reh),
                                    "a0_over_a_transition": float(a0_over_at),
                                },
                            )
                        )

            # Keep only the global top candidates per (target, sector).
            best_by_key: dict[tuple[float, str], _BestCandidate] = {}
            for cand in best:
                key = (float(cand.target_ell), str(cand.sector))
                cur = best_by_key.get(key)
                if cur is None or float(cand.abs_log10_error) < float(cur.abs_log10_error):
                    best_by_key[key] = cand
            best_list = list(best_by_key.values())
            best_list.sort(key=lambda x: (x.sector, x.target_ell))

            # Summary stats (use log10 ℓ quantiles; ℓ itself spans huge ranges).
            log_ell_s = [math.log10(x) for x in ell_s_list if x > 0]
            log_ell_t = [math.log10(x) for x in ell_t_list if x > 0]
            scenario_logs[str(sc.id)] = {"scalar": log_ell_s, "tensor": log_ell_t}
            scenario_counts[str(sc.id)] = {
                "scalar_in_cmb": int(in_cmb_s),
                "tensor_in_cmb": int(in_cmb_t),
                "scalar_total": int(len(ell_s_list)),
                "tensor_total": int(len(ell_t_list)),
            }

            rows = {
                "scenario": sc.__dict__,
                "kept_samples": int(min(len(ell_s_list), len(ell_t_list))),
                "cmb_window": {"min": cmb_window[0], "max": cmb_window[1]},
                "counts": {
                    "scalar_in_cmb": int(in_cmb_s),
                    "tensor_in_cmb": int(in_cmb_t),
                    "scalar_total": int(len(ell_s_list)),
                    "tensor_total": int(len(ell_t_list)),
                },
                "log10_ell_quantiles": {
                    "scalar": _quantiles(log_ell_s, [0.05, 0.5, 0.95]),
                    "tensor": _quantiles(log_ell_t, [0.05, 0.5, 0.95]),
                },
                "ell_minmax": {
                    "scalar": {"min": float(min(ell_s_list)) if ell_s_list else float("nan"), "max": float(max(ell_s_list)) if ell_s_list else float("nan")},
                    "tensor": {"min": float(min(ell_t_list)) if ell_t_list else float("nan"), "max": float(max(ell_t_list)) if ell_t_list else float("nan")},
                },
                "best_candidates": best_list,
            }
            scenario_rows.append(rows)

            # Scenario-level checks:
            frac_s = float(in_cmb_s) / float(len(ell_s_list)) if ell_s_list else 0.0
            frac_t = float(in_cmb_t) / float(len(ell_t_list)) if ell_t_list else 0.0
            all_checks.append(
                Check(
                    check_id=f"{sc.id}_has_any_scalar_in_cmb_window",
                    passed=bool(in_cmb_s > 0),
                    detail=f"fraction={frac_s:.4g} ({in_cmb_s}/{len(ell_s_list)}) in ℓ∈[{cmb_window[0]},{cmb_window[1]}]",
                )
            )
            all_checks.append(
                Check(
                    check_id=f"{sc.id}_has_any_tensor_in_cmb_window",
                    passed=bool(in_cmb_t > 0),
                    detail=f"fraction={frac_t:.4g} ({in_cmb_t}/{len(ell_t_list)}) in ℓ∈[{cmb_window[0]},{cmb_window[1]}]",
                )
            )

        # Global sanity checks
        checks = [
            Check(
                check_id="bounce_source_resolved",
                passed=True,
                detail=f"bounce_source={bounce_source}; k_bounce_s_raw={k_bounce_s_raw}; k_bounce_t_raw={k_bounce_t_raw}",
            ),
            Check(
                check_id="chi_star_computed",
                passed=bool(np.isfinite(chi_star) and chi_star > 0),
                detail=f"chi_star(z*={z_star})={chi_star:.3f} Mpc (flat ΛCDM)",
            ),
            Check(
                check_id="naive_ell_is_huge_as_expected",
                passed=bool(ell_naive_s > 1e40 and ell_naive_t > 1e40),
                detail=f"ell_naive scalar={ell_naive_s:.3e}, tensor={ell_naive_t:.3e}",
            ),
            *all_checks,
        ]

        # Human report (compact; detailed candidates are in results.json).
        lines: list[str] = []
        lines += [
            "Unconventional: cosmology history sampler (k→ℓ feasibility)",
            "",
            "Purpose:",
            "- Explore whether the bounce k-scale can plausibly map into CMB multipoles once an explicit a0/a_transition model is declared.",
            "- This does *not* derive the expansion history; it quantifies feasibility under explicit priors.",
            "",
            f"Inputs:",
            f"- k_calibration config: {kc_path}",
            f"- bounce diagnostics: {bounce_source}",
            f"- k_bounce_s_raw={k_bounce_s_raw:.6g}, k_bounce_t_raw={k_bounce_t_raw:.6g}",
            f"- M = {M_GeV:.3e} GeV (from TFPT M/Mpl and Mpl(reduced)={Mpl_reduced_GeV:.3e} GeV)",
            f"- χ_* = {chi_star:.3f} Mpc (flat ΛCDM; H0={H0}, Ωm={Om}, Ωr={Or}, z*={z_star})",
            "",
            "Naive projection (no absolute normalization):",
            f"- ℓ_bounce_s(naive) ≈ {ell_naive_s:.3e}",
            f"- ℓ_bounce_t(naive) ≈ {ell_naive_t:.3e}",
            "",
            f"CMB window used for diagnostics: ℓ∈[{cmb_window[0]},{cmb_window[1]}]",
            "",
        ]
        for row in scenario_rows:
            sc = row["scenario"]
            counts = row["counts"]
            q = row["log10_ell_quantiles"]
            lines += [
                f"Scenario: {sc['id']}",
                f"- note: {sc['note']}",
                f"- priors: N_infl∈{tuple(sc['N_inflation_range'])}, N_reh∈{tuple(sc['N_reheat_range'])}, log10(Treh/GeV)∈{tuple(sc['log10_T_reheat_GeV_range'])}",
                f"- kept_samples: {row['kept_samples']}",
                f"- scalar in CMB window: {counts['scalar_in_cmb']}/{counts['scalar_total']}",
                f"- tensor in CMB window: {counts['tensor_in_cmb']}/{counts['tensor_total']}",
                f"- log10 ℓ quantiles (scalar): {q['scalar']}",
                f"- log10 ℓ quantiles (tensor): {q['tensor']}",
                "",
            ]

        lines += [
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- If even the extended scenario yields zero CMB hits, TFPT bounce features are likely a small-scale prediction rather than a CMB feature.",
            "- If the plausible scenario yields CMB hits, the next step is not “more tuning”, but encoding a publication-grade expansion history policy (reheating + g*(T) + thresholds) as a deterministic module.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"cosmo_history_sampler_log10ell_png": None, "cosmo_history_sampler_cmb_fraction_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_cosmo_history_sampler(
                out_dir=out_dir,
                scenario_logs=scenario_logs,
                scenario_counts=scenario_counts,
                cmb_window=cmb_window,
                ell_targets=ell_targets,
            )
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "config_file": str(kc_path),
                "bounce_source": bounce_source,
                "bounce": {"k_bounce_s_raw": float(k_bounce_s_raw), "k_bounce_t_raw": float(k_bounce_t_raw)},
                "cosmology": {"H0_km_s_Mpc": H0, "Omega_m": Om, "Omega_r": Or, "z_star": z_star, "chi_star_Mpc": chi_star},
                "constants": {"M_over_Mpl": c.M_over_Mpl, "Mpl_reduced_GeV": Mpl_reduced_GeV, "M_GeV": M_GeV, "GeV_to_Mpc_inv": gev_to_mpc},
                "naive": {"ell_bounce_s": ell_naive_s, "ell_bounce_t": ell_naive_t},
                "diagnostics": {"cmb_window": {"min": cmb_window[0], "max": cmb_window[1]}, "ell_targets": ell_targets, "transition": transition},
                "scenarios": scenario_rows,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

