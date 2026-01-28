from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from mpmath import mp
import numpy as np

from tfpt_suite.alpha_running import AlphaRunningInputs, alpha_bar5_MZ_from_alpha0, alpha_running_inputs_from_pdg
from tfpt_suite.constants import TfptConstants
from tfpt_suite.defect_partition import derive_delta2_from_defect_partition
from tfpt_suite.module_base import (
    Check,
    ModuleResult,
    ModuleSpec,
    TfptModule,
    mk_check_fail,
    mk_check_info,
    mk_check_pass,
    mk_check_warn,
)


def _plot_chi2_contributions(
    *,
    out_dir: Path,
    terms: list[dict[str, Any]],
    r_enabled: bool,
    chi2_r_proxy: mp.mpf,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"chi2_contributions_png": None, "chi2_contributions_no_alpha_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        def _float(x: object) -> float:
            try:
                return float(x)  # mp.mpf, numpy scalar, etc.
            except Exception:
                return float("nan")

        rows: list[dict[str, object]] = []
        for t in terms:
            rows.append(
                {
                    "name": str(t.get("name", "")),
                    "chi2_strict": _float(t.get("chi2_total", 0.0)),
                    "chi2_engineering": _float(t.get("chi2_eff", 0.0)),
                }
            )
        if r_enabled:
            rows.append(
                {
                    "name": "r_upper_95 (proxy)",
                    "chi2_strict": _float(chi2_r_proxy),
                    "chi2_engineering": _float(chi2_r_proxy),
                }
            )

        # Sort for readability by strict chi2
        rows.sort(key=lambda r: float(r.get("chi2_strict", 0.0)), reverse=True)
        labels = [str(r["name"]) for r in rows]
        chi2_strict = [max(float(r["chi2_strict"]), 1e-30) for r in rows]
        chi2_eng = [max(float(r["chi2_engineering"]), 1e-30) for r in rows]

        # Plot 1: strict vs engineering (log x, because alpha can dominate massively)
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 7), sharey=True)
        y = list(range(len(labels)))

        ax1.barh(y, chi2_strict, color="#1f77b4", alpha=0.9)
        ax1.set_xscale("log")
        ax1.set_title("Global consistency χ² contributions (strict)")
        ax1.set_xlabel("χ² (log scale)")
        ax1.grid(True, axis="x", ls=":", alpha=0.4)

        ax2.barh(y, chi2_eng, color="#ff7f0e", alpha=0.9)
        ax2.set_xscale("log")
        ax2.set_title("Global consistency χ² contributions (engineering; sigma_floor)")
        ax2.set_xlabel("χ² (log scale)")
        ax2.grid(True, axis="x", ls=":", alpha=0.4)

        ax2.set_yticks(y)
        ax2.set_yticklabels(labels)
        ax2.invert_yaxis()

        fig.tight_layout()
        p1 = out_dir / "chi2_contributions.png"
        fig.savefig(p1, dpi=160)
        plt.close(fig)
        plot["chi2_contributions_png"] = str(p1)

        # Plot 2: excluding alpha (linear; shows the non-alpha sector structure)
        alpha_keys = {"alpha_inv_0", "alpha_bar5_inv_MZ"}
        rows_wo_alpha = [r for r in rows if str(r.get("name", "")) not in alpha_keys]
        if rows_wo_alpha:
            labels2 = [str(r["name"]) for r in rows_wo_alpha]
            chi2_strict2 = [max(float(r["chi2_strict"]), 0.0) for r in rows_wo_alpha]
            chi2_eng2 = [max(float(r["chi2_engineering"]), 0.0) for r in rows_wo_alpha]
            y2 = list(range(len(labels2)))

            fig2, (bx1, bx2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 6), sharey=True)
            bx1.barh(y2, chi2_strict2, color="#1f77b4", alpha=0.9)
            bx1.set_title("χ² contributions excluding α (strict)")
            bx1.set_xlabel("χ²")
            bx1.grid(True, axis="x", ls=":", alpha=0.4)

            bx2.barh(y2, chi2_eng2, color="#ff7f0e", alpha=0.9)
            bx2.set_title("χ² contributions excluding α (engineering)")
            bx2.set_xlabel("χ²")
            bx2.grid(True, axis="x", ls=":", alpha=0.4)

            bx2.set_yticks(y2)
            bx2.set_yticklabels(labels2)
            bx2.invert_yaxis()

            fig2.tight_layout()
            p2 = out_dir / "chi2_contributions_no_alpha.png"
            fig2.savefig(p2, dpi=160)
            plt.close(fig2)
            plot["chi2_contributions_no_alpha_png"] = str(p2)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


def _plot_global_radar(
    *,
    out_dir: Path,
    terms: list[dict[str, Any]],
    extra_terms: list[dict[str, Any]],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"global_radar_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)
        term_map = {str(t.get("name", "")): t for t in terms}
        extra_map = {str(t.get("name", "")): t for t in extra_terms}
        labels = ["alpha_inv_0", "beta_deg", "cabibbo_lambda", "n_s", "A_s", "ckm_chi2", "pmns_chi2"]

        values: list[float] = []
        for label in labels:
            if label in term_map:
                z = term_map[label].get("z_total", mp.mpf("nan"))
                values.append(float(abs(z)) if mp.isfinite(z) else 0.0)
            elif label == "ckm_chi2" and "ckm_full_pipeline_chi2" in extra_map:
                chi2_val = extra_map["ckm_full_pipeline_chi2"].get("chi2_total", mp.mpf("nan"))
                dof_val = extra_map["ckm_full_pipeline_chi2"].get("dof", 1)
                if mp.isfinite(chi2_val) and int(dof_val) > 0:
                    values.append(float(mp.sqrt(mp.mpf(chi2_val) / mp.mpf(dof_val))))
                else:
                    values.append(0.0)
            elif label == "pmns_chi2" and "pmns_full_pipeline_chi2" in extra_map:
                chi2_val = extra_map["pmns_full_pipeline_chi2"].get("chi2_total", mp.mpf("nan"))
                dof_val = extra_map["pmns_full_pipeline_chi2"].get("dof", 1)
                if mp.isfinite(chi2_val) and int(dof_val) > 0:
                    values.append(float(mp.sqrt(mp.mpf(chi2_val) / mp.mpf(dof_val))))
                else:
                    values.append(0.0)
            else:
                values.append(0.0)

        # Radar plot
        angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
        values_cycle = values + values[:1]
        angles_cycle = angles + angles[:1]

        fig = plt.figure(figsize=(7.0, 6.0))
        ax = fig.add_subplot(111, polar=True)
        ax.plot(angles_cycle, values_cycle, color="#1f77b4", linewidth=2.0)
        ax.fill(angles_cycle, values_cycle, color="#1f77b4", alpha=0.25)
        ax.set_thetagrids([a * 180 / np.pi for a in angles], labels)
        ax.set_rlabel_position(0)
        rmax = max(max(values_cycle), 5.0)
        ax.set_ylim(0, rmax)
        for ring in [2.0, 5.0]:
            ax.plot(angles_cycle, [ring] * len(angles_cycle), color="gray", linestyle="--", linewidth=1.0)
        ax.set_title("Global z-score radar (|z| or √(χ²/dof))")
        fig.tight_layout()

        path = out_dir / "global_radar.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        plot["global_radar_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")
    return plot, warnings


def _chi2_sf(*, chi2: mp.mpf, dof: int) -> mp.mpf:
    """
    Chi-square survival function (upper tail probability):
      p = Q(k/2, chi2/2) = Γ(k/2, chi2/2) / Γ(k/2)
    """
    k = int(dof)
    if k <= 0:
        return mp.mpf("nan")
    x = mp.mpf(chi2) / 2
    s = mp.mpf(k) / 2
    try:
        return mp.gammainc(s, x, mp.inf) / mp.gamma(s)
    except Exception:
        return mp.mpf("nan")


def _read_results_json(*, out_dir: Path, module_id: str) -> dict[str, Any] | None:
    path = out_dir / module_id / "results.json"
    try:
        if not path.is_file():
            return None
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            return payload
    except Exception:
        return None
    return None


class GlobalConsistencyTestModule(TfptModule):
    module_id = "global_consistency_test"
    title = "Global consistency test (multi-observable χ² using a reference table)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["reference table: tfpt_suite/data/global_reference.json"],
            outputs=["per-observable χ² contributions, total χ², watchlist"],
            formulas=[
                "χ² = Σ_i ((pred_i - mean_i)/sigma_i)^2 for Gaussian observables",
                "one-sided r bound handled as: χ²_r = 0 if r<=r_upper else ((r-r_upper)/r_upper)^2 (proxy)",
            ],
            validation=["produces a deterministic score and a sorted contributions table"],
            determinism="Deterministic given the reference table and TFPT invariants.",
            question="Does the current TFPT output set form a self-consistent scorecard once we include the hard sectors (α(0), CKM, PMNS) and enforce explicit FAIL/WARN semantics?",
            objective=[
                "Provide two verification views: engineering (narrow, avoids misleading false negatives) and physics (strict, flags large deviations).",
                "Report per-term χ² contributions plus approximate p-values (dashboard-style; no covariance).",
            ],
            what_was_done=[
                "Compute TFPT predictions for a small set of reference observables and evaluate a reference-table χ².",
                "In physics mode, enable α(0) metrology and import CKM/PMNS χ² from their module outputs to prevent “everything green” optics.",
                "Emit WARN/FAIL checks when |z|>2 or |z|>5 (requested severity semantics).",
            ],
            assumptions=[
                "Reference table is treated as independent Gaussians (no covariance).",
                "r bound is handled by a simple one-sided χ² proxy (not a full likelihood).",
                "Imported CKM/PMNS χ² are diagnostic and depend on scheme/scale policy in the respective pipelines.",
            ],
            gaps=[
                "Not a publication-grade global likelihood fit (no covariance; simplified r proxy).",
                "Hard-sector closure requires finishing finite matching pieces and topology→flavor phase derivations; physics mode intentionally FAILs until that is done.",
            ],
            references=[
                "tfpt_suite/data/global_reference.json (scorecard table; α(0) disabled by default)",
                "tfpt_suite/modules/ckm_full_pipeline.py and pmns_full_pipeline.py (hard-sector χ² imports in physics mode)",
            ],
            maturity="dashboard scorecard (engineering/physics modes; not a full global fit)",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        ref_path = Path(__file__).resolve().parent.parent / "data" / "global_reference.json"
        ref = json.loads(ref_path.read_text(encoding="utf-8"))
        obs: dict[str, Any] = ref.get("observables", {})

        mode = str(getattr(config, "verification_mode", "engineering") or "engineering").strip() or "engineering"
        if mode not in {"engineering", "physics"}:
            mode = "engineering"

        # Physics-mode override: include α(0) metrology and hard sectors in the scorecard.
        # (global_reference.json keeps α(0) disabled by default to avoid false "everything broken" optics.)
        if mode == "physics" and isinstance(obs.get("alpha_inv_0", None), dict):
            obs["alpha_inv_0"] = dict(obs.get("alpha_inv_0", {}))
            obs["alpha_inv_0"]["enabled"] = True

        # --- TFPT predictions (paper v2.4 formulas) ---
        # alpha^{-1} (self-consistent backreaction closure)
        varphi_tree = c.varphi0_tree
        delta_top = c.delta_top
        d2 = derive_delta2_from_defect_partition(delta_top=mp.mpf(delta_top))

        def varphi0_of(alpha: mp.mpf, *, mode: str) -> mp.mpf:
            """
            Backreaction response model for varphi0(α).

            - mode="single_defect" (paper v2.4): varphi = varphi_tree + delta_top * exp(-2α)
            - mode="two_defect": add the next term δ₂ * exp(-4α), where δ₂ is *derived* from a discrete
              defect-partition model (see `defect_partition_derivation` / `tfpt_suite/defect_partition.py`);
              no continuous fit parameter is introduced.
            """
            if mode == "two_defect":
                return varphi_tree + delta_top * mp.e ** (-mp.mpf(2) * alpha) + mp.mpf(d2.delta2) * mp.e ** (-mp.mpf(4) * alpha)
            return varphi_tree + delta_top * mp.e ** (-mp.mpf(2) * alpha)

        def cfe(alpha: mp.mpf, varphi0: mp.mpf) -> mp.mpf:
            return alpha**3 - mp.mpf(2) * (c.c3**3) * alpha**2 - mp.mpf(8) * c.b1 * (c.c3**6) * mp.log(mp.mpf(1) / varphi0)

        def solve_cfe_for(varphi0: mp.mpf) -> mp.mpf:
            f = lambda a: cfe(a, varphi0)
            return mp.findroot(f, (mp.mpf("0.006"), mp.mpf("0.010")))

        def solve_self_consistent_alpha(*, mode: str) -> mp.mpf:
            a = solve_cfe_for(c.varphi0)
            for _ in range(60):
                v = varphi0_of(a, mode=mode)
                nxt = solve_cfe_for(v)
                if abs(nxt - a) < mp.mpf("1e-30"):
                    a = nxt
                    break
                a = nxt
            return a

        alpha_single = solve_self_consistent_alpha(mode="single_defect")
        alpha_two = solve_self_consistent_alpha(mode="two_defect")
        alpha_inv_pred_single = mp.mpf(1) / alpha_single
        alpha_inv_pred_two_defect = mp.mpf(1) / alpha_two

        # α policy:
        # - We keep TFPT's native prediction as α(0) (on-shell, Thomson limit) from the CFE+backreaction closure.
        # - For the primary comparison under the MSbar-at-MZ policy, we map α(0) → ᾱ^(5)(MZ) using an explicit
        #   SM/QED running layer (delta_alpha inputs).
        alpha_inv_0_pred_single = alpha_inv_pred_single
        alpha_inv_0_pred_two_defect = alpha_inv_pred_two_defect

        # Default TFPT α(0) prediction used for downstream comparisons (parameter-free refinement).
        alpha_inv_0_pred = alpha_inv_0_pred_two_defect

        # External running inputs (comparison layer)
        pdg_path = Path(__file__).resolve().parent.parent / "data" / "alpha_running_pdg.json"
        pdg = json.loads(pdg_path.read_text(encoding="utf-8"))
        alpha_inputs = alpha_running_inputs_from_pdg(pdg=pdg)

        alpha0_pred = mp.mpf(1) / alpha_inv_0_pred
        alpha_mz_pred, _dal = alpha_bar5_MZ_from_alpha0(alpha0=alpha0_pred, inputs=alpha_inputs)
        alpha_bar5_inv_MZ_pred = mp.mpf(1) / alpha_mz_pred

        beta_deg_pred = c.beta_deg

        cabibbo_pred = mp.sqrt(c.varphi0) * (mp.mpf(1) - mp.mpf(1) / 2 * c.varphi0)

        # Starobinsky predictions at N=56 (paper Theorem inflation)
        N = mp.mpf(56)
        ns_pred = mp.mpf(1) - mp.mpf(2) / N
        r_pred = mp.mpf(12) / (N**2)
        As_pred = (N**2) / (mp.mpf(24) * mp.pi**2) * (c.M_over_Mpl**2)

        predictions: dict[str, mp.mpf] = {
            "alpha_inv_0": alpha_inv_0_pred,
            "alpha_bar5_inv_MZ": alpha_bar5_inv_MZ_pred,
            "beta_deg": beta_deg_pred,
            "cabibbo_lambda": cabibbo_pred,
            "n_s": ns_pred,
            "r": r_pred,
            "A_s": As_pred,
        }

        # Optional covariance-based gate (alpha-sector only; see global_reference_cov.json).
        cov_gate: dict[str, Any] | None = None
        cov_gate_check: Check | None = None
        cov_path = Path(__file__).resolve().parent.parent / "data" / "global_reference_cov.json"
        if cov_path.is_file():
            try:
                cov_raw = json.loads(cov_path.read_text(encoding="utf-8"))
                labels = cov_raw.get("labels", [])
                cov = cov_raw.get("covariance", [])
                if not isinstance(labels, list) or not isinstance(cov, list) or len(labels) < 1:
                    raise ValueError("invalid labels/covariance")
                n = int(len(labels))
                if not (isinstance(cov, list) and len(cov) == n and all(isinstance(r, list) and len(r) == n for r in cov)):
                    raise ValueError("covariance must be an NxN list")

                # Build residual vector (pred-mean) using the global_reference means.
                rvec: list[mp.mpf] = []
                for name in labels:
                    nm = str(name)
                    if nm not in predictions:
                        raise KeyError(f"covariance label not in predictions: {nm}")
                    cfg = obs.get(nm, None)
                    if not isinstance(cfg, dict) or "mean" not in cfg:
                        raise KeyError(f"covariance label missing mean in global_reference.json: {nm}")
                    mean = mp.mpf(str(cfg["mean"]))
                    rvec.append(mp.mpf(predictions[nm]) - mean)

                # Invert covariance (currently only used for small alpha-sector matrices).
                # Keep it explicit (2x2 analytic inverse) to avoid pulling numpy into the scorecard core.
                if n == 1:
                    s2 = mp.mpf(str(cov[0][0]))
                    if s2 <= 0:
                        raise ValueError("covariance[0][0] must be positive")
                    chi2_cov = (rvec[0] ** 2) / s2
                elif n == 2:
                    a = mp.mpf(str(cov[0][0]))
                    b = mp.mpf(str(cov[0][1]))
                    c2 = mp.mpf(str(cov[1][1]))
                    det = a * c2 - b * b
                    if det <= 0:
                        raise ValueError("covariance matrix must be positive definite")
                    inv00 = c2 / det
                    inv01 = -b / det
                    inv11 = a / det
                    chi2_cov = rvec[0] * (inv00 * rvec[0] + inv01 * rvec[1]) + rvec[1] * (inv01 * rvec[0] + inv11 * rvec[1])
                else:
                    raise ValueError("covariance gate currently supports only 1x1 or 2x2 matrices")

                p_cov = _chi2_sf(chi2=chi2_cov, dof=n)
                cov_gate = {
                    "cov_file": str(cov_path),
                    "labels": [str(x) for x in labels],
                    "chi2": chi2_cov,
                    "dof": n,
                    "p_value": p_cov,
                }

                # Severity: treat p<0.05 as WARN and p<1e-5 as FAIL (mirrors hard-sector χ² checks).
                if not mp.isfinite(p_cov):
                    cov_gate_check = mk_check_warn("alpha_0_with_covariance_gate", f"chi2={chi2_cov}, dof={n}, p={p_cov} (non-finite)")
                elif p_cov < mp.mpf("1e-5"):
                    cov_gate_check = mk_check_fail("alpha_0_with_covariance_gate", f"chi2={chi2_cov}, dof={n}, p={p_cov} (<1e-5)")
                elif p_cov < mp.mpf("0.05"):
                    cov_gate_check = mk_check_warn("alpha_0_with_covariance_gate", f"chi2={chi2_cov}, dof={n}, p={p_cov} (<0.05)")
                else:
                    cov_gate_check = mk_check_pass("alpha_0_with_covariance_gate", f"chi2={chi2_cov}, dof={n}, p={p_cov}")
            except Exception as e:
                cov_gate = {"cov_file": str(cov_path), "error": str(e)}
                cov_gate_check = mk_check_warn("alpha_0_with_covariance_gate", f"covariance gate not computed: {e}")

        # Gaussian terms (strict vs engineering)
        terms: list[dict[str, Any]] = []

        def add_gaussian_term(name: str) -> None:
            cfg = obs.get(name)
            if not isinstance(cfg, dict):
                return
            if cfg.get("enabled", True) is False:
                return
            if "mean" not in cfg or "sigma" not in cfg:
                return

            pred = predictions[name]
            mean = mp.mpf(str(cfg["mean"]))
            sigma_exp = mp.mpf(str(cfg["sigma"]))
            sigma_theory = mp.mpf(str(cfg.get("sigma_theory", 0)))
            sigma_total = mp.sqrt(sigma_exp**2 + sigma_theory**2)
            sigma_floor = mp.mpf(str(cfg.get("sigma_floor", 0)))
            sigma_eff = max(sigma_total, sigma_floor)

            z_exp = (pred - mean) / sigma_exp
            z_total = (pred - mean) / sigma_total if sigma_total != 0 else mp.mpf(0)
            z_eff = (pred - mean) / sigma_eff if sigma_eff != 0 else mp.mpf(0)

            ppm = None
            if mean != 0:
                ppm = (pred - mean) / mean * mp.mpf(1_000_000)

            terms.append(
                {
                    "name": name,
                    "pred": pred,
                    "mean": mean,
                    "sigma_exp": sigma_exp,
                    "sigma_theory": sigma_theory,
                    "sigma_total": sigma_total,
                    "sigma_floor": sigma_floor,
                    "sigma_eff": sigma_eff,
                    "z_exp": z_exp,
                    "z_total": z_total,
                    "z_eff": z_eff,
                    "chi2_total": z_total**2,
                    "chi2_eff": z_eff**2,
                    "ppm": ppm,
                    "source": cfg.get("source"),
                }
            )

        for key in ("alpha_bar5_inv_MZ", "alpha_inv_0", "beta_deg", "cabibbo_lambda", "n_s", "A_s"):
            add_gaussian_term(key)

        # One-sided upper bound (proxy) for r
        r_cfg = obs.get("r_upper_95", None)
        r_upper = mp.mpf(0)
        chi2_r_proxy = mp.mpf(0)
        r_enabled = False
        if isinstance(r_cfg, dict) and r_cfg.get("enabled", True) is not False and "upper" in r_cfg:
            r_enabled = True
            r_upper = mp.mpf(str(r_cfg["upper"]))
            if r_pred > r_upper:
                chi2_r_proxy = ((r_pred - r_upper) / r_upper) ** 2

        chi2_total_strict = mp.fsum([t["chi2_total"] for t in terms]) + (chi2_r_proxy if r_enabled else mp.mpf(0))
        chi2_total_engineering = mp.fsum([t["chi2_eff"] for t in terms]) + (chi2_r_proxy if r_enabled else mp.mpf(0))

        # Excluding alpha: useful diagnostic to avoid single-observable dominance
        alpha_keys = {"alpha_inv_0", "alpha_bar5_inv_MZ"}
        terms_no_alpha = [t for t in terms if t["name"] not in alpha_keys]
        chi2_total_strict_no_alpha = mp.fsum([t["chi2_total"] for t in terms_no_alpha]) + (chi2_r_proxy if r_enabled else mp.mpf(0))
        chi2_total_engineering_no_alpha = mp.fsum([t["chi2_eff"] for t in terms_no_alpha]) + (chi2_r_proxy if r_enabled else mp.mpf(0))

        # Degrees of freedom and p-values (dashboard; no covariance).
        dof_core = int(len(terms) + (1 if r_enabled else 0))
        dof_core_no_alpha = int(len(terms_no_alpha) + (1 if r_enabled else 0))
        p_core_strict = _chi2_sf(chi2=chi2_total_strict, dof=dof_core)
        p_core_engineering = _chi2_sf(chi2=chi2_total_engineering, dof=dof_core)
        p_core_strict_no_alpha = _chi2_sf(chi2=chi2_total_strict_no_alpha, dof=dof_core_no_alpha)
        p_core_engineering_no_alpha = _chi2_sf(chi2=chi2_total_engineering_no_alpha, dof=dof_core_no_alpha)

        # Physics-mode extensions: import hard-sector χ² from other modules (CKM/PMNS).
        extra_terms: list[dict[str, Any]] = []
        chi2_extra = mp.mpf(0)
        dof_extra = 0
        if mode == "physics":
            ckm_payload = _read_results_json(out_dir=Path(config.output_dir), module_id="ckm_full_pipeline")
            if isinstance(ckm_payload, dict):
                ckm_res = ckm_payload.get("results", {}) if isinstance(ckm_payload.get("results", {}), dict) else {}
                chi2_ckm = None
                try:
                    chi2_ckm = mp.mpf(str(ckm_res.get("rg_upward", {}).get("chi2_refscale")))
                except Exception:
                    chi2_ckm = None
                if chi2_ckm is not None and mp.isfinite(chi2_ckm):
                    contrib = ckm_res.get("contributions_best", None)
                    dof_ckm = int(len(contrib)) if isinstance(contrib, list) and len(contrib) > 0 else 9
                    p_ckm = _chi2_sf(chi2=chi2_ckm, dof=dof_ckm)
                    extra_terms.append(
                        {
                            "name": "ckm_full_pipeline_chi2",
                            "chi2_total": chi2_ckm,
                            "chi2_eff": chi2_ckm,
                            "dof": dof_ckm,
                            "p_value": p_ckm,
                            "source": str(Path(config.output_dir) / "ckm_full_pipeline" / "results.json"),
                        }
                    )
                    chi2_extra += chi2_ckm
                    dof_extra += dof_ckm

            pmns_payload = _read_results_json(out_dir=Path(config.output_dir), module_id="pmns_full_pipeline")
            if isinstance(pmns_payload, dict):
                pmns_res = pmns_payload.get("results", {}) if isinstance(pmns_payload.get("results", {}), dict) else {}
                chi2_pmns = None
                try:
                    chi2_pmns = mp.mpf(str(pmns_res.get("pmns_mt", {}).get("best_convention", {}).get("chi2")))
                except Exception:
                    chi2_pmns = None
                if chi2_pmns is not None and mp.isfinite(chi2_pmns):
                    contrib = pmns_res.get("pmns_mt", {}).get("best_convention", {}).get("contributions", None)
                    dof_pmns = int(len(contrib)) if isinstance(contrib, list) and len(contrib) > 0 else 4
                    p_pmns = _chi2_sf(chi2=chi2_pmns, dof=dof_pmns)
                    extra_terms.append(
                        {
                            "name": "pmns_full_pipeline_chi2",
                            "chi2_total": chi2_pmns,
                            "chi2_eff": chi2_pmns,
                            "dof": dof_pmns,
                            "p_value": p_pmns,
                            "source": str(Path(config.output_dir) / "pmns_full_pipeline" / "results.json"),
                        }
                    )
                    chi2_extra += chi2_pmns
                    dof_extra += dof_pmns

        dof_physics = int(dof_core + dof_extra)
        chi2_physics_strict = chi2_total_strict + chi2_extra
        chi2_physics_engineering = chi2_total_engineering + chi2_extra
        p_physics_strict = _chi2_sf(chi2=chi2_physics_strict, dof=dof_physics) if dof_physics > 0 else mp.mpf("nan")
        p_physics_engineering = _chi2_sf(chi2=chi2_physics_engineering, dof=dof_physics) if dof_physics > 0 else mp.mpf("nan")

        watch_strict = sorted((terms + extra_terms), key=lambda t: t["chi2_total"], reverse=True)
        watch_engineering = sorted((terms + extra_terms), key=lambda t: t["chi2_eff"], reverse=True)

        checks: list[Check] = []
        checks.append(
            (mk_check_pass if (isinstance(ref, dict) and isinstance(obs, dict) and len(obs) > 0) else mk_check_fail)(
                "reference_loaded",
                f"loaded {len(obs)} reference observables from {ref_path}",
            )
        )
        checks.append(
            (mk_check_pass if (len(terms) >= 3) else mk_check_fail)(
                "computed_terms",
                f"computed {len(terms)} Gaussian terms (+ r-bound proxy={r_enabled})",
            )
        )
        if cov_gate_check is not None:
            checks.append(cov_gate_check)

        # Physics-facing severity policy (as requested): |z|>5 => FAIL, 2<|z|<=5 => WARN.
        z_warn = mp.mpf(2)
        z_fail = mp.mpf(5)

        def _check_z(*, name: str, z: mp.mpf, extra: str = "") -> Check:
            if not mp.isfinite(z):
                return mk_check_fail(f"{name}_finite", f"z is not finite ({z}) {extra}".strip())
            az = abs(z)
            if az > z_fail:
                return mk_check_fail(f"{name}_within_5sigma", f"|z|={az} > 5 {extra}".strip())
            if az > z_warn:
                return mk_check_warn(f"{name}_within_2sigma", f"|z|={az} > 2 {extra}".strip())
            return mk_check_pass(f"{name}_within_2sigma", f"|z|={az} ≤ 2 {extra}".strip())

        # Per-observable z checks (keep compact: only for enabled terms).
        for t in terms:
            nm = str(t.get("name", ""))
            z = mp.mpf(t.get("z_exp", 0))
            checks.append(_check_z(name=nm, z=z, extra=f"(pred={t.get('pred')}, mean={t.get('mean')}, sigma_exp={t.get('sigma_exp')})"))

        # Hard-sector χ² checks (physics mode): gate by p-value.
        if mode == "physics":
            for t in extra_terms:
                nm = str(t.get("name", ""))
                dof = int(t.get("dof", 0) or 0)
                chi2 = mp.mpf(t.get("chi2_total", 0))
                pval = t.get("p_value", None)
                if pval is None:
                    checks.append(mk_check_warn(f"{nm}_pvalue_missing", f"chi2={chi2}, dof={dof}"))
                    continue
                p = mp.mpf(pval)
                if not mp.isfinite(p):
                    checks.append(mk_check_warn(f"{nm}_pvalue_nonfinite", f"chi2={chi2}, dof={dof}, p={p}"))
                    continue
                if p < mp.mpf("1e-5"):
                    checks.append(mk_check_fail(f"{nm}_pvalue_too_small", f"chi2={chi2}, dof={dof}, p={p}"))
                elif p < mp.mpf("0.05"):
                    checks.append(mk_check_warn(f"{nm}_pvalue_small", f"chi2={chi2}, dof={dof}, p={p}"))
                else:
                    checks.append(mk_check_pass(f"{nm}_pvalue_ok", f"chi2={chi2}, dof={dof}, p={p}"))

        # Total score summary checks (dashboard-style).
        checks.append(mk_check_info("core_score_p_value", f"chi2={chi2_total_strict}, dof={dof_core}, p={p_core_strict} (mode={mode})"))
        if mode == "physics":
            checks.append(mk_check_info("physics_score_p_value", f"chi2={chi2_physics_strict}, dof={dof_physics}, p={p_physics_strict}"))

        lines: list[str] = []
        lines += [
            "Global consistency test (χ²; reference-table driven; strict + engineering views)",
            "",
            f"reference file: {ref_path}",
            f"mode: {mode}",
            "",
            "TFPT predictions used:",
            f"- alpha_inv_0 (self-consistent, single_defect) = {alpha_inv_0_pred_single}",
            f"- alpha_inv_0 (self-consistent, two_defect; used for scoring) = {alpha_inv_0_pred_two_defect}",
            f"- delta2 model (two_defect): {d2.model_id} => delta2={d2.delta2} (delta2/delta_top^2={d2.delta2_over_delta_top2})",
            f"- alpha_bar5_inv_MZ (mapped from alpha(0) via alpha_running_pdg.json; MSbar-at-MZ comparison) = {alpha_bar5_inv_MZ_pred}",
            f"- beta_deg = {beta_deg_pred}",
            f"- cabibbo_lambda = {cabibbo_pred}",
            f"- n_s(N=56) = {ns_pred}",
            f"- r(N=56) = {r_pred}",
            f"- A_s(N=56) = {As_pred}",
            "",
            "Gaussian terms (strict uses sigma_total; engineering uses sigma_eff=max(sigma_total, sigma_floor)):",
        ]
        for t in terms:
            ppm_txt = f", ppm={t['ppm']}" if t["ppm"] is not None else ""
            lines.append(
                f"- {t['name']}: pred={t['pred']} mean={t['mean']} "
                f"sigma_exp={t['sigma_exp']} sigma_theory={t['sigma_theory']} sigma_total={t['sigma_total']} "
                f"sigma_floor={t['sigma_floor']} sigma_eff={t['sigma_eff']} "
                f"z_exp={t['z_exp']} z_total={t['z_total']} z_eff={t['z_eff']} "
                f"chi2_total={t['chi2_total']} chi2_eff={t['chi2_eff']}{ppm_txt}"
            )

        if r_enabled:
            lines.append("")
            lines.append(f"r bound proxy: r_pred={r_pred} <= r_upper_95={r_upper} => chi2_r_proxy={chi2_r_proxy}")

        lines += [
            "",
            f"CORE TOTAL chi2_strict = {chi2_total_strict}  (dof={dof_core}, p={p_core_strict})",
            f"CORE TOTAL chi2_engineering = {chi2_total_engineering}  (dof={dof_core}, p={p_core_engineering})",
            "",
            f"CORE TOTAL chi2_strict (excluding alpha) = {chi2_total_strict_no_alpha}  (dof={dof_core_no_alpha}, p={p_core_strict_no_alpha})",
            f"CORE TOTAL chi2_engineering (excluding alpha) = {chi2_total_engineering_no_alpha}  (dof={dof_core_no_alpha}, p={p_core_engineering_no_alpha})",
            "",
        ]

        if mode == "physics" and extra_terms:
            lines += [
                "Physics-mode extra terms (imported χ² from other modules):",
            ]
            for t in extra_terms:
                lines.append(f"- {t.get('name')}: chi2={t.get('chi2_total')} dof={t.get('dof')} p={t.get('p_value')} (src={t.get('source')})")
            lines += [
                "",
                f"PHYSICS TOTAL chi2_strict = {chi2_physics_strict}  (dof={dof_physics}, p={p_physics_strict})",
                f"PHYSICS TOTAL chi2_engineering = {chi2_physics_engineering}  (dof={dof_physics}, p={p_physics_engineering})",
                "",
            ]

        lines += [
            "Watchlist (largest contributors; strict):",
        ]
        for t in watch_strict[:5]:
            if "z_total" in t:
                lines.append(f"- {t['name']}: chi2_total={t['chi2_total']} (z_total={t['z_total']})")
            else:
                lines.append(f"- {t['name']}: chi2_total={t['chi2_total']} (p={t.get('p_value')}, dof={t.get('dof')})")
        lines += [
            "",
            "Watchlist (largest contributors; engineering):",
        ]
        for t in watch_engineering[:5]:
            if "z_eff" in t:
                lines.append(f"- {t['name']}: chi2_eff={t['chi2_eff']} (z_eff={t['z_eff']})")
            else:
                lines.append(f"- {t['name']}: chi2_eff={t['chi2_eff']} (p={t.get('p_value')}, dof={t.get('dof')})")

        lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module is only as good as the reference table. `global_reference.json` is marked vetted=true and includes per-observable sources/definitions; treat this as a dashboard unless you supply a full likelihood (incl. covariances).",
            "- If `global_reference_cov.json` is present, an auxiliary covariance-based α-sector gate is computed (check_id=alpha_0_with_covariance_gate).",
            "- The r bound handling is a simple one-sided proxy; implement the exact published likelihood if you need rigor.",
            "- If alpha dominates the strict χ² (because σ_exp is tiny), the report also provides an 'excluding alpha' total and an engineering total (sigma_floor).",
            "",
        ]

        warnings: list[str] = []
        if not bool(ref.get("vetted", False)):
            warnings.append("Reference values in global_reference.json are illustrative placeholders unless you update them.")
        if mode == "physics":
            names = {str(t.get("name", "")) for t in extra_terms}
            if "ckm_full_pipeline_chi2" not in names:
                warnings.append("physics_mode_missing_ckm_full_pipeline_output: run ckm_full_pipeline before scorecard (registry order should do this).")
            if "pmns_full_pipeline_chi2" not in names:
                warnings.append("physics_mode_missing_pmns_full_pipeline_output: run pmns_full_pipeline before scorecard (registry order should do this).")

        plot: dict[str, str | None] = {
            "chi2_contributions_png": None,
            "chi2_contributions_no_alpha_png": None,
            "global_radar_png": None,
        }
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_chi2_contributions(
                out_dir=out_dir, terms=(terms + extra_terms), r_enabled=r_enabled, chi2_r_proxy=chi2_r_proxy
            )
            plot_radar, plot_radar_warnings = _plot_global_radar(out_dir=out_dir, terms=terms, extra_terms=extra_terms)
            plot.update(plot_radar)
            warnings.extend(plot_radar_warnings)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "reference_file": str(ref_path),
                "mode": mode,
                "predictions": {
                    "alpha_inv_0": alpha_inv_0_pred,
                    "alpha_inv_0_single_defect": alpha_inv_0_pred_single,
                    "alpha_inv_0_two_defect": alpha_inv_0_pred_two_defect,
                    "delta2_model_id": d2.model_id,
                    "delta2": d2.delta2,
                    "delta2_over_delta_top2": d2.delta2_over_delta_top2,
                    "beta_deg": beta_deg_pred,
                    "cabibbo_lambda": cabibbo_pred,
                    "n_s": ns_pred,
                    "r": r_pred,
                    "A_s": As_pred,
                },
                "terms": terms,
                "extra_terms": extra_terms,
                "r_bound_proxy": {"enabled": r_enabled, "r_upper_95": r_upper if r_enabled else None, "chi2_r_proxy": chi2_r_proxy},
                "covariance_gate": cov_gate,
                "totals": {
                    "chi2_total_strict": chi2_total_strict,
                    "chi2_total_engineering": chi2_total_engineering,
                    "chi2_total_strict_excluding_alpha": chi2_total_strict_no_alpha,
                    "chi2_total_engineering_excluding_alpha": chi2_total_engineering_no_alpha,
                    "dof_core": dof_core,
                    "p_core_strict": p_core_strict,
                    "p_core_engineering": p_core_engineering,
                    "dof_core_excluding_alpha": dof_core_no_alpha,
                    "p_core_strict_excluding_alpha": p_core_strict_no_alpha,
                    "p_core_engineering_excluding_alpha": p_core_engineering_no_alpha,
                    "chi2_physics_strict": chi2_physics_strict if mode == "physics" else None,
                    "chi2_physics_engineering": chi2_physics_engineering if mode == "physics" else None,
                    "dof_physics": dof_physics if mode == "physics" else None,
                    "p_physics_strict": p_physics_strict if mode == "physics" else None,
                    "p_physics_engineering": p_physics_engineering if mode == "physics" else None,
                },
                "watchlist_strict": watch_strict,
                "watchlist_engineering": watch_engineering,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

