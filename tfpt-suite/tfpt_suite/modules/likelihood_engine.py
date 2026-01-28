from __future__ import annotations

import json
import math
import numpy as np
from mpmath import mp
from pathlib import Path
from typing import Any

from tfpt_suite.likelihood_engine import (
    evaluate_likelihood_spec,
    interpolate_grid_chi2,
    load_grid_table,
    load_likelihood_spec,
    predictions_from_global_consistency_results,
)
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass, mk_check_warn

NUFIT_DEFAULT_PARAM_COLUMNS = [
    "sin2_theta12",
    "sin2_theta13",
    "sin2_theta23",
    "delta_cp_deg",
    "dm21_sq_eV2",
    "dm3l_sq_eV2",
]
NUFIT_DEFAULT_CHI2_COLUMN = "chi2"
NUFIT_DEFAULT_INTERP_METHOD = "inverse_distance"
NUFIT_DEFAULT_K_NEAREST = 8
NUFIT_SCALE_SIN2 = 0.02
NUFIT_SCALE_DELTA_CP_DEG = 30.0
NUFIT_SCALE_DM21_SQ = 1.0e-5
NUFIT_SCALE_DM3L_SQ = 1.0e-3
NUFIT_DEFAULT_SCALES = {
    "sin2_theta12": NUFIT_SCALE_SIN2,
    "sin2_theta13": NUFIT_SCALE_SIN2,
    "sin2_theta23": NUFIT_SCALE_SIN2,
    "delta_cp_deg": NUFIT_SCALE_DELTA_CP_DEG,
    "dm21_sq_eV2": NUFIT_SCALE_DM21_SQ,
    "dm3l_sq_eV2": NUFIT_SCALE_DM3L_SQ,
    "dm31_sq_eV2": NUFIT_SCALE_DM3L_SQ,
    "dm32_sq_eV2": NUFIT_SCALE_DM3L_SQ,
}
UNIFIED_SCORE_ALPHA_KEYS = {"alpha_inv_0", "alpha_bar5_inv_MZ"}
UNIFIED_SCORE_CKM_SCALE = "refscale"
UNIFIED_SCORE_PMNS_SCALE = "mt"
UNIFIED_SCORE_DM_DATASET = "Omega_dm_h2_planck2018"


def _chi2_sf(*, chi2: float, dof: int) -> float:
    if not (math.isfinite(chi2) and dof > 0):
        return float("nan")
    try:
        k = mp.mpf(dof) / 2
        x = mp.mpf(chi2) / 2
        return float(mp.gammainc(k, x, mp.inf) / mp.gamma(k))
    except Exception:
        return float("nan")


def _as_float(value: object) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _sin2_from_deg(theta_deg: float) -> float:
    if not math.isfinite(theta_deg):
        return float("nan")
    return float(math.sin(math.radians(theta_deg)) ** 2)


def _pmns_predictions_from_results(payload: dict[str, Any], *, scale_key: str) -> tuple[dict[str, float], dict[str, Any]]:
    res = payload.get("results", {}) if isinstance(payload.get("results", {}), dict) else {}
    pmns_key = "pmns_mt" if scale_key == "mt" else "pmns_mu_uv"
    pmns = res.get(pmns_key, {}) if isinstance(res.get(pmns_key, {}), dict) else {}
    angles = pmns.get("angles_deg", {}) if isinstance(pmns.get("angles_deg", {}), dict) else {}
    best_conv = pmns.get("best_convention", {}) if isinstance(pmns.get("best_convention", {}), dict) else {}
    ordering = str(best_conv.get("ordering", "")).upper() or None

    theta12_deg = _as_float(angles.get("theta12_deg"))
    theta13_deg = _as_float(angles.get("theta13_deg"))
    theta23_deg = _as_float(angles.get("theta23_deg"))
    delta_cp_deg = _as_float(angles.get("delta_cp_deg"))
    masses = res.get("neutrino_masses_proxy_eV", {}) if isinstance(res.get("neutrino_masses_proxy_eV", {}), dict) else {}
    masses_vec = masses.get(scale_key, None)
    dm21_sq = float("nan")
    dm3l_sq = float("nan")
    dm31_sq = float("nan")
    dm32_sq = float("nan")
    if isinstance(masses_vec, list) and len(masses_vec) >= 3:
        m1, m2, m3 = (_as_float(masses_vec[0]), _as_float(masses_vec[1]), _as_float(masses_vec[2]))
        if all(math.isfinite(x) for x in (m1, m2, m3)):
            dm21_sq = float(m2 * m2 - m1 * m1)
            dm31_sq = float(m3 * m3 - m1 * m1)
            dm32_sq = float(m3 * m3 - m2 * m2)
            if ordering == "IO":
                dm3l_sq = dm32_sq
            else:
                dm3l_sq = dm31_sq

    preds = {
        "theta12_deg": theta12_deg,
        "theta13_deg": theta13_deg,
        "theta23_deg": theta23_deg,
        "delta_cp_deg": delta_cp_deg,
        "sin2_theta12": _sin2_from_deg(theta12_deg),
        "sin2_theta13": _sin2_from_deg(theta13_deg),
        "sin2_theta23": _sin2_from_deg(theta23_deg),
        "dm21_sq_eV2": dm21_sq,
        "dm3l_sq_eV2": dm3l_sq,
        "dm31_sq_eV2": dm31_sq,
        "dm32_sq_eV2": dm32_sq,
    }
    meta = {"ordering": ordering, "pmns_key": pmns_key, "masses_key": scale_key}
    return preds, meta


def _scales_for_columns(param_columns: list[str], scale_overrides: dict[str, Any] | None) -> list[float]:
    scales: list[float] = []
    overrides = scale_overrides or {}
    for col in param_columns:
        scale = _as_float(overrides.get(col, NUFIT_DEFAULT_SCALES.get(col, 1.0)))
        if not math.isfinite(scale) or scale <= 0:
            scale = 1.0
        scales.append(scale)
    return scales


class LikelihoodEngineModule(TfptModule):
    module_id = "likelihood_engine"
    title = "Likelihood engine (covariance datasets + nuisance-policy contract; plugin-ready)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "dataset spec: tfpt_suite/data/likelihood_datasets_v1.json",
                "predictions provider: global_consistency_test output (preferred) or direct prediction providers (future)",
            ],
            outputs=[
                "per-dataset chi2/logL summary (explicit covariance/bound handling)",
                "total log-likelihood (sum of enabled datasets)",
                "optional NuFIT PMNS grid logL (non-Gaussian, when enabled)",
            ],
            formulas=[
                "multivariate Gaussian: χ² = rᵀ C⁻¹ r, logL = -1/2 (χ² + log det C + n log 2π)",
                "one-sided upper bound (proxy): χ²=0 if x≤x_max else ((x-x_max)/x_max)²",
                "grid χ²: χ²(θ) interpolated from NuFIT grid; logL = -1/2 χ²",
            ],
            validation=[
                "covariance_enabled_for_reference_tables: PASS if a covariance dataset is evaluated from the unified spec.",
                "nuisance_handling_policy_explicit: PASS if nuisance_policy is present in the spec.",
                "nufit_pmns_grid_plugin_active: INFO/WARN unless a NuFIT grid is enabled and evaluated.",
            ],
            determinism="Deterministic given the dataset spec + prediction providers.",
            question="Can TFPT evaluate explicit covariances/bounds under a declared nuisance policy (upgrade path to Planck/NuFIT plugins)?",
            objective=[
                "Provide a single unified dataset schema that can later host Planck (plugin) and NuFIT (grid chi2) without ad-hoc logic in each physics module.",
            ],
            gaps=[
                "External plugins (Planck likelihood, NuFIT chi2 grids) are not enabled by default; this module defines the contract and evaluates shipped datasets.",
            ],
        )

    def run(self, config) -> ModuleResult:
        data_dir = Path(__file__).resolve().parent.parent / "data"
        spec_path = data_dir / "likelihood_datasets_v1.json"
        spec = load_likelihood_spec(spec_path)

        # Preferred prediction source: global_consistency_test output in the same output_dir.
        gc_path = Path(config.output_dir) / "global_consistency_test" / "results.json"
        preds: dict[str, float] = {}
        preds_source = "none"
        if gc_path.is_file():
            try:
                payload = json.loads(gc_path.read_text(encoding="utf-8"))
                preds = predictions_from_global_consistency_results(payload)
                preds_source = str(gc_path)
            except Exception:
                preds = {}

        # Evaluate datasets (deterministic, no sampling).
        results, logL_datasets = evaluate_likelihood_spec(spec=spec, predictions=preds, include_norm_constant=True)

        checks: list[Check] = []
        checks.append(mk_check_pass("likelihood_engine_runs", f"evaluated {len(results)} dataset(s); logL_datasets={logL_datasets}"))

        nuisance = spec.get("nuisance_policy", {})
        checks.append(
            mk_check_pass("nuisance_handling_policy_explicit", f"policy={nuisance}")
            if isinstance(nuisance, dict) and bool(nuisance.get("kind", ""))
            else mk_check_warn("nuisance_handling_policy_explicit", "missing nuisance_policy.kind in likelihood spec")
        )

        has_cov = any(r.kind == "multivariate_gaussian" and len(r.labels) >= 2 for r in results)
        checks.append(
            mk_check_pass("covariance_enabled_for_reference_tables", "at least one multivariate covariance dataset evaluated")
            if has_cov
            else mk_check_warn("covariance_enabled_for_reference_tables", "no multivariate covariance dataset evaluated")
        )

        # Plugin summary (not executed unless enabled and dependencies exist)
        plugins = spec.get("plugins", [])
        if isinstance(plugins, list) and plugins:
            enabled = [p for p in plugins if isinstance(p, dict) and bool(p.get("enabled", False))]
            checks.append(mk_check_info("plugins_declared", f"plugins={len(plugins)}, enabled={len(enabled)}"))

        # Optional plugin contribution: Planck logL components (provided by boltzmann_transfer when enabled).
        logL_plugins = 0.0
        planck_from_suite: dict[str, Any] | None = None
        bt_path = Path(config.output_dir) / "boltzmann_transfer" / "results.json"
        if bt_path.is_file():
            try:
                payload = json.loads(bt_path.read_text(encoding="utf-8"))
                bt_res = payload.get("results", {}) if isinstance(payload, dict) else {}
                planck_from_suite = {}

                def _include_planck_component(label: str, component: dict[str, Any] | None) -> None:
                    nonlocal logL_plugins
                    if not isinstance(component, dict):
                        return
                    planck_from_suite[label] = component
                    if bool(component.get("enabled", False)) and component.get("logp", None) is not None:
                        logp = float(component["logp"])
                        if math.isfinite(logp):
                            logL_plugins += float(logp)
                            checks.append(mk_check_pass(f"planck_{label}_included", f"logp={logp} (from {bt_path})"))
                        else:
                            checks.append(mk_check_info(f"planck_{label}_included", f"present but non-finite logp (from {bt_path})"))
                    else:
                        checks.append(mk_check_info(f"planck_{label}_included", f"not enabled / no logp (from {bt_path})"))

                _include_planck_component("pliklite", bt_res.get("planck_pliklite", None))
                _include_planck_component("lowl", bt_res.get("planck_lowl", None))
                _include_planck_component("lensing", bt_res.get("planck_lensing", None))

                planck_combined = bt_res.get("planck_combined", None)
                if isinstance(planck_combined, dict):
                    planck_from_suite["combined"] = planck_combined
            except Exception as e:
                checks.append(mk_check_info("planck_likelihood_included", f"failed to read {bt_path}: {e}"))
        else:
            checks.append(mk_check_info("planck_likelihood_included", f"no boltzmann_transfer output found at {bt_path}"))

        # Optional plugin contribution: NuFIT PMNS grid (non-Gaussian, when enabled).
        nufit_result: dict[str, Any] | None = None
        nufit_plugin = next(
            (p for p in plugins if isinstance(p, dict) and p.get("plugin_id") == "nufit_pmns_grid"),
            None,
        )
        if isinstance(nufit_plugin, dict):
            enabled = bool(nufit_plugin.get("enabled", False))
            if not enabled:
                checks.append(mk_check_info("nufit_pmns_grid_plugin_active", "disabled in likelihood spec"))
            else:
                grid_file = str(nufit_plugin.get("grid_file", "")).strip()
                param_columns_raw = nufit_plugin.get("param_columns", None)
                if not isinstance(param_columns_raw, list) or not param_columns_raw:
                    param_columns = list(NUFIT_DEFAULT_PARAM_COLUMNS)
                else:
                    param_columns = [str(x) for x in param_columns_raw]
                chi2_column = str(nufit_plugin.get("chi2_column", NUFIT_DEFAULT_CHI2_COLUMN)).strip() or NUFIT_DEFAULT_CHI2_COLUMN
                interp_cfg = nufit_plugin.get("interpolation", {}) if isinstance(nufit_plugin.get("interpolation", {}), dict) else {}
                method = str(interp_cfg.get("method", NUFIT_DEFAULT_INTERP_METHOD))
                k_nearest = int(interp_cfg.get("k_nearest", NUFIT_DEFAULT_K_NEAREST))
                scale_key = str(nufit_plugin.get("pmns_scale", "mt")).strip() or "mt"
                ordering_policy = str(nufit_plugin.get("ordering_policy", "from_pmns_full_pipeline")).strip()
                chi2_kind = str(nufit_plugin.get("chi2_kind", "delta")).strip()
                normalize_delta = bool(nufit_plugin.get("normalize_delta_chi2", True))
                scale_overrides = nufit_plugin.get("distance_scales", {}) if isinstance(nufit_plugin.get("distance_scales", {}), dict) else {}

                if not grid_file:
                    checks.append(mk_check_warn("nufit_pmns_grid_plugin_active", "enabled but grid_file is missing"))
                else:
                    grid_path = Path(grid_file)
                    if not grid_path.is_absolute():
                        grid_path = data_dir / grid_path
                    if not grid_path.is_file():
                        checks.append(mk_check_warn("nufit_pmns_grid_plugin_active", f"grid file not found: {grid_path}"))
                    else:
                        pmns_path = Path(config.output_dir) / "pmns_full_pipeline" / "results.json"
                        if not pmns_path.is_file():
                            checks.append(mk_check_warn("nufit_pmns_grid_plugin_active", f"pmns_full_pipeline output missing: {pmns_path}"))
                        else:
                            try:
                                pmns_payload = json.loads(pmns_path.read_text(encoding="utf-8"))
                                pmns_preds, pmns_meta = _pmns_predictions_from_results(pmns_payload, scale_key=scale_key)
                                if ordering_policy.lower() == "force_no":
                                    pmns_meta["ordering"] = "NO"
                                elif ordering_policy.lower() == "force_io":
                                    pmns_meta["ordering"] = "IO"

                                target = [_as_float(pmns_preds.get(col, float("nan"))) for col in param_columns]
                                if not all(math.isfinite(x) for x in target):
                                    checks.append(
                                        mk_check_warn("nufit_pmns_grid_plugin_active", f"non-finite PMNS target vector for columns={param_columns}")
                                    )
                                else:
                                    points, chi2_grid, grid_meta = load_grid_table(
                                        grid_path,
                                        param_columns=param_columns,
                                        chi2_column=chi2_column,
                                    )
                                    scales = _scales_for_columns(param_columns, scale_overrides)
                                    chi2_val, interp_meta = interpolate_grid_chi2(
                                        points=points,
                                        chi2=chi2_grid,
                                        target=target,
                                        method=method,
                                        k_nearest=k_nearest,
                                        scales=scales,
                                    )
                                    chi2_min = float(np.nanmin(chi2_grid)) if chi2_grid.size else float("nan")
                                    if chi2_kind == "delta" and normalize_delta and math.isfinite(chi2_min):
                                        chi2_val = float(chi2_val - chi2_min)
                                    logp = float(-0.5 * chi2_val) if math.isfinite(chi2_val) else float("nan")
                                    if math.isfinite(logp):
                                        logL_plugins += logp
                                        checks.append(mk_check_pass("nufit_pmns_grid_plugin_active", f"logp={logp:.6g} (chi2={chi2_val:.6g})"))
                                    else:
                                        checks.append(mk_check_warn("nufit_pmns_grid_plugin_active", f"non-finite logp (chi2={chi2_val})"))

                                    nufit_result = {
                                        "enabled": True,
                                        "grid_file": str(grid_path),
                                        "param_columns": param_columns,
                                        "chi2_column": chi2_column,
                                        "chi2_kind": chi2_kind,
                                        "normalize_delta_chi2": normalize_delta,
                                        "chi2_grid_min": chi2_min,
                                        "pmns_predictions": pmns_preds,
                                        "pmns_meta": pmns_meta,
                                        "target_vector": target,
                                        "interpolation": {"method": method, "k_nearest": k_nearest, "scales": scales, **interp_meta},
                                        "grid_meta": grid_meta,
                                        "chi2": chi2_val,
                                        "logp": logp,
                                    }
                            except Exception as exc:
                                checks.append(mk_check_warn("nufit_pmns_grid_plugin_active", f"failed to evaluate grid: {exc}"))

        # Unified scorecard (α + flavor + CMB + DM).
        unified_components: list[dict[str, Any]] = []
        unified_missing: list[str] = []

        alpha_chi2 = 0.0
        alpha_dof = 0
        for r in results:
            if any(k in UNIFIED_SCORE_ALPHA_KEYS for k in r.labels):
                if math.isfinite(float(r.chi2)) and r.dof > 0:
                    alpha_chi2 += float(r.chi2)
                    alpha_dof += int(r.dof)
        if alpha_dof > 0:
            unified_components.append(
                {"sector": "alpha", "chi2": alpha_chi2, "dof": alpha_dof, "source": "likelihood_datasets_v1.json"}
            )
        else:
            unified_missing.append("alpha")

        # Flavor: CKM chi2 (refscale).
        ckm_path = Path(config.output_dir) / "ckm_full_pipeline" / "results.json"
        if ckm_path.is_file():
            try:
                payload = json.loads(ckm_path.read_text(encoding="utf-8"))
                res = payload.get("results", {}) if isinstance(payload.get("results", {}), dict) else {}
                rg = res.get("rg_upward", {}) if isinstance(res.get("rg_upward", {}), dict) else {}
                chi2_ref = float(rg.get("chi2_refscale", float("nan")))
                chi2_keys = rg.get("chi2_keys", [])
                dof = int(len(chi2_keys)) if isinstance(chi2_keys, list) else 0
                if math.isfinite(chi2_ref) and dof > 0:
                    unified_components.append({"sector": "flavor_ckm", "chi2": chi2_ref, "dof": dof, "source": str(ckm_path)})
                else:
                    unified_missing.append("flavor_ckm")
            except Exception:
                unified_missing.append("flavor_ckm")
        else:
            unified_missing.append("flavor_ckm")

        # Flavor: PMNS chi2 (prefer NuFIT grid if enabled).
        if nufit_result is not None and math.isfinite(float(nufit_result.get("chi2", float("nan")))):
            chi2_val = float(nufit_result["chi2"])
            dof = int(len(nufit_result.get("param_columns", [])))
            if dof > 0:
                unified_components.append(
                    {"sector": "flavor_pmns", "chi2": chi2_val, "dof": dof, "source": str(nufit_result.get("grid_file"))}
                )
            else:
                unified_missing.append("flavor_pmns")
        else:
            pmns_path = Path(config.output_dir) / "pmns_full_pipeline" / "results.json"
            if pmns_path.is_file():
                try:
                    payload = json.loads(pmns_path.read_text(encoding="utf-8"))
                    res = payload.get("results", {}) if isinstance(payload.get("results", {}), dict) else {}
                    pmns = res.get("pmns_mt", {}) if isinstance(res.get("pmns_mt", {}), dict) else {}
                    best = pmns.get("best_convention", {}) if isinstance(pmns.get("best_convention", {}), dict) else {}
                    chi2_val = float(best.get("chi2", float("nan")))
                    contribs = best.get("contributions", [])
                    dof = int(len(contribs)) if isinstance(contribs, list) else 0
                    if math.isfinite(chi2_val) and dof > 0:
                        unified_components.append({"sector": "flavor_pmns", "chi2": chi2_val, "dof": dof, "source": str(pmns_path)})
                    else:
                        unified_missing.append("flavor_pmns")
                except Exception:
                    unified_missing.append("flavor_pmns")
            else:
                unified_missing.append("flavor_pmns")

        # CMB: prefer Planck combined logp (from boltzmann_transfer), else r upper bound proxy.
        cmb_added = False
        if isinstance(planck_from_suite, dict):
            combined = planck_from_suite.get("combined", None)
            if isinstance(combined, dict) and math.isfinite(float(combined.get("chi2", float("nan")))):
                components = combined.get("components", {})
                dof_proxy = len(components) if isinstance(components, dict) and components else 1
                unified_components.append(
                    {"sector": "cmb", "chi2": float(combined["chi2"]), "dof": int(dof_proxy), "source": "planck_combined"}
                )
                cmb_added = True
        if not cmb_added:
            r_proxy = next((r for r in results if r.dataset_id == "r_upper_95_proxy"), None)
            if r_proxy is not None and math.isfinite(float(r_proxy.chi2)):
                unified_components.append(
                    {"sector": "cmb", "chi2": float(r_proxy.chi2), "dof": int(r_proxy.dof), "source": "r_upper_95_proxy"}
                )
            else:
                unified_missing.append("cmb")

        # DM: use axion_dm_pipeline output vs reference sigma from references.json.
        dm_ref_val = float("nan")
        dm_ref_sig = float("nan")
        ref_path = data_dir / "references.json"
        if ref_path.is_file():
            try:
                ref_payload = json.loads(ref_path.read_text(encoding="utf-8"))
                datasets = ref_payload.get("datasets", {}) if isinstance(ref_payload.get("datasets", {}), dict) else {}
                dm_ref = datasets.get(UNIFIED_SCORE_DM_DATASET, {})
                if isinstance(dm_ref, dict):
                    dm_ref_val = float(dm_ref.get("value", float("nan")))
                    dm_ref_sig = float(dm_ref.get("sigma", float("nan")))
            except Exception:
                dm_ref_val = float("nan")
                dm_ref_sig = float("nan")

        dm_path = Path(config.output_dir) / "axion_dm_pipeline" / "results.json"
        if dm_path.is_file() and math.isfinite(dm_ref_val) and math.isfinite(dm_ref_sig) and dm_ref_sig > 0:
            try:
                dm_payload = json.loads(dm_path.read_text(encoding="utf-8"))
                res = dm_payload.get("results", {}) if isinstance(dm_payload.get("results", {}), dict) else {}

                def _omega_a_from_results(results: dict[str, Any]) -> float | None:
                    for key in ("summary", "relic_density"):
                        block = results.get(key, {})
                        if isinstance(block, dict):
                            value = float(block.get("Omega_a_h2", float("nan")))
                            if math.isfinite(value):
                                return float(value)
                    return None

                omega_a = _omega_a_from_results(res)
                if omega_a is not None and math.isfinite(omega_a):
                    chi2_dm = float(((omega_a - dm_ref_val) / dm_ref_sig) ** 2)
                    unified_components.append(
                        {
                            "sector": "dm",
                            "chi2": chi2_dm,
                            "dof": 1,
                            "source": str(dm_path),
                            "omega_a_h2": omega_a,
                            "omega_dm_h2_ref": dm_ref_val,
                            "sigma": dm_ref_sig,
                        }
                    )
                else:
                    unified_missing.append("dm")
            except Exception:
                unified_missing.append("dm")
        else:
            unified_missing.append("dm")

        unified_chi2 = float(sum(c["chi2"] for c in unified_components if math.isfinite(float(c["chi2"]))))
        unified_dof = int(sum(int(c["dof"]) for c in unified_components if isinstance(c.get("dof", None), int)))
        unified_p = _chi2_sf(chi2=unified_chi2, dof=unified_dof)
        all_sectors = {"alpha", "flavor_ckm", "flavor_pmns", "cmb", "dm"}
        missing_unique = sorted({m for m in unified_missing if m in all_sectors})
        if missing_unique:
            checks.append(mk_check_warn("unified_score_all_sectors", f"missing sectors: {missing_unique}"))
        else:
            checks.append(mk_check_pass("unified_score_all_sectors", "all sectors contributing (alpha, flavor, cmb, dm)"))
        if math.isfinite(unified_p):
            checks.append(mk_check_info("unified_score_p_value_reported", f"p≈{unified_p:.6g} (chi2={unified_chi2:.6g}, dof={unified_dof})"))
        else:
            checks.append(mk_check_warn("unified_score_p_value_reported", f"non-finite p-value (chi2={unified_chi2}, dof={unified_dof})"))

        lines: list[str] = []
        lines += [
            "Likelihood engine (unified dataset schema)",
            f"spec: {spec_path}",
            f"predictions source: {preds_source}",
            f"nuisance_policy: {json.dumps(nuisance, ensure_ascii=False)}",
            "",
            "Predictions:",
            json.dumps(preds, indent=2, sort_keys=True),
            "",
            "Datasets:",
        ]
        for r in results:
            lines.append(f"- {r.dataset_id} ({r.kind}): chi2={r.chi2:.6g}, dof={r.dof}, logL={r.loglike:.6g}, labels={r.labels}")
            lines.append(f"  details: {r.details}")
        lines += [
            "",
            f"Total logL (datasets) = {logL_datasets}",
            f"Total logL (plugins) = {logL_plugins}",
            f"Total logL (datasets + plugins) = {logL_datasets + logL_plugins}",
            "",
            "Unified scorecard (χ² proxy):",
            *[
                f"- {c['sector']}: chi2={c['chi2']:.6g}, dof={c['dof']} (source={c.get('source')})"
                for c in unified_components
            ],
            f"- total chi2={unified_chi2:.6g}, dof={unified_dof}, p≈{unified_p:.6g}",
        ]
        if isinstance(nufit_plugin, dict):
            lines += [
                "",
                "NuFIT grid plugin:",
                f"- enabled = {bool(nufit_plugin.get('enabled', False))}",
                f"- grid_file = {nufit_plugin.get('grid_file')}",
                f"- param_columns = {nufit_plugin.get('param_columns', NUFIT_DEFAULT_PARAM_COLUMNS)}",
                f"- chi2_column = {nufit_plugin.get('chi2_column', NUFIT_DEFAULT_CHI2_COLUMN)}",
            ]
            if nufit_result is not None:
                lines.append(f"- chi2 = {nufit_result.get('chi2')}, logp = {nufit_result.get('logp')}")
        lines += [
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "spec_file": str(spec_path),
                "predictions_source": preds_source,
                "predictions": preds,
                "datasets": [
                    {"dataset_id": r.dataset_id, "kind": r.kind, "labels": r.labels, "chi2": r.chi2, "dof": r.dof, "loglike": r.loglike, "details": r.details}
                    for r in results
                ],
                "loglike_total_datasets": logL_datasets,
                "loglike_total_plugins": logL_plugins,
                "loglike_total": float(logL_datasets + logL_plugins),
                "nuisance_policy": nuisance,
                "plugins": plugins,
                "planck_pliklite_from_boltzmann_transfer": planck_from_suite,
                "nufit_pmns_grid": nufit_result,
                "unified_score": {
                    "components": unified_components,
                    "chi2_total": unified_chi2,
                    "dof_total": unified_dof,
                    "p_value": unified_p,
                    "missing_sectors": missing_unique,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

