from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass, mk_check_warn


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


def _chi2_sf(*, chi2: float, dof: int) -> float:
    try:
        from mpmath import mp

        k = int(dof)
        if k <= 0:
            return float("nan")
        x = mp.mpf(str(chi2)) / 2
        s = mp.mpf(k) / 2
        p = mp.gammainc(s, x, mp.inf) / mp.gamma(s)
        return float(p)
    except Exception:
        return float("nan")


class UncertaintyPropagatorModule(TfptModule):
    module_id = "uncertainty_propagator"
    title = "Uncertainty propagator (end-to-end MC summary aggregator; covariance scaffold)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "upstream outputs in out/: msbar_matching_map, ckm_full_pipeline, pmns_full_pipeline (if present)",
            ],
            outputs=[
                "summary of available uncertainty propagation artifacts (means/stds) across modules",
                "joint flavor objective distribution proxy (mean/std) when upstream MC artifacts are available",
            ],
            formulas=[
                "Assume independence: Var(sum) = Σ Var_i (placeholder; replace by covariance matrix when available).",
            ],
            validation=[
                "Reports which uncertainty pieces exist and which are missing (reviewer-proof).",
                "Emits a physics-mode gate check for the joint flavor χ² with uncertainty (proxy).",
            ],
            determinism="Deterministic given upstream artifacts.",
            question="Do we have end-to-end uncertainty propagation artifacts, and what is the joint flavor χ² once uncertainties are accounted for (proxy)?",
            objective=[
                "Centralize uncertainty bookkeeping (avoid scattering MC summaries across modules).",
                "Prepare the suite for a future covariance-based global fit.",
            ],
            gaps=[
                "This module is an aggregator; publication-grade requires explicit covariance support and consistent propagation through all thresholds/policies (incl. QED/EW below MZ).",
            ],
        )

    def run(self, config) -> ModuleResult:
        out_dir = Path(config.output_dir)
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        msbar = _read_results_json(out_dir=out_dir, module_id="msbar_matching_map")
        ckm = _read_results_json(out_dir=out_dir, module_id="ckm_full_pipeline")
        pmns = _read_results_json(out_dir=out_dir, module_id="pmns_full_pipeline")

        checks: list[Check] = []
        checks.append(mk_check_pass("uncertainty_propagator_runs", "module executed"))

        # Extract MC summaries if present.
        def _get(d: dict[str, Any] | None, *keys: str) -> Any:
            cur: Any = d
            for k in keys:
                if not isinstance(cur, dict):
                    return None
                cur = cur.get(k)
            return cur

        # nominal flavor chi2 (no uncertainty)
        chi2_ckm_nom = None
        chi2_pmns_nom = None
        try:
            chi2_ckm_nom = float(_get(ckm, "results", "rg_upward", "chi2_refscale"))
        except Exception:
            chi2_ckm_nom = None
        try:
            chi2_pmns_nom = float(_get(pmns, "results", "pmns_mt", "best_convention", "chi2"))
        except Exception:
            chi2_pmns_nom = None

        # MC summaries
        ckm_mc = _get(ckm, "results", "matching_mc")
        pmns_mc = _get(pmns, "results", "matching_mc")
        msbar_mc = _get(msbar, "results", "matching_mc")

        def _extract_mc_mean_std(mc: Any, mean_key: str, std_key: str) -> tuple[float | None, float | None]:
            if not isinstance(mc, dict):
                return None, None
            try:
                mu = float(mc.get(mean_key))
                sig = float(mc.get(std_key))
                if not (np.isfinite(mu) and np.isfinite(sig) and sig >= 0):
                    return None, None
                return mu, sig
            except Exception:
                return None, None

        ckm_mu, ckm_sig = _extract_mc_mean_std(ckm_mc, "chi2_refscale_mean", "chi2_refscale_std")
        pmns_mu, pmns_sig = _extract_mc_mean_std(pmns_mc, "chi2_uv_mean", "chi2_uv_std")

        # Joint flavor objective proxy (use what we have; keep explicit about scale mismatch for PMNS MC).
        joint_nom = None
        if chi2_ckm_nom is not None and chi2_pmns_nom is not None:
            joint_nom = float(chi2_ckm_nom + chi2_pmns_nom)
            checks.append(mk_check_pass("flavor_joint_objective_nominal_present", f"chi2_joint_nominal={joint_nom}"))
        else:
            checks.append(mk_check_warn("flavor_joint_objective_nominal_present", f"missing nominal chi2 (ckm={chi2_ckm_nom}, pmns={chi2_pmns_nom})"))

        joint_mc = None
        if ckm_mu is not None and ckm_sig is not None and pmns_mu is not None and pmns_sig is not None:
            joint_mu = float(ckm_mu + pmns_mu)
            joint_sig = float(np.sqrt(ckm_sig * ckm_sig + pmns_sig * pmns_sig))
            joint_mc = {"mean": joint_mu, "std": joint_sig, "note": "PMNS MC uses chi2_uv_*; CKM uses chi2_refscale_* (scale mismatch; treat as proxy)."}
            checks.append(mk_check_pass("flavor_joint_objective_with_unc_present", f"mean={joint_mu:.6g}, std={joint_sig:.3g}"))
        else:
            checks.append(mk_check_warn("flavor_joint_objective_with_unc_present", "missing MC summaries for CKM and/or PMNS"))

        # Physics-mode gate (proxy): compute p-value for joint nominal (dof proxy = 9+4).
        dof_proxy = 13
        if joint_nom is not None and np.isfinite(joint_nom):
            p = _chi2_sf(chi2=float(joint_nom), dof=dof_proxy)
            if np.isfinite(p):
                if mode == "physics" and p < 0.05:
                    checks.append(mk_check_fail("flavor_chi2_with_unc", f"FAIL: p≈{p:.3e} (<0.05) chi2≈{joint_nom:.6g} dof={dof_proxy} (nominal proxy)"))
                elif p < 0.05:
                    checks.append(mk_check_warn("flavor_chi2_with_unc", f"WARN: p≈{p:.3e} (<0.05) chi2≈{joint_nom:.6g} dof={dof_proxy} (nominal proxy)"))
                else:
                    checks.append(mk_check_pass("flavor_chi2_with_unc", f"PASS: p≈{p:.3e} chi2≈{joint_nom:.6g} dof={dof_proxy} (nominal proxy)"))
            else:
                checks.append(mk_check_warn("flavor_chi2_with_unc", "p-value non-finite"))
        else:
            checks.append(mk_check_warn("flavor_chi2_with_unc", "chi2_joint_nominal missing/non-finite"))

        lines: list[str] = []
        lines += [
            "Uncertainty propagator (aggregator)",
            "",
            f"mode={mode}",
            "",
            "Upstream artifacts present:",
            f"- msbar_matching_map: {'yes' if msbar is not None else 'no'}",
            f"- ckm_full_pipeline: {'yes' if ckm is not None else 'no'}",
            f"- pmns_full_pipeline: {'yes' if pmns is not None else 'no'}",
            "",
            "Flavor χ² (nominal):",
            f"- chi2_ckm_refscale = {chi2_ckm_nom}",
            f"- chi2_pmns_mt = {chi2_pmns_nom}",
            f"- chi2_joint_nominal = {joint_nom}",
            "",
            "Flavor χ² (MC summaries, if available):",
            f"- ckm chi2_refscale mean±std = {ckm_mu} ± {ckm_sig}",
            f"- pmns chi2_uv mean±std = {pmns_mu} ± {pmns_sig}",
            f"- joint_mc = {joint_mc}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "upstream": {
                    "msbar_matching_map_present": msbar is not None,
                    "ckm_full_pipeline_present": ckm is not None,
                    "pmns_full_pipeline_present": pmns is not None,
                },
                "nominal": {"chi2_ckm_refscale": chi2_ckm_nom, "chi2_pmns_mt": chi2_pmns_nom, "chi2_joint": joint_nom, "dof_proxy": dof_proxy},
                "mc_summaries": {
                    "msbar_matching_map": msbar_mc,
                    "ckm_full_pipeline": ckm_mc,
                    "pmns_full_pipeline": pmns_mc,
                },
                "joint_mc_proxy": joint_mc,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

