from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.constants import TfptConstants
from tfpt_suite.mobius_z3_yukawa_generator import _yukawa_ratios_from_mobius
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass, mk_check_warn

MASS_RATIO_TOL_REL_DEFAULT = 0.05


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


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            return payload
    except Exception:
        return None
    return None


def _mass_ratio_penalty_from_leptons(*, delta_used: float, tol_rel: float) -> tuple[float, dict[str, float]]:
    lep_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
    lep = _read_json(lep_path)
    if not isinstance(lep, dict):
        return float("nan"), {"error": float("nan")}
    try:
        me = float(lep["masses"]["electron"]["mean"])
        mm = float(lep["masses"]["muon"]["mean"])
        mt = float(lep["masses"]["tau"]["mean"])
        meas_tau_mu = float(mt / mm) if mm > 0 else float("nan")
        meas_mu_e = float(mm / me) if me > 0 else float("nan")
        ratios = _yukawa_ratios_from_mobius(delta_used=float(delta_used))
        pred_tau_mu = float(ratios.get("m_tau_over_m_mu", float("nan")))
        pred_mu_e = float(ratios.get("m_mu_over_m_e", float("nan")))
        rel_tau = float((pred_tau_mu / meas_tau_mu) - 1.0) if meas_tau_mu > 0 else float("nan")
        rel_mu = float((pred_mu_e / meas_mu_e) - 1.0) if meas_mu_e > 0 else float("nan")
        tol = float(tol_rel) if tol_rel > 0 else MASS_RATIO_TOL_REL_DEFAULT
        c2_tau = float((rel_tau / tol) ** 2) if np.isfinite(rel_tau) else float("nan")
        c2_mu = float((rel_mu / tol) ** 2) if np.isfinite(rel_mu) else float("nan")
        penalty = float(np.nansum([c2_tau, c2_mu]))
        detail = {
            "rel_tau_over_mu": rel_tau,
            "rel_mu_over_e": rel_mu,
            "chi2_tau_over_mu": c2_tau,
            "chi2_mu_over_e": c2_mu,
            "tolerance_rel": tol,
        }
        return penalty, detail
    except Exception:
        return float("nan"), {"error": float("nan")}


def _chi2_sf(*, chi2: float, dof: int) -> float:
    # Minimal chi2 survival function via scipy is avoided here; we use a tiny approximation:
    # for dof>=1 use numpy's regularized gammaincc via mpmath if needed; keep deterministic.
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


class FlavorJointObjectiveScanModule(TfptModule):
    module_id = "flavor_joint_objective_scan"
    title = "Flavor joint objective scan (CKM + PMNS χ² aggregation; discrete wiring)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "upstream module outputs in out/: ckm_full_pipeline/results.json, pmns_full_pipeline/results.json",
                "topology_phase_map/results.json (discrete candidate set size; docking point)",
                "flavor texture policy: tfpt_suite/data/flavor_texture_v24.json",
                "lepton masses: tfpt_suite/data/lepton_masses_pdg.json",
            ],
            outputs=[
                "joint objective table over discrete CKM variants (χ²_CKM + w1·χ²_PMNS + w2·Mass_Ratio_Penalty)",
                "best candidate under the declared objective",
                "diagnostic joint p-value",
            ],
            formulas=[
                "objective = chi2_ckm_refscale + w_pmns * chi2_pmns_mt + w_mass_ratio * chi2_mass_ratio",
            ],
            validation=[
                "objective is finite and computed from discrete upstream variants (no continuous fitter)",
                "module is explicit about missing upstream artifacts",
            ],
            determinism="Deterministic given upstream artifacts.",
            question="Given the discrete CKM variants and the PMNS best convention, what is the joint objective and which discrete variant wins?",
            objective=[
                "Provide a single joint score used by Physics-mode gates (prevents two separate heuristic narratives).",
                "Keep the search space strictly discrete (no free floats).",
            ],
            gaps=[
                "Full topology→phase map integration requires wiring topology_phase_map candidates into the CKM/PMNS generators (currently this module aggregates the existing discrete scans).",
            ],
        )

    def run(self, config) -> ModuleResult:
        out_dir = Path(config.output_dir)
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        flavor_cfg_path = Path(__file__).resolve().parent.parent / "data" / "flavor_texture_v24.json"
        flavor_cfg = _read_json(flavor_cfg_path)
        policy = flavor_cfg.get("joint_objective_policy", {}) if isinstance(flavor_cfg, dict) else {}
        w_pmns = float(policy.get("w_pmns", 1.0))
        w_mass_ratio = float(policy.get("w_mass_ratio", 1.0))
        mass_ratio_cfg = policy.get("mass_ratio_penalty", {}) if isinstance(policy.get("mass_ratio_penalty", {}), dict) else {}
        tol_rel = float(mass_ratio_cfg.get("tolerance_rel", MASS_RATIO_TOL_REL_DEFAULT))

        ckm_payload = _read_results_json(out_dir=out_dir, module_id="ckm_full_pipeline")
        pmns_payload = _read_results_json(out_dir=out_dir, module_id="pmns_full_pipeline")
        topo_payload = _read_results_json(out_dir=out_dir, module_id="topology_phase_map")

        checks: list[Check] = []
        if isinstance(flavor_cfg, dict):
            checks.append(mk_check_pass("flavor_policy_present", f"policy from {flavor_cfg_path.name}"))
        else:
            checks.append(mk_check_warn("flavor_policy_present", f"missing {flavor_cfg_path.name}; using defaults"))

        if isinstance(flavor_cfg, dict):
            topo_mode = str(flavor_cfg.get("topology_phase_atoms", {}).get("wiring", {}).get("phase_selection_rule_mode", ""))
            ckm_mode = str(flavor_cfg.get("ckm_variants", {}).get("phase_selection_rule_mode", ""))
            if topo_mode == "filter_only" and ckm_mode == "filter_only":
                checks.append(mk_check_pass("phase_selection_rule_filter_only", "topology + CKM phase selection are filter_only"))
            else:
                checks.append(
                    mk_check_fail(
                        "phase_selection_rule_filter_only",
                        f"topology_mode={topo_mode}, ckm_mode={ckm_mode}",
                    )
                )
        else:
            checks.append(mk_check_warn("phase_selection_rule_filter_only", "missing flavor_texture_v24.json"))
        if ckm_payload is None:
            checks.append(mk_check_fail("upstream_ckm_present", "missing ckm_full_pipeline/results.json (run ckm_full_pipeline first)"))
        else:
            checks.append(mk_check_pass("upstream_ckm_present", "ckm_full_pipeline/results.json present"))
        if pmns_payload is None:
            checks.append(mk_check_fail("upstream_pmns_present", "missing pmns_full_pipeline/results.json (run pmns_full_pipeline first)"))
        else:
            checks.append(mk_check_pass("upstream_pmns_present", "pmns_full_pipeline/results.json present"))
        if topo_payload is None:
            checks.append(mk_check_warn("topology_phase_map_present", "missing topology_phase_map/results.json (joint scan still runs; docking info unavailable)"))
        else:
            checks.append(mk_check_pass("topology_phase_map_present", "topology_phase_map/results.json present"))

        if ckm_payload is None or pmns_payload is None:
            return ModuleResult(
                results={"mode": mode, "status": "missing_upstream_outputs"},
                checks=checks,
                report="\n".join(
                    [
                        "Flavor joint objective scan",
                        "",
                        f"mode={mode}",
                        "ERROR: missing required upstream outputs (see checks).",
                        "",
                        "Checks:",
                        *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
                    ]
                ),
                warnings=[],
            )

        ckm_res = ckm_payload.get("results", {}) if isinstance(ckm_payload.get("results", {}), dict) else {}
        pmns_res = pmns_payload.get("results", {}) if isinstance(pmns_payload.get("results", {}), dict) else {}

        # CKM variants are emitted by `ckm_full_pipeline` under `results.flavor_texture.ckm_variant_scan.variants`.
        # Keep a small fallback for older/alternate schema keys.
        variants: list[object] = []
        try:
            if isinstance(ckm_res.get("flavor_texture", {}), dict):
                variants = (
                    ckm_res.get("flavor_texture", {}).get("ckm_variant_scan", {}).get("variants", [])  # type: ignore[assignment]
                )
            elif isinstance(ckm_res.get("ckm_texture", {}), dict):
                variants = ckm_res.get("ckm_texture", {}).get("ckm_variant_scan", {}).get("variants", [])  # type: ignore[assignment]
        except Exception:
            variants = []
        if not isinstance(variants, list) or not variants:
            checks.append(mk_check_fail("ckm_variants_present", "no CKM variants found in ckm_full_pipeline results"))
            variants = []
        else:
            checks.append(mk_check_pass("ckm_variants_present", f"{len(variants)} CKM variants available"))

        topo_present = isinstance(topo_payload, dict)
        if topo_present and variants:
            from_topo = [v for v in variants if isinstance(v, dict) and v.get("source_module") == "topology_phase_map"]
            if len(from_topo) == len(variants):
                checks.append(mk_check_pass("topology_filter_applied", "all CKM variants originate from topology_phase_map"))
            else:
                checks.append(
                    mk_check_warn(
                        "topology_filter_applied",
                        f"{len(from_topo)}/{len(variants)} variants tagged as topology_phase_map",
                    )
                )
        elif not topo_present:
            checks.append(mk_check_warn("topology_filter_applied", "missing topology_phase_map results; cannot verify filter-only topology usage"))

        chi2_pmns = None
        try:
            chi2_pmns = float(pmns_res.get("pmns_mt", {}).get("best_convention", {}).get("chi2"))
        except Exception:
            chi2_pmns = None
        if chi2_pmns is None or not np.isfinite(float(chi2_pmns)):
            checks.append(mk_check_fail("pmns_chi2_present", f"missing/non-finite pmns_mt.best_convention.chi2: {chi2_pmns}"))
            chi2_pmns = float("nan")
        else:
            checks.append(mk_check_pass("pmns_chi2_present", f"chi2_pmns_mt={chi2_pmns}"))

        # Mass ratio penalty (lepton ratios, TFPT delta_used anchor).
        delta_used = None
        try:
            delta_used = float(ckm_res.get("flavor_texture", {}).get("delta_used"))
        except Exception:
            delta_used = None
        if delta_used is None:
            delta_used = float(TfptConstants.compute().delta_star)
            checks.append(mk_check_warn("mass_ratio_anchor", f"delta_used missing; fallback delta_star={delta_used}"))
        else:
            checks.append(mk_check_pass("mass_ratio_anchor", f"delta_used={delta_used}"))

        mass_penalty, mass_detail = _mass_ratio_penalty_from_leptons(delta_used=delta_used, tol_rel=tol_rel)
        if np.isfinite(mass_penalty):
            checks.append(mk_check_pass("mass_ratio_penalty_computed", f"chi2_mass_ratio={mass_penalty:.6g}"))
        else:
            checks.append(mk_check_warn("mass_ratio_penalty_computed", "mass ratio penalty unavailable"))
            mass_penalty = 0.0

        rows: list[dict[str, object]] = []
        best = None
        for v in variants:
            if not isinstance(v, dict):
                continue
            label = str(v.get("label", "variant"))
            chi2_ckm = float(v.get("chi2_refscale", float("nan")))
            obj = float(chi2_ckm + w_pmns * float(chi2_pmns) + w_mass_ratio * float(mass_penalty))
            row = {
                "label": label,
                "s13_mode": str(v.get("s13_mode", "")),
                "delta_mode": str(v.get("delta_mode", "")),
                "chi2_ckm_refscale": chi2_ckm,
                "chi2_pmns_mt": float(chi2_pmns),
                "w_pmns": w_pmns,
                "w_mass_ratio": w_mass_ratio,
                "chi2_mass_ratio": float(mass_penalty),
                "objective": obj,
            }
            rows.append(row)
            if best is None or (np.isfinite(obj) and obj < float(best.get("objective", float("inf")))):
                best = dict(row)

        rows.sort(key=lambda r: float(r.get("objective", float("inf"))))
        checks.append(mk_check_pass("joint_objective_computed", f"{len(rows)} rows; best={best.get('label') if best else None}"))
        checks.append(mk_check_pass("joint_search_space_is_discrete", "objective aggregates discrete upstream scans; no continuous fitter"))

        # Diagnostic p-value (dashboard; dof ~ 9+4 by default).
        dof_ckm = 9
        dof_pmns = 4
        chi2_joint = float(best.get("objective", float("nan"))) if isinstance(best, dict) else float("nan")
        p_joint = _chi2_sf(chi2=chi2_joint, dof=int(dof_ckm + dof_pmns)) if np.isfinite(chi2_joint) else float("nan")
        if np.isfinite(p_joint):
            if mode == "physics" and p_joint < 0.05:
                checks.append(mk_check_fail("joint_pvalue_ok", f"p={p_joint} (<0.05) chi2={chi2_joint} dof={dof_ckm + dof_pmns}"))
            elif p_joint < 0.05:
                checks.append(mk_check_warn("joint_pvalue_ok", f"p={p_joint} (<0.05) chi2={chi2_joint} dof={dof_ckm + dof_pmns}"))
            else:
                checks.append(mk_check_pass("joint_pvalue_ok", f"p={p_joint} chi2={chi2_joint} dof={dof_ckm + dof_pmns}"))
        else:
            checks.append(mk_check_warn("joint_pvalue_ok", f"p-value not computed (chi2={chi2_joint})"))

        # Search-space sizing info (topology docking point)
        topo_pairs = None
        if isinstance(topo_payload, dict):
            topo_res = topo_payload.get("results", {}) if isinstance(topo_payload.get("results", {}), dict) else {}
            try:
                topo_pairs = int(len(topo_res.get("pairs", [])))
            except Exception:
                topo_pairs = None

        lines: list[str] = []
        lines += [
            "Flavor joint objective scan (discrete aggregation)",
            "",
            f"mode={mode}",
            "",
            f"PMNS: chi2_mt(best_convention) = {chi2_pmns}",
            f"Mass ratio penalty: chi2_mass_ratio = {mass_penalty}",
            f"CKM variants considered: {len(rows)}",
            f"Objective: chi2_ckm_refscale + {w_pmns} * chi2_pmns_mt + {w_mass_ratio} * chi2_mass_ratio",
            "",
            f"Best: {best}",
            f"Joint diagnostic: chi2={chi2_joint}, dof~{dof_ckm + dof_pmns}, p≈{p_joint}",
            f"Topology docking (pairs) = {topo_pairs}",
            "",
            "Top rows:",
        ]
        for r in rows[:5]:
            lines.append(f"- {r['label']}: objective={r['objective']:.6g} (chi2_ckm={r['chi2_ckm_refscale']:.6g})")
        lines += [
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "weights": {"w_pmns": w_pmns, "w_mass_ratio": w_mass_ratio},
                "pmns": {"chi2_mt": float(chi2_pmns)},
                "mass_ratio": {"chi2": float(mass_penalty), "detail": mass_detail, "tolerance_rel": tol_rel},
                "ckm": {"variants": rows},
                "best": best,
                "joint": {"chi2": chi2_joint, "dof_proxy": int(dof_ckm + dof_pmns), "p_value": p_joint},
                "topology_docking": {"pairs": topo_pairs},
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

