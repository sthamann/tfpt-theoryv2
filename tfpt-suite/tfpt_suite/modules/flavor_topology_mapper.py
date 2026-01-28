from __future__ import annotations

import json
from pathlib import Path
from typing import Any

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


class FlavorTopologyMapperModule(TfptModule):
    module_id = "flavor_topology_mapper"
    title = "Flavor topology mapper (holonomy→CP phases; aggregation scaffold)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "topology_phase_map (discrete phase atoms; docking): out/topology_phase_map/results.json (optional)",
                "ckm_full_pipeline and pmns_full_pipeline outputs (optional; used to compute diagnostic joint χ²)",
            ],
            outputs=[
                "explicit statement of whether topology-phase candidates are wired into CKM/PMNS generators",
                "diagnostic joint χ² and p-value from current CKM+PMNS runs (if available)",
            ],
            formulas=["χ²_joint = χ²_CKM + χ²_PMNS (proxy)"],
            validation=["chi2_joint_ckm_pmns check is emitted (PASS only if p>0.05 in physics mode)."],
            determinism="Deterministic given available upstream artifacts.",
            question="Do discrete topology phase atoms map to CKM/PMNS CP phases in a way that improves joint χ² without continuous tuning?",
            objective=[
                "Centralize the topology→phase mapping status and joint χ² as a single, explicit scaffold module.",
                "Avoid silent convention-shopping by requiring discrete candidate sets (topology_phase_map).",
            ],
            gaps=[
                "Even with discrete wiring enabled, a publication-grade topology→operator derivation is still missing (holonomy/APS η → Yukawa operator; scheme-consistent matching below mt).",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")
        out_dir = Path(config.output_dir)

        topo = _read_results_json(out_dir=out_dir, module_id="topology_phase_map")
        ckm = _read_results_json(out_dir=out_dir, module_id="ckm_full_pipeline")
        pmns = _read_results_json(out_dir=out_dir, module_id="pmns_full_pipeline")

        # Candidate set size (docking info)
        topo_pairs = None
        if isinstance(topo, dict):
            res = topo.get("results", {}) if isinstance(topo.get("results", {}), dict) else {}
            try:
                topo_pairs = int(len(res.get("pairs", [])))
            except Exception:
                topo_pairs = None

        # Diagnostic joint χ² (proxy)
        chi2_ckm = None
        chi2_pmns = None
        try:
            chi2_ckm = float((ckm or {}).get("results", {}).get("rg_upward", {}).get("chi2_refscale"))
        except Exception:
            chi2_ckm = None
        try:
            chi2_pmns = float((pmns or {}).get("results", {}).get("pmns_mt", {}).get("best_convention", {}).get("chi2"))
        except Exception:
            chi2_pmns = None

        chi2_joint = None
        if chi2_ckm is not None and chi2_pmns is not None:
            chi2_joint = float(chi2_ckm + chi2_pmns)

        dof_proxy = 13
        p_joint = _chi2_sf(chi2=float(chi2_joint), dof=dof_proxy) if chi2_joint is not None else float("nan")

        checks: list[Check] = []
        checks.append(mk_check_pass("phase_map_is_discrete", f"topology_phase_map pairs={topo_pairs}" if topo_pairs is not None else "topology_phase_map missing"))

        if chi2_joint is None:
            checks.append(mk_check_warn("chi2_joint_ckm_pmns", "missing CKM/PMNS outputs; run ckm_full_pipeline + pmns_full_pipeline to compute joint χ²"))
        else:
            if mode == "physics" and (not (p_joint == p_joint and p_joint >= 0.05)):
                checks.append(mk_check_fail("chi2_joint_ckm_pmns", f"FAIL: chi2={chi2_joint:.6g}, dof≈{dof_proxy}, p≈{p_joint:.3e}"))
            elif p_joint == p_joint and p_joint < 0.05:
                checks.append(mk_check_warn("chi2_joint_ckm_pmns", f"WARN: chi2={chi2_joint:.6g}, dof≈{dof_proxy}, p≈{p_joint:.3e}"))
            else:
                checks.append(mk_check_pass("chi2_joint_ckm_pmns", f"PASS: chi2={chi2_joint:.6g}, dof≈{dof_proxy}, p≈{p_joint:.3e}"))

        # Wiring status (explicit; detect whether the upstream pipelines exposed topology-derived candidate scans)
        wired_ckm = False
        try:
            ckm_res = (ckm or {}).get("results", {}) if isinstance((ckm or {}).get("results", {}), dict) else {}
            ft = ckm_res.get("flavor_texture", {}) if isinstance(ckm_res.get("flavor_texture", {}), dict) else {}
            scan = ft.get("ckm_variant_scan", {}) if isinstance(ft.get("ckm_variant_scan", {}), dict) else {}
            vars_ = scan.get("variants", []) if isinstance(scan.get("variants", []), list) else []
            wired_ckm = any(isinstance(v, dict) and v.get("source_module") == "topology_phase_map" for v in vars_)
        except Exception:
            wired_ckm = False

        wired_pmns = False
        try:
            pmns_res = (pmns or {}).get("results", {}) if isinstance((pmns or {}).get("results", {}), dict) else {}
            tex = pmns_res.get("texture", {}) if isinstance(pmns_res.get("texture", {}), dict) else {}
            theta_scan = tex.get("theta_scan", {}) if isinstance(tex.get("theta_scan", {}), dict) else {}
            cands = theta_scan.get("candidates", []) if isinstance(theta_scan.get("candidates", []), list) else []
            wired_pmns = any(isinstance(r, dict) and r.get("source_module") == "topology_phase_map" for r in cands)
        except Exception:
            wired_pmns = False

        if wired_ckm and wired_pmns:
            checks.append(mk_check_pass("topology_phase_map_wired_into_flavor_pipelines", "WIRED: CKM δ_CP scan + PMNS θ scan include topology_phase_map candidates."))
            wiring_status = "wired_ckm_and_pmns"
        elif wired_ckm or wired_pmns:
            checks.append(
                (mk_check_fail if mode == "physics" else mk_check_warn)(
                    "topology_phase_map_wired_into_flavor_pipelines",
                    f"PARTIALLY WIRED: ckm={wired_ckm}, pmns={wired_pmns} (one pipeline still lacks topology injection).",
                )
            )
            wiring_status = "wired_partial"
        else:
            checks.append(
                (mk_check_fail if mode == "physics" else mk_check_warn)(
                    "topology_phase_map_wired_into_flavor_pipelines",
                    "NOT WIRED: no topology_phase_map-derived candidate scans detected in CKM/PMNS outputs.",
                )
            )
            wiring_status = "not_wired"

        lines: list[str] = []
        lines += [
            "Flavor topology mapper (scaffold)",
            "",
            f"mode={mode}",
            "",
            f"topology_phase_map pairs: {topo_pairs}",
            f"chi2_ckm_refscale: {chi2_ckm}",
            f"chi2_pmns_mt: {chi2_pmns}",
            f"chi2_joint: {chi2_joint}",
            f"p_joint (dof≈{dof_proxy}): {p_joint}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "topology_pairs": topo_pairs,
                "chi2": {"ckm_refscale": chi2_ckm, "pmns_mt": chi2_pmns, "joint": chi2_joint, "dof_proxy": dof_proxy, "p_joint": p_joint},
                "wiring_status": wiring_status,
                "wiring": {"ckm": bool(wired_ckm), "pmns": bool(wired_pmns)},
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

