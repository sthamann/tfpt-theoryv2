from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

from tfpt_suite.constants import TfptConstants
from tfpt_suite.axion_inputs import resolve_axion_claim
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_info, mk_check_pass, mk_check_warn
from tfpt_suite.modules.axion_dm_pipeline import _omega_misalignment_h2, _parse_float_or_fraction


class TorsionDmPipelineModule(TfptModule):
    module_id = "torsion_dm_pipeline"
    title = "Torsion DM pipeline (optional; explicit placeholder)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["(optional) torsion regimes: tfpt_suite/data/torsion_regimes.json", "(optional) torsion bounds: tfpt_suite/data/torsion_bounds_vetted.json"],
            outputs=[
                "explicit status record (implemented vs not implemented)",
                "dependency ledger for a future torsion-as-DM closure module",
                "explicit constraints checklist (Ω_DM target, direct detection, torsion bounds)",
            ],
            formulas=[
                r"\rho_T \sim \frac12 m_T^2 S^2 + \cdots \;\; (model-dependent; requires a torsion excitation spectrum)",
            ],
            validation=[
                "Physics mode must not allow this to be silently green if DM closure depends on it.",
                "Constraints are documented explicitly (Ω_DM target, direct detection, torsion bounds).",
            ],
            determinism="Deterministic (no computation).",
            question="Is a torsion-excitation dark-matter channel implemented and constrained as part of TFPT’s DM closure?",
            objective=[
                "Make the optional torsion-as-DM branch explicit and testable (avoid narrative-only scope creep).",
            ],
            gaps=[
                "Requires a torsion excitation spectrum (mass/couplings), production mechanism, relic density, and comparison to bounds.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        # This module is explicitly optional: it should only become a physics-mode FAIL if the primary DM
        # closure path (axion sector) under-produces Ω_DM in the declared discrete policy.
        data_dir = Path(__file__).resolve().parent.parent / "data"
        ax_path = data_dir / "axion_tfpt_v106.json"
        ax = json.loads(ax_path.read_text(encoding="utf-8")) if ax_path.is_file() else {}

        cst = TfptConstants.compute()
        pol = ax.get("cosmo_policy", {}) if isinstance(ax.get("cosmo_policy", {}), dict) else {}
        scen_pol = pol.get("scenario_policy", {}) if isinstance(pol.get("scenario_policy", {}), dict) else {}
        scen_selected = str(scen_pol.get("selected", "")).strip()
        post = pol.get("post_inflation", {}) if isinstance(pol.get("post_inflation", {}), dict) else {}

        claim_resolved = resolve_axion_claim(ax_raw=ax, output_dir=Path(config.output_dir))
        f_a_GeV = float(claim_resolved.get("f_a_GeV", float("nan")))
        fa_source = str(claim_resolved.get("source", "quoted"))

        mis = pol.get("misalignment_model", {}) if isinstance(pol.get("misalignment_model", {}), dict) else {}
        omega_norm = float(mis.get("omega_norm_at_fa_5e11", 0.12))
        f_norm = float(mis.get("fa_norm_GeV", 5.0e11))
        p = float(mis.get("power_law_index", 1.165))
        omega_dm_ref = float(pol.get("omega_dm_h2_ref", 0.12))

        theta_varphi0 = float(cst.varphi0)
        theta_rms = float(post.get("theta_rms", math.pi / math.sqrt(3.0)))

        strings_dw_factor = _parse_float_or_fraction(post.get("strings_domain_walls_factor", 1.0), 1.0)
        strings_dw_policy = post.get("strings_domain_walls_factor_policy", {})
        if isinstance(strings_dw_policy, dict) and str(strings_dw_policy.get("kind", "")).strip() == "mobius_cusp_charge_sum":
            charges = strings_dw_policy.get("charges", [])
            if isinstance(charges, list) and charges:
                try:
                    strings_dw_factor = float(sum(Fraction(str(c)) for c in charges))
                except Exception:
                    pass
        strings_dw_factor = float(max(strings_dw_factor, 1.0))

        scenario_id = scen_selected or "pre_inflation_single_theta_varphi0"
        if scenario_id == "post_inflation_theta_rms_no_strings":
            theta_eff = theta_rms
            strings_factor = 1.0
        elif scenario_id == "post_inflation_theta_rms_with_strings_dw_factor":
            theta_eff = theta_rms
            strings_factor = strings_dw_factor
        else:
            theta_eff = theta_varphi0
            strings_factor = 1.0

        omega_a_h2 = (
            _omega_misalignment_h2(theta=theta_eff, f_a_GeV=f_a_GeV, omega_norm=omega_norm, f_norm_GeV=f_norm, p=p) * float(strings_factor)
        )
        frac_dm = omega_a_h2 / omega_dm_ref if omega_dm_ref > 0 else float("nan")

        tol_rel = 0.20
        axion_closes_dm = bool(math.isfinite(frac_dm) and abs(frac_dm - 1.0) <= tol_rel)

        if axion_closes_dm:
            chk = mk_check_pass(
                "torsion_dm_not_required",
                f"axion-first DM closure passes within {tol_rel:.0%}; torsion DM optional (Omega_a h^2≈{omega_a_h2:.3g}, fraction≈{frac_dm:.3g})",
            )
        else:
            chk = (
                mk_check_fail(
                    "torsion_dm_not_implemented",
                    f"axion-first DM closure leaves a gap (fraction≈{frac_dm:.3g}); implement torsion DM channel if ToE requires full Ω_DM closure",
                )
                if mode == "physics"
                else mk_check_warn(
                    "torsion_dm_not_implemented",
                    f"axion-first DM closure leaves a gap (fraction≈{frac_dm:.3g}); torsion DM remains an optional extension",
                )
            )

        lines = [
            "Torsion DM pipeline (optional)",
            "",
            f"mode={mode}",
            "",
            "Status: OPTIONAL BRANCH (only required if axion sector does not close Ω_DM)",
            "",
            "If implemented, this module should:",
            "- define a torsion excitation spectrum (mass m_T, degrees of freedom, coupling to SM)",
            "- compute a relic density (misalignment / freeze-in / freeze-out) under an explicit cosmology history",
            "- compare to Ω_DM target and report the remaining DM fraction",
            "- compare to torsion bounds (lab + astro) and direct-detection limits",
            "- document the observable channel(s) used for constraints (timing residuals, polarimetry, PSD, etc.)",
            "- emit a PASS/FAIL scorecard and a scenario matrix (torsion-only vs mixed DM)",
            "",
            "Axion-first closure diagnostic (same scaling law as axion_dm_pipeline; policy-driven):",
            f"- axion config: {ax_path}",
            f"- scenario_policy.selected = {scen_selected or '(missing)'}",
            f"- strings/domain-walls factor (effective) = {strings_factor:g}",
            f"- Ω_a h^2 ≈ {omega_a_h2:.3g} (fraction≈{frac_dm:.3g} of Ω_DM h^2≈{omega_dm_ref})",
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "axion_first_closure": {
                    "config_file": str(ax_path),
                    "scenario_selected": scen_selected or None,
                    "strings_domain_walls_factor_effective": float(strings_factor),
                    "fraction_of_dm": float(frac_dm) if math.isfinite(frac_dm) else None,
                    "tolerance_rel": float(tol_rel),
                    "closes_dm_within_tolerance": bool(axion_closes_dm),
                },
                "status": "optional_branch" if axion_closes_dm else "required_if_toe_requires_full_dm",
            },
            checks=[
                chk,
                mk_check_info("dm_branching_policy", "Torsion DM is optional unless axion-first closure fails under the declared discrete policy."),
                mk_check_info("torsion_dm_constraints_documented", "constraints checklist recorded in report"),
            ],
            report="\n".join(lines),
            warnings=[],
        )

