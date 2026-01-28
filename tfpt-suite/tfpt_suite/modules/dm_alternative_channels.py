from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

from tfpt_suite.axion_inputs import resolve_axion_claim
from tfpt_suite.constants import TfptConstants
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
from tfpt_suite.modules.axion_dm_pipeline import _omega_misalignment_h2, _parse_float_or_fraction


class DmAlternativeChannelsModule(TfptModule):
    module_id = "dm_alternative_channels"
    title = "DM alternative channels (torsion excitations; explicit placeholder)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[],
            outputs=["explicit status record (implemented vs not implemented)", "pointer to torsion_dm_pipeline placeholder"],
            formulas=[],
            validation=["Physics mode must not allow this to be silently green if DM closure depends on it."],
            determinism="Deterministic (no computation).",
            question="Are non-axion TFPT-native dark matter channels (e.g. torsion excitations) implemented and constrained?",
            objective=[
                "Track alternative DM channels explicitly to avoid hidden scope creep and narrative-only claims.",
            ],
            gaps=[
                "Implement a concrete torsion excitation spectrum + production mechanism + relic density + bounds (see torsion_dm_pipeline).",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        data_dir = Path(__file__).resolve().parent.parent / "data"
        ax_path = data_dir / "axion_tfpt_v106.json"
        ax = json.loads(ax_path.read_text(encoding="utf-8")) if ax_path.is_file() else {}

        # Reuse the same (engineering-level) misalignment scaling as axion_dm_pipeline to decide whether
        # *additional* TFPT-native DM channels are required.
        cst = TfptConstants.compute()
        pol = ax.get("cosmo_policy", {}) if isinstance(ax.get("cosmo_policy", {}), dict) else {}
        scen_pol = pol.get("scenario_policy", {}) if isinstance(pol.get("scenario_policy", {}), dict) else {}
        scen_selected = str(scen_pol.get("selected", "")).strip()
        post = pol.get("post_inflation", {}) if isinstance(pol.get("post_inflation", {}), dict) else {}

        claim_resolved = resolve_axion_claim(ax_raw=ax, output_dir=Path(config.output_dir))
        f_a_GeV = float(claim_resolved.get("f_a_GeV", float("nan")))

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
        omega_missing_h2 = max(0.0, float(omega_dm_ref) - float(omega_a_h2)) if math.isfinite(omega_a_h2) else float("nan")

        # Conservative diagnostic acceptance: we want a falsifiable closure path without continuous tuning.
        # Treat 20% relative agreement as a "PASS" threshold for this engineering-level relic scaling.
        tol_rel = 0.20
        ok = bool(math.isfinite(frac_dm) and abs(frac_dm - 1.0) <= tol_rel)
        if not math.isfinite(frac_dm):
            chk = (mk_check_fail if mode == "physics" else mk_check_warn)("omega_dm_match", f"non-finite fraction_of_dm={frac_dm}")
        elif ok:
            chk = mk_check_pass("omega_dm_match", f"Omega_a h^2≈{omega_a_h2:.3g} matches Omega_DM h^2≈{omega_dm_ref} within {tol_rel:.0%} (fraction≈{frac_dm:.3g})")
        else:
            # If axion under-produces in this deterministic policy, the torsion DM channel becomes a required extension.
            chk = (
                mk_check_fail("omega_dm_match", f"axion leaves Ω_missing h^2≈{omega_missing_h2:.3g} (fraction≈{frac_dm:.3g}); implement torsion_dm_pipeline")
                if mode == "physics"
                else mk_check_warn(
                    "omega_dm_match",
                    f"axion leaves Ω_missing h^2≈{omega_missing_h2:.3g} (fraction≈{frac_dm:.3g}); torsion DM remains an optional extension",
                )
            )

        lines = [
            "DM alternative channels (axion-first closure + optional torsion branch)",
            "",
            f"mode={mode}",
            "",
            "Primary DM path: axion (post-inflation policy), checked against Ω_DM h^2 target.",
            "",
            "Notes:",
            f"- axion config: {ax_path}",
            f"- scenario_policy.selected = {scen_selected or '(missing)'}",
            f"- strings/domain-walls factor (effective) = {strings_factor:g}",
            f"- Ω_a h^2 ≈ {omega_a_h2:.3g} (fraction≈{frac_dm:.3g} of Ω_DM h^2≈{omega_dm_ref})",
            f"- Ω_missing h^2 ≈ {omega_missing_h2:.3g}",
            "- If Ω_missing>0 in physics mode, the torsion DM branch becomes required (implement torsion_dm_pipeline).",
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "axion_policy": {
                    "config_file": str(ax_path),
                    "scenario_selected": scen_selected or None,
                    "theta_eff": float(theta_eff),
                    "strings_domain_walls_factor_effective": float(strings_factor),
                    "strings_domain_walls_factor_policy": strings_dw_policy if isinstance(strings_dw_policy, dict) else None,
                },
                "dm_closure": {
                    "Omega_DM_h2_ref": float(omega_dm_ref),
                    "Omega_a_h2": float(omega_a_h2),
                    "fraction_of_dm": float(frac_dm) if math.isfinite(frac_dm) else None,
                    "Omega_missing_h2": float(omega_missing_h2) if math.isfinite(omega_missing_h2) else None,
                    "tolerance_rel": float(tol_rel),
                },
                "recommended_module_if_missing": "torsion_dm_pipeline" if (math.isfinite(omega_missing_h2) and omega_missing_h2 > 0) else None,
            },
            checks=[chk, mk_check_info("dm_branching_policy", "Axion-first closure; torsion DM is only required if Ω_missing>0 in physics mode.")],
            report="\n".join(lines),
            warnings=[],
        )

