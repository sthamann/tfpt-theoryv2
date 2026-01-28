from __future__ import annotations

from pathlib import Path

from mpmath import mp

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass, mk_check_warn
from tfpt_suite.reference_ledger import get_dataset
from theoryv3_suite.utils import ensure_ascii, read_json_if_exists


Z_TOL_TM1 = 3.0


class PmnsTm1AuditModule(TfptModule):
    module_id = "pmns_tm1_audit"
    title = "PMNS TM1 audit (theta12 from theta13)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "sin2_theta13 from flavor_pattern_audit",
                "NuFIT sin2_theta12 reference (normal ordering)",
            ],
            outputs=[
                "sin2_theta12 from TM1 sum rule",
                "z score vs NuFIT reference",
            ],
            formulas=[
                "sin2_theta12 = (1/3) * (1 - 2 * sin2_theta13)",
            ],
            validation=[
                "sin2_theta12 within conservative band (3 sigma)",
            ],
            determinism="Deterministic (closed form).",
            question="Does the TM1 sum rule yield a viable sin2(theta12) given sin2(theta13)?",
            objective=[
                "Audit the TM1 sum rule with explicit inputs and a conservative band.",
            ],
        )

    def run(self, config) -> ModuleResult:
        exp_path = Path(getattr(config, "output_dir", "")) / "flavor_pattern_audit" / "results.json"
        payload = read_json_if_exists(exp_path) or {}
        sin2_theta13_val = payload.get("results", {}).get("sin2_theta13", None)
        sin2_theta13 = mp.mpf(str(sin2_theta13_val)) if sin2_theta13_val is not None else mp.mpf("nan")

        sin2_theta12 = (mp.mpf(1) / 3) * (mp.mpf(1) - mp.mpf(2) * sin2_theta13)

        ref = get_dataset("pmns_sin2_theta12_nufit53_no")
        sin2_ref = mp.mpf(str(ref["value"]))
        sin2_sigma = mp.mpf(str(ref["sigma"]))
        z = (sin2_theta12 - sin2_ref) / sin2_sigma if sin2_sigma != 0 else mp.mpf("nan")

        checks: list[Check] = []
        checks.append(
            mk_check_pass("sin2_theta12_within_band", f"z={z}")
            if mp.isfinite(z) and abs(z) <= Z_TOL_TM1
            else mk_check_warn("sin2_theta12_within_band", f"z={z}")
        )
        checks.append(
            mk_check_info(
                "theta23_placeholder",
                "theta23 and delta_CP are placeholders until Z3-breaking is wired into this audit",
            )
        )

        lines = [
            "PMNS TM1 audit",
            "",
            f"sin2_theta13 (input) = {sin2_theta13}",
            f"sin2_theta12 (TM1) = {sin2_theta12}",
            f"sin2_theta12 ref = {sin2_ref} ± {sin2_sigma} ({ref.get('version')})",
            f"z = {z}",
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        return ModuleResult(
            results={
                "reference": ref,
                "sin2_theta13": str(sin2_theta13),
                "sin2_theta12": str(sin2_theta12),
                "z": str(z),
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )
