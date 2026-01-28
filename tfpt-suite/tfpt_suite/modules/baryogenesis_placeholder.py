from __future__ import annotations

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass
from tfpt_suite.modules.baryogenesis_mechanism import BaryogenesisMechanismModule


class BaryogenesisPlaceholderModule(TfptModule):
    module_id = "baryogenesis_placeholder"
    title = "Baryogenesis placeholder (explicitly not implemented)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[],
            outputs=[
                "explicit status record (implemented vs not implemented)",
                "dependency ledger for a future baryogenesis mechanism module",
            ],
            formulas=[],
            validation=[
                "Physics mode must not allow this to be silently green.",
            ],
            determinism="Deterministic (no computation).",
            question="Is a baryogenesis mechanism implemented as a falsifiable module?",
            objective=[
                "Avoid narrative-only claims by making the absence of a baryogenesis mechanism explicit and testable.",
            ],
            gaps=[
                "Requires a concrete mechanism module (e.g. anomaly/inflow + out-of-equilibrium dynamics) with a falsifiable η_b output.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        # Backward-compatible alias: run the concrete mechanism module and expose a simple gate.
        mech = BaryogenesisMechanismModule()
        mech_res = mech.run(config)
        ok = any((c.check_id == "eta_b_match" and c.passed) for c in mech_res.checks)
        chk: Check = mk_check_pass("baryogenesis_mechanism_present", "baryogenesis_mechanism implemented and eta_b_match gate evaluated") if ok else mk_check_info(
            "baryogenesis_mechanism_present", "baryogenesis_mechanism ran; see eta_b_match for pass/fail"
        )

        lines = [
            "Baryogenesis placeholder (alias to baryogenesis_mechanism)",
            "",
            f"mode={mode}",
            "",
            "Status: implemented via baryogenesis_mechanism",
            "",
            "Delegated result (from baryogenesis_mechanism):",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in mech_res.checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "delegated_module": "baryogenesis_mechanism",
                "delegated": mech_res.results,
            },
            checks=[chk, *mech_res.checks],
            report="\n".join(lines),
            warnings=[],
        )

