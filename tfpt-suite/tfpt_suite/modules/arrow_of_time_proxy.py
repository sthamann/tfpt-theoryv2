from __future__ import annotations

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass
from tfpt_suite.modules.arrow_mechanism import ArrowMechanismModule


class ArrowOfTimeProxyModule(TfptModule):
    module_id = "arrow_of_time_proxy"
    title = "Arrow of time proxy (explicitly not implemented; ToE gap marker)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[],
            outputs=[
                "explicit status record (implemented vs not implemented)",
                "dependency ledger for a future irreversible torsion-flux / topology module",
            ],
            formulas=[
                r"S_{\rm prod} = \int d^4x \, \sigma(x) \ge 0  (entropy production proxy; requires a non-invertible mechanism)",
            ],
            validation=["Physics mode must not allow this to be silently green."],
            determinism="Deterministic (no computation).",
            question="Is an explicit, testable arrow-of-time mechanism implemented as a falsifiable module?",
            objective=[
                "Avoid narrative-only claims by making the missing arrow-of-time mechanism explicit and testable.",
                "Provide a concrete checklist for the minimal 'irreversibility' closure work.",
            ],
            gaps=[
                "Requires a concrete irreversible mechanism (e.g. torsion flux with non-invertible transitions) and an observable proxy (entropy production, CPT-odd signal, etc.).",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        mech = ArrowMechanismModule()
        mech_res = mech.run(config)
        ok = any((c.check_id == "entropy_production_testable" and c.passed) for c in mech_res.checks)
        chk: Check = mk_check_pass("arrow_mechanism_present", "arrow_mechanism implemented and entropy_production_testable evaluated") if ok else mk_check_info(
            "arrow_mechanism_present", "arrow_mechanism ran; see entropy_production_testable for pass/fail"
        )

        lines = [
            "Arrow of time proxy (alias to arrow_mechanism)",
            "",
            f"mode={mode}",
            "",
            "Delegated result (from arrow_mechanism):",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in mech_res.checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "delegated_module": "arrow_mechanism",
                "delegated": mech_res.results,
            },
            checks=[chk, *mech_res.checks],
            report="\n".join(lines),
            warnings=[],
        )

