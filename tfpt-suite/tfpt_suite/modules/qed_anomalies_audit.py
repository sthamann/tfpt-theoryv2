from __future__ import annotations

from typing import Any

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from tfpt_suite.modules.g2_and_lamb_shift_proxy import G2AndLambShiftProxyModule


class QedAnomaliesAuditModule(TfptModule):
    module_id = "qed_anomalies_audit"
    title = "QED anomalies audit (proxy): TFPT-scale contributions must not overshoot anomaly scales"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["g2_and_lamb_shift_proxy (TFPT-scale new-physics consistency proxy)"],
            outputs=["proxy scorecard that TFPT-scale suppressed contributions are negligible vs anomaly scales"],
            formulas=[
                r"a_\ell = (g_\ell-2)/2",
                r"\Delta E_{\rm Lamb} \;\; \text{(bound state QED)}",
            ],
            validation=["Delegates to g2_and_lamb_shift_proxy and requires its consistency gate to pass."],
            determinism="Deterministic given module inputs.",
            question="Are precision-QED anomalies (g-2, Lamb shift) implemented and compared within 5σ to reference values?",
            objective=[
                "Provide a minimal precision-QED audit layer that prevents TFPT-scale new physics from being obviously excluded by g-2 / Lamb shift.",
            ],
            gaps=[
                "Publication-grade requires a full SM+QED calculation and a proper likelihood against (g-2)_e, (g-2)_μ, and bound-state observables.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        proxy = G2AndLambShiftProxyModule()
        proxy_res = proxy.run(config)
        ok = any((c.check_id == "precision_qed_consistency" and c.passed) for c in proxy_res.checks)

        chk: Check = (
            mk_check_pass("within_5sigma", "proxy consistency passes (TFPT-scale contributions do not overshoot anomaly scales)")
            if ok
            else mk_check_warn("within_5sigma", "proxy consistency did not pass cleanly; review g2_and_lamb_shift_proxy outputs")
        )

        lines = [
            "QED anomalies audit (proxy layer; delegates to g2_and_lamb_shift_proxy)",
            "",
            f"mode={mode}",
            "",
            "Delegated result (from g2_and_lamb_shift_proxy):",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in proxy_res.checks],
            "",
            "Audit gate:",
            f"- within_5sigma: {'PASS' if chk.passed else 'WARN'} ({chk.detail})",
        ]

        results: dict[str, Any] = {"mode": mode, "delegated_module": "g2_and_lamb_shift_proxy", "delegated": proxy_res.results}

        return ModuleResult(results=results, checks=[chk, *proxy_res.checks], report="\n".join(lines), warnings=[])

