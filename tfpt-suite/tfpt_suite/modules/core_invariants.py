from __future__ import annotations

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


class CoreInvariantsModule(TfptModule):
    module_id = "core_invariants"
    title = "Core invariants (c3, varphi0, g_aγγ, beta)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["(none)"],
            outputs=["c3", "varphi0", "g_a_gamma_gamma", "beta_rad", "beta_deg"],
            formulas=[
                "c3 = 1/(8π)",
                "varphi0 = 1/(6π) + 3/(256π^4)",
                "g_{aγγ} = -4 c3 = -1/(2π)",
                "beta_rad = varphi0/(4π)",
                "beta_deg = (180/π) beta_rad",
            ],
            validation=[
                "algebraic identities hold to numerical precision",
                "beta_deg matches paper value (~0.2424°) within tight tolerance",
            ],
            determinism="Deterministic (mpmath precision controlled by config.mp_dps).",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        # validations (mpmath arithmetic should make these exact given enough dps)
        checks: list[Check] = []

        checks.append(
            Check(
                check_id="c3_definition",
                passed=mp.almosteq(c.c3, mp.mpf(1) / (8 * c.pi)),
                detail="c3 == 1/(8π)",
            )
        )
        checks.append(
            Check(
                check_id="varphi0_definition",
                passed=mp.almosteq(c.varphi0, mp.mpf(1) / (6 * c.pi) + mp.mpf(3) / (256 * c.pi**4)),
                detail="varphi0 == 1/(6π) + 3/(256π^4)",
            )
        )
        checks.append(
            Check(
                check_id="gagg_definition",
                passed=mp.almosteq(c.g_a_gamma_gamma, -mp.mpf(1) / (2 * c.pi)),
                detail="g_aγγ == -1/(2π)",
            )
        )
        checks.append(
            Check(
                check_id="beta_deg_value",
                passed=abs(c.beta_deg - mp.mpf("0.2424")) < mp.mpf("5e-4"),
                detail="beta_deg ~ 0.2424° (paper v2.4)",
            )
        )

        report = "\n".join(
            [
                "TFPT core invariants (paper v2.4 conventions)",
                "",
                f"mp.dps = {mp.dps}",
                "",
                f"c3                = {c.c3}",
                f"varphi0_tree      = {c.varphi0_tree}",
                f"delta_top         = {c.delta_top}",
                f"varphi0           = {c.varphi0}",
                f"b1                = {c.b1}",
                f"g_aγγ              = {c.g_a_gamma_gamma}",
                f"beta_rad          = {c.beta_rad}",
                f"beta_deg          = {c.beta_deg}",
                f"M/Mpl (R2 scale)   = {c.M_over_Mpl}",
                "",
                "Checks:",
                *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
                "",
            ]
        )

        return ModuleResult(
            results={"constants": c.as_dict()},
            checks=checks,
            report=report,
            warnings=[],
        )

