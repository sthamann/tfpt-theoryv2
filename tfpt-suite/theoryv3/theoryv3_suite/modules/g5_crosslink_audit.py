from __future__ import annotations

import json
from pathlib import Path

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.defect_partition import derive_delta2_from_defect_partition
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from theoryv3_suite.utils import ensure_ascii


def _plot_g5_links(*, out_dir: Path, g: float) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"g5_links_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        values = [g / 4.0, g / 2.0, g / (g + 1.0)]
        labels = ["g/4 (delta2)", "g/2 (delta_b3)", "g/(g+1) (gamma0)"]

        fig, ax = plt.subplots(figsize=(7.5, 3.4))
        ax.bar(labels, values, color="#718096")
        ax.set_ylabel("value")
        ax.set_title("g=5 crosslink ratios")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.tick_params(axis="x", rotation=15)

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "g5_links.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["g5_links_png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class G5CrosslinkAuditModule(TfptModule):
    module_id = "g5_crosslink_audit"
    title = "g=5 crosslink audit (delta2, unification patch, gamma0)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "defect partition g",
                "unification_gate_policy delta_b3 candidates",
                "TFPT gamma0",
            ],
            outputs=[
                "g/4, g/2, g/(g+1) ratios",
                "consistency checks for g=5 usage",
            ],
            formulas=[
                "delta2 factor = g/4",
                "unification patch candidate = g/2",
                "gamma0 = g/(g+1)",
            ],
            validation=[
                "g resolves to 5",
                "delta_b3 candidates include g/2",
                "gamma0 equals g/(g+1)",
            ],
            determinism="Deterministic (finite comparisons).",
            question="Is g=5 consistently threaded across multiple sectors?",
            objective=[
                "Summarize the discrete g=5 signature across delta2, unification, and E8 ladder.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        d2 = derive_delta2_from_defect_partition(delta_top=mp.mpf(c.delta_top))
        g_value = d2.justification.get("effective_multiplicity") if isinstance(d2.justification, dict) else None
        g_float = float(g_value) if g_value is not None else float("nan")

        policy_path = Path(__file__).resolve().parents[3] / "tfpt_suite" / "data" / "unification_gate_policy.json"
        policy = json.loads(policy_path.read_text(encoding="utf-8"))
        candidates = policy.get("delta_b3_candidates", [])
        has_g_over_2 = any(abs(float(x) - (g_float / 2.0)) < 1e-12 for x in candidates if isinstance(x, (int, float)))

        gamma0_pred = g_float / (g_float + 1.0) if g_float == g_float else float("nan")
        gamma0_ref = float(c.gamma0)

        checks: list[Check] = []
        checks.append(
            mk_check_pass("g_equals_5", f"g={g_value}")
            if g_value == 5
            else mk_check_warn("g_equals_5", f"g={g_value} (expected 5)")
        )
        checks.append(
            mk_check_pass("delta_b3_has_g_over_2", f"delta_b3 includes {g_float/2.0}")
            if has_g_over_2
            else mk_check_warn("delta_b3_has_g_over_2", "delta_b3 candidates missing g/2")
        )
        checks.append(
            mk_check_pass("gamma0_matches_g_over_g_plus_1", f"gamma0={gamma0_ref}")
            if abs(gamma0_pred - gamma0_ref) < 1e-12
            else mk_check_warn("gamma0_matches_g_over_g_plus_1", f"gamma0={gamma0_ref}, g/(g+1)={gamma0_pred}")
        )

        lines = [
            "g=5 crosslink audit",
            "",
            f"g = {g_value}",
            f"g/4 = {g_float/4.0}",
            f"g/2 = {g_float/2.0}",
            f"g/(g+1) = {gamma0_pred}",
            f"gamma0 (TFPT) = {gamma0_ref}",
            "",
            f"delta_b3 candidates include g/2: {has_g_over_2}",
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"g5_links_png": None}
        if getattr(config, "plot", True) and g_float == g_float:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_g5_links(out_dir=out_dir, g=g_float)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "g_value": g_value,
                "g_over_4": g_float / 4.0 if g_float == g_float else float("nan"),
                "g_over_2": g_float / 2.0 if g_float == g_float else float("nan"),
                "g_over_g_plus_1": gamma0_pred,
                "gamma0_ref": gamma0_ref,
                "delta_b3_candidates": candidates,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
