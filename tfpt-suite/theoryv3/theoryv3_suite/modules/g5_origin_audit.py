from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from tfpt_suite.modules.mobius_cusp_classification import _su5_hypercharge_spectrum
from theoryv3_suite.utils import ensure_ascii


def _plot_g5_origin(*, out_dir: Path, counts: dict[str, int]) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"g5_origin.png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels = list(counts.keys())
        values = [counts[k] for k in labels]

        fig, ax = plt.subplots(figsize=(6.5, 3.2))
        ax.bar(labels, values, color="#4a5568")
        ax.set_ylabel("degeneracy")
        ax.set_title("SU(5) holonomy degeneracies (origin for g)")
        ax.grid(True, axis="y", ls=":", alpha=0.4)

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "g5_origin.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["g5_origin.png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class G5OriginAuditModule(TfptModule):
    module_id = "g5_origin_audit"
    title = "g=5 origin audit (single-source SU(5) holonomy degeneracy)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["SU(5) hypercharge holonomy spectrum (fundamental)"],
            outputs=["degeneracy counts", "g from degeneracy sum"],
            formulas=[
                "g := sum of degeneracies in SU(5) fundamental holonomy spectrum",
                "degeneracies inferred from repeated eigenvalues",
            ],
            validation=[
                "g equals 5 (3 color + 2 weak holonomy channels)",
            ],
            determinism="Deterministic (no fitting, no scanning).",
            question="Can g be derived from a single discrete origin (SU(5) holonomy degeneracy)?",
            objective=[
                "Fix g from one source only, not by crosslinking multiple sectors.",
            ],
        )

    def run(self, config) -> ModuleResult:
        hol = _su5_hypercharge_spectrum()
        eigenvalues = list(hol.eigenvalues_fund)
        counts_raw = Counter([str(v) for v in eigenvalues])
        counts = {k: int(v) for k, v in counts_raw.items()}
        g_value = sum(counts.values())

        degens_sorted = sorted(counts.values(), reverse=True)
        expected_degens = [3, 2]
        matches_degeneracy = degens_sorted == expected_degens

        checks: list[Check] = []
        checks.append(
            mk_check_pass("g_equals_5_from_holonomy", f"g={g_value}")
            if g_value == 5
            else mk_check_warn("g_equals_5_from_holonomy", f"g={g_value} (expected 5)")
        )
        checks.append(
            mk_check_pass("degeneracy_pattern", f"degeneracies={degens_sorted}")
            if matches_degeneracy
            else mk_check_warn("degeneracy_pattern", f"degeneracies={degens_sorted} (expected {expected_degens})")
        )

        lines = [
            "g=5 origin audit (SU(5) holonomy degeneracy)",
            "",
            f"eigenvalues_fund = {eigenvalues}",
            f"degeneracies = {counts}",
            f"g = {g_value}",
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"g5_origin.png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_g5_origin(out_dir=out_dir, counts=counts)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "eigenvalues_fund": [str(v) for v in eigenvalues],
                "degeneracies": counts,
                "g": g_value,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
