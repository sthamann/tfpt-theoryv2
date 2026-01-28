from __future__ import annotations

from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass, mk_check_warn
from tfpt_suite.modules.alpha_precision_audit import _fixed_point_alpha, _ppm
from tfpt_suite.reference_ledger import get_dataset
from theoryv3_suite.utils import ensure_ascii


K_GRID = [mp.mpf(0), mp.mpf(1), mp.mpf("1.5"), mp.mpf(2), mp.mpf("2.5"), mp.mpf(3)]


def _plot_ppm_vs_k(*, out_dir: Path, series: list[dict[str, Any]]) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"alpha_backreaction_ppm.png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        ks = [float(row["k"]) for row in series]
        ppm = [float(row["ppm_vs_codata"]) for row in series]

        fig, ax = plt.subplots(figsize=(7.5, 3.6))
        ax.plot(ks, ppm, marker="o", lw=1.8, color="#2c5282")
        ax.axhline(0.0, color="black", lw=1.0, alpha=0.7)
        ax.axvline(2.0, color="black", lw=1.0, ls="--", alpha=0.7, label="k=2")
        ax.set_xlabel("k (backreaction exponent)")
        ax.set_ylabel("ppm vs CODATA")
        ax.set_title("alpha backreaction sensitivity")
        ax.grid(True, ls=":", alpha=0.4)
        ax.legend(loc="best")

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "alpha_backreaction_ppm.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["alpha_backreaction_ppm.png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class AlphaBackreactionSensitivityAuditModule(TfptModule):
    module_id = "alpha_backreaction_sensitivity_audit"
    title = "Alpha backreaction sensitivity audit (k sweep)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants (c3, varphi0_tree, delta_top, b1)",
                "CODATA alpha_inv_0 reference",
            ],
            outputs=[
                "alpha_inv_0 for k in {0,1,1.5,2,2.5,3}",
                "ppm deviation vs CODATA",
            ],
            formulas=[
                "varphi(alpha)=varphi_tree + delta_top * exp(-k alpha)",
                "alpha from CFE fixed point (self-consistent)",
                "ppm = 1e6*(pred-ref)/ref",
            ],
            validation=[
                "k=2 is near the minimum |ppm| in the grid",
            ],
            determinism="Deterministic (fixed-point iteration).",
            question="Is k=2 a natural choice from alpha backreaction sensitivity?",
            objective=[
                "Quantify how alpha_inv_0 depends on k and show the k=2 near-crossing.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        ref = get_dataset("alpha_inv_codata_2022")
        alpha_ref = mp.mpf(str(ref["value"]))

        series: list[dict[str, Any]] = []
        for k in K_GRID:
            sol = _fixed_point_alpha(
                c3=c.c3,
                b1=c.b1,
                varphi_tree=c.varphi0_tree,
                delta_top=c.delta_top,
                k=k,
            )
            ppm = _ppm(sol.alpha_inv, alpha_ref)
            series.append(
                {
                    "k": str(k),
                    "alpha_inv": sol.alpha_inv,
                    "ppm_vs_codata": ppm,
                    "iterations": sol.iterations,
                    "converged": sol.converged,
                }
            )

        abs_ppm = {float(row["k"]): abs(float(row["ppm_vs_codata"])) for row in series}
        min_k = min(abs_ppm, key=abs_ppm.get)
        k2_is_best = abs_ppm.get(2.0) == abs_ppm.get(min_k)

        diffs = [float(series[i + 1]["ppm_vs_codata"]) - float(series[i]["ppm_vs_codata"]) for i in range(len(series) - 1)]
        monotonic = all(d >= 0 for d in diffs) or all(d <= 0 for d in diffs)

        checks: list[Check] = []
        checks.append(
            mk_check_pass("k2_is_near_zero_crossing", f"min_k={min_k}, abs_ppm={abs_ppm}")
            if k2_is_best
            else mk_check_warn("k2_is_near_zero_crossing", f"min_k={min_k}, abs_ppm={abs_ppm}")
        )
        checks.append(mk_check_info("monotonicity_check", f"monotonic={monotonic}, diffs={diffs}"))

        lines = [
            "Alpha backreaction sensitivity audit",
            "",
            f"alpha_ref = {alpha_ref} ({ref.get('version')})",
            "k grid results:",
            *[f"- k={row['k']}: alpha_inv={row['alpha_inv']}, ppm={row['ppm_vs_codata']}" for row in series],
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"alpha_backreaction_ppm.png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_ppm_vs_k(out_dir=out_dir, series=series)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "reference": ref,
                "series": series,
                "abs_ppm": abs_ppm,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
