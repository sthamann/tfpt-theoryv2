from __future__ import annotations

from pathlib import Path

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.defect_partition import alpha_inv_0_from_delta2, derive_delta2_from_defect_partition
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from tfpt_suite.reference_ledger import get_dataset
from theoryv3_suite.utils import Z_TOL_DEFAULT, ensure_ascii


DELTA2_FACTORS = [
    ("half", mp.mpf("0.5")),
    ("one", mp.mpf("1.0")),
    ("g_over_4", mp.mpf("1.25")),
]
G_CANDIDATES = [4, 5, 6]


def _plot_alpha_defect(*, out_dir: Path, series: list[tuple[str, float]], alpha_ref: float) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"alpha_defect_series_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels = [s[0] for s in series]
        values = [s[1] for s in series]

        fig, ax = plt.subplots(figsize=(7.5, 4.0))
        ax.bar(labels, values, color="#2c7a7b")
        ax.axhline(float(alpha_ref), color="black", lw=1.0, ls="--", alpha=0.8, label="CODATA")
        ax.set_ylabel("alpha_inv_0")
        ax.set_title("Defect partition candidates for alpha_inv_0")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.legend(loc="best")

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "alpha_defect_series.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["alpha_defect_series_png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


def _plot_z_by_g(*, out_dir: Path, zs: list[tuple[int, float]]) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"g_negative_control.png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels = [str(g) for g, _ in zs]
        values = [z for _, z in zs]

        fig, ax = plt.subplots(figsize=(6.5, 3.4))
        ax.plot(labels, values, marker="o", lw=1.8, color="#2c7a7b")
        ax.axhline(0.0, color="black", lw=1.0, alpha=0.7)
        ax.set_xlabel("g")
        ax.set_ylabel("z (CODATA)")
        ax.set_title("Negative control: z vs g")
        ax.grid(True, ls=":", alpha=0.4)

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "g_negative_control.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["g_negative_control.png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class DefectPartitionG5AuditModule(TfptModule):
    module_id = "defect_partition_g5_audit"
    title = "Defect partition g=5 audit (delta2 -> alpha_inv_0)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT constants (delta_top)",
                "defect partition g from discrete enumeration",
                "CODATA alpha_inv_0 reference",
            ],
            outputs=[
                "delta2 = (g/4) delta_top^2",
                "alpha_inv_0 from delta2",
                "z score vs CODATA",
                "negative control: z(g) for g=4,5,6",
            ],
            formulas=[
                "delta2 = (g/4) * delta_top^2",
                "alpha_inv_0 from CFE + backreaction with delta2",
                "z = (pred - mean)/sigma",
            ],
            validation=[
                "g is discrete and equals 5 under current enumeration",
                "alpha_inv_0 is within 2 sigma of CODATA",
                "g=5 minimizes |z| among g in {4,5,6}",
            ],
            determinism="Deterministic (finite enumeration + root finding).",
            question="Does the g=5 defect partition close alpha(0) without a fit parameter?",
            objective=[
                "Lock delta2 to a discrete g=5 multiplicity.",
                "Show the resulting alpha_inv_0 is within CODATA sigma.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        d2 = derive_delta2_from_defect_partition(delta_top=mp.mpf(c.delta_top))
        alpha_inv = alpha_inv_0_from_delta2(delta2=d2.delta2, mp_dps=int(getattr(config, "mp_dps", 80)))

        ref = get_dataset("alpha_inv_codata_2022")
        alpha_ref = mp.mpf(str(ref["value"]))
        alpha_sigma = mp.mpf(str(ref["sigma"]))
        z = (alpha_inv - alpha_ref) / alpha_sigma if alpha_sigma != 0 else mp.mpf("nan")

        checks: list[Check] = []
        g_value = d2.justification.get("effective_multiplicity") if isinstance(d2.justification, dict) else None
        checks.append(
            mk_check_pass("g_equals_5", f"g={g_value}")
            if g_value == 5
            else mk_check_warn("g_equals_5", f"g={g_value} (expected 5)")
        )
        checks.append(
            mk_check_pass("alpha_within_2sigma", f"z={z}")
            if mp.isfinite(z) and abs(z) <= Z_TOL_DEFAULT
            else mk_check_warn("alpha_within_2sigma", f"z={z}")
        )

        series: list[tuple[str, float]] = []
        for label, factor in DELTA2_FACTORS:
            delta2 = factor * (mp.mpf(c.delta_top) ** 2)
            alpha_val = alpha_inv_0_from_delta2(delta2=delta2, mp_dps=int(getattr(config, "mp_dps", 80)))
            series.append((label, float(alpha_val)))

        z_by_g: list[tuple[int, float]] = []
        z_abs_by_g: dict[int, float] = {}
        for g in G_CANDIDATES:
            delta2_g = (mp.mpf(g) / mp.mpf(4)) * (mp.mpf(c.delta_top) ** 2)
            alpha_val = alpha_inv_0_from_delta2(delta2=delta2_g, mp_dps=int(getattr(config, "mp_dps", 80)))
            z_g = float((alpha_val - alpha_ref) / alpha_sigma) if alpha_sigma != 0 else float("nan")
            z_by_g.append((g, z_g))
            z_abs_by_g[g] = abs(z_g)

        z5 = z_abs_by_g.get(5, float("nan"))
        next_best = min([z_abs_by_g[g] for g in G_CANDIDATES if g != 5]) if len(G_CANDIDATES) > 1 else float("nan")
        gap = (next_best - z5) if (z5 == z5 and next_best == next_best) else float("nan")
        checks.append(
            mk_check_pass("g5_is_best_by_z", f"|z5|={z5}, others={z_abs_by_g}")
            if z5 == min(z_abs_by_g.values())
            else mk_check_warn("g5_is_best_by_z", f"|z5|={z5}, others={z_abs_by_g}")
        )
        checks.append(
            mk_check_pass("gap_to_next_best", f"gap={gap}")
            if gap == gap and gap >= float(Z_TOL_DEFAULT)
            else mk_check_warn("gap_to_next_best", f"gap={gap}")
        )

        lines = [
            "Defect partition g=5 audit",
            "",
            f"delta_top = {c.delta_top}",
            f"delta2 (g/4) = {d2.delta2} (model_id={d2.model_id})",
            f"alpha_inv_0 (g/4) = {alpha_inv}",
            f"alpha_inv_0 ref (CODATA) = {alpha_ref} ± {alpha_sigma} ({ref.get('version')})",
            f"z = {z}",
            "",
            "Candidate series:",
            *[f"- {label}: alpha_inv_0={val}" for label, val in series],
            "",
            "Negative control (z by g):",
            *[f"- g={g}: z={zv}" for g, zv in z_by_g],
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"alpha_defect_series_png": None, "g_negative_control.png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_alpha_defect(out_dir=out_dir, series=series, alpha_ref=float(alpha_ref))
            plot2, plot_warnings2 = _plot_z_by_g(out_dir=out_dir, zs=z_by_g)
            plot.update(plot2)
            plot_warnings.extend(plot_warnings2)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "reference": ref,
                "delta2": {
                    "model_id": d2.model_id,
                    "delta2": str(d2.delta2),
                    "delta2_over_delta_top2": str(d2.delta2_over_delta_top2),
                    "g_value": g_value,
                },
                "alpha_inv_0": {
                    "pred": str(alpha_inv),
                    "ref": str(alpha_ref),
                    "sigma": str(alpha_sigma),
                    "z": str(z),
                },
                "series": [{"label": label, "alpha_inv_0": val} for label, val in series],
                "negative_control": {
                    "g_candidates": G_CANDIDATES,
                    "z_by_g": [{"g": g, "z": zval} for g, zval in z_by_g],
                    "gap_to_next_best": gap,
                },
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
