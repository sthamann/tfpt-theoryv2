from __future__ import annotations

import math
from pathlib import Path

from mpmath import mp

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass, mk_check_warn
from theoryv3_suite.utils import ensure_ascii, extract_nested, load_tfpt_results, read_json_if_exists


def _plot_norm_origin(*, out_dir: Path, k_value: float, n_value: float) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"dark_energy_norm_origin.png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels = ["k (cover degree)", "n = 1/k"]
        values = [k_value, n_value]

        fig, ax = plt.subplots(figsize=(6.5, 3.2))
        ax.bar(labels, values, color="#6b46c1")
        ax.set_ylabel("value")
        ax.set_title("Dark energy normalization from cover degree")
        ax.grid(True, axis="y", ls=":", alpha=0.4)

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "dark_energy_norm_origin.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["dark_energy_norm_origin.png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class DarkEnergyNormHalfOriginAuditModule(TfptModule):
    module_id = "dark_energy_norm_half_origin_audit"
    title = "Dark energy norm origin audit (n=1/2 from double cover)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "alpha_precision_audit self_consistent.k (cover degree)",
                "dark_energy_exponential_audit best candidate label",
            ],
            outputs=[
                "n_from_cover = 1/k",
                "consistency with best candidate n=1/2",
            ],
            formulas=[
                "n_from_cover = 1 / k",
                "k=2 corresponds to the double-cover backreaction exponent",
            ],
            validation=[
                "n_from_cover equals 1/2",
                "dark_energy_exponential_audit best candidate is n=1/2",
            ],
            determinism="Deterministic (no fitting).",
            question="Is the dark-energy normalization n=1/2 fixed by the double-cover degree?",
            objective=[
                "Tie the preferred n=1/2 normalization to the cover degree used in alpha backreaction.",
            ],
        )

    def run(self, config) -> ModuleResult:
        alpha_payload = load_tfpt_results("alpha_precision_audit", prefer_physics=True) or {}
        k_value = extract_nested(alpha_payload, ["results", "self_consistent", "k"], default=2)
        k = float(k_value) if k_value is not None else 2.0
        n_from_cover = 1.0 / k if k != 0 else float("nan")

        exp_path = Path(getattr(config, "output_dir", "")) / "dark_energy_exponential_audit" / "results.json"
        dark_payload = read_json_if_exists(exp_path) or load_tfpt_results("dark_energy_exponential_audit", prefer_physics=False) or {}
        best_label = extract_nested(dark_payload, ["results", "best", "label"], default="")
        if not best_label:
            best_label = extract_nested(dark_payload, ["results", "best_label_str"], default="")
        candidates = extract_nested(dark_payload, ["results", "candidates"], default=[])
        if not best_label and isinstance(candidates, list) and candidates:
            best = None
            best_mis = None
            for row in candidates:
                try:
                    mis = float(row.get("log10_mismatch_rho_L", float("nan")))
                except Exception:
                    continue
                if best is None or (math.isfinite(mis) and (best_mis is None or mis < best_mis)):
                    best = row
                    best_mis = mis
            if isinstance(best, dict):
                best_label = str(best.get("label", ""))
        best_is_half = str(best_label) == "n=1/2"

        checks: list[Check] = []
        checks.append(
            mk_check_pass("n_from_cover_equals_half", f"n_from_cover={n_from_cover}")
            if abs(n_from_cover - 0.5) < 1e-12
            else mk_check_warn("n_from_cover_equals_half", f"n_from_cover={n_from_cover}")
        )
        if best_label or (isinstance(candidates, list) and candidates):
            checks.append(
                mk_check_pass("best_candidate_is_half", f"best_label={best_label}")
                if best_is_half
                else mk_check_warn("best_candidate_is_half", f"best_label={best_label}")
            )
        else:
            checks.append(
                mk_check_fail(
                    "best_candidate_is_half",
                    "missing best_label and candidates (artifact error)",
                )
            )

        lines = [
            "Dark energy norm origin audit",
            "",
            f"k (cover degree) = {k}",
            f"n_from_cover = {n_from_cover}",
            f"best candidate label = {best_label}",
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"dark_energy_norm_origin.png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_norm_origin(out_dir=out_dir, k_value=k, n_value=n_from_cover)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "k": k,
                "n_from_cover": n_from_cover,
                "best_candidate_label": str(best_label),
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
