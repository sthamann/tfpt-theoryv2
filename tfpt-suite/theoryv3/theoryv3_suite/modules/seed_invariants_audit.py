from __future__ import annotations

from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from theoryv3_suite.utils import coerce_mpf, ensure_ascii, load_tfpt_results


TOL_COMPARE = mp.mpf("1e-24")


def _plot_seed_invariants(*, out_dir: Path, values: dict[str, mp.mpf]) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"seed_invariants_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        keys = ["c3", "varphi0_tree", "delta_top", "varphi0", "beta_rad"]
        labels = ["c3", "varphi0_tree", "delta_top", "varphi0", "beta_rad"]
        vals = [float(values[k]) for k in keys]

        fig, ax = plt.subplots(figsize=(8.5, 4.0))
        ax.bar(labels, vals, color="#2b6cb0")
        ax.set_yscale("log")
        ax.set_ylabel("value (log scale)")
        ax.set_title("Seed invariants from pi")
        ax.grid(True, axis="y", ls=":", alpha=0.4)

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "seed_invariants.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["seed_invariants_png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class SeedInvariantsAuditModule(TfptModule):
    module_id = "seed_invariants_audit"
    title = "Seed invariants audit (pi -> c3, varphi0, beta_rad)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["pi"],
            outputs=["c3", "varphi0_tree", "delta_top", "varphi0", "beta_rad"],
            formulas=[
                "c3 = 1/(8*pi)",
                "varphi0_tree = 1/(6*pi)",
                "delta_top = 3/(256*pi^4)",
                "varphi0 = varphi0_tree + delta_top",
                "beta_rad = varphi0/(4*pi)",
            ],
            validation=[
                "internal identities are satisfied",
                "optional match to core_invariants output if present",
            ],
            determinism="Deterministic (pure algebra).",
            question="Does pi alone fix the core invariants used across TFPT?",
            objective=[
                "Make the pi-seed dependency explicit and auditable.",
                "Provide a compact reference table for downstream audits.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        values = {
            "pi": mp.mpf(c.pi),
            "c3": mp.mpf(c.c3),
            "varphi0_tree": mp.mpf(c.varphi0_tree),
            "delta_top": mp.mpf(c.delta_top),
            "varphi0": mp.mpf(c.varphi0),
            "beta_rad": mp.mpf(c.beta_rad),
        }

        # Optional comparison to existing suite output.
        core_payload = load_tfpt_results("core_invariants", prefer_physics=True)
        compare: dict[str, Any] = {}
        checks: list[Check] = []
        if core_payload:
            core_constants = core_payload.get("results", {}).get("constants", {})
            for key in ("c3", "varphi0_tree", "delta_top", "varphi0", "beta_rad"):
                ref = coerce_mpf(core_constants.get(key))
                diff = values[key] - ref
                ok = bool(mp.isfinite(ref) and abs(diff) <= TOL_COMPARE)
                compare[key] = {"ref": ref, "diff": diff}
                checks.append(
                    mk_check_pass(
                        f"match_{key}",
                        f"{key} matches core_invariants (diff={diff})",
                    )
                    if ok
                    else mk_check_warn(
                        f"match_{key}",
                        f"{key} mismatch vs core_invariants (diff={diff})",
                    )
                )
        else:
            checks.append(mk_check_warn("core_invariants_missing", "core_invariants results.json not found (comparison skipped)"))

        # Internal identities.
        checks.append(
            mk_check_pass("varphi0_identity", f"varphi0 = varphi0_tree + delta_top ({values['varphi0']})")
            if abs(values["varphi0"] - (values["varphi0_tree"] + values["delta_top"])) <= TOL_COMPARE
            else mk_check_warn("varphi0_identity", "varphi0 identity mismatch")
        )
        checks.append(
            mk_check_pass("beta_rad_identity", f"beta_rad = varphi0/(4*pi) ({values['beta_rad']})")
            if abs(values["beta_rad"] - (values["varphi0"] / (mp.mpf(4) * values["pi"]))) <= TOL_COMPARE
            else mk_check_warn("beta_rad_identity", "beta_rad identity mismatch")
        )

        lines = [
            "Seed invariants audit (pi -> c3, varphi0, beta_rad)",
            "",
            f"pi = {values['pi']}",
            f"c3 = {values['c3']}",
            f"varphi0_tree = {values['varphi0_tree']}",
            f"delta_top = {values['delta_top']}",
            f"varphi0 = {values['varphi0']}",
            f"beta_rad = {values['beta_rad']}",
            "",
        ]
        if compare:
            lines.append("Comparison to core_invariants (if present):")
            for key, row in compare.items():
                lines.append(f"- {key}: ref={row['ref']}, diff={row['diff']}")
            lines.append("")

        lines.append("Checks:")
        lines.extend(
            [
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ]
        )

        warnings: list[str] = []
        plot: dict[str, str | None] = {"seed_invariants_png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_seed_invariants(out_dir=out_dir, values=values)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "values": {k: str(v) for k, v in values.items()},
                "comparison": {k: {"ref": str(v["ref"]), "diff": str(v["diff"])} for k, v in compare.items()},
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
