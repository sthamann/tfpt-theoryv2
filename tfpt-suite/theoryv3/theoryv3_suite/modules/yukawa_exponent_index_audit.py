from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from tfpt_suite.reference_ledger import get_dataset
from theoryv3_suite.utils import REL_ERROR_DEFAULT, coerce_float, ensure_ascii, load_tfpt_results, safe_log10


DENOMINATORS = (2, 3, 4, 5, 6, 8, 10, 12, 16, 24)
MAX_NUMERATOR = 200
REL_ERR_TOL = float(REL_ERROR_DEFAULT)
RATIO_DATASET_MAP = {
    "m_mu_over_m_e": "mass_ratio_mu_over_e_pdg",
    "m_tau_over_m_mu": "mass_ratio_tau_over_mu_pdg",
    "m_s_over_m_d": "mass_ratio_ms_over_md",
    "m_b_over_m_s": "mass_ratio_mb_over_ms",
    "m_c_over_m_u": "mass_ratio_mc_over_mu_quark",
    "m_t_over_m_c": "mass_ratio_mt_over_mc",
}


def _best_rational(value: float) -> Fraction:
    assert DENOMINATORS, "denominator set must not be empty"
    best = Fraction(0, 1)
    best_err = float("inf")
    for d in DENOMINATORS:
        for n in range(1, MAX_NUMERATOR + 1):
            cand = Fraction(n, d)
            err = abs(float(value) - float(cand))
            if err < best_err:
                best = cand
                best_err = err
    return best


def _plot_q_errors(
    *, out_dir: Path, rows: list[dict[str, Any]]
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"yukawa_q_errors_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels = [r["ratio"] for r in rows]
        values = [abs(float(r["rel_error"])) for r in rows]

        fig, ax = plt.subplots(figsize=(9.0, 3.6))
        ax.bar(labels, values, color="#dd6b20")
        ax.axhline(REL_ERR_TOL, color="black", lw=1.0, ls="--", alpha=0.8, label="rel error tol")
        ax.set_ylabel("abs relative error")
        ax.set_title("Yukawa exponent index rationalization")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.tick_params(axis="x", rotation=20)
        ax.legend(loc="best")

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "yukawa_q_errors.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["yukawa_q_errors_png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class YukawaExponentIndexAuditModule(TfptModule):
    module_id = "yukawa_exponent_index_audit"
    title = "Yukawa exponent index audit (q_ij rationalization)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "mass ratios from mass_spectrum_deriver (preferred)",
                "TFPT invariants (c3, varphi0)",
            ],
            outputs=[
                "q_ij = (c3/varphi0) * ln(m_i/m_j)",
                "best rational approximation for q_ij",
                "reconstruction error for each ratio",
            ],
            formulas=[
                "q_ij = (c3/varphi0) * ln(ratio)",
                "ratio_recon = exp(q_rational * varphi0 / c3)",
            ],
            validation=[
                "relative reconstruction errors are below 2% for the selected ratios",
            ],
            determinism="Deterministic (finite rational scan).",
            question="Do the Mobius-derived ratios align with rational exponent indices?",
            objective=[
                "Quantize the exponent indices with a fixed denominator grammar.",
                "Report which ratios are most stable under the rational map.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        varphi0 = mp.mpf(c.varphi0)
        c3 = mp.mpf(c.c3)

        payload = load_tfpt_results("mass_spectrum_deriver", prefer_physics=True)
        ratios: dict[str, float] = {}
        if payload:
            ratios_raw = payload.get("results", {}).get("predicted_from_delta_star", {}).get("ratios", {})
            for key, value in ratios_raw.items():
                ratios[key] = coerce_float(value)

        if not ratios:
            # Fallback to a minimal set (avoid empty run).
            ratios = {"m_mu_over_m_e": 197.84536441898265}

        rows: list[dict[str, Any]] = []
        for ratio_name, ratio_val in ratios.items():
            if not math.isfinite(ratio_val) or ratio_val <= 0:
                continue
            q_float = (c3 / varphi0) * mp.log(mp.mpf(ratio_val))
            q_best = _best_rational(float(q_float))
            ratio_recon = mp.e ** (mp.mpf(q_best.numerator) / mp.mpf(q_best.denominator) * (varphi0 / c3))
            rel_error = float((ratio_recon / mp.mpf(ratio_val)) - 1)
            rows.append(
                {
                    "ratio": ratio_name,
                    "value": ratio_val,
                    "q_float": float(q_float),
                    "q_rational": f"{q_best.numerator}/{q_best.denominator}",
                    "ratio_recon": float(ratio_recon),
                    "rel_error": rel_error,
                    "log10_ratio": safe_log10(float(ratio_val)),
                }
            )

        rows.sort(key=lambda r: abs(float(r["rel_error"])))
        ok = all(abs(float(r["rel_error"])) <= REL_ERR_TOL for r in rows)

        checks: list[Check] = [
            mk_check_pass("rationalization_within_tol", f"max_rel_error={max(abs(float(r['rel_error'])) for r in rows):.4g}")
            if ok and rows
            else mk_check_warn("rationalization_within_tol", "some ratios exceed rel error tolerance")
        ]

        ratio_policies: dict[str, dict[str, str]] = {}
        for row in rows:
            ratio = row["ratio"]
            dataset_id = RATIO_DATASET_MAP.get(ratio, "")
            if not dataset_id:
                continue
            ds = get_dataset(dataset_id)
            status = str(ds.get("status", "direct"))
            ratio_policies[ratio] = {
                "dataset_id": dataset_id,
                "status": status,
                "scheme": str(ds.get("scheme", "n/a")),
                "scale": str(ds.get("scale", "n/a")),
            }
            if status == "needs_rg":
                checks.append(mk_check_warn(f"ratio_needs_rg:{ratio}", f"scheme={ds.get('scheme')}, scale={ds.get('scale')}"))
            else:
                checks.append(mk_check_pass(f"ratio_direct:{ratio}", "scale-independent reference"))

        lines = [
            "Yukawa exponent index audit",
            "",
            f"c3/varphi0 = {c3/varphi0}",
            f"denominators = {DENOMINATORS}, max_numerator = {MAX_NUMERATOR}",
            "",
            "Ratios:",
        ]
        for row in rows:
            lines.append(
                f"- {row['ratio']}: q={row['q_float']:.6f}, q_rational={row['q_rational']}, rel_error={row['rel_error']:.4g}"
            )
        lines += [
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"yukawa_q_errors_png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_q_errors(out_dir=out_dir, rows=rows)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "ratios": ratios,
                "rows": rows,
                "ratio_policies": ratio_policies,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
