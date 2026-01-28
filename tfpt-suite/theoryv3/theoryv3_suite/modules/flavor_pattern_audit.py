from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from tfpt_suite.modules.mobius_cusp_classification import _derive_cusp_set_from_su5_hypercharge
from tfpt_suite.reference_ledger import get_dataset
from theoryv3_suite.utils import Z_TOL_DEFAULT, coerce_float, ensure_ascii, load_tfpt_results, safe_log10


PMNS_ORDERING = "normal_ordering"


def _mobius_ratio(y: mp.mpf, delta: mp.mpf) -> mp.mpf:
    return ((y + delta) / (y - delta)) ** 2


def _plot_flavor_anchors(
    *,
    out_dir: Path,
    lambda_pred: float,
    lambda_ref: float,
    sin2_pred: float,
    sin2_ref: float,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"flavor_anchors_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels = ["lambda", "sin2_theta13"]
        pred = [lambda_pred, sin2_pred]
        ref = [lambda_ref, sin2_ref]

        fig, ax = plt.subplots(figsize=(7.5, 3.6))
        x = range(len(labels))
        ax.bar([i - 0.15 for i in x], pred, width=0.3, label="pred", color="#2f855a")
        ax.bar([i + 0.15 for i in x], ref, width=0.3, label="ref", color="#c05621")
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels)
        ax.set_ylabel("value")
        ax.set_title("Flavor anchors from varphi0")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.legend(loc="best")

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "flavor_anchors.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["flavor_anchors_png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


def _plot_mobius_ratios(
    *,
    out_dir: Path,
    ratios: dict[str, float],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"mobius_ratios_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels = list(ratios.keys())
        values = [safe_log10(float(ratios[k])) for k in labels]

        fig, ax = plt.subplots(figsize=(8.5, 3.6))
        ax.bar(labels, values, color="#3182ce")
        ax.set_ylabel("log10 ratio")
        ax.set_title("Mobius mass ratios (delta_star)")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.tick_params(axis="x", rotation=20)

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "mobius_ratios.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["mobius_ratios_png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class FlavorPatternAuditModule(TfptModule):
    module_id = "flavor_pattern_audit"
    title = "Flavor pattern audit (lambda, delta_star, delta_cp, PMNS theta13)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT constants (varphi0, delta_star)",
                "CKM/PMNS reference tables",
                "Mobius cusp set from SU(5) hypercharge",
            ],
            outputs=[
                "lambda (Cabibbo) from varphi0",
                "delta_star and delta_cp",
                "PMNS sin^2 theta13 from varphi0",
                "Mobius mass ratios for quark/lepton hierarchy",
            ],
            formulas=[
                "lambda = sqrt(varphi0) * (1 - varphi0/2)",
                "delta_star = 3/5 + varphi0/6",
                "delta_cp = pi * (1 - delta_star)",
                "sin2_theta13 = varphi0 * exp(-5/6)",
                "M_y(delta) = (y + delta)/(y - delta)",
            ],
            validation=[
                "lambda within 2 sigma of reference",
                "sin2_theta13 within 2 sigma of NuFIT reference",
                "cusp set equals {1, 1/3, 2/3}",
            ],
            determinism="Deterministic given constants and reference tables.",
            question="Do the compact flavor formulas reproduce the anchor observables?",
            objective=[
                "Make the flavor anchor formulas explicit and testable.",
                "Expose Mobius ratio predictions in one place.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        varphi0 = mp.mpf(c.varphi0)
        delta_star = mp.mpf(c.delta_star)

        lambda_pred = mp.sqrt(varphi0) * (mp.mpf(1) - varphi0 / 2)
        delta_cp_rad = mp.pi * (mp.mpf(1) - delta_star)
        delta_cp_deg = float(delta_cp_rad * 180 / mp.pi)

        sin2_theta13 = varphi0 * mp.e ** (-(mp.mpf(5) / 6))
        theta13_deg = float(mp.asin(mp.sqrt(sin2_theta13)) * 180 / mp.pi)

        # References
        cab_ref = get_dataset("cabibbo_lambda_pdg2024")
        lambda_ref = mp.mpf(str(cab_ref["value"]))
        lambda_sigma = mp.mpf(str(cab_ref["sigma"]))
        pmns_ref = get_dataset("pmns_sin2_theta13_nufit53_no")
        sin2_ref = mp.mpf(str(pmns_ref["value"]))
        sin2_sigma = mp.mpf(str(pmns_ref["sigma"]))

        z_lambda = (lambda_pred - lambda_ref) / lambda_sigma if lambda_sigma != 0 else mp.mpf("nan")
        z_sin2 = (sin2_theta13 - sin2_ref) / sin2_sigma if sin2_sigma != 0 else mp.mpf("nan")

        # Mobius ratios
        ratios: dict[str, float] = {}
        payload = load_tfpt_results("mass_spectrum_deriver", prefer_physics=True)
        if payload:
            ratios_raw = payload.get("results", {}).get("predicted_from_delta_star", {}).get("ratios", {})
            for key, value in ratios_raw.items():
                ratios[key] = coerce_float(value)
        if not ratios:
            y_1 = mp.mpf(1)
            y_23 = mp.mpf(2) / 3
            ratios = {
                "m_s_over_m_d": float(_mobius_ratio(y_1, delta_star)),
                "m_b_over_m_s": float(_mobius_ratio(y_1, delta_star) * (1 + delta_star) ** 2),
                "m_c_over_m_u": float(_mobius_ratio(y_23, delta_star)),
                "m_t_over_m_c": float(((y_23) / (y_23 - delta_star)) ** 2),
            }

        cusp_set = _derive_cusp_set_from_su5_hypercharge()
        cusp_expected = {Fraction(1, 1), Fraction(1, 3), Fraction(2, 3)}

        checks: list[Check] = []
        checks.append(
            mk_check_pass("lambda_within_2sigma", f"z={z_lambda}")
            if mp.isfinite(z_lambda) and abs(z_lambda) <= Z_TOL_DEFAULT
            else mk_check_warn("lambda_within_2sigma", f"z={z_lambda}")
        )
        checks.append(
            mk_check_pass("sin2_theta13_within_2sigma", f"z={z_sin2}")
            if mp.isfinite(z_sin2) and abs(z_sin2) <= Z_TOL_DEFAULT
            else mk_check_warn("sin2_theta13_within_2sigma", f"z={z_sin2}")
        )
        checks.append(
            mk_check_pass("cusp_set_matches", f"cusp_set={sorted(cusp_set)}")
            if cusp_set == cusp_expected
            else mk_check_warn("cusp_set_matches", f"cusp_set={sorted(cusp_set)} (expected {sorted(cusp_expected)})")
        )

        lines = [
            "Flavor pattern audit",
            "",
            f"lambda_pred = {lambda_pred} (ref={lambda_ref}, z={z_lambda})",
            f"delta_star = {delta_star}",
            f"delta_cp_rad = {delta_cp_rad} (deg={delta_cp_deg})",
            f"sin2_theta13 = {sin2_theta13} (ref={sin2_ref}, z={z_sin2})",
            f"theta13_deg = {theta13_deg}",
            "",
            f"cusp_set = {sorted(cusp_set)}",
            "",
            "Mobius ratios:",
            *[f"- {k}: {v}" for k, v in ratios.items()],
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"flavor_anchors_png": None, "mobius_ratios_png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot1, warn1 = _plot_flavor_anchors(
                out_dir=out_dir,
                lambda_pred=float(lambda_pred),
                lambda_ref=float(lambda_ref),
                sin2_pred=float(sin2_theta13),
                sin2_ref=float(sin2_ref),
            )
            plot2, warn2 = _plot_mobius_ratios(out_dir=out_dir, ratios=ratios)
            plot.update(plot1)
            plot.update(plot2)
            warnings.extend(warn1 + warn2)

        return ModuleResult(
            results={
                "references": {"cabibbo_lambda": cab_ref, "pmns_sin2_theta13": pmns_ref},
                "lambda_pred": str(lambda_pred),
                "lambda_ref": str(lambda_ref),
                "lambda_z": str(z_lambda),
                "delta_star": str(delta_star),
                "delta_cp_rad": str(delta_cp_rad),
                "delta_cp_deg": delta_cp_deg,
                "sin2_theta13": str(sin2_theta13),
                "sin2_theta13_ref": str(sin2_ref),
                "sin2_theta13_z": str(z_sin2),
                "theta13_deg": theta13_deg,
                "cusp_set": [str(x) for x in sorted(cusp_set)],
                "mobius_ratios": ratios,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
