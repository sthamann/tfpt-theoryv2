from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from tfpt_suite.matching import match_alpha3_quark_threshold_2loop, match_gauge
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _plot_matching_metamorphic(
    *,
    out_dir: Path,
    alpha3_rel_errors: list[float],
    gauge_rel_errors: list[float],
    alpha3_tol_rel: float,
    gauge_tol_rel: float,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"matching_metamorphic_errors_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        def to_log10(xs: list[float]) -> list[float]:
            out: list[float] = []
            for x in xs:
                v = float(x)
                if not np.isfinite(v):
                    continue
                out.append(float(np.log10(max(1e-20, v))))
            return out

        a = to_log10(alpha3_rel_errors)
        g = to_log10(gauge_rel_errors)

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(11.5, 4.2))

        if a:
            ax1.hist(a, bins=60, color="#1f77b4", alpha=0.85)
        ax1.axvline(float(np.log10(alpha3_tol_rel)), color="black", lw=1.0, ls="--", alpha=0.8, label="tolerance")
        ax1.set_title(r"$\alpha_s$ threshold roundtrip errors")
        ax1.set_xlabel(r"$\log_{10}(\mathrm{rel\_err})$")
        ax1.set_ylabel("count")
        ax1.grid(True, ls=":", alpha=0.35)
        ax1.legend(loc="best")

        if g:
            ax2.hist(g, bins=60, color="#ff7f0e", alpha=0.85)
        ax2.axvline(float(np.log10(gauge_tol_rel)), color="black", lw=1.0, ls="--", alpha=0.8, label="tolerance")
        ax2.set_title(r"Gauge finite matching roundtrip errors")
        ax2.set_xlabel(r"$\log_{10}(\mathrm{rel\_err})$")
        ax2.set_ylabel("count")
        ax2.grid(True, ls=":", alpha=0.35)
        ax2.legend(loc="best")

        fig.tight_layout()
        path = out_dir / "matching_metamorphic_errors.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["matching_metamorphic_errors_png"] = str(path)

    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class _WorstCase:
    label: str
    inputs: dict[str, float]
    outputs: dict[str, float]
    abs_error: float
    rel_error: float


def _relerr(a: float, b: float) -> float:
    denom = max(1e-30, abs(a), abs(b))
    return abs(a - b) / denom


class MatchingMetamorphicAuditModule(TfptModule):
    """
    Metamorphic / “property-style” checks for the matching layer.

    This does not add new physics. It hardens the suite by checking that
    small matching primitives satisfy basic invertibility / stability properties.
    """

    module_id = "ux_matching_metamorphic_audit"
    title = "Unconventional: matching metamorphic audit (invertibility + stability)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "`tfpt_suite/matching.py` primitives",
                "deterministic RNG seed (SuiteConfig.seed)",
            ],
            outputs=[
                "max abs/rel errors for metamorphic properties",
                "worst-case counterexamples (if any)",
            ],
            formulas=[
                "Metamorphic property examples:",
                "- α3 threshold: down(up(α3)) ≈ α3 (within truncation order)",
                "- finite gauge matching: down(up(g_i, δα)) = g_i (should be exact up to float rounding)",
            ],
            validation=[
                "report max observed errors and fail explicitly if they exceed conservative tolerances",
            ],
            determinism="Deterministic given config seed + precision.",
        )

    def run(self, config) -> ModuleResult:
        # Keep the runtime small and deterministic.
        n_alpha3 = 500
        n_gauge = 400

        worst_alpha3: _WorstCase | None = None
        max_abs_alpha3 = 0.0
        max_rel_alpha3 = 0.0
        alpha3_rel_errors: list[float] = []

        for _ in range(n_alpha3):
            a0 = float(0.01 + 0.29 * random.random())  # [0.01, 0.30]
            a_up = float(match_alpha3_quark_threshold_2loop(alpha3=a0, direction="up"))
            a_back = float(match_alpha3_quark_threshold_2loop(alpha3=a_up, direction="down"))
            abs_err = float(abs(a_back - a0))
            rel_err = float(_relerr(a_back, a0))
            alpha3_rel_errors.append(float(rel_err))
            if abs_err > max_abs_alpha3 or worst_alpha3 is None:
                max_abs_alpha3 = abs_err
                max_rel_alpha3 = max(max_rel_alpha3, rel_err)
                worst_alpha3 = _WorstCase(
                    label="alpha3_threshold_up_then_down",
                    inputs={"alpha3": a0},
                    outputs={"alpha3_up": a_up, "alpha3_back": a_back},
                    abs_error=abs_err,
                    rel_error=rel_err,
                )
            else:
                max_rel_alpha3 = max(max_rel_alpha3, rel_err)

        worst_gauge: _WorstCase | None = None
        max_abs_gauge = 0.0
        max_rel_gauge = 0.0
        gauge_rel_errors: list[float] = []

        for _ in range(n_gauge):
            # Sample around physically sensible magnitudes, but this is a math property test, not a fit.
            gY = float(0.2 + 0.4 * random.random())  # ~[0.2, 0.6]
            g2 = float(0.4 + 0.6 * random.random())  # ~[0.4, 1.0]
            g3 = float(0.9 + 0.7 * random.random())  # ~[0.9, 1.6]

            # Finite α-shifts (small, to keep α >= 0 robustly).
            dY = float((2.0 * random.random() - 1.0) * 1e-5)
            d2 = float((2.0 * random.random() - 1.0) * 1e-5)
            d3 = float((2.0 * random.random() - 1.0) * 1e-5)

            below = {"gY": gY, "g2": g2, "g3": g3}
            above, out_up = match_gauge(
                threshold_id="ux_finite_delta_smoke",
                mu_thr_GeV=173.0,
                direction="up",
                couplings_below=below,
                scheme="MSbar",
                loop_order=2,
                finite_delta_alpha={"alphaY": dY, "alpha2": d2, "alpha3": d3},
            )
            back, out_down = match_gauge(
                threshold_id="ux_finite_delta_smoke",
                mu_thr_GeV=173.0,
                direction="down",
                couplings_below=above,
                scheme="MSbar",
                loop_order=2,
                finite_delta_alpha={"alphaY": dY, "alpha2": d2, "alpha3": d3},
            )

            # Aggregate error over the three couplings (L-infinity).
            abs_err = max(abs(float(back["gY"]) - gY), abs(float(back["g2"]) - g2), abs(float(back["g3"]) - g3))
            rel_err = max(_relerr(float(back["gY"]), gY), _relerr(float(back["g2"]), g2), _relerr(float(back["g3"]), g3))
            gauge_rel_errors.append(float(rel_err))

            if abs_err > max_abs_gauge or worst_gauge is None:
                max_abs_gauge = float(abs_err)
                max_rel_gauge = max(max_rel_gauge, float(rel_err))
                worst_gauge = _WorstCase(
                    label="match_gauge_finite_delta_up_then_down",
                    inputs={"gY": gY, "g2": g2, "g3": g3, "dalphaY": dY, "dalpha2": d2, "dalpha3": d3},
                    outputs={"gY_up": float(above["gY"]), "g2_up": float(above["g2"]), "g3_up": float(above["g3"]), "gY_back": float(back["gY"]), "g2_back": float(back["g2"]), "g3_back": float(back["g3"])},
                    abs_error=float(abs_err),
                    rel_error=float(rel_err),
                )
            else:
                max_rel_gauge = max(max_rel_gauge, float(rel_err))

            # Sanity: both calls should declare matching active.
            assert out_up.matching_active and out_down.matching_active

        # Conservative tolerances:
        #
        # - α3 decoupling is a truncated series; exact inversion is not expected.
        # - finite matching deltas in match_gauge are applied explicitly and should be nearly exact.
        alpha3_tol_rel = 1e-5
        gauge_tol_rel = 1e-12

        checks = [
            Check(
                check_id="alpha3_threshold_up_down_near_identity",
                passed=bool(max_rel_alpha3 <= alpha3_tol_rel),
                detail=f"max_rel_err={max_rel_alpha3:.3e} (tol={alpha3_tol_rel:.1e}); max_abs_err={max_abs_alpha3:.3e}; samples={n_alpha3}",
            ),
            Check(
                check_id="match_gauge_finite_delta_up_down_identity",
                passed=bool(max_rel_gauge <= gauge_tol_rel),
                detail=f"max_rel_err={max_rel_gauge:.3e} (tol={gauge_tol_rel:.1e}); max_abs_err={max_abs_gauge:.3e}; samples={n_gauge}",
            ),
        ]

        report_lines: list[str] = []
        report_lines += [
            "Unconventional: matching metamorphic audit",
            "",
            "Goal:",
            "- Harden the matching layer by verifying metamorphic properties (invertibility/stability).",
            "- This is a math/engineering audit; it does not add new physics content.",
            "",
            f"Samples: alpha3={n_alpha3}, gauge={n_gauge}",
            "",
            "Results:",
            f"- alpha3 threshold: max_abs_err={max_abs_alpha3:.6e}, max_rel_err={max_rel_alpha3:.6e}",
            f"- gauge finite matching: max_abs_err={max_abs_gauge:.6e}, max_rel_err={max_rel_gauge:.6e}",
            "",
            "Worst cases:",
            f"- alpha3: {worst_alpha3}" if worst_alpha3 is not None else "- alpha3: n/a",
            f"- gauge:  {worst_gauge}" if worst_gauge is not None else "- gauge:  n/a",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- `match_alpha3_quark_threshold_2loop` is a truncated series: exact invertibility is not expected.",
            "- `match_gauge` with `finite_delta_alpha` is explicit by construction: up+down should invert within float rounding.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"matching_metamorphic_errors_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_matching_metamorphic(
                out_dir=out_dir,
                alpha3_rel_errors=alpha3_rel_errors,
                gauge_rel_errors=gauge_rel_errors,
                alpha3_tol_rel=alpha3_tol_rel,
                gauge_tol_rel=gauge_tol_rel,
            )
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "samples": {"alpha3": n_alpha3, "gauge": n_gauge},
                "alpha3_threshold": {
                    "max_abs_error": max_abs_alpha3,
                    "max_rel_error": max_rel_alpha3,
                    "tolerance_rel": alpha3_tol_rel,
                    "worst_case": worst_alpha3,
                },
                "gauge_finite_matching": {
                    "max_abs_error": max_abs_gauge,
                    "max_rel_error": max_rel_gauge,
                    "tolerance_rel": gauge_tol_rel,
                    "worst_case": worst_gauge,
                },
                "plot": plot,
            },
            checks=checks,
            report="\n".join(report_lines),
            warnings=warnings,
        )

