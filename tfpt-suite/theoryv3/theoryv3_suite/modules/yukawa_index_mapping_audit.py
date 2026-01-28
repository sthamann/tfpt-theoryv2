from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from tfpt_suite.reference_ledger import get_dataset
from theoryv3_suite.utils import REL_ERROR_DEFAULT, coerce_float, ensure_ascii, load_tfpt_results


DENOMINATORS = (2, 3, 4, 5, 6, 8, 10, 12, 16, 24)
MAX_NUMERATOR = 200
MAX_COEFF = 12
REL_ERR_TOL = float(REL_ERROR_DEFAULT)
RATIO_DATASET_MAP = {
    "m_mu_over_m_e": "mass_ratio_mu_over_e_pdg",
    "m_tau_over_m_mu": "mass_ratio_tau_over_mu_pdg",
    "m_s_over_m_d": "mass_ratio_ms_over_md",
    "m_b_over_m_s": "mass_ratio_mb_over_ms",
    "m_c_over_m_u": "mass_ratio_mc_over_mu_quark",
    "m_t_over_m_c": "mass_ratio_mt_over_mc",
}


@dataclass(frozen=True)
class IndexBasis:
    name: str
    value: mp.mpf


def _best_rational(value: float) -> float:
    best = 0.0
    best_err = float("inf")
    for d in DENOMINATORS:
        for n in range(1, MAX_NUMERATOR + 1):
            cand = float(n) / float(d)
            err = abs(float(value) - cand)
            if err < best_err:
                best = cand
                best_err = err
    return best


def _load_index_basis() -> list[IndexBasis]:
    data_path = Path(__file__).resolve().parents[3] / "tfpt_suite" / "data" / "microscopic_action_tfpt_v25.json"
    payload = json.loads(data_path.read_text(encoding="utf-8"))
    sm = payload.get("fields", {}).get("sm", {})
    fermions = sm.get("fermions_left_handed_weyl", []) if isinstance(sm, dict) else []
    basis: list[IndexBasis] = []
    for f in fermions:
        try:
            name = str(f.get("name"))
            reps = f.get("reps", {})
            dims = f.get("dims", {})
            y = mp.mpf(str(reps.get("U1Y")))
            dim_su3 = mp.mpf(str(dims.get("SU3c")))
            dim_su2 = mp.mpf(str(dims.get("SU2L")))
            if not (mp.isfinite(y) and mp.isfinite(dim_su3) and mp.isfinite(dim_su2)):
                continue
            # Charge-squared sum proxy: I = dim(SU3) * dim(SU2) * Y^2
            index_val = dim_su3 * dim_su2 * (y**2)
            if index_val == 0:
                continue
            basis.append(IndexBasis(name=name, value=index_val))
        except Exception:
            continue
    return basis


def _best_index_sum(target: float, basis: list[IndexBasis]) -> tuple[dict[str, int], float, float]:
    best_combo: dict[str, int] = {}
    best_value = 0.0
    best_error = float("inf")
    basis_vals = [float(b.value) for b in basis]
    basis_names = [b.name for b in basis]

    # Brute-force coefficient search (bounded, deterministic).
    for a in range(0, MAX_COEFF + 1):
        for b in range(0, MAX_COEFF + 1):
            for c in range(0, MAX_COEFF + 1):
                for d in range(0, MAX_COEFF + 1):
                    for e in range(0, MAX_COEFF + 1):
                        coeffs = [a, b, c, d, e][: len(basis_vals)]
                        val = 0.0
                        for coeff, v in zip(coeffs, basis_vals):
                            if coeff:
                                val += float(coeff) * v
                        err = abs(val - float(target))
                        if err < best_error:
                            best_error = err
                            best_value = val
                            best_combo = {name: coeff for name, coeff in zip(basis_names, coeffs) if coeff}
                        # Early exit if exact match
                        if err == 0:
                            return best_combo, best_value, best_error
    return best_combo, best_value, best_error


def _plot_mapping_errors(*, out_dir: Path, rows: list[dict[str, Any]]) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"yukawa_index_mapping.png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        labels = [r["ratio"] for r in rows]
        values = [abs(float(r["rel_error"])) for r in rows]

        fig, ax = plt.subplots(figsize=(9.0, 3.6))
        ax.bar(labels, values, color="#dd6b20")
        ax.axhline(REL_ERR_TOL, color="black", lw=1.0, ls="--", alpha=0.8, label="rel error tol")
        ax.set_ylabel("abs relative error")
        ax.set_title("Yukawa index mapping (charge-squared sums)")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.tick_params(axis="x", rotation=20)
        ax.legend(loc="best")

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "yukawa_index_mapping.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["yukawa_index_mapping.png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class YukawaIndexMappingAuditModule(TfptModule):
    module_id = "yukawa_index_mapping_audit"
    title = "Yukawa index mapping audit (q_ij -> charge-squared sums)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "mass ratios from mass_spectrum_deriver (preferred)",
                "TFPT invariants (c3, varphi0)",
                "microscopic_action_tfpt_v25.json (charge-squared sums)",
            ],
            outputs=[
                "q_ij rational targets",
                "best charge-squared sum mapping",
                "relative mapping errors",
            ],
            formulas=[
                "q_ij = (c3/varphi0) * ln(ratio)",
                "I_field = dim(SU3) * dim(SU2) * Y^2",
                "q_target ~ sum_i n_i * I_i (integer coefficients)",
            ],
            validation=[
                "relative mapping errors are <= 2% for the selected ratios",
            ],
            determinism="Deterministic (finite bounded search).",
            question="Can q_ij be mapped to discrete charge-squared index sums from the microscopic action?",
            objective=[
                "Tie exponent indices to explicit charge-squared sums (no continuous tuning).",
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
            ratios = {"m_mu_over_m_e": 197.84536441898265}

        basis = _load_index_basis()
        rows: list[dict[str, Any]] = []
        for ratio_name, ratio_val in ratios.items():
            if not math.isfinite(ratio_val) or ratio_val <= 0:
                continue
            q_float = (c3 / varphi0) * mp.log(mp.mpf(ratio_val))
            q_target = _best_rational(float(q_float))
            combo, best_value, err = _best_index_sum(float(q_target), basis)
            rel_error = float(err / q_target) if q_target != 0 else float("nan")
            rows.append(
                {
                    "ratio": ratio_name,
                    "value": ratio_val,
                    "q_float": float(q_float),
                    "q_target": str(q_target),
                    "combo": combo,
                    "combo_value": float(best_value),
                    "abs_error": float(err),
                    "rel_error": rel_error,
                }
            )

        rows.sort(key=lambda r: abs(float(r["rel_error"])))
        ok = all(abs(float(r["rel_error"])) <= REL_ERR_TOL for r in rows)

        checks: list[Check] = [
            mk_check_pass("mapping_within_tol", f"max_rel_error={max(abs(float(r['rel_error'])) for r in rows):.4g}")
            if ok and rows
            else mk_check_warn("mapping_within_tol", "some ratios exceed rel error tolerance")
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
            "Yukawa index mapping audit",
            "",
            f"c3/varphi0 = {c3/varphi0}",
            f"basis (I = dim(SU3)*dim(SU2)*Y^2): {[f'{b.name}={b.value}' for b in basis]}",
            "",
            "Mappings:",
        ]
        for row in rows:
            lines.append(
                f"- {row['ratio']}: q_target={row['q_target']}, combo={row['combo']}, rel_error={row['rel_error']:.4g}"
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
        plot: dict[str, str | None] = {"yukawa_index_mapping.png": None}
        if getattr(config, "plot", True):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_mapping_errors(out_dir=out_dir, rows=rows)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "ratios": ratios,
                "basis": {b.name: str(b.value) for b in basis},
                "rows": rows,
                "ratio_policies": ratio_policies,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
