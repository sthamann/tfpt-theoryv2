from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


GEV_TO_HZ = 2.417989242e23  # E/h
K_MINIMAL_COUPLING = 0.75  # minimal coupling mapping (|b_mu| ≃ (3/4)|A_mu|)


def _plot_torsion_regime_designer(
    *,
    out_dir: Path,
    tight_bound_GeV: float,
    anchor_rows: list[dict[str, Any]],
    candidates: list[dict[str, float | str]],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"torsion_regime_designer_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        # Prepare candidate distributions of log10(c_spin_max).
        def logs(kind: str) -> list[float]:
            xs: list[float] = []
            for r in candidates:
                if str(r.get("coupling_kind", "")) != kind:
                    continue
                cs = float(r.get("c_spin_max_to_saturate_tightest_bound", float("nan")))
                if not (np.isfinite(cs) and cs > 0):
                    continue
                xs.append(float(np.log10(cs)))
            return xs

        log_mpl = logs("Mpl_reduced")
        log_tfpt = logs("TFPT_M")

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12.2, 4.6))

        # Panel 1: histogram of log10(c_spin_max).
        bins = 60
        if log_mpl:
            ax1.hist(log_mpl, bins=bins, alpha=0.7, label="Mpl_reduced", color="#1f77b4")
        if log_tfpt:
            ax1.hist(log_tfpt, bins=bins, alpha=0.7, label="TFPT_M", color="#ff7f0e")
        ax1.axvline(0.0, color="black", lw=1.0, ls="--", alpha=0.8)
        ax1.set_xlabel(r"$\log_{10}(c_{\rm spin,max})$ needed to saturate tightest bound")
        ax1.set_ylabel("count")
        ax1.set_title("How strong can the coupling be? (distribution)")
        ax1.grid(True, ls=":", alpha=0.35)
        ax1.legend(loc="best")

        # Panel 2: anchor scenarios vs tightest bound (S_pred with c_spin=1).
        labels = [str(r.get("label", "?")) for r in anchor_rows]
        S1 = [float(r.get("S_mu_pred_cspin_1_GeV", float("nan"))) for r in anchor_rows]
        x = np.arange(len(labels), dtype=float)
        ax2.bar(x, S1, color="#2ca02c", alpha=0.85, label=r"$|S_\mu|$ (c_spin=1)")
        if np.isfinite(float(tight_bound_GeV)) and float(tight_bound_GeV) > 0:
            ax2.axhline(float(tight_bound_GeV), color="black", lw=1.2, alpha=0.85, label="tightest bound")
        ax2.set_yscale("log")
        ax2.set_xticks(x)
        ax2.set_xticklabels(labels, rotation=25, ha="right")
        ax2.set_ylabel(r"$|S_\mu|$ [GeV] (log)")
        ax2.set_title("Anchor scenarios (toy model)")
        ax2.grid(True, axis="y", ls=":", alpha=0.35)
        ax2.legend(loc="best")

        fig.tight_layout()
        path = out_dir / "torsion_regime_designer.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["torsion_regime_designer_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class BoundSummary:
    label: str
    component: str
    abs_max_GeV: float
    source: str


@dataclass(frozen=True)
class SpinMediumScenario:
    label: str
    number_density_value: float
    number_density_unit: str  # "fm^-3" | "cm^-3"
    polarization: float
    spin_per_particle: float
    coupling_kind: str  # "Mpl_reduced" | "TFPT_M"


def _tfpt_suite_dir() -> Path:
    # .../tfpt-suite/unconventional/tfpt_unconventional/modules/<this_file>
    return Path(__file__).resolve().parents[3]


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _number_density_to_GeV3(value: float, unit: str) -> float:
    v = float(value)
    if v <= 0:
        raise ValueError("number density must be positive")
    u = str(unit).strip()
    if u == "cm^-3":
        # 1 cm = 5.0677307e13 GeV^-1
        cm_in_GeV_inv = 5.0677307e13
        return float(v * (1.0 / cm_in_GeV_inv) ** 3)
    if u == "fm^-3":
        # 1 fm = 5.0677307 GeV^-1
        fm_in_GeV_inv = 5.0677307
        return float(v * (1.0 / fm_in_GeV_inv) ** 3)
    raise ValueError(f"unsupported number density unit: {unit!r}")


def _M_eff_GeV(*, kind: str) -> float:
    kind = str(kind).strip()
    Mpl_reduced_GeV = 2.435e18
    if kind == "Mpl_reduced":
        return float(Mpl_reduced_GeV)
    if kind == "TFPT_M":
        c = TfptConstants.compute()
        return float(float(c.M_over_Mpl) * float(Mpl_reduced_GeV))
    raise ValueError(f"unsupported coupling kind: {kind!r}")


def _S_mu_pred_GeV(
    *,
    n_value: float,
    n_unit: str,
    polarization: float,
    spin_per_particle: float,
    c_spin: float,
    coupling_kind: str,
) -> float:
    n_GeV3 = _number_density_to_GeV3(n_value, n_unit)
    pol = float(polarization)
    spin = float(spin_per_particle)
    cs = float(c_spin)
    if not (0.0 <= pol <= 1.0):
        raise ValueError("polarization must be in [0,1]")
    if spin <= 0:
        raise ValueError("spin_per_particle must be positive")
    if cs <= 0:
        raise ValueError("c_spin must be positive")
    M = _M_eff_GeV(kind=coupling_kind)
    spin_density = float(pol * spin * n_GeV3)
    return float(abs(cs) * abs(spin_density) / (M * M))


class TorsionRegimeDesignerModule(TfptModule):
    """
    Unconventional module (F): propose nontrivial present-day torsion regimes.

    The main suite already supports regime evaluation in `torsion_bounds_mapping` with:
      model = "spin_polarized_medium"

    This module helps *design* candidate regimes and makes the “needed suppression” explicit:
    it computes, for chosen matter scenarios, what dimensionless factor c_spin would saturate
    the tightest vetted SME torsion bounds.
    """

    module_id = "ux_torsion_regime_designer"
    title = "Unconventional: torsion regime designer (spin-polarized medium + bound saturation targets)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "vetted bounds: tfpt_suite/data/torsion_bounds_vetted.json (A_mu component-wise limits)",
                "TFPT scale M from constants.py (for coupling_kind='TFPT_M')",
            ],
            outputs=[
                "tightest bound summary",
                "scenario table (S_pred for c_spin=1, and c_spin needed to saturate bounds)",
                "JSON regime proposals compatible with tfpt_suite/data/torsion_regimes.json",
            ],
            formulas=[
                "toy spin-medium model: |S_mu| ~ c_spin * (pol * spin_per_particle * n) / M_eff^2",
                "compute c_spin_max := bound / S_pred(c_spin=1) for saturation targets",
            ],
            validation=[
                "produce at least one explicit regime proposal (even if conservative)",
                "make all assumptions explicit (toy model; c_spin is a derivation target)",
            ],
            determinism="Deterministic given config seed (random sampling for scenario exploration).",
        )

    def run(self, config) -> ModuleResult:
        suite_dir = _tfpt_suite_dir()
        bounds_path = suite_dir / "tfpt_suite" / "data" / "torsion_bounds_vetted.json"
        bounds_raw = _read_json(bounds_path)
        bounds_list = bounds_raw.get("bounds", []) if isinstance(bounds_raw.get("bounds", []), list) else []

        vetted: list[BoundSummary] = []
        for b in bounds_list:
            if not isinstance(b, dict):
                continue
            if str(b.get("coefficient")) not in {"A_mu", "S_mu"}:
                continue
            if "abs_max" not in b:
                continue
            vetted.append(
                BoundSummary(
                    label=str(b.get("label", "?")),
                    component=str(b.get("component", "?")),
                    abs_max_GeV=float(b.get("abs_max")),
                    source=str(b.get("source", "")),
                )
            )
        vetted.sort(key=lambda x: x.abs_max_GeV)
        tight = vetted[0] if vetted else BoundSummary(label="n/a", component="?", abs_max_GeV=float("nan"), source="")

        # A small set of physically motivated “anchors” (numbers are order-of-magnitude placeholders).
        anchors: list[SpinMediumScenario] = [
            SpinMediumScenario(
                label="laboratory_electron_gas",
                number_density_value=1e23,  # ~solid density electrons, cm^-3
                number_density_unit="cm^-3",
                polarization=0.5,
                spin_per_particle=0.5,
                coupling_kind="Mpl_reduced",
            ),
            SpinMediumScenario(
                label="nuclear_matter_core",
                number_density_value=0.16,  # fm^-3 (nuclear saturation density)
                number_density_unit="fm^-3",
                polarization=0.5,
                spin_per_particle=0.5,
                coupling_kind="TFPT_M",
            ),
            SpinMediumScenario(
                label="extreme_dense_spin_fluid",
                number_density_value=1.0,  # fm^-3 (very dense)
                number_density_unit="fm^-3",
                polarization=1.0,
                spin_per_particle=0.5,
                coupling_kind="TFPT_M",
            ),
        ]

        # Random exploration around the anchors (to find “interesting” regimes close to bounds).
        samples = 1200
        candidates: list[dict[str, float | str]] = []

        for _ in range(samples):
            coupling_kind = random.choice(["Mpl_reduced", "TFPT_M"])
            if coupling_kind == "Mpl_reduced":
                # laboratory-ish number densities, log-uniform 1e18..1e26 cm^-3
                logn = random.uniform(18.0, 26.0)
                n_val = 10.0**logn
                n_unit = "cm^-3"
            else:
                # compact objects: log-uniform 1e-3..1 fm^-3
                logn = random.uniform(-3.0, 0.0)
                n_val = 10.0**logn
                n_unit = "fm^-3"

            pol = random.uniform(0.01, 1.0)
            spin = 0.5
            S1 = _S_mu_pred_GeV(
                n_value=float(n_val),
                n_unit=n_unit,
                polarization=float(pol),
                spin_per_particle=float(spin),
                c_spin=1.0,
                coupling_kind=str(coupling_kind),
            )
            if not math.isfinite(S1) or S1 <= 0:
                continue

            c_spin_max = float(tight.abs_max_GeV / S1) if (math.isfinite(tight.abs_max_GeV) and tight.abs_max_GeV > 0) else float("nan")
            b1 = float(K_MINIMAL_COUPLING * float(S1)) if math.isfinite(float(S1)) else float("nan")
            dnu1 = float(2.0 * b1 * GEV_TO_HZ) if math.isfinite(b1) else float("nan")
            candidates.append(
                {
                    "coupling_kind": coupling_kind,
                    "n_value": float(n_val),
                    "n_unit": n_unit,
                    "polarization": float(pol),
                    "spin_per_particle": float(spin),
                    "S_mu_pred_cspin_1_GeV": float(S1),
                    "delta_nu_Hz_cspin_1": float(dnu1),
                    "c_spin_max_to_saturate_tightest_bound": float(c_spin_max),
                }
            )

        # Find “best” candidates per coupling kind: maximize S_pred but still allow c_spin<=1.
        def best_for(kind: str) -> dict[str, float | str] | None:
            rows = [r for r in candidates if r["coupling_kind"] == kind and float(r["c_spin_max_to_saturate_tightest_bound"]) > 0]
            if not rows:
                return None
            # Want c_spin_max close to 1 (no suppression) and large S.
            def key(r: dict[str, float | str]) -> tuple[float, float]:
                cs = float(r["c_spin_max_to_saturate_tightest_bound"])
                # Prefer cs>=1 (allowed even with c_spin=1), otherwise closer to 1 is better.
                penalty = 0.0 if cs >= 1.0 else (1.0 / max(1e-30, cs))
                return (penalty, -float(r["S_mu_pred_cspin_1_GeV"]))

            rows.sort(key=key)
            return rows[0]

        best_mpl = best_for("Mpl_reduced")
        best_tfpt = best_for("TFPT_M")

        # Build explicit proposals in torsion_regimes.json schema (model="spin_polarized_medium").
        proposals: list[dict[str, Any]] = []
        for idx, row in enumerate([best_mpl, best_tfpt], start=1):
            if row is None:
                continue
            cs_max = float(row["c_spin_max_to_saturate_tightest_bound"])
            # Choose a conservative c_spin: half of cs_max, capped at 1.
            cs_use = float(min(1.0, 0.5 * cs_max))
            proposals.append(
                {
                    "id": f"spin_medium_candidate_{idx}",
                    "label": f"Spin-polarized medium candidate {idx} ({row['coupling_kind']})",
                    "model": "spin_polarized_medium",
                    "compare_to_bounds": True,
                    "number_density": {"value": float(row["n_value"]), "unit": str(row["n_unit"])},
                    "polarization": float(row["polarization"]),
                    "spin_per_particle": float(row["spin_per_particle"]),
                    "c_spin": cs_use,
                    "coupling_scale": {"kind": str(row["coupling_kind"])},
                    "note": (
                        "Generated by ux_torsion_regime_designer. "
                        "Toy model: |S_mu| ~ c_spin*(pol*spin*n)/M_eff^2. "
                        "c_spin is not derived here; it is a target for the microscopic operator derivation. "
                        f"Chosen c_spin={cs_use:.3g} (≈0.5*c_spin_max) to stay safely below the tightest bound."
                    ),
                }
            )

        # Anchor table: compute S_pred and c_spin_max for the named scenarios.
        anchor_rows: list[dict[str, Any]] = []
        for sc in anchors:
            S1 = _S_mu_pred_GeV(
                n_value=sc.number_density_value,
                n_unit=sc.number_density_unit,
                polarization=sc.polarization,
                spin_per_particle=sc.spin_per_particle,
                c_spin=1.0,
                coupling_kind=sc.coupling_kind,
            )
            c_spin_max = float(tight.abs_max_GeV / S1) if (math.isfinite(tight.abs_max_GeV) and tight.abs_max_GeV > 0 and S1 > 0) else float("nan")
            b1 = float(K_MINIMAL_COUPLING * float(S1)) if math.isfinite(float(S1)) else float("nan")
            dnu1 = float(2.0 * b1 * GEV_TO_HZ) if math.isfinite(b1) else float("nan")
            anchor_rows.append(
                {
                    "label": sc.label,
                    "coupling_kind": sc.coupling_kind,
                    "n_value": sc.number_density_value,
                    "n_unit": sc.number_density_unit,
                    "polarization": sc.polarization,
                    "spin_per_particle": sc.spin_per_particle,
                    "S_mu_pred_cspin_1_GeV": S1,
                    "delta_nu_Hz_cspin_1": dnu1,
                    "c_spin_max_to_saturate_tightest_bound": c_spin_max,
                }
            )

        checks = [
            Check(
                check_id="vetted_bounds_loaded",
                passed=bool(vetted),
                detail=f"bounds={len(vetted)}; tightest={tight.label}:{tight.abs_max_GeV:.3e} GeV",
            ),
            Check(
                check_id="proposals_generated",
                passed=bool(len(proposals) >= 1),
                detail=f"proposals={len(proposals)} (compatible with torsion_regimes.json schema)",
            ),
        ]

        lines: list[str] = []
        lines += [
            "Unconventional: torsion regime designer",
            "",
            "Purpose:",
            "- Provide concrete, nontrivial candidate regimes for present-day torsion falsifiability.",
            "- Make explicit what suppression/normalization factor c_spin would be required to saturate vetted bounds.",
            "",
            f"Bounds file: {bounds_path}",
            f"Tightest component-wise bound: {tight.label} (component={tight.component}) <= {tight.abs_max_GeV:.3e} GeV",
            "",
            "Toy model used (explicit; not derived here):",
            "- |S_mu| ~ c_spin * (polarization * spin_per_particle * n) / M_eff^2",
            f"- observable proxy: Δν ≈ 2|b|·(GeV→Hz), with |b|≈(3/4)|S| and GeV→Hz={GEV_TO_HZ:.3e}",
            "- coupling_kind=Mpl_reduced uses M_eff = Mpl(reduced)",
            "- coupling_kind=TFPT_M uses M_eff = TFPT M scale (from effective_action_r2 / constants)",
            "",
            "Named anchor scenarios (c_spin=1):",
        ]
        for r in anchor_rows:
            lines.append(f"- {r}")
        lines += [
            "",
            "Proposed regime JSON entries (copy into tfpt_suite/data/torsion_regimes.json if desired):",
        ]
        for p in proposals:
            lines.append(json.dumps(p, indent=2, sort_keys=True))
        lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module is a *design tool*. It does not claim an operator-level derivation of the spin→torsion coupling.",
            "- The main suite's `torsion_bounds_mapping` already supports evaluating `spin_polarized_medium` regimes; these proposals are meant to populate that config.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"torsion_regime_designer_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_torsion_regime_designer(
                out_dir=out_dir,
                tight_bound_GeV=float(tight.abs_max_GeV),
                anchor_rows=anchor_rows,
                candidates=candidates,
            )
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "bounds": {"file": str(bounds_path), "tightest": tight.__dict__, "all": [b.__dict__ for b in vetted]},
                "anchors": anchor_rows,
                "sampling": {"samples": samples, "candidates_kept": len(candidates), "best_mpl_reduced": best_mpl, "best_tfpt_M": best_tfpt},
                "proposals": proposals,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

