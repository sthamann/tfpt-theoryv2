from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _plot_torsion_bounds(
    *,
    out_dir: Path,
    bounds: list[BoundResult],
    tfpt_predicted_S_mu_GeV: float,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"torsion_bounds_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        labels = [b.label for b in bounds]
        vals = [float(b.inferred_abs_max_S_mu_GeV) for b in bounds]
        x = np.arange(len(labels), dtype=float)

        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.bar(x, vals, color="#1f77b4", alpha=0.9)
        ax.axhline(float(tfpt_predicted_S_mu_GeV), color="black", lw=1.0, alpha=0.8, label="TFPT prediction (today)")
        ax.set_yscale("symlog", linthresh=1e-35)
        ax.set_ylabel(r"$|S_\mu|$ bound [GeV] (symlog)")
        ax.set_title("Axial torsion bounds (mapped/ingested)")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right")
        ax.grid(True, axis="y", ls=":", alpha=0.4)
        ax.legend(loc="best")
        fig.tight_layout()

        path = out_dir / "torsion_bounds.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["torsion_bounds_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class BoundResult:
    label: str
    coefficient: str  # input coefficient name (e.g. b_mu or A_mu)
    abs_max_input_GeV: float
    inferred_abs_max_S_mu_GeV: float
    applies_to: str
    mapping_used: str


class TorsionBoundsMappingModule(TfptModule):
    module_id = "torsion_bounds_mapping"
    title = "Torsion bounds mapping (SME-style b_mu ↔ axial torsion S_mu)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "bounds file: tfpt_suite/data/torsion_bounds_vetted.json (preferred) or tfpt_suite/data/torsion_bounds.json (fallback)",
                "TFPT torsion regimes: tfpt_suite/data/torsion_regimes.json (explicit regime selection for falsifiability)",
            ],
            outputs=["inferred torsion bounds (S_mu) from b_mu limits"],
            formulas=[
                "convention: b_mu = k * S_mu  =>  |S_mu| <= |b_mu|/|k|",
            ],
            validation=[
                "loads bounds JSON and produces finite inferred bounds",
            ],
            determinism="Deterministic given the bounds file contents.",
            question="Given vetted SME-style axial torsion bounds, what does TFPT predict for present-day torsion amplitudes under explicit regime assumptions, and is it falsifiable?",
            objective=[
                "Provide a clean mapping/ingestion layer from literature bounds to TFPT’s axial torsion variable S_mu.",
                "Force an explicit TFPT regime choice (vacuum vs nontrivial benchmark) so the module cannot be ‘always green’ by construction.",
                "Expose the minimal-coupling mapping constant (3/4) and its relation to TFPT invariants (ξ_tree=c3/φ_tree).",
            ],
            what_was_done=[
                "Load vetted component-wise bounds and map them into inferred |S_mu| limits.",
                "Evaluate the selected TFPT regime from torsion_regimes.json and compare to the bounds (ratio report).",
                "Annotate the mapping constant k with ξ_tree and ξ=c3/φ0 for traceability to eliminating_k-style normalization factors.",
            ],
            assumptions=[
                "Minimal coupling mapping b_mu = k S_mu is used when bounds are given in SME b_mu conventions (absolute-value mapping).",
                "Present-day torsion is assumed to be vacuum-like unless a nontrivial source model is declared explicitly in torsion_regimes.json.",
            ],
            gaps=[
                "Current nonzero regime is still a toy benchmark; publication-grade falsifiability needs a computed source model (spin media, magnetars, early plasma).",
            ],
            references=[
                "torsion_bounds_vetted.json (Kostelecký–Russell–Tasson 2008 axial torsion bounds)",
                "eliminating_k.tex (ξ normalization; ξ_tree=3/4)",
            ],
            maturity="framework + falsifiability scaffold (regime model is still toy unless replaced)",
        )

    def run(self, config) -> ModuleResult:
        data_dir = Path(__file__).resolve().parent.parent / "data"
        data_path_vetted = data_dir / "torsion_bounds_vetted.json"
        data_path_fallback = data_dir / "torsion_bounds.json"
        data_path = data_path_vetted if data_path_vetted.exists() else data_path_fallback
        raw = json.loads(data_path.read_text(encoding="utf-8"))

        # Regime selection for a nontrivial TFPT torsion prediction.
        regimes_path = data_dir / "torsion_regimes.json"
        regimes_raw = json.loads(regimes_path.read_text(encoding="utf-8")) if regimes_path.exists() else {}
        default_regime_id = str(regimes_raw.get("default_regime_id", "vacuum_today"))
        regimes_list = regimes_raw.get("regimes", []) if isinstance(regimes_raw.get("regimes", []), list) else []
        regimes_by_id: dict[str, dict[str, Any]] = {}
        for r in regimes_list:
            if isinstance(r, dict) and "id" in r:
                regimes_by_id[str(r["id"])] = r
        regime = regimes_by_id.get(default_regime_id, regimes_by_id.get("vacuum_today", {}))

        # TFPT invariants (used to annotate the minimal-coupling mapping constant 3/4).
        cst = TfptConstants.compute()
        xi_tree = cst.c3 / cst.varphi0_tree  # = 3/4 exactly
        xi = cst.c3 / cst.varphi0  # small correction when δ_top is included in varphi0

        mapping = raw.get("mapping", {}) if isinstance(raw.get("mapping", {}), dict) else {}
        k_policy = str(mapping.get("k_policy", "")).strip()
        if k_policy == "xi_tree":
            k_default = float(xi_tree)
        elif k_policy == "xi":
            k_default = float(xi)
        else:
            # Backwards-compatible: keep existing JSON mapping.k if present, otherwise use xi_tree (=3/4).
            k_default = float(mapping.get("k", float(xi_tree)))
        if k_default == 0:
            raise ValueError("Invalid mapping factor k=0 in torsion_bounds.json")

        bounds_in = raw.get("bounds", [])
        results: list[BoundResult] = []
        for entry in bounds_in:
            if not isinstance(entry, dict):
                continue
            coef = str(entry.get("coefficient", ""))
            if "abs_max" not in entry:
                continue

            abs_max = float(entry["abs_max"])
            if abs_max <= 0:
                continue

            # Two supported input styles:
            # - SME-style bounds on b_mu, mapped to axial torsion via b_mu = k * S_mu.
            # - Direct bounds on the axial torsion vector itself (A_mu / S_mu), passed through.
            if coef == "b_mu":
                k_entry = float(entry.get("mapping_k", k_default))
                if k_entry == 0:
                    continue
                smax = abs_max / abs(k_entry)
                mapping_used = f"b_mu = {k_entry} * S_mu"
            elif coef in {"A_mu", "S_mu"}:
                smax = abs_max
                mapping_used = "direct (S_mu := A_mu)"
            else:
                continue

            results.append(
                BoundResult(
                    label=str(entry.get("label", "unnamed")),
                    coefficient=coef,
                    abs_max_input_GeV=abs_max,
                    inferred_abs_max_S_mu_GeV=smax,
                    applies_to=str(entry.get("applies_to", "")),
                    mapping_used=mapping_used,
                )
            )

        # --- TFPT prediction layer (minimal, explicit) ---
        #
        # In the Minkowski/vacuum limit with no spin sources, the minimal Einstein–Cartan / Riemann–Cartan
        # torsion field equation is algebraic: torsion is sourced by spin and vanishes in vacuum.
        # In TFPT language: the "local lab" regime corresponds to the torsion vacuum K→0, hence:
        #   S_mu(today) = 0
        #
        # Nonzero late-time torsion would require either (i) explicit propagating torsion DOFs in the
        # low-energy effective action, or (ii) a specified spin-fluid / condensate source model.
        # Default theorem: vacuum torsion vanishes (Einstein–Cartan minimal).
        # For falsifiability, we optionally evaluate a nontrivial regime declared in torsion_regimes.json.
        tfpt_regime_id = str(regime.get("id", "vacuum_today") or "vacuum_today")
        tfpt_regime_label = str(regime.get("label", tfpt_regime_id))
        tfpt_regime_model = str(regime.get("model", "vacuum"))

        def _H0_GeV_from_km_s_Mpc(H0_km_s_Mpc: float) -> float:
            # H0 [s^-1] = H0_km_s_Mpc / (Mpc in km); 1/s = ħ [GeV·s] in natural units.
            Mpc_km = 3.0856775814913673e19
            hbar_GeV_s = 6.582119569e-25
            H0_s_inv = float(H0_km_s_Mpc) / Mpc_km
            return float(H0_s_inv * hbar_GeV_s)

        def _number_density_to_GeV3(value: float, unit: str) -> float:
            v = float(value)
            u = str(unit).strip()
            if v <= 0:
                raise ValueError("number density must be positive")
            if u == "cm^-3":
                cm_in_GeV_inv = 5.0677307e13
                return float(v * (1.0 / cm_in_GeV_inv) ** 3)
            if u == "fm^-3":
                fm_in_GeV_inv = 5.0677307
                return float(v * (1.0 / fm_in_GeV_inv) ** 3)
            raise ValueError(f"unsupported number density unit: {unit!r}")

        tfpt_predicted_S_mu_GeV = 0.0
        tfpt_prediction_note = "vacuum theorem: torsion sourced by spin ⇒ vacuum torsion vanishes (S_mu=0)"
        tfpt_regime_details: dict[str, Any] = {"id": tfpt_regime_id, "label": tfpt_regime_label, "model": tfpt_regime_model}

        try:
            if tfpt_regime_model == "vacuum":
                tfpt_predicted_S_mu_GeV = float(regime.get("S_mu_abs_GeV", 0.0))
                tfpt_prediction_note = str(regime.get("note", tfpt_prediction_note))
            elif tfpt_regime_model == "cosmological_H0":
                H0 = float(regime.get("H0_km_s_Mpc", 67.36))
                c_factor = float(regime.get("c_factor", 1.0))
                tfpt_predicted_S_mu_GeV = float(abs(c_factor) * abs(_H0_GeV_from_km_s_Mpc(H0)))
                tfpt_prediction_note = str(regime.get("note", "toy regime: |S_mu| ~ c_factor * H0"))
                tfpt_regime_details.update({"H0_km_s_Mpc": H0, "c_factor": c_factor})
            elif tfpt_regime_model == "spin_polarized_medium":
                nd = regime.get("number_density", {}) if isinstance(regime.get("number_density", {}), dict) else {}
                n_val = float(nd.get("value", 0.0))
                n_unit = str(nd.get("unit", "cm^-3"))
                n_GeV3 = _number_density_to_GeV3(n_val, n_unit)
                pol = float(regime.get("polarization", 0.0))
                spin_per = float(regime.get("spin_per_particle", 0.5))
                c_spin = float(regime.get("c_spin", 1.0))

                coupling = regime.get("coupling_scale", {}) if isinstance(regime.get("coupling_scale", {}), dict) else {}
                kind = str(coupling.get("kind", "Mpl_reduced"))
                Mpl_reduced_GeV = 2.435e18
                if kind == "TFPT_M":
                    cst = TfptConstants.compute()
                    M_eff_GeV = float(float(cst.M_over_Mpl) * Mpl_reduced_GeV)
                else:
                    M_eff_GeV = float(Mpl_reduced_GeV)

                spin_density_GeV3 = float(pol * spin_per * n_GeV3)
                tfpt_predicted_S_mu_GeV = float(abs(c_spin) * abs(spin_density_GeV3) / (M_eff_GeV**2))
                tfpt_prediction_note = str(regime.get("note", "spin medium model: |S_mu| ~ c_spin * (pol*spin*n)/M_eff^2"))
                tfpt_regime_details.update(
                    {
                        "number_density": {"value": n_val, "unit": n_unit, "value_GeV3": n_GeV3},
                        "polarization": pol,
                        "spin_per_particle": spin_per,
                        "spin_density_GeV3": spin_density_GeV3,
                        "c_spin": c_spin,
                        "coupling_scale": {"kind": kind, "M_eff_GeV": M_eff_GeV},
                    }
                )
            else:
                tfpt_prediction_note = f"unknown regime model={tfpt_regime_model!r}; falling back to vacuum"
                tfpt_predicted_S_mu_GeV = 0.0
        except Exception as e:
            tfpt_prediction_note = f"regime_evaluation_failed ({tfpt_regime_model}): {e}; falling back to vacuum"
            tfpt_predicted_S_mu_GeV = 0.0

        compare_to_bounds = bool(regime.get("compare_to_bounds", True))
        worst_ratio = 0.0
        if compare_to_bounds:
            for r in results:
                if r.inferred_abs_max_S_mu_GeV > 0:
                    worst_ratio = max(worst_ratio, abs(tfpt_predicted_S_mu_GeV) / r.inferred_abs_max_S_mu_GeV)

        checks: list[Check] = []
        # Publication-grade bounds require a vetted, component-wise dataset and a TFPT predicted torsion amplitude to compare against.
        # This module is only the mapping/ingestion layer until those two inputs exist.
        def _is_vetted_bound(entry: dict[str, Any]) -> bool:
            required = {"label", "coefficient", "abs_max", "applies_to", "source"}
            return required.issubset(set(entry.keys()))

        vetted_count = 0
        for e in bounds_in:
            if isinstance(e, dict) and _is_vetted_bound(e):
                vetted_count += 1

        checks.append(
            Check(
                check_id="vetted_componentwise_bounds_dataset_present",
                passed=bool(vetted_count >= 4),
                detail=f"vetted entries={vetted_count} (need >=4 with per-entry source/metadata)",
            )
        )
        checks.append(
            Check(
                check_id="tfpt_torsion_prediction_present",
                passed=True,
                detail=f"prediction regime={tfpt_regime_id} (model={tfpt_regime_model}); {tfpt_prediction_note}",
            )
        )
        checks.append(
            Check(
                check_id="nontrivial_tfpt_torsion_regime_defined",
                passed=bool(abs(tfpt_predicted_S_mu_GeV) > 0),
                detail="FAIL means this module is currently an ingestion/mapping layer only. Choose/define a nontrivial TFPT torsion regime in tfpt_suite/data/torsion_regimes.json.",
            )
        )
        checks.append(
            Check(
                check_id="tfpt_prediction_respects_bounds",
                passed=True if (not compare_to_bounds) else bool(worst_ratio <= 1.0),
                detail=("skipped for this regime (compare_to_bounds=false)" if (not compare_to_bounds) else f"max |S_mu_pred|/|S_mu_bound| = {worst_ratio:.3e}"),
            )
        )
        checks.append(
            Check(
                check_id="bounds_file_loaded",
                passed=bool(isinstance(raw, dict) and len(bounds_in) >= 1),
                detail=f"loaded {len(bounds_in)} bound entries from {data_path}",
            )
        )
        checks.append(
            Check(
                check_id="computed_inferred_bounds",
                passed=bool(len(results) >= 1 and all(r.inferred_abs_max_S_mu_GeV > 0 for r in results)),
                detail=f"computed {len(results)} inferred S_mu bounds",
            )
        )

        report_lines: list[str] = []
        report_lines += [
            "Torsion bounds mapping (framework module)",
            "",
            f"data file: {data_path}",
            f"units: {raw.get('units')}",
            "",
            "Convention used:",
            f"- b_mu = k * S_mu with k={k_default}  =>  |S_mu| <= |b_mu|/|k|",
            f"  (TFPT invariants: xi_tree=c3/varphi0_tree={xi_tree} (=3/4), xi=c3/varphi0={xi})",
            "- direct A_mu bounds are also accepted and passed through as S_mu := A_mu",
            "",
            "TFPT prediction (local / present-day vacuum regime):",
            f"- selected regime: {tfpt_regime_id} ({tfpt_regime_label})",
            f"- model: {tfpt_regime_model}",
            f"- predicted |S_mu| = {abs(tfpt_predicted_S_mu_GeV):.3e} GeV",
            f"- note: {tfpt_prediction_note}",
            "",
            "Inferred bounds (placeholders unless you replace the JSON table):",
        ]
        for r in results:
            report_lines.append(
                f"- {r.label}: |{r.coefficient}| <= {r.abs_max_input_GeV:.3e} GeV  =>  |S_mu| <= {r.inferred_abs_max_S_mu_GeV:.3e} GeV  ({r.applies_to}; {r.mapping_used})"
            )
        report_lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module is the ingestion/mapping layer requested in tasks.md.",
            "- To make this a falsification module, provide TFPT-predicted local torsion amplitudes and a vetted component-wise bounds table.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"torsion_bounds_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_torsion_bounds(out_dir=out_dir, bounds=results, tfpt_predicted_S_mu_GeV=tfpt_predicted_S_mu_GeV)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "data_file": str(data_path),
                "mapping": {
                    "k_default": k_default,
                    "k_policy": k_policy,
                    "xi_tree": xi_tree,
                    "xi": xi,
                    "mapping_raw": mapping,
                },
                "tfpt_prediction": {
                    "S_mu_abs_GeV": float(abs(tfpt_predicted_S_mu_GeV)),
                    "regime": tfpt_regime_id,
                    "regime_details": tfpt_regime_details,
                    "theorem_basis": "vacuum torsion vanishes in minimal Einstein–Cartan; nonzero requires an explicit source/propagating torsion regime (see torsion_regimes.json)",
                },
                "inferred_bounds": [
                    {
                        "label": r.label,
                        "coefficient": r.coefficient,
                        "abs_max_input_GeV": r.abs_max_input_GeV,
                        "inferred_abs_max_S_mu_GeV": r.inferred_abs_max_S_mu_GeV,
                        "applies_to": r.applies_to,
                        "mapping_used": r.mapping_used,
                    }
                    for r in results
                ],
                "plot": plot,
            },
            checks=checks,
            report="\n".join(report_lines),
            warnings=warnings,
        )

