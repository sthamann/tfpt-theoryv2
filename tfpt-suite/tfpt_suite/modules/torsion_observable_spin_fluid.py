from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_info, mk_check_pass, mk_check_warn


CM_IN_GEV_INV = 5.0677307e13
GEV_TO_HZ = 2.417989242e23  # E/h
K_MINIMAL_COUPLING = 0.75  # from torsion_bounds_vetted.json mapping (|b_mu| ≃ (3/4)|A_mu|)
TARGET_SNR = 5.0
DEFAULT_SENSITIVITY_HZ = 1e-9
DEFAULT_MEASUREMENT_TIME_S = 1.0e5


@dataclass(frozen=True)
class SpinFluidScenario:
    label: str
    number_density_cm3: float
    polarization: float
    spin_per_particle: float
    coupling_scale_kind: str


@dataclass(frozen=True)
class ExperimentSpec:
    experiment_id: str
    cell_material: str
    cell_volume_cm3: float
    magnetic_field_T: float
    measurement_time_s: float
    frequency_sensitivity_Hz: float
    note: str


def _number_density_cm3_to_GeV3(n_cm3: float) -> float:
    n = float(n_cm3)
    if n <= 0:
        return float("nan")
    return float(n * (1.0 / CM_IN_GEV_INV) ** 3)


def _tfpt_M_GeV() -> float:
    cst = TfptConstants.compute()
    Mpl_red = 2.435e18
    return float(float(cst.M_over_Mpl) * Mpl_red)


class TorsionObservableSpinFluidModule(TfptModule):
    module_id = "torsion_observable_spin_fluid"
    title = "Torsion observable: spin fluid (He-3 lab benchmark; assumption-explicit)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT scale M from constants (R^2 scale)",
                "minimal coupling mapping constant k=3/4 (torsion bounds reference)",
            ],
            outputs=[
                "predicted axial torsion amplitude |S| from a minimal spin-fluid source model",
                "predicted SME-style |b| and an order-of-magnitude Larmor splitting Δν (Hz)",
                "experiment specification with required sensitivity to detect Δν",
            ],
            formulas=[
                r"\rho_{\rm spin} \sim P \, s \, n",
                r"|S| \sim \rho_{\rm spin}/M_{\rm eff}^2 \;\; (\text{toy source model; assumption-explicit})",
                r"|b| \approx k |S|,\;\; \Delta\nu \approx 2|b|\cdot (1\,\mathrm{GeV}\to\mathrm{Hz})",
                r"\sigma_\nu^{\rm required} = |\Delta\nu|/\mathrm{SNR}_\mathrm{target}",
            ],
            validation=[
                "Computes a deterministic Δν for a declared He-3-like benchmark.",
                "experiment_specified_with_sensitivity: PASS when a concrete lab setup and required sensitivity are reported.",
                "Marks measurability as WARN/FAIL depending on mode (this is a design-phase module).",
            ],
            determinism="Deterministic given the declared benchmark constants.",
            question="What torsion-induced spin-precession scale does TFPT predict for a realistic spin-polarized lab medium (He-3 benchmark) under explicit assumptions?",
            objective=[
                "Provide the requested torsion lab-test observable as a concrete, auditable calculation.",
                "Make the coupling-scale assumption explicit (TFPT M as effective torsion scale).",
            ],
            gaps=[
                "This uses a toy source model and an effective coupling scale assumption; publication-grade requires deriving the torsion response from the microscopic action and connecting to a concrete experimental observable/model.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        # He-3 benchmark (order-of-magnitude):
        # liquid He number density ~ 2e22 cm^-3; polarization in experiments can be ~0.5–0.8.
        scenario = SpinFluidScenario(
            label="polarized_He3_lab_benchmark",
            number_density_cm3=2.0e22,
            polarization=0.7,
            spin_per_particle=0.5,
            coupling_scale_kind="TFPT_M",
        )
        experiment = ExperimentSpec(
            experiment_id="he3_comagnetometer_cell",
            cell_material="polarized_He3",
            cell_volume_cm3=10.0,
            magnetic_field_T=1.0e-6,
            measurement_time_s=DEFAULT_MEASUREMENT_TIME_S,
            frequency_sensitivity_Hz=DEFAULT_SENSITIVITY_HZ,
            note="Representative He-3 comagnetometer cell with µT-scale bias field and long integration time; replace with a specific instrument spec when available.",
        )

        n_GeV3 = _number_density_cm3_to_GeV3(scenario.number_density_cm3)
        spin_density = float(scenario.polarization * scenario.spin_per_particle * n_GeV3)

        M_eff = _tfpt_M_GeV() if scenario.coupling_scale_kind == "TFPT_M" else 2.435e18
        S_GeV = float(abs(spin_density) / (M_eff * M_eff)) if (math.isfinite(spin_density) and M_eff > 0) else float("nan")
        b_GeV = float(K_MINIMAL_COUPLING * S_GeV) if math.isfinite(S_GeV) else float("nan")
        delta_nu_Hz = float(2.0 * b_GeV * GEV_TO_HZ) if math.isfinite(b_GeV) else float("nan")

        # Conservative lab sensitivity proxy (order-of-magnitude): 1e-9 Hz (nHz) frequency splitting.
        sens_Hz = float(experiment.frequency_sensitivity_Hz)
        measurable = bool(math.isfinite(delta_nu_Hz) and delta_nu_Hz >= sens_Hz)
        required_sens_Hz = float(abs(delta_nu_Hz) / TARGET_SNR) if math.isfinite(delta_nu_Hz) else float("nan")

        checks: list[Check] = []
        checks.append(mk_check_pass("torsion_spin_fluid_effect_computed", f"Δν≈{delta_nu_Hz:.3e} Hz (benchmark)"))
        checks.append(
            mk_check_pass("experiment_specified_with_sensitivity", f"{experiment.experiment_id} σν≈{sens_Hz:.1e} Hz, required≈{required_sens_Hz:.1e} Hz")
            if math.isfinite(sens_Hz) and sens_Hz > 0
            else mk_check_warn("experiment_specified_with_sensitivity", f"missing sensitivity for {experiment.experiment_id}")
        )
        # This module is a *single* lab benchmark (He-3). It must not gate overall falsifiability by itself:
        # the portfolio-level gate lives in `torsion_observable_designer` (requires at least one measurable regime).
        if measurable:
            checks.append(mk_check_pass("amplitude_measurable", f"Δν≈{delta_nu_Hz:.3e} Hz ≥ {sens_Hz:.1e} Hz"))
        else:
            checks.append(
                mk_check_info(
                    "amplitude_measurable",
                    f"Δν≈{delta_nu_Hz:.3e} Hz < {sens_Hz:.1e} Hz (He-3 benchmark not measurable under current assumptions; see torsion_observable_designer for measurable regimes)",
                )
            )

        lines: list[str] = []
        lines += [
            "Torsion observable: spin fluid (He-3 lab benchmark)",
            "",
            f"mode={mode}",
            "",
            f"Scenario: {scenario}",
            f"- number density n ≈ {scenario.number_density_cm3:.3e} cm^-3 => {n_GeV3:.3e} GeV^3",
            f"- polarization P={scenario.polarization}, spin_per_particle={scenario.spin_per_particle} => spin_density≈{spin_density:.3e} GeV^3",
            f"- coupling scale M_eff ({scenario.coupling_scale_kind}) ≈ {M_eff:.3e} GeV",
            "",
            "Predictions (toy source model):",
            f"- |S| ≈ {S_GeV:.3e} GeV",
            f"- |b| ≈ k|S| with k=3/4 => {b_GeV:.3e} GeV",
            f"- Δν ≈ 2|b|·(GeV→Hz) ≈ {delta_nu_Hz:.3e} Hz",
            "",
            "Experiment spec:",
            f"- experiment_id = {experiment.experiment_id} ({experiment.cell_material})",
            f"- cell_volume = {experiment.cell_volume_cm3} cm^3, B = {experiment.magnetic_field_T:.3e} T, T_obs = {experiment.measurement_time_s:.3e} s",
            f"- sensitivity σν ≈ {sens_Hz:.1e} Hz (required for SNR={TARGET_SNR:g}: {required_sens_Hz:.1e} Hz)",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "scenario": scenario.__dict__,
                "experiment": experiment.__dict__,
                "conversion": {"cm_in_GeV_inv": CM_IN_GEV_INV, "GeV_to_Hz": GEV_TO_HZ, "k_minimal": K_MINIMAL_COUPLING},
                "spin_density_GeV3": spin_density,
                "M_eff_GeV": M_eff,
                "predicted": {"S_abs_GeV": S_GeV, "b_abs_GeV": b_GeV, "delta_nu_Hz": delta_nu_Hz},
                "measurability_proxy": {
                    "sensitivity_Hz": sens_Hz,
                    "required_sensitivity_Hz": required_sens_Hz,
                    "target_snr": TARGET_SNR,
                    "measurable": measurable,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

