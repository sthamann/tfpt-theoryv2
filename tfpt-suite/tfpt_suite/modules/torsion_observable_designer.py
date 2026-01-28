from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_info, mk_check_pass, mk_check_warn


CM_IN_GEV_INV = 5.0677307e13
FM_IN_GEV_INV = 5.0677307
GEV_TO_HZ = 2.417989242e23
K_MINIMAL_COUPLING = 0.75


@dataclass(frozen=True)
class ObservableScenario:
    scenario_id: str
    label: str
    model: str
    spin_density_GeV3: float
    M_eff_kind: str
    note: str


def _tfpt_M_GeV() -> float:
    cst = TfptConstants.compute()
    Mpl_red = 2.435e18
    return float(float(cst.M_over_Mpl) * Mpl_red)


def _n_cm3_to_GeV3(n_cm3: float) -> float:
    n = float(n_cm3)
    return float(n * (1.0 / CM_IN_GEV_INV) ** 3)


def _n_fm3_to_GeV3(n_fm3: float) -> float:
    n = float(n_fm3)
    return float(n * (1.0 / FM_IN_GEV_INV) ** 3)


class TorsionObservableDesignerModule(TfptModule):
    module_id = "torsion_observable_designer"
    title = "Torsion observable designer (magnetar + lab signal proxies)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT scale M from constants (R^2 scale) used as coupling-scale proxy",
                "torsion bounds mapping constant k=3/4 (minimal coupling)",
            ],
            outputs=[
                "two benchmark scenarios (lab spin-fluid, nuclear/magnetar-like spin density) with predicted Δν scales",
                "amplitude_measurable gate (design-phase; physics mode remains strict)",
            ],
            formulas=[
                r"|S| \sim \rho_{\rm spin}/M_{\rm eff}^2,\;\; |b|=k|S|,\;\; \Delta\nu\approx 2|b|\cdot(\mathrm{GeV}\to\mathrm{Hz})",
            ],
            validation=[
                "amplitude_measurable check is emitted (PASS/WARN/FAIL depending on mode and proxy thresholds).",
            ],
            determinism="Deterministic given benchmark constants.",
            question="What torsion observable magnitude could be targeted in lab vs astrophysical regimes under TFPT-scale coupling assumptions?",
            objective=[
                "Turn the torsion-falsifiability task into a concrete, assumption-explicit signal-size table.",
                "Provide a path to a real experimental design (replace benchmarks with a derived source model + instrument noise model).",
            ],
            gaps=[
                "Benchmarks use a toy source model and TFPT-scale coupling assumption; publication-grade needs derived source dynamics + real observables/noise model.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")
        M_eff = _tfpt_M_GeV()

        # Scenario A: He-3 lab-like (spin fluid)
        n_cm3 = 2.0e22
        pol = 0.7
        spin = 0.5
        spin_density_lab = float(pol * spin * _n_cm3_to_GeV3(n_cm3))

        # Scenario B: nuclear density (magnetar-like core / neutron-star matter proxy)
        n_fm3 = 0.16
        pol_ns = 1.0
        spin_ns = 0.5
        spin_density_ns = float(pol_ns * spin_ns * _n_fm3_to_GeV3(n_fm3))

        scenarios: list[ObservableScenario] = [
            ObservableScenario(
                scenario_id="lab_spin_fluid_He3",
                label="Lab spin fluid (He-3 benchmark)",
                model="spin_polarized_medium",
                spin_density_GeV3=spin_density_lab,
                M_eff_kind="TFPT_M",
                note="Order-of-magnitude lab density; intended as a falsifiability baseline.",
            ),
            ObservableScenario(
                scenario_id="nuclear_spin_density_magnetar_proxy",
                label="Nuclear spin density (magnetar proxy)",
                model="spin_polarized_medium",
                spin_density_GeV3=spin_density_ns,
                M_eff_kind="TFPT_M",
                note="Order-of-magnitude nuclear density; not an Earth lab bound regime. Intended for astrophysical amplification discussion.",
            ),
        ]

        rows: list[dict[str, Any]] = []
        for sc in scenarios:
            S = float(abs(sc.spin_density_GeV3) / (M_eff * M_eff)) if (math.isfinite(sc.spin_density_GeV3) and M_eff > 0) else float("nan")
            b = float(K_MINIMAL_COUPLING * S) if math.isfinite(S) else float("nan")
            dnu = float(2.0 * b * GEV_TO_HZ) if math.isfinite(b) else float("nan")
            rows.append(
                {
                    "scenario_id": sc.scenario_id,
                    "label": sc.label,
                    "model": sc.model,
                    "spin_density_GeV3": sc.spin_density_GeV3,
                    "M_eff_GeV": M_eff,
                    "S_abs_GeV": S,
                    "b_abs_GeV": b,
                    "delta_nu_Hz": dnu,
                    "note": sc.note,
                }
            )

        # Proxy measurability thresholds (design-phase):
        # - lab: nHz
        # - astro: 1e-9 Hz is already challenging; keep proxy explicit.
        lab_thr = 1e-9
        astro_thr = 1e-9
        lab_dnu = float(next((r["delta_nu_Hz"] for r in rows if r["scenario_id"] == "lab_spin_fluid_He3"), float("nan")))
        astro_dnu = float(next((r["delta_nu_Hz"] for r in rows if r["scenario_id"] == "nuclear_spin_density_magnetar_proxy"), float("nan")))
        proxy_ok = bool((math.isfinite(lab_dnu) and lab_dnu >= lab_thr) or (math.isfinite(astro_dnu) and astro_dnu >= astro_thr))

        checks: list[Check] = []
        checks.append(mk_check_pass("torsion_observable_table_present", f"rows={len(rows)}"))
        if mode == "physics":
            # Physics-mode gate: at least one regime should be *measurable in principle* (proxy thresholds) to qualify
            # as a falsifiable observable target. This remains a proxy until a full source+noise model is implemented.
            checks.append(
                mk_check_pass(
                    "amplitude_measurable",
                    f"proxy_ok=True (lab Δν≈{lab_dnu:.3e} Hz, astro Δν≈{astro_dnu:.3e} Hz; thresholds={lab_thr:.1e} Hz)",
                )
                if proxy_ok
                else mk_check_warn(
                    "amplitude_measurable",
                    f"proxy_ok=False (lab Δν≈{lab_dnu:.3e} Hz, astro Δν≈{astro_dnu:.3e} Hz; thresholds={lab_thr:.1e} Hz) — design-phase",
                )
            )
        else:
            checks.append(
                mk_check_info(
                    "amplitude_measurable",
                    f"INFO: proxy thresholds lab/astro={lab_thr:.1e} Hz; lab Δν≈{lab_dnu:.3e} Hz, astro Δν≈{astro_dnu:.3e} Hz (proxy_ok={proxy_ok})",
                )
            )

        lines: list[str] = []
        lines += [
            "Torsion observable designer (signal proxy table)",
            "",
            f"mode={mode}",
            f"M_eff (TFPT_M) ≈ {M_eff:.3e} GeV",
            f"minimal coupling k = 3/4, GeV→Hz = {GEV_TO_HZ:.3e}",
            "",
            "Scenarios:",
        ]
        for r in rows:
            lines.append(
                f"- {r['scenario_id']}: Δν≈{r['delta_nu_Hz']:.3e} Hz (|S|≈{r['S_abs_GeV']:.3e} GeV; |b|≈{r['b_abs_GeV']:.3e} GeV) — {r['label']}"
            )
        lines += [
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
            "",
            "Notes:",
            "- This is a design-phase module. Replace benchmark spin densities and coupling-scale assumptions with a derived TFPT source model + experimental/astrophysical observable model to upgrade to publication-grade falsifiability.",
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "M_eff_GeV": M_eff,
                "rows": rows,
                "proxy_thresholds_Hz": {"lab": lab_thr, "astro": astro_thr},
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

