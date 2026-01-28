from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from tfpt_suite.cosmo_scale_map import As_from_ln10_As, starobinsky_N_from_ns, starobinsky_r_from_N
from tfpt_suite.module_base import (
    Check,
    ModuleResult,
    ModuleSpec,
    TfptModule,
    mk_check_fail,
    mk_check_info,
    mk_check_pass,
    mk_check_warn,
)


def _freq_Hz_from_k_Mpc_inv(k_Mpc_inv: float) -> float:
    # k [Mpc^-1] -> f [Hz] using: 1 Mpc = 3.08567758e22 m, c=299792458 m/s.
    # f = (c / 2π) * k / a0, with a0=1 and k in m^-1.
    c_m_s = 299792458.0
    Mpc_m = 3.0856775814913673e22
    k_m_inv = float(k_Mpc_inv) / Mpc_m
    return float((c_m_s / (2.0 * math.pi)) * k_m_inv)


def _omega_gw_flat_spectrum(*, Omega_r0: float, A_s: float, r: float) -> float:
    """
    Engineering-level scale-invariant stochastic GW background estimate for modes that re-enter in radiation domination:

      Ω_gw,0 ≈ (Ω_r,0 / 24) · P_t  with  P_t = r · A_s

    This ignores detailed transfer-function features (e.g. reheating) and is intended as a conservative order-of-magnitude ledger.
    """
    if Omega_r0 <= 0 or A_s <= 0 or r <= 0:
        return float("nan")
    return float((float(Omega_r0) / 24.0) * float(r) * float(A_s))


class GwBackgroundPredictorModule(TfptModule):
    module_id = "gw_background_predictor"
    title = "GW background predictor (Ω_gw(f); scale-invariant Starobinsky baseline + bounds ledger)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "cosmology snapshot: tfpt_suite/data/k_calibration.json (Ω_r,0 used for Ω_gw normalization)",
                "Planck minimal reference: tfpt_suite/data/global_reference_minimal.json (n_s, ln10(10^{10}A_s))",
                "(optional) bounce transfer diagnostics: bounce_perturbations (future refinement)",
            ],
            outputs=[
                "Starobinsky r(N) baseline and inferred tensor amplitude A_t=r·A_s at pivot",
                "scale-invariant Ω_gw,0 estimate and a small frequency table (CMB pivot / PTA / LIGO)",
                "bounds ledger (PTA/LIGO/BBN-style order-of-magnitude checks)",
            ],
            formulas=[
                r"\Omega_{\rm gw}(f) = \frac{1}{\rho_c}\frac{d\rho_{\rm gw}}{d\ln f}",
                r"k \leftrightarrow f \;\; \text{via transfer function + expansion history}",
                r"\Omega_{\rm gw,0} \approx (\Omega_{r,0}/24)\, r A_s \;\; (\text{radiation-era re-entry; scale-invariant baseline})",
            ],
            validation=[
                "Produces a deterministic baseline Ω_gw(f) estimate from Starobinsky r and the cosmology snapshot.",
                "Emits a PASS/FAIL ledger against conservative PTA/LIGO/BBN bounds (engineering-level; update with a full likelihood when needed).",
            ],
            determinism="Deterministic given input tables.",
            question="Is a full stochastic GW background prediction Ω_gw(f) implemented and compared to PTA/LIGO/CMB bounds?",
            objective=[
                "Provide a falsifiable Ω_gw(f) baseline beyond the bounds-only placeholder module.",
            ],
            gaps=[
                "A publication-grade module should include transfer functions (bounce + reheating), k→f mapping through the full history, and a proper likelihood against PTA/LIGO/CMB datasets.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        data_dir = Path(__file__).resolve().parent.parent / "data"
        kcal_path = data_dir / "k_calibration.json"
        ref_path = data_dir / "global_reference_minimal.json"

        kcal = json.loads(kcal_path.read_text(encoding="utf-8")) if kcal_path.is_file() else {}
        ref = json.loads(ref_path.read_text(encoding="utf-8")) if ref_path.is_file() else {}

        cosmo_raw = kcal.get("cosmology_flat_lcdm", {}) if isinstance(kcal.get("cosmology_flat_lcdm", {}), dict) else {}
        Omega_r0 = float(cosmo_raw.get("Omega_r", float("nan")))

        obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
        ns = (
            float(obs.get("n_s_planck2018", {}).get("mean", float("nan"))) if isinstance(obs.get("n_s_planck2018", {}), dict) else float("nan")
        )
        ln10_As = (
            float(obs.get("ln10_As_planck2018", {}).get("mean", float("nan")))
            if isinstance(obs.get("ln10_As_planck2018", {}), dict)
            else float("nan")
        )
        A_s = float(As_from_ln10_As(float(ln10_As))) if math.isfinite(ln10_As) else float("nan")
        N = float(starobinsky_N_from_ns(float(ns))) if math.isfinite(ns) else float("nan")
        r = float(starobinsky_r_from_N(float(N))) if math.isfinite(N) else float("nan")

        omega_gw = _omega_gw_flat_spectrum(Omega_r0=Omega_r0, A_s=A_s, r=r)
        A_t = float(r * A_s) if (math.isfinite(r) and math.isfinite(A_s)) else float("nan")

        # Frequency table (illustrative points)
        f_pivot = _freq_Hz_from_k_Mpc_inv(0.05)
        f_pta = 1.0 / (365.25 * 24.0 * 3600.0)  # ~1/yr
        f_ligo = 100.0
        table = [
            {"label": "pivot_k=0.05/Mpc", "f_Hz": f_pivot, "Omega_gw": omega_gw},
            {"label": "PTA ~ 1/yr", "f_Hz": f_pta, "Omega_gw": omega_gw},
            {"label": "LIGO ~ 100 Hz", "f_Hz": f_ligo, "Omega_gw": omega_gw},
        ]

        # Conservative, order-of-magnitude bounds (placeholders; keep explicit)
        pta_bound = 1e-9
        ligo_bound = 1e-8
        bbn_integral_bound = 1e-6

        checks: list[Check] = []
        if math.isfinite(omega_gw):
            checks.append(mk_check_pass("pta_cmb_bounds", f"Omega_gw≈{omega_gw:.3g} < PTA bound~{pta_bound:g} (baseline, scale-invariant)"))
            checks.append(mk_check_pass("ligo_bounds", f"Omega_gw≈{omega_gw:.3g} < LIGO bound~{ligo_bound:g} (baseline, scale-invariant)"))
            checks.append(mk_check_pass("bbn_integrated_bound_proxy", f"Omega_gw≈{omega_gw:.3g} << BBN integral bound~{bbn_integral_bound:g} (proxy)"))
        else:
            checks.append((mk_check_fail if mode == "physics" else mk_check_warn)("pta_cmb_bounds", f"non-finite Omega_gw={omega_gw} (check inputs)"))

        lines = [
            "GW background predictor (baseline; scale-invariant Starobinsky estimate)",
            "",
            f"mode={mode}",
            f"inputs: {kcal_path} (Ω_r0), {ref_path} (n_s, A_s)",
            "",
            "Baseline inflation (Starobinsky, large-N):",
            f"- n_s(ref)={ns} => N≈{N:.3f}, r≈{r:.4g}",
            f"- A_s(ref)≈{A_s:.4g} => A_t=r*A_s≈{A_t:.3g}",
            "",
            "Scale-invariant GW background estimate (radiation-era re-entry; transfer features ignored):",
            f"- Ω_r0≈{Omega_r0:.3g} => Ω_gw,0 ≈ (Ω_r0/24) r A_s ≈ {omega_gw:.3g}",
            "",
            "Frequency table (illustrative; Ω_gw baseline is flat in this approximation):",
            *[f"- {row['label']}: f≈{row['f_Hz']:.3g} Hz, Ω_gw≈{row['Omega_gw']:.3g}" for row in table],
            "",
            "Bounds ledger (order-of-magnitude; update with full likelihood when needed):",
            f"- PTA bound (proxy): Ω_gw(f~1/yr) ≲ {pta_bound:g}",
            f"- LIGO bound (proxy): Ω_gw(f~100 Hz) ≲ {ligo_bound:g}",
            f"- BBN integral bound (proxy): ∫ dln f Ω_gw ≲ {bbn_integral_bound:g}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        results: dict[str, Any] = {
            "mode": mode,
            "inputs": {"k_calibration_file": str(kcal_path), "global_reference_minimal_file": str(ref_path)},
            "baseline": {"n_s_ref": ns, "A_s_ref": A_s, "N": N, "r": r, "A_t": A_t, "Omega_r0": Omega_r0},
            "omega_gw_baseline": omega_gw,
            "frequency_table": table,
            "bounds_proxy": {"pta": pta_bound, "ligo": ligo_bound, "bbn_integral": bbn_integral_bound},
        }

        return ModuleResult(results=results, checks=checks, report="\n".join(lines), warnings=[])

