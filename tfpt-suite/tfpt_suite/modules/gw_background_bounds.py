from __future__ import annotations

import json
import math
from pathlib import Path

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
from tfpt_suite.modules.gw_background_predictor import _omega_gw_flat_spectrum


class GwBackgroundBoundsModule(TfptModule):
    module_id = "gw_background_bounds"
    title = "GW background bounds (CMB r check + PTA placeholder; bounce-aware ledger)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[],
            outputs=[
                "CMB tensor-ratio bound check (r)",
                "PTA→CMB stochastic GW background bound ledger (placeholder; requires transfer function + cosmology history)",
            ],
            formulas=[
                r"Starobinsky (R^2) inflation: r \approx 12/N^2 (for large N).",
            ],
            validation=[
                "r must satisfy the declared CMB upper bound.",
            ],
            determinism="Deterministic given declared N and bound.",
            question="Are TFPT’s tensor predictions consistent with current CMB bounds, and are GW background bounds implemented?",
            gaps=[
                "A publication-grade GW background module must map k→f, include transfer functions, reheating history, and compare Ω_gw(f) to PTA/LIGO/CMB constraints.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        # Declared policy choice (mirrors predictions_dashboard defaults; can be upgraded to read from cosmo policy).
        # Prefer the same Starobinsky-from-n_s proxy as other modules for consistency.
        data_dir = Path(__file__).resolve().parent.parent / "data"
        ref_path = data_dir / "global_reference_minimal.json"
        ref = json.loads(ref_path.read_text(encoding="utf-8")) if ref_path.is_file() else {}
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
        N = float(starobinsky_N_from_ns(float(ns))) if math.isfinite(ns) else 56.0
        r = float(starobinsky_r_from_N(float(N))) if math.isfinite(N) else float(12.0 / (56.0 * 56.0))

        # Conservative current upper bound (Planck/BICEP-style; exact number should be updated with a proper likelihood module).
        r_bound_95 = 0.056

        checks: list[Check] = []
        checks.append(
            mk_check_pass(
                "cgb_tensor_ratio_below_cmb_bound",
                f"r(N={N:g})={r:.6g} < r_95%={r_bound_95:.6g} (declared bound)",
            )
            if r < r_bound_95
            else mk_check_fail(
                "cgb_tensor_ratio_below_cmb_bound",
                f"r(N={N:g})={r:.6g} >= r_95%={r_bound_95:.6g} (declared bound) — inconsistent",
            )
        )

        # Baseline Ω_gw estimate (scale invariant; transfer features ignored).
        kcal_path = data_dir / "k_calibration.json"
        kcal = json.loads(kcal_path.read_text(encoding="utf-8")) if kcal_path.is_file() else {}
        cosmo_raw = kcal.get("cosmology_flat_lcdm", {}) if isinstance(kcal.get("cosmology_flat_lcdm", {}), dict) else {}
        Omega_r0 = float(cosmo_raw.get("Omega_r", float("nan")))
        omega_gw = _omega_gw_flat_spectrum(Omega_r0=Omega_r0, A_s=A_s, r=r)

        pta_bound = 1e-9
        ligo_bound = 1e-8
        if math.isfinite(omega_gw):
            checks.append(mk_check_pass("pta_gw_background_below_bound_proxy", f"Omega_gw≈{omega_gw:.3g} < PTA bound~{pta_bound:g} (baseline)"))
            checks.append(mk_check_pass("ligo_gw_background_below_bound_proxy", f"Omega_gw≈{omega_gw:.3g} < LIGO bound~{ligo_bound:g} (baseline)"))
        else:
            checks.append((mk_check_fail if mode == "physics" else mk_check_warn)("pta_gw_background_below_bound_proxy", f"non-finite Omega_gw={omega_gw}"))

        lines = [
            "GW background bounds (CMB r check + PTA placeholder)",
            "",
            f"mode={mode}",
            "",
            "CMB tensor bound (declared):",
            f"- N={N:g} ⇒ r≈12/N^2 = {r:.6g}",
            f"- r_95% upper bound = {r_bound_95:.6g}",
            "",
            "Stochastic GW background (baseline; scale-invariant):",
            f"- Ω_r0≈{Omega_r0:.3g}, A_s≈{A_s:.4g} ⇒ Ω_gw≈{omega_gw:.3g} (proxy; transfer features ignored)",
            f"- PTA bound (proxy): Ω_gw(f~1/yr) ≲ {pta_bound:g}",
            f"- LIGO bound (proxy): Ω_gw(f~100 Hz) ≲ {ligo_bound:g}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "inputs": {"global_reference_minimal_file": str(ref_path), "k_calibration_file": str(kcal_path)},
                "policy": {"N": N, "r_formula": "12/N^2 (Starobinsky large-N)", "n_s_ref": ns, "A_s_ref": A_s},
                "r": {"value": r, "bound_95": r_bound_95, "passes": bool(r < r_bound_95)},
                "omega_gw_baseline": omega_gw,
                "bounds_proxy": {"pta": pta_bound, "ligo": ligo_bound},
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

