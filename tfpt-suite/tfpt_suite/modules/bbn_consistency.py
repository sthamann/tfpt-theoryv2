from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

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


def _eta10_from_omega_b_h2(omega_b_h2: float) -> float:
    # Standard conversion: η10 := 10^10 (n_b/n_γ) ≈ 273.9 · Ω_b h^2 (for T_CMB≈2.725 K).
    return float(273.9 * float(omega_b_h2))


def _bbn_fit_Yp(*, eta10: float, N_eff: float) -> float:
    """
    Engineering-level fit for primordial helium fraction.

    Form used as a compact proxy (not a full network):
      Y_p ≈ 0.2471 + 0.0016(η10 − 6) + 0.013(ΔN_eff)
    """
    dNeff = float(N_eff) - 3.046
    return float(0.2471 + 0.0016 * (float(eta10) - 6.0) + 0.013 * dNeff)


def _bbn_fit_D_over_H(*, eta10: float, N_eff: float) -> float:
    """
    Engineering-level fit for primordial deuterium abundance.

    Proxy scaling:
      (D/H) ≈ 2.54e-5 · (6/η10)^1.6 · (1 + 0.135 ΔN_eff)
    """
    dNeff = float(N_eff) - 3.046
    if eta10 <= 0:
        return float("nan")
    base = float(2.54e-5) * (6.0 / float(eta10)) ** 1.6
    return float(base * (1.0 + 0.135 * dNeff))


def _z_score(pred: float, mean: float, sigma_exp: float, sigma_theory: float) -> float:
    s2 = float(sigma_exp) ** 2 + float(sigma_theory) ** 2
    if s2 <= 0 or not math.isfinite(s2):
        return float("nan")
    return float((float(pred) - float(mean)) / math.sqrt(s2))


class BbnConsistencyModule(TfptModule):
    module_id = "bbn_consistency"
    title = "BBN consistency (light elements; engineering-level fit + reference-table check)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "Planck reference: tfpt_suite/data/global_reference_minimal.json (Ω_b h^2)",
                "BBN reference: tfpt_suite/data/bbn_reference.json (Y_p, D/H, N_eff anchor values)",
                "(optional) thermodynamic ledger: bbn_neff_sanity (g_*(T), N_eff marker)",
            ],
            outputs=[
                "η10 from Ω_b h^2",
                "engineering-level BBN fit predictions for (Y_p, D/H) and a reference-table z-score ledger",
                "light_elements_match gate (PASS/WARN/FAIL by mode)",
            ],
            formulas=[
                r"\eta_{10} \approx 273.9\,\Omega_b h^2",
                r"Y_p \approx 0.2471 + 0.0016(\eta_{10}-6) + 0.013(\Delta N_{\rm eff}) \;\; (\text{proxy fit})",
                r"(D/H) \approx 2.54\times10^{-5}\,(6/\eta_{10})^{1.6}\,(1+0.135\Delta N_{\rm eff}) \;\; (\text{proxy fit})",
            ],
            validation=[
                "Computes finite η10, Y_p, and D/H under the declared reference inputs.",
                "Provides a falsifiable reference-table check (with an explicit theory-floor for the proxy-fit uncertainty).",
            ],
            determinism="Deterministic given inputs.",
            question="Are BBN light-element abundances consistent with the declared Ω_b h^2 and N_eff under a transparent proxy model?",
            objective=[
                "Close the 'BBN missing' placeholder by providing an explicit, auditable light-element consistency module.",
                "Keep publication-grade scope explicit: this is a proxy-fit, not a full nuclear network.",
            ],
            gaps=[
                "A publication-grade BBN module would integrate weak rates and the full nuclear network and compare to a dedicated likelihood for (D/H, Y_p, N_eff).",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        data_dir = Path(__file__).resolve().parent.parent / "data"
        ref_min_path = data_dir / "global_reference_minimal.json"
        bbn_ref_path = data_dir / "bbn_reference.json"

        ref_min = json.loads(ref_min_path.read_text(encoding="utf-8")) if ref_min_path.is_file() else {}
        bbn_ref = json.loads(bbn_ref_path.read_text(encoding="utf-8")) if bbn_ref_path.is_file() else {}

        obs_min = ref_min.get("observables", {}) if isinstance(ref_min.get("observables", {}), dict) else {}
        omega_b = obs_min.get("omega_b_h2_planck2018", {}) if isinstance(obs_min.get("omega_b_h2_planck2018", {}), dict) else {}
        omega_b_mean = float(omega_b.get("mean", float("nan")))
        omega_b_sigma = float(omega_b.get("sigma", float("nan")))

        obs_bbn = bbn_ref.get("observables", {}) if isinstance(bbn_ref.get("observables", {}), dict) else {}
        Yp_ref = obs_bbn.get("Y_p_primordial", {}) if isinstance(obs_bbn.get("Y_p_primordial", {}), dict) else {}
        DH_ref = obs_bbn.get("D_over_H_primordial", {}) if isinstance(obs_bbn.get("D_over_H_primordial", {}), dict) else {}
        Neff_ref = obs_bbn.get("N_eff_BBN", {}) if isinstance(obs_bbn.get("N_eff_BBN", {}), dict) else {}

        # Use SM marker for N_eff unless bbn_neff_sanity output exists in the run directory.
        neff_marker = 3.046
        neff_out = Path(config.output_dir) / "bbn_neff_sanity" / "results.json"
        if neff_out.is_file():
            try:
                payload = json.loads(neff_out.read_text(encoding="utf-8"))
                neff_marker = float(payload.get("results", {}).get("neff_marker", neff_marker))
            except Exception:
                neff_marker = 3.046

        eta10 = _eta10_from_omega_b_h2(omega_b_mean) if math.isfinite(omega_b_mean) else float("nan")
        Yp_pred = _bbn_fit_Yp(eta10=eta10, N_eff=neff_marker) if math.isfinite(eta10) else float("nan")
        DH_pred = _bbn_fit_D_over_H(eta10=eta10, N_eff=neff_marker) if math.isfinite(eta10) else float("nan")

        # Proxy-model theory floors (explicit): these reflect the fact that we are not running a full BBN network here.
        sigma_theory_Yp = 0.002
        sigma_theory_DH = 1.0e-6
        sigma_theory_Neff = 0.0

        # z-scores vs reference table
        z_Yp = _z_score(Yp_pred, float(Yp_ref.get("mean", float("nan"))), float(Yp_ref.get("sigma_exp", float("nan"))), sigma_theory_Yp)
        z_DH = _z_score(DH_pred, float(DH_ref.get("mean", float("nan"))), float(DH_ref.get("sigma_exp", float("nan"))), sigma_theory_DH)
        z_Neff = _z_score(neff_marker, float(Neff_ref.get("mean", float("nan"))), float(Neff_ref.get("sigma_exp", float("nan"))), sigma_theory_Neff)

        checks: list[Check] = []
        checks.append(mk_check_pass("eta10_computed", f"eta10≈{eta10:.3g} from Ω_b h^2={omega_b_mean:.6g}"))
        checks.append(mk_check_info("neff_used", f"N_eff={neff_marker} (bbn_neff_sanity marker if available)"))

        def _gate(check_id: str, z: float, pred: float, mean: float, sigma_exp: float, sigma_theory: float) -> Check:
            if not math.isfinite(z):
                return (mk_check_fail if mode == "physics" else mk_check_warn)(
                    check_id, f"non-finite z (pred={pred}, mean={mean}, sigma_exp={sigma_exp}, sigma_theory={sigma_theory})"
                )
            if abs(z) <= 2.0:
                return mk_check_pass(check_id, f"|z|={abs(z):.3g} (pred={pred}, mean={mean}, sigma_eff≈{math.sqrt(sigma_exp**2+sigma_theory**2):.3g})")
            if abs(z) <= 5.0:
                return mk_check_warn(check_id, f"|z|={abs(z):.3g} (pred={pred}, mean={mean}, sigma_eff≈{math.sqrt(sigma_exp**2+sigma_theory**2):.3g})")
            return (mk_check_fail if mode == "physics" else mk_check_warn)(
                check_id, f"|z|={abs(z):.3g} (pred={pred}, mean={mean}, sigma_eff≈{math.sqrt(sigma_exp**2+sigma_theory**2):.3g})"
            )

        checks.append(
            _gate(
                "Yp_within_2sigma",
                z_Yp,
                Yp_pred,
                float(Yp_ref.get("mean", float("nan"))),
                float(Yp_ref.get("sigma_exp", float("nan"))),
                sigma_theory_Yp,
            )
        )
        checks.append(
            _gate(
                "D_over_H_within_2sigma",
                z_DH,
                DH_pred,
                float(DH_ref.get("mean", float("nan"))),
                float(DH_ref.get("sigma_exp", float("nan"))),
                sigma_theory_DH,
            )
        )
        checks.append(
            _gate(
                "Neff_within_2sigma",
                z_Neff,
                neff_marker,
                float(Neff_ref.get("mean", float("nan"))),
                float(Neff_ref.get("sigma_exp", float("nan"))),
                sigma_theory_Neff,
            )
        )

        # Aggregate gate
        ok = all(c.passed for c in checks if c.check_id in {"Yp_within_2sigma", "D_over_H_within_2sigma", "Neff_within_2sigma"})
        checks.append(
            mk_check_pass("light_elements_match", "Yp, D/H, N_eff all within the declared proxy tolerances")
            if ok
            else (mk_check_fail if mode == "physics" else mk_check_warn)("light_elements_match", "one or more BBN observables outside proxy tolerances")
        )

        lines = [
            "BBN consistency (engineering-level fit; not a full nuclear network)",
            "",
            f"mode={mode}",
            f"reference: {ref_min_path} (Ω_b h^2), {bbn_ref_path} (Yp, D/H, N_eff)",
            "",
            "Inputs:",
            f"- Ω_b h^2 (Planck 2018): {omega_b_mean} ± {omega_b_sigma}",
            f"- N_eff (marker): {neff_marker}",
            "",
            "Derived:",
            f"- η10 ≈ {eta10:.6g}",
            "",
            "Proxy BBN predictions:",
            f"- Y_p(pred) ≈ {Yp_pred:.6g}",
            f"- D/H(pred) ≈ {DH_pred:.6g}",
            "",
            "Theory floors (proxy-fit uncertainty; explicit):",
            f"- sigma_theory(Y_p) = {sigma_theory_Yp}",
            f"- sigma_theory(D/H) = {sigma_theory_DH}",
            f"- sigma_theory(N_eff) = {sigma_theory_Neff}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        results: dict[str, Any] = {
            "mode": mode,
            "inputs": {
                "global_reference_minimal_file": str(ref_min_path),
                "bbn_reference_file": str(bbn_ref_path),
                "omega_b_h2": {"mean": omega_b_mean, "sigma": omega_b_sigma},
                "N_eff": neff_marker,
            },
            "derived": {"eta10": eta10},
            "predictions": {"Y_p": Yp_pred, "D_over_H": DH_pred},
            "theory_floors": {"Y_p": sigma_theory_Yp, "D_over_H": sigma_theory_DH, "N_eff": sigma_theory_Neff},
            "z_scores": {"Y_p": z_Yp, "D_over_H": z_DH, "N_eff": z_Neff},
        }

        return ModuleResult(results=results, checks=checks, report="\n".join(lines), warnings=[])

