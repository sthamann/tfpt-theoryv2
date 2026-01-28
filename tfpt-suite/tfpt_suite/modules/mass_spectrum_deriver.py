from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.mobius_z3_yukawa_generator import _yukawa_ratios_from_mobius
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass, mk_check_warn


class MassSpectrumDeriverModule(TfptModule):
    module_id = "mass_spectrum_deriver"
    title = "Mass spectrum deriver (Mobius/Z3 ratios from δ⋆; lepton ratio checks)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants: tfpt_suite/constants.py (δ⋆ anchor)",
                "lepton masses: tfpt_suite/data/lepton_masses_pdg.json",
            ],
            outputs=[
                "Mobius-predicted mass hierarchy ratios from δ⋆ (dimensionless)",
                "measured lepton ratios and comparison to predictions",
            ],
            formulas=[
                r"δ⋆ = 3/5 + φ0/6",
                r"M_y(δ)=(y+δ)/(y-δ)",
                r"m_τ/m_μ ≈ M_1(δ)^2,  m_μ/m_e ≈ (M_1(δ)|M_{1/3}(δ)|)^2  (suite convention; see generator)",
            ],
            validation=[
                "Lepton ratio predictions from δ⋆ are within a declared tolerance of PDG ratios (diagnostic gate).",
                "No continuous fit parameters are introduced.",
            ],
            determinism="Deterministic given TFPT constants and PDG lepton masses table.",
            question="Do Mobius/Z3 mass-ratio formulas anchored at δ⋆ reproduce observed lepton hierarchy ratios at the level of order-of-magnitude closure?",
            objective=[
                "Turn mass-ratio claims into a machine-checkable module (lepton sector first).",
                "Keep scheme-dependence explicit: quark ratios are emitted but not gated (mass scheme/scale matters).",
            ],
            gaps=[
                "Full publication-grade mass derivation requires a unique holonomy/operator selection rule and a quark-mass scheme/scale reference table.",
            ],
        )

    def run(self, config) -> ModuleResult:
        lep_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        lep = json.loads(lep_path.read_text(encoding="utf-8"))

        me = float(lep["masses"]["electron"]["mean"])
        mm = float(lep["masses"]["muon"]["mean"])
        mt = float(lep["masses"]["tau"]["mean"])

        meas_tau_mu = float(mt / mm) if mm > 0 else float("nan")
        meas_mu_e = float(mm / me) if me > 0 else float("nan")

        c = TfptConstants.compute()
        delta_star = float(c.delta_star)
        ratios = _yukawa_ratios_from_mobius(delta_used=float(delta_star))

        pred_tau_mu = float(ratios.get("m_tau_over_m_mu", float("nan")))
        pred_mu_e = float(ratios.get("m_mu_over_m_e", float("nan")))

        def _rel_err(pred: float, meas: float) -> float:
            if not (pred == pred and meas == meas) or meas == 0:
                return float("nan")
            return float((pred / meas) - 1.0)

        rel_tau_mu = _rel_err(pred_tau_mu, meas_tau_mu)
        rel_mu_e = _rel_err(pred_mu_e, meas_mu_e)

        checks: list[Check] = []
        tol = 0.05  # 5% relative tolerance (diagnostic; not a publication-grade fit)
        ok_tau = (abs(rel_tau_mu) <= tol) if rel_tau_mu == rel_tau_mu else False
        ok_mu = (abs(rel_mu_e) <= tol) if rel_mu_e == rel_mu_e else False

        if ok_tau and ok_mu:
            checks.append(mk_check_pass("ratios_match_pdg", f"PASS: rel(tau/mu)={rel_tau_mu:.3e}, rel(mu/e)={rel_mu_e:.3e} (tol={tol})"))
        else:
            checks.append(
                (mk_check_warn if (ok_tau or ok_mu) else mk_check_fail)(
                    "ratios_match_pdg",
                    f"rel(tau/mu)={rel_tau_mu}, rel(mu/e)={rel_mu_e} (tol={tol})",
                )
            )

        # Always record the full ratio dictionary (includes quark ratios used by the generator).
        checks.append(mk_check_pass("mass_ratio_dictionary_present", f"{sorted(ratios.keys())}"))

        lines: list[str] = []
        lines += [
            "Mass spectrum deriver (Mobius/Z3 ratios from δ⋆)",
            "",
            f"delta_star (TFPT): {delta_star}",
            "",
            "Measured lepton ratios (from lepton_masses_pdg.json):",
            f"- m_tau/m_mu = {meas_tau_mu:.12g}",
            f"- m_mu/m_e   = {meas_mu_e:.12g}",
            "",
            "Predicted ratios from δ⋆ via Mobius map (suite convention):",
            f"- m_tau/m_mu(pred) = {pred_tau_mu:.12g}  (rel={rel_tau_mu:.3e})",
            f"- m_mu/m_e(pred)   = {pred_mu_e:.12g}  (rel={rel_mu_e:.3e})",
            "",
            "Additional ratios emitted (scheme/scale dependent; not gated):",
            f"- m_s/m_d(pred) = {float(ratios.get('m_s_over_m_d', float('nan'))):.12g}",
            f"- m_b/m_s(pred) = {float(ratios.get('m_b_over_m_s', float('nan'))):.12g}",
            f"- m_c/m_u(pred) = {float(ratios.get('m_c_over_m_u', float('nan'))):.12g}",
            f"- m_t/m_c(pred) = {float(ratios.get('m_t_over_m_c', float('nan'))):.12g}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "inputs": {"lepton_masses_file": str(lep_path)},
                "tfpt": {"delta_star": delta_star},
                "measured": {"m_tau_over_m_mu": meas_tau_mu, "m_mu_over_m_e": meas_mu_e},
                "predicted_from_delta_star": {
                    "ratios": ratios,
                    "m_tau_over_m_mu": pred_tau_mu,
                    "m_mu_over_m_e": pred_mu_e,
                    "rel_errors": {"tau_over_mu": rel_tau_mu, "mu_over_e": rel_mu_e},
                    "tolerance_rel": tol,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

