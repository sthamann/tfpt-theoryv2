from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from tfpt_suite.constants import TfptConstants
from tfpt_suite.cosmo_scale_map import MPL_REDUCED_GEV
from tfpt_suite.module_base import (
    Check,
    ModuleResult,
    ModuleSpec,
    TfptModule,
    mk_check_info,
    mk_check_pass,
    mk_check_warn,
)


class G2AndLambShiftProxyModule(TfptModule):
    module_id = "g2_and_lamb_shift_proxy"
    title = "g-2 / Lamb shift proxy (TFPT-scale new-physics contribution is negligible; consistency ledger)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT scale M from constants (R^2 scale M/Mpl)",
                "lepton masses: tfpt_suite/data/lepton_masses_pdg.json",
            ],
            outputs=[
                "dimensionless new-physics scaling proxy δa_ℓ ~ (m_ℓ/M)^2 for ℓ=e,μ",
                "Lamb-shift scale proxy δν ~ (m_e^3/M^2)·(GeV→Hz) (order-of-magnitude)",
                "consistency gates against conservative anomaly scales (does not claim to *explain* anomalies)",
            ],
            formulas=[
                r"\delta a_\ell \sim (m_\ell/M)^2 \;\; (\text{dimensional-analysis proxy})",
                r"\delta E \sim m_e^3/M^2,\;\; \delta\nu = \delta E\cdot(\mathrm{GeV}\to\mathrm{Hz}) \;\; (\text{Lamb-shift proxy})",
            ],
            validation=[
                "Shows that TFPT-scale suppressed contributions are far below current g-2/Lamb-shift anomaly scales (consistency; not an explanation).",
            ],
            determinism="Deterministic given inputs.",
            question="Are high-precision QED observables (g-2, Lamb shift) audited as part of the ToE closure?",
            objective=[
                "Provide a minimal, falsifiable consistency check: TFPT-scale suppressed contributions must not overshoot known anomaly scales.",
            ],
            gaps=[
                "This is not a full QED calculation. Publication-grade requires bound-state QED + hadronic inputs + a declared likelihood.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        data_dir = Path(__file__).resolve().parent.parent / "data"
        lep_path = data_dir / "lepton_masses_pdg.json"
        lep = json.loads(lep_path.read_text(encoding="utf-8")) if lep_path.is_file() else {}
        masses = lep.get("masses", {}) if isinstance(lep.get("masses", {}), dict) else {}
        m_e = float((masses.get("electron", {}) if isinstance(masses.get("electron", {}), dict) else {}).get("mean", float("nan")))
        m_mu = float((masses.get("muon", {}) if isinstance(masses.get("muon", {}), dict) else {}).get("mean", float("nan")))

        c = TfptConstants.compute()
        M_eff = float(c.M_over_Mpl) * float(MPL_REDUCED_GEV)

        # New-physics scaling proxies (dimension analysis; assumes O(1) coupling at scale M_eff).
        da_e = float((m_e / M_eff) ** 2) if (math.isfinite(m_e) and math.isfinite(M_eff) and M_eff > 0) else float("nan")
        da_mu = float((m_mu / M_eff) ** 2) if (math.isfinite(m_mu) and math.isfinite(M_eff) and M_eff > 0) else float("nan")

        # Lamb-shift frequency proxy: δE ~ m_e^3/M^2, δν = δE·(GeV→Hz)
        GeV_to_Hz = 2.417989242e23
        dE_lamb_GeV = float((m_e**3) / (M_eff**2)) if (math.isfinite(m_e) and math.isfinite(M_eff) and M_eff > 0) else float("nan")
        dnu_lamb_Hz = float(dE_lamb_GeV * GeV_to_Hz) if math.isfinite(dE_lamb_GeV) else float("nan")

        # Conservative anomaly scales (order-of-magnitude): used only as "do not overshoot" gates.
        da_mu_scale = 1e-9
        da_e_scale = 1e-12
        dnu_lamb_scale_Hz = 1.0  # 1 Hz is an extremely conservative floor for atomic precision

        checks: list[Check] = []
        checks.append(mk_check_info("tfpt_scale", f"M_eff≈{M_eff:.3g} GeV (from M/Mpl={c.M_over_Mpl})"))
        checks.append(mk_check_info("np_scaling", f"δa_e≈{da_e:.3g}, δa_μ≈{da_mu:.3g}, δν_Lamb≈{dnu_lamb_Hz:.3g} Hz (proxy)"))

        ok = (
            math.isfinite(da_e)
            and math.isfinite(da_mu)
            and math.isfinite(dnu_lamb_Hz)
            and (da_e < da_e_scale)
            and (da_mu < da_mu_scale)
            and (dnu_lamb_Hz < dnu_lamb_scale_Hz)
        )
        if ok:
            checks.append(
                mk_check_pass(
                    "precision_qed_consistency",
                    f"TFPT-scale NP is negligible: δa_e<{da_e_scale:g}, δa_μ<{da_mu_scale:g}, δν_Lamb<{dnu_lamb_scale_Hz:g} Hz",
                )
            )
        else:
            # In engineering mode keep this as WARN (this is a proxy); in physics mode treat as WARN as well
            # because the failure would indicate an internal inconsistency rather than a precise QED claim.
            checks.append(
                mk_check_warn(
                    "precision_qed_consistency",
                    f"non-finite or unexpectedly large proxy (δa_e={da_e}, δa_mu={da_mu}, δν={dnu_lamb_Hz})",
                )
            )

        lines = [
            "g-2 / Lamb shift proxy",
            "",
            f"mode={mode}",
            "",
            "TFPT-scale suppressed new-physics proxy (consistency only; no full QED calculation):",
            f"- M_eff ≈ {M_eff:.6g} GeV",
            f"- δa_e ~ (m_e/M)^2 ≈ {da_e:.3e}",
            f"- δa_μ ~ (m_μ/M)^2 ≈ {da_mu:.3e}",
            f"- δν_Lamb ~ (m_e^3/M^2)·(GeV→Hz) ≈ {dnu_lamb_Hz:.3e} Hz",
            "",
            "Conservative anomaly-scale gates (do-not-overshoot):",
            f"- δa_e < {da_e_scale:g}",
            f"- δa_μ < {da_mu_scale:g}",
            f"- δν_Lamb < {dnu_lamb_scale_Hz:g} Hz",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "inputs": {"lepton_masses_file": str(lep_path)},
                "tfpt": {"M_eff_GeV": M_eff},
                "proxies": {"delta_a_e": da_e, "delta_a_mu": da_mu, "delta_nu_lamb_Hz": dnu_lamb_Hz},
                "bounds_proxy": {"delta_a_e_scale": da_e_scale, "delta_a_mu_scale": da_mu_scale, "delta_nu_lamb_scale_Hz": dnu_lamb_scale_Hz},
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

