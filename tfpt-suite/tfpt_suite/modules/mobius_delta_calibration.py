from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _mobius_map(y: mp.mpf, delta: mp.mpf) -> mp.mpf:
    """
    Möbius ladder map used in the TFPT flavor notes / Paper v1.06:
      M_y(δ) = (y + δ) / (y - δ)
    """
    denom = y - delta
    if denom == 0:
        return mp.mpf("inf")
    return (y + delta) / denom


class MobiusDeltaCalibrationModule(TfptModule):
    module_id = "mobius_delta_calibration"
    title = "Möbius δ calibration (τ/μ) vs geometric anchor δ⋆ = 3/5 + varphi0/6"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariant varphi0 (for δ⋆)",
                "lepton masses (τ, μ): tfpt_suite/data/lepton_masses_pdg.json",
            ],
            outputs=[
                "δ_M from τ/μ",
                "δ⋆ from geometry (varphi0)",
                "percent deviation (δ⋆-δ_M)/δ_M",
                "Möbius map values at the canonical cusp set y ∈ {1, 1/3, 2/3}",
            ],
            formulas=[
                "δ_M = (sqrt(m_tau/m_mu) - 1) / (sqrt(m_tau/m_mu) + 1)",
                "δ⋆ = 3/5 + varphi0/6",
                "M_y(δ) = (y + δ)/(y - δ), y ∈ {1, 1/3, 2/3}",
            ],
            validation=[
                "δ_M and δ⋆ are in (0,1) and agree at the sub-percent level (historically ~0.1–0.2%).",
            ],
            determinism="Deterministic (fixed data inputs).",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        data_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        raw = json.loads(data_path.read_text(encoding="utf-8"))
        masses = raw.get("masses", {})
        mtau = mp.mpf(str(masses["tau"]["mean"]))
        mmu = mp.mpf(str(masses["muon"]["mean"]))

        R = mp.sqrt(mtau / mmu)
        delta_M = (R - 1) / (R + 1)
        delta_star = mp.mpf(c.delta_star)

        deviation_percent = (delta_star - delta_M) / delta_M * 100

        cusps = [mp.mpf(1), mp.mpf(1) / 3, mp.mpf(2) / 3]
        mobius_vals = {str(y): _mobius_map(y, delta_M) for y in cusps}

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="delta_M_in_unit_interval",
                passed=bool(delta_M > 0 and delta_M < 1),
                detail=f"delta_M={delta_M}",
            )
        )
        checks.append(
            Check(
                check_id="delta_star_in_unit_interval",
                passed=bool(delta_star > 0 and delta_star < 1),
                detail=f"delta_star={delta_star}",
            )
        )
        checks.append(
            Check(
                check_id="delta_star_close_to_delta_M",
                passed=bool(mp.fabs(deviation_percent) < 1),
                detail=f"deviation_percent={deviation_percent}",
            )
        )

        lines: list[str] = []
        lines += [
            "Möbius δ calibration (τ/μ) vs geometric δ⋆",
            "",
            f"inputs file: {data_path}",
            "",
            "Definitions:",
            "- δ_M = (sqrt(m_tau/m_mu) - 1) / (sqrt(m_tau/m_mu) + 1)",
            "- δ⋆  = 3/5 + varphi0/6",
            "- M_y(δ) = (y+δ)/(y-δ)",
            "",
            "Lepton masses (GeV):",
            f"- m_tau = {mtau}",
            f"- m_mu  = {mmu}",
            "",
            f"R = sqrt(m_tau/m_mu) = {R}",
            f"δ_M   = {delta_M}",
            f"δ⋆    = {delta_star}",
            f"(δ⋆-δ_M)/δ_M = {deviation_percent} %",
            "",
            "Canonical cusp set (from hypercharge holonomy): y ∈ {1, 1/3, 2/3}",
            "Möbius map values at δ_M:",
        ]
        for y in cusps:
            lines.append(f"- M_{y}(δ_M) = {mobius_vals[str(y)]}")

        lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
        ]

        return ModuleResult(
            results={
                "inputs": {
                    "lepton_masses_file": str(data_path),
                    "mtau_GeV": mtau,
                    "mmu_GeV": mmu,
                },
                "delta": {
                    "delta_M_from_tau_mu": delta_M,
                    "delta_star_from_varphi0": delta_star,
                    "deviation_percent": deviation_percent,
                },
                "cusps": [str(y) for y in cusps],
                "mobius_map_at_delta_M": {k: v for k, v in mobius_vals.items()},
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

