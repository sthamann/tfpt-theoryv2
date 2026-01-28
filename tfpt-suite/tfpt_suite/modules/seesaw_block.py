from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.e8_ladder import e8_phi_n
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


@dataclass(frozen=True)
class SeesawReference:
    MR_GeV: mp.mpf
    mnu3_eV: mp.mpf
    dm31_sq_eV2: mp.mpf


class SeesawBlockModule(TfptModule):
    module_id = "seesaw_block"
    title = "Seesaw block (paper v1.06 anchor): MR scale + m_nu3 order-of-magnitude"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants (varphi0, gamma(0), lambda) from tfpt_suite/constants.py",
                "SM v reference from tfpt_suite/data/sm_inputs_mz.json",
                "paper v1.06 anchor numbers (MR ~ 1.311e15 GeV, mnu3 ~ 0.048 eV)",
            ],
            outputs=[
                "dimensionless ladder step varphi_5",
                "inferred ζ_NR (from MR = ζ_NR M_Pl varphi_5)",
                "m_nu3 and Δm31^2 estimates for y_nu3 ~ 1 using v^2/MR",
            ],
            formulas=[
                "varphi_5 = varphi0 * exp(-gamma(0)) * (D_5/D_1)^lambda, D_n=60-2n, D_1=58",
                "MR = ζ_NR * M_Pl * varphi_5",
                "m_nu3 ~ v^2 / MR  (paper v1.06 uses this normalization for y_{ν3}~1)",
                "Δm31^2 ~ m_nu3^2 (order-of-magnitude anchor)",
            ],
            validation=[
                "Produces the paper-scale MR and m_nu3 order-of-magnitude without any fit parameters (beyond the stated anchors).",
            ],
            determinism="Deterministic (fixed inputs).",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        v_sm = mp.mpf(str(sm_raw.get("v_ev_GeV", 246.0)))

        # Paper v1.06 anchors (Seesaw block n=5):
        # - MR ≈ 1.311×10^15 GeV
        # - mν3 ≈ 0.04807 eV (using v_H≈251 GeV in that document)
        ref = SeesawReference(
            MR_GeV=mp.mpf("1.311e15"),
            mnu3_eV=mp.mpf("0.04807"),
            dm31_sq_eV2=mp.mpf("2.31e-3"),
        )

        # We use the (unreduced) Planck mass as in the paper's PyR@TE block.
        MPl_GeV = mp.mpf("1.221e19")

        varphi5 = e8_phi_n(n=5, c=c)
        zeta_nr = ref.MR_GeV / (MPl_GeV * varphi5)

        # Convert GeV -> eV
        GeV_to_eV = mp.mpf("1e9")

        # neutrino mass estimate using the suite's v reference
        mnu3_sm_eV = (v_sm**2) / ref.MR_GeV * GeV_to_eV
        dm31_sm_eV2 = mnu3_sm_eV**2

        # Optional: reproduce the paper's stated mν3 by using its implicit v_H ≈ sqrt(mν3 * MR) (in GeV).
        v_from_paper_GeV = mp.sqrt((ref.mnu3_eV / GeV_to_eV) * ref.MR_GeV)

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="varphi5_positive",
                passed=bool(varphi5 > 0),
                detail=f"varphi5={varphi5}",
            )
        )
        checks.append(
            Check(
                check_id="zeta_nr_finite",
                passed=bool(mp.isfinite(zeta_nr)),
                detail=f"zeta_NR={zeta_nr}",
            )
        )
        checks.append(
            Check(
                check_id="mnu3_order_of_magnitude_ok",
                passed=bool(mnu3_sm_eV > mp.mpf("0.01") and mnu3_sm_eV < mp.mpf("0.2")),
                detail=f"mnu3_sm_eV={mnu3_sm_eV}",
            )
        )

        lines: list[str] = []
        lines += [
            "Seesaw block (paper v1.06 anchor, n=5)",
            "",
            "Inputs:",
            f"- SM v reference (from {sm_path}): v = {v_sm} GeV",
            f"- TFPT invariants: varphi0={c.varphi0}, gamma(0)={c.gamma0}, lambda={c.e8_lambda}",
            "",
            "E8 ladder step (dimensionless):",
            f"- varphi_5 = {varphi5}",
            "",
            "Paper v1.06 anchors:",
            f"- MR (target) = {ref.MR_GeV} GeV",
            f"- mnu3 (paper) = {ref.mnu3_eV} eV, Δm31^2 (paper) = {ref.dm31_sq_eV2} eV^2",
            "",
            "Inferred block calibration (from MR = ζ_NR M_Pl varphi_5):",
            f"- M_Pl = {MPl_GeV} GeV (unreduced)",
            f"- ζ_NR = {zeta_nr}",
            "",
            "Neutrino mass estimate (y_{ν3}~1, using mν3 ~ v^2/MR):",
            f"- mnu3(v={v_sm} GeV) = {mnu3_sm_eV} eV",
            f"- Δm31^2 ≈ mnu3^2 = {dm31_sm_eV2} eV^2",
            "",
            "Implied Higgs VEV from the paper's mnu3 and MR (v = sqrt(mnu3 * MR)):",
            f"- v_from_paper = {v_from_paper_GeV} GeV",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module is an anchor-level consistency block: it does not yet implement a full 3×3 Y_ν texture, RG running, or threshold matching.",
            "- Once Y_ν and thresholds are implemented (mt→UV-only policy), this block should be upgraded to a full PMNS+mass-spectrum pipeline.",
            "",
        ]

        return ModuleResult(
            results={
                "inputs": {
                    "sm_inputs_mz_file": str(sm_path),
                    "v_sm_GeV": v_sm,
                    "MPl_GeV": MPl_GeV,
                },
                "e8_ladder": {
                    "varphi5": varphi5,
                    "gamma0": c.gamma0,
                    "lambda": c.e8_lambda,
                    "definition": "varphi_n = varphi0 * exp(-gamma0) * (D_n/D_1)^lambda for n>=1; D_n=60-2n, D_1=58",
                },
                "paper_v1_06_anchor": {
                    "MR_GeV": ref.MR_GeV,
                    "mnu3_eV": ref.mnu3_eV,
                    "dm31_sq_eV2": ref.dm31_sq_eV2,
                },
                "derived": {
                    "zeta_NR": zeta_nr,
                    "mnu3_from_v_sm_eV": mnu3_sm_eV,
                    "dm31_sq_from_v_sm_eV2": dm31_sm_eV2,
                    "v_from_paper_GeV": v_from_paper_GeV,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

