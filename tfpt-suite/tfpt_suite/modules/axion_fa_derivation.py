from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.e8_ladder import e8_phi_n
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass, mk_check_warn


PQ_BLOCK_N = 10
PQ_BLOCK_RANK = mp.mpf(1)
PQ_BLOCK_I1 = mp.mpf(1) / mp.mpf(3)
MPL_UNREDUCED_GEV = mp.mpf("1.221e19")
FA_REL_TOL = mp.mpf("5e-3")
MA_REL_TOL = mp.mpf("5e-2")
MICROEV_PER_GEV = mp.mpf("1e15")


def _zeta_block(*, c3: mp.mpf, r_B: mp.mpf, I1: mp.mpf) -> mp.mpf:
    """
    Block factor (paper v2.4/v2.5 block-constants definition):
      zeta_B = (pi c3) * exp(-beta_B pi c3) * exp(-k_B / c3)
      beta_B = (8 - r_B)/8,  k_B = (3/2) I1
    """
    beta_B = (mp.mpf(8) - r_B) / mp.mpf(8)
    k_B = mp.mpf(3) / mp.mpf(2) * I1
    return (mp.pi * c3) * mp.e ** (-beta_B * mp.pi * c3) * mp.e ** (-k_B / c3)


def _axion_mass_micro_eV(*, f_a_GeV: mp.mpf, m_pi_GeV: mp.mpf, f_pi_GeV: mp.mpf, m_u_GeV: mp.mpf, m_d_GeV: mp.mpf) -> mp.mpf:
    """
    QCD axion mass (standard):
      m_a = (m_pi f_pi / f_a) * sqrt(m_u m_d / (m_u + m_d)^2)
    """
    if f_a_GeV <= 0:
        return mp.mpf("nan")
    mass = (m_pi_GeV * f_pi_GeV / f_a_GeV) * mp.sqrt((m_u_GeV * m_d_GeV) / ((m_u_GeV + m_d_GeV) ** 2))
    return mass * MICROEV_PER_GEV


class AxionFaDerivationModule(TfptModule):
    module_id = "axion_fa_derivation"
    title = "Axion PQ block: derive f_a from E8 ladder + block constants"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants: tfpt_suite/constants.py (c3, varphi0, gamma0, lambda)",
                "axion claim: tfpt_suite/data/axion_tfpt_v106.json (block stage n=10, quoted f_a/m_a for cross-check)",
            ],
            outputs=[
                "E8 ladder varphi_n at n=10",
                "block factor zeta_PQ from (r_B, I1) and c3",
                "derived f_a and implied m_a",
            ],
            formulas=[
                "varphi_n = varphi0 * exp(-gamma0) * (D_n/D_1)^lambda, D_n=60-2n",
                "zeta_B = (pi c3) * exp(-beta_B pi c3) * exp(-k_B/c3), beta_B=(8-r_B)/8, k_B=(3/2)I1",
                "f_a = zeta_PQ * M_Pl * varphi_10 (unreduced Planck mass)",
                "m_a = (m_pi f_pi / f_a) * sqrt(m_u m_d / (m_u + m_d)^2)",
            ],
            validation=[
                "f_a_derived_not_quoted: derived f_a matches the quoted TFPT benchmark within tolerance.",
                "m_a_follows_from_f_a: derived m_a matches the quoted micro-eV value within tolerance.",
            ],
            determinism="Deterministic given the canonical TFPT inputs.",
            question="Can the PQ breaking scale be derived from the E8 ladder + block constants without quoting f_a?",
            objective=[
                "Compute f_a from the ladder/block rule at n=10 (no ad-hoc fit parameters).",
                "Validate that the implied m_a matches the quoted TFPT benchmark.",
            ],
            references=[
                "tfpt-theory-fullv24.tex (Block constants definition; PQ block n=10)",
                "tfpt_suite/data/axion_tfpt_v106.json (quoted benchmark for cross-check)",
            ],
            maturity="derivation (deterministic; cross-checked against quoted benchmark)",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        data_dir = Path(__file__).resolve().parent.parent / "data"
        ax_path = data_dir / "axion_tfpt_v106.json"
        ax = json.loads(ax_path.read_text(encoding="utf-8")) if ax_path.is_file() else {}

        claim = ax.get("tfpt_claim", {}) if isinstance(ax.get("tfpt_claim", {}), dict) else {}
        claim_n = int(claim.get("block_stage_n", PQ_BLOCK_N))
        claim_fa = mp.mpf(str(claim.get("f_a_GeV", "nan")))
        claim_ma = mp.mpf(str(claim.get("m_a_micro_eV", "nan")))

        qcd = ax.get("qcd_inputs", {}) if isinstance(ax.get("qcd_inputs", {}), dict) else {}
        m_pi = mp.mpf(str(qcd.get("m_pi_GeV", "0.1349768")))
        f_pi = mp.mpf(str(qcd.get("f_pi_GeV", "0.09207")))
        m_u = mp.mpf(str(qcd.get("m_u_GeV", "0.00216")))
        m_d = mp.mpf(str(qcd.get("m_d_GeV", "0.00467")))

        varphi_n = e8_phi_n(n=PQ_BLOCK_N, c=c)
        zeta_pq = _zeta_block(c3=mp.mpf(c.c3), r_B=PQ_BLOCK_RANK, I1=PQ_BLOCK_I1)
        fa_derived = zeta_pq * MPL_UNREDUCED_GEV * varphi_n
        ma_derived = _axion_mass_micro_eV(f_a_GeV=fa_derived, m_pi_GeV=m_pi, f_pi_GeV=f_pi, m_u_GeV=m_u, m_d_GeV=m_d)

        checks: list[Check] = []
        checks.append(mk_check_pass("pq_block_stage_n", f"n={PQ_BLOCK_N} (claim_n={claim_n})") if claim_n == PQ_BLOCK_N else mk_check_warn("pq_block_stage_n", f"n={PQ_BLOCK_N} (claim_n={claim_n})"))

        if mp.isfinite(claim_fa) and claim_fa > 0:
            rel = abs(fa_derived - claim_fa) / claim_fa
            checks.append(
                mk_check_pass("f_a_derived_not_quoted", f"f_a={fa_derived} GeV (rel={rel})")
                if rel <= FA_REL_TOL
                else mk_check_warn("f_a_derived_not_quoted", f"f_a={fa_derived} GeV (rel={rel})")
            )
        else:
            checks.append(mk_check_warn("f_a_derived_not_quoted", "quoted f_a missing; cannot compare"))

        if mp.isfinite(claim_ma) and claim_ma > 0:
            rel = abs(ma_derived - claim_ma) / claim_ma
            checks.append(
                mk_check_pass("m_a_follows_from_f_a", f"m_a={ma_derived} µeV (rel={rel})")
                if rel <= MA_REL_TOL
                else mk_check_warn("m_a_follows_from_f_a", f"m_a={ma_derived} µeV (rel={rel})")
            )
        else:
            checks.append(mk_check_warn("m_a_follows_from_f_a", "quoted m_a missing; cannot compare"))

        lines: list[str] = []
        lines += [
            "Axion PQ block derivation (n=10)",
            "",
            "Inputs:",
            f"- axion claim file: {ax_path}",
            f"- TFPT invariants: c3={c.c3}, varphi0={c.varphi0}, gamma0={c.gamma0}, lambda={c.e8_lambda}",
            f"- PQ block: n={PQ_BLOCK_N}, r_B={PQ_BLOCK_RANK}, I1={PQ_BLOCK_I1}",
            "",
            "E8 ladder:",
            f"- varphi_{PQ_BLOCK_N} = {varphi_n}",
            "",
            "Block factor:",
            f"- zeta_PQ = {zeta_pq}",
            "",
            "Derived scales:",
            f"- f_a = {fa_derived} GeV (M_Pl={MPL_UNREDUCED_GEV} GeV)",
            f"- m_a = {ma_derived} µeV (QCD inputs: m_pi={m_pi}, f_pi={f_pi}, m_u={m_u}, m_d={m_d})",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
        ]

        return ModuleResult(
            results={
                "inputs": {
                    "axion_claim_file": str(ax_path),
                    "block_stage_n": PQ_BLOCK_N,
                    "block_rank_rB": float(PQ_BLOCK_RANK),
                    "block_I1": float(PQ_BLOCK_I1),
                    "Mpl_unreduced_GeV": float(MPL_UNREDUCED_GEV),
                },
                "e8_ladder": {
                    "varphi_n": varphi_n,
                    "gamma0": c.gamma0,
                    "lambda": c.e8_lambda,
                    "definition": "varphi_n = varphi0 * exp(-gamma0) * (D_n/D_1)^lambda for n>=1; D_n=60-2n, D_1=58",
                },
                "block_factor": {
                    "zeta_PQ": zeta_pq,
                    "definition": "zeta_B = (pi c3) * exp(-beta_B pi c3) * exp(-k_B/c3), beta_B=(8-r_B)/8, k_B=(3/2)I1",
                },
                "derived": {
                    "f_a_GeV": fa_derived,
                    "m_a_micro_eV": ma_derived,
                },
                "claim": {
                    "block_stage_n": claim_n,
                    "f_a_GeV": claim_fa,
                    "m_a_micro_eV": claim_ma,
                },
                "qcd_inputs": {
                    "m_pi_GeV": m_pi,
                    "f_pi_GeV": f_pi,
                    "m_u_GeV": m_u,
                    "m_d_GeV": m_d,
                },
                "publication_grade_gap": False,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )
