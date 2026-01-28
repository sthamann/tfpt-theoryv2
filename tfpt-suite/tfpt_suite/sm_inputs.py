from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tfpt_suite.conventions import g1_gut_from_gY


@dataclass(frozen=True)
class SmMzInputs:
    """
    Minimal SM MSbar-ish inputs at μ = M_Z to initialize gauge couplings.

    Conventions:
    - Use GUT-normalized g1 ≡ sqrt(5/3) g_Y (so g_Y is the SM hypercharge coupling g′).
    """

    mu_GeV: float
    alpha_em_inv: float
    sin2_thetaW: float
    alpha_s: float


def gauge_couplings_from_mz_inputs(inp: SmMzInputs) -> tuple[float, float, float]:
    """
    Compute (g1_GUT, g2, g3) at μ=MZ from (α_em(MZ), sin^2 θW(MZ), α_s(MZ)).

    Notes:
    - This is an algebraic conversion, not RG running.
    - `g1_GUT` is the GUT-normalized hypercharge coupling.
    """
    alpha_em = 1.0 / float(inp.alpha_em_inv)
    e = float(np.sqrt(4.0 * np.pi * alpha_em))

    s2 = float(inp.sin2_thetaW)
    s = float(np.sqrt(s2))
    c = float(np.sqrt(max(0.0, 1.0 - s2)))
    if s == 0 or c == 0:
        raise ValueError("Invalid sin^2(thetaW) input")

    g2 = e / s
    gY = e / c
    g1_gut = g1_gut_from_gY(gY)

    g3 = float(np.sqrt(4.0 * np.pi * float(inp.alpha_s)))
    return float(g1_gut), float(g2), float(g3)

