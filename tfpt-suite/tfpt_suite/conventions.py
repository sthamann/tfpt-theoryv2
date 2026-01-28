from __future__ import annotations

"""
Central conventions and conversions used across the TFPT suite.

This file exists to make "scheme + normalization" choices explicit and reusable.

Gauge / hypercharge conventions:
- Hypercharge definition is the SM one: Q = T3 + Y (see `data/microscopic_action_tfpt_v25.json`).
- Some RG sources use SM-normalized hypercharge coupling g' (= g_Y).
- Others use the GUT-normalized g1 := sqrt(5/3) * g_Y.
"""

import math
from dataclasses import dataclass


_G1_GUT_OVER_GY = math.sqrt(5.0 / 3.0)


@dataclass(frozen=True)
class GaugeConvention:
    """
    Minimal gauge-coupling convention descriptor.

    `g1_is_gut_normalized` tells you what a variable named `g1` means:
    - True  -> g1 = sqrt(5/3) g_Y (GUT normalization)
    - False -> g1 = g_Y = g'       (SM normalization)
    """

    g1_is_gut_normalized: bool
    hypercharge_definition: str = "Q = T3 + Y"

    @staticmethod
    def sm_hypercharge() -> "GaugeConvention":
        return GaugeConvention(g1_is_gut_normalized=False)

    @staticmethod
    def gut_hypercharge() -> "GaugeConvention":
        return GaugeConvention(g1_is_gut_normalized=True)


def g1_gut_over_gY() -> float:
    """
    Return the fixed normalization factor g1_GUT / gY = sqrt(5/3).
    """

    return float(_G1_GUT_OVER_GY)


def g1_gut_from_gY(gY: float) -> float:
    return float(g1_gut_over_gY() * float(gY))


def gY_from_g1_gut(g1_gut: float) -> float:
    return float(float(g1_gut) / g1_gut_over_gY())


def alpha_from_g(g: float) -> float:
    """
    α = g^2 / (4π)
    """

    return float((float(g) ** 2) / (4.0 * math.pi))


def g_from_alpha(alpha: float) -> float:
    """
    g = sqrt(4π α)
    """

    a = float(alpha)
    if a < 0:
        raise ValueError("alpha must be non-negative")
    return float(math.sqrt(4.0 * math.pi * a))


def convert_g1(*, g1: float, from_conv: GaugeConvention, to_conv: GaugeConvention) -> float:
    """
    Convert a hypercharge coupling between SM and GUT normalizations.
    """

    if from_conv.hypercharge_definition != to_conv.hypercharge_definition:
        raise ValueError(
            f"Hypercharge definition mismatch: {from_conv.hypercharge_definition!r} vs {to_conv.hypercharge_definition!r}"
        )
    if from_conv.g1_is_gut_normalized == to_conv.g1_is_gut_normalized:
        return float(g1)
    if from_conv.g1_is_gut_normalized:
        return gY_from_g1_gut(g1)
    return g1_gut_from_gY(g1)

