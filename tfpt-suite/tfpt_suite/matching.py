from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Literal, Optional

import numpy as np

from tfpt_suite.conventions import alpha_from_g, g_from_alpha


Scheme = Literal["MSbar"]
Direction = Literal["up", "down"]

PI = math.pi
SQRT2 = math.sqrt(2.0)
ZETA2 = math.pi**2 / 6.0
ZETA3 = 1.202056903159594
ZETA5 = 1.0369277551433698
B0_GAUSS_LEGENDRE_N = 160
LOG_ARG_FLOOR = 1e-300
DEFAULT_NC = 3


@lru_cache(maxsize=8)
def _gauss_legendre_nodes_weights(n_points: int) -> tuple[np.ndarray, np.ndarray]:
    if n_points < 2:
        raise ValueError("n_points must be >= 2 for Gauss-Legendre integration")
    nodes, weights = np.polynomial.legendre.leggauss(int(n_points))
    x = 0.5 * (nodes + 1.0)
    w = 0.5 * weights
    return x, w


def a0_msbar(*, mass_GeV: float, mu_GeV: float) -> float:
    """
    Buttazzo et al. (2013), Appendix A: A0(M) = M^2 [1 - ln(M^2 / mu^2)].
    """
    m = float(mass_GeV)
    mu = float(mu_GeV)
    if mu <= 0:
        raise ValueError("mu_GeV must be positive")
    if m < 0:
        raise ValueError("mass_GeV must be non-negative")
    if m == 0.0:
        return 0.0
    m2 = m * m
    mu2 = mu * mu
    return float(m2 * (1.0 - math.log(m2 / mu2)))


def b0_msbar(
    *,
    p_GeV: float,
    m1_GeV: float,
    m2_GeV: float,
    mu_GeV: float,
    n_points: int = B0_GAUSS_LEGENDRE_N,
) -> float:
    """
    Buttazzo et al. (2013), Appendix A: B0(p; M1, M2) = -∫_0^1 ln[(x M1^2 + (1-x) M2^2 - x(1-x) p^2)/mu^2] dx.
    Returns the real part (log |arg|) for above-threshold kinematics.
    """
    p = float(p_GeV)
    m1 = float(m1_GeV)
    m2 = float(m2_GeV)
    mu = float(mu_GeV)
    if mu <= 0:
        raise ValueError("mu_GeV must be positive")
    if p < 0 or m1 < 0 or m2 < 0:
        raise ValueError("p_GeV, m1_GeV, m2_GeV must be non-negative")
    p2 = p * p
    m1_2 = m1 * m1
    m2_2 = m2 * m2
    mu2 = mu * mu
    x, w = _gauss_legendre_nodes_weights(int(n_points))
    arg = (x * m1_2) + ((1.0 - x) * m2_2) - (x * (1.0 - x) * p2)
    arg_abs = np.abs(arg)
    arg_floor = float(LOG_ARG_FLOOR) * mu2
    if arg_floor <= 0.0:
        raise ValueError("LOG_ARG_FLOOR must be positive")
    arg_abs = np.maximum(arg_abs, arg_floor)
    integrand = np.log(arg_abs / mu2)
    value = -float(np.sum(w * integrand))
    if not math.isfinite(value):
        raise RuntimeError("B0 integral returned a non-finite value")
    return value


def delta_y_t_ew_hempfling(
    *,
    m_t_GeV: float,
    m_h_GeV: float,
    mu_GeV: float,
    v_ev_GeV: float,
    n_c: int = DEFAULT_NC,
) -> float:
    """
    Hempfling & Kniehl (hep-ph/9408313), Eq. (2.6):
    δ_t^w(μ) = (G_F m_t^2 / (8 π^2 √2)) [
        -(N_c + 3/2) ln(m_t^2/μ^2) + N_c/2 + 4 - r + 2 r (2 r - 3) ln(4 r)
        - 8 r^2 (1 - 1/r)^{3/2} arcosh(√r) ]  for r ≥ 1,
    with r = M_H^2 / (4 m_t^2).
    For r < 1, replace (1 - 1/r)^{3/2} arcosh(√r) by (1/r - 1)^{3/2} arccos(√r).
    """
    mt = float(m_t_GeV)
    mh = float(m_h_GeV)
    mu = float(mu_GeV)
    v = float(v_ev_GeV)
    nc = int(n_c)
    if mt <= 0 or mh <= 0 or mu <= 0 or v <= 0:
        raise ValueError("m_t_GeV, m_h_GeV, mu_GeV, v_ev_GeV must be positive")
    if nc <= 0:
        raise ValueError("n_c must be positive")
    gf = 1.0 / (SQRT2 * v * v)
    r = (mh * mh) / (4.0 * mt * mt)
    if r <= 0:
        raise ValueError("r must be positive")
    log_mt = math.log((mt * mt) / (mu * mu))
    if r >= 1.0:
        sqrt_r = math.sqrt(r)
        factor = math.acosh(sqrt_r)
        pow_term = (1.0 - 1.0 / r) ** 1.5
    else:
        sqrt_r = math.sqrt(r)
        factor = math.acos(sqrt_r)
        pow_term = (1.0 / r - 1.0) ** 1.5
    bracket = (
        -(nc + 1.5) * log_mt
        + (nc / 2.0)
        + 4.0
        - r
        + 2.0 * r * (2.0 * r - 3.0) * math.log(4.0 * r)
        - 8.0 * r * r * pow_term * factor
    )
    pref = (gf * mt * mt) / (8.0 * PI * PI * SQRT2)
    return float(pref * bracket)


def delta_lambda_h_1loop_buttazzo(
    *,
    m_h_GeV: float,
    m_t_GeV: float,
    m_w_GeV: float,
    m_z_GeV: float,
    mu_GeV: float,
    v_ev_GeV: float,
) -> float:
    """
    Buttazzo et al. (arXiv:1307.3536), Appendix A.1:
    λ(1)(μ) = (1/(4π)^2 V^4) Re[ ... ],
    with A0 and B0 defined by Eq. (89).
    """
    mh = float(m_h_GeV)
    mt = float(m_t_GeV)
    mw = float(m_w_GeV)
    mz = float(m_z_GeV)
    mu = float(mu_GeV)
    v = float(v_ev_GeV)
    if mh <= 0 or mt <= 0 or mw <= 0 or mz <= 0 or mu <= 0 or v <= 0:
        raise ValueError("m_h_GeV, m_t_GeV, m_w_GeV, m_z_GeV, mu_GeV, v_ev_GeV must be positive")
    mh2 = mh * mh
    mt2 = mt * mt
    mw2 = mw * mw
    mz2 = mz * mz
    mh4 = mh2 * mh2
    mw4 = mw2 * mw2
    mz4 = mz2 * mz2
    v4 = v * v * v * v
    pref = 1.0 / ((4.0 * PI) ** 2 * v4)

    term = 0.0
    term += 3.0 * mt2 * (mh2 - 4.0 * mt2) * b0_msbar(p_GeV=mh, m1_GeV=mt, m2_GeV=mt, mu_GeV=mu)
    term += 3.0 * mh2 * a0_msbar(mass_GeV=mt, mu_GeV=mu)
    term += 0.25 * (mh4 - 4.0 * mh2 * mz2 + 12.0 * mz4) * b0_msbar(p_GeV=mh, m1_GeV=mz, m2_GeV=mz, mu_GeV=mu)
    term += (mh2 * (7.0 * mw2 - 4.0 * mz2) / (2.0 * (mz2 - mw2))) * a0_msbar(mass_GeV=mz, mu_GeV=mu)
    term += 0.5 * (mh4 - 4.0 * mh2 * mw2 + 12.0 * mw4) * b0_msbar(p_GeV=mh, m1_GeV=mw, m2_GeV=mw, mu_GeV=mu)
    term += -(3.0 * mh2 * mw2 / (2.0 * (mh2 - mw2))) * a0_msbar(mass_GeV=mh, mu_GeV=mu)
    term += (mh2 / 2.0) * (
        -11.0 + (3.0 * mh2 / (mh2 - mw2)) - (3.0 * mw2 / (mz2 - mw2))
    ) * a0_msbar(mass_GeV=mw, mu_GeV=mu)
    term += 2.25 * mh4 * b0_msbar(p_GeV=mh, m1_GeV=mh, m2_GeV=mh, mu_GeV=mu)
    term += 0.25 * (mh4 + mh2 * (mz2 + 2.0 * mw2 - 6.0 * mt2) - 8.0 * (mz4 + 2.0 * mw4))

    return float(pref * term)


def delta_alpha_msbar_on_shell_shift_pdg(*, alpha0: float, alpha_s_MZ: float, MZ_GeV: float, MW_GeV: float) -> float:
    """
    PDG (2024) Eq. 10.13: Δα̂(MZ) − Δα(MZ) (finite MSbar↔on-shell shift).

    Δα̂(MZ) − Δα(MZ) = (α/π) [100/27 - 1/6 - (7/4) ln(MZ^2/MW^2)
        + (α_s/π)(605/108 - (44/9) ζ(3))
        + (α_s/π)^2 (976481/23328 - (253/36) ζ(2) - (781/18) ζ(3) + (275/27) ζ(5))
        + (α_s/π)^3 (-349280239/559872 + (1351/14) ζ(2) + (18229/72) ζ(3)
          + (331/2) ζ(2) ζ(3) - (6125/16) ζ(5)) ].
    """
    a0 = float(alpha0)
    a_s = float(alpha_s_MZ)
    mz = float(MZ_GeV)
    mw = float(MW_GeV)
    if a0 <= 0:
        raise ValueError("alpha0 must be positive")
    if a_s < 0:
        raise ValueError("alpha_s_MZ must be non-negative")
    if mz <= 0 or mw <= 0:
        raise ValueError("MZ_GeV and MW_GeV must be positive")

    pi = math.pi
    log_term = math.log((mz * mz) / (mw * mw))
    a_s_pi = a_s / pi
    term0 = (100.0 / 27.0) - (1.0 / 6.0) - (7.0 / 4.0) * log_term
    term1 = a_s_pi * ((605.0 / 108.0) - (44.0 / 9.0) * ZETA3)
    term2 = (a_s_pi**2) * (
        (976481.0 / 23328.0) - (253.0 / 36.0) * ZETA2 - (781.0 / 18.0) * ZETA3 + (275.0 / 27.0) * ZETA5
    )
    term3 = (a_s_pi**3) * (
        (-349280239.0 / 559872.0)
        + (1351.0 / 14.0) * ZETA2
        + (18229.0 / 72.0) * ZETA3
        + (331.0 / 2.0) * ZETA2 * ZETA3
        - (6125.0 / 16.0) * ZETA5
    )
    return float((a0 / pi) * (term0 + term1 + term2 + term3))


@dataclass(frozen=True)
class MatchingOutcome:
    """
    Generic matching outcome container.
    """

    status: str
    note: str
    matching_active: bool
    deltas: dict[str, float]
    details: dict[str, object]


def match_alpha3_quark_threshold_2loop(*, alpha3: float, direction: Direction) -> float:
    """
    Two-loop MSbar decoupling constant for α_s at a heavy-quark threshold matched at μ = m_Q.

    Used in `rge_sm.run_sm_gauge_only_2loop_thresholds` and exposed here as a reusable primitive.

    Conventions:
    - direction="up": nf → nf+1 (integrating *in* a heavy quark when running upward)
    - direction="down": nf+1 → nf (integrating *out* a heavy quark when running downward)
    """
    a = float(alpha3)
    if a < 0:
        raise ValueError("alpha3 must be non-negative")
    c2 = 11.0 / 72.0
    if a == 0.0:
        return 0.0

    # Forward (down) map:
    #   α_nf = α_{nf+1} [1 + c2 (α_{nf+1}/π)^2]
    # The upward direction must invert this relation (not merely flip the sign),
    # otherwise down(up(α)) differs from α by O(α^5) at typical SM α_s values.
    if direction == "up":
        # Solve for x >= 0:
        #   a = x * (1 + c2 (x/pi)^2)
        # i.e. (c2/pi^2) x^3 + x - a = 0
        pi = math.pi
        k = c2 / (pi * pi)

        # Newton iteration with a conservative starting guess.
        x = a
        for _ in range(30):
            fx = (k * x * x * x) + x - a
            dfx = (3.0 * k * x * x) + 1.0
            if dfx == 0.0:
                break
            step = fx / dfx
            x_new = x - step
            # Keep positivity; the physical solution is positive for a>0.
            if x_new <= 0.0:
                x_new = 0.5 * x
            if abs(step) <= 1e-16 * max(1.0, abs(x_new)):
                x = x_new
                break
            x = x_new

        # Fallback: solve cubic via roots if Newton produced a non-finite result.
        if not math.isfinite(x) or x <= 0.0:
            roots = np.roots([k, 0.0, 1.0, -a])
            real_pos = [float(r.real) for r in roots if abs(float(r.imag)) < 1e-10 and float(r.real) > 0.0]
            if not real_pos:
                raise RuntimeError("Failed to invert 2-loop α_s threshold matching (no positive real root)")
            x = min(real_pos)
        return float(x)
    if direction == "down":
        return float(a * (1.0 + c2 * (a / math.pi) ** 2))
    raise ValueError(f"Invalid direction: {direction!r}")


def match_gauge(
    *,
    threshold_id: str,
    mu_thr_GeV: float,
    direction: Direction,
    couplings_below: dict[str, float],
    scheme: Scheme = "MSbar",
    loop_order: int = 1,
    active_fields_before: Optional[list[str]] = None,
    active_fields_after: Optional[list[str]] = None,
    finite_delta_alpha: Optional[dict[str, float]] = None,
) -> tuple[dict[str, float], MatchingOutcome]:
    """
    Gauge-coupling matching across a threshold μ=mu_thr.

    Phase-1 (minimal) implementation:
    - Deterministic, explicit, and conservative.
    - By default returns the identity map (continuous matching at μ=mu_thr).
    - Optional `finite_delta_alpha` can encode *known* finite pieces in terms of δα_i,
      keyed by "alphaY", "alpha2", "alpha3".

    NOTE:
    One-loop MSbar decoupling for gauge couplings is log-only; if you match at μ=M_heavy,
    the one-loop step is zero. Finite steps start at higher loop order.
    """
    if scheme != "MSbar":
        raise ValueError(f"Unsupported scheme: {scheme!r}")
    if not str(threshold_id).strip():
        raise ValueError("threshold_id must be a non-empty string")
    if float(mu_thr_GeV) <= 0:
        raise ValueError("mu_thr_GeV must be positive")
    if loop_order < 1:
        raise ValueError("loop_order must be >= 1")
    if direction not in ("up", "down"):
        raise ValueError(f"Invalid direction: {direction!r}")

    out = {k: float(v) for k, v in dict(couplings_below).items()}
    deltas = dict(finite_delta_alpha or {})
    fields_before = [str(x) for x in (active_fields_before or [])]
    fields_after = [str(x) for x in (active_fields_after or [])]
    delta_out: dict[str, float] = {}
    delta_alpha_applied: dict[str, float] = {}

    # Apply finite shifts in α variables (if provided).
    # This keeps the transformation explicit and numerically safe (α>=0).
    for key_g, key_a in (("gY", "alphaY"), ("g2", "alpha2"), ("g3", "alpha3")):
        if key_g not in out:
            continue
        if key_a not in deltas:
            continue
        a0 = float(alpha_from_g(out[key_g]))
        # Convention: finite_delta_alpha encodes (alpha_above - alpha_below) for direction="up".
        # For direction="down" we invert the step.
        da = float(deltas[key_a])
        if direction == "down":
            da = -da
        a1 = float(a0 + da)
        if a1 < 0:
            raise ValueError(f"Finite matching would make {key_a}<0 at μ={mu_thr_GeV}: {a1}")
        g0 = float(out[key_g])
        g1 = float(g_from_alpha(a1))
        out[key_g] = g1
        delta_alpha_applied[key_a] = da
        delta_out[f"delta_{key_g}"] = float(g1 - g0)
        delta_out[f"delta_{key_a}"] = float(da)

    has_finite = any(abs(float(v)) > 0 for v in delta_alpha_applied.values())
    if has_finite:
        status = "matched_with_finite_pieces"
        note = "finite α-shifts applied (see deltas/details); sign convention is (above-below) for direction='up'"
    else:
        status = "matched_1loop_log_only_identity" if loop_order == 1 else f"matched_{loop_order}loop_identity"
        note = (
            "identity matching at μ=threshold (1-loop decoupling is log-only, so the step vanishes at μ=M); "
            "finite pieces require explicit inputs or higher-loop decoupling constants"
        )

    details: dict[str, Any] = {
        "threshold_id": str(threshold_id),
        "mu_thr_GeV": float(mu_thr_GeV),
        "direction": str(direction),
        "scheme": str(scheme),
        "loop_order": int(loop_order),
        "active_fields_before": fields_before,
        "active_fields_after": fields_after,
        "finite_delta_alpha_input": {str(k): float(v) for k, v in deltas.items()},
        "finite_delta_alpha_applied": {str(k): float(v) for k, v in delta_alpha_applied.items()},
    }
    return out, MatchingOutcome(
        status=status,
        note=note,
        matching_active=True,
        deltas=delta_out,
        details=details,
    )


def match_yukawa(
    *,
    threshold_id: str,
    mu_thr_GeV: float,
    direction: Direction,
    yukawas_below: dict[str, object],
    scheme: Scheme = "MSbar",
    loop_order: int = 1,
    active_fields_before: Optional[list[str]] = None,
    active_fields_after: Optional[list[str]] = None,
    finite_delta_yukawa: Optional[dict[str, float]] = None,
) -> tuple[dict[str, object], MatchingOutcome]:
    """
    Yukawa matching across a threshold.

    Phase-1 placeholder:
    - Deterministic identity mapping (continuous matching).
    - The top-threshold pole→MSbar conversion is handled in `rge_sm.sm_boundary_conditions_at_mt`.
    - Optional `finite_delta_yukawa` encodes relative shifts (y_above = y_below * (1 + δ) for direction="up").
    """
    if scheme != "MSbar":
        raise ValueError(f"Unsupported scheme: {scheme!r}")
    if not str(threshold_id).strip():
        raise ValueError("threshold_id must be a non-empty string")
    if float(mu_thr_GeV) <= 0:
        raise ValueError("mu_thr_GeV must be positive")
    if loop_order < 1:
        raise ValueError("loop_order must be >= 1")
    if direction not in ("up", "down"):
        raise ValueError(f"Invalid direction: {direction!r}")

    fields_before = [str(x) for x in (active_fields_before or [])]
    fields_after = [str(x) for x in (active_fields_after or [])]

    out: dict[str, object] = {}
    max_abs_in: dict[str, float] = {}
    for k, v in dict(yukawas_below).items():
        key = str(k)
        if isinstance(v, (int, float)):
            out[key] = float(v)
            max_abs_in[key] = abs(float(v))
        elif isinstance(v, np.ndarray):
            arr = np.array(v, dtype=complex, copy=True)
            out[key] = arr
            max_abs_in[key] = float(np.max(np.abs(arr))) if arr.size else 0.0
        else:
            raise TypeError(f"Unsupported Yukawa value type for {key}: {type(v)}")

    deltas = dict(finite_delta_yukawa or {})
    delta_out: dict[str, float] = {}
    delta_applied: dict[str, float] = {}
    for key, delta in deltas.items():
        if key not in out:
            continue
        d = float(delta)
        if not math.isfinite(d):
            raise ValueError(f"finite_delta_yukawa[{key}] must be finite")
        factor = 1.0 + d
        if factor == 0.0:
            raise ValueError(f"finite_delta_yukawa[{key}] would zero out the Yukawa (1+delta=0)")
        if direction == "down":
            factor = 1.0 / factor
        if isinstance(out[key], np.ndarray):
            out[key] = out[key] * factor
        else:
            out[key] = float(out[key]) * factor
        applied = float(factor - 1.0)
        delta_applied[key] = applied
        delta_out[f"delta_{key}_rel"] = applied

    has_finite = any(abs(float(v)) > 0 for v in delta_applied.values())
    if has_finite:
        status = "matched_with_finite_pieces"
        note = "finite Yukawa shifts applied (relative deltas; direction='up' uses y_above=y_below*(1+δ))"
    else:
        note = (
            "identity matching at μ=threshold (1-loop decoupling is log-only, so the step vanishes at μ=M); "
            "top pole→MSbar(mt) is handled in rge_sm.sm_boundary_conditions_at_mt"
        )
        status = "matched_1loop_log_only_identity" if loop_order == 1 else f"matched_{loop_order}loop_identity"
    details: dict[str, Any] = {
        "threshold_id": str(threshold_id),
        "mu_thr_GeV": float(mu_thr_GeV),
        "direction": str(direction),
        "scheme": str(scheme),
        "loop_order": int(loop_order),
        "active_fields_before": fields_before,
        "active_fields_after": fields_after,
        "yukawa_maxabs_input": max_abs_in,
        "finite_delta_yukawa_input": {str(k): float(v) for k, v in deltas.items()},
        "finite_delta_yukawa_applied": {str(k): float(v) for k, v in delta_applied.items()},
    }
    return out, MatchingOutcome(
        status=status,
        note=note,
        matching_active=True,
        deltas=delta_out,
        details=details,
    )


def match_quartic(
    *,
    threshold_id: str,
    mu_thr_GeV: float,
    direction: Direction,
    quartics_below: dict[str, float],
    scheme: Scheme = "MSbar",
    loop_order: int = 1,
    active_fields_before: Optional[list[str]] = None,
    active_fields_after: Optional[list[str]] = None,
    finite_delta_quartic: Optional[dict[str, float]] = None,
) -> tuple[dict[str, float], MatchingOutcome]:
    """
    Quartic matching across a threshold.

    Phase-1 placeholder:
    - Deterministic identity mapping (continuous matching).
    """
    if scheme != "MSbar":
        raise ValueError(f"Unsupported scheme: {scheme!r}")
    if not str(threshold_id).strip():
        raise ValueError("threshold_id must be a non-empty string")
    if float(mu_thr_GeV) <= 0:
        raise ValueError("mu_thr_GeV must be positive")
    if loop_order < 1:
        raise ValueError("loop_order must be >= 1")
    if direction not in ("up", "down"):
        raise ValueError(f"Invalid direction: {direction!r}")
    fields_before = [str(x) for x in (active_fields_before or [])]
    fields_after = [str(x) for x in (active_fields_after or [])]
    out = {k: float(v) for k, v in dict(quartics_below).items()}
    deltas = dict(finite_delta_quartic or {})
    delta_out: dict[str, float] = {}
    delta_applied: dict[str, float] = {}
    for key, delta in deltas.items():
        if key not in out:
            continue
        d = float(delta)
        if not math.isfinite(d):
            raise ValueError(f"finite_delta_quartic[{key}] must be finite")
        if direction == "down":
            d = -d
        out[key] = float(out[key] + d)
        delta_applied[key] = float(d)
        delta_out[f"delta_{key}"] = float(d)

    has_finite = any(abs(float(v)) > 0 for v in delta_applied.values())
    if has_finite:
        status = "matched_with_finite_pieces"
        note = "finite quartic shifts applied (additive deltas; direction='up' uses λ_above=λ_below+Δ)"
    else:
        status = "matched_1loop_log_only_identity" if loop_order == 1 else f"matched_{loop_order}loop_identity"
        note = "identity matching at μ=threshold (explicit); finite quartic matching not implemented yet"
    details: dict[str, Any] = {
        "threshold_id": str(threshold_id),
        "mu_thr_GeV": float(mu_thr_GeV),
        "direction": str(direction),
        "scheme": str(scheme),
        "loop_order": int(loop_order),
        "active_fields_before": fields_before,
        "active_fields_after": fields_after,
        "finite_delta_quartic_input": {str(k): float(v) for k, v in deltas.items()},
        "finite_delta_quartic_applied": {str(k): float(v) for k, v in delta_applied.items()},
    }
    return out, MatchingOutcome(status=status, note=note, matching_active=True, deltas=delta_out, details=details)

