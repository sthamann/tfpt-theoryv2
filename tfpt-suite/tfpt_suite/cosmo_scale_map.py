from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class CosmoScaleInputs:
    """
    Minimal cosmology mapping inputs used by `k_calibration`.

    This is *not* a full reheating solver; it is a deterministic mapping layer that turns
    explicitly stated assumptions into a scale-factor ratio a0/a_transition.
    """

    transition: str  # semantic label only (bounce_end | inflation_start | horizon_exit_of_pivot)
    N_inflation_from_transition: float
    # optional extra expansion between inflation end and reheating (a_reh/a_end = exp(N_reheat))
    N_reheat: float = 0.0
    # instantaneous reheating temperature (used for entropy-conserving a0/a_reh)
    T_reheat_GeV: float = 0.0
    g_star_s_reheat: float = 106.75
    g_star_s_today: float = 3.91
    T0_K: float = 2.7255


def gev_to_mpc_inv() -> float:
    """
    1 GeV in natural units corresponds to 1/(ħc) in inverse length.
      ħc = 0.1973269804 GeV·fm  ⇒  1 GeV = 5.0677307 fm^{-1}
      1 fm^{-1} = 1e15 m^{-1}
      1 Mpc = 3.085677581491367e22 m
    """
    hbarc_GeV_fm = 0.1973269804
    inv_fm_per_GeV = 1.0 / hbarc_GeV_fm
    inv_m_per_GeV = inv_fm_per_GeV * 1e15
    m_per_Mpc = 3.085677581491367e22
    return float(inv_m_per_GeV * m_per_Mpc)


def a0_over_a_reheat_from_entropy(
    *,
    T_reheat_GeV: float,
    g_star_s_reheat: float,
    g_star_s_today: float,
    T0_K: float,
) -> float:
    """
    Entropy-conserving scaling after reheating (instantaneous):
      a0/a_reh = (T_reh/T0) * (g*_s,reheat / g*_s,0)^{1/3}
    """
    kB_eV_per_K = 8.617333262145e-5
    T0_eV = float(T0_K) * kB_eV_per_K
    T0_GeV = T0_eV * 1e-9
    if T_reheat_GeV <= 0 or T0_GeV <= 0:
        return float("nan")
    if g_star_s_today <= 0:
        return float("nan")
    g_ratio = (float(g_star_s_reheat) / float(g_star_s_today)) ** (1.0 / 3.0)
    return float((float(T_reheat_GeV) / float(T0_GeV)) * g_ratio)


def a0_over_a_transition(inp: CosmoScaleInputs) -> tuple[float | None, str]:
    """
    Compute a0/a_transition under explicit assumptions.

    Returns:
      (value_or_None, note)
    """
    # Inflation expansion from transition to end:
    a_end_over_a_tr = math.exp(float(inp.N_inflation_from_transition)) if float(inp.N_inflation_from_transition) != 0.0 else 1.0

    # Reheating expansion (optional): a_reh/a_end
    a_reh_over_a_end = math.exp(float(inp.N_reheat)) if float(inp.N_reheat) != 0.0 else 1.0

    if float(inp.T_reheat_GeV) <= 0.0:
        return None, "T_reheat_GeV not provided; cannot fix absolute a0/a_transition (needs reheating+entropy model)"

    a0_over_areh = a0_over_a_reheat_from_entropy(
        T_reheat_GeV=float(inp.T_reheat_GeV),
        g_star_s_reheat=float(inp.g_star_s_reheat),
        g_star_s_today=float(inp.g_star_s_today),
        T0_K=float(inp.T0_K),
    )
    if not math.isfinite(a0_over_areh) or a0_over_areh <= 0:
        return None, "invalid a0/a_reh from entropy mapping"

    a0_over_a_tr = float(a0_over_areh * a_reh_over_a_end * a_end_over_a_tr)
    return a0_over_a_tr, f"a0/a_transition = (a0/a_reh)*(a_reh/a_end)*exp(N_inflation_from_transition); transition={inp.transition}"


def k_mpc_inv_from_k_hat(*, k_hat: float, M_GeV: float, a0_over_a_tr: float) -> float:
    """
    Map dimensionless k_hat ≡ k/M (used in the bounce solver) to a comoving k in Mpc^{-1}
    under a specified scale-factor normalization a0/a_transition.
    """
    if a0_over_a_tr <= 0:
        return float("nan")
    return float(float(k_hat) * float(M_GeV) * gev_to_mpc_inv() / float(a0_over_a_tr))


def ell_from_k_hat(*, k_hat: float, M_GeV: float, a0_over_a_tr: float, chi_star_Mpc: float) -> float:
    """
    Flat-sky proxy: ℓ ≈ k(Mpc^{-1}) · χ_*(Mpc).
    """
    k_mpc = k_mpc_inv_from_k_hat(k_hat=k_hat, M_GeV=M_GeV, a0_over_a_tr=a0_over_a_tr)
    return float(k_mpc * float(chi_star_Mpc))


# ---------------------------------------------------------------------------
# v1.06 reheating-window helpers (deterministic policy layer)
# ---------------------------------------------------------------------------


MPL_REDUCED_GEV = 2.435e18


def As_from_ln10_As(ln10_As: float) -> float:
    """
    Convert Planck-style ln(10^{10} A_s) to A_s.
    """
    return float(math.exp(float(ln10_As)) * 1e-10)


def starobinsky_N_from_ns(ns: float) -> float:
    """
    Starobinsky / α-attractor plateau relation used in v1.06:
      n_s ≈ 1 - 2/N  =>  N ≈ 2/(1-n_s)
    """
    ns_f = float(ns)
    if not (0.0 < ns_f < 1.0):
        return float("nan")
    denom = 1.0 - ns_f
    if denom <= 0:
        return float("nan")
    return float(2.0 / denom)


def starobinsky_r_from_N(N: float) -> float:
    """
    Starobinsky relation:
      r ≈ 12/N^2
    """
    Nf = float(N)
    if Nf <= 0:
        return float("nan")
    return float(12.0 / (Nf * Nf))


def V_star_from_As_r(As: float, r: float) -> float:
    """
    Slow-roll relation (reduced Planck units M̄_P=1):
      A_s = V / (24π^2 ε), ε = r/16  =>  V = (3/2) π^2 A_s r
    Returns V in units of M̄_P^4.
    """
    Asf = float(As)
    rf = float(r)
    if Asf <= 0 or rf <= 0:
        return float("nan")
    return float(1.5 * (math.pi**2) * Asf * rf)


def rho_reheat_GeV4(*, T_reheat_GeV: float, g_star: float) -> float:
    """
    Radiation energy density at reheating:
      ρ_reh = (π^2/30) g* T^4
    """
    T = float(T_reheat_GeV)
    if T <= 0:
        return float("nan")
    return float((math.pi**2 / 30.0) * float(g_star) * (T**4))


def rho_end_GeV4_from_As_r(*, As: float, r: float, c_end: float, Mpl_reduced_GeV: float = MPL_REDUCED_GEV) -> float:
    """
    v1.06-style proxy for end-of-inflation energy density:
      ρ_end ≈ c_end · V0, with V0 fixed by A_s and r at the pivot (plateau proxy).

    We implement:
      V0 ≈ V_*(A_s,r) in reduced Planck units,
      then ρ_end[GeV^4] = c_end · V0 · M̄_P^4.
    """
    V0 = V_star_from_As_r(As, r)
    if not math.isfinite(V0) or V0 <= 0:
        return float("nan")
    M = float(Mpl_reduced_GeV)
    return float(float(c_end) * V0 * (M**4))


def N_reheat_from_rho_ratio(*, w_reh: float, rho_reh: float, rho_end: float) -> float:
    """
    Reheating e-folds from ρ scaling:
      ρ ∝ a^{-3(1+w)}  =>  a_reh/a_end = (ρ_end/ρ_reh)^{1/(3(1+w))}
      N_reheat := ln(a_reh/a_end)
    """
    w = float(w_reh)
    if rho_reh <= 0 or rho_end <= 0:
        return float("nan")
    denom = 3.0 * (1.0 + w)
    if denom <= 0:
        return float("nan")
    return float((1.0 / denom) * math.log(float(rho_end) / float(rho_reh)))


def deltaN_from_rho_ratio(*, w_reh: float, rho_reh: float, rho_end: float) -> float:
    """
    v1.06-style ΔN reheating shift:
      ΔN ≈ (1-3w)/(12(1+w)) ln(ρ_reh/ρ_end)
    """
    w = float(w_reh)
    if rho_reh <= 0 or rho_end <= 0:
        return float("nan")
    denom = 12.0 * (1.0 + w)
    if denom == 0:
        return float("nan")
    return float(((1.0 - 3.0 * w) / denom) * math.log(float(rho_reh) / float(rho_end)))

