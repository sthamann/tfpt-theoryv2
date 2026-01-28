from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from tfpt_suite.conventions import alpha_from_g, g1_gut_over_gY, gY_from_g1_gut
from tfpt_suite.rg_authority import enforce_no_legacy_sm_rge_above_mz
from tfpt_suite.sm_inputs import SmMzInputs, gauge_couplings_from_mz_inputs


def _beta_gauge_1loop(*, g1: float, g2: float, g3: float) -> tuple[float, float, float]:
    """
    1-loop SM beta functions for gauge couplings in GUT normalization.

    (16π^2) dg_i/dt = b_i g_i^3
      b1 = 41/10, b2 = -19/6, b3 = -7
    """
    b1 = 41.0 / 10.0
    b2 = -19.0 / 6.0
    b3 = -7.0
    fac = 1.0 / (16.0 * np.pi**2)
    return fac * b1 * (g1**3), fac * b2 * (g2**3), fac * b3 * (g3**3)


def _trace_yukawas(*, Yu: np.ndarray, Yd: np.ndarray, Ye: np.ndarray) -> float:
    """
    T = Tr(3 Yu^†Yu + 3 Yd^†Yd + Ye^†Ye)
    """
    Hu = Yu.conj().T @ Yu
    Hd = Yd.conj().T @ Yd
    He = Ye.conj().T @ Ye
    return float(3.0 * np.trace(Hu).real + 3.0 * np.trace(Hd).real + np.trace(He).real)


def _beta_yukawas_1loop(*, g1: float, g2: float, g3: float, Yu: np.ndarray, Yd: np.ndarray, Ye: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    1-loop SM matrix Yukawa RGEs (GUT-normalized g1).

    (16π^2) dY_u/dt = Y_u [ 3/2 (Y_u^†Y_u - Y_d^†Y_d) + T - (17/20 g1^2 + 9/4 g2^2 + 8 g3^2) ]
    (16π^2) dY_d/dt = Y_d [ 3/2 (Y_d^†Y_d - Y_u^†Y_u) + T - ( 1/4 g1^2 + 9/4 g2^2 + 8 g3^2) ]
    (16π^2) dY_e/dt = Y_e [ 3/2 (Y_e^†Y_e)              + T - ( 9/4 g1^2 + 9/4 g2^2) ]
    """
    Hu = Yu.conj().T @ Yu
    Hd = Yd.conj().T @ Yd
    He = Ye.conj().T @ Ye
    T = _trace_yukawas(Yu=Yu, Yd=Yd, Ye=Ye)

    I3 = np.eye(3, dtype=complex)
    gu = (17.0 / 20.0) * (g1**2) + (9.0 / 4.0) * (g2**2) + 8.0 * (g3**2)
    gd = (1.0 / 4.0) * (g1**2) + (9.0 / 4.0) * (g2**2) + 8.0 * (g3**2)
    ge = (9.0 / 4.0) * (g1**2) + (9.0 / 4.0) * (g2**2)

    fac = 1.0 / (16.0 * np.pi**2)
    dYu = fac * (Yu @ ((1.5 * (Hu - Hd)) + (T - gu) * I3))
    dYd = fac * (Yd @ ((1.5 * (Hd - Hu)) + (T - gd) * I3))
    dYe = fac * (Ye @ ((1.5 * (He)) + (T - ge) * I3))
    return dYu, dYd, dYe


def _pack_state(*, g1: float, g2: float, g3: float, Yu: np.ndarray, Yd: np.ndarray, Ye: np.ndarray) -> np.ndarray:
    """
    Pack (g1,g2,g3,Yu,Yd,Ye) into a real vector.
    """
    parts = [np.array([g1, g2, g3], dtype=float)]
    for M in (Yu, Yd, Ye):
        parts.append(M.real.reshape(-1))
        parts.append(M.imag.reshape(-1))
    return np.concatenate(parts, axis=0)


def _unpack_state(y: np.ndarray) -> tuple[float, float, float, np.ndarray, np.ndarray, np.ndarray]:
    """
    Unpack vector into (g1,g2,g3,Yu,Yd,Ye).
    """
    g1, g2, g3 = (float(y[0]), float(y[1]), float(y[2]))
    idx = 3

    def read_mat() -> np.ndarray:
        nonlocal idx
        re = y[idx : idx + 9].reshape(3, 3)
        idx += 9
        im = y[idx : idx + 9].reshape(3, 3)
        idx += 9
        return (re + 1j * im).astype(complex)

    Yu = read_mat()
    Yd = read_mat()
    Ye = read_mat()
    return g1, g2, g3, Yu, Yd, Ye


def run_sm_rge_1loop(
    *,
    mu_start_GeV: float,
    mu_end_GeV: float,
    g_start: tuple[float, float, float],
    Yu_start: np.ndarray,
    Yd_start: np.ndarray,
    Ye_start: np.ndarray,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    method: str = "DOP853",
) -> dict[str, object]:
    """
    Run 1-loop SM RGEs from mu_start to mu_end (both in GeV) in t = ln(mu) time.

    Returns end-state and diagnostics.
    """
    mu_start = float(mu_start_GeV)
    mu_end = float(mu_end_GeV)
    if mu_start <= 0 or mu_end <= 0:
        raise ValueError("mu_start_GeV and mu_end_GeV must be positive")
    enforce_no_legacy_sm_rge_above_mz(mu_start_GeV=mu_start, mu_end_GeV=mu_end, caller="rge_sm.run_sm_rge_1loop")

    t0 = float(np.log(mu_start))
    t1 = float(np.log(mu_end))

    g1_0, g2_0, g3_0 = g_start
    y0 = _pack_state(g1=g1_0, g2=g2_0, g3=g3_0, Yu=Yu_start, Yd=Yd_start, Ye=Ye_start)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        g1, g2, g3, Yu, Yd, Ye = _unpack_state(y)
        dg1, dg2, dg3 = _beta_gauge_1loop(g1=g1, g2=g2, g3=g3)
        dYu, dYd, dYe = _beta_yukawas_1loop(g1=g1, g2=g2, g3=g3, Yu=Yu, Yd=Yd, Ye=Ye)
        return _pack_state(g1=dg1, g2=dg2, g3=dg3, Yu=dYu, Yd=dYd, Ye=dYe)

    sol = solve_ivp(
        rhs,
        t_span=(t0, t1),
        y0=y0,
        method=method,
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(f"SM RGE integration failed: {sol.message}")

    g1_1, g2_1, g3_1, Yu_1, Yd_1, Ye_1 = _unpack_state(sol.y[:, -1])
    return {
        "mu_start_GeV": mu_start,
        "mu_end_GeV": mu_end,
        "g_start": {"g1": g1_0, "g2": g2_0, "g3": g3_0},
        "g_end": {"g1": g1_1, "g2": g2_1, "g3": g3_1},
        "Yu_end": Yu_1,
        "Yd_end": Yd_1,
        "Ye_end": Ye_1,
        "n_steps": int(sol.t.size),
    }


def run_sm_gauge_only_1loop(
    *,
    mu_start_GeV: float,
    mu_end_GeV: float,
    g_start: tuple[float, float, float],
    rtol: float = 1e-10,
    atol: float = 1e-12,
    method: str = "DOP853",
) -> tuple[float, float, float]:
    """
    Convenience: integrate only gauge couplings at 1-loop.
    """
    mu_start = float(mu_start_GeV)
    mu_end = float(mu_end_GeV)
    if mu_start <= 0 or mu_end <= 0:
        raise ValueError("mu_start_GeV and mu_end_GeV must be positive")
    enforce_no_legacy_sm_rge_above_mz(mu_start_GeV=mu_start, mu_end_GeV=mu_end, caller="rge_sm.run_sm_gauge_only_1loop")
    t0 = float(np.log(mu_start))
    t1 = float(np.log(mu_end))
    g1_0, g2_0, g3_0 = g_start
    y0 = np.array([g1_0, g2_0, g3_0], dtype=float)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        g1, g2, g3 = float(y[0]), float(y[1]), float(y[2])
        dg1, dg2, dg3 = _beta_gauge_1loop(g1=g1, g2=g2, g3=g3)
        return np.array([dg1, dg2, dg3], dtype=float)

    sol = solve_ivp(rhs, (t0, t1), y0=y0, method=method, rtol=rtol, atol=atol, dense_output=False)
    if not sol.success:
        raise RuntimeError(f"SM gauge-only integration failed: {sol.message}")
    return float(sol.y[0, -1]), float(sol.y[1, -1]), float(sol.y[2, -1])


def _qcd_nf(*, mu_GeV: float, mc_GeV: float, mb_GeV: float, mt_GeV: float) -> int:
    """
    Effective number of active quark flavors for QCD beta function coefficient.
    """
    if mu_GeV < mc_GeV:
        return 3
    if mu_GeV < mb_GeV:
        return 4
    if mu_GeV < mt_GeV:
        return 5
    return 6


def _b3_from_nf(nf: int) -> float:
    """
    1-loop QCD coefficient for dg3/dt in our convention:
      (16π^2) dg3/dt = b3 g3^3
    with b3 = -11 + 2/3 nf.
    """
    return float(-11.0 + (2.0 / 3.0) * float(nf))


def run_sm_gauge_only_1loop_thresholds(
    *,
    mu_start_GeV: float,
    mu_end_GeV: float,
    g_start: tuple[float, float, float],
    mc_GeV: float = 1.27,
    mb_GeV: float = 4.18,
    mt_GeV: float = 173.0,
    rtol: float = 1e-10,
    atol: float = 1e-12,
    method: str = "DOP853",
) -> tuple[float, float, float]:
    """
    1-loop gauge-only running with an nf-threshold model for b3 (QCD).

    Conventions:
    - g1 is GUT-normalized (sqrt(5/3) gY)
    - t = ln(mu)

    Notes:
    - At strict 1-loop, αs matching at quark thresholds is continuous. We therefore do not
      apply a finite matching step here (that starts at 2-loop).
    """
    mu0 = float(mu_start_GeV)
    mu1 = float(mu_end_GeV)
    if mu0 <= 0 or mu1 <= 0:
        raise ValueError("mu_start_GeV and mu_end_GeV must be positive")
    enforce_no_legacy_sm_rge_above_mz(
        mu_start_GeV=mu0, mu_end_GeV=mu1, caller="rge_sm.run_sm_gauge_only_1loop_thresholds"
    )

    # Constant SM b1,b2 (full SM field content).
    b1 = 41.0 / 10.0
    b2 = -19.0 / 6.0

    # Segment boundaries for QCD thresholds that lie between mu0 and mu1
    thresholds = [mc_GeV, mb_GeV, mt_GeV]
    lo = min(mu0, mu1)
    hi = max(mu0, mu1)
    cuts = [x for x in thresholds if lo < x < hi]
    cuts.sort()
    if mu1 < mu0:
        cuts = cuts[::-1]

    points = [mu0] + cuts + [mu1]
    y = np.array(list(g_start), dtype=float)

    fac = 1.0 / (16.0 * np.pi**2)

    def integrate_segment(mu_a: float, mu_b: float, y0: np.ndarray) -> np.ndarray:
        t0 = float(np.log(mu_a))
        t1 = float(np.log(mu_b))
        # Use a probe scale inside the segment to avoid boundary ambiguities at thresholds.
        mu_probe = float(np.sqrt(mu_a * mu_b))
        nf = _qcd_nf(mu_GeV=mu_probe, mc_GeV=mc_GeV, mb_GeV=mb_GeV, mt_GeV=mt_GeV)
        b3 = _b3_from_nf(nf)

        def rhs(t: float, yy: np.ndarray) -> np.ndarray:
            g1, g2, g3 = float(yy[0]), float(yy[1]), float(yy[2])
            dg1 = fac * b1 * (g1**3)
            dg2 = fac * b2 * (g2**3)
            dg3 = fac * b3 * (g3**3)
            return np.array([dg1, dg2, dg3], dtype=float)

        sol = solve_ivp(rhs, (t0, t1), y0=y0, method=method, rtol=rtol, atol=atol, dense_output=False)
        if not sol.success:
            raise RuntimeError(f"SM gauge-only 1-loop integration failed: {sol.message}")
        return sol.y[:, -1]

    for i in range(len(points) - 1):
        mu_a = float(points[i])
        mu_b = float(points[i + 1])
        y = integrate_segment(mu_a, mu_b, y)

    return float(y[0]), float(y[1]), float(y[2])


def _beta_gauge_2loop(*, g1: float, g2: float, g3: float, b1: float, b2: float, b3: float) -> tuple[float, float, float]:
    """
    Minimal 2-loop SM gauge RGEs (GUT-normalized g1), with an nf-dependent b3 in the 1-loop term.

    We use the standard SM 2-loop gauge matrix B (neglecting Yukawa-trace pieces here; those are included
    once the full 2-loop Yukawa+gauge system is implemented).
    """
    # 2-loop gauge matrix (SM, GUT-normalized g1)
    B = np.array([[199 / 50, 27 / 10, 44 / 5], [9 / 10, 35 / 6, 12], [11 / 10, 9 / 2, -26]], dtype=float)
    fac1 = 1.0 / (16.0 * np.pi**2)
    fac2 = fac1**2

    g = np.array([g1, g2, g3], dtype=float)
    b = np.array([b1, b2, b3], dtype=float)

    one = fac1 * b * (g**3)
    two = fac2 * (g**3) * (B @ (g**2))
    return float(one[0] + two[0]), float(one[1] + two[1]), float(one[2] + two[2])


def run_sm_gauge_only_2loop_thresholds(
    *,
    mu_start_GeV: float,
    mu_end_GeV: float,
    g_start: tuple[float, float, float],
    mc_GeV: float = 1.27,
    mb_GeV: float = 4.18,
    mt_GeV: float = 173.0,
    apply_alpha3_matching: bool = True,
    rtol: float = 1e-10,
    atol: float = 1e-12,
    method: str = "DOP853",
) -> tuple[float, float, float]:
    """
    Two-loop gauge-only running with a simple QCD flavor-threshold model for b3 and optional 2-loop α_s matching.

    This is intended as an intermediate step toward the full 2-loop + threshold-matching pipeline
    requested for the CKM/PMNS RG-dressed modules.

    Conventions:
    - g1 is GUT-normalized (sqrt(5/3) gY)
    - t = ln(mu)
    """
    mu0 = float(mu_start_GeV)
    mu1 = float(mu_end_GeV)
    if mu0 <= 0 or mu1 <= 0:
        raise ValueError("mu_start_GeV and mu_end_GeV must be positive")
    enforce_no_legacy_sm_rge_above_mz(
        mu_start_GeV=mu0, mu_end_GeV=mu1, caller="rge_sm.run_sm_gauge_only_2loop_thresholds"
    )

    # Constant SM b1,b2 (full SM field content); QCD b3 uses nf-effective coefficient below mt.
    b1 = 41.0 / 10.0
    b2 = -19.0 / 6.0

    # Segment boundaries for QCD thresholds that lie between mu0 and mu1
    thresholds = [mc_GeV, mb_GeV, mt_GeV]
    lo = min(mu0, mu1)
    hi = max(mu0, mu1)
    cuts = [x for x in thresholds if lo < x < hi]
    cuts.sort()
    if mu1 < mu0:
        cuts = cuts[::-1]

    # Build breakpoints
    points = [mu0] + cuts + [mu1]

    y = np.array(list(g_start), dtype=float)

    def integrate_segment(mu_a: float, mu_b: float, y0: np.ndarray) -> np.ndarray:
        t0 = float(np.log(mu_a))
        t1 = float(np.log(mu_b))
        # Use a probe scale inside the segment to avoid boundary ambiguities at thresholds.
        mu_probe = float(np.sqrt(mu_a * mu_b))
        nf = _qcd_nf(mu_GeV=mu_probe, mc_GeV=mc_GeV, mb_GeV=mb_GeV, mt_GeV=mt_GeV)
        b3 = _b3_from_nf(nf)

        def rhs(t: float, yy: np.ndarray) -> np.ndarray:
            g1, g2, g3 = float(yy[0]), float(yy[1]), float(yy[2])
            dg1, dg2, dg3 = _beta_gauge_2loop(g1=g1, g2=g2, g3=g3, b1=b1, b2=b2, b3=b3)
            return np.array([dg1, dg2, dg3], dtype=float)

        sol = solve_ivp(rhs, (t0, t1), y0=y0, method=method, rtol=rtol, atol=atol, dense_output=False)
        if not sol.success:
            raise RuntimeError(f"SM gauge-only 2-loop integration failed: {sol.message}")
        return sol.y[:, -1]

    def match_alpha3_up(alpha_below: float) -> float:
        # Two-loop decoupling at μ=m_Q: α_nf = α_nf+1 * [1 + (11/72)*(α/π)^2]
        c2 = 11.0 / 72.0
        a = float(alpha_below)
        return float(a * (1.0 - c2 * (a / np.pi) ** 2))

    def match_alpha3_down(alpha_above: float) -> float:
        c2 = 11.0 / 72.0
        a = float(alpha_above)
        return float(a * (1.0 + c2 * (a / np.pi) ** 2))

    for i in range(len(points) - 1):
        mu_a = float(points[i])
        mu_b = float(points[i + 1])
        y = integrate_segment(mu_a, mu_b, y)

        # Apply matching at heavy-quark thresholds (only affects α3 / g3)
        if apply_alpha3_matching and mu_b in thresholds:
            alpha3 = float((y[2] ** 2) / (4.0 * np.pi))
            if mu1 > mu0:
                alpha3_new = match_alpha3_up(alpha3)
            else:
                alpha3_new = match_alpha3_down(alpha3)
            y[2] = float(np.sqrt(4.0 * np.pi * alpha3_new))

    return float(y[0]), float(y[1]), float(y[2])


def pole_mass_to_msbar_mass_qcd_2loop(*, pole_mass_GeV: float, alpha_s_mu: float, n_light: int) -> float:
    """
    Convert a quark pole mass M to an MSbar running mass m(μ) at μ = M, using QCD up to 2 loops.

    We use the standard expansion (Chetyrkin/Steinhauser):

      M / m(M) = 1 + (4/3) a + (13.4434 - 1.0414 n_l) a^2 + O(a^3),  a = αs(M)/π

    where n_l is the number of *lighter* active flavors (excluding the heavy quark itself).

    For the top quark: n_l = 5.
    """
    M = float(pole_mass_GeV)
    a_s = float(alpha_s_mu)
    if M <= 0:
        raise ValueError("pole_mass_GeV must be positive")
    if not (a_s >= 0.0):
        raise ValueError("alpha_s_mu must be non-negative")
    nl = int(n_light)
    if nl < 0:
        raise ValueError("n_light must be >= 0")

    a = a_s / math.pi
    c1 = 4.0 / 3.0
    c2 = 13.4434 - 1.0414 * float(nl)
    return float(M / (1.0 + c1 * a + c2 * a * a))


def msbar_mass_run_qcd_1loop(*, m_mu0_GeV: float, alpha_s_mu0: float, alpha_s_mu1: float, nf: int) -> float:
    """
    1-loop QCD running of a MSbar quark mass:
      m(μ1) = m(μ0) * (αs(μ1)/αs(μ0))^(12/(33-2 nf))

    Notes:
    - This is the standard LO relation used as a deterministic primitive in the suite.
    - Upgrade path (publication-grade matching): replace by higher-loop running + matching.
    """
    m0 = float(m_mu0_GeV)
    if m0 <= 0:
        raise ValueError("m_mu0_GeV must be positive")
    a0 = float(alpha_s_mu0)
    a1 = float(alpha_s_mu1)
    if a0 <= 0 or a1 <= 0:
        raise ValueError("alpha_s inputs must be positive for QCD mass running")
    nf_i = int(nf)
    if nf_i < 3:
        raise ValueError("nf must be >=3")
    expo = 12.0 / (33.0 - 2.0 * float(nf_i))
    return float(m0 * (a1 / a0) ** expo)


@dataclass(frozen=True)
class SmBoundaryAtMt:
    """
    Boundary conditions for SM running at μ = m_t (mt→UV policy).
    """

    mu_mt_GeV: float
    scheme: str
    conventions: dict[str, str]
    route_2loop: dict[str, float]
    route_1loop: dict[str, float]
    diffs: dict[str, float]


def sm_boundary_conditions_at_mt(*, sm_inputs_mz: dict[str, Any]) -> SmBoundaryAtMt:
    """
    Construct a consistent set of MSbar boundary conditions at μ = m_t for mt→UV RG pipelines.

    Inputs (expected in `sm_inputs_mz`):
    - α_em_inv, sin2_thetaW, α_s at μ=MZ (interpreted as MSbar inputs at MZ)
    - pole/threshold masses: mu_GeV (=MZ), mt_GeV, mc_GeV, mb_GeV, mH_GeV (mW optional)

    Outputs:
    - (gY, g2, g3) at μ=mt (SM normalization for gY ≡ g′)
    - yt(mt) from mt_pole → mt_MSbar(mt) using QCD 2-loop matching
    - yb, ytau, lambda, v as explicit inputs (or deterministic fallbacks if not provided)

    Two routes are produced for reproducibility diagnostics:
    - route_2loop: 2-loop gauge-only running with QCD nf thresholds + 2-loop αs matching (via `run_sm_gauge_only_2loop_thresholds`)
    - route_1loop: 1-loop gauge-only running with QCD nf thresholds (via `run_sm_gauge_only_1loop_thresholds`)
    """
    mu_mz = float(sm_inputs_mz["mu_GeV"])
    mt_pole = float(sm_inputs_mz.get("mt_GeV", 172.76))
    enforce_no_legacy_sm_rge_above_mz(
        mu_start_GeV=mu_mz, mu_end_GeV=mt_pole, caller="rge_sm.sm_boundary_conditions_at_mt"
    )
    mc = float(sm_inputs_mz.get("mc_GeV", 1.27))
    mb = float(sm_inputs_mz.get("mb_GeV", 4.18))
    mH = float(sm_inputs_mz.get("mH_GeV", 125.25))

    alpha_em_inv = float(sm_inputs_mz["alpha_em_inv"])
    sin2 = float(sm_inputs_mz["sin2_thetaW"])
    alpha_s_mz = float(sm_inputs_mz["alpha_s"])

    # Optional explicit MSbar-ish boundary knobs (kept here as explicit inputs, not fitted).
    #
    # If yb_mt / ytau_mt are not provided, we derive them deterministically from PDG-style masses:
    # - mb(mt) from mb(mb) via LO QCD running (nf=5), using αs computed in the same gauge-running policy.
    # - yτ(mt) from mτ/v (no QED running in this helper; see msbar_matching_map for explicit policy notes).
    v_ev = float(sm_inputs_mz.get("v_ev_GeV", 246.0))
    yb_mt_in = sm_inputs_mz.get("yb_mt", None)
    ytau_mt_in = sm_inputs_mz.get("ytau_mt", None)

    if "lambda_mt" in sm_inputs_mz:
        lam_mt = float(sm_inputs_mz["lambda_mt"])
        lam_source = "input"
    else:
        lam_mt = float((mH * mH) / (2.0 * v_ev * v_ev))
        lam_source = "tree_from_mH_and_v"

    # MZ inputs object (for gauge-coupling initialization)
    inp = SmMzInputs(mu_GeV=mu_mz, alpha_em_inv=alpha_em_inv, sin2_thetaW=sin2, alpha_s=alpha_s_mz)
    g_mz_gut = gauge_couplings_from_mz_inputs(inp)

    # Route A: 2-loop gauge-only (with QCD thresholds + 2-loop αs matching)
    g_mt_gut_2l = run_sm_gauge_only_2loop_thresholds(
        mu_start_GeV=mu_mz,
        mu_end_GeV=mt_pole,
        g_start=g_mz_gut,
        mc_GeV=mc,
        mb_GeV=mb,
        mt_GeV=mt_pole,
        apply_alpha3_matching=True,
    )
    gY_mt_2l = gY_from_g1_gut(g_mt_gut_2l[0])
    g2_mt_2l = float(g_mt_gut_2l[1])
    g3_mt_2l = float(g_mt_gut_2l[2])
    alpha_s_mt_2l = alpha_from_g(g3_mt_2l)
    mt_msbar_mt_2l = pole_mass_to_msbar_mass_qcd_2loop(pole_mass_GeV=mt_pole, alpha_s_mu=alpha_s_mt_2l, n_light=5)
    yt_mt_2l = float(math.sqrt(2.0) * mt_msbar_mt_2l / v_ev)

    # --- yb(mt), ytau(mt) (explicit inputs or deterministic derivations) ---
    if yb_mt_in is not None:
        yb_mt = float(yb_mt_in)
        yb_source = "input"
    else:
        # nf=5 αs for mb running: use the same 2-loop gauge-only policy but *without* top-threshold matching.
        # (Between MZ and mt, nf=5; this yields αs(mt) in the nf=5 effective theory.)
        g_mt_nf5 = run_sm_gauge_only_2loop_thresholds(
            mu_start_GeV=mu_mz,
            mu_end_GeV=mt_pole,
            g_start=g_mz_gut,
            mc_GeV=mc,
            mb_GeV=mb,
            mt_GeV=mt_pole,
            apply_alpha3_matching=False,
        )
        g_mb_nf5 = run_sm_gauge_only_2loop_thresholds(
            mu_start_GeV=mu_mz,
            mu_end_GeV=mb,
            g_start=g_mz_gut,
            mc_GeV=mc,
            mb_GeV=mb,
            mt_GeV=mt_pole,
            apply_alpha3_matching=False,
        )
        alpha_s_mt_nf5 = alpha_from_g(float(g_mt_nf5[2]))
        alpha_s_mb_nf5 = alpha_from_g(float(g_mb_nf5[2]))
        mb_mt = msbar_mass_run_qcd_1loop(m_mu0_GeV=mb, alpha_s_mu0=alpha_s_mb_nf5, alpha_s_mu1=alpha_s_mt_nf5, nf=5)
        yb_mt = float(math.sqrt(2.0) * mb_mt / v_ev)
        yb_source = "derived_from_mb_mb_via_1loop_qcd_nf5"

    if ytau_mt_in is not None:
        ytau_mt = float(ytau_mt_in)
        ytau_source = "input"
    else:
        # Deterministic proxy: yτ ≈ √2 mτ/v (no QED running here).
        try:
            from pathlib import Path
            import json

            lep_path = Path(__file__).resolve().parent / "data" / "lepton_masses_pdg.json"
            lep = json.loads(lep_path.read_text(encoding="utf-8")) if lep_path.exists() else {}
            mtau = float(lep.get("masses", {}).get("tau", {}).get("mean", 1.77686))
        except Exception:
            mtau = 1.77686
        ytau_mt = float(math.sqrt(2.0) * mtau / v_ev)
        ytau_source = "derived_from_mtau_over_v_tree"

    # Route B: 1-loop gauge-only (with QCD thresholds)
    g_mt_gut_1l = run_sm_gauge_only_1loop_thresholds(
        mu_start_GeV=mu_mz,
        mu_end_GeV=mt_pole,
        g_start=g_mz_gut,
        mc_GeV=mc,
        mb_GeV=mb,
        mt_GeV=mt_pole,
    )
    gY_mt_1l = gY_from_g1_gut(g_mt_gut_1l[0])
    g2_mt_1l = float(g_mt_gut_1l[1])
    g3_mt_1l = float(g_mt_gut_1l[2])
    alpha_s_mt_1l = alpha_from_g(g3_mt_1l)
    mt_msbar_mt_1l = pole_mass_to_msbar_mass_qcd_2loop(pole_mass_GeV=mt_pole, alpha_s_mu=alpha_s_mt_1l, n_light=5)
    yt_mt_1l = float(math.sqrt(2.0) * mt_msbar_mt_1l / v_ev)

    route_2loop = {
        "mu_mz_GeV": mu_mz,
        "mu_mt_GeV": mt_pole,
        "g1_gut_mt": float(g_mt_gut_2l[0]),
        "gY_mt": gY_mt_2l,
        "g2_mt": g2_mt_2l,
        "g3_mt": g3_mt_2l,
        "alpha_s_mt": alpha_s_mt_2l,
        "mt_msbar_mt_GeV": mt_msbar_mt_2l,
        "yt_mt": yt_mt_2l,
        "yb_mt": yb_mt,
        "yb_source": yb_source,
        "ytau_mt": ytau_mt,
        "ytau_source": ytau_source,
        "lambda_mt": lam_mt,
        "v_ev_GeV": v_ev,
        "lambda_source": lam_source,
    }
    route_1loop = {
        "mu_mz_GeV": mu_mz,
        "mu_mt_GeV": mt_pole,
        "g1_gut_mt": float(g_mt_gut_1l[0]),
        "gY_mt": gY_mt_1l,
        "g2_mt": g2_mt_1l,
        "g3_mt": g3_mt_1l,
        "alpha_s_mt": alpha_s_mt_1l,
        "mt_msbar_mt_GeV": mt_msbar_mt_1l,
        "yt_mt": yt_mt_1l,
        "yb_mt": yb_mt,
        "yb_source": yb_source,
        "ytau_mt": ytau_mt,
        "ytau_source": ytau_source,
        "lambda_mt": lam_mt,
        "v_ev_GeV": v_ev,
        "lambda_source": lam_source,
    }

    diffs = {
        "gY_mt": float(route_2loop["gY_mt"] - route_1loop["gY_mt"]),
        "g2_mt": float(route_2loop["g2_mt"] - route_1loop["g2_mt"]),
        "g3_mt": float(route_2loop["g3_mt"] - route_1loop["g3_mt"]),
        "yt_mt": float(route_2loop["yt_mt"] - route_1loop["yt_mt"]),
    }

    return SmBoundaryAtMt(
        mu_mt_GeV=float(mt_pole),
        scheme="MSbar (gauge MSbar-at-MZ inputs; mt pole→MSbar(mt) via QCD 2-loop)",
        conventions={
            "hypercharge_definition": "Q = T3 + Y",
            "g1_gut_over_gY": f"{g1_gut_over_gY():.15g}",
            "g1_gut_definition": "g1_GUT = sqrt(5/3) * gY",
            "gY_definition": "gY ≡ g′ (SM hypercharge coupling)",
            "lambda_source": lam_source,
        },
        route_2loop=route_2loop,
        route_1loop=route_1loop,
        diffs=diffs,
    )

