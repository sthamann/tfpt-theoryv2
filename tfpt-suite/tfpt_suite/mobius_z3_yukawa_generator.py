from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np

from tfpt_suite.flavor_textures import PhaseMode, left_unitary_from_yukawa, theta_of_delta


S13Mode = Literal["A_lam3_over_3", "A_lam3_times_1_minus_delta"]
DeltaCkmMode = Literal["pi_times_delta", "pi_times_1_minus_delta", "2pi_times_delta", "koide_pi_over_12", "external_delta_cp_override"]


@dataclass(frozen=True)
class CkmConstruction:
    """
    Deterministic CKM construction used by the Möbius/Z3 Yukawa generator.

    This encodes the v1.07SM “cold pass” idea:
      λ from φ0, A from Z3 cusp slopes, s23 ~ A λ^2, s13 ~ O(λ^3),
    and an explicit CP-phase convention derived from the Möbius deformation δ.
    """

    lam: float
    A: float
    s12: float
    s23: float
    s13: float
    delta_cp_rad: float
    s13_mode: S13Mode
    delta_mode: DeltaCkmMode


@dataclass(frozen=True)
class YukawaGenerationMeta:
    scheme: str
    reference_scale_GeV: float
    delta_used: float
    phase_mode: str
    ckm: CkmConstruction
    ratios: dict[str, float]


@dataclass(frozen=True)
class QuarkPhaseMap:
    """
    Explicit naming map to avoid conflating distinct "phases" in the flavor pipeline.

    - delta_mobius: δ extracted from lepton masses (τ/μ) via Möbius map inversion.
      (Paper v2.5: Sec. 7.1 "Universal Phase δ", Eq. δ_M; and Sec. 7.2 "Möbius Mass Relations".)
    - theta_texture_rad: θ entering the Z3 circulant C(δ) via ζ = exp(i θ).
      This is a texture-layer phase and depends on the explicit `phase_mode`.
    - delta_ckm_rad: PDG CKM Dirac CP phase δ_CKM used in V_CKM.
      The suite currently uses an explicit convention `delta_mode` (v1.07SM-style “cold pass”);
      a full topology→phase map derivation is still pending and tracked by the QFT ledger.
    """

    delta_mobius: float
    delta_used: float
    delta_source: str
    theta_texture_rad: float
    phase_mode: PhaseMode
    delta_ckm_rad: float
    delta_mode: DeltaCkmMode
    source: str
    notes: list[str]


def quark_phase_map(
    *,
    delta_mobius: float,
    delta_used: float,
    delta_source: str,
    phase_mode: PhaseMode,
    delta_mode: DeltaCkmMode,
    delta_ckm_override_rad: float | None = None,
) -> QuarkPhaseMap:
    """
    Compute the explicit phase map (δ_Mobius, θ_texture, δ_CKM) used by the suite.

    Note: `delta_used` is the δ actually fed into the texture/generator; it may differ from
    `delta_mobius` if `delta_source` chooses δ_* instead of δ_M.
    """
    d_used = float(delta_used)

    # Texture phase: ζ = exp(i θ(δ_used))
    theta = float(theta_of_delta(d_used, phase_mode=phase_mode))

    # CKM Dirac phase convention (explicit; tracked as a postulate until formalized).
    #
    # If an explicit override is provided (e.g. from topology_phase_map), prefer it and mark the mode accordingly.
    if delta_ckm_override_rad is not None:
        delta_ckm = float(delta_ckm_override_rad)
        delta_mode = "external_delta_cp_override"
    else:
        if delta_mode == "pi_times_delta":
            delta_ckm = float(np.pi * d_used)
        elif delta_mode == "pi_times_1_minus_delta":
            delta_ckm = float(np.pi * (1.0 - d_used))
        elif delta_mode == "2pi_times_delta":
            delta_ckm = float(2.0 * np.pi * d_used)
        elif delta_mode == "koide_pi_over_12":
            delta_ckm = float(np.pi / 12.0)
        elif delta_mode == "external_delta_cp_override":
            # Backwards-compatible: allow callers to explicitly declare override mode.
            delta_ckm = float(np.pi * (1.0 - d_used))
        else:
            raise ValueError(f"Unsupported delta_mode: {delta_mode}")
    delta_ckm = float(delta_ckm % (2.0 * np.pi))

    src = "paper v2.5 Sec. 7.1–7.3 (Flavor Structure) + update_tfptv1_07sm.tex (CKM cold-pass phase convention)"
    notes = [
        "delta_mobius is the τ/μ-extracted Möbius deformation parameter (δ_M).",
        "theta_texture_rad is a Z3-texture layer parameter (ζ = exp(i θ)); it is not the CKM Dirac phase.",
        "delta_ckm_rad is the PDG CKM Dirac CP phase δ_CKM used in V_CKM; current mapping is an explicit suite convention.",
        f"delta_source={delta_source} (delta_used may differ from delta_mobius if δ_* is chosen).",
    ]
    return QuarkPhaseMap(
        delta_mobius=float(delta_mobius),
        delta_used=d_used,
        delta_source=str(delta_source),
        theta_texture_rad=theta,
        phase_mode=phase_mode,
        delta_ckm_rad=delta_ckm,
        delta_mode=delta_mode,
        source=src,
        notes=notes,
    )


def mobius_map(y: float, delta: float) -> float:
    """
    Möbius map used in the TFPT flavor notes:
      M_y(δ) = (y + δ)/(y - δ).
    """
    y = float(y)
    delta = float(delta)
    if abs(y - delta) < 1e-15:
        raise ZeroDivisionError("Möbius map singular at y=delta")
    return float((y + delta) / (y - delta))


def _ckm_from_angles(*, s12: float, s23: float, s13: float, delta_cp_rad: float) -> np.ndarray:
    """
    PDG convention:
      V = R23 * U13(δ) * R12
    """
    s12 = float(s12)
    s23 = float(s23)
    s13 = float(s13)
    d = float(delta_cp_rad)

    c12 = float(np.sqrt(max(0.0, 1.0 - s12 * s12)))
    c23 = float(np.sqrt(max(0.0, 1.0 - s23 * s23)))
    c13 = float(np.sqrt(max(0.0, 1.0 - s13 * s13)))

    e_minus = np.exp(-1j * d)
    e_plus = np.exp(1j * d)

    R23 = np.array([[1.0, 0.0, 0.0], [0.0, c23, s23], [0.0, -s23, c23]], dtype=complex)
    U13 = np.array([[c13, 0.0, s13 * e_minus], [0.0, 1.0, 0.0], [-s13 * e_plus, 0.0, c13]], dtype=complex)
    R12 = np.array([[c12, s12, 0.0], [-s12, c12, 0.0], [0.0, 0.0, 1.0]], dtype=complex)
    return (R23 @ U13 @ R12).astype(complex)


def _tfpt_lambda_from_varphi0(varphi0: float) -> float:
    """
    v1.07SM: λ = sqrt(varphi0) * (1 - varphi0/2).
    """
    v = float(varphi0)
    return float(np.sqrt(v) * (1.0 - 0.5 * v))


def _ckm_construction_v107sm(
    *,
    varphi0: float,
    delta_used: float,
    s13_mode: S13Mode = "A_lam3_times_1_minus_delta",
    delta_mode: DeltaCkmMode = "pi_times_1_minus_delta",
    delta_cp_override_rad: float | None = None,
) -> CkmConstruction:
    lam = _tfpt_lambda_from_varphi0(varphi0)

    # v1.07SM: sector slopes a_u=2, a_d=1 => A=(a_u+a_d)/6 = 5/6
    A = float(5.0 / 6.0)

    s12 = float(lam)
    s23 = float(A * (lam**2))

    if s13_mode == "A_lam3_over_3":
        s13 = float(A * (lam**3) / 3.0)
    elif s13_mode == "A_lam3_times_1_minus_delta":
        # Empirically (and neatly) matches |Vub| scale if δ is fixed from τ/μ:
        #   |Vub| ~ A λ^3 (1-δ)
        s13 = float(A * (lam**3) * (1.0 - float(delta_used)))
    else:
        raise ValueError(f"Unsupported s13_mode: {s13_mode}")

    if delta_cp_override_rad is not None:
        delta_cp = float(delta_cp_override_rad)
        delta_mode = "external_delta_cp_override"
    else:
        if delta_mode == "pi_times_delta":
            delta_cp = float(np.pi * float(delta_used))
        elif delta_mode == "pi_times_1_minus_delta":
            # Gives δ_CP ≈ 70° for δ≈0.608 (τ/μ), close to observed CKM δ_CP.
            delta_cp = float(np.pi * (1.0 - float(delta_used)))
        elif delta_mode == "2pi_times_delta":
            delta_cp = float(2.0 * np.pi * float(delta_used))
        elif delta_mode == "koide_pi_over_12":
            delta_cp = float(np.pi / 12.0)
        elif delta_mode == "external_delta_cp_override":
            # Backwards-compatible: allow callers to explicitly declare override mode.
            delta_cp = float(np.pi * (1.0 - float(delta_used)))
        else:
            raise ValueError(f"Unsupported delta_mode: {delta_mode}")

    # wrap to [0,2π)
    delta_cp = float(delta_cp % (2.0 * np.pi))

    return CkmConstruction(
        lam=lam,
        A=A,
        s12=s12,
        s23=s23,
        s13=s13,
        delta_cp_rad=delta_cp,
        s13_mode=s13_mode,
        delta_mode=delta_mode,
    )


def _yukawa_ratios_from_mobius(*, delta_used: float) -> dict[str, float]:
    """
    Build the *dimensionless* hierarchy ratios referenced in update_tfptv1_07sm.

    Returns ratios for:
    - down: ms/md, mb/ms
    - up: mc/mu, mt/mc
    - charged leptons: mtau/mmu (by construction of δ), mmu/me (uses 1/3 cusp)
    """
    d = float(delta_used)
    M1 = mobius_map(1.0, d)
    M1_3 = mobius_map(1.0 / 3.0, d)
    M2_3 = mobius_map(2.0 / 3.0, d)

    r_tau_mu = float(M1**2)  # √(mτ/mμ)=M1(δ)
    r_s_d = float(M1**2)  # √(ms/md)=M1(δ)
    r_b_s = float((M1 * (1.0 + d)) ** 2)

    r_c_u = float(M2_3**2)
    r_t_c = float(((2.0 / 3.0) / ((2.0 / 3.0) - d)) ** 2)

    r_mu_e = float((M1 * abs(M1_3)) ** 2)

    return {
        "m_tau_over_m_mu": r_tau_mu,
        "m_mu_over_m_e": r_mu_e,
        "m_s_over_m_d": r_s_d,
        "m_b_over_m_s": r_b_s,
        "m_c_over_m_u": r_c_u,
        "m_t_over_m_c": r_t_c,
    }


def _ckm_from_yukawas(Yu: np.ndarray, Yd: np.ndarray) -> np.ndarray:
    Uu = left_unitary_from_yukawa(Yu)
    Ud = left_unitary_from_yukawa(Yd)
    return (Uu.conj().T @ Ud).astype(complex)


def generate_quark_yukawas_mt(
    *,
    varphi0: float,
    delta_used: float,
    yt_mt: float,
    yb_mt: float,
    scheme: str = "MSbar",
    reference_scale_GeV: float = 172.76,
    s13_mode: S13Mode = "A_lam3_times_1_minus_delta",
    delta_mode: DeltaCkmMode = "pi_times_1_minus_delta",
    delta_cp_override_rad: float | None = None,
) -> tuple[np.ndarray, np.ndarray, YukawaGenerationMeta]:
    """
    Construct (Yu,Yd) at the mt boundary using:
    - eigenvalue hierarchies from Möbius relations (δ fixed from τ/μ),
    - mixing angles from v1.07SM-style “cold pass” CKM,
    - and a deterministic choice of basis:
        Yu = diag(y_u,y_c,y_t),
        Yd = V_ckm * diag(y_d,y_s,y_b).
    """
    ratios = _yukawa_ratios_from_mobius(delta_used=delta_used)

    yt = float(yt_mt)
    yb = float(yb_mt)
    if yt <= 0 or yb <= 0:
        raise ValueError("yt_mt and yb_mt must be positive")

    # Up-type: y_t/y_c = m_t/m_c, y_c/y_u = m_c/m_u
    y_c = yt / float(ratios["m_t_over_m_c"])
    y_u = y_c / float(ratios["m_c_over_m_u"])

    # Down-type: y_b/y_s = m_b/m_s, y_s/y_d = m_s/m_d
    y_s = yb / float(ratios["m_b_over_m_s"])
    y_d = y_s / float(ratios["m_s_over_m_d"])

    Yu = np.diag([y_u, y_c, yt]).astype(complex)
    Yd_diag = np.diag([y_d, y_s, yb]).astype(complex)

    ckm = _ckm_construction_v107sm(
        varphi0=varphi0,
        delta_used=delta_used,
        s13_mode=s13_mode,
        delta_mode=delta_mode,
        delta_cp_override_rad=delta_cp_override_rad,
    )
    V = _ckm_from_angles(s12=ckm.s12, s23=ckm.s23, s13=ckm.s13, delta_cp_rad=ckm.delta_cp_rad)
    Yd = (V @ Yd_diag).astype(complex)

    meta = YukawaGenerationMeta(
        scheme=str(scheme),
        reference_scale_GeV=float(reference_scale_GeV),
        delta_used=float(delta_used),
        phase_mode="CKM(PDG) from TFPT invariants; see CkmConstruction fields",
        ckm=ckm,
        ratios=ratios,
    )
    return Yu, Yd, meta

