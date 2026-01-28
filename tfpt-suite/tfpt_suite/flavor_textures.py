from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np


PhaseMode = Literal["2pi_delta", "delta_rad", "koide_pi_over_12"]


@dataclass(frozen=True)
class TextureConfig:
    """
    Minimal configuration for the Z3-flavor Yukawa texture layer.

    The current paper notes specify a *circulant + diagonal* form:
      Y^(y) = y_* [ C(δ) + a_y varphi0 D + b c3 I ],
      D = diag(1,0,-1),
    but do not pin down all sector coefficients. Those are treated as explicit inputs here.
    """

    phase_mode: PhaseMode
    b: float
    a_u: float
    a_d: float
    a_e: float
    a_nu: float


@dataclass(frozen=True)
class TextureCoefficients:
    """
    Fixed sector coefficients for the Z3 texture layer.

    These values are treated as *theory inputs* for v1.07SM/v1.06-style TFPT flavor:
    - a_u = 2/3 (up-sector cusp)
    - a_d = 1   (down-sector slope used in the notes)
    - a_e = 1   (charged-lepton slope)
    - a_nu ≈ 1.103 (neutrino slope quoted in update_tfptv1_07sm)
    - b = 1     (minimal topological offset multiplying c3·I)

    They are intentionally not optimized/fitted in the suite.
    """

    a_u: float
    a_d: float
    a_e: float
    a_nu: float
    b: float
    source: str


def coefficients_v107sm() -> TextureCoefficients:
    """
    Return the fixed texture coefficients as stated/used in the TFPT v1.06/v1.07SM notes.
    """
    return TextureCoefficients(
        a_u=2.0 / 3.0,
        a_d=1.0,
        a_e=1.0,
        a_nu=1.103,
        b=1.0,
        source="update_tfptv1_07sm.tex (v1.07SM) + paper_v1_06_01_09_2025.tex (v1.06) conventions",
    )


def circulant(c0: complex, c1: complex, c2: complex) -> np.ndarray:
    """
    3x3 circulant matrix from its first row (c0,c1,c2):

      circ(c0,c1,c2) =
        [[c0, c1, c2],
         [c2, c0, c1],
         [c1, c2, c0]]
    """
    return np.array([[c0, c1, c2], [c2, c0, c1], [c1, c2, c0]], dtype=complex)


def theta_of_delta(delta: float, *, phase_mode: PhaseMode) -> float:
    """
    Map the Möbius deformation parameter δ to a complex phase θ used in C(δ).

    The docs are ambiguous on whether δ is used as a literal angle or a fractional phase.
    We therefore expose explicit modes:
    - 2pi_delta: θ = 2π δ   (recommended default)
    - delta_rad: θ = δ
    - koide_pi_over_12: θ = π/12   (Koide alignment phase; ignores δ in the circulant phase)
    """
    if phase_mode == "2pi_delta":
        return float(2.0 * np.pi * float(delta))
    if phase_mode == "delta_rad":
        return float(delta)
    if phase_mode == "koide_pi_over_12":
        return float(np.pi / 12.0)
    raise ValueError(f"Unsupported phase_mode: {phase_mode}")


def C_delta(delta: float, *, phase_mode: PhaseMode) -> np.ndarray:
    """
    Minimal Hermitian Z3-circulant matrix built from a single phase:

      ζ(δ) = exp(i θ(δ))
      C(δ) = circ(1, ζ, ζ*) =
        [[1, ζ, ζ*],
         [ζ*, 1, ζ],
         [ζ, ζ*, 1]]
    """
    theta = theta_of_delta(delta, phase_mode=phase_mode)
    zeta = complex(np.cos(theta), np.sin(theta))
    return circulant(1.0, zeta, np.conj(zeta))


def D_diag() -> np.ndarray:
    """
    D = diag(1,0,-1) (the diagonal Z3-breaking direction used in the v1.07SM texture formula).
    """
    return np.diag([1.0, 0.0, -1.0]).astype(complex)


def scale_y_star_to_match_sigma_max(*, target_y3: float, base: np.ndarray) -> float:
    """
    Choose y_* such that sigma_max(y_* base) = target_y3.
    """
    if target_y3 <= 0:
        raise ValueError("target_y3 must be positive")
    s = np.linalg.svd(base, compute_uv=False)
    smax = float(np.max(s)) if s.size else 0.0
    if smax <= 0:
        raise ValueError("base matrix has zero singular values; cannot scale")
    return float(target_y3) / smax


def yukawa_texture_matrix(
    *,
    delta: float,
    varphi0: float,
    c3: float,
    a_y: float,
    b: float,
    y_star: float,
    phase_mode: PhaseMode,
    theta_override_rad: float | None = None,
) -> np.ndarray:
    """
    Construct a Yukawa matrix from the stated "circulant + diagonal + identity" texture.

    Optional wiring hook:
    - If `theta_override_rad` is provided, the circulant phase is taken as ζ = exp(i·theta_override_rad)
      instead of ζ(delta; phase_mode). This allows discrete topology-phase atoms (Wilson-line holonomies)
      to be injected downstream without introducing any continuous fitting parameter.
    """
    if theta_override_rad is None:
        C = C_delta(delta, phase_mode=phase_mode)
    else:
        theta = float(theta_override_rad)
        zeta = complex(np.cos(theta), np.sin(theta))
        C = circulant(1.0, zeta, np.conj(zeta))
    D = D_diag()
    I = np.eye(3, dtype=complex)
    return float(y_star) * (C + float(a_y) * float(varphi0) * D + float(b) * float(c3) * I)


def left_unitary_from_yukawa(Y: np.ndarray) -> np.ndarray:
    """
    Deterministic left-unitary diagonalizer of H = Y Y†.

    We sort eigenvalues ascending (light→heavy) and fix column phases so that the
    largest-magnitude component of each eigenvector is real-positive.
    """
    H = Y @ Y.conj().T
    w, U = np.linalg.eigh(H)
    idx = np.argsort(w.real)
    U = U[:, idx]

    # Deterministic rephasing: for each column, make the largest component real-positive.
    U2 = U.copy()
    for k in range(U2.shape[1]):
        col = U2[:, k]
        i = int(np.argmax(np.abs(col)))
        a = col[i]
        if abs(a) == 0:
            continue
        phase = a / abs(a)
        U2[:, k] = col / phase
        # ensure positive real
        if U2[i, k].real < 0:
            U2[:, k] *= -1
    return U2.astype(complex)

