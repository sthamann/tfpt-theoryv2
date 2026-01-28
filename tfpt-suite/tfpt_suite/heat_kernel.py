from __future__ import annotations

from dataclasses import dataclass

from mpmath import mp


@dataclass(frozen=True)
class LaplaceTypeBlock:
    """
    One Laplace-type operator block Δ = -∇^2 + E acting on a bundle of rank `rank`.

    This is a minimal, machine-readable input contract for heat-kernel computations.
    For constant-curvature bookkeeping we parameterize:

      E = (E_over_R) * R
      Ω_{μν}Ω^{μν} = (Omega_sq_over_R2) * R^2

    where E_over_R and Omega_sq_over_R2 are endomorphisms on the bundle.
    They may be scalars (mp.mpf) or mp.matrix.
    """

    name: str
    rank: int
    statistics: str  # "boson" | "fermion" | "ghost"
    E_over_R: object  # mp.mpf or mp.matrix
    Omega_sq_over_R2: object  # mp.mpf or mp.matrix


def _is_matrix(x: object) -> bool:
    return isinstance(x, mp.matrix)


def _eye(n: int) -> mp.matrix:
    return mp.matrix([[mp.mpf(1) if i == j else mp.mpf(0) for j in range(n)] for i in range(n)])


def _trace(x: object, *, rank: int) -> mp.mpf:
    if _is_matrix(x):
        M = x  # type: ignore[assignment]
        return mp.fsum([M[i, i] for i in range(M.rows)])
    # scalar interpreted as scalar*I
    return mp.mpf(x) * mp.mpf(rank)


def _square(x: object, *, rank: int) -> object:
    if _is_matrix(x):
        M = x  # type: ignore[assignment]
        return M * M
    return mp.mpf(x) ** 2


def a2_R2_coeff_constant_curvature_4d(*, block: LaplaceTypeBlock) -> mp.mpf:
    """
    Return the *curly-brace* coefficient of R^2 in the local a2 integrand in 4D,
    on a maximally symmetric background, for a Laplace-type block.

    Conventions (dropping total derivatives):

      a2 ⊃ (1/360) tr( 5 R^2 - 2 Ric^2 + 2 Riem^2 + 60 R E + 180 E^2 + 30 Ω^2 )

    For 4D maximally symmetric (constant curvature):
      Ric^2 = R^2/4,  Riem^2 = R^2/6
    hence
      (1/360)(5R^2 - 2Ric^2 + 2Riem^2) = (29/2160) R^2.

    With E = (E_over_R) R and Ω^2 = (Omega_sq_over_R2) R^2, the R^2 coefficient is:

      tr( (29/2160) I + (1/6) E_over_R + (1/2) E_over_R^2 + (1/12) Omega_sq_over_R2 )
    """
    n = int(block.rank)
    I = _eye(n)

    # Normalize inputs to endomorphisms (mp.matrix or scalar)
    EoR = block.E_over_R
    OoR2 = block.Omega_sq_over_R2

    # Build the endomorphism inside tr(...)
    term_geom = (mp.mpf(29) / mp.mpf(2160)) * I
    term_RE = (mp.mpf(1) / mp.mpf(6)) * (EoR if _is_matrix(EoR) else mp.mpf(EoR))  # type: ignore[operator]
    term_E2 = (mp.mpf(1) / mp.mpf(2)) * (_square(EoR, rank=n) if _is_matrix(EoR) else mp.mpf(EoR) ** 2)  # type: ignore[operator]
    term_O2 = (mp.mpf(1) / mp.mpf(12)) * (OoR2 if _is_matrix(OoR2) else mp.mpf(OoR2))  # type: ignore[operator]

    if _is_matrix(term_RE):
        endo = term_geom + term_RE + term_E2 + term_O2  # type: ignore[operator]
        return _trace(endo, rank=n)

    # scalar case: endo is scalar, trace multiplies by rank
    endo_scalar = (mp.mpf(29) / mp.mpf(2160)) + mp.mpf(term_RE) + mp.mpf(term_E2) + mp.mpf(term_O2)
    return endo_scalar * mp.mpf(n)


def beta_R2_from_a2_R2_coeff_4d(*, a2_R2_coeff_curly: mp.mpf, prefactor: mp.mpf = mp.mpf(1) / 2) -> mp.mpf:
    """
    Map a *curly-brace* a2 coefficient to the R^2 coefficient β_R2 in:

      Γ_1 = prefactor * Tr ln Δ  ⇒  Γ_1 ⊃ ∫ √g β_R2 R^2

    Using:
      Γ_1(local) ⊃ prefactor * (4π)^(-2) * a2

    So:
      β_R2 = prefactor * (4π)^(-2) * a2_R2_coeff_curly
           = prefactor * a2_R2_coeff_curly / (16 π^2).
    """
    return mp.mpf(prefactor) * a2_R2_coeff_curly / (mp.mpf(16) * (mp.pi**2))

