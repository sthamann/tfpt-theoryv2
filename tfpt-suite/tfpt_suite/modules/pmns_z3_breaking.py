from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _plot_pmns_z3_breaking_variants(
    *,
    out_dir: Path,
    variants: list[VariantResult],
    eps_deg: float,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"pmns_z3_breaking_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        names = [v.variant for v in variants]
        dth = [float(v.delta_theta23_deg) for v in variants]
        dd = [float(v.delta_delta_deg) for v in variants]

        x = np.arange(len(names), dtype=float)
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 6), sharex=True)

        ax1.bar(x, dth, color="#1f77b4", alpha=0.9)
        ax1.axhline(0.0, color="black", lw=1.0, alpha=0.7)
        ax1.set_ylabel(r"$\Delta\theta_{23}$ [deg]")
        ax1.set_title(rf"PMNS Z3-breaking variants (ε≈{eps_deg:.3f}°)")
        ax1.grid(True, axis="y", ls=":", alpha=0.4)

        ax2.bar(x, dd, color="#ff7f0e", alpha=0.9)
        ax2.axhline(0.0, color="black", lw=1.0, alpha=0.7)
        ax2.set_ylabel(r"$\Delta\delta_{CP}$ [deg]")
        ax2.grid(True, axis="y", ls=":", alpha=0.4)

        ax2.set_xticks(x)
        ax2.set_xticklabels(names, rotation=25, ha="right")
        fig.tight_layout()

        path = out_dir / "pmns_z3_breaking.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["pmns_z3_breaking_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


def _rot12(theta: float) -> np.ndarray:
    c = float(np.cos(theta))
    s = float(np.sin(theta))
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]], dtype=complex)


def _rot23(theta: float) -> np.ndarray:
    c = float(np.cos(theta))
    s = float(np.sin(theta))
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]], dtype=complex)


def _u13(theta: float, delta: float) -> np.ndarray:
    c = float(np.cos(theta))
    s = float(np.sin(theta))
    e_minus = np.exp(-1j * delta)
    e_plus = np.exp(1j * delta)
    return np.array(
        [[c, 0.0, s * e_minus], [0.0, 1.0, 0.0], [-s * e_plus, 0.0, c]],
        dtype=complex,
    )


def _pmns_pdg(theta12: float, theta13: float, theta23: float, delta: float) -> np.ndarray:
    return _rot23(theta23) @ _u13(theta13, delta) @ _rot12(theta12)


@dataclass(frozen=True)
class MixingAngles:
    theta12_deg: float
    theta13_deg: float
    theta23_deg: float
    delta_cp_deg: float


def _extract_pdg_angles(U: np.ndarray) -> MixingAngles:
    """
    Extract (θ12, θ13, θ23, δ) from a unitary PMNS matrix in the PDG convention using
    rephasing-invariant relations (|U_e3|, |U_e2|, |U_μ3|, Jarlskog).

    This ignores Majorana phases (not needed for the angles/δ extraction).
    """
    # s13
    s13 = float(np.abs(U[0, 2]))
    s13 = min(max(s13, 0.0), 1.0)
    c13 = float(np.sqrt(max(0.0, 1.0 - s13 * s13)))

    # s12 from |U_e2| = s12 c13
    s12 = float(np.abs(U[0, 1]) / c13) if c13 > 0 else 0.0
    s12 = min(max(s12, 0.0), 1.0)
    c12 = float(np.sqrt(max(0.0, 1.0 - s12 * s12)))

    # s23 from |U_mu3| = s23 c13
    s23 = float(np.abs(U[1, 2]) / c13) if c13 > 0 else 0.0
    s23 = min(max(s23, 0.0), 1.0)
    c23 = float(np.sqrt(max(0.0, 1.0 - s23 * s23)))

    # Jarlskog invariant
    J = float(np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])))

    denom_sin = c12 * c23 * (c13**2) * s12 * s23 * s13
    sin_delta = J / denom_sin if denom_sin != 0 else 0.0
    sin_delta = float(np.clip(sin_delta, -1.0, 1.0))

    # cos δ from |U_{μ1}|^2 relation:
    # |U_{μ1}|^2 = s12^2 c23^2 + c12^2 s23^2 s13^2 + 2 s12 c12 c23 s23 s13 cos δ
    Um1_sq = float(np.abs(U[1, 0]) ** 2)
    denom_cos = 2.0 * s12 * c12 * c23 * s23 * s13
    if denom_cos != 0:
        cos_delta = (Um1_sq - (s12**2) * (c23**2) - (c12**2) * (s23**2) * (s13**2)) / denom_cos
        cos_delta = float(np.clip(cos_delta, -1.0, 1.0))
    else:
        cos_delta = 1.0

    delta = float(np.arctan2(sin_delta, cos_delta))
    if delta < 0:
        delta += 2.0 * np.pi

    return MixingAngles(
        theta12_deg=float(np.degrees(np.arcsin(s12))),
        theta13_deg=float(np.degrees(np.arcsin(s13))),
        theta23_deg=float(np.degrees(np.arcsin(s23))),
        delta_cp_deg=float(np.degrees(delta)),
    )


@dataclass(frozen=True)
class VariantResult:
    variant: str
    angles: MixingAngles
    delta_theta23_deg: float
    delta_delta_deg: float


def _z3_permutation_matrix() -> np.ndarray:
    """
    Cyclic Z3 permutation acting on family space:
      (1,2,3) -> (2,3,1)
    """
    P = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=complex)
    return P


def _z3_project_majorana(M: np.ndarray) -> np.ndarray:
    """
    Z3 projector for a Majorana mass matrix (family space):
      M_inv = (1/3) * (M + P^T M P + P^{2T} M P^2)
    """
    P = _z3_permutation_matrix()
    P2 = P @ P
    return (M + (P.T @ M @ P) + (P2.T @ M @ P2)) / 3.0


def _symmetric_basis_3x3() -> list[np.ndarray]:
    """
    Real symmetric 3x3 basis (6 elements) spanning all Majorana mass matrices.
    """
    basis: list[np.ndarray] = []
    for i in range(3):
        M = np.zeros((3, 3), dtype=complex)
        M[i, i] = 1.0
        basis.append(M)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        M = np.zeros((3, 3), dtype=complex)
        M[i, j] = 1.0
        M[j, i] = 1.0
        basis.append(M)
    return basis


def _svd_basis(mats: list[np.ndarray], *, tol: float = 1e-12) -> list[np.ndarray]:
    """
    Return an orthonormal basis for the span of mats (flattened), via SVD.
    """
    if not mats:
        return []
    X = np.stack([m.reshape(-1) for m in mats], axis=1)  # (9, n)
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    if s.size == 0:
        return []
    s0 = float(s[0])
    thr = tol * s0 if s0 > 0 else tol
    r = int(np.sum(s > thr))
    out: list[np.ndarray] = []
    for k in range(r):
        v = U[:, k].reshape(3, 3)
        v = 0.5 * (v + v.T)  # enforce symmetry numerically
        out.append(v)
    return out


def _max_abs_entry(M: np.ndarray) -> float:
    return float(np.max(np.abs(M)))


class PmnsZ3BreakingModule(TfptModule):
    module_id = "pmns_z3_breaking"
    title = "PMNS Z3-breaking scan (ε = varphi0/6; discrete operator variants)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["TFPT varphi0, gamma(0)=5/6, epsilon=varphi0/6 (paper v2.4)"],
            outputs=["baseline TM1 angles and a discrete Z3-breaking variant table (θ23, δCP shifts)"],
            formulas=[
                "sin^2(theta13) = varphi0 * exp(-gamma(0)) with gamma(0)=5/6 (paper identity table)",
                "TM1 sum rule: sin^2(theta12) = (1/3)*(1 - 2 sin^2(theta13))",
                "leading order: theta23=45°, delta=90° (Z3 symmetry)",
                "Z3-breaking scale: epsilon = varphi0/6 (fixed; no new continuous parameter)",
            ],
            validation=[
                "perturbative shifts |Δθ23| and |Δδ| are O(epsilon) in radians for the scanned variants",
            ],
            determinism="Deterministic (finite variant list).",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        gamma0 = mp.mpf(5) / 6
        sin2_theta13 = c.varphi0 * mp.e ** (-gamma0)
        theta13 = float(mp.asin(mp.sqrt(sin2_theta13)))

        sin2_theta12 = (mp.mpf(1) / 3) * (mp.mpf(1) - mp.mpf(2) * sin2_theta13)
        theta12 = float(mp.asin(mp.sqrt(sin2_theta12)))

        theta23_0 = float(np.deg2rad(45.0))
        delta0 = float(np.deg2rad(90.0))

        U0 = _pmns_pdg(theta12, theta13, theta23_0, delta0)
        base_angles = _extract_pdg_angles(U0)

        eps = float(c.varphi0 / 6)  # dimensionless; used here as a small angle proxy

        # --- Derive the Z3-invariant operator basis (and the breaking complement) ---
        P = _z3_permutation_matrix()
        basis6 = _symmetric_basis_3x3()
        inv_candidates = [_z3_project_majorana(B) for B in basis6]
        brk_candidates = [B - _z3_project_majorana(B) for B in basis6]
        inv_basis = _svd_basis(inv_candidates)
        brk_basis = _svd_basis(brk_candidates)

        def _is_z3_invariant(M: np.ndarray, *, atol: float = 1e-12) -> bool:
            return _max_abs_entry(P.T @ M @ P - M) < atol

        inv_ok = all(_is_z3_invariant(M) for M in inv_basis)

        # Discrete operator variants (left vs right action; sign choices; optional Z3 phase insertion)
        phase_Z3 = np.exp(2j * np.pi / 3)
        variants: list[tuple[str, Callable[[np.ndarray], np.ndarray]]] = [
            ("L_R23(+eps)", lambda U: _rot23(+eps) @ U),
            ("L_R23(-eps)", lambda U: _rot23(-eps) @ U),
            ("R_R23(+eps)", lambda U: U @ _rot23(+eps)),
            ("R_R23(-eps)", lambda U: U @ _rot23(-eps)),
            ("L_R23(+eps)*Z3phase", lambda U: (_rot23(+eps) @ U) @ np.diag([1.0, phase_Z3, 1.0])),
            ("L_R23(-eps)*Z3phase", lambda U: (_rot23(-eps) @ U) @ np.diag([1.0, phase_Z3, 1.0])),
        ]

        results: list[VariantResult] = []
        for name, op in variants:
            U = op(U0)
            ang = _extract_pdg_angles(U)
            results.append(
                VariantResult(
                    variant=name,
                    angles=ang,
                    delta_theta23_deg=ang.theta23_deg - base_angles.theta23_deg,
                    delta_delta_deg=((ang.delta_cp_deg - base_angles.delta_cp_deg + 540.0) % 360.0) - 180.0,
                )
            )

        # Basic perturbativity check: shifts should be ~eps (converted to degrees)
        eps_deg = float(np.degrees(eps))
        max_shift_theta23 = max(abs(r.delta_theta23_deg) for r in results)
        max_shift_delta = max(abs(r.delta_delta_deg) for r in results)

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="tfpt_z3_breaking_operator_basis_derived",
                passed=bool(inv_ok and len(inv_basis) == 2 and len(brk_basis) == 4),
                detail=f"derived Z3-invariant subspace dim={len(inv_basis)} and breaking subspace dim={len(brk_basis)} (expected 2 and 4); invariance checks={inv_ok}",
            )
        )
        checks.append(
            Check(
                check_id="epsilon_value",
                passed=bool(eps > 0 and eps < 0.02),
                detail=f"epsilon=varphi0/6 = {eps} (rad proxy), ~{eps_deg:.3f}°",
            )
        )
        checks.append(
            Check(
                check_id="perturbative_shifts",
                passed=bool(max_shift_theta23 < 5 * eps_deg and max_shift_delta < 50 * eps_deg),
                detail=f"max Δθ23={max_shift_theta23:.3f}°, max Δδ={max_shift_delta:.3f}°",
            )
        )

        lines: list[str] = []
        lines += [
            "PMNS Z3-breaking scan (ε fixed = varphi0/6)",
            "",
            "Derived operator basis (family-space Majorana mass matrices; Z3 permutation P):",
            f"- dim(invariant) = {len(inv_basis)}, dim(breaking) = {len(brk_basis)}",
            "",
            "Baseline (TM1 + Z3 leading order):",
            f"- sin^2(theta13) = varphi0*exp(-5/6) = {sin2_theta13}",
            f"- theta13 = {base_angles.theta13_deg:.4f}°",
            f"- theta12 (TM1 sum rule) = {base_angles.theta12_deg:.4f}°",
            f"- theta23 (LO) = {base_angles.theta23_deg:.4f}°",
            f"- delta_CP (LO) = {base_angles.delta_cp_deg:.4f}°",
            "",
            f"epsilon = varphi0/6 = {eps:.6g} (rad proxy) ~ {eps_deg:.3f}°",
            "",
            "Variants:",
            "variant                 theta23[deg]   deltaCP[deg]   Δtheta23[deg]   Δdelta[deg]",
        ]
        for r in results:
            lines.append(
                f"{r.variant:<22s}  {r.angles.theta23_deg:>10.4f}  {r.angles.delta_cp_deg:>12.4f}  {r.delta_theta23_deg:>12.4f}  {r.delta_delta_deg:>10.4f}"
            )

        lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This implements a finite discrete-operator scan with ε fixed by TFPT (no new continuous parameters).",
            "- The Z3-invariant vs Z3-breaking operator subspaces are derived from the family-space permutation symmetry and exposed in results.json for downstream restrictions.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"pmns_z3_breaking_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_pmns_z3_breaking_variants(out_dir=out_dir, variants=results, eps_deg=eps_deg)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "epsilon": {"eps": eps, "eps_deg": eps_deg},
                "z3_operator_basis": {
                    "P": [[str(x) for x in row] for row in P.tolist()],
                    "dim_invariant": len(inv_basis),
                    "dim_breaking": len(brk_basis),
                    "invariant_basis_max_abs_entry": [_max_abs_entry(M) for M in inv_basis],
                    "breaking_basis_max_abs_entry": [_max_abs_entry(M) for M in brk_basis],
                },
                "baseline": base_angles.__dict__,
                "variants": [
                    {
                        "variant": r.variant,
                        "angles": r.angles.__dict__,
                        "delta_theta23_deg": r.delta_theta23_deg,
                        "delta_delta_deg": r.delta_delta_deg,
                    }
                    for r in results
                ],
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

