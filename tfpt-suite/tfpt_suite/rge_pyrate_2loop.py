from __future__ import annotations

import importlib
import importlib.util
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Optional

import numpy as np
from scipy.integrate import solve_ivp

from tfpt_suite.conventions import g1_gut_over_gY
from tfpt_suite.matching import match_gauge, match_quartic, match_yukawa


ModelKind = Literal["sm_tfpt_2loop_v25", "e8_sigma_yN_2loop"]


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _workspace_root() -> Path:
    # .../wolfram_latex_attachments/tfpt-suite/tfpt_suite/<file>.py -> parents[2] is workspace root
    return Path(__file__).resolve().parents[2]


def _relpath(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(_workspace_root()))
    except Exception:
        return str(path)


def _kappa(loop: int) -> float:
    """
    PyR@TE convention:
      dX/d ln(mu) = Σ_{n>=1} kappa(n) * beta_X(nLoop=n)
    where
      kappa(n) = 1/(4π)^(2n) = 1/(16π^2)^n.
    """
    if loop < 1:
        raise ValueError("loop must be >= 1")
    return float(1.0 / (4.0 * np.pi) ** (2 * loop))


@dataclass(frozen=True)
class PyrateBetaModule:
    """
    Container for a PyR@TE-generated `RGEs.py` module imported from a PythonOutput directory.
    """

    kind: ModelKind
    pythonoutput_dir: Path
    pythonoutput_module_file: Path
    pythonoutput_module_sha256: str
    model_name_expected: str
    yaml_source: Optional[Path]
    yaml_source_sha256: Optional[str]
    rges: Any


@dataclass(frozen=True)
class ThresholdRule:
    """
    Declarative threshold bookkeeping.

    This makes every segment switch explicit, even if we currently assume continuity (no finite matching).
    """

    threshold_id: str
    scale_GeV: float
    action: str
    status: str
    affected_parameters: list[str]
    note: str


def load_pyrate_beta_module(
    *,
    kind: ModelKind,
    pythonoutput_dir: Path,
    model_name_expected: str,
    yaml_source: Optional[Path] = None,
) -> PyrateBetaModule:
    """
    Import the PyR@TE-generated `RGEs.py` from a PythonOutput directory.

    IMPORTANT:
    PyR@TE emits many different `RGEs.py` files with the *same module name*. We must load them
    via file location under a unique module name to avoid `sys.modules['RGEs']` collisions
    between SM and E8 models.
    """
    if not pythonoutput_dir.is_dir():
        raise FileNotFoundError(f"PythonOutput dir not found: {pythonoutput_dir}")
    exp = str(model_name_expected).strip()
    if not exp:
        raise ValueError("model_name_expected must be a non-empty string")
    if exp.lower() not in str(pythonoutput_dir).lower():
        raise RuntimeError(
            "Model mismatch (fail-fast): "
            f"model_name_expected={exp!r} not found in pythonoutput_dir={str(pythonoutput_dir)!r}"
        )

    rges_path = pythonoutput_dir / "RGEs.py"
    if not rges_path.is_file():
        raise FileNotFoundError(f"RGEs.py not found in: {pythonoutput_dir}")
    rges_sha = _sha256_file(rges_path)
    yaml_sha = _sha256_file(yaml_source) if (yaml_source is not None and yaml_source.exists()) else None

    uniq = hashlib.sha256(str(rges_path.resolve()).encode("utf-8")).hexdigest()[:12]
    module_name = f"tfpt_pyrate_rges_{kind}_{uniq}"
    spec = importlib.util.spec_from_file_location(module_name, str(rges_path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Failed to load module spec for {rges_path}")
    rges = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rges)  # type: ignore[union-attr]

    return PyrateBetaModule(
        kind=kind,
        pythonoutput_dir=pythonoutput_dir,
        pythonoutput_module_file=rges_path,
        pythonoutput_module_sha256=rges_sha,
        model_name_expected=exp,
        yaml_source=yaml_source if (yaml_source is not None and yaml_source.exists()) else None,
        yaml_source_sha256=yaml_sha,
        rges=rges,
    )


def _as_matrix(x: np.ndarray) -> np.matrix:
    # PyR@TE RGEs expect `np.matrix` to support `.H` (Hermitian transpose) via `adjoint(x) = x.H`.
    return np.matrix(x, dtype=complex)


def _pack_complex(*parts: np.ndarray) -> np.ndarray:
    return np.concatenate([p.reshape(-1).astype(complex) for p in parts], axis=0)


def _unpack_complex(y: np.ndarray, *, shapes: list[tuple[int, int] | None]) -> list[Any]:
    """
    Unpack a complex state vector into scalars/matrices given `shapes`:
    - shape=None means a scalar (complex)
    - shape=(m,n) means an (m,n) complex matrix
    """
    out: list[Any] = []
    idx = 0
    for sh in shapes:
        if sh is None:
            out.append(complex(y[idx]))
            idx += 1
        else:
            m, n = sh
            size = m * n
            out.append(np.array(y[idx : idx + size]).reshape(m, n).astype(complex))
            idx += size
    if idx != y.size:
        raise ValueError("state vector length mismatch")
    return out


def _solve_segment(
    *,
    rhs,
    mu_start_GeV: float,
    mu_end_GeV: float,
    y0: np.ndarray,
    rtol: float,
    atol: float,
    method: str,
) -> np.ndarray:
    t0 = float(np.log(mu_start_GeV))
    t1 = float(np.log(mu_end_GeV))
    sol = solve_ivp(rhs, t_span=(t0, t1), y0=y0, method=method, rtol=rtol, atol=atol, dense_output=False)
    if not sol.success:
        raise RuntimeError(f"2-loop PyR@TE-driven RGE integration failed: {sol.message}")
    return sol.y[:, -1].astype(complex)


def run_flavor_rge_2loop_thresholds(
    *,
    mu_start_GeV: float,
    mu_end_GeV: float,
    thresholds_GeV: dict[str, float],
    g_start: tuple[float, float, float],
    Yu_start: np.ndarray,
    Yd_start: np.ndarray,
    Ye_start: np.ndarray,
    lambda_start: float,
    yN_start: np.ndarray | None = None,
    # model beta sources
    beta_sm: PyrateBetaModule,
    beta_e8: PyrateBetaModule,
    # threshold toggles
    apply_sigma_threshold: bool = True,
    apply_g8_delta_b3: bool = True,
    delta_b3_g8: float = 2.0,
    apply_g8_delta_b3_2loop: bool = False,
    g8_delta_b3_2loop: float | None = None,
    # threshold matching (optional; default is continuous-by-assumption)
    apply_matching: bool = False,
    matching_loop_order: int = 1,
    matching_finite_delta_alpha: dict[str, dict[str, float]] | None = None,
    # optional gravity α^3 patch (runner-level; not part of PyR@TE model)
    apply_gravity_alpha3: bool = False,
    gravity_kappa: tuple[float, float, float] = (0.0, 0.0, 0.0),
    gravity_c3: float = float(1.0 / (8.0 * np.pi)),
    alpha_definition_for_u1: str = "alphaY",
    # loop order (PyR@TE provides beta_X(nLoop) coefficients; we can truncate to 1-loop for diagnostics)
    max_loop: int = 2,
    # numerical
    rtol: float = 1e-8,
    atol: float = 1e-10,
    method: str = "DOP853",
    # optional sampling of the trajectory (needed by diagnostic gates such as unification)
    return_trajectory: bool = False,
    trajectory_n_points: int = 260,
) -> dict[str, Any]:
    """
    2-loop RGE engine for SM Yukawas (and optionally yN) using PyR@TE-generated beta functions,
    with explicit threshold switching:

    - below MSigma: use SM betas (SM_TFPT_2Loop_v25)
    - above MSigma: use E8 sigma+yN betas (E8Cascade2LoopGravityV2 / E8Cascade2LoopGravity)
    - above MG8: patch the 1-loop QCD coefficient by Δb3 (paper v1.06 note), i.e. b3:-7→-7+Δb3
      (optional) and apply a 2-loop g3^5 patch consistent with a Weyl adjoint G8 assumption.

    Notes:
    - This integrates in t = ln(mu); PyR@TE’s `RGEs.py` provides `beta_X(nLoop)` coefficients and
      we apply kappa(n)=1/(4π)^(2n) explicitly.
    - Gauge couplings use the PyR@TE convention internally: g1 ≡ gY (SM hypercharge, not GUT-normalized).
    - We keep lHphi ≡ 0 in the E8 beta system (it stays 0 under its own beta function), which is consistent
      with using the minimal coupling subset for flavor running.
    """
    mu_start = float(mu_start_GeV)
    mu_end = float(mu_end_GeV)
    if mu_start <= 0 or mu_end <= 0:
        raise ValueError("mu_start_GeV and mu_end_GeV must be positive")
    if mu_start == mu_end:
        raise ValueError("mu_start_GeV and mu_end_GeV must differ")
    direction = "up" if mu_end > mu_start else "down"
    lo, hi = (mu_start, mu_end) if mu_start < mu_end else (mu_end, mu_start)

    MSigma = float(thresholds_GeV.get("MSigma", 1.0e3))
    MG8 = float(thresholds_GeV.get("MG8", 1.8e10))
    MNR1 = float(thresholds_GeV.get("MNR1", 1.0e14))
    MNR2 = float(thresholds_GeV.get("MNR2", 3.0e14))
    MNR3 = float(thresholds_GeV.get("MNR3", 8.0e14))
    if int(matching_loop_order) < 1:
        raise ValueError("matching_loop_order must be >= 1")
    max_loop_i = int(max_loop)
    if max_loop_i not in (1, 2):
        raise ValueError("max_loop must be 1 or 2")

    delta_b3_g8_2loop_value: float | None = None
    g8_delta_b3_2loop_mode = "disabled"
    if apply_g8_delta_b3_2loop and float(delta_b3_g8) != 0.0:
        if g8_delta_b3_2loop is not None:
            delta_b3_g8_2loop_value = float(g8_delta_b3_2loop)
            g8_delta_b3_2loop_mode = "explicit_delta_b3_2loop"
        else:
            # Weyl adjoint assumption (G8): Δb3(1-loop)=(2/3)C_A, Δb3(2-loop)=(16/3)C_A^2.
            ca_su3 = 3.0
            delta_b3_per_weyl_adj = (2.0 / 3.0) * ca_su3
            delta_b3_2loop_per_weyl_adj = (16.0 / 3.0) * (ca_su3**2)
            if delta_b3_per_weyl_adj != 0.0:
                n_adj = float(delta_b3_g8) / delta_b3_per_weyl_adj
                delta_b3_g8_2loop_value = n_adj * delta_b3_2loop_per_weyl_adj
                g8_delta_b3_2loop_mode = "weyl_adjoint_match"

    # Decide if we evolve yN at all
    evolve_yN = yN_start is not None
    if not evolve_yN:
        yN_start = np.zeros((3, 3), dtype=complex)

    # Declarative threshold rules (used for reporting + future matching upgrades).
    threshold_rules: list[ThresholdRule] = []
    thr_status = "matched_1loop_log_only_identity" if apply_matching and int(matching_loop_order) == 1 else "matched_identity"
    if apply_sigma_threshold and (lo < MSigma < hi):
        threshold_rules.append(
            ThresholdRule(
                threshold_id="MSigma",
                scale_GeV=float(MSigma),
                action="beta_source_switch(sm↔e8)",
                status=thr_status if apply_matching else "continuous_by_assumption",
                affected_parameters=["g", "Yu", "Yd", "Ye", "lambda"] + (["yN"] if evolve_yN else []),
                note=(
                    "Switch beta source SM→E8 at μ=MSigma (Sigma integrated in). "
                    + (
                        f"Matching enabled (loop_order={int(matching_loop_order)}; 1-loop is identity at μ=threshold)."
                        if apply_matching
                        else "No finite matching implemented yet."
                    )
                ),
            )
        )
    if apply_g8_delta_b3 and (lo < MG8 < hi):
        threshold_rules.append(
            ThresholdRule(
                threshold_id="MG8",
                scale_GeV=float(MG8),
                action="beta_patch(Δb3)",
                status=thr_status if apply_matching else "continuous_by_assumption",
                affected_parameters=["g3"],
                note=(
                    "Apply Δb3 patch above MG8 (paper v1.06 note). "
                    + (
                        "2-loop adjoint patch enabled (Weyl adjoint assumption). "
                        if (apply_g8_delta_b3_2loop and max_loop_i >= 2 and delta_b3_g8_2loop_value is not None)
                        else ""
                    )
                    + (
                        f"Matching enabled (loop_order={int(matching_loop_order)}; 1-loop is identity at μ=threshold)."
                        if apply_matching
                        else "No finite matching implemented yet."
                    )
                ),
            )
        )
    if evolve_yN:
        for i, m in enumerate([MNR1, MNR2, MNR3], start=1):
            if lo < m < hi:
                threshold_rules.append(
                    ThresholdRule(
                        threshold_id=f"MNR{i}",
                        scale_GeV=float(m),
                        action="activate_yN_column",
                        status=thr_status if apply_matching else "continuous_by_assumption",
                        affected_parameters=["yN"],
                        note=(
                            "Activate yN column above μ=MNRi via projector. "
                            + (
                                f"Matching enabled (loop_order={int(matching_loop_order)}; identity at μ=threshold)."
                                if apply_matching
                                else "No finite matching implemented in this engine (PMNS EFT does explicit tree-level κ matching)."
                            )
                        ),
                    )
                )

    # segment boundaries (always include endpoints)
    cuts: list[float] = [mu_start, mu_end]
    if apply_sigma_threshold and (lo < MSigma < hi):
        cuts.append(MSigma)
    if apply_g8_delta_b3 and (lo < MG8 < hi):
        cuts.append(MG8)
    if evolve_yN:
        for m in (MNR1, MNR2, MNR3):
            if lo < m < hi:
                cuts.append(m)
    cuts = sorted(set(cuts))
    if direction == "down":
        cuts = list(reversed(cuts))

    # Optional: sample a gauge-coupling trajectory on a log grid (plus exact segment boundaries).
    want_traj = bool(return_trajectory)
    traj_mu: list[float] = []
    traj_gY: list[float] = []
    traj_g2: list[float] = []
    traj_g3: list[float] = []
    mu_samples: list[float] = []
    if want_traj:
        n = int(trajectory_n_points)
        if n < 2:
            raise ValueError("trajectory_n_points must be >= 2")
        lo_s, hi_s = (mu_start, mu_end) if mu_start < mu_end else (mu_end, mu_start)
        base = np.logspace(np.log10(lo_s), np.log10(hi_s), n)
        base[0] = lo_s
        base[-1] = hi_s
        extra = np.array([float(x) for x in cuts], dtype=float)
        all_pts = np.concatenate([base, extra], axis=0)
        # Uniquify by ln(mu) key to avoid float duplicates, while preserving exact boundaries.
        uniq: dict[float, float] = {}
        for x in all_pts:
            xf = float(x)
            if xf <= 0:
                continue
            if xf < lo_s or xf > hi_s:
                continue
            key = float(round(float(np.log(xf)), 12))
            uniq[key] = xf
        mu_samples = sorted(uniq.values())
        if direction == "down":
            mu_samples = list(reversed(mu_samples))

    # state layout (scalars + matrices)
    shapes: list[tuple[int, int] | None] = [None, None, None, None, (3, 3), (3, 3), (3, 3)]
    if evolve_yN:
        shapes.append((3, 3))

    # initial state
    g1_0, g2_0, g3_0 = (float(g_start[0]), float(g_start[1]), float(g_start[2]))
    y_parts = [np.array([g1_0, g2_0, g3_0, float(lambda_start)], dtype=complex), Yu_start, Yd_start, Ye_start]
    if evolve_yN:
        y_parts.append(yN_start)
    y = _pack_complex(*y_parts)

    # helpers
    def rhs_for(kind: ModelKind, *, delta_b3_active: bool):
        if kind == "sm_tfpt_2loop_v25":
            rges = beta_sm.rges

            def rhs(t: float, yvec: np.ndarray) -> np.ndarray:
                if evolve_yN:
                    g1, g2, g3, lam, Yu, Yd, Ye, yN = _unpack_complex(
                        yvec, shapes=[None, None, None, None, (3, 3), (3, 3), (3, 3), (3, 3)]
                    )
                else:
                    g1, g2, g3, lam, Yu, Yd, Ye = _unpack_complex(
                        yvec, shapes=[None, None, None, None, (3, 3), (3, 3), (3, 3)]
                    )
                    yN = np.zeros((3, 3), dtype=complex)
                g1r, g2r, g3r, lamr = float(np.real(g1)), float(np.real(g2)), float(np.real(g3)), float(np.real(lam))
                Yu_m, Yd_m, Ye_m = _as_matrix(Yu), _as_matrix(Yd), _as_matrix(Ye)

                dg1 = _kappa(1) * rges.beta_g1(1, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m)
                dg2 = _kappa(1) * rges.beta_g2(1, g2r, g1r, g3r, Yu_m, Yd_m, Ye_m)
                dg3_1 = rges.beta_g3(1, g3r, g1r, g2r, Yu_m, Yd_m)
                if delta_b3_active:
                    # Patch only the 1-loop coefficient as per paper v1.06 note: Δb3 = +2 above MG8.
                    dg3_1 = dg3_1 + float(delta_b3_g8) * (g3r**3)
                dg3 = _kappa(1) * dg3_1
                if max_loop_i >= 2:
                    dg1 = dg1 + _kappa(2) * rges.beta_g1(2, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m)
                    dg2 = dg2 + _kappa(2) * rges.beta_g2(2, g2r, g1r, g3r, Yu_m, Yd_m, Ye_m)
                    dg3 = dg3 + _kappa(2) * rges.beta_g3(2, g3r, g1r, g2r, Yu_m, Yd_m)
                    if delta_b3_active and apply_g8_delta_b3_2loop and delta_b3_g8_2loop_value is not None:
                        dg3 = dg3 + _kappa(2) * float(delta_b3_g8_2loop_value) * (g3r**5)

                # Optional TFPT gravity α^3 correction (runner-level):
                #   Δ(dα_i)/d ln μ = κ_i * c3 * α_i^3  ⇒  Δ(dg_i)/d ln μ = κ_i * c3 * g_i^5 / (32 π^2)
                if apply_gravity_alpha3:
                    if str(alpha_definition_for_u1).strip().lower() not in ("alphay", "alpha_y", "alpha-y"):
                        raise ValueError(f"Unsupported alpha_definition_for_u1={alpha_definition_for_u1!r}; expected 'alphaY'")
                    k1, k2, k3 = [float(x) for x in gravity_kappa]
                    c3v = float(gravity_c3)
                    if c3v != 0.0:
                        dg1 = dg1 + k1 * c3v * (g1r**5) / (32.0 * (np.pi**2))
                        dg2 = dg2 + k2 * c3v * (g2r**5) / (32.0 * (np.pi**2))
                        dg3 = dg3 + k3 * c3v * (g3r**5) / (32.0 * (np.pi**2))

                dYu = _kappa(1) * rges.beta_Yu(1, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, lamr)
                dYd = _kappa(1) * rges.beta_Yd(1, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, lamr)
                dYe = _kappa(1) * rges.beta_Ye(1, g1r, g2r, Yu_m, Yd_m, Ye_m, g3r, lamr)
                dlam = _kappa(1) * rges.beta_lambda_(1, g1r, g2r, Yu_m, Yd_m, Ye_m, lamr, g3r)
                if max_loop_i >= 2:
                    dYu = dYu + _kappa(2) * rges.beta_Yu(2, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, lamr)
                    dYd = dYd + _kappa(2) * rges.beta_Yd(2, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, lamr)
                    dYe = dYe + _kappa(2) * rges.beta_Ye(2, g1r, g2r, Yu_m, Yd_m, Ye_m, g3r, lamr)
                    dlam = dlam + _kappa(2) * rges.beta_lambda_(2, g1r, g2r, Yu_m, Yd_m, Ye_m, lamr, g3r)

                parts = [
                    np.array([dg1, dg2, dg3, dlam], dtype=complex),
                    np.array(dYu, dtype=complex),
                    np.array(dYd, dtype=complex),
                    np.array(dYe, dtype=complex),
                ]
                if evolve_yN:
                    # yN is not part of the SM beta system; keep it frozen below MSigma.
                    parts.append(np.zeros((3, 3), dtype=complex))
                return _pack_complex(*parts)

            return rhs

        # E8 sigma + yN model
        rges = beta_e8.rges

        def rhs(t: float, yvec: np.ndarray) -> np.ndarray:
            if evolve_yN:
                g1, g2, g3, lam, Yu, Yd, Ye, yN = _unpack_complex(
                    yvec, shapes=[None, None, None, None, (3, 3), (3, 3), (3, 3), (3, 3)]
                )
            else:
                g1, g2, g3, lam, Yu, Yd, Ye = _unpack_complex(yvec, shapes=[None, None, None, None, (3, 3), (3, 3), (3, 3)])
                yN = np.zeros((3, 3), dtype=complex)

            g1r, g2r, g3r, lamr = float(np.real(g1)), float(np.real(g2)), float(np.real(g3)), float(np.real(lam))
            Yu_m, Yd_m, Ye_m = _as_matrix(Yu), _as_matrix(Yd), _as_matrix(Ye)

            # RH neutrino thresholds: activate yN columns sequentially above (MNR1,MNR2,MNR3).
            # We model this as a right-multiplication projector on yN (columns correspond to NR1..3).
            if evolve_yN:
                mu = float(np.exp(t))
                P = np.diag(
                    [
                        1.0 if mu >= MNR1 else 0.0,
                        1.0 if mu >= MNR2 else 0.0,
                        1.0 if mu >= MNR3 else 0.0,
                    ]
                ).astype(complex)
                yN_eff = (yN @ P).astype(complex)
            else:
                P = np.eye(3, dtype=complex)
                yN_eff = yN
            yN_m = _as_matrix(yN_eff)

            dg1 = _kappa(1) * rges.beta_g1(1, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, yN_m)
            dg2 = _kappa(1) * rges.beta_g2(1, g2r, g1r, g3r, Yu_m, Yd_m, Ye_m, yN_m)
            dg3_1 = rges.beta_g3(1, g3r, g1r, g2r, Yu_m, Yd_m)
            if delta_b3_active:
                dg3_1 = dg3_1 + float(delta_b3_g8) * (g3r**3)
            dg3 = _kappa(1) * dg3_1
            if max_loop_i >= 2:
                dg1 = dg1 + _kappa(2) * rges.beta_g1(2, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, yN_m)
                dg2 = dg2 + _kappa(2) * rges.beta_g2(2, g2r, g1r, g3r, Yu_m, Yd_m, Ye_m, yN_m)
                dg3 = dg3 + _kappa(2) * rges.beta_g3(2, g3r, g1r, g2r, Yu_m, Yd_m)
                if delta_b3_active and apply_g8_delta_b3_2loop and delta_b3_g8_2loop_value is not None:
                    dg3 = dg3 + _kappa(2) * float(delta_b3_g8_2loop_value) * (g3r**5)

            if apply_gravity_alpha3:
                if str(alpha_definition_for_u1).strip().lower() not in ("alphay", "alpha_y", "alpha-y"):
                    raise ValueError(f"Unsupported alpha_definition_for_u1={alpha_definition_for_u1!r}; expected 'alphaY'")
                k1, k2, k3 = [float(x) for x in gravity_kappa]
                c3v = float(gravity_c3)
                if c3v != 0.0:
                    dg1 = dg1 + k1 * c3v * (g1r**5) / (32.0 * (np.pi**2))
                    dg2 = dg2 + k2 * c3v * (g2r**5) / (32.0 * (np.pi**2))
                    dg3 = dg3 + k3 * c3v * (g3r**5) / (32.0 * (np.pi**2))

            # keep lHphi=0 (it stays 0 if initialized 0); beta functions accept it as a scalar
            lHphi = 0.0
            dYu = _kappa(1) * rges.beta_Yu(1, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, yN_m, lamr, lHphi)
            dYd = _kappa(1) * rges.beta_Yd(1, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, yN_m, lamr, lHphi)
            dYe = _kappa(1) * rges.beta_Ye(1, g1r, g2r, Yu_m, Yd_m, Ye_m, yN_m, g3r, lamr, lHphi)
            dlam = _kappa(1) * rges.beta_lambda_(1, g1r, g2r, Yu_m, Yd_m, Ye_m, yN_m, lamr, lHphi, g3r)
            if max_loop_i >= 2:
                dYu = dYu + _kappa(2) * rges.beta_Yu(2, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, yN_m, lamr, lHphi)
                dYd = dYd + _kappa(2) * rges.beta_Yd(2, g1r, g2r, g3r, Yu_m, Yd_m, Ye_m, yN_m, lamr, lHphi)
                dYe = dYe + _kappa(2) * rges.beta_Ye(2, g1r, g2r, Yu_m, Yd_m, Ye_m, yN_m, g3r, lamr, lHphi)
                dlam = dlam + _kappa(2) * rges.beta_lambda_(2, g1r, g2r, Yu_m, Yd_m, Ye_m, yN_m, lamr, lHphi, g3r)

            parts = [np.array([dg1, dg2, dg3, dlam], dtype=complex), np.array(dYu, dtype=complex), np.array(dYd, dtype=complex), np.array(dYe, dtype=complex)]
            if evolve_yN:
                dyN = _kappa(1) * rges.beta_yN(1, g1r, g2r, Yu_m, Yd_m, Ye_m, yN_m, g3r, lamr, lHphi)
                if max_loop_i >= 2:
                    dyN = dyN + _kappa(2) * rges.beta_yN(2, g1r, g2r, Yu_m, Yd_m, Ye_m, yN_m, g3r, lamr, lHphi)
                # Keep inactive columns frozen below their thresholds.
                parts.append(np.array(dyN, dtype=complex) @ P)

            return _pack_complex(*parts)

        return rhs

    # Map threshold scale → rule (exact float match; thresholds originate from the same dict we used for cuts)
    rules_by_mu: dict[float, ThresholdRule] = {float(r.scale_GeV): r for r in threshold_rules}

    def active_fields_at(mu_GeV: float) -> list[str]:
        """
        Minimal “active fields” bookkeeping for audit trails.

        This is intentionally conservative and only encodes what the runner actually toggles:
        - SM↔E8 beta source switch at MSigma (Sigma integrated in)
        - Δb3 patch above MG8 (G8 bridge note)
        - sequential activation of yN columns above MNRi (proxy for integrating in N_Ri)
        """
        mu = float(mu_GeV)
        fields: list[str] = ["SM"]
        if apply_sigma_threshold and mu >= MSigma:
            fields.append("Sigma")
        if apply_g8_delta_b3 and mu >= MG8:
            fields.append("G8")
        if evolve_yN:
            if mu >= MNR1:
                fields.append("NR1")
            if mu >= MNR2:
                fields.append("NR2")
            if mu >= MNR3:
                fields.append("NR3")
        return fields

    # integrate piecewise
    events: list[dict[str, Any]] = []
    for a, bnd in zip(cuts[:-1], cuts[1:]):
        mu_probe = float(np.sqrt(float(a) * float(bnd)))
        # Apply matching exactly at the segment start (after any threshold transition).
        rule_at_start = rules_by_mu.get(float(a))
        threshold_match: dict[str, Any] | None = None
        if apply_matching and rule_at_start is not None and float(a) != float(mu_start):
            # Unpack, match (gauge, Yukawa blocks, quartic) explicitly, repack.
            if evolve_yN:
                g1, g2, g3, lam, Yu, Yd, Ye, yN = _unpack_complex(y, shapes=[None, None, None, None, (3, 3), (3, 3), (3, 3), (3, 3)])
            else:
                g1, g2, g3, lam, Yu, Yd, Ye = _unpack_complex(y, shapes=[None, None, None, None, (3, 3), (3, 3), (3, 3)])
                yN = np.zeros((3, 3), dtype=complex)
            g1r, g2r, g3r = float(np.real(g1)), float(np.real(g2)), float(np.real(g3))
            finite = (matching_finite_delta_alpha or {}).get(str(rule_at_start.threshold_id), None)
            eps = 1e-12
            mu_below = float(a) * (1.0 - eps)
            mu_above = float(a) * (1.0 + eps)
            fields_below = active_fields_at(mu_below)
            fields_above = active_fields_at(mu_above)
            if direction == "up":
                fields_before, fields_after = fields_below, fields_above
            else:
                fields_before, fields_after = fields_above, fields_below

            g_above, out_g = match_gauge(
                threshold_id=str(rule_at_start.threshold_id),
                mu_thr_GeV=float(a),
                direction=direction,
                couplings_below={"gY": g1r, "g2": g2r, "g3": g3r},
                scheme="MSbar",
                loop_order=int(matching_loop_order),
                active_fields_before=fields_before,
                active_fields_after=fields_after,
                finite_delta_alpha=finite,
            )
            # Yukawa and quartic matching are currently identity (but explicit + audited).
            yuk_in: dict[str, object] = {"Yu": np.array(Yu, dtype=complex), "Yd": np.array(Yd, dtype=complex), "Ye": np.array(Ye, dtype=complex)}
            if evolve_yN:
                yuk_in["yN"] = np.array(yN, dtype=complex)
            yuk_out, out_y = match_yukawa(
                threshold_id=str(rule_at_start.threshold_id),
                mu_thr_GeV=float(a),
                direction=direction,
                yukawas_below=yuk_in,
                scheme="MSbar",
                loop_order=int(matching_loop_order),
                active_fields_before=fields_before,
                active_fields_after=fields_after,
            )
            q_out, out_q = match_quartic(
                threshold_id=str(rule_at_start.threshold_id),
                mu_thr_GeV=float(a),
                direction=direction,
                quartics_below={"lambda": float(np.real(lam))},
                scheme="MSbar",
                loop_order=int(matching_loop_order),
                active_fields_before=fields_before,
                active_fields_after=fields_after,
            )

            g1 = complex(float(g_above.get("gY", g1r)))
            g2 = complex(float(g_above.get("g2", g2r)))
            g3 = complex(float(g_above.get("g3", g3r)))
            lam = complex(float(q_out.get("lambda", float(np.real(lam)))))
            Yu = np.array(yuk_out.get("Yu", Yu), dtype=complex)
            Yd = np.array(yuk_out.get("Yd", Yd), dtype=complex)
            Ye = np.array(yuk_out.get("Ye", Ye), dtype=complex)
            if evolve_yN:
                yN = np.array(yuk_out.get("yN", yN), dtype=complex)

            parts = [np.array([g1, g2, g3, lam], dtype=complex), Yu, Yd, Ye]
            if evolve_yN:
                parts.append(yN)
            y = _pack_complex(*parts)

            threshold_match = {
                "threshold_id": str(rule_at_start.threshold_id),
                "matching_active": True,
                "scheme": "MSbar",
                "loop_order": int(matching_loop_order),
                "gauge": {"status": out_g.status, "note": out_g.note, "deltas": out_g.deltas, "details": out_g.details},
                "yukawa": {"status": out_y.status, "note": out_y.note, "deltas": out_y.deltas, "details": out_y.details},
                "quartic": {"status": out_q.status, "note": out_q.note, "deltas": out_q.deltas, "details": out_q.details},
            }

        # choose model per segment
        use_sm = (mu_probe < MSigma) if apply_sigma_threshold else False
        kind: ModelKind = "sm_tfpt_2loop_v25" if use_sm else "e8_sigma_yN_2loop"
        delta_b3_active = apply_g8_delta_b3 and (mu_probe >= MG8)
        rhs = rhs_for(kind, delta_b3_active=delta_b3_active)
        if want_traj:
            # Evaluate the segment on the global mu grid restricted to [a,bnd].
            seg_lo, seg_hi = (float(a), float(bnd)) if float(a) < float(bnd) else (float(bnd), float(a))
            eps = 1e-14
            mu_eval = [m for m in mu_samples if (seg_lo * (1.0 - eps) <= float(m) <= seg_hi * (1.0 + eps))]
            if not mu_eval:
                mu_eval = [float(a), float(bnd)]
            mu_eval = sorted(set([float(x) for x in mu_eval]))
            if direction == "down":
                mu_eval = list(reversed(mu_eval))

            t0 = float(np.log(float(a)))
            t1 = float(np.log(float(bnd)))
            t_eval = np.array([float(np.log(float(m))) for m in mu_eval], dtype=float)
            sol = solve_ivp(rhs, t_span=(t0, t1), y0=y, t_eval=t_eval, method=method, rtol=rtol, atol=atol, dense_output=False)
            if not sol.success:
                raise RuntimeError(f"2-loop PyR@TE-driven RGE integration failed: {sol.message}")

            # Record gauge couplings (SM convention: g1 ≡ gY) at sampled points.
            for j, mu_j in enumerate(mu_eval):
                gY_j = float(np.real(sol.y[0, j]))
                g2_j = float(np.real(sol.y[1, j]))
                g3_j = float(np.real(sol.y[2, j]))
                if traj_mu and abs(float(mu_j) - float(traj_mu[-1])) / max(float(mu_j), float(traj_mu[-1]), 1.0) < 1e-13:
                    # Replace the last point (segment boundary) with the post-matching value.
                    traj_mu[-1] = float(mu_j)
                    traj_gY[-1] = float(gY_j)
                    traj_g2[-1] = float(g2_j)
                    traj_g3[-1] = float(g3_j)
                else:
                    traj_mu.append(float(mu_j))
                    traj_gY.append(float(gY_j))
                    traj_g2.append(float(g2_j))
                    traj_g3.append(float(g3_j))

            y = sol.y[:, -1].astype(complex)
        else:
            y = _solve_segment(rhs=rhs, mu_start_GeV=float(a), mu_end_GeV=float(bnd), y0=y, rtol=rtol, atol=atol, method=method)
        event: dict[str, Any] = {
            "mu_start_GeV": float(a),
            "mu_end_GeV": float(bnd),
            # Public/portable aliases (preferred going forward).
            "mu_start": float(a),
            "mu_end": float(bnd),
            "model": kind,
            "delta_b3_active": bool(delta_b3_active),
        }
        patches_active: list[str] = []
        if delta_b3_active:
            patches_active.append("delta_b3_g8")
        if delta_b3_active and apply_g8_delta_b3_2loop and max_loop_i >= 2 and delta_b3_g8_2loop_value is not None:
            patches_active.append("delta_b3_g8_2loop")
        if apply_gravity_alpha3:
            patches_active.append("gravity_alpha3")
        event["patches_active"] = patches_active
        # Public/portable alias (preferred going forward).
        event["patches"] = list(patches_active)
        # Always present for schema stability (None if no threshold at start).
        event["threshold_match"] = None
        if evolve_yN:
            # With cuts including MNR thresholds, yN column activation is constant within each segment.
            event["yN_active_cols"] = [bool(mu_probe >= MNR1), bool(mu_probe >= MNR2), bool(mu_probe >= MNR3)]
        rule = rules_by_mu.get(float(a))
        if rule is not None and float(a) != float(mu_start):
            event["threshold_transition_at_start"] = rule
            event["threshold_actions_at_start"] = [rule]
            if apply_matching and threshold_match is not None:
                event["threshold_match"] = threshold_match
            else:
                # Explicitly mark that this segment start is “by assumption” (not finite-matched).
                event["threshold_match"] = {
                    "threshold_id": str(getattr(rule, "threshold_id", "?")),
                    "matching_active": False,
                    "status": str(getattr(rule, "status", "continuous_by_assumption")),
                    "note": str(getattr(rule, "note", "")),
                }
        else:
            event["threshold_actions_at_start"] = []
        events.append(event)

    # unpack end state
    if evolve_yN:
        g1, g2, g3, lam, Yu, Yd, Ye, yN = _unpack_complex(y, shapes=[None, None, None, None, (3, 3), (3, 3), (3, 3), (3, 3)])
    else:
        g1, g2, g3, lam, Yu, Yd, Ye = _unpack_complex(y, shapes=[None, None, None, None, (3, 3), (3, 3), (3, 3)])
        yN = np.zeros((3, 3), dtype=complex)

    gY_end = float(np.real(g1))
    g2_end = float(np.real(g2))
    g3_end = float(np.real(g3))
    g1_gut_end = float(g1_gut_over_gY() * gY_end)
    blocked_thresholds = [str(r.threshold_id) for r in threshold_rules] if (threshold_rules and not apply_matching) else []
    threshold_matching_ok = (len(threshold_rules) == 0) or bool(apply_matching)

    return {
        "direction": direction,
        "mu_start_GeV": mu_start,
        "mu_end_GeV": mu_end,
        "thresholds_GeV": {"MSigma": MSigma, "MG8": MG8, "MNR1": MNR1, "MNR2": MNR2, "MNR3": MNR3},
        "u1_convention": {
            "internal_g1_is": "gY (SM hypercharge, PyR@TE convention)",
            "output_alpha1_gut_relation": "alpha1_GUT = (5/3) * alphaY",
            "g1_gut_over_gY": float(g1_gut_over_gY()),
        },
        "matching": {
            "enabled": bool(apply_matching),
            "loop_order": int(matching_loop_order),
            "finite_delta_alpha_by_threshold": matching_finite_delta_alpha,
        },
        "publication_grade": {
            "threshold_matching_ok": bool(threshold_matching_ok),
            "blocked_thresholds": blocked_thresholds,
            "note": "Publication-grade threshold bookkeeping requires matching_active=True at each threshold boundary; "
            "if matching is disabled, segments are labeled 'continuous_by_assumption'.",
        },
        "gravity_alpha3_patch": {
            "enabled": bool(apply_gravity_alpha3),
            "kappa_vector": [float(x) for x in gravity_kappa],
            "c3": float(gravity_c3),
            "alpha_definition_for_U1": str(alpha_definition_for_u1),
        },
        "g8_patch": {
            "delta_b3_g8": float(delta_b3_g8),
            "apply_g8_delta_b3_2loop": bool(apply_g8_delta_b3_2loop),
            "delta_b3_g8_2loop": float(delta_b3_g8_2loop_value) if delta_b3_g8_2loop_value is not None else None,
            "delta_b3_2loop_mode": str(g8_delta_b3_2loop_mode),
            "note": "2-loop Δb3 patch assumes Weyl adjoint if delta_b3_g8_2loop is not explicitly provided.",
        },
        "beta_sources": {
            "sm": {
                "kind": beta_sm.kind,
                "model_name_expected": beta_sm.model_name_expected,
                "pythonoutput_dir": _relpath(beta_sm.pythonoutput_dir),
                "pythonoutput_module_file": _relpath(beta_sm.pythonoutput_module_file),
                "pythonoutput_module_sha256": beta_sm.pythonoutput_module_sha256,
                "yaml_source": _relpath(beta_sm.yaml_source) if beta_sm.yaml_source is not None else None,
                "yaml_source_sha256": beta_sm.yaml_source_sha256,
            },
            "e8": {
                "kind": beta_e8.kind,
                "model_name_expected": beta_e8.model_name_expected,
                "pythonoutput_dir": _relpath(beta_e8.pythonoutput_dir),
                "pythonoutput_module_file": _relpath(beta_e8.pythonoutput_module_file),
                "pythonoutput_module_sha256": beta_e8.pythonoutput_module_sha256,
                "yaml_source": _relpath(beta_e8.yaml_source) if beta_e8.yaml_source is not None else None,
                "yaml_source_sha256": beta_e8.yaml_source_sha256,
            },
        },
        "threshold_rules": threshold_rules,
        "segments": events,
        **(
            {
                "trajectory": {
                    "mu_GeV": [float(x) for x in traj_mu],
                    "log10_mu": [float(np.log10(float(x))) for x in traj_mu],
                    "g_sm": {"gY": [float(x) for x in traj_gY], "g2": [float(x) for x in traj_g2], "g3": [float(x) for x in traj_g3]},
                    "note": "Sampled gauge-coupling trajectory after applying any explicit matching at segment starts; g1 is reported as gY (SM hypercharge, PyR@TE convention).",
                }
            }
            if want_traj
            else {}
        ),
        # Back-compat key: g1 here means gY (SM hypercharge normalization).
        "g_end": {"g1": gY_end, "g2": g2_end, "g3": g3_end},
        "g_end_sm": {"gY": gY_end, "g2": g2_end, "g3": g3_end},
        "g_end_gut": {"g1_gut": g1_gut_end, "g2": g2_end, "g3": g3_end},
        "lambda_end": float(np.real(lam)),
        "Yu_end": Yu,
        "Yd_end": Yd,
        "Ye_end": Ye,
        "yN_end": yN,
    }

