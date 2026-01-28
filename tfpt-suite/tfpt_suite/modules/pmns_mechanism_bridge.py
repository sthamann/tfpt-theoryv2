from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from itertools import permutations
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.flavor_textures import (
    coefficients_v107sm,
    left_unitary_from_yukawa,
    scale_y_star_to_match_sigma_max,
    theta_of_delta,
    yukawa_texture_matrix,
)
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.pyrate_boundary_runner import sm_boundary_conditions_at_mt


def _workspace_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _relpath(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(_workspace_root()))
    except Exception:
        return str(path)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


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
    e_minus = np.exp(-1j * float(delta))
    e_plus = np.exp(1j * float(delta))
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
    """
    s13 = float(np.abs(U[0, 2]))
    s13 = min(max(s13, 0.0), 1.0)
    c13 = float(np.sqrt(max(0.0, 1.0 - s13 * s13)))

    s12 = float(np.abs(U[0, 1]) / c13) if c13 > 0 else 0.0
    s12 = min(max(s12, 0.0), 1.0)
    c12 = float(np.sqrt(max(0.0, 1.0 - s12 * s12)))

    s23 = float(np.abs(U[1, 2]) / c13) if c13 > 0 else 0.0
    s23 = min(max(s23, 0.0), 1.0)
    c23 = float(np.sqrt(max(0.0, 1.0 - s23 * s23)))

    J = float(np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])))
    denom_sin = c12 * c23 * (c13**2) * s12 * s23 * s13
    sin_delta = J / denom_sin if denom_sin != 0 else 0.0
    sin_delta = float(np.clip(sin_delta, -1.0, 1.0))

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


def _pmns_from_ye_and_kappa(*, Ye: np.ndarray, kappa: np.ndarray, v_ev_GeV: float, tol: float = 1e-12) -> tuple[np.ndarray, np.ndarray]:
    """
    Build PMNS from charged lepton Yukawa Ye and Weinberg operator κ (Majorana).

    Convention (matches pmns_full_pipeline):
      mν = v^2 κ, and we use a Takagi proxy via SVD for the symmetric complex matrix mν.
    """
    if Ye.shape != (3, 3) or kappa.shape != (3, 3):
        raise ValueError("Ye and kappa must be 3x3")
    if not np.allclose(kappa, kappa.T, atol=tol):
        raise ValueError("kappa must be symmetric (Majorana)")

    Ue = left_unitary_from_yukawa(Ye)
    mnu = (float(v_ev_GeV) ** 2) * kappa

    # Takagi proxy via SVD (for symmetric matrices, V ≈ U* and mnu ≈ U Σ U^T).
    #
    # IMPORTANT: np.linalg.svd returns singular values in descending order; for PDG angle extraction
    # we want columns ordered by increasing neutrino masses (normal ordering labels).
    Uu, s, _vh = np.linalg.svd(mnu, full_matrices=True)
    idx = np.argsort(s.real)  # ascending masses
    s = s[idx]
    U = Uu[:, idx].copy()
    for k in range(3):
        col = U[:, k]
        i = int(np.argmax(np.abs(col)))
        a = col[i]
        if abs(a) == 0:
            continue
        phase = a / abs(a)
        U[:, k] = col / phase
        if U[i, k].real < 0:
            U[:, k] *= -1

    U_pmns = (Ue.conj().T @ U).astype(complex)
    mnu_eV = np.sort((s.real * 1.0e9).astype(float))  # singular values in eV
    return U_pmns, mnu_eV


def _beta_kappa_1loop(
    *,
    g2: float,
    lambda_h: float,
    Yu: np.ndarray,
    Yd: np.ndarray,
    Ye: np.ndarray,
    yN: np.ndarray,
    kappa: np.ndarray,
) -> np.ndarray:
    """
    1-loop EFT beta for the SM Weinberg operator κ (dimension-5) in MSbar (minimal form).

    16π^2 dκ/dt = [ -3 g2^2 + 6 λ + 2 Tr(3 Yu†Yu + 3 Yd†Yd + Ye†Ye + yN†yN) ] κ
                 - 3/2 [ (Ye†Ye)^T κ + κ (Ye†Ye) ].
    """
    Hu = Yu.conj().T @ Yu
    Hd = Yd.conj().T @ Yd
    He = Ye.conj().T @ Ye
    Hn = yN.conj().T @ yN
    T = float(np.real(3.0 * np.trace(Hu) + 3.0 * np.trace(Hd) + np.trace(He) + np.trace(Hn)))

    alpha = (-3.0 * (g2**2)) + (6.0 * lambda_h) + (2.0 * T)
    term = alpha * kappa - 1.5 * ((He.T @ kappa) + (kappa @ He))
    term = 0.5 * (term + term.T)  # enforce symmetry
    return (1.0 / (16.0 * np.pi**2)) * term


def _mat_real_imag(M: np.ndarray) -> dict[str, list[list[float]]]:
    M = np.array(M, dtype=complex)
    return {"real": M.real.astype(float).tolist(), "imag": M.imag.astype(float).tolist()}


class PmnsMechanismBridgeModule(TfptModule):
    module_id = "pmns_mechanism_bridge"
    title = "PMNS mechanism bridge (Z3-breaking angles → κ(mt) → yN+MNR reconstruction)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT varphi0, gamma(0)=5/6, epsilon=varphi0/6 (Z3-breaking scan constraints)",
                "SM boundary at mt from tfpt_suite/data/sm_inputs_mz.json (v(mt), ytau(mt), g2(mt), lambda(mt))",
                "charged-lepton Yukawa texture Ye(mt) from tfpt_suite/data/flavor_texture_v24.json (for Ue basis)",
                "heavy-neutrino thresholds MNR1..3 from tfpt_suite/data/rge_thresholds_v25.json",
            ],
            outputs=[
                "selected Z3-breaking variant angles (θ12, θ13, θ23, δCP)",
                "reconstructed κ(mt) consistent with those angles (and a minimal neutrino-mass spectrum choice)",
                "a consistent (yN, MNR) factorization with κ = yN M^{-1} yN^T and explicit κ decoupling above MNR3",
                "EFT κ running (mt→0.99 MNR1) and angle stability diagnostic",
            ],
            formulas=[
                "TM1 baseline: sin^2(theta13)=varphi0*exp(-5/6), sin^2(theta12)=(1/3)(1-2 sin^2 theta13), theta23=45°, δ=90°",
                "Z3-breaking scale: ε=varphi0/6 (fixed; no new continuous parameter); discrete operator variants act on U_PMNS",
                "Majorana EFT: mν=v^2 κ; (Takagi) mν ≈ U diag(m_i) U^T",
                "Tree-level seesaw matching: κ = Σ_i (y_i y_i^T)/M_i, with diagonal M_i at thresholds",
            ],
            validation=[
                "Reconstructed κ reproduces the selected angles at mt (numerical tolerance).",
                "κ is symmetric and decouples to ~0 above MNR3 under tree-level matching inversion.",
            ],
            determinism="Deterministic (no scan; fixed selection rule for the Z3-breaking variant).",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        flavor_cfg_path = Path(__file__).resolve().parent.parent / "data" / "flavor_texture_v24.json"
        thresholds_path = Path(__file__).resolve().parent.parent / "data" / "rge_thresholds_v25.json"
        lepton_masses_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"

        # --- Z3-breaking scan constraints (mirror pmns_z3_breaking; choose a "best" discrete variant) ---
        gamma0 = 5.0 / 6.0
        sin2_theta13 = float(c.varphi0) * float(np.exp(-gamma0))
        theta13 = float(np.arcsin(np.sqrt(max(0.0, min(1.0, sin2_theta13)))))
        sin2_theta12 = (1.0 / 3.0) * (1.0 - 2.0 * sin2_theta13)
        theta12 = float(np.arcsin(np.sqrt(max(0.0, min(1.0, sin2_theta12)))))
        theta23_0 = float(np.deg2rad(45.0))
        delta0 = float(np.deg2rad(90.0))
        U0 = _pmns_pdg(theta12, theta13, theta23_0, delta0)
        base_angles = _extract_pdg_angles(U0)

        eps = float(c.varphi0 / 6.0)
        phase_Z3 = np.exp(2j * np.pi / 3)
        variants = [
            ("L_R23(+eps)", lambda U: _rot23(+eps) @ U),
            ("L_R23(-eps)", lambda U: _rot23(-eps) @ U),
            ("R_R23(+eps)", lambda U: U @ _rot23(+eps)),
            ("R_R23(-eps)", lambda U: U @ _rot23(-eps)),
            ("L_R23(+eps)*Z3phase", lambda U: (_rot23(+eps) @ U) @ np.diag([1.0, phase_Z3, 1.0])),
            ("L_R23(-eps)*Z3phase", lambda U: (_rot23(-eps) @ U) @ np.diag([1.0, phase_Z3, 1.0])),
        ]

        @dataclass(frozen=True)
        class _Variant:
            name: str
            angles: MixingAngles
            delta_theta23_deg: float
            delta_delta_deg: float

        var_rows: list[_Variant] = []
        for name, op in variants:
            Uv = op(U0)
            ang = _extract_pdg_angles(Uv)
            dth = ang.theta23_deg - base_angles.theta23_deg
            dd = ((ang.delta_cp_deg - base_angles.delta_cp_deg + 540.0) % 360.0) - 180.0
            var_rows.append(_Variant(name=name, angles=ang, delta_theta23_deg=dth, delta_delta_deg=dd))

        # Selection rule ("best variant"):
        #   1) minimize |Δδ| (keep δCP ≈ 90°),
        #   2) prefer Δθ23 > 0 (second octant direction),
        #   3) maximize Δθ23 among ties.
        def score(v: _Variant) -> tuple[float, int, float]:
            return (abs(v.delta_delta_deg), 0 if v.delta_theta23_deg > 0 else 1, -v.delta_theta23_deg)

        best = min(var_rows, key=score)

        # --- Inputs for Ye(mt) and v(mt) from the central SM boundary layer ---
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        sm_bc_mt = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_raw)
        v_ev_GeV = float(sm_bc_mt.route_2loop["v_ev_GeV"])
        ytau_mt_target = float(sm_bc_mt.route_2loop["ytau_mt"])
        yt_mt = float(sm_bc_mt.route_2loop["yt_mt"])
        yb_mt = float(sm_bc_mt.route_2loop["yb_mt"])
        g2_mt = float(sm_bc_mt.route_2loop["g2_mt"])
        lambda_mt = float(sm_bc_mt.route_2loop["lambda_mt"])

        # --- Charged-lepton basis (Ye texture) ---
        flavor_cfg = json.loads(flavor_cfg_path.read_text(encoding="utf-8"))
        tex = flavor_cfg.get("yukawa_texture", {})
        delta_source = str(tex.get("delta_source", "tau_mu")).strip()
        phase_mode = str(tex.get("phase_mode", "2pi_delta")).strip()
        if phase_mode not in ("2pi_delta", "delta_rad", "koide_pi_over_12"):
            raise ValueError(f"Unsupported phase_mode in {flavor_cfg_path}: {phase_mode}")

        lep = json.loads(lepton_masses_path.read_text(encoding="utf-8"))
        mtau = float(lep["masses"]["tau"]["mean"])
        mmu = float(lep["masses"]["muon"]["mean"])
        R = float(np.sqrt(mtau / mmu))
        delta_M = float((R - 1.0) / (R + 1.0))
        delta_star = float(c.delta_star)
        if delta_source == "tau_mu":
            delta_used = delta_M
        elif delta_source in ("delta_star", "delta_*", "delta_star_varphi0"):
            delta_used = delta_star
        else:
            raise ValueError(f"Unsupported delta_source in {flavor_cfg_path}: {delta_source}")
        theta = float(theta_of_delta(delta_used, phase_mode=phase_mode))
        zeta = complex(np.cos(theta), np.sin(theta))

        coeff = coefficients_v107sm()
        Ye_base = yukawa_texture_matrix(
            delta=float(delta_used),
            varphi0=float(c.varphi0),
            c3=float(c.c3),
            a_y=float(coeff.a_e),
            b=float(coeff.b),
            y_star=1.0,
            phase_mode=phase_mode,
        )
        ystar_e = scale_y_star_to_match_sigma_max(target_y3=ytau_mt_target, base=Ye_base)
        Ye_mt = (ystar_e * Ye_base).astype(complex)

        Ue = left_unitary_from_yukawa(Ye_mt)

        # --- Build κ(mt) from the selected angles + a minimal neutrino spectrum choice ---
        th12 = float(np.deg2rad(best.angles.theta12_deg))
        th13 = float(np.deg2rad(best.angles.theta13_deg))
        th23 = float(np.deg2rad(best.angles.theta23_deg))
        dcp = float(np.deg2rad(best.angles.delta_cp_deg))
        U_pmns_target = _pmns_pdg(th12, th13, th23, dcp)

        # Neutrino mass spectrum choice (normal ordering, minimal m1=0).
        dm21_sq_eV2 = 7.4e-5
        dm31_sq_eV2 = 2.5e-3
        m1_eV = 0.0
        m2_eV = float(np.sqrt(dm21_sq_eV2))
        m3_eV = float(np.sqrt(dm31_sq_eV2))
        m_eV = np.array([m1_eV, m2_eV, m3_eV], dtype=float)
        m_GeV = (m_eV * 1.0e-9).astype(float)

        # Enforce PMNS = Ue† Uν by choosing Uν = Ue U_PMNS, and Takagi mν = Uν diag(m) Uν^T.
        U_nu = (Ue @ U_pmns_target).astype(complex)
        mnu_GeV = (U_nu @ np.diag(m_GeV) @ U_nu.T).astype(complex)
        kappa_mt = (mnu_GeV / (float(v_ev_GeV) ** 2)).astype(complex)
        kappa_mt = 0.5 * (kappa_mt + kappa_mt.T)

        # Verify PMNS angles from (Ye, κ) match the target angles at mt.
        U_mt_check, mnu_mt_eV = _pmns_from_ye_and_kappa(Ye=Ye_mt, kappa=kappa_mt, v_ev_GeV=v_ev_GeV)
        angles_mt = _extract_pdg_angles(U_mt_check)

        # --- EFT RG: run κ only (mt → 0.99 MNR1) with fixed mt couplings as a stability diagnostic ---
        thr_raw = json.loads(thresholds_path.read_text(encoding="utf-8"))
        thr = thr_raw.get("thresholds_GeV", {})
        MNR = [float(thr.get("MNR1", 1.0e14)), float(thr.get("MNR2", 3.0e14)), float(thr.get("MNR3", 8.0e14))]
        mu_start = float(sm_bc_mt.mu_mt_GeV)
        mu_end = 0.99 * float(MNR[0])

        Yu_diag = np.diag([0.0, 0.0, float(yt_mt)]).astype(complex)
        Yd_diag = np.diag([0.0, 0.0, float(yb_mt)]).astype(complex)
        yN_off = np.zeros((3, 3), dtype=complex)

        def rhs(t: float, y: np.ndarray) -> np.ndarray:
            # t = ln(mu); y stores [Re(kappa_flat), Im(kappa_flat)]
            k_re = y[:9].reshape(3, 3)
            k_im = y[9:].reshape(3, 3)
            k = (k_re + 1j * k_im).astype(complex)
            dk = _beta_kappa_1loop(g2=g2_mt, lambda_h=lambda_mt, Yu=Yu_diag, Yd=Yd_diag, Ye=Ye_mt, yN=yN_off, kappa=k)
            dk = 0.5 * (dk + dk.T)
            return np.concatenate([dk.real.reshape(-1), dk.imag.reshape(-1)], axis=0).astype(float)

        y0 = np.concatenate([kappa_mt.real.reshape(-1), kappa_mt.imag.reshape(-1)], axis=0).astype(float)
        sol = solve_ivp(
            rhs,
            t_span=(float(np.log(mu_start)), float(np.log(mu_end))),
            y0=y0,
            method="DOP853",
            rtol=1e-9,
            atol=1e-12,
        )
        y_end = sol.y[:, -1]
        kappa_end = (y_end[:9].reshape(3, 3) + 1j * y_end[9:].reshape(3, 3)).astype(complex)
        kappa_end = 0.5 * (kappa_end + kappa_end.T)
        U_end, _mnu_end_eV = _pmns_from_ye_and_kappa(Ye=Ye_mt, kappa=kappa_end, v_ev_GeV=v_ev_GeV)
        angles_end = _extract_pdg_angles(U_end)

        # --- Reconstruct (yN, MNR) from κ via κ = yN M^{-1} yN^T (tree-level) ---
        #
        # Since κ was constructed explicitly from (Uν, m_i), we can factorize it exactly as
        #   κ = (Uν diag(sqrt(m_i))/v) (Uν diag(sqrt(m_i))/v)^T
        # and then assign heavy masses MNRi by column-wise rescaling.
        sqrt_m_GeV = np.sqrt(np.clip(m_GeV, 0.0, float("inf")))
        Ytilde = (U_nu @ np.diag(sqrt_m_GeV / float(v_ev_GeV))).astype(complex)  # κ = Ytilde Ytilde^T

        sqrtM = np.sqrt(np.array(MNR, dtype=float))
        invM = np.diag(1.0 / np.array(MNR, dtype=float))

        best_perm = (0, 1, 2)
        best_maxabs = float("inf")
        yN_best = np.zeros((3, 3), dtype=complex)
        for perm in permutations((0, 1, 2)):
            yN = (Ytilde[:, perm] @ np.diag(sqrtM)).astype(complex)
            mabs = float(np.max(np.abs(yN)))
            if mabs < best_maxabs:
                best_maxabs = mabs
                best_perm = perm
                yN_best = yN

        kappa_rec = (yN_best @ invM @ yN_best.T).astype(complex)
        kappa_rec = 0.5 * (kappa_rec + kappa_rec.T)
        kappa_rec_dev = float(np.max(np.abs(kappa_rec - kappa_mt)))

        # Explicit κ decoupling ladder (integrate-in N_Ri upward)
        cols = [yN_best[:, i].reshape(3, 1) for i in range(3)]
        k_above_1 = kappa_mt - (cols[0] @ cols[0].T) / float(MNR[0])
        k_above_2 = k_above_1 - (cols[1] @ cols[1].T) / float(MNR[1])
        k_above_3 = k_above_2 - (cols[2] @ cols[2].T) / float(MNR[2])
        k_above_3 = 0.5 * (k_above_3 + k_above_3.T)
        kappa_decouple_resid = float(np.max(np.abs(k_above_3)))

        # --- Checks ---
        angle_mismatch_deg = max(
            abs(angles_mt.theta12_deg - best.angles.theta12_deg),
            abs(angles_mt.theta13_deg - best.angles.theta13_deg),
            abs(angles_mt.theta23_deg - best.angles.theta23_deg),
            abs(((angles_mt.delta_cp_deg - best.angles.delta_cp_deg + 540.0) % 360.0) - 180.0),
        )
        drift_deg = max(
            abs(angles_end.theta12_deg - angles_mt.theta12_deg),
            abs(angles_end.theta13_deg - angles_mt.theta13_deg),
            abs(angles_end.theta23_deg - angles_mt.theta23_deg),
            abs(((angles_end.delta_cp_deg - angles_mt.delta_cp_deg + 540.0) % 360.0) - 180.0),
        )

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="selected_variant_defined",
                passed=bool(best.name in [v[0] for v in variants]),
                detail=f"selected={best.name} (selection: min|Δδ|, prefer Δθ23>0, then max Δθ23)",
            )
        )
        checks.append(
            Check(
                check_id="kappa_symmetric",
                passed=bool(float(np.max(np.abs(kappa_mt - kappa_mt.T))) < 1e-12),
                detail=f"max|kappa-kappa^T|(mt)={float(np.max(np.abs(kappa_mt - kappa_mt.T))):.3e}",
            )
        )
        checks.append(
            Check(
                check_id="pmns_angles_reproduced_at_mt",
                passed=bool(angle_mismatch_deg < 1e-6),
                detail=f"max angle mismatch at mt = {angle_mismatch_deg:.3e} deg",
            )
        )
        checks.append(
            Check(
                check_id="kappa_factorization_recovers_kappa",
                passed=bool(kappa_rec_dev < 1e-12 * float(np.max(np.abs(kappa_mt)))),
                detail=f"max|kappa_rec-kappa|(mt)={kappa_rec_dev:.3e} GeV^-1 (perm={best_perm}, max|yN|={best_maxabs:.3e})",
            )
        )
        checks.append(
            Check(
                check_id="kappa_decouples_above_MNR3_tree_level",
                passed=bool(kappa_decouple_resid < 1e-12 * float(np.max(np.abs(kappa_mt)))),
                detail=f"max|kappa_above_MNR3|={kappa_decouple_resid:.3e} GeV^-1",
            )
        )
        checks.append(
            Check(
                check_id="pmns_angles_stable_under_kappa_running",
                passed=bool(drift_deg < 0.05),
                detail=f"max drift (mt→0.99 MNR1) = {drift_deg:.3e} deg",
            )
        )

        # --- Report ---
        lines: list[str] = []
        lines += [
            "PMNS mechanism bridge (Z3-breaking angles → κ(mt) → yN+MNR reconstruction)",
            "",
            f"sm inputs file: {_relpath(sm_path)} (sha256={_sha256_file(sm_path)})",
            f"flavor texture config: {_relpath(flavor_cfg_path)} (sha256={_sha256_file(flavor_cfg_path)})",
            f"lepton masses file: {_relpath(lepton_masses_path)} (sha256={_sha256_file(lepton_masses_path)})",
            f"thresholds file: {_relpath(thresholds_path)} (sha256={_sha256_file(thresholds_path)})",
            "",
            "Z3-breaking constraint (pmns_z3_breaking logic, deterministic selection):",
            f"- baseline: θ12={base_angles.theta12_deg:.4f}°, θ13={base_angles.theta13_deg:.4f}°, θ23={base_angles.theta23_deg:.4f}°, δ={base_angles.delta_cp_deg:.2f}°",
            f"- epsilon = varphi0/6 = {eps:.6g} rad = {float(np.degrees(eps)):.3f}°",
            f"- selected variant: {best.name}",
            f"  angles: θ12={best.angles.theta12_deg:.4f}°, θ13={best.angles.theta13_deg:.4f}°, θ23={best.angles.theta23_deg:.4f}°, δ={best.angles.delta_cp_deg:.2f}°",
            "",
            "Charged-lepton basis (Ye texture, mt):",
            f"- delta_source={delta_source}, delta_M={delta_M:.9f}, delta_star={delta_star:.9f}, delta_used={delta_used:.9f}",
            f"- phase_mode={phase_mode}, theta(delta_used)={theta:.6g} rad, zeta=exp(i theta)={zeta}",
            f"- ytau(mt) target={ytau_mt_target:.6g} => y*_e={ystar_e:.6g}",
            "",
            "SM boundary @ mt (matching layer):",
            f"- v(mt)={v_ev_GeV:.6g} GeV, g2(mt)={g2_mt:.6g}, lambda(mt)={lambda_mt:.6g}",
            f"- (diag approx for κ beta traces) yt(mt)={yt_mt:.6g}, yb(mt)={yb_mt:.6g}",
            "",
            "Reconstructed κ(mt) (Majorana EFT; mν=v^2 κ):",
            f"- neutrino masses used (normal ordering, m1=0): m(eV)={m_eV.tolist()}",
            f"- max|kappa|(mt)={float(np.max(np.abs(kappa_mt))):.3e} GeV^-1",
            f"- PMNS from (Ye,κ) @ mt: θ12={angles_mt.theta12_deg:.4f}°, θ13={angles_mt.theta13_deg:.4f}°, θ23={angles_mt.theta23_deg:.4f}°, δ={angles_mt.delta_cp_deg:.2f}°",
            f"- neutrino masses proxy from SVD(mν) (eV, sorted): {mnu_mt_eV.tolist()}",
            "",
            "κ running stability (EFT below MNR1; couplings frozen at mt):",
            f"- run: mu_start=mt={mu_start:.6g} GeV → mu_end=0.99*MNR1={mu_end:.6g} GeV",
            f"- PMNS from (Ye,κ) @ mu_end: θ12={angles_end.theta12_deg:.4f}°, θ13={angles_end.theta13_deg:.4f}°, θ23={angles_end.theta23_deg:.4f}°, δ={angles_end.delta_cp_deg:.2f}°",
            "",
            "Seesaw inversion / reconstruction (tree-level):",
            f"- MNR (GeV) = {MNR}",
            f"- reconstructed yN (perm={best_perm}) max|yN|={best_maxabs:.3e}",
            f"- kappa_rec dev: max|kappa_rec-kappa|={kappa_rec_dev:.3e} GeV^-1",
            f"- kappa decoupling residual above MNR3: max|kappa_above_MNR3|={kappa_decouple_resid:.3e} GeV^-1",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This is a bridge layer: it takes Z3-breaking angle constraints and produces a consistent κ(mt) and a (yN,MNR) factorization.",
            "- The (yN,MNR) reconstruction is not unique (Casas–Ibarra freedom); here we choose a minimal Takagi/SVD factorization and a permutation that minimizes max|yN|.",
        ]

        return ModuleResult(
            results={
                "selected_variant": {
                    "name": best.name,
                    "angles_deg": best.angles.__dict__,
                    "delta_theta23_deg": best.delta_theta23_deg,
                    "delta_delta_deg": best.delta_delta_deg,
                },
                "kappa_mt": _mat_real_imag(kappa_mt),
                "kappa_running": {
                    "mu_start_GeV": mu_start,
                    "mu_end_GeV": mu_end,
                    "angles_deg_mt": angles_mt.__dict__,
                    "angles_deg_mu_end": angles_end.__dict__,
                    "max_angle_drift_deg": drift_deg,
                },
                "ye_texture": {
                    "delta_source": delta_source,
                    "phase_mode": phase_mode,
                    "delta_M": delta_M,
                    "delta_star": delta_star,
                    "delta_used": delta_used,
                    "theta_rad": theta,
                    "y_star_e": ystar_e,
                },
                "neutrino_masses_input_eV": m_eV.tolist(),
                "thresholds": {"MNR": MNR},
                "reconstruction": {
                    "yN": _mat_real_imag(yN_best),
                    "perm": list(best_perm),
                    "maxabs_yN": best_maxabs,
                    "kappa_rec_dev_GeVinv": kappa_rec_dev,
                    "kappa_above_MNR1_maxabs_GeVinv": float(np.max(np.abs(k_above_1))),
                    "kappa_above_MNR2_maxabs_GeVinv": float(np.max(np.abs(k_above_2))),
                    "kappa_above_MNR3_maxabs_GeVinv": kappa_decouple_resid,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

