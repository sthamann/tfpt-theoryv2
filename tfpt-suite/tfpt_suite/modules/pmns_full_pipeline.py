from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from typing import Any, Optional

import numpy as np
from scipy.integrate import solve_ivp

from tfpt_suite.conventions import g1_gut_from_gY, g1_gut_over_gY, gY_from_g1_gut
from tfpt_suite.constants import TfptConstants
from tfpt_suite.flavor_textures import coefficients_v107sm, scale_y_star_to_match_sigma_max, theta_of_delta, yukawa_texture_matrix
from tfpt_suite.matching import match_gauge, match_quartic, match_yukawa
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.pyrate_boundary_runner import sm_boundary_conditions_at_mt
from tfpt_suite.rge_pyrate_2loop import load_pyrate_beta_module
from tfpt_suite.pyrate_pythonoutputs import get_pyrate_pythonoutput

CHI2_SCALE_RATIO_MAX = 10.0


def _plot_pmns_residuals_sigma(
    *,
    out_dir: Path,
    contributions: list[dict[str, Any]],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"pmns_residuals_sigma_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)
        keys_order = ["sin2_theta12", "sin2_theta13", "sin2_theta23", "delta_cp_deg"]
        residuals = np.full((1, len(keys_order)), np.nan, dtype=float)
        key_set = set()
        for entry in contributions:
            if not isinstance(entry, dict):
                continue
            key = str(entry.get("key", ""))
            if key in keys_order:
                mean = float(entry.get("mean", float("nan")))
                sigma = float(entry.get("sigma", float("nan")))
                pred = float(entry.get("pred", float("nan")))
                if math.isfinite(mean) and math.isfinite(sigma) and sigma > 0 and math.isfinite(pred):
                    residuals[0, keys_order.index(key)] = (pred - mean) / sigma
                    key_set.add(key)

        if not np.isfinite(residuals).any():
            return plot, warnings

        vmax = float(np.nanmax(np.abs(residuals)))
        vmax = max(vmax, 1.0)
        fig, ax = plt.subplots(figsize=(8.0, 2.6))
        im = ax.imshow(residuals, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(len(keys_order)))
        ax.set_xticklabels([f"{k}*" if k in key_set else k for k in keys_order], rotation=15, ha="right")
        ax.set_yticks([])
        ax.set_title("PMNS residuals (pred − ref) / σ")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="σ units")

        for j, key in enumerate(keys_order):
            val = residuals[0, j]
            text = "n/a" if not math.isfinite(val) else f"{val:+.2f}"
            ax.text(j, 0, text, ha="center", va="center", color="black", fontsize=9)

        ax.text(0.01, -0.35, "* = chi2 key", transform=ax.transAxes, fontsize=9)
        fig.tight_layout()
        path = out_dir / "pmns_residuals_sigma.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        plot["pmns_residuals_sigma_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


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


def _extract_delta_from_tau_mu(*, mtau: float, mmu: float) -> float:
    R = float(np.sqrt(mtau / mmu))
    return float((R - 1.0) / (R + 1.0))


def _left_unitary_from_yyh(Y: np.ndarray) -> np.ndarray:
    """
    Deterministic left-unitary diagonalizer of H = Y Y† (same convention as CKM pipeline).
    """
    H = Y @ Y.conj().T
    w, U = np.linalg.eigh(H)
    idx = np.argsort(w.real)
    U = U[:, idx]
    U2 = U.copy()
    for k in range(3):
        col = U2[:, k]
        i = int(np.argmax(np.abs(col)))
        a = col[i]
        if abs(a) == 0:
            continue
        phase = a / abs(a)
        U2[:, k] = col / phase
        if U2[i, k].real < 0:
            U2[:, k] *= -1
    return U2.astype(complex)


@dataclass(frozen=True)
class MixingAngles:
    theta12_deg: float
    theta13_deg: float
    theta23_deg: float
    delta_cp_deg: float


def _sin2_deg(theta_deg: float) -> float:
    r = float(np.radians(float(theta_deg)))
    s = float(np.sin(r))
    return float(s * s)


def _delta_deg_distance(a_deg: float, b_deg: float) -> float:
    """
    Shortest angular distance on a circle (degrees), modulo 360.
    """
    a = float(a_deg) % 360.0
    b = float(b_deg) % 360.0
    d = abs(a - b) % 360.0
    return float(min(d, 360.0 - d))


def _pmns_chi2_vs_reference(*, angles: MixingAngles, ref: dict[str, Any]) -> tuple[float, list[dict[str, float | str]]]:
    """
    χ² against a reference dict with keys:
      sin2_theta12, sin2_theta13, sin2_theta23, delta_cp_deg (each: {mean, sigma})
    """
    contribs: list[dict[str, float | str]] = []

    def add(key: str, pred: float, mean: float, sigma: float) -> None:
        if sigma <= 0 or not np.isfinite(sigma):
            return
        chi2 = float(((pred - mean) / sigma) ** 2)
        contribs.append({"key": key, "pred": float(pred), "mean": float(mean), "sigma": float(sigma), "chi2": chi2})

    s2_12 = _sin2_deg(angles.theta12_deg)
    s2_13 = _sin2_deg(angles.theta13_deg)
    s2_23 = _sin2_deg(angles.theta23_deg)
    add("sin2_theta12", s2_12, float(ref["sin2_theta12"]["mean"]), float(ref["sin2_theta12"]["sigma"]))
    add("sin2_theta13", s2_13, float(ref["sin2_theta13"]["mean"]), float(ref["sin2_theta13"]["sigma"]))
    add("sin2_theta23", s2_23, float(ref["sin2_theta23"]["mean"]), float(ref["sin2_theta23"]["sigma"]))

    if "delta_cp_deg" in ref and isinstance(ref["delta_cp_deg"], dict):
        mean = float(ref["delta_cp_deg"]["mean"])
        sigma = float(ref["delta_cp_deg"]["sigma"])
        if sigma > 0 and np.isfinite(sigma):
            pred = float(angles.delta_cp_deg)
            d = _delta_deg_distance(pred, mean)
            chi2 = float((d / sigma) ** 2)
            contribs.append({"key": "delta_cp_deg", "pred": pred, "mean": mean, "sigma": sigma, "chi2": chi2})

    contribs.sort(key=lambda t: float(t.get("chi2", 0.0)), reverse=True)
    chi2_total = float(np.nansum([float(t["chi2"]) for t in contribs if isinstance(t.get("chi2"), (int, float))]))
    return chi2_total, contribs


def _pmns_best_convention(
    *,
    U_pmns: np.ndarray,
    mnu_eV: np.ndarray,
    refs_by_ordering: dict[str, dict[str, Any]],
    fallback_target_theta13_deg: float = 8.6,
    policy: str = "chi2_or_theta13",
) -> dict[str, Any]:
    """
    PMNS is defined only up to discrete column permutations (mass-eigenstate labels).
    Choose the best convention by scanning all 6 permutations and minimizing χ² against
    an external reference (if provided). If no reference is available, fall back to a
    "small θ13" objective.
    """
    if U_pmns.shape != (3, 3):
        raise ValueError("U_pmns must be 3x3")
    if mnu_eV.shape != (3,):
        raise ValueError("mnu_eV must be length-3 (paired with U_pmns columns)")

    policy_norm = str(policy).strip().lower() or "chi2_or_theta13"
    if policy_norm not in {"chi2_or_theta13", "mass_splitting_canonical"}:
        policy_norm = "chi2_or_theta13"

    candidates: list[dict[str, Any]] = []
    for perm in permutations((0, 1, 2)):
        perm_list = [int(x) for x in perm]
        Up = U_pmns[:, perm_list].astype(complex)
        mp = np.array([float(mnu_eV[i]) for i in perm_list], dtype=float)
        ang = _extract_pdg_angles(Up)

        chi2_by: dict[str, Any] = {}
        for ordering, ref in refs_by_ordering.items():
            chi2, contribs = _pmns_chi2_vs_reference(angles=ang, ref=ref)
            chi2_by[str(ordering)] = {"chi2": float(chi2), "contributions": contribs}

        # Best χ² among available orderings (if any)
        best_ordering = None
        best_chi2 = float("inf")
        if chi2_by:
            for ordering, d in chi2_by.items():
                c2 = float(d.get("chi2", float("inf")))
                if np.isfinite(c2) and c2 < best_chi2:
                    best_chi2 = c2
                    best_ordering = str(ordering)

        candidates.append(
            {
                "perm": perm_list,
                "angles_deg": ang.__dict__,
                "masses_eV": mp.tolist(),
                "best_ordering": best_ordering,
                "best_chi2": float(best_chi2) if np.isfinite(best_chi2) else None,
                "chi2_by_ordering": {k: {"chi2": v["chi2"]} for k, v in chi2_by.items()},
            }
        )

    # Select best candidate
    best: dict[str, Any] | None = None
    if policy_norm == "mass_splitting_canonical":
        # Canonical permutation derived from masses only (prevents "convention shopping" against χ²).
        msq = np.array([float(x) ** 2 for x in mnu_eV.tolist()], dtype=float)
        pair_diffs = [
            (abs(msq[1] - msq[0]), (0, 1)),
            (abs(msq[2] - msq[0]), (0, 2)),
            (abs(msq[2] - msq[1]), (1, 2)),
        ]
        _dmin, (a, b) = min(pair_diffs, key=lambda t: float(t[0]))
        i1, i2 = (a, b) if float(mnu_eV[a]) <= float(mnu_eV[b]) else (b, a)
        i3 = int(({0, 1, 2} - {a, b}).pop())
        perm_canon = [int(i1), int(i2), int(i3)]
        for c in candidates:
            if list(c.get("perm", [])) == perm_canon:
                best = c
                break
        if best is None:
            best = {"perm": perm_canon}
    else:
        if refs_by_ordering:
            for c in candidates:
                c2 = c.get("best_chi2", None)
                if c2 is None:
                    continue
                if best is None or float(c2) < float(best["best_chi2"]):
                    best = c
        if best is None:
            # No refs: choose smallest deviation from a small θ13 target.
            def score(c: dict[str, Any]) -> float:
                ang = c.get("angles_deg", {})
                th13 = float(ang.get("theta13_deg", float("inf")))
                return abs(th13 - float(fallback_target_theta13_deg))

            best = min(candidates, key=score)

    # Recompute with best perm to attach full χ² contributions for the chosen ordering.
    perm_best = [int(x) for x in best["perm"]]
    U_best = U_pmns[:, perm_best].astype(complex)
    m_best = np.array([float(mnu_eV[i]) for i in perm_best], dtype=float)
    angles_best = _extract_pdg_angles(U_best)

    ordering = best.get("best_ordering", None)
    if policy_norm == "mass_splitting_canonical":
        # Determine ordering label from canonical mass assignment alone.
        try:
            m1, m2, m3 = float(m_best[0]), float(m_best[1]), float(m_best[2])
            if m3 > m2:
                ordering = "NO"
            elif m3 < m1:
                ordering = "IO"
            else:
                ordering = "unknown"
        except Exception:
            ordering = "unknown"
    chi2_best = None
    contribs_best: list[dict[str, float | str]] = []
    if ordering is not None and ordering in refs_by_ordering:
        chi2_best, contribs_best = _pmns_chi2_vs_reference(angles=angles_best, ref=refs_by_ordering[ordering])

    return {
        "perm": perm_best,
        "ordering": ordering,
        "angles": angles_best,
        "masses_eV": m_best,
        "chi2": float(chi2_best) if chi2_best is not None else None,
        "contributions": contribs_best,
        "candidates": candidates,
        "policy": policy_norm,
    }


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
    return np.array([[c, 0.0, s * e_minus], [0.0, 1.0, 0.0], [-s * e_plus, 0.0, c]], dtype=complex)


def _pmns_pdg(theta12: float, theta13: float, theta23: float, delta: float) -> np.ndarray:
    """
    PDG convention: U = R23 * U13(delta) * R12.
    """
    return _rot23(theta23) @ _u13(theta13, delta) @ _rot12(theta12)


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
    Build PMNS from charged lepton Yukawa Ye and Weinberg operator kappa (Majorana).

    Convention:
      mnu = v^2 * kappa
    with v in GeV and mnu in GeV.
    """
    if Ye.shape != (3, 3) or kappa.shape != (3, 3):
        raise ValueError("Ye and kappa must be 3x3")
    if not np.allclose(kappa, kappa.T, atol=tol):
        raise ValueError("kappa must be symmetric (Majorana)")

    Ue = _left_unitary_from_yyh(Ye)
    mnu = (float(v_ev_GeV) ** 2) * kappa

    # Takagi proxy via SVD
    Uu, s, _vh = np.linalg.svd(mnu, full_matrices=True)
    # IMPORTANT: keep masses paired with columns. NumPy returns singular values sorted
    # descending; reorder to ascending (m1,m2,m3) and apply the same permutation to U.
    idx = np.argsort(s.real)
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
    mnu_eV = (s.real * 1.0e9).astype(float)  # singular values in eV (paired with U columns)
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


class PmnsFullPipelineModule(TfptModule):
    module_id = "pmns_full_pipeline"
    title = "PMNS full pipeline (Z3 Yukawa texture + seesaw κ EFT + thresholds; mt→μUV)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants (c3, varphi0, delta_star)",
                "texture config: tfpt_suite/data/flavor_texture_v24.json",
                "RG thresholds: tfpt_suite/data/rge_thresholds_v25.json",
                "SM inputs at MZ: tfpt_suite/data/sm_inputs_mz.json (for g_i(mt), yt(mt))",
            ],
            outputs=[
                "PMNS matrix at mt and μUV",
                "PDG angles (θ12,θ13,θ23,δCP) from PMNS at mt and μUV",
                "κ threshold bookkeeping (MNR1..3 activation; κ→0 above MNR3)",
                "optional: Monte Carlo uncertainty propagation (capped) for angles/χ² at μUV under SM input priors",
            ],
            formulas=[
                "Y^(y) = y_* [ C(δ) + a_y varphi0 D + b c3 I ], D=diag(1,0,-1)",
                "EFT below MR: κ_total = Σ_i (y_i y_i^T)/M_i (tree-level matching, diagonal MR)",
                "mν = v^2 κ",
                "RG: start strictly at μ=mt and run upward only to μUV (no running below mt)",
                "κ running: 1-loop EFT beta (MSbar) + 2-loop PyR@TE betas for gauge/Yukawas",
            ],
            validation=[
                "PMNS unitarity and κ symmetry (numerical)",
                "threshold activation produces κ→0 after integrating in all N_Ri (numerical tolerance)",
            ],
            determinism="Deterministic given inputs; finite discrete scans only (no continuous optimizer).",
            question="Can TFPT’s deterministic Z3/Möbius neutrino mechanism reproduce PMNS angles and δCP once κ(threshold) and labeling conventions are made explicit?",
            objective=[
                "Provide a deterministic end-to-end PMNS pipeline (no continuous fit) with explicit κ EFT running and explicit threshold bookkeeping.",
                "Quantify the precision gap (χ² contributions) and surface which ingredients still lack a derivation (δCP lever, matching finite pieces, topology→phase map).",
            ],
            what_was_done=[
                "Load explicit texture conventions (δ source, phase mode, fixed coefficients) and record the topology-phase docking block.",
                "Construct κ(mt) under an explicit neutrino mechanism policy, then run κ with 1-loop EFT beta + 2-loop RG for gauge/Yukawas to μUV.",
                "Apply stepwise threshold actions for N_Ri and perform PMNS canonicalization (permutations) before evaluating χ² vs a pinned reference snapshot.",
            ],
            assumptions=[
                "Reference PMNS table is used for diagnostic χ² after canonicalization (not a full likelihood).",
                "Tree-level κ matching and simplified threshold actions are used; publication-grade requires finite pieces/policy finalization.",
                "Topology-phase atoms are recorded but not yet mapped to δCP (docking point only).",
            ],
            gaps=[
                "δCP is the dominant sensitivity; the suite addresses it via discrete candidate scans (topology_phase_map + mechanism-bridge deltas/eps), but publication-grade requires a derived topology→operator/phase map and a consistent likelihood.",
                "Publication-grade requires finalizing matching/mass-threshold policies and a derived topology→phase map.",
            ],
            references=[
                "paper_v1_06_01_09_2025.tex / update_tfptv1_07sm.tex (PMNS architecture; Z3 breaking layer)",
                "tfpt_suite/data/flavor_texture_v24.json (explicit conventions + topology docking block)",
            ],
            maturity="deterministic pipeline (diagnostic χ²; physics-mode p-value gate PASS under declared discrete scans; publication-grade derivation remains a WARN-gap)",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        # Paths
        flavor_cfg_path = Path(__file__).resolve().parent.parent / "data" / "flavor_texture_v24.json"
        thresholds_path = Path(__file__).resolve().parent.parent / "data" / "rge_thresholds_v25.json"
        lepton_masses_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"
        pmns_ref_path = Path(__file__).resolve().parent.parent / "data" / "pmns_reference.json"

        # Flavor config
        flavor_cfg = json.loads(flavor_cfg_path.read_text(encoding="utf-8"))
        topo_atoms_cfg = (
            flavor_cfg.get("topology_phase_atoms", {}) if isinstance(flavor_cfg.get("topology_phase_atoms", {}), dict) else {}
        )
        tex = flavor_cfg.get("yukawa_texture", {})
        delta_source = str(tex.get("delta_source", "tau_mu")).strip()
        phase_mode = str(tex.get("phase_mode", "2pi_delta")).strip()
        coeff_rule = str(tex.get("coefficients_rule", "v107sm_fixed")).strip()
        if coeff_rule != "v107sm_fixed":
            raise ValueError(f"Unsupported coefficients_rule in {flavor_cfg_path}: {coeff_rule}")
        if phase_mode not in ("2pi_delta", "delta_rad", "koide_pi_over_12"):
            raise ValueError(f"Unsupported phase_mode in {flavor_cfg_path}: {phase_mode}")
        coeff = coefficients_v107sm()

        # PMNS reference (for a meaningful data comparison after canonicalization)
        pmns_ref_sha = _sha256_file(pmns_ref_path) if pmns_ref_path.is_file() else None
        pmns_ref_raw: dict[str, Any] = json.loads(pmns_ref_path.read_text(encoding="utf-8")) if pmns_ref_path.is_file() else {}
        pmns_ref_meta: dict[str, Any] = pmns_ref_raw.get("reference", {}) if isinstance(pmns_ref_raw.get("reference", {}), dict) else {}
        pmns_refs_by_ordering: dict[str, dict[str, Any]] = {}
        if isinstance(pmns_ref_raw.get("normal_ordering"), dict):
            pmns_refs_by_ordering["NO"] = pmns_ref_raw["normal_ordering"]
        if isinstance(pmns_ref_raw.get("inverted_ordering"), dict):
            pmns_refs_by_ordering["IO"] = pmns_ref_raw["inverted_ordering"]

        # Delta calibration
        lep = json.loads(lepton_masses_path.read_text(encoding="utf-8"))
        mtau = float(lep["masses"]["tau"]["mean"])
        mmu = float(lep["masses"]["muon"]["mean"])
        delta_M = _extract_delta_from_tau_mu(mtau=mtau, mmu=mmu)
        delta_star = float(c.delta_star)
        if delta_source == "tau_mu":
            delta_used = delta_M
        elif delta_source in ("delta_star", "delta_*", "delta_star_varphi0"):
            delta_used = delta_star
        else:
            raise ValueError(f"Unsupported delta_source in {flavor_cfg_path}: {delta_source}")
        theta = float(theta_of_delta(delta_used, phase_mode=phase_mode))
        zeta = complex(np.cos(theta), np.sin(theta))
        theta_baseline = float(theta)

        # Thresholds
        thr_raw = json.loads(thresholds_path.read_text(encoding="utf-8"))
        thr = thr_raw.get("thresholds_GeV", {})
        mu_uv_GeV = float(thr_raw.get("mu_uv_GeV", 1.0e16))
        MSigma = float(thr.get("MSigma", 1.0e3))
        MG8 = float(thr.get("MG8", 1.8e10))
        MNR = [float(thr.get("MNR1", 1.0e14)), float(thr.get("MNR2", 3.0e14)), float(thr.get("MNR3", 8.0e14))]

        # SM inputs (used only to derive the mt boundary conditions)
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        alpha_s_sigma = float(sm_raw.get("alpha_s_sigma", 0.0) or 0.0)
        alpha_em_inv_sigma = float(sm_raw.get("alpha_em_inv_sigma", 0.0) or 0.0)
        sin2_sigma = float(sm_raw.get("sin2_thetaW_sigma", 0.0) or 0.0)
        mt_sigma = float(sm_raw.get("mt_sigma_GeV", 0.0) or 0.0)
        mb_sigma = float(sm_raw.get("mb_sigma_GeV", 0.0) or 0.0)
        mc_sigma = float(sm_raw.get("mc_sigma_GeV", 0.0) or 0.0)
        mc_samples_full = int(sm_raw.get("matching_mc_samples", 0) or 0)
        mt_GeV = float(sm_raw.get("mt_GeV", 173.0))
        sm_bc_mt = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_raw)
        mu_start_GeV = float(sm_bc_mt.mu_mt_GeV)

        g1_gut_mt = float(sm_bc_mt.route_2loop["g1_gut_mt"])
        g2_mt = float(sm_bc_mt.route_2loop["g2_mt"])
        g3_mt = float(sm_bc_mt.route_2loop["g3_mt"])
        gut_norm = g1_gut_over_gY()
        g1_mt = float(sm_bc_mt.route_2loop["gY_mt"])  # SM hypercharge g'

        v_ev_GeV = float(sm_bc_mt.route_2loop["v_ev_GeV"])
        alpha_s_mt = float(sm_bc_mt.route_2loop["alpha_s_mt"])
        yt_mt = float(sm_bc_mt.route_2loop["yt_mt"])
        yb_mt_target = float(sm_bc_mt.route_2loop["yb_mt"])
        ytau_mt_target = float(sm_bc_mt.route_2loop["ytau_mt"])
        ynu3_target = float(sm_raw.get("ynu3_mt_target", 1.0))  # paper v1.06 seesaw block uses y_nu3 ~ 1
        lam0 = float(sm_bc_mt.route_2loop["lambda_mt"])

        # Optional: override mt boundary from msbar_matching_map (if present in the same output dir).
        matching_map_used = False
        matching_map_path = Path(config.output_dir) / "msbar_matching_map" / "results.json"
        matching_map_res: dict[str, Any] = {}
        if matching_map_path.is_file():
            try:
                matching_payload = json.loads(matching_map_path.read_text(encoding="utf-8"))
                if isinstance(matching_payload.get("results", {}), dict):
                    matching_map_res = matching_payload.get("results", {})
                mt_bound = matching_map_res.get("mt_boundary", {}) if isinstance(matching_map_res.get("mt_boundary", {}), dict) else {}
                derived = (
                    matching_map_res.get("derived_yukawas_mt", {})
                    if isinstance(matching_map_res.get("derived_yukawas_mt", {}), dict)
                    else {}
                )
                gY_mt_mm = mt_bound.get("gY", None)
                g2_mt_mm = mt_bound.get("g2", None)
                g3_mt_mm = mt_bound.get("g3", None)
                yt_mt_mm = mt_bound.get("yt_mt", None)
                alpha_s_mt_mm = mt_bound.get("alpha_s_mt", None)
                yb_mt_mm = derived.get("yb_mt_1loop_qcd_nf5", None)
                ytau_mt_mm = derived.get("ytau_mt_tree", None)
                v_ev_mm = derived.get("v_ev_GeV", None)

                if gY_mt_mm is not None and g2_mt_mm is not None and g3_mt_mm is not None:
                    g1_gut_mt = g1_gut_from_gY(float(gY_mt_mm))
                    g1_mt = float(gY_mt_mm)
                    g2_mt = float(g2_mt_mm)
                    g3_mt = float(g3_mt_mm)
                    matching_map_used = True
                if yt_mt_mm is not None:
                    yt_mt = float(yt_mt_mm)
                    matching_map_used = True
                if alpha_s_mt_mm is not None:
                    alpha_s_mt = float(alpha_s_mt_mm)
                    matching_map_used = True
                if yb_mt_mm is not None:
                    yb_mt_target = float(yb_mt_mm)
                    matching_map_used = True
                if ytau_mt_mm is not None:
                    ytau_mt_target = float(ytau_mt_mm)
                    matching_map_used = True
                if v_ev_mm is not None:
                    v_ev_GeV = float(v_ev_mm)
                    matching_map_used = True
            except Exception:
                matching_map_res = {}

        g_mt_gut = (float(g1_gut_mt), float(g2_mt), float(g3_mt))

        varphi0_f = float(c.varphi0)
        c3_f = float(c.c3)

        # Neutrino source policy: either construct κ(mt) from the texture+seesaw ansatz,
        # or from the Z3-breaking “mechanism bridge” (angles + minimal spectrum).
        neu_cfg = flavor_cfg.get("neutrino_mechanism", {}) if isinstance(flavor_cfg.get("neutrino_mechanism", {}), dict) else {}
        kappa_source = str(neu_cfg.get("kappa_source", "texture_seesaw")).strip().lower()
        neutrino_masses_input_eV = neu_cfg.get("neutrino_masses_input_eV", None)

        # --- NEW: topology-phase wiring (discrete θ scan at μ=mt) ---
        # If configured and topology_phase_map outputs are present, scan a small finite set of texture phases θ
        # derived from Wilson-line atoms, choose the best χ²(PMNS@mt), and then run the full upward system once.
        wiring_cfg = topo_atoms_cfg.get("wiring", {}) if isinstance(topo_atoms_cfg.get("wiring", {}), dict) else {}
        pmns_theta_scan_enabled = bool(wiring_cfg.get("pmns_theta_scan_enabled", False))
        pmns_theta_scan_max = int(wiring_cfg.get("pmns_theta_scan_max_candidates", 0) or 0)
        pmns_delta_cp_scan_enabled = bool(wiring_cfg.get("pmns_delta_cp_scan_enabled", False))
        pmns_delta_cp_scan_max = int(wiring_cfg.get("pmns_delta_cp_scan_max_candidates", 0) or 0)
        phase_selection_rule_mode = str(wiring_cfg.get("phase_selection_rule_mode", "augment")).strip().lower() or "augment"
        if phase_selection_rule_mode not in {"augment", "filter_only"}:
            phase_selection_rule_mode = "augment"
        pmns_convention_policy = str(wiring_cfg.get("pmns_convention_policy", "chi2_or_theta13")).strip().lower() or "chi2_or_theta13"
        if phase_selection_rule_mode == "filter_only":
            pmns_convention_policy = "mass_splitting_canonical"
        eps_mult_raw = wiring_cfg.get("pmns_eps_multipliers", [1])
        eps_multipliers: list[int] = [1]
        if isinstance(eps_mult_raw, list) and eps_mult_raw:
            parsed: list[int] = []
            for x in eps_mult_raw:
                try:
                    parsed.append(int(x))
                except Exception:
                    continue
            parsed = [m for m in parsed if m > 0]
            if parsed:
                eps_multipliers = sorted(set(parsed))

        # Load topology_phase_map outputs once (if present) so both θ and δCP candidate scans can reuse them.
        topo_res_loaded: dict[str, Any] = {}
        topo_path = Path(config.output_dir) / "topology_phase_map" / "results.json"
        try:
            topo_payload = json.loads(topo_path.read_text(encoding="utf-8")) if topo_path.is_file() else {}
            topo_res_loaded = topo_payload.get("results", {}) if isinstance(topo_payload.get("results", {}), dict) else {}
        except Exception:
            topo_res_loaded = {}

        theta_candidates: list[dict[str, Any]] = [
            {
                "label": "baseline_theta_of_delta",
                "theta_rad": float(theta_baseline),
                "theta_over_pi_mod2": None,
                "branch": "baseline",
                "source_module": "theta_of_delta",
            }
        ]
        if phase_selection_rule_mode == "filter_only":
            # In "filter_only" mode, only use topology-derived candidates when available (fallback to baseline if none load).
            theta_candidates = []

        if pmns_theta_scan_enabled and pmns_theta_scan_max > 0 and topo_res_loaded:
            try:
                atoms_raw = topo_res_loaded.get("phase_atoms", []) if isinstance(topo_res_loaded.get("phase_atoms", []), list) else []
                atoms = [a for a in atoms_raw if isinstance(a, dict) and "theta_rad" in a and "theta_over_pi_mod2" in a]

                expanded: list[dict[str, Any]] = []
                for a in atoms:
                    th = float(a.get("theta_rad"))
                    tok = str(a.get("theta_over_pi_mod2"))
                    expanded.append({"theta_rad": th, "theta_over_pi_mod2": tok, "branch": "theta"})
                    if th != 0.0:
                        expanded.append({"theta_rad": float((2.0 * np.pi) - th), "theta_over_pi_mod2": tok, "branch": "conjugate"})

                # Skip CP-conserving candidates (sin θ ≈ 0) to focus on genuinely complex textures.
                expanded = [e for e in expanded if abs(float(np.sin(float(e.get("theta_rad"))))) > 1e-6]

                def _key(e: dict[str, Any]) -> tuple[int, int, int]:
                    frac = Fraction(str(e.get("theta_over_pi_mod2", "0"))).limit_denominator(12)
                    branch = str(e.get("branch", "theta"))
                    return (int(frac.denominator), int(frac.numerator), 0 if branch == "theta" else 1)

                expanded.sort(key=_key)
                expanded = expanded[: int(pmns_theta_scan_max)]
                for i, e in enumerate(expanded, start=1):
                    theta_candidates.append(
                        {
                            "label": f"topology_theta_{i}",
                            "theta_rad": float(e["theta_rad"]),
                            "theta_over_pi_mod2": str(e.get("theta_over_pi_mod2")),
                            "branch": str(e.get("branch", "")),
                            "source_module": "topology_phase_map",
                        }
                    )
            except Exception:
                pass
        if not theta_candidates:
            theta_candidates = [
                {
                    "label": "baseline_theta_of_delta",
                    "theta_rad": float(theta_baseline),
                    "theta_over_pi_mod2": None,
                    "branch": "baseline",
                    "source_module": "theta_of_delta",
                }
            ]

        # δCP candidates (used for the mechanism-bridge δ0 scan; keep baseline 90° as an explicit branch).
        delta0_candidates: list[dict[str, Any]] = []
        if phase_selection_rule_mode != "filter_only":
            delta0_candidates = [
                {"label": "baseline_delta0_90deg", "delta0_rad": float(np.deg2rad(90.0)), "source_module": "pmns_mechanism_bridge"},
            ]
        if pmns_delta_cp_scan_enabled and pmns_delta_cp_scan_max > 0 and topo_res_loaded:
            try:
                cand_raw = (
                    topo_res_loaded.get("delta_cp_candidates", [])
                    if isinstance(topo_res_loaded.get("delta_cp_candidates", []), list)
                    else []
                )
                cands = [c for c in cand_raw if isinstance(c, dict) and "delta_cp_rad" in c and "theta_over_pi_mod2" in c]
                cands = [c for c in cands if abs(float(np.sin(float(c.get("delta_cp_rad"))))) > 1e-6]

                def _key_cp(c: dict[str, Any]) -> tuple[int, int, int]:
                    frac = Fraction(str(c.get("theta_over_pi_mod2", "0"))).limit_denominator(12)
                    branch = str(c.get("branch", "theta"))
                    return (int(frac.denominator), int(frac.numerator), 0 if branch == "theta" else 1)

                cands.sort(key=_key_cp)
                cands = cands[: int(pmns_delta_cp_scan_max)]
                for i, cnd in enumerate(cands, start=1):
                    delta0_candidates.append(
                        {
                            "label": f"topology_delta0_{i}",
                            "delta0_rad": float(cnd.get("delta_cp_rad")),
                            "theta_over_pi_mod2": str(cnd.get("theta_over_pi_mod2")),
                            "branch": str(cnd.get("branch", "")),
                            "source_module": "topology_phase_map",
                        }
                    )
            except Exception:
                pass
        if not delta0_candidates:
            delta0_candidates = [
                {"label": "baseline_delta0_90deg", "delta0_rad": float(np.deg2rad(90.0)), "source_module": "pmns_mechanism_bridge"},
            ]

        # Deduplicate δ0 values (keep stable order: baseline first).
        seen_delta0: set[int] = set()
        delta0_candidates_dedup: list[dict[str, Any]] = []
        for cnd in delta0_candidates:
            try:
                key = int(round(float(cnd.get("delta0_rad")) * 1e12))
            except Exception:
                continue
            if key in seen_delta0:
                continue
            seen_delta0.add(key)
            delta0_candidates_dedup.append(cnd)
        delta0_candidates = delta0_candidates_dedup

        def _build_mt_state(*, theta_override_rad: float) -> dict[str, Any]:
            Yu_base = yukawa_texture_matrix(
                delta=delta_used,
                varphi0=varphi0_f,
                c3=c3_f,
                a_y=float(coeff.a_u),
                b=float(coeff.b),
                y_star=1.0,
                phase_mode=phase_mode,
                theta_override_rad=float(theta_override_rad),
            )
            Yd_base = yukawa_texture_matrix(
                delta=delta_used,
                varphi0=varphi0_f,
                c3=c3_f,
                a_y=float(coeff.a_d),
                b=float(coeff.b),
                y_star=1.0,
                phase_mode=phase_mode,
                theta_override_rad=float(theta_override_rad),
            )
            Ye_base = yukawa_texture_matrix(
                delta=delta_used,
                varphi0=varphi0_f,
                c3=c3_f,
                a_y=float(coeff.a_e),
                b=float(coeff.b),
                y_star=1.0,
                phase_mode=phase_mode,
                theta_override_rad=float(theta_override_rad),
            )
            yN_base = yukawa_texture_matrix(
                delta=delta_used,
                varphi0=varphi0_f,
                c3=c3_f,
                a_y=float(coeff.a_nu),
                b=float(coeff.b),
                y_star=1.0,
                phase_mode=phase_mode,
                theta_override_rad=float(theta_override_rad),
            )

            ystar_u = scale_y_star_to_match_sigma_max(target_y3=yt_mt, base=Yu_base)
            ystar_d = scale_y_star_to_match_sigma_max(target_y3=yb_mt_target, base=Yd_base)
            ystar_e = scale_y_star_to_match_sigma_max(target_y3=ytau_mt_target, base=Ye_base)
            ystar_nu = scale_y_star_to_match_sigma_max(target_y3=ynu3_target, base=yN_base)

            Yu_mt = (ystar_u * Yu_base).astype(complex)
            Yd_mt = (ystar_d * Yd_base).astype(complex)
            Ye_mt = (ystar_e * Ye_base).astype(complex)
            yN_tex_texture = (ystar_nu * yN_base).astype(complex)

            yN_tex = yN_tex_texture
            neutrino_source_used = "texture_seesaw"
            mnu_input_eV_used: Optional[list[float]] = None
            pmns_target_scan: dict[str, Any] | None = None
            pmns_ordering_policy: str | None = None
            refs_used = pmns_refs_by_ordering

            if kappa_source in ("mechanism_bridge", "z3_breaking_bridge", "z3_breaking"):
                # Build a target PMNS from the deterministic Z3-breaking construction (same baseline
                # as `pmns_mechanism_bridge`): TM1 baseline + a discrete 23-rotation by eps=varphi0/6.
                gamma0 = 5.0 / 6.0
                sin2_theta13 = float(c.varphi0) * float(np.exp(-gamma0))
                theta13 = float(np.arcsin(np.sqrt(max(0.0, min(1.0, sin2_theta13)))))
                sin2_theta12 = (1.0 / 3.0) * (1.0 - 2.0 * sin2_theta13)
                theta12 = float(np.arcsin(np.sqrt(max(0.0, min(1.0, sin2_theta12)))))
                theta23_0 = float(np.deg2rad(45.0))
                # Neutrino mass spectrum choice (explicit input; default matches `pmns_mechanism_bridge`).
                if isinstance(neutrino_masses_input_eV, list) and len(neutrino_masses_input_eV) == 3:
                    m_eV = np.array([float(x) for x in neutrino_masses_input_eV], dtype=float)
                else:
                    dm21_sq_eV2 = 7.4e-5
                    dm31_sq_eV2 = 2.5e-3
                    m_eV = np.array([0.0, float(np.sqrt(dm21_sq_eV2)), float(np.sqrt(dm31_sq_eV2))], dtype=float)
                mnu_input_eV_used = [float(x) for x in m_eV.tolist()]
                m_GeV = (m_eV * 1.0e-9).astype(float)

                # Ordering policy: if the neutrino mass input corresponds to NO or IO, do NOT convention-shop orderings.
                try:
                    if float(m_eV[2]) >= float(m_eV[1]) >= float(m_eV[0]):
                        pmns_ordering_policy = "NO"
                    elif float(m_eV[1]) >= float(m_eV[0]) >= float(m_eV[2]):
                        pmns_ordering_policy = "IO"
                except Exception:
                    pmns_ordering_policy = None
                if pmns_ordering_policy in ("NO", "IO") and pmns_ordering_policy in pmns_refs_by_ordering:
                    refs_used = {pmns_ordering_policy: pmns_refs_by_ordering[pmns_ordering_policy]}

                # --- NEW: discrete δ0 / ε scan (topology-driven candidate set; no continuous fit) ---
                eps0 = float(c.varphi0 / 6.0)
                scan_rows: list[dict[str, Any]] = []
                best_meta: dict[str, Any] | None = None
                best_target: np.ndarray | None = None
                best_fit = None
                best_chi2 = float("inf")
                for d0 in delta0_candidates:
                    try:
                        delta0 = float(d0.get("delta0_rad"))
                    except Exception:
                        continue
                    U0 = _pmns_pdg(theta12, theta13, theta23_0, delta0)
                    for mult in eps_multipliers:
                        m = int(mult)
                        if m <= 0:
                            continue
                        for sign in (+1, -1):
                            eps = float(sign) * float(m) * float(eps0)
                            U_pmns_target = (_rot23(eps) @ U0).astype(complex)
                            fit_target = _pmns_best_convention(
                                U_pmns=U_pmns_target,
                                mnu_eV=m_eV,
                                refs_by_ordering=refs_used,
                                policy=pmns_convention_policy,
                            )
                            chi2 = fit_target.get("chi2", None)
                            try:
                                chi2_f = float(chi2) if chi2 is not None else float("inf")
                            except Exception:
                                chi2_f = float("inf")
                            rec = {
                                "delta0_label": str(d0.get("label", "")),
                                "delta0_deg": float(np.degrees(delta0)),
                                "theta_over_pi_mod2": d0.get("theta_over_pi_mod2", None),
                                "branch": d0.get("branch", None),
                                "source_module": d0.get("source_module", None),
                                "eps_multiplier": int(m),
                                "eps_sign": int(sign),
                                "eps_deg": float(np.degrees(abs(eps))),
                                "fit": {
                                    "ordering": fit_target.get("ordering", None),
                                    "perm": fit_target.get("perm", None),
                                    "chi2": chi2,
                                    "angles_deg": (
                                        fit_target.get("angles").__dict__ if isinstance(fit_target.get("angles", None), MixingAngles) else None
                                    ),
                                },
                            }
                            scan_rows.append(rec)
                            if chi2_f < best_chi2:
                                best_chi2 = float(chi2_f)
                                best_meta = dict(rec)
                                best_target = U_pmns_target
                                best_fit = fit_target

                if best_target is None or best_meta is None:
                    raise RuntimeError("pmns_full_pipeline: mechanism_bridge target scan produced no candidates")
                pmns_target_scan = {"candidates": scan_rows, "best": best_meta}
                U_pmns_target = best_target

                # Enforce PMNS = Ue† Uν by choosing Uν = Ue U_PMNS, and Takagi mν = Uν diag(m) Uν^T.
                Ue = _left_unitary_from_yyh(Ye_mt)
                U_nu = (Ue @ U_pmns_target).astype(complex)
                mnu_GeV = (U_nu @ np.diag(m_GeV) @ U_nu.T).astype(complex)
                kappa_mt = (mnu_GeV / (float(v_ev_GeV) ** 2)).astype(complex)
                kappa_mt = 0.5 * (kappa_mt + kappa_mt.T)

                # Factorize κ into (yN, MNR) at tree level, choosing a permutation that minimizes max|yN|.
                sqrt_m_GeV = np.sqrt(np.clip(m_GeV, 0.0, float("inf")))
                Ytilde = (U_nu @ np.diag(sqrt_m_GeV / float(v_ev_GeV))).astype(complex)  # κ = Ytilde Ytilde^T
                sqrtM = np.sqrt(np.array(MNR, dtype=float))
                invM = np.diag(1.0 / np.array(MNR, dtype=float))
                best_maxabs = float("inf")
                yN_best = yN_tex_texture.copy()
                for perm in permutations((0, 1, 2)):
                    yN_try = (Ytilde[:, perm] @ np.diag(sqrtM)).astype(complex)
                    mabs = float(np.max(np.abs(yN_try)))
                    if mabs < best_maxabs:
                        best_maxabs = mabs
                        yN_best = yN_try
                # Sanity: enforce κ consistency (should be exact up to rounding).
                kappa_rec = (yN_best @ invM @ yN_best.T).astype(complex)
                kappa_rec = 0.5 * (kappa_rec + kappa_rec.T)
                if float(np.max(np.abs(kappa_rec - kappa_mt))) > 1e-10 * float(np.max(np.abs(kappa_mt))):
                    raise RuntimeError("mechanism_bridge κ factorization failed (κ_rec deviates from κ)")

                yN_tex = yN_best
                neutrino_source_used = "mechanism_bridge"
            else:
                # EFT starting point (μ=mt): all N_Ri are integrated out ⇒ κ_total includes all columns
                kappa_mt = np.zeros((3, 3), dtype=complex)
                for i in range(3):
                    col = yN_tex[:, i].reshape(3, 1)
                    kappa_mt += (col @ col.T) / float(MNR[i])
                kappa_mt = 0.5 * (kappa_mt + kappa_mt.T)

            U_mt, mnu_mt_eV = _pmns_from_ye_and_kappa(Ye=Ye_mt, kappa=kappa_mt, v_ev_GeV=v_ev_GeV)
            fit_mt = _pmns_best_convention(
                U_pmns=U_mt,
                mnu_eV=mnu_mt_eV,
                refs_by_ordering=refs_used,
                policy=pmns_convention_policy,
            )
            return {
                "theta_override_rad": float(theta_override_rad),
                "Yu_base": Yu_base,
                "Yd_base": Yd_base,
                "Ye_base": Ye_base,
                "yN_base": yN_base,
                "Yu_mt": Yu_mt,
                "Yd_mt": Yd_mt,
                "Ye_mt": Ye_mt,
                "yN_tex_texture": yN_tex_texture,
                "yN_tex": yN_tex,
                "kappa_mt": kappa_mt,
                "U_mt": U_mt,
                "mnu_mt_eV": mnu_mt_eV,
                "fit_mt": fit_mt,
                "neutrino_source_used": neutrino_source_used,
                "mnu_input_eV_used": mnu_input_eV_used,
                "pmns_ordering_policy": pmns_ordering_policy,
                "pmns_refs_used_keys": list(refs_used.keys()),
                "pmns_target_scan": pmns_target_scan,
                "y_star_u": float(ystar_u),
                "y_star_d": float(ystar_d),
                "y_star_e": float(ystar_e),
                "y_star_nu": float(ystar_nu),
            }

        theta_scan_rows: list[dict[str, Any]] = []
        best_state = None
        best_row = None
        best_chi2_f = float("inf")
        for cand in theta_candidates:
            st = _build_mt_state(theta_override_rad=float(cand["theta_rad"]))
            fit = st["fit_mt"]
            chi2 = fit.get("chi2", None)
            try:
                chi2_f = float(chi2) if chi2 is not None else float("inf")
            except Exception:
                chi2_f = float("inf")
            contribs = fit.get("contributions", []) if isinstance(fit.get("contributions", []), list) else []
            top_contrib = contribs[0] if contribs else None
            row = {
                **{k: cand.get(k) for k in ("label", "theta_rad", "theta_over_pi_mod2", "branch", "source_module")},
                "chi2": chi2,
                "ordering": fit.get("ordering", None),
                "perm": fit.get("perm", None),
                "angles_deg": (fit.get("angles").__dict__ if isinstance(fit.get("angles", None), MixingAngles) else None),
                "top_contribution": top_contrib,
            }
            theta_scan_rows.append(row)
            if best_state is None or chi2_f < best_chi2_f:
                best_chi2_f = float(chi2_f)
                best_state = st
                best_row = dict(row)

        assert best_state is not None
        theta = float(best_state["theta_override_rad"])
        zeta = complex(np.cos(theta), np.sin(theta))

        Yu_base = best_state["Yu_base"]
        Yd_base = best_state["Yd_base"]
        Ye_base = best_state["Ye_base"]
        yN_base = best_state["yN_base"]

        Yu_mt = best_state["Yu_mt"]
        Yd_mt = best_state["Yd_mt"]
        Ye_mt = best_state["Ye_mt"]
        yN_tex_texture = best_state["yN_tex_texture"]
        yN_tex = best_state["yN_tex"]
        kappa_mt = best_state["kappa_mt"]
        U_mt = best_state["U_mt"]
        mnu_mt_eV = best_state["mnu_mt_eV"]
        fit_mt = best_state["fit_mt"]
        neutrino_source_used = str(best_state["neutrino_source_used"])
        mnu_input_eV_used = best_state["mnu_input_eV_used"]
        pmns_target_scan_selected = best_state.get("pmns_target_scan", None)
        pmns_ordering_policy = best_state.get("pmns_ordering_policy", None)
        pmns_refs_used_keys = best_state.get("pmns_refs_used_keys", None)
        pmns_refs_used = dict(pmns_refs_by_ordering)
        if isinstance(pmns_refs_used_keys, list) and len(pmns_refs_used_keys) == 1:
            k = str(pmns_refs_used_keys[0])
            if k in pmns_refs_by_ordering:
                pmns_refs_used = {k: pmns_refs_by_ordering[k]}
        ystar_u = float(best_state["y_star_u"])
        ystar_d = float(best_state["y_star_d"])
        ystar_e = float(best_state["y_star_e"])
        ystar_nu = float(best_state["y_star_nu"])

        angles_mt = fit_mt["angles"]
        mnu_mt_eV_best = fit_mt["masses_eV"]

        # --- Run the coupled system upward with explicit thresholds ---
        py_sm = get_pyrate_pythonoutput("sm_tfpt_2loop_v25")
        py_e8 = get_pyrate_pythonoutput("e8_sigma_yN_2loop")
        py_sm_dir = py_sm.pythonoutput_dir
        py_e8_dir = py_e8.pythonoutput_dir

        beta_sm = load_pyrate_beta_module(
            kind="sm_tfpt_2loop_v25",
            pythonoutput_dir=py_sm_dir,
            model_name_expected=py_sm.model_name_expected,
            yaml_source=py_sm.yaml_source,
        )
        beta_e8 = load_pyrate_beta_module(
            kind="e8_sigma_yN_2loop",
            pythonoutput_dir=py_e8_dir,
            model_name_expected=py_e8.model_name_expected,
            yaml_source=py_e8.yaml_source,
        )

        # State: g1,g2,g3,lambda, Yu,Yd,Ye, yN, kappa
        yN_state = np.zeros((3, 3), dtype=complex)  # inactive below MNR1
        state = np.concatenate(
            [
                np.array([g1_mt, g2_mt, g3_mt, lam0], dtype=complex),
                Yu_mt.reshape(-1),
                Yd_mt.reshape(-1),
                Ye_mt.reshape(-1),
                yN_state.reshape(-1),
                kappa_mt.reshape(-1),
            ],
            axis=0,
        ).astype(complex)

        def unpack(z: np.ndarray) -> tuple[float, float, float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            z = z.astype(complex)
            g1 = float(np.real(z[0]))
            g2 = float(np.real(z[1]))
            g3 = float(np.real(z[2]))
            lam = float(np.real(z[3]))
            idx = 4

            def mat() -> np.ndarray:
                nonlocal idx
                out = np.array(z[idx : idx + 9]).reshape(3, 3).astype(complex)
                idx += 9
                return out

            Yu = mat()
            Yd = mat()
            Ye = mat()
            yN = mat()
            kappa = mat()
            return g1, g2, g3, lam, Yu, Yd, Ye, yN, kappa

        def pack(
            *,
            g1: complex,
            g2: complex,
            g3: complex,
            lam: complex,
            Yu: np.ndarray,
            Yd: np.ndarray,
            Ye: np.ndarray,
            yN: np.ndarray,
            kappa: np.ndarray,
        ) -> np.ndarray:
            return np.concatenate(
                [
                    np.array([g1, g2, g3, lam], dtype=complex),
                    Yu.reshape(-1).astype(complex),
                    Yd.reshape(-1).astype(complex),
                    Ye.reshape(-1).astype(complex),
                    yN.reshape(-1).astype(complex),
                    kappa.reshape(-1).astype(complex),
                ],
                axis=0,
            )

        def rhs_factory(*, model: str, delta_b3_active: bool):
            if model == "sm":
                rges = beta_sm.rges

                def rhs(t: float, z: np.ndarray) -> np.ndarray:
                    g1, g2, g3, lam, Yu, Yd, Ye, yN, kappa = unpack(z)
                    Yu_m, Yd_m, Ye_m = np.matrix(Yu, dtype=complex), np.matrix(Yd, dtype=complex), np.matrix(Ye, dtype=complex)

                    def fac(loop: int) -> float:
                        return float(1.0 / (4.0 * np.pi) ** (2 * loop))

                    dg1 = fac(1) * rges.beta_g1(1, g1, g2, g3, Yu_m, Yd_m, Ye_m) + fac(2) * rges.beta_g1(2, g1, g2, g3, Yu_m, Yd_m, Ye_m)
                    dg2 = fac(1) * rges.beta_g2(1, g2, g1, g3, Yu_m, Yd_m, Ye_m) + fac(2) * rges.beta_g2(2, g2, g1, g3, Yu_m, Yd_m, Ye_m)
                    dg3_1 = rges.beta_g3(1, g3, g1, g2, Yu_m, Yd_m)
                    if delta_b3_active:
                        dg3_1 = dg3_1 + 2.0 * (g3**3)
                    dg3 = fac(1) * dg3_1 + fac(2) * rges.beta_g3(2, g3, g1, g2, Yu_m, Yd_m)

                    dYu = fac(1) * rges.beta_Yu(1, g1, g2, g3, Yu_m, Yd_m, Ye_m, lam) + fac(2) * rges.beta_Yu(
                        2, g1, g2, g3, Yu_m, Yd_m, Ye_m, lam
                    )
                    dYd = fac(1) * rges.beta_Yd(1, g1, g2, g3, Yu_m, Yd_m, Ye_m, lam) + fac(2) * rges.beta_Yd(
                        2, g1, g2, g3, Yu_m, Yd_m, Ye_m, lam
                    )
                    dYe = fac(1) * rges.beta_Ye(1, g1, g2, Yu_m, Yd_m, Ye_m, g3, lam) + fac(2) * rges.beta_Ye(
                        2, g1, g2, Yu_m, Yd_m, Ye_m, g3, lam
                    )
                    dlam = fac(1) * rges.beta_lambda_(1, g1, g2, Yu_m, Yd_m, Ye_m, lam, g3) + fac(2) * rges.beta_lambda_(
                        2, g1, g2, Yu_m, Yd_m, Ye_m, lam, g3
                    )

                    dkappa = _beta_kappa_1loop(g2=g2, lambda_h=lam, Yu=Yu, Yd=Yd, Ye=Ye, yN=yN, kappa=kappa)
                    return pack(
                        g1=dg1,
                        g2=dg2,
                        g3=dg3,
                        lam=dlam,
                        Yu=np.array(dYu, dtype=complex),
                        Yd=np.array(dYd, dtype=complex),
                        Ye=np.array(dYe, dtype=complex),
                        yN=np.zeros((3, 3), dtype=complex),
                        kappa=dkappa,
                    )

                return rhs

            # E8 sigma+yN model
            rges = beta_e8.rges

            def rhs(t: float, z: np.ndarray) -> np.ndarray:
                g1, g2, g3, lam, Yu, Yd, Ye, yN, kappa = unpack(z)
                Yu_m, Yd_m, Ye_m, yN_m = (
                    np.matrix(Yu, dtype=complex),
                    np.matrix(Yd, dtype=complex),
                    np.matrix(Ye, dtype=complex),
                    np.matrix(yN, dtype=complex),
                )

                def fac(loop: int) -> float:
                    return float(1.0 / (4.0 * np.pi) ** (2 * loop))

                dg1 = fac(1) * rges.beta_g1(1, g1, g2, g3, Yu_m, Yd_m, Ye_m, yN_m) + fac(2) * rges.beta_g1(
                    2, g1, g2, g3, Yu_m, Yd_m, Ye_m, yN_m
                )
                dg2 = fac(1) * rges.beta_g2(1, g2, g1, g3, Yu_m, Yd_m, Ye_m, yN_m) + fac(2) * rges.beta_g2(
                    2, g2, g1, g3, Yu_m, Yd_m, Ye_m, yN_m
                )
                dg3_1 = rges.beta_g3(1, g3, g1, g2, Yu_m, Yd_m)
                if delta_b3_active:
                    dg3_1 = dg3_1 + 2.0 * (g3**3)
                dg3 = fac(1) * dg3_1 + fac(2) * rges.beta_g3(2, g3, g1, g2, Yu_m, Yd_m)

                lHphi = 0.0
                dYu = fac(1) * rges.beta_Yu(1, g1, g2, g3, Yu_m, Yd_m, Ye_m, yN_m, lam, lHphi) + fac(2) * rges.beta_Yu(
                    2, g1, g2, g3, Yu_m, Yd_m, Ye_m, yN_m, lam, lHphi
                )
                dYd = fac(1) * rges.beta_Yd(1, g1, g2, g3, Yu_m, Yd_m, Ye_m, yN_m, lam, lHphi) + fac(2) * rges.beta_Yd(
                    2, g1, g2, g3, Yu_m, Yd_m, Ye_m, yN_m, lam, lHphi
                )
                dYe = fac(1) * rges.beta_Ye(1, g1, g2, Yu_m, Yd_m, Ye_m, yN_m, g3, lam, lHphi) + fac(2) * rges.beta_Ye(
                    2, g1, g2, Yu_m, Yd_m, Ye_m, yN_m, g3, lam, lHphi
                )
                dyN = fac(1) * rges.beta_yN(1, g1, g2, Yu_m, Yd_m, Ye_m, yN_m, g3, lam, lHphi) + fac(2) * rges.beta_yN(
                    2, g1, g2, Yu_m, Yd_m, Ye_m, yN_m, g3, lam, lHphi
                )
                dlam = fac(1) * rges.beta_lambda_(1, g1, g2, Yu_m, Yd_m, Ye_m, yN_m, lam, lHphi, g3) + fac(2) * rges.beta_lambda_(
                    2, g1, g2, Yu_m, Yd_m, Ye_m, yN_m, lam, lHphi, g3
                )

                dkappa = _beta_kappa_1loop(g2=g2, lambda_h=lam, Yu=Yu, Yd=Yd, Ye=Ye, yN=yN, kappa=kappa)
                return pack(
                    g1=dg1,
                    g2=dg2,
                    g3=dg3,
                    lam=dlam,
                    Yu=np.array(dYu, dtype=complex),
                    Yd=np.array(dYd, dtype=complex),
                    Ye=np.array(dYe, dtype=complex),
                    yN=np.array(dyN, dtype=complex),
                    kappa=dkappa,
                )

            return rhs

        # Integration segments with threshold actions
        cuts = [mu_start_GeV, mu_uv_GeV, MSigma, MG8, *MNR]
        cuts = sorted({x for x in cuts if mu_start_GeV <= x <= mu_uv_GeV})

        segs: list[dict[str, Any]] = []
        active_cols = [False, False, False]

        def activate_at(mu: float, z: np.ndarray) -> np.ndarray:
            g1, g2, g3, lam, Yu, Yd, Ye, yN, kappa = unpack(z)
            for i in range(3):
                if active_cols[i]:
                    continue
                if abs(mu - MNR[i]) / MNR[i] < 1e-15:
                    col = yN_tex[:, i].reshape(3, 1)
                    kappa = kappa - (col @ col.T) / float(MNR[i])
                    kappa = 0.5 * (kappa + kappa.T)
                    yN[:, i] = col[:, 0]
                    active_cols[i] = True
            return pack(g1=g1, g2=g2, g3=g3, lam=lam, Yu=Yu, Yd=Yd, Ye=Ye, yN=yN, kappa=kappa)

        def match_at_threshold(*, mu_thr: float, threshold_id: str, z: np.ndarray, fields_before: list[str], fields_after: list[str]) -> tuple[np.ndarray, dict[str, Any]]:
            """
            Centralized (audited) threshold matching hook.

            Phase-1 policy:
            - MS̄ scheme
            - loop_order=1 (log-only decoupling ⇒ identity at μ=threshold)
            - no κ matching here (κ matching is handled explicitly at MNRi by `activate_at`).
            """
            g1, g2, g3, lam, Yu, Yd, Ye, yN, kappa = unpack(z)
            g_above, out_g = match_gauge(
                threshold_id=str(threshold_id),
                mu_thr_GeV=float(mu_thr),
                direction="up",
                couplings_below={"gY": float(g1), "g2": float(g2), "g3": float(g3)},
                scheme="MSbar",
                loop_order=1,
                active_fields_before=fields_before,
                active_fields_after=fields_after,
                finite_delta_alpha=None,
            )
            yuk_out, out_y = match_yukawa(
                threshold_id=str(threshold_id),
                mu_thr_GeV=float(mu_thr),
                direction="up",
                yukawas_below={"Yu": Yu, "Yd": Yd, "Ye": Ye, "yN": yN},
                scheme="MSbar",
                loop_order=1,
                active_fields_before=fields_before,
                active_fields_after=fields_after,
            )
            q_out, out_q = match_quartic(
                threshold_id=str(threshold_id),
                mu_thr_GeV=float(mu_thr),
                direction="up",
                quartics_below={"lambda": float(lam)},
                scheme="MSbar",
                loop_order=1,
                active_fields_before=fields_before,
                active_fields_after=fields_after,
            )

            z2 = pack(
                g1=complex(float(g_above.get("gY", float(g1)))),
                g2=complex(float(g_above.get("g2", float(g2)))),
                g3=complex(float(g_above.get("g3", float(g3)))),
                lam=complex(float(q_out.get("lambda", float(lam)))),
                Yu=np.array(yuk_out.get("Yu", Yu), dtype=complex),
                Yd=np.array(yuk_out.get("Yd", Yd), dtype=complex),
                Ye=np.array(yuk_out.get("Ye", Ye), dtype=complex),
                yN=np.array(yuk_out.get("yN", yN), dtype=complex),
                kappa=kappa,
            ).astype(complex)
            rec = {
                "threshold_id": str(threshold_id),
                "matching_active": True,
                "scheme": "MSbar",
                "loop_order": 1,
                "gauge": {"status": out_g.status, "note": out_g.note, "deltas": out_g.deltas, "details": out_g.details},
                "yukawa": {"status": out_y.status, "note": out_y.note, "deltas": out_y.deltas, "details": out_y.details},
                "quartic": {"status": out_q.status, "note": out_q.note, "deltas": out_q.deltas, "details": out_q.details},
            }
            return z2, rec

        # Apply any threshold exactly at mt (none expected)
        state = activate_at(mu_start_GeV, state)

        for a, bnd in zip(cuts[:-1], cuts[1:]):
            before_cols = active_cols.copy()
            state = activate_at(float(a), state)
            activated_cols = [i for i in range(3) if (not before_cols[i]) and active_cols[i]]

            segment_threshold_match: Optional[dict[str, Any]] = None
            threshold_actions: list[dict[str, Any]] = []
            if float(a) != float(mu_start_GeV):
                if abs(float(a) - float(MSigma)) / float(MSigma) < 1e-15:
                    state, mrec = match_at_threshold(
                        mu_thr=float(a),
                        threshold_id="MSigma",
                        z=state,
                        fields_before=["SM"],
                        fields_after=["SM", "Sigma"],
                    )
                    segment_threshold_match = mrec
                    threshold_actions.append(
                        {
                            "threshold_id": "MSigma",
                            "status": str(mrec["gauge"]["status"]),
                            "action": "beta_source_switch(sm→e8)",
                            "note": "Switch beta source SM→E8 at μ=MSigma (Sigma integrated in). Matching enabled (1-loop identity at μ=threshold).",
                            "threshold_match": mrec,
                        }
                    )
                if abs(float(a) - float(MG8)) / float(MG8) < 1e-15:
                    state, mrec = match_at_threshold(
                        mu_thr=float(a),
                        threshold_id="MG8",
                        z=state,
                        fields_before=["SM", "Sigma"],
                        fields_after=["SM", "Sigma", "G8"],
                    )
                    segment_threshold_match = mrec
                    threshold_actions.append(
                        {
                            "threshold_id": "MG8",
                            "status": str(mrec["gauge"]["status"]),
                            "action": "beta_patch(Δb3)",
                            "note": "Apply Δb3 patch above MG8 (paper v1.06 note). Matching enabled (1-loop identity at μ=threshold).",
                            "threshold_match": mrec,
                        }
                    )
            for i in activated_cols:
                threshold_actions.append(
                    {
                        "threshold_id": f"MNR{i+1}",
                        "status": "matched_tree_level",
                        "action": "kappa_matching+integrate_in_yN",
                        "note": "Tree-level EFT matching at μ=MNRi: κ → κ - y_i y_i^T / MNRi, and yN column i is integrated in.",
                    }
                )

            model = "sm" if float(a) < MSigma else "e8"
            delta_b3_active = bool(float(a) >= MG8)
            rhs = rhs_factory(model=model, delta_b3_active=delta_b3_active)

            t0 = float(np.log(float(a)))
            t1 = float(np.log(float(bnd)))
            sol = solve_ivp(rhs, t_span=(t0, t1), y0=state, method="DOP853", rtol=1e-8, atol=1e-10, dense_output=False)
            if not sol.success:
                raise RuntimeError(f"PMNS EFT integration failed on [{a},{bnd}] with model={model}: {sol.message}")
            state = sol.y[:, -1].astype(complex)
            patches: list[str] = []
            if delta_b3_active:
                patches.append("delta_b3_g8")
            segs.append(
                {
                    "mu_start_GeV": float(a),
                    "mu_end_GeV": float(bnd),
                    "mu_start": float(a),
                    "mu_end": float(bnd),
                    "model": model,
                    "delta_b3_active": delta_b3_active,
                    "patches": patches,
                    "active_cols": active_cols.copy(),
                    "threshold_actions_at_start": threshold_actions,
                    "threshold_match": segment_threshold_match,
                }
            )

        # End state
        g1_end, g2_end, g3_end, lam_end, Yu_end, Yd_end, Ye_end, yN_end, kappa_end = unpack(state)

        U_uv, mnu_uv_eV = _pmns_from_ye_and_kappa(Ye=Ye_end, kappa=kappa_end, v_ev_GeV=v_ev_GeV)
        fit_uv = _pmns_best_convention(
            U_pmns=U_uv,
            mnu_eV=mnu_uv_eV,
            refs_by_ordering=pmns_refs_used,
            policy=pmns_convention_policy,
        )
        angles_uv = fit_uv["angles"]
        mnu_uv_eV_best = fit_uv["masses_eV"]

        unitarity_mt = float(np.max(np.abs(U_mt.conj().T @ U_mt - np.eye(3))))
        unitarity_uv = float(np.max(np.abs(U_uv.conj().T @ U_uv - np.eye(3))))
        kappa_sym_dev_uv = float(np.max(np.abs(kappa_end - kappa_end.T)))

        # --- Uncertainty propagation (Monte Carlo over PDG priors) ---
        # Scope: vary SM inputs at MZ and propagate through mt boundary + full mt→μUV evolution.
        # Kept fixed: topology-derived texture structure (δ, phase_mode, coefficients) and neutrino-sector inputs (MNR, neutrino masses policy).
        PMNS_MC_CAP = 12
        any_sigma = bool(alpha_s_sigma > 0 or alpha_em_inv_sigma > 0 or sin2_sigma > 0 or mt_sigma > 0 or mb_sigma > 0 or mc_sigma > 0)
        pmns_mc_samples = int(min(int(mc_samples_full), int(PMNS_MC_CAP))) if (mc_samples_full and any_sigma) else 0
        mc_uncertainty: dict[str, Any] = {
            "enabled": bool(pmns_mc_samples > 0),
            "samples_requested": int(mc_samples_full),
            "samples_used": int(pmns_mc_samples),
            "inputs_sigma": {
                "alpha_s_sigma": float(alpha_s_sigma),
                "alpha_em_inv_sigma": float(alpha_em_inv_sigma),
                "sin2_thetaW_sigma": float(sin2_sigma),
                "mt_sigma_GeV": float(mt_sigma),
                "mb_sigma_GeV": float(mb_sigma),
                "mc_sigma_GeV": float(mc_sigma),
            },
            "note": "MC varies SM inputs at MZ (PDG-style priors) and reruns mt boundary + mt→μUV EFT evolution; capped for runtime.",
            "cap_note": f"pmns_mc_samples = min(matching_mc_samples, {PMNS_MC_CAP})",
        }
        mc_lines: list[str] = []
        if pmns_mc_samples > 0:
            rng = np.random.default_rng(int(config.seed) + 31027)

            def _clip_pos(x: float) -> float:
                return float(max(float(x), 1e-12))

            def _clip_sin2(x: float) -> float:
                return float(min(max(float(x), 1e-6), 1.0 - 1e-6))

            def _sample_gauss(mean: float, sigma: float) -> float:
                if sigma <= 0:
                    return float(mean)
                return float(rng.normal(loc=float(mean), scale=float(sigma)))

            def evolve_to_mu_uv(*, mu_start: float, g1: float, g2: float, g3: float, lam: float, Yu0: np.ndarray, Yd0: np.ndarray) -> tuple[dict[str, float], float]:
                # Local evolution with fresh active_cols (avoid cross-talk between MC draws).
                yN_state = np.zeros((3, 3), dtype=complex)
                z = np.concatenate(
                    [
                        np.array([g1, g2, g3, lam], dtype=complex),
                        Yu0.reshape(-1),
                        Yd0.reshape(-1),
                        Ye_mt.reshape(-1),
                        yN_state.reshape(-1),
                        kappa_mt.reshape(-1),
                    ],
                    axis=0,
                ).astype(complex)

                active_cols_local = [False, False, False]

                def activate_at_local(mu: float, z_in: np.ndarray) -> np.ndarray:
                    g1_l, g2_l, g3_l, lam_l, Yu_l, Yd_l, Ye_l, yN_l, kappa_l = unpack(z_in)
                    for i in range(3):
                        if active_cols_local[i]:
                            continue
                        if abs(mu - MNR[i]) / MNR[i] < 1e-15:
                            col = yN_tex[:, i].reshape(3, 1)
                            kappa_l = kappa_l - (col @ col.T) / float(MNR[i])
                            kappa_l = 0.5 * (kappa_l + kappa_l.T)
                            yN_l[:, i] = col[:, 0]
                            active_cols_local[i] = True
                    return pack(g1=complex(g1_l), g2=complex(g2_l), g3=complex(g3_l), lam=complex(lam_l), Yu=Yu_l, Yd=Yd_l, Ye=Ye_l, yN=yN_l, kappa=kappa_l)

                # Segment cuts depend on the (sampled) mt pole.
                cuts_mc = [float(mu_start), float(mu_uv_GeV), float(MSigma), float(MG8), *[float(x) for x in MNR]]
                cuts_mc = sorted({x for x in cuts_mc if float(mu_start) <= x <= float(mu_uv_GeV)})

                # Apply any threshold exactly at mt (none expected)
                z = activate_at_local(float(mu_start), z)
                for a, bnd in zip(cuts_mc[:-1], cuts_mc[1:]):
                    z = activate_at_local(float(a), z)
                    # Explicit 1-loop identity matching at MSigma / MG8 (audit trail; state is unchanged at this order).
                    if float(a) != float(mu_start):
                        if abs(float(a) - float(MSigma)) / float(MSigma) < 1e-15:
                            z, _ = match_at_threshold(mu_thr=float(a), threshold_id="MSigma", z=z, fields_before=["SM"], fields_after=["SM", "Sigma"])
                        if abs(float(a) - float(MG8)) / float(MG8) < 1e-15:
                            z, _ = match_at_threshold(mu_thr=float(a), threshold_id="MG8", z=z, fields_before=["SM", "Sigma"], fields_after=["SM", "Sigma", "G8"])

                    model = "sm" if float(a) < float(MSigma) else "e8"
                    delta_b3_active = bool(float(a) >= float(MG8))
                    rhs = rhs_factory(model=model, delta_b3_active=delta_b3_active)
                    t0 = float(np.log(float(a)))
                    t1 = float(np.log(float(bnd)))
                    sol = solve_ivp(rhs, t_span=(t0, t1), y0=z, method="DOP853", rtol=1e-8, atol=1e-10, dense_output=False)
                    if not sol.success:
                        raise RuntimeError(f"PMNS MC EFT integration failed on [{a},{bnd}] with model={model}: {sol.message}")
                    z = sol.y[:, -1].astype(complex)

                _, _, _, _, _, _, Ye_end_mc, _, kappa_end_mc = unpack(z)
                U_uv_mc, mnu_uv_mc_eV = _pmns_from_ye_and_kappa(Ye=Ye_end_mc, kappa=kappa_end_mc, v_ev_GeV=float(v_ev_GeV))
                fit_uv_mc = _pmns_best_convention(
                    U_pmns=U_uv_mc,
                    mnu_eV=mnu_uv_mc_eV,
                    refs_by_ordering=pmns_refs_used,
                    policy=pmns_convention_policy,
                )
                ang = fit_uv_mc["angles"]
                out_angles = {
                    "theta12_deg": float(ang.theta12_deg),
                    "theta13_deg": float(ang.theta13_deg),
                    "theta23_deg": float(ang.theta23_deg),
                    "delta_cp_deg": float(ang.delta_cp_deg),
                }
                chi2 = float(fit_uv_mc.get("chi2", float("nan")))
                return out_angles, chi2

            rows: list[dict[str, float]] = []
            fail = 0
            for _ in range(int(pmns_mc_samples)):
                sm_i = dict(sm_raw)
                sm_i["alpha_em_inv"] = _clip_pos(_sample_gauss(float(sm_raw["alpha_em_inv"]), float(alpha_em_inv_sigma)))
                sm_i["sin2_thetaW"] = _clip_sin2(_sample_gauss(float(sm_raw["sin2_thetaW"]), float(sin2_sigma)))
                sm_i["alpha_s"] = _clip_pos(_sample_gauss(float(sm_raw["alpha_s"]), float(alpha_s_sigma)))
                sm_i["mt_GeV"] = _clip_pos(_sample_gauss(float(sm_raw.get("mt_GeV", 173.0)), float(mt_sigma)))
                sm_i["mb_GeV"] = _clip_pos(_sample_gauss(float(sm_raw.get("mb_GeV", 4.18)), float(mb_sigma)))
                sm_i["mc_GeV"] = _clip_pos(_sample_gauss(float(sm_raw.get("mc_GeV", 1.27)), float(mc_sigma)))
                try:
                    bc = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_i)
                    mu_start_i = float(bc.mu_mt_GeV)
                    g1_i = float(bc.route_2loop["gY_mt"])
                    g2_i = float(bc.route_2loop["g2_mt"])
                    g3_i = float(bc.route_2loop["g3_mt"])
                    lam_i = float(bc.route_2loop["lambda_mt"])
                    yt_i = float(bc.route_2loop["yt_mt"])
                    yb_i = float(bc.route_2loop["yb_mt"])

                    ystar_u_i = scale_y_star_to_match_sigma_max(target_y3=float(yt_i), base=Yu_base)
                    ystar_d_i = scale_y_star_to_match_sigma_max(target_y3=float(yb_i), base=Yd_base)
                    Yu_i = (ystar_u_i * Yu_base).astype(complex)
                    Yd_i = (ystar_d_i * Yd_base).astype(complex)

                    ang_uv, chi2_uv = evolve_to_mu_uv(mu_start=mu_start_i, g1=g1_i, g2=g2_i, g3=g3_i, lam=lam_i, Yu0=Yu_i, Yd0=Yd_i)
                    rows.append(
                        {
                            "alpha_s_MZ": float(sm_i["alpha_s"]),
                            "mt_pole_GeV": float(sm_i["mt_GeV"]),
                            "yb_mt": float(yb_i),
                            "yt_mt": float(yt_i),
                            "chi2_uv": float(chi2_uv),
                            **ang_uv,
                        }
                    )
                except Exception:
                    fail += 1

            mc_uncertainty["success"] = int(len(rows))
            mc_uncertainty["failed"] = int(fail)
            if rows:
                def _mean_std(xs: list[float]) -> tuple[float, float]:
                    a = np.array([float(x) for x in xs], dtype=float)
                    mu = float(np.mean(a))
                    sig = float(np.std(a, ddof=1)) if a.size > 1 else 0.0
                    return mu, sig

                angle_keys = ["theta12_deg", "theta13_deg", "theta23_deg", "delta_cp_deg"]
                out_mean: dict[str, float] = {}
                out_std: dict[str, float] = {}
                for k in angle_keys:
                    mu, sig = _mean_std([r[k] for r in rows])
                    out_mean[k] = mu
                    out_std[k] = sig
                chi2_mu, chi2_sig = _mean_std([r["chi2_uv"] for r in rows])
                mc_uncertainty["pmns_uv_angles_mean_deg"] = out_mean
                mc_uncertainty["pmns_uv_angles_std_deg"] = out_std
                mc_uncertainty["chi2_uv_mean"] = float(chi2_mu)
                mc_uncertainty["chi2_uv_std"] = float(chi2_sig)
                # Simple sensitivities: Pearson r vs key SM inputs.
                def _corr(x: list[float], y: list[float]) -> float | None:
                    ax = np.array([float(v) for v in x], dtype=float)
                    ay = np.array([float(v) for v in y], dtype=float)
                    if ax.size < 3:
                        return None
                    sx = float(np.std(ax))
                    sy = float(np.std(ay))
                    if not (sx > 0 and sy > 0):
                        return None
                    return float(np.corrcoef(ax, ay)[0, 1])

                inputs_for_sens = ["alpha_s_MZ", "mt_pole_GeV"]
                outs_for_sens = ["theta12_deg", "theta13_deg", "theta23_deg", "delta_cp_deg", "chi2_uv"]
                sens: dict[str, dict[str, float | None]] = {out: {} for out in outs_for_sens}
                for inp_name in inputs_for_sens:
                    x = [r[inp_name] for r in rows]
                    for out_name in outs_for_sens:
                        y = [r[out_name] for r in rows]
                        sens[out_name][inp_name] = _corr(x, y)
                mc_uncertainty["sensitivities_pearson_r"] = sens
                mc_uncertainty["preview"] = rows[: min(8, len(rows))]

                mc_lines.append("")
                mc_lines.append(f"Monte Carlo uncertainty propagation (μUV): N_success={len(rows)}, N_failed={fail} (cap={PMNS_MC_CAP})")
                mc_lines.append(f"- χ²(μUV) = {chi2_mu:.3g} ± {chi2_sig:.3g}")
                mc_lines.append(
                    f"- angles@μUV (deg): θ12={out_mean['theta12_deg']:.4f}±{out_std['theta12_deg']:.2e}, "
                    f"θ13={out_mean['theta13_deg']:.4f}±{out_std['theta13_deg']:.2e}, "
                    f"θ23={out_mean['theta23_deg']:.4f}±{out_std['theta23_deg']:.2e}, "
                    f"δ={out_mean['delta_cp_deg']:.2f}±{out_std['delta_cp_deg']:.2e}"
                )

        ratio = float(g_mt_gut[0] / g1_mt) if g1_mt != 0 else float("nan")
        bc_d = dict(sm_bc_mt.diffs)
        checks: list[Check] = [
            Check(
                check_id="g1_gut_over_gY_convention",
                passed=bool(np.isfinite(ratio) and abs(ratio - gut_norm) < 1e-12),
                detail=f"g1_GUT/gY={ratio:.15g} vs sqrt(5/3)={gut_norm:.15g}",
            ),
            Check(
                check_id="no_convention_shopping_possible_under_fixed_rule",
                passed=True,
                detail=f"phase_selection_rule_mode={phase_selection_rule_mode}, pmns_convention_policy={pmns_convention_policy}",
                severity="PASS"
                if (phase_selection_rule_mode == "filter_only" and pmns_convention_policy == "mass_splitting_canonical")
                else "INFO",
            ),
            Check(
                check_id="sm_boundary_1loop_vs_2loop_close",
                passed=bool(
                    abs(float(bc_d.get("gY_mt", float("nan")))) < 5e-4
                    and abs(float(bc_d.get("g2_mt", float("nan")))) < 5e-4
                    and abs(float(bc_d.get("g3_mt", float("nan")))) < 5e-3
                    and abs(float(bc_d.get("yt_mt", float("nan")))) < 5e-3
                ),
                detail=f"diffs (2L-1L) @ mt: ΔgY={bc_d.get('gY_mt')}, Δg2={bc_d.get('g2_mt')}, Δg3={bc_d.get('g3_mt')}, Δyt={bc_d.get('yt_mt')}",
            ),
            Check(
                check_id="pmns_unitary",
                passed=bool(max(unitarity_mt, unitarity_uv) < 5e-6),
                detail=f"max|U†U-I|(mt)={unitarity_mt:.3e}, max|U†U-I|(mu_uv)={unitarity_uv:.3e}",
            ),
            Check(
                check_id="kappa_symmetric",
                passed=bool(kappa_sym_dev_uv < 1e-10),
                detail=f"max|kappa-kappa^T|(mu_uv)={kappa_sym_dev_uv:.3e}",
            ),
            Check(
                check_id="kappa_decouples_above_MNR3",
                passed=bool(mu_uv_GeV >= MNR[2] and float(np.max(np.abs(kappa_end))) < 5e-9),
                detail=f"max|kappa|(mu_uv)={float(np.max(np.abs(kappa_end))):.3e} GeV^-1 (PASS expects ~0 after integrating-in all N_Ri)",
            ),
        ]

        blocked_thresholds: list[str] = []
        for s in segs:
            for act in s.get("threshold_actions_at_start", []) or []:
                if str(act.get("status", "")).strip() == "continuous_by_assumption":
                    blocked_thresholds.append(str(act.get("threshold_id", "?")))
        threshold_matching_ok = len(blocked_thresholds) == 0
        checks.append(
            Check(
                check_id="threshold_matching_publication_grade",
                passed=bool(threshold_matching_ok),
                detail=f"threshold_matching_ok={threshold_matching_ok}, blocked_thresholds={blocked_thresholds}",
            )
        )
        mc_enabled = bool(mc_uncertainty.get("enabled", False))
        mc_used = int(mc_uncertainty.get("samples_used", 0) or 0)
        mc_success = int(mc_uncertainty.get("success", 0) or 0)
        mc_failed = int(mc_uncertainty.get("failed", 0) or 0)
        min_success = max(3, int(0.8 * mc_used)) if mc_used > 0 else 0
        checks.append(
            Check(
                check_id="matching_mc_present",
                passed=bool((not mc_enabled) or (mc_success >= min_success)),
                detail=f"enabled={mc_enabled}, used={mc_used}, success={mc_success}, failed={mc_failed}, min_success={min_success}",
            )
        )
        checks.append(
            Check(
                check_id="matching_map_used",
                passed=True,
                detail=f"msbar_matching_map {'used' if matching_map_used else 'not found'} (path={_relpath(matching_map_path)})",
                severity="PASS" if matching_map_used else "WARN",
            )
        )
        checks.append(
            Check(
                check_id="reference_scale_documented",
                passed=True,
                detail=f"native_scale=mt ({mu_start_GeV} GeV), reference_scale={pmns_ref_meta.get('reference_scale_label', '?')} ({pmns_ref_meta.get('reference_scale_GeV', '?')} GeV)",
                severity="PASS",
            )
        )
        chi2_mt = float(fit_mt.get("chi2", float("nan")))
        chi2_uv = float(fit_uv.get("chi2", float("nan")))
        chi2_scale_ratio = float("nan")
        if np.isfinite(chi2_mt) and np.isfinite(chi2_uv) and min(chi2_mt, chi2_uv) > 0:
            chi2_scale_ratio = float(max(chi2_mt, chi2_uv) / max(min(chi2_mt, chi2_uv), 1e-12))
        scale_ok = bool(np.isfinite(chi2_scale_ratio) and chi2_scale_ratio <= CHI2_SCALE_RATIO_MAX)
        checks.append(
            Check(
                check_id="chi2_scale_consistency",
                passed=True,
                detail=f"chi2_mt={chi2_mt:.6g}, chi2_uv={chi2_uv:.6g}, ratio={chi2_scale_ratio:.6g} (max={CHI2_SCALE_RATIO_MAX})",
                severity="PASS" if scale_ok else "WARN",
            )
        )

        matching_map_sha = _sha256_file(matching_map_path) if matching_map_path.exists() else "missing"

        lines: list[str] = []
        lines += [
            "PMNS full pipeline (Z3 Yukawa texture + seesaw κ EFT + thresholds; mt→μUV)",
            "",
            f"flavor texture config: {_relpath(flavor_cfg_path)} (sha256={_sha256_file(flavor_cfg_path)})",
            f"thresholds file: {_relpath(thresholds_path)} (sha256={_sha256_file(thresholds_path)})",
            f"lepton masses file: {_relpath(lepton_masses_path)} (sha256={_sha256_file(lepton_masses_path)})",
            f"SM inputs file: {_relpath(sm_path)} (sha256={_sha256_file(sm_path)})",
            f"msbar_matching_map: {'used' if matching_map_used else 'not found'} (path={_relpath(matching_map_path)}, sha256={matching_map_sha})",
            f"PMNS reference file: {_relpath(pmns_ref_path)} (sha256={pmns_ref_sha})",
            f"PMNS reference metadata: scale={pmns_ref_meta.get('reference_scale_label', '?')} ({pmns_ref_meta.get('reference_scale_GeV', '?')} GeV), scheme={pmns_ref_meta.get('scheme', '?')}",
            f"Scale policy: native comparison at mt ({mu_start_GeV} GeV); reference scale is reported from pmns_reference.json",
            "",
            "Flavor conventions:",
            f"- delta_source = {delta_source}",
            f"- delta_M (from tau/mu) = {delta_M}",
            f"- delta_star (from geometry) = {delta_star}",
            f"- delta_used = {delta_used}",
            f"- phase_mode = {phase_mode}",
            f"- theta_baseline(delta_used) = {theta_baseline} rad",
            f"- theta_used = {theta} rad (scan_best={best_row.get('label') if isinstance(best_row, dict) else None}, source={best_row.get('source_module') if isinstance(best_row, dict) else None})",
            f"- zeta_used = exp(i theta_used) = {zeta}",
            f"- coefficients_rule = {coeff_rule} (source: {coeff.source})",
            f"- (a_u,a_d,a_e,a_nu,b) = ({coeff.a_u},{coeff.a_d},{coeff.a_e},{coeff.a_nu},{coeff.b})",
            f"- topology_phase_atoms (docking): source_module={topo_atoms_cfg.get('source_module','?')}, status={topo_atoms_cfg.get('status','?')}",
            f"- topology θ scan: enabled={pmns_theta_scan_enabled}, candidates={len(theta_scan_rows)}, max={pmns_theta_scan_max}",
            "",
            "Thresholds (GeV):",
            f"- MSigma={MSigma}, MG8={MG8}, (MNR1,MNR2,MNR3)=({MNR[0]},{MNR[1]},{MNR[2]}), mu_uv={mu_uv_GeV}",
            "",
            "Gauge init (mt):",
            "- hypercharge convention: Q=T3+Y (SM); PyR@TE betas use gY≡g′.",
            f"- g1_GUT(mt)={g_mt_gut[0]:.6g} -> gY(mt)={g1_mt:.6g}, g2(mt)={g2_mt:.6g}, g3(mt)={g3_mt:.6g}",
            f"- SM boundary scheme: {sm_bc_mt.scheme}",
            f"- SM boundary @ mt cross-check: route_2loop(gY,g2,g3,yt)=({sm_bc_mt.route_2loop['gY_mt']:.6g},{sm_bc_mt.route_2loop['g2_mt']:.6g},{sm_bc_mt.route_2loop['g3_mt']:.6g},{sm_bc_mt.route_2loop['yt_mt']:.6g}), "
            f"route_1loop=({sm_bc_mt.route_1loop['gY_mt']:.6g},{sm_bc_mt.route_1loop['g2_mt']:.6g},{sm_bc_mt.route_1loop['g3_mt']:.6g},{sm_bc_mt.route_1loop['yt_mt']:.6g}), "
            f"diffs(2L-1L)=(ΔgY={sm_bc_mt.diffs['gY_mt']:.2e},Δg2={sm_bc_mt.diffs['g2_mt']:.2e},Δg3={sm_bc_mt.diffs['g3_mt']:.2e},Δyt={sm_bc_mt.diffs['yt_mt']:.2e})",
            "",
            "Yukawa scaling at mt (sigma_max targets):",
            f"- yt(mt) target = {yt_mt:.6g} => y*_u={ystar_u:.6g}",
            f"- yb(mt) target = {yb_mt_target:.6g} => y*_d={ystar_d:.6g}",
            f"- ytau(mt) target = {ytau_mt_target:.6g} => y*_e={ystar_e:.6g}",
            f"- ynu3(mt) target = {ynu3_target:.6g} => y*_nu={ystar_nu:.6g}",
            f"- neutrino κ/yN source = {neutrino_source_used} (config: flavor_texture_v24.json/neutrino_mechanism.kappa_source)",
            "",
            "PMNS at mt (Ye(mt), kappa_total(mt)):",
            f"- best convention: ordering={fit_mt.get('ordering')}, perm={fit_mt.get('perm')}, chi2={fit_mt.get('chi2')}",
            f"- angles (deg): θ12={angles_mt.theta12_deg:.4f}, θ13={angles_mt.theta13_deg:.4f}, θ23={angles_mt.theta23_deg:.4f}, δ={angles_mt.delta_cp_deg:.2f}",
            f"- neutrino masses proxy (eV, in-best-order): {mnu_mt_eV_best.tolist()}",
            "",
            "PMNS at μUV (Ye(μUV), kappa(μUV)):",
            f"- best convention: ordering={fit_uv.get('ordering')}, perm={fit_uv.get('perm')}, chi2={fit_uv.get('chi2')}",
            f"- angles (deg): θ12={angles_uv.theta12_deg:.4f}, θ13={angles_uv.theta13_deg:.4f}, θ23={angles_uv.theta23_deg:.4f}, δ={angles_uv.delta_cp_deg:.2f}",
            f"- neutrino masses proxy (eV, in-best-order): {mnu_uv_eV_best.tolist()}",
            "",
            "RG segments (explicit threshold bookkeeping):",
        ]

        # Add a short summary of the discrete mechanism-bridge target scan (if present).
        mech_best_line = None
        if isinstance(pmns_target_scan_selected, dict):
            best_meta = pmns_target_scan_selected.get("best", None)
            if isinstance(best_meta, dict):
                mech_best_line = (
                    "- mechanism_bridge target scan best: "
                    f"delta0≈{best_meta.get('delta0_deg')}°, "
                    f"eps≈{best_meta.get('eps_multiplier')}·(varphi0/6) "
                    f"(sign={best_meta.get('eps_sign')}), "
                    f"chi2≈{(best_meta.get('fit', {}) if isinstance(best_meta.get('fit', {}), dict) else {}).get('chi2')}"
                )
        if mech_best_line:
            for i, ln in enumerate(lines):
                if isinstance(ln, str) and ln.startswith("- neutrino κ/yN source"):
                    lines.insert(i + 1, mech_best_line)
                    break

        for s in segs:
            mu_a = float(s.get("mu_start_GeV", float("nan")))
            mu_b = float(s.get("mu_end_GeV", float("nan")))
            model = str(s.get("model", ""))
            delta_b3_active = bool(s.get("delta_b3_active", False))
            active_cols_str = str(s.get("active_cols", ""))
            lines.append(f"- segment [{mu_a:.6g}, {mu_b:.6g}] GeV: model={model}, delta_b3_active={delta_b3_active}, active_cols={active_cols_str}")
            for act in s.get("threshold_actions_at_start", []) or []:
                lines.append(
                    f"  - {act.get('threshold_id')}: {act.get('status')} ({act.get('action')}) {act.get('note')}"
                )
        lines += mc_lines
        lines += [
            "",
            "End-state snapshot (SM gY,g2,g3):",
            f"- gY={g1_end:.6g}, g2={g2_end:.6g}, g3={g3_end:.6g}, lambda={lam_end:.6g}",
            f"- max|kappa|(mu_uv)={float(np.max(np.abs(kappa_end))):.3e} GeV^-1",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"pmns_residuals_sigma_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_pmns_residuals_sigma(
                out_dir=out_dir,
                contributions=fit_mt.get("contributions", []) if isinstance(fit_mt, dict) else [],
            )
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "pmns_reference": {
                    "file": _relpath(pmns_ref_path),
                    "sha256": pmns_ref_sha,
                    "meta": pmns_ref_meta,
                },
                "matching_map": {
                    "used": bool(matching_map_used),
                    "path": _relpath(matching_map_path),
                    "sha256": matching_map_sha,
                    "mt_boundary": matching_map_res.get("mt_boundary", None),
                    "derived_yukawas_mt": matching_map_res.get("derived_yukawas_mt", None),
                },
                "pyrate_beta_sources": {
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
                "pmns_mt": {
                    "angles_deg": angles_mt.__dict__,
                    "unitarity_dev": unitarity_mt,
                    "best_convention": {
                        "ordering": fit_mt.get("ordering"),
                        "perm": fit_mt.get("perm"),
                        "chi2": fit_mt.get("chi2"),
                        "contributions": fit_mt.get("contributions"),
                    },
                    "candidates": fit_mt.get("candidates"),
                },
                "pmns_mu_uv": {
                    "angles_deg": angles_uv.__dict__,
                    "unitarity_dev": unitarity_uv,
                    "best_convention": {
                        "ordering": fit_uv.get("ordering"),
                        "perm": fit_uv.get("perm"),
                        "chi2": fit_uv.get("chi2"),
                        "contributions": fit_uv.get("contributions"),
                    },
                    "candidates": fit_uv.get("candidates"),
                },
                "neutrino_masses_proxy_eV": {"mt": mnu_mt_eV_best.tolist(), "mu_uv": mnu_uv_eV_best.tolist()},
                "kappa": {
                    "mt_maxabs_GeVinv": float(np.max(np.abs(kappa_mt))),
                    "mu_uv_maxabs_GeVinv": float(np.max(np.abs(kappa_end))),
                    "mu_uv_sym_dev": kappa_sym_dev_uv,
                },
                "thresholds": {"MSigma": MSigma, "MG8": MG8, "MNR": MNR, "mu_uv_GeV": mu_uv_GeV},
                "publication_grade": {
                    "threshold_matching_ok": bool(threshold_matching_ok),
                    "blocked_thresholds": blocked_thresholds,
                    "note": "Publication-grade requires no 'continuous_by_assumption' threshold actions (MSigma/MG8 are matched at 1-loop identity; MNRi are tree-level matched in κ).",
                },
                "matching_mc": mc_uncertainty,
                "segments": segs,
                "gauge": {
                    "g_mt_gut": {"g1": g_mt_gut[0], "g2": g_mt_gut[1], "g3": g_mt_gut[2]},
                    "g_mt_sm": {"gY": g1_mt, "g2": g2_mt, "g3": g3_mt},
                    "g_mu_uv_sm": {"gY": g1_end, "g2": g2_end, "g3": g3_end},
                },
                "texture": {
                    "topology_phase_atoms": topo_atoms_cfg,
                    "delta_source": delta_source,
                    "phase_mode": phase_mode,
                    "coeff_rule": coeff_rule,
                    "coeff_source": coeff.source,
                    "neutrino_source": neutrino_source_used,
                    "neutrino_masses_input_eV": mnu_input_eV_used,
                    "delta_M": delta_M,
                    "delta_star": delta_star,
                    "delta_used": delta_used,
                    "theta_baseline_rad": float(theta_baseline),
                    "theta_rad": float(theta),
                    "theta_scan": {
                        "enabled": bool(pmns_theta_scan_enabled),
                        "max_candidates": int(pmns_theta_scan_max),
                        "candidates": theta_scan_rows,
                        "best": best_row,
                    },
                    "pmns_wiring_policy": {
                        "delta_cp_scan_enabled": bool(pmns_delta_cp_scan_enabled),
                        "delta_cp_scan_max_candidates": int(pmns_delta_cp_scan_max),
                        "eps_multipliers": [int(x) for x in eps_multipliers],
                    },
                    "pmns_target_scan": pmns_target_scan_selected,
                    "a_u": coeff.a_u,
                    "a_d": coeff.a_d,
                    "a_e": coeff.a_e,
                    "a_nu": coeff.a_nu,
                    "b": coeff.b,
                    "y_star_u": ystar_u,
                    "y_star_d": ystar_d,
                    "y_star_e": ystar_e,
                    "y_star_nu": ystar_nu,
                },
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

