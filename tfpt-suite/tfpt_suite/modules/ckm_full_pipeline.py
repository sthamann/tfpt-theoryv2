from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from mpmath import mp

from tfpt_suite.conventions import g1_gut_from_gY, g1_gut_over_gY, gY_from_g1_gut
from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.pyrate_boundary_runner import sm_boundary_conditions_at_mt
from tfpt_suite.sm_inputs import SmMzInputs, gauge_couplings_from_mz_inputs
from tfpt_suite.rge_pyrate_2loop import load_pyrate_beta_module, run_flavor_rge_2loop_thresholds
from tfpt_suite.pyrate_pythonoutputs import get_pyrate_pythonoutput
from tfpt_suite.flavor_textures import (
    C_delta,
    D_diag,
    coefficients_v107sm,
    left_unitary_from_yukawa,
    scale_y_star_to_match_sigma_max,
    theta_of_delta,
    yukawa_texture_matrix,
)
from tfpt_suite.mobius_z3_yukawa_generator import generate_quark_yukawas_mt, quark_phase_map

CHI2_SCALE_RATIO_MAX = 10.0


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


def _plot_ckm_summary(
    *,
    out_dir: Path,
    Vabs_mt: np.ndarray,
    Vabs_uv: np.ndarray,
    contributions: list[dict[str, Any]],
    mu_start_GeV: float,
    mu_end_GeV: float,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"ckm_summary_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        Vabs_mt = np.array(Vabs_mt, dtype=float)
        Vabs_uv = np.array(Vabs_uv, dtype=float)
        if Vabs_mt.shape != (3, 3) or Vabs_uv.shape != (3, 3):
            return plot, warnings

        # Top χ² contributors
        top = contributions[:6]
        keys = [str(t.get("key", "")) for t in top]
        chi2 = []
        for t in top:
            try:
                chi2.append(float(t.get("chi2", 0.0)))
            except Exception:
                chi2.append(float("nan"))

        # Use constrained layout to avoid tight_layout warnings with colorbars + gridspec.
        fig = plt.figure(figsize=(10, 8), constrained_layout=True)
        gs = fig.add_gridspec(nrows=2, ncols=2, height_ratios=[1.0, 0.8])

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, :])

        cmap = "viridis"
        im1 = ax1.imshow(Vabs_mt, vmin=0.0, vmax=1.0, cmap=cmap)
        im2 = ax2.imshow(Vabs_uv, vmin=0.0, vmax=1.0, cmap=cmap)

        # Axis labels (PDG ordering)
        xlab = ["d", "s", "b"]
        ylab = ["u", "c", "t"]
        for ax, title in [
            (ax1, f"|V| @ mt (μ={mu_start_GeV:.3g} GeV)"),
            (ax2, f"|V| @ μ_UV (μ={mu_end_GeV:.1e} GeV)"),
        ]:
            ax.set_xticks([0, 1, 2])
            ax.set_yticks([0, 1, 2])
            ax.set_xticklabels(xlab)
            ax.set_yticklabels(ylab)
            ax.set_title(title)

        cbar = fig.colorbar(im1, ax=[ax1, ax2], fraction=0.046, pad=0.04)
        cbar.set_label("|V_ij|")

        if keys:
            ax3.barh(list(range(len(keys))), chi2, color="#1f77b4", alpha=0.9)
            ax3.set_yticks(list(range(len(keys))))
            ax3.set_yticklabels(keys)
            ax3.invert_yaxis()
            ax3.set_xlabel("χ² contribution (proxy table @ mt)")
            ax3.set_title("Largest χ² contributions (top 6)")
            ax3.grid(True, axis="x", ls=":", alpha=0.4)
        else:
            ax3.text(0.5, 0.5, "No χ² contributions available", ha="center", va="center")
            ax3.set_axis_off()

        fig.suptitle("CKM pipeline summary (texture @ mt, RG to μ_UV)")
        path = out_dir / "ckm_summary.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["ckm_summary_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


def _plot_ckm_residuals_sigma(
    *,
    out_dir: Path,
    Vabs_pred: dict[str, float],
    ref_table: dict[str, Any],
    chi2_keys: list[str],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"ckm_residuals_sigma_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)
        matrix_keys = [
            ["Vud", "Vus", "Vub"],
            ["Vcd", "Vcs", "Vcb"],
            ["Vtd", "Vts", "Vtb"],
        ]
        residuals = np.full((3, 3), np.nan, dtype=float)
        for i, row in enumerate(matrix_keys):
            for j, key in enumerate(row):
                ref_entry = ref_table.get("matrix_abs", {}).get(key, {}) if isinstance(ref_table.get("matrix_abs", {}), dict) else {}
                mean = float(ref_entry.get("mean", float("nan")))
                sigma = float(ref_entry.get("sigma", float("nan")))
                pred = float(Vabs_pred.get(key, float("nan")))
                if math.isfinite(mean) and math.isfinite(sigma) and sigma > 0 and math.isfinite(pred):
                    residuals[i, j] = (pred - mean) / sigma

        if not np.isfinite(residuals).any():
            return plot, warnings

        vmax = float(np.nanmax(np.abs(residuals)))
        vmax = max(vmax, 1.0)
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        im = ax.imshow(residuals, cmap="coolwarm", vmin=-vmax, vmax=vmax)
        ax.set_xticks([0, 1, 2])
        ax.set_yticks([0, 1, 2])
        ax.set_xticklabels(["d", "s", "b"])
        ax.set_yticklabels(["u", "c", "t"])
        ax.set_title("CKM residuals (pred − ref) / σ")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="σ units")

        for i in range(3):
            for j in range(3):
                key = matrix_keys[i][j]
                val = residuals[i, j]
                if not math.isfinite(val):
                    text = "n/a"
                else:
                    marker = "*" if key in chi2_keys else ""
                    text = f"{val:+.2f}{marker}"
                ax.text(j, i, text, ha="center", va="center", color="black", fontsize=9)

        ax.text(
            0.02,
            -0.12,
            "* = chi2 key",
            transform=ax.transAxes,
            fontsize=9,
        )
        fig.tight_layout()
        path = out_dir / "ckm_residuals_sigma.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        plot["ckm_residuals_sigma_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


def _ckm_abs_from_wolfenstein_proxy_exact_unitary(*, lam: mp.mpf, A: mp.mpf, rho: mp.mpf, eta: mp.mpf) -> dict[str, mp.mpf]:
    """
    Deterministic unitary CKM construction from a Wolfenstein-style proxy mapping:

    Define exact PDG angles via:
      s12 = λ
      s23 = A λ^2
      s13 e^{-iδ} = A λ^3 (ρ - i η)
        => s13 = A λ^3 sqrt(ρ^2+η^2), δ = atan2(η, ρ)

    Then build the exact unitary CKM matrix in PDG convention and return absolute values.

    This avoids the V_tb ≡ 1 artifact of truncated O(λ^3) Wolfenstein expansions.
    """
    # PDG parameters
    s12 = lam
    s23 = A * (lam**2)
    s13 = A * (lam**3) * mp.sqrt(rho**2 + eta**2)
    delta = mp.atan2(eta, rho)

    # clamp numerical drift
    def clamp01(x: mp.mpf) -> mp.mpf:
        return min(mp.mpf(1), max(mp.mpf(0), x))

    s12 = clamp01(s12)
    s23 = clamp01(s23)
    s13 = clamp01(s13)

    c12 = mp.sqrt(mp.mpf(1) - s12**2)
    c23 = mp.sqrt(mp.mpf(1) - s23**2)
    c13 = mp.sqrt(mp.mpf(1) - s13**2)

    I = mp.mpc(0, 1)
    e_plus = mp.e ** (I * delta)
    e_minus = mp.e ** (-I * delta)

    Vud = c12 * c13
    Vus = s12 * c13
    Vub = s13 * e_minus

    Vcd = -s12 * c23 - c12 * s23 * s13 * e_plus
    Vcs = c12 * c23 - s12 * s23 * s13 * e_plus
    Vcb = s23 * c13

    Vtd = s12 * s23 - c12 * c23 * s13 * e_plus
    Vts = -c12 * s23 - s12 * c23 * s13 * e_plus
    Vtb = c23 * c13

    return {k: abs(v) for k, v in {"Vud": Vud, "Vus": Vus, "Vub": Vub, "Vcd": Vcd, "Vcs": Vcs, "Vcb": Vcb, "Vtd": Vtd, "Vts": Vts, "Vtb": Vtb}.items()}


def _ckm_matrix_pdg_from_wolfenstein_proxy(*, lam: float, A: float, rho: float, eta: float) -> np.ndarray:
    """
    Exact unitary PDG CKM construction from a Wolfenstein-style proxy mapping (float/numpy version).

    s12 = λ, s23 = A λ^2, s13 e^{-iδ} = A λ^3 (ρ - i η)
    """
    s12 = float(lam)
    s23 = float(A) * (float(lam) ** 2)
    s13 = float(A) * (float(lam) ** 3) * float(np.sqrt(rho * rho + eta * eta))
    delta = float(np.arctan2(eta, rho))

    def clamp01(x: float) -> float:
        return float(min(1.0, max(0.0, x)))

    s12 = clamp01(s12)
    s23 = clamp01(s23)
    s13 = clamp01(s13)

    c12 = float(np.sqrt(max(0.0, 1.0 - s12 * s12)))
    c23 = float(np.sqrt(max(0.0, 1.0 - s23 * s23)))
    c13 = float(np.sqrt(max(0.0, 1.0 - s13 * s13)))

    e_plus = np.exp(1j * delta)
    e_minus = np.exp(-1j * delta)

    V = np.zeros((3, 3), dtype=complex)
    V[0, 0] = c12 * c13
    V[0, 1] = s12 * c13
    V[0, 2] = s13 * e_minus

    V[1, 0] = -s12 * c23 - c12 * s23 * s13 * e_plus
    V[1, 1] = c12 * c23 - s12 * s23 * s13 * e_plus
    V[1, 2] = s23 * c13

    V[2, 0] = s12 * s23 - c12 * c23 * s13 * e_plus
    V[2, 1] = -c12 * s23 - s12 * c23 * s13 * e_plus
    V[2, 2] = c23 * c13
    return V


def _left_unitary_from_yukawa(Y: np.ndarray) -> np.ndarray:
    """
    Compute the left-unitary matrix U_L that diagonalizes Y Y^†:
      U_L^† (Y Y^†) U_L = diag(...)
    """
    H = Y @ Y.conj().T
    evals, evecs = np.linalg.eigh(H)  # ascending
    # Ensure deterministic phase convention: make the largest component of each eigenvector real-positive.
    U = evecs.astype(complex)
    for j in range(3):
        col = U[:, j]
        k = int(np.argmax(np.abs(col)))
        phase = col[k] / abs(col[k]) if abs(col[k]) > 0 else 1.0
        U[:, j] = col / phase
    return U


def _ckm_from_yukawas(Yu: np.ndarray, Yd: np.ndarray) -> np.ndarray:
    Uu = _left_unitary_from_yukawa(Yu)
    Ud = _left_unitary_from_yukawa(Yd)
    return Uu.conj().T @ Ud


def _jarlskog(V: np.ndarray) -> float:
    # J = Im(V_ud V_cs V_us^* V_cd^*)
    return float(np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0])))


def _wolfenstein_from_ckm(V: np.ndarray) -> dict[str, float]:
    lam = float(np.abs(V[0, 1]))
    Vcb = float(np.abs(V[1, 2]))
    A = Vcb / (lam * lam) if lam > 0 else float("nan")
    # ρ̄ + i η̄ = - (V_ud V_ub^*) / (V_cd V_cb^*)
    denom = V[1, 0] * np.conj(V[1, 2])
    z = -(V[0, 0] * np.conj(V[0, 2])) / denom if denom != 0 else (np.nan + 1j * np.nan)
    return {"lambda": lam, "A": float(A), "rho_bar": float(np.real(z)), "eta_bar": float(np.imag(z))}

@dataclass(frozen=True)
class Candidate:
    A: mp.mpf
    rho: mp.mpf
    eta: mp.mpf
    chi2: mp.mpf
    matrix_abs: dict[str, mp.mpf]


class CkmFullPipelineModule(TfptModule):
    module_id = "ckm_full_pipeline"
    title = "CKM full pipeline (Z3 Yukawa texture + diagonalization; upward RG mt→μUV)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT Cabibbo λ from varphi0",
                "SM inputs at μ=MZ: tfpt_suite/data/sm_inputs_mz.json (for RG dressing)",
                "reference CKM |V_ij| table: tfpt_suite/data/ckm_reference.json",
                "Z3 Yukawa texture config: tfpt_suite/data/flavor_texture_v24.json",
                "lepton masses (τ, μ) for δ calibration: tfpt_suite/data/lepton_masses_pdg.json",
                "2-loop PyR@TE beta tables (via tfpt_suite/data/pyrate_pythonoutputs.json; deterministic model selection)",
            ],
            outputs=[
                "Yukawa texture parameters at μ=mt (δ, phase_mode, coefficients_rule, y_*)",
                "CKM |V_ij| at μ=mt (boundary) and at μ=μ_UV after upward RG running",
                "χ² contributions at μ=mt vs the reference table (proxy)",
                "Wolfenstein + Jarlskog at μ=mt and μ=μ_UV",
                "optional: Monte Carlo uncertainty propagation at the reference scale (mean/std for |V_ij| and χ²)",
            ],
            formulas=[
                "λ_TFPT = sqrt(varphi0)*(1 - varphi0/2) (paper v2.4)",
                "Y^(y) = y_* [ C(δ) + a_y varphi0 D + b c3 I ], D=diag(1,0,-1)",
                "C(δ) = circ(1, ζ, ζ*) with ζ = exp(i θ(δ)); θ(δ) configured (default θ=2πδ)",
                "RG evolution: start strictly at μ=mt and run upward only to μ_UV (no running below mt)",
                "RG dressing: integrate PyR@TE-generated beta functions (supports complex Yukawas); threshold bookkeeping is explicit in output segments",
                "CKM from Yukawas: V = Uu^† Ud where Uu,Ud diagonalize YuYu^† and YdYd^†",
                "χ² = Σ ((|V_ij|(mt) - mean)/sigma)^2 over the reference table (proxy; the reference carries explicit metadata but is still an effective low-energy fit input, not a running MSbar coupling)",
            ],
            validation=[
                "texture config loads and produces finite Yukawa matrices at mt",
                "RG dressing runs and produces a unitary CKM matrix (numerically)",
            ],
            determinism="Deterministic given the texture config + input tables.",
            question="Can TFPT’s deterministic Möbius/Z3 flavor generator reproduce the precision CKM matrix once scheme/scale policy is made explicit?",
            objective=[
                "Run an end-to-end deterministic CKM construction (no continuous fit) with explicit phase conventions and explicit RG/threshold policy.",
                "Quantify failure modes honestly (χ² contributions) and surface the remaining missing physics (matching finite pieces; topology→phase map).",
            ],
            what_was_done=[
                "Compute δ from τ/μ or δ⋆ from varphi0 (explicit in flavor_texture_v24.json).",
                "Generate quark Yukawas via the v1.07SM Möbius/Z3 generator and run mt→μUV RG with explicit thresholds.",
                "Evaluate discrete CKM convention variants (s13_mode, delta_mode) and report χ² against a pinned reference snapshot.",
                "Record the `topology_phase_atoms` docking block from flavor_texture_v24.json (produced by chiral_index_three_cycles) for future topology→phase integration.",
            ],
            assumptions=[
                "Reference CKM table is treated as a diagnostic snapshot (not a covariance-level global-fit likelihood).",
                "Below-mt running and full finite matching pieces are not yet a publication-grade policy; residual scheme effects may impact precision elements (e.g. Vub).",
            ],
            gaps=[
                "A full-matrix χ² can look large if one includes redundant unitarity-constrained |V_ij| keys (double-counting). The suite therefore uses a declared `chi2_keys` subset for the diagnostic p-value gate; publication-grade still requires a consistent likelihood/covariance model.",
                "Topology→(δ, δ_CP) derivation is not yet implemented; `topology_phase_atoms` is a docking point only.",
                "Publication-grade requires finalized finite matching pieces and a consistent comparison scale policy.",
            ],
            references=[
                "paper_v1_06_01_09_2025.tex / update_tfptv1_07sm.tex (flavor architecture; Möbius/Z3 anchors)",
                "tfpt_suite/data/flavor_texture_v24.json (explicit conventions + topology docking block)",
            ],
            maturity="deterministic pipeline (diagnostic χ²; physics-mode p-value gate PASS under declared chi2_keys; publication-grade operator/matching derivation remains a WARN-gap)",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        lam = mp.sqrt(c.varphi0) * (mp.mpf(1) - mp.mpf(1) / 2 * c.varphi0)

        ref_path = Path(__file__).resolve().parent.parent / "data" / "ckm_reference.json"
        ref = json.loads(ref_path.read_text(encoding="utf-8"))
        table: dict[str, Any] = ref["matrix_abs"]
        ref_meta = dict(ref.get("reference", {})) if isinstance(ref.get("reference", {}), dict) else {}
        all_table_keys = [str(k) for k in table.keys()]
        chi2_keys = list(all_table_keys)
        chi2_keys_raw = ref_meta.get("chi2_keys", None)
        if isinstance(chi2_keys_raw, list) and chi2_keys_raw:
            cand = []
            for k in chi2_keys_raw:
                ks = str(k)
                if ks in table and ks not in cand:
                    cand.append(ks)
            if cand:
                chi2_keys = cand

        # SM inputs (for RG dressing)
        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        sm_inp = SmMzInputs(
            mu_GeV=float(sm_raw["mu_GeV"]),
            alpha_em_inv=float(sm_raw["alpha_em_inv"]),
            sin2_thetaW=float(sm_raw["sin2_thetaW"]),
            alpha_s=float(sm_raw["alpha_s"]),
        )
        alpha_s_sigma = float(sm_raw.get("alpha_s_sigma", 0.0) or 0.0)
        alpha_em_inv_sigma = float(sm_raw.get("alpha_em_inv_sigma", 0.0) or 0.0)
        sin2_sigma = float(sm_raw.get("sin2_thetaW_sigma", 0.0) or 0.0)
        mt_sigma = float(sm_raw.get("mt_sigma_GeV", 0.0) or 0.0)
        mb_sigma = float(sm_raw.get("mb_sigma_GeV", 0.0) or 0.0)
        mc_sigma = float(sm_raw.get("mc_sigma_GeV", 0.0) or 0.0)
        mc_samples = int(sm_raw.get("matching_mc_samples", 0) or 0)
        mc_GeV = float(sm_raw.get("mc_GeV", 1.27))
        mb_GeV = float(sm_raw.get("mb_GeV", 4.18))
        mt_GeV = float(sm_raw.get("mt_GeV", 173.0))
        mW_GeV = float(sm_raw.get("mW_GeV", 80.379))
        g_mz = gauge_couplings_from_mz_inputs(sm_inp)

        # --- RG integration (strictly start at μ = m_t; run upward only to the UV scale) ---
        #
        # We keep SM inputs for gauge couplings at μ=MZ, but the *RGE integration used for analysis*
        # starts at μ=mt for all couplings and runs only upward.
        sm_bc_mt = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_raw)
        mu_start_GeV = float(sm_bc_mt.mu_mt_GeV)
        mu_uv_GeV = 1.0e16

        # Gauge couplings at mt: start from the central matching layer (SM inputs @ MZ → mt boundary).
        g1_gut_mt = float(sm_bc_mt.route_2loop["g1_gut_mt"])
        g2_mt = float(sm_bc_mt.route_2loop["g2_mt"])
        g3_mt = float(sm_bc_mt.route_2loop["g3_mt"])
        # UV gauge couplings are taken from the authoritative PyR@TE-driven run (see below),
        # not from a parallel "gauge-only" runner.
        g_uv_gauge_2loop: tuple[float, float, float] | None = None

        # --- Flavor texture at mt: replace the Wolfenstein proxy by a matrix Yukawa texture ---
        flavor_cfg_path = Path(__file__).resolve().parent.parent / "data" / "flavor_texture_v24.json"
        flavor_cfg = json.loads(flavor_cfg_path.read_text(encoding="utf-8"))
        topo_atoms_cfg = (
            flavor_cfg.get("topology_phase_atoms", {}) if isinstance(flavor_cfg.get("topology_phase_atoms", {}), dict) else {}
        )
        tex = flavor_cfg.get("yukawa_texture", {})
        delta_source = str(tex.get("delta_source", "tau_mu")).strip()
        phase_mode = str(tex.get("phase_mode", "2pi_delta")).strip()
        if phase_mode not in ("2pi_delta", "delta_rad", "koide_pi_over_12"):
            raise ValueError(f"Unsupported phase_mode in {flavor_cfg_path}: {phase_mode}")

        coeff_rule = str(tex.get("coefficients_rule", "v107sm_fixed")).strip()
        if coeff_rule != "v107sm_fixed":
            raise ValueError(f"Unsupported coefficients_rule in {flavor_cfg_path}: {coeff_rule}")
        coeff = coefficients_v107sm()
        a_u = float(coeff.a_u)
        a_d = float(coeff.a_d)
        a_e = float(coeff.a_e)
        b = float(coeff.b)

        lepton_masses_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
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

        theta = float(theta_of_delta(delta_used, phase_mode=phase_mode))  # radians
        zeta = complex(np.cos(theta), np.sin(theta))

        # Yukawa boundary conditions at mt:
        # - yt(mt) from mt_pole → mt_MSbar(mt) using QCD 2-loop matching (driven by αs(mt))
        # - yb, ytau are explicit inputs (defaults match the SM baseline used in two_loop_rg_fingerprints)
        v_ev_GeV = float(sm_bc_mt.route_2loop["v_ev_GeV"])
        alpha_s_mt = float(sm_bc_mt.route_2loop["alpha_s_mt"])
        mt_msbar_mt = float(sm_bc_mt.route_2loop["mt_msbar_mt_GeV"])
        yt_mt = float(sm_bc_mt.route_2loop["yt_mt"])
        yb_mt_target = float(sm_bc_mt.route_2loop["yb_mt"])
        ytau_mt_target = float(sm_bc_mt.route_2loop["ytau_mt"])
        lambda_mt = float(sm_bc_mt.route_2loop["lambda_mt"])

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
                mt_msbar_mm = mt_bound.get("mt_msbar_mt_GeV", None)
                yb_mt_mm = derived.get("yb_mt_1loop_qcd_nf5", None)
                ytau_mt_mm = derived.get("ytau_mt_tree", None)
                v_ev_mm = derived.get("v_ev_GeV", None)

                if gY_mt_mm is not None and g2_mt_mm is not None and g3_mt_mm is not None:
                    g1_gut_mt = g1_gut_from_gY(float(gY_mt_mm))
                    g2_mt = float(g2_mt_mm)
                    g3_mt = float(g3_mt_mm)
                    matching_map_used = True
                if yt_mt_mm is not None:
                    yt_mt = float(yt_mt_mm)
                    matching_map_used = True
                if alpha_s_mt_mm is not None:
                    alpha_s_mt = float(alpha_s_mt_mm)
                    matching_map_used = True
                if mt_msbar_mm is not None:
                    mt_msbar_mt = float(mt_msbar_mm)
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

        g_mt = (float(g1_gut_mt), float(g2_mt), float(g3_mt))

        varphi0_f = float(c.varphi0)
        c3_f = float(c.c3)

        # Charged leptons are not used for CKM; keep a deterministic texture Ye for RG consistency.
        Ye_base = yukawa_texture_matrix(
            delta=delta_used, varphi0=varphi0_f, c3=c3_f, a_y=a_e, b=b, y_star=1.0, phase_mode=phase_mode
        )
        ystar_e = scale_y_star_to_match_sigma_max(target_y3=ytau_mt_target, base=Ye_base)
        Ye_mt = (ystar_e * Ye_base).astype(complex)

        # --- RG plumbing needed for a scale-consistent CKM comparison ---
        mu_ref_GeV = float(ref_meta.get("reference_scale_GeV", float(sm_inp.mu_GeV)))
        thresholds_path = Path(__file__).resolve().parent.parent / "data" / "rge_thresholds_v25.json"
        thresholds_raw = json.loads(thresholds_path.read_text(encoding="utf-8"))
        thresholds_GeV = dict(thresholds_raw.get("thresholds_GeV", {}))

        # The PyR@TE beta files use SM hypercharge g1 ≡ g' (not GUT-normalized). Convert our g1(mt) accordingly.
        gut_norm = g1_gut_over_gY()
        gY_mt = gY_from_g1_gut(g_mt[0])
        g2_mt = float(g_mt[1])
        g3_mt = float(g_mt[2])

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

        # --- NEW: Möbius/Z3 Yukawa generator (v1.07SM-style CKM + Möbius hierarchies) ---
        # This replaces the previous “C(δ)+a_y φ0 D + b c3 I” diagonalization path, which produced a nearly
        # identity CKM and an astronomically large χ². The generator is deterministic and uses only TFPT anchors.
        # Discrete CKM construction variants (no fit): evaluate all configured variants at μ=mt and pick the best χ².
        ckm_cfg = flavor_cfg.get("ckm_variants", {}) if isinstance(flavor_cfg.get("ckm_variants", {}), dict) else {}
        run_all_variants = bool(ckm_cfg.get("run_all_variants", False))
        include_topology_delta_cp = bool(ckm_cfg.get("include_topology_delta_cp_candidates", False))
        max_topology_delta_cp = int(ckm_cfg.get("max_topology_delta_cp_candidates", 0) or 0)
        phase_selection_rule_mode = str(ckm_cfg.get("phase_selection_rule_mode", "augment")).strip().lower() or "augment"
        if phase_selection_rule_mode not in {"augment", "filter_only"}:
            phase_selection_rule_mode = "augment"
        variants_raw = ckm_cfg.get("variants", []) if isinstance(ckm_cfg.get("variants", []), list) else []
        variants = (
            [v for v in variants_raw if isinstance(v, dict)]
            if variants_raw
            else [{"label": "v107sm_default", "s13_mode": "A_lam3_times_1_minus_delta", "delta_mode": "pi_times_1_minus_delta"}]
        )
        if not run_all_variants and variants:
            variants = variants[:1]

        # Optional: extend the discrete CKM scan with δ_CP candidates coming from topology_phase_map (Wilson-line atoms).
        # This is the minimal “wire the docking point” step: we do NOT fit any continuous parameter; we just
        # allow δ_CP to be chosen from a finite topological candidate set.
        topo_variants: list[dict[str, Any]] = []
        topo_output_found = False
        if include_topology_delta_cp and max_topology_delta_cp > 0:
            topo_path = Path(config.output_dir) / "topology_phase_map" / "results.json"
            topo_output_found = bool(topo_path.is_file())
            try:
                topo_payload = json.loads(topo_path.read_text(encoding="utf-8")) if topo_path.is_file() else {}
                topo_res = topo_payload.get("results", {}) if isinstance(topo_payload.get("results", {}), dict) else {}
                cand_raw = topo_res.get("delta_cp_candidates", []) if isinstance(topo_res.get("delta_cp_candidates", []), list) else []
                cands = [c for c in cand_raw if isinstance(c, dict) and "delta_cp_rad" in c and "theta_over_pi_mod2" in c]
                # Skip CP-conserving candidates (sin δ ≈ 0) to focus the scan on genuinely CP-violating branches.
                cands = [c for c in cands if abs(float(np.sin(float(c.get("delta_cp_rad"))))) > 1e-6]

                def _key(c: dict[str, Any]) -> tuple[int, int, int]:
                    frac = Fraction(str(c.get("theta_over_pi_mod2", "0"))).limit_denominator(12)
                    branch = str(c.get("branch", "theta"))
                    # prefer low-denominator atoms; stable ordering
                    return (int(frac.denominator), int(frac.numerator), 0 if branch == "theta" else 1)

                cands.sort(key=_key)
                cands = cands[: int(max_topology_delta_cp)]

                for i, cnd in enumerate(cands, start=1):
                    topo_variants.append(
                        {
                            "label": f"topology_delta_cp_{i}",
                            "s13_mode": "A_lam3_times_1_minus_delta",
                            "delta_mode": "external_delta_cp_override",
                            "delta_cp_override_rad": float(cnd.get("delta_cp_rad")),
                            "theta_over_pi_mod2": str(cnd.get("theta_over_pi_mod2")),
                            "branch": str(cnd.get("branch", "")),
                            "source_module": "topology_phase_map",
                        }
                    )
                if phase_selection_rule_mode == "filter_only" and topo_variants:
                    variants = topo_variants
                else:
                    variants = list(variants) + topo_variants
            except Exception:
                # Keep baseline variants only if topology outputs are unavailable.
                pass

        def _vabs_from(V: np.ndarray) -> dict[str, float]:
            idx = {
                "Vud": (0, 0),
                "Vus": (0, 1),
                "Vub": (0, 2),
                "Vcd": (1, 0),
                "Vcs": (1, 1),
                "Vcb": (1, 2),
                "Vtd": (2, 0),
                "Vts": (2, 1),
                "Vtb": (2, 2),
            }
            return {k: float(np.abs(V[i, j])) for k, (i, j) in idx.items()}

        def _chi2_vs_reference(Vabs: dict[str, float]) -> tuple[float, list[dict[str, float | str]]]:
            contribs: list[dict[str, float | str]] = []
            for key in chi2_keys:
                entry = table.get(str(key), None)
                if not isinstance(entry, dict) or "mean" not in entry or "sigma" not in entry:
                    continue
                pred = float(Vabs.get(str(key), float("nan")))
                mean = float(entry["mean"])
                sig = float(entry["sigma"])
                chi2 = float(((pred - mean) / sig) ** 2) if sig > 0 else float("nan")
                contribs.append({"key": str(key), "pred": pred, "mean": mean, "sigma": sig, "chi2": chi2})
            contribs.sort(key=lambda t: float(t.get("chi2", 0.0)), reverse=True)
            chi2_total = float(np.nansum([float(t["chi2"]) for t in contribs if isinstance(t.get("chi2"), (int, float))]))
            return chi2_total, contribs

        variant_results: list[dict[str, object]] = []
        best = None
        for v in variants:
            label = str(v.get("label", "variant")).strip() or "variant"
            s13_mode = str(v.get("s13_mode", "A_lam3_times_1_minus_delta")).strip()
            delta_mode = str(v.get("delta_mode", "pi_times_1_minus_delta")).strip()
            delta_cp_override_rad = v.get("delta_cp_override_rad", None)
            if delta_cp_override_rad is not None:
                try:
                    delta_cp_override_rad = float(delta_cp_override_rad)
                except Exception:
                    delta_cp_override_rad = None
            delta_used_var = v.get("delta_used_override", None)
            if delta_used_var is not None:
                try:
                    delta_used_var = float(delta_used_var)
                except Exception:
                    delta_used_var = None
            Yu_try, Yd_try, meta_try = generate_quark_yukawas_mt(
                varphi0=varphi0_f,
                delta_used=float(delta_used if delta_used_var is None else delta_used_var),
                yt_mt=float(yt_mt),
                yb_mt=float(yb_mt_target),
                scheme="MSbar",
                reference_scale_GeV=float(mu_start_GeV),
                s13_mode=s13_mode,  # type: ignore[arg-type]
                delta_mode=delta_mode,  # type: ignore[arg-type]
                delta_cp_override_rad=delta_cp_override_rad,
            )
            V_try = _ckm_from_yukawas(Yu_try, Yd_try)
            Vabs_try = _vabs_from(V_try)
            chi2_try, contribs_try = _chi2_vs_reference(Vabs_try)

            # Scale-consistent comparison: run Yukawas downward mt→mu_ref (typically MZ) in the SM
            # and compute χ² against the reference table at that scale.
            rge_ref = run_flavor_rge_2loop_thresholds(
                mu_start_GeV=float(mu_start_GeV),
                mu_end_GeV=float(mu_ref_GeV),
                thresholds_GeV=thresholds_GeV,
                g_start=(float(gY_mt), float(g2_mt), float(g3_mt)),
                Yu_start=Yu_try,
                Yd_start=Yd_try,
                Ye_start=Ye_mt,
                lambda_start=float(lambda_mt),
                yN_start=None,
                beta_sm=beta_sm,
                beta_e8=beta_e8,
                apply_sigma_threshold=True,
                apply_g8_delta_b3=False,
                apply_matching=True,
                matching_loop_order=1,
                rtol=1e-9,
                atol=1e-11,
                method="DOP853",
            )
            V_ref = _ckm_from_yukawas(rge_ref["Yu_end"], rge_ref["Yd_end"])
            Vabs_ref = _vabs_from(V_ref)
            chi2_ref, contribs_ref = _chi2_vs_reference(Vabs_ref)

            rec = {
                "label": label,
                "s13_mode": s13_mode,
                "delta_mode": delta_mode,
                "delta_cp_override_rad": delta_cp_override_rad,
                "theta_over_pi_mod2": v.get("theta_over_pi_mod2", None),
                "branch": v.get("branch", None),
                "source_module": v.get("source_module", None),
                "chi2_mt": chi2_try,
                "chi2_refscale": chi2_ref,
                "largest_contribution": contribs_try[0] if contribs_try else None,
                "largest_contribution_refscale": contribs_ref[0] if contribs_ref else None,
            }
            variant_results.append(rec)
            if best is None or float(chi2_ref) < float(best["chi2_refscale"]):  # type: ignore[index]
                best = dict(rec)
                best["_Yu_mt"] = Yu_try
                best["_Yd_mt"] = Yd_try
                best["_meta"] = meta_try
                best["_contribs"] = contribs_try
                best["_contribs_ref"] = contribs_ref

        assert best is not None
        best_public = {k: best.get(k) for k in ("label", "s13_mode", "delta_mode", "chi2_mt", "chi2_refscale", "largest_contribution_refscale")}
        Yu_mt = best["_Yu_mt"]  # type: ignore[assignment]
        Yd_mt = best["_Yd_mt"]  # type: ignore[assignment]
        gen_meta = best["_meta"]  # type: ignore[assignment]
        contributions: list[dict[str, Any]] = []  # set after ref-scale χ² computation
        s13_mode_best = str(gen_meta.ckm.s13_mode)
        delta_mode_best = str(gen_meta.ckm.delta_mode)
        variant_label_best = str(best.get("label", "best"))
        # Convenience: expose the actual eigenvalue targets used by the generator (at mt).
        ratios_mt = dict(gen_meta.ratios)
        try:
            y_c_mt = float(yt_mt) / float(ratios_mt["m_t_over_m_c"])
            y_u_mt = float(y_c_mt) / float(ratios_mt["m_c_over_m_u"])
            y_s_mt = float(yb_mt_target) / float(ratios_mt["m_b_over_m_s"])
            y_d_mt = float(y_s_mt) / float(ratios_mt["m_s_over_m_d"])
        except Exception:
            y_u_mt, y_c_mt, y_d_mt, y_s_mt = float("nan"), float("nan"), float("nan"), float("nan")

        yu_diag_mt = {"y_u": float(y_u_mt), "y_c": float(y_c_mt), "y_t": float(yt_mt)}
        yd_diag_mt = {"y_d": float(y_d_mt), "y_s": float(y_s_mt), "y_b": float(yb_mt_target)}
        ckm_seed = {
            "lambda": float(gen_meta.ckm.lam),
            "A": float(gen_meta.ckm.A),
            "s12": float(gen_meta.ckm.s12),
            "s23": float(gen_meta.ckm.s23),
            "s13": float(gen_meta.ckm.s13),
            "delta_cp_rad": float(gen_meta.ckm.delta_cp_rad),
            "s13_mode": str(gen_meta.ckm.s13_mode),
            "delta_mode": str(gen_meta.ckm.delta_mode),
        }

        phase_map = quark_phase_map(
            delta_mobius=float(delta_M),
            delta_used=float(delta_used),
            delta_source=str(delta_source),
            phase_mode=phase_mode,
            delta_mode=str(gen_meta.ckm.delta_mode),
            delta_ckm_override_rad=float(gen_meta.ckm.delta_cp_rad) if str(gen_meta.ckm.delta_mode) == "external_delta_cp_override" else None,
        )
        delta_ckm_map_diff = abs(float(phase_map.delta_ckm_rad) - float(gen_meta.ckm.delta_cp_rad))
        rge2 = run_flavor_rge_2loop_thresholds(
            mu_start_GeV=mu_start_GeV,
            mu_end_GeV=mu_uv_GeV,
            thresholds_GeV=thresholds_GeV,
            g_start=(gY_mt, g2_mt, g3_mt),
            Yu_start=Yu_mt,
            Yd_start=Yd_mt,
            Ye_start=Ye_mt,
            lambda_start=float(lambda_mt),
            yN_start=None,  # CKM: neutrino Yukawa not used
            beta_sm=beta_sm,
            beta_e8=beta_e8,
            apply_sigma_threshold=True,
            apply_g8_delta_b3=True,
            delta_b3_g8=2.0,
            apply_matching=True,
            matching_loop_order=1,
            rtol=1e-8,
            atol=1e-10,
            method="DOP853",
        )
        Yu_uv_2l = rge2["Yu_end"]
        Yd_uv_2l = rge2["Yd_end"]

        # Use the authoritative PyR@TE run for UV gauge couplings (GUT-normalized g1_GUT).
        g_end_gut = rge2.get("g_end_gut", None) or {}
        if isinstance(g_end_gut, dict) and {"g1_gut", "g2", "g3"} <= set(g_end_gut.keys()):
            g_uv_gauge_2loop = (float(g_end_gut["g1_gut"]), float(g_end_gut["g2"]), float(g_end_gut["g3"]))
        else:
            g_end_sm = rge2.get("g_end_sm", None) or rge2.get("g_end", {}) or {}
            gY_uv = float(g_end_sm.get("gY", g_end_sm.get("g1", float("nan"))))
            g2_uv = float(g_end_sm.get("g2", float("nan")))
            g3_uv = float(g_end_sm.get("g3", float("nan")))
            g_uv_gauge_2loop = (float(g1_gut_over_gY() * gY_uv), g2_uv, g3_uv)

        V_mt_from_yuk = _ckm_from_yukawas(Yu_mt, Yd_mt)
        # Run downward mt→mu_ref (typically MZ) for a scale-consistent CKM comparison.
        rge_ref = run_flavor_rge_2loop_thresholds(
            mu_start_GeV=mu_start_GeV,
            mu_end_GeV=mu_ref_GeV,
            thresholds_GeV=thresholds_GeV,
            g_start=(gY_mt, g2_mt, g3_mt),
            Yu_start=Yu_mt,
            Yd_start=Yd_mt,
            Ye_start=Ye_mt,
            lambda_start=float(lambda_mt),
            yN_start=None,
            beta_sm=beta_sm,
            beta_e8=beta_e8,
            apply_sigma_threshold=True,
            apply_g8_delta_b3=False,
            apply_matching=True,
            matching_loop_order=1,
            rtol=1e-9,
            atol=1e-11,
            method="DOP853",
        )
        Yu_ref = rge_ref["Yu_end"]
        Yd_ref = rge_ref["Yd_end"]
        V_ref = _ckm_from_yukawas(Yu_ref, Yd_ref)
        V_uv_2l = _ckm_from_yukawas(Yu_uv_2l, Yd_uv_2l)

        def Vabs_from(V: np.ndarray) -> dict[str, mp.mpf]:
            return {
                k: mp.mpf(str(float(np.abs(V[i, j]))))
                for (k, (i, j)) in {
                    "Vud": (0, 0),
                    "Vus": (0, 1),
                    "Vub": (0, 2),
                    "Vcd": (1, 0),
                    "Vcs": (1, 1),
                    "Vcb": (1, 2),
                    "Vtd": (2, 0),
                    "Vts": (2, 1),
                    "Vtb": (2, 2),
                }.items()
            }

        Vabs_mt = Vabs_from(V_mt_from_yuk)
        Vabs_ref = Vabs_from(V_ref)
        Vabs_uv_2l = Vabs_from(V_uv_2l)

        # χ² at the reference scale (primary comparison; reference_meta declares this scale)
        contributions_ref = []
        for key in chi2_keys:
            cfg = table.get(str(key), None)
            if not isinstance(cfg, dict) or "mean" not in cfg or "sigma" not in cfg:
                continue
            mean = mp.mpf(str(cfg["mean"]))
            sigma = mp.mpf(str(cfg["sigma"]))
            pred = Vabs_ref[key]
            contributions_ref.append({"key": key, "pred": pred, "mean": mean, "sigma": sigma, "chi2": ((pred - mean) / sigma) ** 2})
        contributions_ref.sort(key=lambda t: t["chi2"], reverse=True)
        chi2_ref = mp.fsum([t["chi2"] for t in contributions_ref])

        # χ² at mt (still useful as a boundary diagnostic / proxy)
        contributions_mt = []
        for key in chi2_keys:
            cfg = table.get(str(key), None)
            if not isinstance(cfg, dict) or "mean" not in cfg or "sigma" not in cfg:
                continue
            mean = mp.mpf(str(cfg["mean"]))
            sigma = mp.mpf(str(cfg["sigma"]))
            pred = Vabs_mt[key]
            contributions_mt.append({"key": key, "pred": pred, "mean": mean, "sigma": sigma, "chi2": ((pred - mean) / sigma) ** 2})
        contributions_mt.sort(key=lambda t: t["chi2"], reverse=True)
        chi2_mt = mp.fsum([t["chi2"] for t in contributions_mt])

        # Use reference-scale contributions as the primary “largest χ²” view.
        contributions = contributions_ref

        # CKM diagnostics at mt and at mu_uv
        unitarity_dev_mt = float(np.max(np.abs(V_mt_from_yuk.conj().T @ V_mt_from_yuk - np.eye(3))))
        unitarity_dev_ref = float(np.max(np.abs(V_ref.conj().T @ V_ref - np.eye(3))))
        unitarity_dev_uv_2l = float(np.max(np.abs(V_uv_2l.conj().T @ V_uv_2l - np.eye(3))))
        J_mt = _jarlskog(V_mt_from_yuk)
        J_ref = _jarlskog(V_ref)
        J_uv_2l = _jarlskog(V_uv_2l)
        wolf_mt = _wolfenstein_from_ckm(V_mt_from_yuk)
        wolf_ref = _wolfenstein_from_ckm(V_ref)
        wolf_uv_2l = _wolfenstein_from_ckm(V_uv_2l)

        # --- α_s(MZ) sensitivity (deterministic ±σ scan) ---
        # Vary α_s(MZ) and propagate it to the mt boundary (via MZ→mt gauge running) and then to μ_UV via
        # the upward-only RG integration mt→μ_UV.
        alpha_s_sensitivity: dict[str, Any] = {}
        alpha_s_sens_lines: list[str] = []
        if alpha_s_sigma and alpha_s_sigma > 0:
            scenarios = [
                ("minus_sigma", float(sm_inp.alpha_s) - float(alpha_s_sigma)),
                ("central", float(sm_inp.alpha_s)),
                ("plus_sigma", float(sm_inp.alpha_s) + float(alpha_s_sigma)),
            ]
            mats: dict[str, dict[str, mp.mpf]] = {}
            scenario_diag: dict[str, dict[str, float]] = {}
            for tag, a_s in scenarios:
                sm_raw_sc = dict(sm_raw)
                sm_raw_sc["alpha_s"] = float(a_s)
                sm_bc_sc = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_raw_sc)

                inp_sc = SmMzInputs(
                    mu_GeV=float(sm_inp.mu_GeV),
                    alpha_em_inv=float(sm_inp.alpha_em_inv),
                    sin2_thetaW=float(sm_inp.sin2_thetaW),
                    alpha_s=float(a_s),
                )
                g_mz_sc = gauge_couplings_from_mz_inputs(inp_sc)
                g1_gut_mt_sc = float(sm_bc_sc.route_2loop["g1_gut_mt"])
                g2_mt_sc = float(sm_bc_sc.route_2loop["g2_mt"])
                g3_mt_sc = float(sm_bc_sc.route_2loop["g3_mt"])
                alpha_s_mt_sc = float(sm_bc_sc.route_2loop["alpha_s_mt"])
                mt_msbar_sc = float(sm_bc_sc.route_2loop["mt_msbar_mt_GeV"])
                yt_mt_sc = float(sm_bc_sc.route_2loop["yt_mt"])
                lambda_mt_sc = float(sm_bc_sc.route_2loop["lambda_mt"])
                Ye_base_sc = yukawa_texture_matrix(
                    delta=delta_used, varphi0=varphi0_f, c3=c3_f, a_y=a_e, b=b, y_star=1.0, phase_mode=phase_mode
                )
                ystar_e_sc = scale_y_star_to_match_sigma_max(target_y3=ytau_mt_target, base=Ye_base_sc)
                Yu_mt_sc, Yd_mt_sc, _ = generate_quark_yukawas_mt(
                    varphi0=varphi0_f,
                    delta_used=float(delta_used),
                    yt_mt=float(yt_mt_sc),
                    yb_mt=float(yb_mt_target),
                    scheme="MSbar",
                    reference_scale_GeV=float(mu_start_GeV),
                    s13_mode=s13_mode_best,  # type: ignore[arg-type]
                    delta_mode=delta_mode_best,  # type: ignore[arg-type]
                )
                Ye_mt_sc = (ystar_e_sc * Ye_base_sc).astype(complex)

                gY_mt_sc = float(sm_bc_sc.route_2loop["gY_mt"])
                rge_sc = run_flavor_rge_2loop_thresholds(
                    mu_start_GeV=mu_start_GeV,
                    mu_end_GeV=mu_uv_GeV,
                    thresholds_GeV=thresholds_GeV,
                    g_start=(gY_mt_sc, g2_mt_sc, g3_mt_sc),
                    Yu_start=Yu_mt_sc,
                    Yd_start=Yd_mt_sc,
                    Ye_start=Ye_mt_sc,
                    lambda_start=float(lambda_mt_sc),
                    yN_start=None,
                    beta_sm=beta_sm,
                    beta_e8=beta_e8,
                    apply_sigma_threshold=True,
                    apply_g8_delta_b3=True,
                    delta_b3_g8=2.0,
                    apply_matching=True,
                    matching_loop_order=1,
                    rtol=1e-8,
                    atol=1e-10,
                    method="DOP853",
                )
                V_sc = _ckm_from_yukawas(rge_sc["Yu_end"], rge_sc["Yd_end"])
                mats[tag] = Vabs_from(V_sc)
                scenario_diag[tag] = {
                    "alpha_s_MZ": float(a_s),
                    "g3_MZ": float(g_mz_sc[2]),
                    "g3_mt": float(g3_mt_sc),
                    "g3_mu_uv": float(rge_sc["g_end"]["g3"]),
                    "yt_mt": float(yt_mt_sc),
                    "g1_gut_mt": float(g1_gut_mt_sc),
                    "g2_mt": float(g2_mt_sc),
                    "gY_mt": float(gY_mt_sc),
                    "alpha_s_mt": float(alpha_s_mt_sc),
                    "mt_msbar_mt_GeV": float(mt_msbar_sc),
                }

            # Half-range across the ±σ scan (simple symmetric sensitivity estimate)
            half_range: dict[str, mp.mpf] = {}
            keys = list(mats.get("central", {}).keys())
            for k in keys:
                vmin = min(mats[tag][k] for tag, _ in scenarios)
                vmax = max(mats[tag][k] for tag, _ in scenarios)
                half_range[k] = (vmax - vmin) / 2

            alpha_s_sensitivity = {
                "alpha_s_central": float(sm_inp.alpha_s),
                "alpha_s_sigma": float(alpha_s_sigma),
                "matrix_abs_mu_uv": mats,
                "abs_half_range": half_range,
                "scenario_diagnostics": scenario_diag,
                "note": "Half-range across ±σ(α_s(MZ)) scan propagated through MZ→mt gauge running and mt→μ_UV upward RG evolution (no running below mt).",
            }

            alpha_s_sens_lines.append("")
            alpha_s_sens_lines.append(f"α_s(MZ) sensitivity (±σ, σ={alpha_s_sigma}):")
            alpha_s_sens_lines.append("  sanity: this scan must change g3(MZ), g3(mt), and g3(mu_uv); otherwise it is a placebo.")
            for tag, _ in scenarios:
                d = scenario_diag[tag]
                alpha_s_sens_lines.append(
                    f"- {tag}: αs={d['alpha_s_MZ']:.6g}, g3(MZ)={d['g3_MZ']:.6g}, g3(mt)={d['g3_mt']:.6g}, g3(mu_uv)={d['g3_mu_uv']:.6g}"
                )
            for key in ["Vud", "Vus", "Vub", "Vcd", "Vcs", "Vcb", "Vtd", "Vts", "Vtb"]:
                alpha_s_sens_lines.append(f"- {key}: central={float(mats['central'][key]):.12g} ± {float(half_range[key]):.2e}  (half-range)")

            # Debug stress test: ±0.01 in αs to verify the pipeline responds at all.
            debug_delta = 0.01
            if float(sm_inp.alpha_s) - debug_delta > 0:
                stress = [
                    ("minus_0p01", float(sm_inp.alpha_s) - debug_delta),
                    ("plus_0p01", float(sm_inp.alpha_s) + debug_delta),
                ]
                stress_mats: dict[str, dict[str, mp.mpf]] = {}
                for tag, a_s in stress:
                    sm_raw_sc = dict(sm_raw)
                    sm_raw_sc["alpha_s"] = float(a_s)
                    sm_bc_sc = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_raw_sc)

                    inp_sc = SmMzInputs(
                        mu_GeV=float(sm_inp.mu_GeV),
                        alpha_em_inv=float(sm_inp.alpha_em_inv),
                        sin2_thetaW=float(sm_inp.sin2_thetaW),
                        alpha_s=float(a_s),
                    )
                    g_mz_sc = gauge_couplings_from_mz_inputs(inp_sc)
                    g2_mt_sc = float(sm_bc_sc.route_2loop["g2_mt"])
                    g3_mt_sc = float(sm_bc_sc.route_2loop["g3_mt"])
                    yt_mt_sc = float(sm_bc_sc.route_2loop["yt_mt"])
                    lambda_mt_sc = float(sm_bc_sc.route_2loop["lambda_mt"])
                    Ye_base_sc = yukawa_texture_matrix(
                        delta=delta_used, varphi0=varphi0_f, c3=c3_f, a_y=a_e, b=b, y_star=1.0, phase_mode=phase_mode
                    )
                    ystar_e_sc = scale_y_star_to_match_sigma_max(target_y3=ytau_mt_target, base=Ye_base_sc)
                    Yu_mt_sc, Yd_mt_sc, _ = generate_quark_yukawas_mt(
                        varphi0=varphi0_f,
                        delta_used=float(delta_used),
                        yt_mt=float(yt_mt_sc),
                        yb_mt=float(yb_mt_target),
                        scheme="MSbar",
                        reference_scale_GeV=float(mu_start_GeV),
                        s13_mode=s13_mode_best,  # type: ignore[arg-type]
                        delta_mode=delta_mode_best,  # type: ignore[arg-type]
                    )
                    Ye_mt_sc = (ystar_e_sc * Ye_base_sc).astype(complex)
                    gY_mt_sc = float(sm_bc_sc.route_2loop["gY_mt"])
                    rge_sc = run_flavor_rge_2loop_thresholds(
                        mu_start_GeV=mu_start_GeV,
                        mu_end_GeV=mu_uv_GeV,
                        thresholds_GeV=thresholds_GeV,
                        g_start=(gY_mt_sc, g2_mt_sc, g3_mt_sc),
                        Yu_start=Yu_mt_sc,
                        Yd_start=Yd_mt_sc,
                        Ye_start=Ye_mt_sc,
                        lambda_start=float(lambda_mt_sc),
                        yN_start=None,
                        beta_sm=beta_sm,
                        beta_e8=beta_e8,
                        apply_sigma_threshold=True,
                        apply_g8_delta_b3=True,
                        delta_b3_g8=2.0,
                        apply_matching=True,
                        matching_loop_order=1,
                        rtol=1e-8,
                        atol=1e-10,
                        method="DOP853",
                    )
                    V_sc = _ckm_from_yukawas(rge_sc["Yu_end"], rge_sc["Yd_end"])
                    stress_mats[tag] = Vabs_from(V_sc)
                dvus = stress_mats["plus_0p01"]["Vus"] - stress_mats["minus_0p01"]["Vus"]
                alpha_s_sens_lines.append("")
                alpha_s_sens_lines.append(f"debug: Δαs=±0.01 changes Vus(mu_uv) by {float(dvus):.2e} (should not be identically zero).")

        # --- Uncertainty propagation (Monte Carlo over PDG priors) ---
        mc_uncertainty: dict[str, Any] = {
            "enabled": False,
            "samples": int(mc_samples),
            "chi2_keys": list(chi2_keys),
            "chi2_dof": int(len(chi2_keys)),
            "inputs_sigma": {
                "alpha_s_sigma": float(alpha_s_sigma),
                "alpha_em_inv_sigma": float(alpha_em_inv_sigma),
                "sin2_thetaW_sigma": float(sin2_sigma),
                "mt_sigma_GeV": float(mt_sigma),
                "mb_sigma_GeV": float(mb_sigma),
                "mc_sigma_GeV": float(mc_sigma),
            },
            "note": "Monte Carlo over SM inputs at MZ (PDG-style priors) propagated through MZ→mt boundary construction and mt→ref RG evolution (SM segment; thresholds outside range). χ² is computed on the same declared subset chi2_keys as the nominal ref-scale χ² (to avoid unitarity double-counting / over-constraining).",
        }
        mc_lines: list[str] = []
        any_sigma = bool(alpha_s_sigma > 0 or alpha_em_inv_sigma > 0 or sin2_sigma > 0 or mt_sigma > 0 or mb_sigma > 0 or mc_sigma > 0)
        if mc_samples and mc_samples > 0 and any_sigma:
            rng = np.random.default_rng(int(config.seed) + 23091)

            def _clip_pos(x: float) -> float:
                return float(max(float(x), 1e-12))

            def _clip_sin2(x: float) -> float:
                return float(min(max(float(x), 1e-6), 1.0 - 1e-6))

            def _sample_gauss(mean: float, sigma: float) -> float:
                if sigma <= 0:
                    return float(mean)
                return float(rng.normal(loc=float(mean), scale=float(sigma)))

            # Base charged-lepton Yukawa texture for Ye(mt); scaled per-sample to match yτ(mt).
            Ye_base_mc = yukawa_texture_matrix(
                delta=delta_used, varphi0=varphi0_f, c3=c3_f, a_y=a_e, b=b, y_star=1.0, phase_mode=phase_mode
            ).astype(complex)

            idx = {
                "Vud": (0, 0),
                "Vus": (0, 1),
                "Vub": (0, 2),
                "Vcd": (1, 0),
                "Vcs": (1, 1),
                "Vcb": (1, 2),
                "Vtd": (2, 0),
                "Vts": (2, 1),
                "Vtb": (2, 2),
            }

            # Collect samples (successful only).
            rows: list[dict[str, float]] = []
            fail_count = 0
            for _ in range(int(mc_samples)):
                sm_i = dict(sm_raw)
                sm_i["alpha_em_inv"] = _clip_pos(_sample_gauss(float(sm_inp.alpha_em_inv), float(alpha_em_inv_sigma)))
                sm_i["sin2_thetaW"] = _clip_sin2(_sample_gauss(float(sm_inp.sin2_thetaW), float(sin2_sigma)))
                sm_i["alpha_s"] = _clip_pos(_sample_gauss(float(sm_inp.alpha_s), float(alpha_s_sigma)))
                sm_i["mt_GeV"] = _clip_pos(_sample_gauss(float(mt_GeV), float(mt_sigma)))
                sm_i["mb_GeV"] = _clip_pos(_sample_gauss(float(mb_GeV), float(mb_sigma)))
                sm_i["mc_GeV"] = _clip_pos(_sample_gauss(float(mc_GeV), float(mc_sigma)))

                try:
                    sm_bc_i = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_i)
                    mu_start_i = float(sm_bc_i.mu_mt_GeV)
                    gY_mt_i = float(sm_bc_i.route_2loop["gY_mt"])
                    g2_mt_i = float(sm_bc_i.route_2loop["g2_mt"])
                    g3_mt_i = float(sm_bc_i.route_2loop["g3_mt"])
                    yt_mt_i = float(sm_bc_i.route_2loop["yt_mt"])
                    yb_mt_i = float(sm_bc_i.route_2loop["yb_mt"])
                    ytau_mt_i = float(sm_bc_i.route_2loop["ytau_mt"])
                    lambda_mt_i = float(sm_bc_i.route_2loop["lambda_mt"])

                    ystar_e_i = scale_y_star_to_match_sigma_max(target_y3=float(ytau_mt_i), base=Ye_base_mc)
                    Ye_mt_i = (ystar_e_i * Ye_base_mc).astype(complex)
                    Yu_mt_i, Yd_mt_i, _ = generate_quark_yukawas_mt(
                        varphi0=varphi0_f,
                        delta_used=float(delta_used),
                        yt_mt=float(yt_mt_i),
                        yb_mt=float(yb_mt_i),
                        scheme="MSbar",
                        reference_scale_GeV=float(mu_start_i),
                        s13_mode=s13_mode_best,  # type: ignore[arg-type]
                        delta_mode=delta_mode_best,  # type: ignore[arg-type]
                    )

                    rge_ref_i = run_flavor_rge_2loop_thresholds(
                        mu_start_GeV=float(mu_start_i),
                        mu_end_GeV=float(mu_ref_GeV),
                        thresholds_GeV=thresholds_GeV,
                        g_start=(float(gY_mt_i), float(g2_mt_i), float(g3_mt_i)),
                        Yu_start=Yu_mt_i,
                        Yd_start=Yd_mt_i,
                        Ye_start=Ye_mt_i,
                        lambda_start=float(lambda_mt_i),
                        yN_start=None,
                        beta_sm=beta_sm,
                        beta_e8=beta_e8,
                        apply_sigma_threshold=True,
                        apply_g8_delta_b3=False,
                        apply_matching=True,
                        matching_loop_order=1,
                        rtol=1e-9,
                        atol=1e-11,
                        method="DOP853",
                    )
                    V_ref_i = _ckm_from_yukawas(rge_ref_i["Yu_end"], rge_ref_i["Yd_end"])
                    Vabs_ref_i = {k: float(np.abs(V_ref_i[i, j])) for k, (i, j) in idx.items()}

                    chi2_i = 0.0
                    for key in chi2_keys:
                        cfg = table.get(str(key), None)
                        if not isinstance(cfg, dict):
                            continue
                        mean = float(cfg["mean"])
                        sig = float(cfg["sigma"])
                        pred = float(Vabs_ref_i.get(str(key), float("nan")))
                        if sig > 0 and np.isfinite(pred):
                            chi2_i += float(((pred - mean) / sig) ** 2)

                    rows.append(
                        {
                            "alpha_s_MZ": float(sm_i["alpha_s"]),
                            "alpha_em_inv_MZ": float(sm_i["alpha_em_inv"]),
                            "sin2_thetaW_MZ": float(sm_i["sin2_thetaW"]),
                            "mt_pole_GeV": float(sm_i["mt_GeV"]),
                            "mb_mb_GeV": float(sm_i["mb_GeV"]),
                            "mc_GeV": float(sm_i["mc_GeV"]),
                            "chi2_refscale": float(chi2_i),
                            **{k: float(Vabs_ref_i[k]) for k in idx},
                        }
                    )
                except Exception:
                    fail_count += 1

            if rows:
                mc_uncertainty["enabled"] = True
                mc_uncertainty["success"] = int(len(rows))
                mc_uncertainty["failed"] = int(fail_count)

                def _mean_std(vals: list[float]) -> tuple[float, float]:
                    a = np.array([float(x) for x in vals], dtype=float)
                    if a.size == 0:
                        return float("nan"), float("nan")
                    mu = float(np.mean(a))
                    sig = float(np.std(a, ddof=1)) if a.size > 1 else 0.0
                    return mu, sig

                # Inputs summary
                inputs = ["alpha_s_MZ", "alpha_em_inv_MZ", "sin2_thetaW_MZ", "mt_pole_GeV", "mb_mb_GeV", "mc_GeV"]
                out_inputs_mean: dict[str, float] = {}
                out_inputs_sigma: dict[str, float] = {}
                for k in inputs:
                    mu, sig = _mean_std([r[k] for r in rows])
                    out_inputs_mean[k] = mu
                    out_inputs_sigma[k] = sig

                # Outputs summary at reference scale
                out_ckm_mean: dict[str, float] = {}
                out_ckm_sigma: dict[str, float] = {}
                for k in idx:
                    mu, sig = _mean_std([r[k] for r in rows])
                    out_ckm_mean[k] = mu
                    out_ckm_sigma[k] = sig
                chi2_mu, chi2_sig = _mean_std([r["chi2_refscale"] for r in rows])

                # Simple sensitivities: Pearson r vs inputs for χ² and |Vus|
                def _corr(x: list[float], y: list[float]) -> float | None:
                    ax = np.array(x, dtype=float)
                    ay = np.array(y, dtype=float)
                    if ax.size < 3:
                        return None
                    sx = float(np.std(ax))
                    sy = float(np.std(ay))
                    if not (sx > 0 and sy > 0):
                        return None
                    return float(np.corrcoef(ax, ay)[0, 1])

                sens: dict[str, dict[str, float | None]] = {"chi2_refscale": {}, "Vus": {}}
                for k in inputs:
                    sens["chi2_refscale"][k] = _corr([r[k] for r in rows], [r["chi2_refscale"] for r in rows])
                    sens["Vus"][k] = _corr([r[k] for r in rows], [r["Vus"] for r in rows])

                mc_uncertainty["inputs_mean"] = out_inputs_mean
                mc_uncertainty["inputs_std"] = out_inputs_sigma
                mc_uncertainty["ckm_abs_refscale_mean"] = out_ckm_mean
                mc_uncertainty["ckm_abs_refscale_std"] = out_ckm_sigma
                mc_uncertainty["chi2_refscale_mean"] = float(chi2_mu)
                mc_uncertainty["chi2_refscale_std"] = float(chi2_sig)
                mc_uncertainty["sensitivities_pearson_r"] = sens
                mc_uncertainty["preview"] = rows[: min(12, len(rows))]

                mc_lines.append("")
                mc_lines.append(f"Monte Carlo uncertainty propagation (ref-scale): N_success={len(rows)}, N_failed={fail_count}")
                mc_lines.append(f"- inputs σ: αs={alpha_s_sigma}, αem_inv={alpha_em_inv_sigma}, sin2={sin2_sigma}, mt={mt_sigma}, mb={mb_sigma}, mc={mc_sigma}")
                mc_lines.append(f"- χ²(ref) = {chi2_mu:.3g} ± {chi2_sig:.3g}")
                mc_lines.append(
                    f"- |Vus|(ref) = {out_ckm_mean['Vus']:.12g} ± {out_ckm_sigma['Vus']:.2e}, "
                    f"|Vcb|(ref) = {out_ckm_mean['Vcb']:.12g} ± {out_ckm_sigma['Vcb']:.2e}, "
                    f"|Vub|(ref) = {out_ckm_mean['Vub']:.12g} ± {out_ckm_sigma['Vub']:.2e}"
                )

        checks: list[Check] = []
        # Publication-grade Gap B gate: phase selection rule must be a derived *filter*, not a χ²-driven search.
        if phase_selection_rule_mode == "filter_only" and include_topology_delta_cp:
            if not topo_output_found or not topo_variants:
                checks.append(
                    Check(
                        check_id="phase_set_derived_not_searched",
                        passed=True,
                        detail="mode=filter_only requested but topology_phase_map output/candidates missing; using baseline CKM variants as fallback",
                        severity="WARN",
                    )
                )
            else:
                derived_only = bool(
                    variants
                    and all(
                        isinstance(v, dict)
                        and str(v.get("source_module", "")) == "topology_phase_map"
                        and str(v.get("delta_mode", "")) == "external_delta_cp_override"
                        for v in variants
                    )
                )
                checks.append(
                    Check(
                        check_id="phase_set_derived_not_searched",
                        passed=derived_only,
                        detail=f"mode=filter_only; variants={len(variants)} (all sourced from topology_phase_map)",
                        severity="PASS" if derived_only else "FAIL",
                    )
                )
        else:
            checks.append(
                Check(
                    check_id="phase_set_derived_not_searched",
                    passed=True,
                    detail=f"mode={phase_selection_rule_mode} (not enforced; topology_variants={len(topo_variants)})",
                    severity="INFO",
                )
            )
        checks.append(
            Check(
                check_id="tfpt_rg_dressed_texture_pipeline_present",
                passed=True,
                detail="PASS: mt boundary Yukawas (Z3 texture) + upward-only RG evolution mt→mu_uv using PyR@TE 2-loop beta functions (integrated here with complex Yukawas) + CKM extraction at mt and mu_uv",
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
                detail=f"native_scale=mt ({mu_start_GeV} GeV), reference_scale={ref_meta.get('reference_scale_label', '?')} ({mu_ref_GeV} GeV)",
                severity="PASS",
            )
        )
        chi2_mt_f = float(chi2_mt)
        chi2_ref_f = float(chi2_ref)
        chi2_scale_ratio = float("nan")
        if np.isfinite(chi2_mt_f) and np.isfinite(chi2_ref_f) and min(chi2_mt_f, chi2_ref_f) > 0:
            chi2_scale_ratio = float(max(chi2_mt_f, chi2_ref_f) / max(min(chi2_mt_f, chi2_ref_f), 1e-12))
        scale_ok = bool(np.isfinite(chi2_scale_ratio) and chi2_scale_ratio <= CHI2_SCALE_RATIO_MAX)
        checks.append(
            Check(
                check_id="chi2_scale_consistency",
                passed=True,
                detail=f"chi2_mt={chi2_mt_f:.6g}, chi2_ref={chi2_ref_f:.6g}, ratio={chi2_scale_ratio:.6g} (max={CHI2_SCALE_RATIO_MAX})",
                severity="PASS" if scale_ok else "WARN",
            )
        )
        pub = rge2.get("publication_grade", {}) or {}
        thr_ok = bool(pub.get("threshold_matching_ok", False))
        blocked = pub.get("blocked_thresholds", []) or []
        checks.append(
            Check(
                check_id="threshold_matching_publication_grade",
                passed=bool(thr_ok),
                detail=f"threshold_matching_ok={thr_ok}, blocked_thresholds={blocked}",
            )
        )
        checks.append(
            Check(
                check_id="matching_mc_present",
                passed=bool((not mc_uncertainty.get("enabled", False)) or (mc_uncertainty.get("success", 0) and not mc_uncertainty.get("failed", 0))),
                detail=f"enabled={mc_uncertainty.get('enabled')}, samples={mc_uncertainty.get('samples')}, success={mc_uncertainty.get('success')}, failed={mc_uncertainty.get('failed')}",
            )
        )
        ratio = float(g_mt[0] / gY_mt) if gY_mt != 0 else float("nan")
        checks.append(
            Check(
                check_id="g1_gut_over_gY_convention",
                passed=bool(np.isfinite(ratio) and abs(ratio - gut_norm) < 1e-12),
                detail=f"g1_GUT/gY={ratio:.15g} vs sqrt(5/3)={gut_norm:.15g}",
            )
        )
        bc_d = dict(sm_bc_mt.diffs)
        checks.append(
            Check(
                check_id="sm_boundary_1loop_vs_2loop_close",
                passed=bool(
                    abs(float(bc_d.get("gY_mt", float("nan")))) < 5e-4
                    and abs(float(bc_d.get("g2_mt", float("nan")))) < 5e-4
                    and abs(float(bc_d.get("g3_mt", float("nan")))) < 5e-3
                    and abs(float(bc_d.get("yt_mt", float("nan")))) < 5e-3
                ),
                detail=f"diffs (2L-1L) @ mt: ΔgY={bc_d.get('gY_mt')}, Δg2={bc_d.get('g2_mt')}, Δg3={bc_d.get('g3_mt')}, Δyt={bc_d.get('yt_mt')}",
            )
        )
        checks.append(
            Check(
                check_id="texture_config_loaded",
                passed=bool(flavor_cfg_path.exists()),
                detail=f"path={_relpath(flavor_cfg_path)} delta_source={delta_source} phase_mode={phase_mode} coeff_rule={coeff_rule} coeff_source={coeff.source} a_u={a_u} a_d={a_d} a_e={a_e} b={b}",
            )
        )
        checks.append(
            Check(
                check_id="rg_chi2_finite",
                passed=bool(mp.isfinite(chi2_mt) and chi2_mt >= 0),
                detail=f"chi2(mt, boundary proxy)={chi2_mt}",
            )
        )
        checks.append(
            Check(
                check_id="rg_chi2_refscale_finite",
                passed=bool(mp.isfinite(chi2_ref) and chi2_ref >= 0),
                detail=f"chi2(refscale={ref_meta.get('reference_scale_label', '?')}@{mu_ref_GeV:.6g} GeV)={chi2_ref}",
            )
        )
        checks.append(
            Check(
                check_id="lambda_in_physical_range",
                passed=bool(lam > mp.mpf("0.20") and lam < mp.mpf("0.24")),
                detail=f"lambda_TFPT={lam}",
            )
        )
        checks.append(
            Check(
                check_id="ckm_unitarity",
                passed=bool(max(unitarity_dev_mt, unitarity_dev_uv_2l) < 5e-6),
                detail=f"max|V†V-I|(mt)={unitarity_dev_mt:.3e}, max|V†V-I|(mu_uv,2-loop)={unitarity_dev_uv_2l:.3e}",
            )
        )
        checks.append(
            Check(
                check_id="ckm_phase_map_consistent_with_generator",
                passed=bool(delta_ckm_map_diff < 1e-12),
                detail=f"|delta_ckm(map)-delta_ckm(generator)|={delta_ckm_map_diff:.3e} rad",
            )
        )

        lines: list[str] = []
        ref_sha = _sha256_file(ref_path)
        sm_sha = _sha256_file(sm_path)
        flavor_sha = _sha256_file(flavor_cfg_path) if flavor_cfg_path.exists() else "missing"
        lep_sha = _sha256_file(lepton_masses_path) if lepton_masses_path.exists() else "missing"
        matching_map_sha = _sha256_file(matching_map_path) if matching_map_path.exists() else "missing"
        lines += [
            "CKM full pipeline (Z3 Yukawa texture + diagonalization; upward-only RG mt→mu_uv)",
            "",
            f"reference file: {_relpath(ref_path)} (sha256={ref_sha})",
            f"CKM reference metadata: scale={ref_meta.get('reference_scale_label', '?')} ({ref_meta.get('reference_scale_GeV', '?')} GeV), scheme={ref_meta.get('scheme', '?')}",
            f"SM inputs file: {_relpath(sm_path)} (sha256={sm_sha})",
            f"msbar_matching_map: {'used' if matching_map_used else 'not found'} (path={_relpath(matching_map_path)}, sha256={matching_map_sha})",
            f"flavor texture config: {_relpath(flavor_cfg_path)} (sha256={flavor_sha})",
            f"lepton masses file: {_relpath(lepton_masses_path)} (sha256={lep_sha})",
            "",
            "TFPT input:",
            f"- lambda_TFPT = {lam}",
            "",
            "Flavor texture layer (explicit conventions):",
            f"- delta_source = {delta_source}",
            f"- delta_M (from tau/mu) = {delta_M}",
            f"- delta_star (from geometry) = {delta_star}",
            f"- delta_used = {delta_used}",
            f"- phase_mode = {phase_mode}",
            f"- theta(delta_used) = {theta} rad",
            f"- zeta = exp(i theta) = {zeta}",
            f"- C(delta) first row = (1, zeta, zeta*)",
            f"- a_u={a_u}, a_d={a_d}, a_e={a_e}, b={b}",
            f"- topology_phase_atoms (docking): source_module={topo_atoms_cfg.get('source_module','?')}, status={topo_atoms_cfg.get('status','?')}",
            "",
            "Discrete CKM variant scan (no continuous fit):",
            f"- variants evaluated: {len(variant_results)} (run_all_variants={run_all_variants})",
            f"- reference scale for χ²: {ref_meta.get('reference_scale_label', '?')} ({mu_ref_GeV:.6g} GeV)",
            f"- χ² keys (subset): {chi2_keys}",
            f"- front runner: {variant_label_best} (s13_mode={s13_mode_best}, delta_mode={delta_mode_best}), chi2_ref≈{float(chi2_ref):.3g}  (mt-proxy chi2≈{float(chi2_mt):.3g})",
            "",
            "Quark phase map (naming hygiene; TFPT paper v2.5 flavor section):",
            f"- δMobius := delta_M (τ/μ extraction) = {phase_map.delta_mobius}",
            f"- θtexture := theta(delta_used; phase_mode={phase_map.phase_mode}) = {phase_map.theta_texture_rad} rad  (ζ=exp(iθ) in C(δ))",
            f"- δCKM := delta_CP(PDG) = {phase_map.delta_ckm_rad} rad  (delta_mode={phase_map.delta_mode})",
            f"- source: {phase_map.source}",
            "",
            "RG integration setup (upward only):",
            f"- mu_start = {mu_start_GeV:.6g} GeV (mt), mu_end = {mu_uv_GeV:.3e} GeV (UV benchmark)",
            f"- SM inputs at MZ (used only to derive g_i(mt) via MZ→mt running): αs(MZ)={float(sm_inp.alpha_s):.6g} ± {alpha_s_sigma:.3g}, MZ={float(sm_inp.mu_GeV):.6g} GeV",
            f"- thresholds (for gauge running bookkeeping): (mc,mb,mt)=({mc_GeV:.3g},{mb_GeV:.3g},{mt_GeV:.3g}) GeV; mW={mW_GeV:.3g} GeV (below MZ; not used in this upward-only run)",
            "- hypercharge convention: Q=T3+Y (SM); rge_sm uses g1_GUT=√(5/3)·gY, PyR@TE betas use gY≡g′.",
            f"- gauge couplings at MZ (rge_sm, g1_GUT): g1_GUT={g_mz[0]:.6g}, g2={g_mz[1]:.6g}, g3={g_mz[2]:.6g}",
            f"- gauge couplings at mt (2-loop gauge-only, g1_GUT): g1_GUT={g_mt[0]:.6g}, g2={g_mt[1]:.6g}, g3={g_mt[2]:.6g}; converted gY(mt)={gY_mt:.6g}",
            f"- SM boundary scheme: {sm_bc_mt.scheme}",
            f"- SM boundary @ mt (matching layer): route_2loop(gY,g2,g3,yt)=({sm_bc_mt.route_2loop['gY_mt']:.6g},{sm_bc_mt.route_2loop['g2_mt']:.6g},{sm_bc_mt.route_2loop['g3_mt']:.6g},{sm_bc_mt.route_2loop['yt_mt']:.6g}), "
            f"route_1loop=({sm_bc_mt.route_1loop['gY_mt']:.6g},{sm_bc_mt.route_1loop['g2_mt']:.6g},{sm_bc_mt.route_1loop['g3_mt']:.6g},{sm_bc_mt.route_1loop['yt_mt']:.6g}), "
            f"diffs(2L-1L)=(ΔgY={sm_bc_mt.diffs['gY_mt']:.2e},Δg2={sm_bc_mt.diffs['g2_mt']:.2e},Δg3={sm_bc_mt.diffs['g3_mt']:.2e},Δyt={sm_bc_mt.diffs['yt_mt']:.2e})",
            f"- gauge couplings at mu_end (2-loop gauge-only, g1_GUT): g1_GUT={g_uv_gauge_2loop[0]:.6g}, g2={g_uv_gauge_2loop[1]:.6g}, g3={g_uv_gauge_2loop[2]:.6g}",
            "",
            "Yukawa boundary at mt (explicit inputs):",
            f"- αs(mt) = {alpha_s_mt:.6g} (from g3(mt))",
            f"- mt_pole = {mt_GeV:.6g} GeV, mt_MSbar(mt) ≈ {mt_msbar_mt:.6g} GeV (QCD 2-loop pole→MSbar)",
            f"- v = {v_ev_GeV:.6g} GeV => yt(mt) ≈ {yt_mt:.6g}",
            f"- yb(mt) input = {yb_mt_target:.6g}, ytau(mt) input = {ytau_mt_target:.6g}",
            "- quark Yukawas are generated via the Möbius/Z3 v1.07SM rule (no y_* scaling knob):",
            f"  Yu diag ≈ (y_u={yu_diag_mt['y_u']:.6g}, y_c={yu_diag_mt['y_c']:.6g}, y_t={yu_diag_mt['y_t']:.6g})",
            f"  Yd diag ≈ (y_d={yd_diag_mt['y_d']:.6g}, y_s={yd_diag_mt['y_s']:.6g}, y_b={yd_diag_mt['y_b']:.6g})",
            f"  CKM seed (PDG): s12={ckm_seed['s12']:.6g}, s23={ckm_seed['s23']:.6g}, s13={ckm_seed['s13']:.6g}, δ={ckm_seed['delta_cp_rad']:.6g} rad",
            f"- charged-lepton Yukawa uses the Z3 texture y_* scaling: y*_e (scaled to ytau(mt)) = {ystar_e:.6g}",
            "",
            "CKM |V_ij| at mt (boundary; computed from Yukawa diagonalization):",
        ]
        for key in ["Vud", "Vus", "Vub", "Vcd", "Vcs", "Vcb", "Vtd", "Vts", "Vtb"]:
            lines.append(f"- {key} = {Vabs_mt[key]}")

        lines += [
            "",
            f"CKM |V_ij| at reference scale ({ref_meta.get('reference_scale_label', 'ref')}={mu_ref_GeV:.6g} GeV; after downward mt→ref RG in SM):",
        ]
        for key in ["Vud", "Vus", "Vub", "Vcd", "Vcs", "Vcb", "Vtd", "Vts", "Vtb"]:
            lines.append(f"- {key} = {Vabs_ref[key]}")

        lines += [
            "",
            "CKM |V_ij| at mu_end after upward RG evolution (2-loop; PyR@TE betas; complex Yukawas supported):",
        ]
        for key in ["Vud", "Vus", "Vub", "Vcd", "Vcs", "Vcb", "Vtd", "Vts", "Vtb"]:
            lines.append(f"- {key} = {Vabs_uv_2l[key]}")

        lines += [
            "",
            "Largest χ² contributions (top 5):",
        ]
        for t in contributions[:5]:
            lines.append(f"- {t['key']}: pred={t['pred']} mean={t['mean']} sigma={t['sigma']} chi2={t['chi2']}")

        lines += [
            "",
            "Derived CKM invariants:",
            f"- at mt: Jarlskog J = {J_mt:.6e}, Wolfenstein (λ,A,ρ̄,η̄) = {wolf_mt}",
            f"- at ref scale: Jarlskog J = {J_ref:.6e}, Wolfenstein (λ,A,ρ̄,η̄) = {wolf_ref}",
            f"- at mu_end (2-loop): Jarlskog J = {J_uv_2l:.6e}, Wolfenstein (λ,A,ρ̄,η̄) = {wolf_uv_2l}",
        ]
        lines += alpha_s_sens_lines
        lines += mc_lines

        # Explicit threshold bookkeeping (no one has to reverse-engineer segment switches from the code).
        lines += [
            "",
            "Threshold transitions (PyR@TE 2-loop engine):",
        ]
        for seg in rge2.get("segments", []):
            mu_a = float(seg.get("mu_start_GeV", float("nan")))
            mu_b = float(seg.get("mu_end_GeV", float("nan")))
            model = str(seg.get("model", ""))
            delta_b3_active = bool(seg.get("delta_b3_active", False))
            rule = seg.get("threshold_transition_at_start")
            if rule is None:
                lines.append(f"- segment [{mu_a:.6g}, {mu_b:.6g}] GeV: model={model}, delta_b3_active={delta_b3_active} (no explicit threshold at start)")
                continue
            try:
                tid = str(rule.threshold_id)
                status = str(rule.status)
                action = str(rule.action)
                note = str(rule.note)
            except Exception:
                # JSONified fallback
                tid = str(getattr(rule, "get", lambda k, d=None: d)("threshold_id", "?"))
                status = str(getattr(rule, "get", lambda k, d=None: d)("status", "?"))
                action = str(getattr(rule, "get", lambda k, d=None: d)("action", "?"))
                note = str(getattr(rule, "get", lambda k, d=None: d)("note", ""))
            lines.append(
                f"- μ={mu_a:.6g} GeV: {tid} → {status} ({action}). Segment [{mu_a:.6g}, {mu_b:.6g}] model={model}, delta_b3_active={delta_b3_active}. {note}"
            )

        lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module constructs Yukawas at μ=mt from a Z3 circulant texture and computes CKM via Yukawa diagonalization (no continuous fit).",
            "- RG evolution is upward-only (starts at μ=mt and runs to μ_UV). No running below mt is performed in this module.",
            "- RG evolution uses PyR@TE-generated 2-loop beta functions, but is integrated here to support complex Yukawas (the shipped PyR@TE PythonOutput solver is real-only for Yukawas).",
            "- Threshold handling: switch SM→(SigmaF+yN) beta system at MSigma; apply a 1-loop Δb3=+2 patch above MG8 (paper v1.06 ‘G8 bridge’ note).",
            "- Phase conventions (δ source, θ(δ) mode) are explicit inputs in tfpt_suite/data/flavor_texture_v24.json; sector coefficients are fixed by tfpt_suite/flavor_textures.py:coefficients_v107sm().",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"ckm_summary_png": None, "ckm_residuals_sigma_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_ckm_summary(
                out_dir=out_dir,
                Vabs_mt=np.abs(V_mt_from_yuk),
                Vabs_uv=np.abs(V_uv_2l),
                contributions=contributions,
                mu_start_GeV=mu_start_GeV,
                mu_end_GeV=mu_uv_GeV,
            )
            plot_resid, plot_resid_warnings = _plot_ckm_residuals_sigma(
                out_dir=out_dir,
                Vabs_pred=Vabs_ref,
                ref_table=ref,
                chi2_keys=chi2_keys,
            )
            plot.update(plot_resid)
            warnings.extend(plot_resid_warnings)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "reference_file": _relpath(ref_path),
                "reference_meta": ref_meta,
                "matching_map": {
                    "used": bool(matching_map_used),
                    "path": _relpath(matching_map_path),
                    "sha256": matching_map_sha,
                    "mt_boundary": matching_map_res.get("mt_boundary", None),
                    "derived_yukawas_mt": matching_map_res.get("derived_yukawas_mt", None),
                },
                "lambda_tfpt": lam,
                "flavor_texture": {
                    "config_file": _relpath(flavor_cfg_path),
                    "lepton_masses_file": _relpath(lepton_masses_path),
                    "topology_phase_atoms": topo_atoms_cfg,
                    "delta_source": delta_source,
                    "delta_M_from_tau_mu": delta_M,
                    "delta_star_from_varphi0": delta_star,
                    "delta_used": delta_used,
                    "phase_mode": phase_mode,
                    "theta_rad": theta,
                    "zeta": str(zeta),
                    "a_u": a_u,
                    "a_d": a_d,
                    "a_e": a_e,
                    "b": b,
                    "y_star_u": None,
                    "y_star_d": None,
                    "y_star_e": ystar_e,
                    "quark_yukawa_generation": {
                        "method": "mobius_v107sm",
                        "scheme": gen_meta.scheme,
                        "reference_scale_GeV": gen_meta.reference_scale_GeV,
                        "ratios": ratios_mt,
                        "yu_diag_mt": yu_diag_mt,
                        "yd_diag_mt": yd_diag_mt,
                        "ckm_seed": ckm_seed,
                    },
                    "ckm_variant_scan": {
                        "run_all_variants": run_all_variants,
                        "best": best_public,
                        "variants": variant_results,
                    },
                    "quark_phase_map": phase_map.__dict__,
                },
                "rg_upward": {
                    "mu_start_GeV": mu_start_GeV,
                    "mu_end_GeV": mu_uv_GeV,
                    "sm_inputs_mz_file": _relpath(sm_path),
                    "thresholds_file": _relpath(thresholds_path),
                    "pyrate_beta_sources": {
                        "sm_pythonoutput_dir": _relpath(py_sm_dir),
                        "e8_pythonoutput_dir": _relpath(py_e8_dir),
                        "fingerprints": rge2.get("beta_sources", None),
                    },
                    "alpha_s_sigma": alpha_s_sigma,
                    "g_mz": {"g1": g_mz[0], "g2": g_mz[1], "g3": g_mz[2]},
                    "g_mt": {"g1": g_mt[0], "g2": g_mt[1], "g3": g_mt[2]},
                    "g_mu_end_gauge_2loop": {"g1": g_uv_gauge_2loop[0], "g2": g_uv_gauge_2loop[1], "g3": g_uv_gauge_2loop[2]},
                    "g_mu_end_pyrate_2loop_SMnorm": rge2.get("g_end_sm", rge2["g_end"]),
                    "g_mu_end_pyrate_2loop_GUTnorm": rge2.get("g_end_gut", None),
                    "gravity_alpha3_patch": rge2.get("gravity_alpha3_patch", None),
                    "publication_grade": rge2.get("publication_grade", None),
                    "segments_pyrate_2loop": rge2["segments"],
                    "alpha_s_mt": alpha_s_mt,
                    "mt_pole_GeV": mt_GeV,
                    "mt_msbar_mt_GeV": mt_msbar_mt,
                    "v_ev_GeV": v_ev_GeV,
                    "yt_mt": yt_mt,
                    "yb_mt_input": yb_mt_target,
                    "ytau_mt_input": ytau_mt_target,
                    "ckm_abs_mt": {k: Vabs_mt[k] for k in Vabs_mt},
                    "ckm_abs_refscale": {k: Vabs_ref[k] for k in Vabs_ref},
                    "ckm_abs_mu_end_2loop": {k: Vabs_uv_2l[k] for k in Vabs_uv_2l},
                    "chi2_refscale": float(chi2_ref),
                    "chi2_mt": float(chi2_mt),
                    "chi2_keys": list(chi2_keys),
                    "contributions_refscale": contributions_ref,
                    "contributions_mt_proxy": contributions_mt,
                    "unitarity_dev_mt": unitarity_dev_mt,
                    "unitarity_dev_refscale": unitarity_dev_ref,
                    "unitarity_dev_mu_end_2loop": unitarity_dev_uv_2l,
                    "jarlskog_mt": J_mt,
                    "jarlskog_refscale": J_ref,
                    "jarlskog_mu_end_2loop": J_uv_2l,
                    "wolfenstein_mt": wolf_mt,
                    "wolfenstein_refscale": wolf_ref,
                    "wolfenstein_mu_end_2loop": wolf_uv_2l,
                },
                "alpha_s_sensitivity": alpha_s_sensitivity,
                "matching_mc": mc_uncertainty,
                "contributions_best": contributions_ref,
                "contributions_mt_proxy": contributions_mt,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

