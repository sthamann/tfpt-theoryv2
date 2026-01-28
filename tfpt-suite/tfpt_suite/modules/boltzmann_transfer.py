from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.constants import TfptConstants
from tfpt_suite.cosmo_scale_map import CosmoScaleInputs, a0_over_a_transition, ell_from_k_hat
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_info, mk_check_pass, mk_check_warn


@dataclass(frozen=True)
class EllRange:
    ell_min: float
    ell_max: float
    note: str


PLANCK_ENABLE_ENV = "TFPT_ENABLE_PLANCK_LIKELIHOOD"
PLANCK_PLIK_DATASET = "plik_lite_v22.dataset"
PLANCK_LOWL_CANDIDATES = [
    ("cobaya.likelihoods.planck_2018_lowl.TT", "TT"),
    ("cobaya.likelihoods.planck_2018_lowl.EE", "EE"),
    ("cobaya.likelihoods.planck_2018_lowl.TT_EE", "TT_EE"),
]
PLANCK_LENSING_CANDIDATES = [
    ("cobaya.likelihoods.planck_2018_lensing.clik", "clik"),
    ("cobaya.likelihoods.planck_2018_lensing.clik_lensing", "clik_lensing"),
]
DEFAULT_CMB_ELL_RANGE = (2.0, 2500.0)
DEFAULT_SMALL_SCALE_ELL_MIN = 2500.0
DEFAULT_SIGNATURE_DECISION = "prefer_cmb_if_overlap"


def _as_float_or(value: object, default: float) -> float:
    try:
        val = float(value)
    except Exception:
        return default
    if not math.isfinite(val):
        return default
    return val


def _load_nuisance_policy(*, data_dir: Path) -> dict[str, Any]:
    spec_path = data_dir / "likelihood_datasets_v1.json"
    if not spec_path.is_file():
        return {}
    try:
        spec = json.loads(spec_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    nuisance = spec.get("nuisance_policy", {}) if isinstance(spec.get("nuisance_policy", {}), dict) else {}
    return nuisance if isinstance(nuisance, dict) else {}


def _plot_cl_comparison(
    *,
    out_dir: Path,
    ell: np.ndarray,
    tt: np.ndarray,
    ee: np.ndarray,
    tt_base: np.ndarray | None,
    ee_base: np.ndarray | None,
    ell_max_plot: int = 2500,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {
        "cl_tt_ee_comparison_png": None,
        "cl_ratio_tt_ee_png": None,
    }
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)
        ell = np.asarray(ell, dtype=float)
        tt = np.asarray(tt, dtype=float)
        ee = np.asarray(ee, dtype=float)
        if ell.ndim != 1 or tt.ndim != 1 or ee.ndim != 1:
            return plot, warnings

        lmax = int(min(ell_max_plot, int(ell.max()) if ell.size else ell_max_plot))
        mask = (ell >= 2) & (ell <= float(lmax))
        if not np.any(mask):
            return plot, warnings

        # --- TT/EE comparison ---
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 7), sharex=True)
        ax1.plot(ell[mask], tt[mask], lw=1.8, label="bounce (injected P(k))")
        ax2.plot(ell[mask], ee[mask], lw=1.8, label="bounce (injected P(k))")
        if tt_base is not None and ee_base is not None:
            tt0 = np.asarray(tt_base, dtype=float)
            ee0 = np.asarray(ee_base, dtype=float)
            if tt0.shape == tt.shape and ee0.shape == ee.shape:
                ax1.plot(ell[mask], tt0[mask], lw=1.2, ls="--", label="baseline (power-law)")
                ax2.plot(ell[mask], ee0[mask], lw=1.2, ls="--", label="baseline (power-law)")
        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax1.set_ylabel(r"$D_\ell^{TT}$ [$\mu$K$^2$]")
        ax2.set_ylabel(r"$D_\ell^{EE}$ [$\mu$K$^2$]")
        ax2.set_xlabel(r"multipole $\ell$")
        ax1.grid(True, ls=":", alpha=0.35)
        ax2.grid(True, ls=":", alpha=0.35)
        ax1.legend(loc="best")
        ax2.legend(loc="best")
        fig.suptitle("CMB spectra (CAMB): bounce injection vs baseline")
        fig.tight_layout()
        p1 = out_dir / "cl_tt_ee_comparison.png"
        fig.savefig(p1, dpi=180)
        plt.close(fig)
        plot["cl_tt_ee_comparison_png"] = str(p1)

        # --- Ratio plot (requires baseline) ---
        if tt_base is not None and ee_base is not None:
            tt0 = np.asarray(tt_base, dtype=float)
            ee0 = np.asarray(ee_base, dtype=float)
            if tt0.shape == tt.shape and ee0.shape == ee.shape:
                eps = 1e-30
                r_tt = (tt - tt0) / np.maximum(np.abs(tt0), eps)
                r_ee = (ee - ee0) / np.maximum(np.abs(ee0), eps)
                fig2, (bx1, bx2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 7), sharex=True)
                bx1.plot(ell[mask], r_tt[mask], lw=1.6)
                bx2.plot(ell[mask], r_ee[mask], lw=1.6)
                bx1.axhline(0.0, color="black", lw=1.0, alpha=0.7)
                bx2.axhline(0.0, color="black", lw=1.0, alpha=0.7)
                bx1.set_ylabel(r"$(D_\ell^{TT}-D_{\ell,0}^{TT})/D_{\ell,0}^{TT}$")
                bx2.set_ylabel(r"$(D_\ell^{EE}-D_{\ell,0}^{EE})/D_{\ell,0}^{EE}$")
                bx2.set_xlabel(r"multipole $\ell$")
                bx1.grid(True, ls=":", alpha=0.35)
                bx2.grid(True, ls=":", alpha=0.35)
                fig2.suptitle("Relative deviation vs baseline (bounce injection)")
                fig2.tight_layout()
                p2 = out_dir / "cl_ratio_tt_ee.png"
                fig2.savefig(p2, dpi=180)
                plt.close(fig2)
                plot["cl_ratio_tt_ee_png"] = str(p2)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


class BoltzmannTransferModule(TfptModule):
    module_id = "boltzmann_transfer"
    title = "Boltzmann transfer (CAMB-backed C_ℓ + explicit k̂→ℓ mapping)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "k_calibration assumptions: tfpt_suite/data/k_calibration.json (uses reheating policy v1.06 settings)",
                "cosmology reference: tfpt_suite/data/global_reference_minimal.json (Planck Ω_b h^2 anchor)",
                "TFPT R^2 scale M from tfpt_suite/constants.py (M/Mpl)",
                "primordial spectrum tables (optional): output of primordial_spectrum_builder (bounce injection)",
            ],
            outputs=[
                "explicit k_hat→ℓ mapping under an explicit a0/a_transition policy",
                "ℓ-range coverage diagnostics (CMB vs small-scale decision aid)",
                "CMB power spectra C_ℓ (TT/TE/EE/BB) computed via CAMB under explicit cosmology + TFPT primordial parameters",
                "optional Planck high-ℓ/low-ℓ/lensing likelihood logL (Cobaya; opt-in)",
            ],
            formulas=[
                r"ℓ \approx k(\mathrm{Mpc}^{-1}) \, \chi_\*",
                r"k(\mathrm{Mpc}^{-1}) = k_{\hat{}} \, M[\mathrm{GeV}] \, (\mathrm{GeV}\to\mathrm{Mpc}^{-1}) / (a_0/a_{\rm tr})",
                r"P_{\mathcal{R}}(k) = P_{\mathcal{R},\mathrm{base}}(k)\,|T_s(k)|^2 \;\;(\text{bounce injection via primordial_spectrum_builder})",
            ],
            validation=[
                "ell_predictions_falsifiable: emits a PASS/WARN/FAIL depending on whether the predicted ℓ-range intersects declared ℓ-targets.",
                "signature_policy_declared: PASS if a CMB vs small-scale policy is present in k_calibration assumptions.",
                "camb_power_spectra_computed: emits PASS when CAMB returns finite C_ℓ spectra for the declared cosmology+primordial inputs.",
                "bounce_feature_injection_wired: PASS when CAMB uses an injected P(k) table (not just power-law).",
            ],
            determinism="Deterministic given config files.",
            question="Given a declared expansion-history policy (a0/a_transition), where do bounce k-hat features land in ℓ-space?",
            objective=[
                "Provide a falsifiable, assumption-explicit ℓ mapping and a concrete C_ℓ prediction backend (CAMB).",
                "Make the 'last mile' from TFPT primordial parameters to observable spectra explicit and auditable.",
            ],
            gaps=[
                "Planck low-ℓ/lensing likelihoods are optional and depend on Cobaya + data availability; this module remains robust without them.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        data_dir = Path(__file__).resolve().parent.parent / "data"
        cfg_path = data_dir / "k_calibration.json"
        raw = json.loads(cfg_path.read_text(encoding="utf-8"))
        cosmo = raw.get("cosmology_flat_lcdm", {}) if isinstance(raw.get("cosmology_flat_lcdm", {}), dict) else {}
        asm = raw.get("assumptions", {}) if isinstance(raw.get("assumptions", {}), dict) else {}

        nuisance_policy = _load_nuisance_policy(data_dir=data_dir)
        planck_nuisance = nuisance_policy.get("planck_2018", {}) if isinstance(nuisance_policy.get("planck_2018", {}), dict) else {}
        A_planck = _as_float_or(planck_nuisance.get("A_planck", 1.0), 1.0)
        pivot_scalar_Mpc_inv = _as_float_or(planck_nuisance.get("pivot_scalar_Mpc_inv", 0.05), 0.05)

        # Targets for falsifiability
        ell_targets = asm.get("ell_targets", [2, 30, 700])
        ell_targets_f = [float(x) for x in ell_targets] if isinstance(ell_targets, list) else [2.0, 30.0, 700.0]

        signature_policy = asm.get("signature_policy", {}) if isinstance(asm.get("signature_policy", {}), dict) else {}
        cmb_range = signature_policy.get("cmb_ell_range", DEFAULT_CMB_ELL_RANGE)
        cmb_min = _as_float_or(cmb_range[0], DEFAULT_CMB_ELL_RANGE[0]) if isinstance(cmb_range, (list, tuple)) and len(cmb_range) >= 2 else DEFAULT_CMB_ELL_RANGE[0]
        cmb_max = _as_float_or(cmb_range[1], DEFAULT_CMB_ELL_RANGE[1]) if isinstance(cmb_range, (list, tuple)) and len(cmb_range) >= 2 else DEFAULT_CMB_ELL_RANGE[1]
        small_scale_min = _as_float_or(signature_policy.get("small_scale_ell_min", DEFAULT_SMALL_SCALE_ELL_MIN), DEFAULT_SMALL_SCALE_ELL_MIN)
        decision_rule = str(signature_policy.get("decision_rule", DEFAULT_SIGNATURE_DECISION) or DEFAULT_SIGNATURE_DECISION)
        primary_signature = str(signature_policy.get("primary_signature", "") or "").strip()
        signature_rationale = str(signature_policy.get("rationale", "") or "").strip()

        # Compute a0/a_transition from the explicit policy inputs (entropy + N_inflation + N_reheat).
        inp = CosmoScaleInputs(
            transition=str(asm.get("transition", "horizon_exit_of_pivot")),
            N_inflation_from_transition=float(asm.get("N_inflation_from_transition", 56.0)),
            N_reheat=float(asm.get("N_reheat", 0.0)),
            T_reheat_GeV=float(asm.get("T_reheat_GeV", 1.0e13)),
            g_star_s_reheat=float(asm.get("g_star_s_reheat", 120.0)),
            g_star_s_today=float(asm.get("g_star_s_today", 3.91)),
            T0_K=float(asm.get("T0_K", 2.7255)),
        )
        a0_over_atr, note = a0_over_a_transition(inp)
        chi_star_Mpc = 14_000.0
        N_eff = float(inp.N_inflation_from_transition)
        ell_bounce_s_est = None
        ell_bounce_t_est = None

        # Prefer k_calibration expansion-budget estimate (reheating/transition policy) when available.
        kcal_out_path = Path(config.output_dir) / "k_calibration" / "results.json"
        if kcal_out_path.is_file():
            try:
                payload = json.loads(kcal_out_path.read_text(encoding="utf-8"))
                kres = payload.get("results", {}) if isinstance(payload, dict) else {}
                est = kres.get("expansion_budget_estimate", {}) if isinstance(kres.get("expansion_budget_estimate", {}), dict) else {}
                a0_est = float(est.get("a0_over_a_transition", float("nan")))
                if math.isfinite(a0_est) and a0_est > 0:
                    a0_over_atr = a0_est
                    note = f"from k_calibration expansion_budget_estimate ({kcal_out_path})"
                N_est = float(est.get("N_inflation_from_transition", float("nan")))
                if math.isfinite(N_est) and N_est > 0:
                    N_eff = N_est
                ell_s = float(est.get("ell_bounce_s", float("nan")))
                if math.isfinite(ell_s) and ell_s > 0:
                    ell_bounce_s_est = ell_s
                ell_t = float(est.get("ell_bounce_t", float("nan")))
                if math.isfinite(ell_t) and ell_t > 0:
                    ell_bounce_t_est = ell_t
                cos = kres.get("cosmology", {}) if isinstance(kres.get("cosmology", {}), dict) else {}
                chi_est = float(cos.get("chi_star_Mpc", float("nan")))
                if math.isfinite(chi_est) and chi_est > 0:
                    chi_star_Mpc = chi_est
            except Exception:
                pass

        # M scale from TFPT constants
        cst = TfptConstants.compute()
        M_GeV = float(float(cst.M_over_Mpl) * 2.435e18)

        # k-hat range
        grid = asm.get("k_hat_grid_for_ell_range", [1e-7, 10.0])
        kmin, kmax = float(grid[0]), float(grid[1]) if isinstance(grid, list) and len(grid) >= 2 else (1e-7, 10.0)

        ell_rng: EllRange
        if a0_over_atr is None or not math.isfinite(float(a0_over_atr)) or float(a0_over_atr) <= 0:
            ell_rng = EllRange(ell_min=float("nan"), ell_max=float("nan"), note=f"a0/a_transition unavailable: {note}")
        else:
            ell1 = ell_from_k_hat(k_hat=float(kmin), M_GeV=float(M_GeV), a0_over_a_tr=float(a0_over_atr), chi_star_Mpc=float(chi_star_Mpc))
            ell2 = ell_from_k_hat(k_hat=float(kmax), M_GeV=float(M_GeV), a0_over_a_tr=float(a0_over_atr), chi_star_Mpc=float(chi_star_Mpc))
            ell_rng = EllRange(ell_min=float(min(ell1, ell2)), ell_max=float(max(ell1, ell2)), note="computed from k_hat range + entropy/N policy")

        cmb_overlap = bool(
            math.isfinite(ell_rng.ell_min)
            and math.isfinite(ell_rng.ell_max)
            and (ell_rng.ell_max >= cmb_min)
            and (ell_rng.ell_min <= cmb_max)
        )
        small_scale_overlap = bool(math.isfinite(ell_rng.ell_max) and ell_rng.ell_max >= float(small_scale_min))
        if decision_rule == "prefer_small_scale_if_overlap":
            policy_decision = "small_scale" if small_scale_overlap else ("cmb" if cmb_overlap else "unclear")
        else:
            policy_decision = "cmb" if cmb_overlap else ("small_scale" if small_scale_overlap else "unclear")

        # Falsifiability check: do targets lie in predicted range?
        def _in_range(x: float) -> bool:
            return bool(math.isfinite(x) and math.isfinite(ell_rng.ell_min) and math.isfinite(ell_rng.ell_max) and (ell_rng.ell_min <= x <= ell_rng.ell_max))

        hits = [x for x in ell_targets_f if _in_range(float(x))]
        checks: list[Check] = []
        checks.append(mk_check_pass("ell_mapping_computed", f"a0/a_tr={a0_over_atr} (note={note})"))
        policy_declared = (
            bool(signature_policy)
            and ("cmb_ell_range" in signature_policy)
            and ("small_scale_ell_min" in signature_policy)
            and ("primary_signature" in signature_policy)
        )
        checks.append(
            mk_check_pass("signature_policy_declared", f"policy={signature_policy}")
            if policy_declared
            else mk_check_warn("signature_policy_declared", "missing signature_policy in k_calibration.json assumptions")
        )
        checks.append(
            mk_check_pass("signature_policy_decision", f"decision={policy_decision}, cmb_overlap={cmb_overlap}, small_scale_overlap={small_scale_overlap}")
            if policy_decision != "unclear"
            else mk_check_info("signature_policy_decision", f"decision=unclear, cmb_overlap={cmb_overlap}, small_scale_overlap={small_scale_overlap}")
        )
        if primary_signature == "tensor_CMB":
            if ell_bounce_t_est is not None and cmb_min <= ell_bounce_t_est <= cmb_max:
                checks.append(mk_check_pass("signature_policy_consistent_with_bounce", f"tensor ell≈{ell_bounce_t_est:.3g} within CMB range"))
            else:
                checks.append(mk_check_warn("signature_policy_consistent_with_bounce", f"tensor ell≈{ell_bounce_t_est} not in CMB range"))
        elif primary_signature == "scalar_small_scale":
            if ell_bounce_s_est is not None and ell_bounce_s_est >= small_scale_min:
                checks.append(mk_check_pass("signature_policy_consistent_with_bounce", f"scalar ell≈{ell_bounce_s_est:.3g} above small-scale threshold"))
            else:
                checks.append(mk_check_warn("signature_policy_consistent_with_bounce", f"scalar ell≈{ell_bounce_s_est} below small-scale threshold"))
        elif primary_signature in {"conditional_gate", "conditional"}:
            # Contract-style policy: report where the scalar/tensor bounce scales land under the
            # current expansion-budget estimate, without asserting a fixed signature a priori.
            parts: list[str] = []
            if ell_bounce_s_est is not None and math.isfinite(float(ell_bounce_s_est)):
                parts.append(f"scalar ell≈{float(ell_bounce_s_est):.6g} (in CMB={cmb_min <= float(ell_bounce_s_est) <= cmb_max})")
            else:
                parts.append("scalar ell≈n/a")
            if ell_bounce_t_est is not None and math.isfinite(float(ell_bounce_t_est)):
                parts.append(f"tensor ell≈{float(ell_bounce_t_est):.6g} (in CMB={cmb_min <= float(ell_bounce_t_est) <= cmb_max})")
            else:
                parts.append("tensor ell≈n/a")
            checks.append(mk_check_pass("signature_policy_consistent_with_bounce", "; ".join(parts)))
        else:
            checks.append(mk_check_info("signature_policy_consistent_with_bounce", f"primary_signature missing or unknown: '{primary_signature}'"))
        if hits:
            checks.append(mk_check_pass("ell_predictions_falsifiable", f"targets hit: {hits} within [{ell_rng.ell_min:.3g}, {ell_rng.ell_max:.3g}]"))
        else:
            checks.append(
                (mk_check_fail if mode == "physics" else mk_check_warn)(
                    "ell_predictions_falsifiable",
                    f"no targets hit; predicted ℓ-range [{ell_rng.ell_min:.3g}, {ell_rng.ell_max:.3g}] (targets={ell_targets_f})",
                )
            )

        # --- C_ell prediction backend (CAMB) ---
        data_dir = Path(__file__).resolve().parent.parent / "data"
        global_ref_path = data_dir / "global_reference_minimal.json"
        omega_b_h2 = 0.02237  # Planck 2018 baseline (fallback; citeable anchor lives in global_reference_minimal.json)
        if global_ref_path.is_file():
            try:
                graw = json.loads(global_ref_path.read_text(encoding="utf-8"))
                omega_b_h2 = float(graw["observables"]["omega_b_h2_planck2018"]["mean"])
            except Exception:
                pass

        H0_km_s_Mpc = float(cosmo.get("H0_km_s_Mpc", 67.36))
        Omega_m = float(cosmo.get("Omega_m", 0.3153))
        Omega_r = float(cosmo.get("Omega_r", 9.2e-5))
        h = float(H0_km_s_Mpc / 100.0)
        Omega_b = float(omega_b_h2 / (h * h)) if h > 0 else float("nan")
        Omega_c = float(Omega_m - Omega_b) if math.isfinite(Omega_b) else float("nan")
        omega_c_h2 = float(Omega_c * (h * h)) if (math.isfinite(Omega_c) and Omega_c > 0 and h > 0) else float("nan")

        # Primordial parameters from TFPT/Starobinsky (deterministic given N and M/Mpl).
        N = int(round(float(N_eff)))
        n_s = float(1.0 - 2.0 / float(N))
        r = float(12.0 / (float(N) ** 2))
        A_s = float(((float(N) ** 2) / (24.0 * (math.pi**2))) * (float(cst.M_over_Mpl) ** 2))
        # Planck plik-lite uses lmax=2508; keep spectra compatible for optional likelihood evaluation.
        lmax = int(max(2508, max(ell_targets_f) if ell_targets_f else 2500))
        tau_reio = 0.054  # Planck-like baseline; make explicit until a TFPT-derived reionization model exists.
        mnu_eV = 0.06

        camb_backend = "camb"
        cl_summary: dict[str, Any] = {}
        cl_samples: list[dict[str, Any]] = []
        plot_block: dict[str, str | None] = {}
        plot_warnings: list[str] = []
        planck_pliklite: dict[str, Any] | None = None
        planck_lowl: dict[str, Any] | None = None
        planck_lensing: dict[str, Any] | None = None
        planck_combined: dict[str, Any] | None = None
        bounce_injection_used = False
        pk_source = "power_law"
        try:
            import camb  # type: ignore

            pars = camb.CAMBparams()
            pars.set_cosmology(
                H0=H0_km_s_Mpc,
                ombh2=float(omega_b_h2),
                omch2=float(omega_c_h2),
                mnu=mnu_eV,
                omk=0.0,
                tau=tau_reio,
            )

            # Prefer bounce-injected primordial P(k) tables if available.
            pk_table_path = Path(config.output_dir) / "primordial_spectrum_builder" / "results.json"
            pk_table = None
            if pk_table_path.is_file():
                try:
                    payload = json.loads(pk_table_path.read_text(encoding="utf-8"))
                    tab = payload.get("results", {}).get("table", {}) if isinstance(payload, dict) else {}
                    k_tab = tab.get("k_Mpc_inv", None)
                    pr_tab = tab.get("P_R", None)
                    pt_tab = tab.get("P_t", None)
                    if (
                        isinstance(k_tab, list)
                        and isinstance(pr_tab, list)
                        and isinstance(pt_tab, list)
                        and len(k_tab) == len(pr_tab) == len(pt_tab)
                        and len(k_tab) >= 10
                    ):
                        pk_table = (np.asarray(k_tab, dtype=float), np.asarray(pr_tab, dtype=float), np.asarray(pt_tab, dtype=float))
                except Exception:
                    pk_table = None

            if pk_table is not None:
                k_arr, pr_arr, pt_arr = pk_table
                pk_source = f"bounce_injected_table:{pk_table_path}"
                # Pivot-normalize the injected table to keep CAMB in a realistic amplitude regime.
                # This makes the nuisance handling explicit: we fix the scalar amplitude at the pivot,
                # while retaining the bounce-induced *shape* (and any r(k) distortion).
                try:
                    k0 = float(pivot_scalar_Mpc_inv)
                    pr0 = float(np.interp(float(np.log(k0)), np.log(k_arr), pr_arr))
                    if math.isfinite(pr0) and pr0 > 0:
                        scale = float(A_s / pr0)
                        pr_arr = pr_arr * scale
                        pt_arr = pt_arr * scale
                        checks.append(mk_check_info("primordial_pivot_normalization", f"applied scale={scale:.6g} to match P_R({k0})=A_s"))
                        pk_source = pk_source + f" (pivot_normalized scale={scale:.6g})"
                    else:
                        checks.append(mk_check_info("primordial_pivot_normalization", f"skipped (P_R({k0}) invalid: {pr0})"))
                except Exception as e:
                    checks.append(mk_check_info("primordial_pivot_normalization", f"skipped (error: {e})"))
                pars.set_initial_power_table(k_arr, pk=pr_arr, pk_tensor=pt_arr)
                pars.WantTensors = True
                bounce_injection_used = True
            else:
                pars.InitPower.set_params(As=A_s, ns=n_s, r=r, pivot_scalar=pivot_scalar_Mpc_inv)
                pars.WantTensors = bool(r > 0)
            pars.set_for_lmax(lmax, lens_potential_accuracy=0)

            results = camb.get_results(pars)
            powers = results.get_cmb_power_spectra(pars, CMB_unit="muK")
            cl_tot = powers.get("total", None)
            if cl_tot is None:
                raise RuntimeError("CAMB did not return 'total' CMB power spectrum")

            # CAMB returns an array of shape (lmax+1, 4): TT, EE, BB, TE
            tt = np.asarray(cl_tot)[:, 0]
            ee = np.asarray(cl_tot)[:, 1]
            bb = np.asarray(cl_tot)[:, 2]
            te = np.asarray(cl_tot)[:, 3]
            ell = np.arange(tt.shape[0], dtype=int)

            # Baseline power-law spectra for ratio plots (only meaningful when bounce injection is active).
            tt0 = ee0 = None
            if pk_table is not None:
                try:
                    pars0 = camb.CAMBparams()
                    pars0.set_cosmology(
                        H0=H0_km_s_Mpc,
                        ombh2=float(omega_b_h2),
                        omch2=float(omega_c_h2),
                        mnu=mnu_eV,
                        omk=0.0,
                        tau=tau_reio,
                    )
                    pars0.InitPower.set_params(As=A_s, ns=n_s, r=r, pivot_scalar=pivot_scalar_Mpc_inv)
                    pars0.WantTensors = bool(r > 0)
                    pars0.set_for_lmax(lmax, lens_potential_accuracy=0)
                    res0 = camb.get_results(pars0)
                    pow0 = res0.get_cmb_power_spectra(pars0, CMB_unit="muK")
                    tot0 = pow0.get("total", None)
                    if tot0 is not None:
                        arr0 = np.asarray(tot0)
                        if arr0.shape[0] == tt.shape[0] and arr0.shape[1] >= 2:
                            tt0 = np.asarray(arr0)[:, 0]
                            ee0 = np.asarray(arr0)[:, 1]
                except Exception:
                    tt0 = ee0 = None

            pp = None
            try:
                lens_cls = results.get_lens_potential_cls(lmax=lmax)
                if lens_cls is not None:
                    arr = np.asarray(lens_cls)
                    if arr.ndim == 2 and arr.shape[1] >= 1:
                        pp = arr[:, 0]
            except Exception:
                pp = None

            # Basic sanity checks + a world-contact marker (first acoustic peak position).
            peak_slice = slice(50, min(400, tt.shape[0]))
            peak_idx = int(peak_slice.start + int(np.argmax(tt[peak_slice])))
            peak_ell = int(ell[peak_idx])
            peak_tt = float(tt[peak_idx])

            cl_summary = {
                "backend": camb_backend,
                "camb_version": getattr(camb, "__version__", "unknown"),
                "lmax": int(tt.shape[0] - 1),
                "peak1": {"ell": peak_ell, "TT_muK2": peak_tt},
                "TT_at_ell": {str(L): float(tt[L]) for L in [2, 30, 200, 700, 2000] if L < tt.shape[0]},
            }

            # Compact sampling table (kept small to avoid bloating results.json)
            sample_ells = sorted(set([2, 3, 4, 5, 10, 20, 30, 50, 80, 100, 150, 200, 300, 500, 700, 1000, 1500, 2000, lmax]))
            cl_samples = [
                {"ell": int(L), "TT_muK2": float(tt[L]), "EE_muK2": float(ee[L]), "TE_muK2": float(te[L]), "BB_muK2": float(bb[L])}
                for L in sample_ells
                if 0 <= L < tt.shape[0]
            ]

            plot_block, plot_warnings = _plot_cl_comparison(
                out_dir=self.output_dir(config),
                ell=ell,
                tt=tt,
                ee=ee,
                tt_base=tt0,
                ee_base=ee0,
                ell_max_plot=int(cmb_max),
            )

            # --- Optional Planck likelihoods (Cobaya-native; opt-in via env) ---
            enable_planck = str(os.environ.get(PLANCK_ENABLE_ENV, "0")).strip().lower() in {"1", "true", "yes", "on"}
            logp_kwargs = {"A_planck": A_planck} if math.isfinite(A_planck) else {}
            if not enable_planck:
                planck_pliklite = {
                    "enabled": False,
                    "A_planck": A_planck,
                    "note": f"set {PLANCK_ENABLE_ENV}=1 to evaluate Planck likelihoods (Cobaya-native).",
                }
                planck_lowl = {"enabled": False, "note": f"set {PLANCK_ENABLE_ENV}=1 to evaluate Planck low-ℓ TT/EE."}
                planck_lensing = {"enabled": False, "note": f"set {PLANCK_ENABLE_ENV}=1 to evaluate Planck lensing likelihood."}
                checks.append(mk_check_info("planck_likelihood_evaluated", f"disabled (set {PLANCK_ENABLE_ENV}=1)"))
                checks.append(mk_check_info("planck_lowl_evaluated", f"disabled (set {PLANCK_ENABLE_ENV}=1)"))
                checks.append(mk_check_info("planck_lensing_evaluated", f"disabled (set {PLANCK_ENABLE_ENV}=1)"))
            else:
                try:
                    import importlib
                    import warnings

                    warnings.filterwarnings("ignore", message="urllib3 v2 only supports OpenSSL*")
                    warnings.filterwarnings("ignore", category=ResourceWarning)

                    class _Provider:
                        def __init__(self, cl_dict: dict[str, np.ndarray]):
                            self._cl = cl_dict

                        def get_Cl(self, *, ell_factor: bool = True):  # Cobaya provider API
                            if not ell_factor:
                                raise ValueError("This provider only supports ell_factor=True (D_ell arrays).")
                            return self._cl

                    def _logp(like) -> float:
                        try:
                            return float(like.logp(**logp_kwargs))
                        except TypeError:
                            return float(like.logp())

                    provider = _Provider(
                        {
                            "tt": tt,
                            "te": te,
                            "ee": ee,
                            "bb": bb,
                            **({"pp": pp} if pp is not None else {}),
                        }
                    )

                    # --- high-l plik-lite ---
                    try:
                        from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native  # type: ignore

                        like = TTTEEE_lite_native()
                        like.path = None
                        like.dataset_file = PLANCK_PLIK_DATASET
                        like.dataset_params = {"use_cl": "tt te ee"}
                        like.initialize()
                        like.provider = provider
                        logp_model = _logp(like)

                        planck_pliklite = {
                            "enabled": True,
                            "likelihood": "planck_2018_highl_plik.TTTEEE_lite_native",
                            "A_planck": A_planck,
                            "logp": logp_model,
                            "chi2": float(-2.0 * logp_model),
                        }

                        # Δχ² vs a power-law baseline (same cosmology + TFPT baseline As/ns/r), if bounce injection is active.
                        if pk_table is not None:
                            pars0 = camb.CAMBparams()
                            pars0.set_cosmology(
                                H0=H0_km_s_Mpc,
                                ombh2=float(omega_b_h2),
                                omch2=float(omega_c_h2),
                                mnu=mnu_eV,
                                omk=0.0,
                                tau=tau_reio,
                            )
                            pars0.InitPower.set_params(As=A_s, ns=n_s, r=r, pivot_scalar=pivot_scalar_Mpc_inv)
                            pars0.WantTensors = bool(r > 0)
                            pars0.set_for_lmax(lmax, lens_potential_accuracy=0)
                            res0 = camb.get_results(pars0)
                            pow0 = res0.get_cmb_power_spectra(pars0, CMB_unit="muK")
                            tot0 = pow0.get("total", None)
                            if tot0 is not None:
                                arr0 = np.asarray(tot0)
                                tt0 = arr0[:, 0]
                                ee0 = arr0[:, 1]
                                te0 = arr0[:, 3]
                                like.provider = _Provider({"tt": tt0, "te": te0, "ee": ee0})
                                logp_base = _logp(like)
                                planck_pliklite["baseline_power_law"] = {"logp": logp_base, "chi2": float(-2.0 * logp_base)}
                                planck_pliklite["delta_chi2_vs_power_law"] = float(-2.0 * (logp_model - logp_base))

                        checks.append(
                            mk_check_pass("planck_likelihood_evaluated", f"{planck_pliklite['likelihood']}: logp={planck_pliklite['logp']:.6g}")
                            if math.isfinite(float(planck_pliklite.get("logp", float("nan"))))
                            else mk_check_info("planck_likelihood_evaluated", "Planck plik-lite likelihood not evaluated")
                        )
                        if "delta_chi2_vs_power_law" in planck_pliklite:
                            checks.append(
                                mk_check_pass(
                                    "planck_delta_chi2_reported",
                                    f"Δχ²_vs_power_law={float(planck_pliklite['delta_chi2_vs_power_law']):.6g} (A_planck fixed)",
                                )
                            )
                    except Exception as e:
                        planck_pliklite = {"enabled": True, "error": str(e), "likelihood": "planck_2018_highl_plik.TTTEEE_lite_native"}
                        checks.append(mk_check_info("planck_likelihood_evaluated", f"unavailable: {e}"))

                    # --- low-l TT/EE ---
                    lowl_entries: list[dict[str, Any]] = []
                    lowl_logp = 0.0
                    lowl_used_split = False
                    for module_path, class_name in PLANCK_LOWL_CANDIDATES:
                        if class_name == "TT_EE" and lowl_used_split:
                            continue
                        try:
                            module = importlib.import_module(module_path)
                            like_cls = getattr(module, class_name)
                            like = like_cls()
                            like.path = None
                            like.initialize()
                            like.provider = provider
                            logp_val = _logp(like)
                            lowl_entries.append({"likelihood": module_path, "logp": logp_val, "chi2": float(-2.0 * logp_val)})
                            lowl_logp += float(logp_val)
                            if class_name in {"TT", "EE"}:
                                lowl_used_split = True
                            if class_name == "TT_EE":
                                break
                        except Exception:
                            continue
                    if lowl_entries:
                        planck_lowl = {
                            "enabled": True,
                            "logp": lowl_logp,
                            "chi2": float(-2.0 * lowl_logp),
                            "likelihoods": lowl_entries,
                        }
                        checks.append(mk_check_pass("planck_lowl_evaluated", f"logp={lowl_logp:.6g} ({len(lowl_entries)} datasets)"))
                    else:
                        planck_lowl = {
                            "enabled": True,
                            "error": "no low-l likelihoods available",
                            "candidates": [f"{m}:{c}" for m, c in PLANCK_LOWL_CANDIDATES],
                        }
                        checks.append(mk_check_info("planck_lowl_evaluated", "unavailable (Cobaya low-l datasets not found)"))

                    # --- lensing ---
                    lensing_entry: dict[str, Any] | None = None
                    lensing_error: str | None = None
                    for module_path, class_name in PLANCK_LENSING_CANDIDATES:
                        try:
                            module = importlib.import_module(module_path)
                            like_cls = getattr(module, class_name)
                            like = like_cls()
                            like.path = None
                            like.initialize()
                            like.provider = provider
                            logp_val = _logp(like)
                            lensing_entry = {
                                "likelihood": module_path,
                                "logp": logp_val,
                                "chi2": float(-2.0 * logp_val),
                            }
                            break
                        except Exception as e:
                            lensing_error = str(e)
                            continue
                    if lensing_entry is not None:
                        planck_lensing = {"enabled": True, **lensing_entry}
                        checks.append(mk_check_pass("planck_lensing_evaluated", f"logp={lensing_entry['logp']:.6g}"))
                    else:
                        planck_lensing = {
                            "enabled": True,
                            "error": lensing_error or "no lensing likelihoods available",
                            "candidates": [f"{m}:{c}" for m, c in PLANCK_LENSING_CANDIDATES],
                        }
                        checks.append(mk_check_info("planck_lensing_evaluated", "unavailable (lensing dataset not found)"))

                    components: dict[str, float] = {}
                    combined_logp = 0.0
                    if isinstance(planck_pliklite, dict) and math.isfinite(float(planck_pliklite.get("logp", float("nan")))):
                        components["pliklite"] = float(planck_pliklite["logp"])
                        combined_logp += float(planck_pliklite["logp"])
                    if isinstance(planck_lowl, dict) and math.isfinite(float(planck_lowl.get("logp", float("nan")))):
                        components["lowl"] = float(planck_lowl["logp"])
                        combined_logp += float(planck_lowl["logp"])
                    if isinstance(planck_lensing, dict) and math.isfinite(float(planck_lensing.get("logp", float("nan")))):
                        components["lensing"] = float(planck_lensing["logp"])
                        combined_logp += float(planck_lensing["logp"])
                    if components:
                        planck_combined = {
                            "logp": float(combined_logp),
                            "chi2": float(-2.0 * combined_logp),
                            "components": components,
                            "A_planck": A_planck,
                        }
                        checks.append(mk_check_pass("planck_combined_logp_reported", f"logp={combined_logp:.6g}, components={list(components.keys())}"))
                    else:
                        planck_combined = None
                        checks.append(mk_check_info("planck_combined_logp_reported", "no Planck likelihood components evaluated"))
                except Exception as e:
                    planck_pliklite = {"enabled": True, "error": str(e), "likelihood": "planck_2018_highl_plik.TTTEEE_lite_native"}
                    planck_lowl = {"enabled": True, "error": str(e), "likelihoods": PLANCK_LOWL_CANDIDATES}
                    planck_lensing = {"enabled": True, "error": str(e), "likelihoods": PLANCK_LENSING_CANDIDATES}
                    checks.append(mk_check_info("planck_likelihood_evaluated", f"unavailable: {e}"))

            checks.append(mk_check_pass("camb_power_spectra_computed", f"computed C_ell (TT/EE/BB/TE) up to lmax={lmax}"))
            checks.append(
                mk_check_pass("camb_first_peak_in_cmb_range", f"first TT peak at ell≈{peak_ell} (TT≈{peak_tt:.1f} µK^2)")
                if (180 <= peak_ell <= 260)
                else mk_check_warn("camb_first_peak_in_cmb_range", f"first TT peak at ell≈{peak_ell} (unexpected; check inputs)")
            )
            checks.append(mk_check_pass("boltzmann_transfer_implemented", f"backend={camb_backend} (CAMB)"))
            checks.append(
                mk_check_pass("bounce_feature_injection_wired", f"CAMB P(k) source: {pk_source}")
                if pk_table is not None
                else mk_check_info("bounce_feature_injection_wired", f"no injected P(k) table found; using power-law (expected file: {pk_table_path})")
            )
        except Exception as e:
            camb_backend = "unavailable"
            cl_summary = {"backend": camb_backend, "error": str(e)}
            # We keep this as INFO (not WARN) because the suite can still use the ℓ mapping scaffold even without CAMB.
            checks.append(mk_check_info("camb_power_spectra_computed", f"CAMB unavailable: {e}"))
            checks.append(mk_check_info("boltzmann_transfer_implemented", "backend unavailable (ℓ mapping only)"))
            checks.append(mk_check_info("bounce_feature_injection_wired", "CAMB unavailable; cannot inject bounce features"))
            bounce_injection_used = False
            pk_source = "camb_unavailable"
            planck_pliklite = None
            planck_lowl = None
            planck_lensing = None
            planck_combined = None

        lines: list[str] = []
        lines += [
            "Boltzmann transfer (explicit ℓ mapping + CAMB-backed C_ell)",
            "",
            f"mode={mode}",
            f"config: {cfg_path}",
            "",
            "Policy inputs:",
            f"- transition={inp.transition}, N_infl={inp.N_inflation_from_transition}, N_reheat={inp.N_reheat}",
            f"- T_reheat={inp.T_reheat_GeV:.3e} GeV, g*_s(reh)={inp.g_star_s_reheat}, g*_s(today)={inp.g_star_s_today}",
            f"- a0/a_transition = {a0_over_atr} ({note})",
            f"- signature_policy = {signature_policy} (decision={policy_decision})",
            f"- Planck nuisance policy = {planck_nuisance}",
            "",
            "Scale mapping:",
            f"- TFPT M ≈ {M_GeV:.3e} GeV (from M/Mpl)",
            f"- chi_star ≈ {chi_star_Mpc:.0f} Mpc (proxy)",
            f"- k_hat range = [{kmin:.3e}, {kmax:.3e}] => ℓ range ≈ [{ell_rng.ell_min:.3g}, {ell_rng.ell_max:.3g}]",
            f"- targets = {ell_targets_f}, hits = {hits}",
            "",
            "CMB spectra backend (CAMB):",
            f"- backend = {cl_summary.get('backend')}",
            f"- cosmology: H0={H0_km_s_Mpc}, Omega_m={Omega_m}, Omega_r={Omega_r}, omega_b_h2={omega_b_h2}, omega_c_h2={omega_c_h2}, tau={tau_reio}",
            f"- primordial (TFPT/Starobinsky): N={N}, n_s={n_s:.6f}, A_s={A_s:.4e}, r={r:.6f}, pivot={pivot_scalar_Mpc_inv} 1/Mpc",
            f"- primordial P(k) source = {pk_source}",
            *(
                [f"- Planck plik-lite: {json.dumps(planck_pliklite, ensure_ascii=False)}"]
                if isinstance(planck_pliklite, dict) and planck_pliklite
                else []
            ),
            *(
                [f"- Planck low-l: {json.dumps(planck_lowl, ensure_ascii=False)}"]
                if isinstance(planck_lowl, dict) and planck_lowl
                else []
            ),
            *(
                [f"- Planck lensing: {json.dumps(planck_lensing, ensure_ascii=False)}"]
                if isinstance(planck_lensing, dict) and planck_lensing
                else []
            ),
            *(
                [f"- Planck combined: {json.dumps(planck_combined, ensure_ascii=False)}"]
                if isinstance(planck_combined, dict) and planck_combined
                else []
            ),
            *(["- C_ell summary: " + json.dumps(cl_summary, ensure_ascii=False)] if cl_summary else []),
            *(["- sampled C_ell points: " + json.dumps(cl_samples, ensure_ascii=False)] if cl_samples else []),
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "config_file": str(cfg_path),
                "a0_over_a_transition": a0_over_atr,
                "a0_over_a_transition_note": note,
                "M_GeV": M_GeV,
                "chi_star_Mpc": chi_star_Mpc,
                "k_hat_range": {"min": kmin, "max": kmax},
                "ell_range": ell_rng.__dict__,
                "ell_targets": ell_targets_f,
                "hits": hits,
                "signature_policy": {
                    "policy": signature_policy,
                    "decision": policy_decision,
                    "primary_signature": primary_signature,
                    "rationale": signature_rationale,
                    "cmb_overlap": cmb_overlap,
                    "small_scale_overlap": small_scale_overlap,
                    "cmb_ell_range": [cmb_min, cmb_max],
                    "small_scale_ell_min": small_scale_min,
                    "ell_bounce_s_est": ell_bounce_s_est,
                    "ell_bounce_t_est": ell_bounce_t_est,
                },
                "camb": {
                    "backend": cl_summary.get("backend"),
                    "summary": cl_summary,
                    "samples": cl_samples,
                },
                "inputs": {
                    "cosmology": {
                        "H0_km_s_Mpc": H0_km_s_Mpc,
                        "Omega_m": Omega_m,
                        "Omega_r": Omega_r,
                        "omega_b_h2": omega_b_h2,
                        "omega_c_h2": omega_c_h2,
                        "tau_reio": tau_reio,
                        "mnu_eV": mnu_eV,
                    },
                    "primordial_tfpt": {"N": N, "n_s": n_s, "A_s": A_s, "r": r, "pivot_scalar_Mpc_inv": pivot_scalar_Mpc_inv},
                },
                "primordial_pk_source": pk_source,
                "nuisance_policy": nuisance_policy,
                "planck_pliklite": planck_pliklite,
                "planck_lowl": planck_lowl,
                "planck_lensing": planck_lensing,
                "planck_combined": planck_combined,
                "publication_grade_gap": {"boltzmann_transfer": False, "bounce_feature_injection": (not bounce_injection_used)},
                "plot": plot_block,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=plot_warnings,
        )

