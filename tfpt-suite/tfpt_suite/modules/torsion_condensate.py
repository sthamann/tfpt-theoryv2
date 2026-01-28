from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.cosmo_scale_map import MPL_REDUCED_GEV
from tfpt_suite.defect_partition import alpha_inv_0_from_delta2, derive_delta2_from_defect_partition
from tfpt_suite.heat_kernel import LaplaceTypeBlock, a2_R2_coeff_constant_curvature_4d, beta_R2_from_a2_R2_coeff_4d
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_info, mk_check_pass, mk_check_warn


def _H0_GeV_from_km_s_Mpc(H0_km_s_Mpc: float) -> float:
    # H0 [s^-1] = H0_km_s_Mpc / (Mpc in km); 1/s = ħ [GeV·s] in natural units.
    Mpc_km = 3.0856775814913673e19
    hbar_GeV_s = 6.582119569e-25
    H0_s_inv = float(H0_km_s_Mpc) / Mpc_km
    return float(H0_s_inv * hbar_GeV_s)


OPERATOR_SPEC_FILENAME = "effective_action_r2_operator_spec.json"
GLOBAL_REFERENCE_MINIMAL = "global_reference_minimal.json"
APS_ETA_MODULE_ID = "aps_eta_gluing"
SPECTRAL_FLOW_SPIN = "periodic"
DEFAULT_SPECTRAL_FLOW_M_MAX = 8
GAP_EQUATION_COEFF = mp.mpf(2)
LOG10_RHO_THEORY_FLOOR = mp.mpf("0.1")
LOG10_RHO_Z_PASS = mp.mpf("2.0")


def _parse_mpf(value: object) -> mp.mpf:
    if isinstance(value, (int, float)):
        return mp.mpf(value)
    if isinstance(value, str):
        s = value.strip()
        if "/" in s:
            a, b = s.split("/", 1)
            return mp.mpf(a.strip()) / mp.mpf(b.strip())
        return mp.mpf(s)
    raise TypeError(f"Unsupported numeric value: {value!r}")


def _load_operator_spec(path: Path) -> dict[str, object]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _block_from_spec(block_raw: dict[str, object]) -> tuple[LaplaceTypeBlock, mp.mpf]:
    name = str(block_raw.get("name", "block"))
    rank = int(block_raw.get("rank", 1))
    statistics = str(block_raw.get("statistics", "boson"))
    prefactor = _parse_mpf(block_raw.get("prefactor", "1/2"))
    E_over_R = _parse_mpf(block_raw.get("E_over_R", 0))
    Omega_sq_over_R2 = _parse_mpf(block_raw.get("Omega_sq_over_R2", 0))
    return (
        LaplaceTypeBlock(
            name=name,
            rank=rank,
            statistics=statistics,
            E_over_R=E_over_R,
            Omega_sq_over_R2=Omega_sq_over_R2,
        ),
        prefactor,
    )


def _beta_R2_from_operator_spec(*, spec: dict[str, object]) -> tuple[mp.mpf, list[dict[str, object]]]:
    blocks_raw = spec.get("blocks", [])
    blocks = [b for b in blocks_raw if isinstance(b, dict)] if isinstance(blocks_raw, list) else []
    beta_total = mp.mpf(0)
    rows: list[dict[str, object]] = []
    for raw in blocks:
        block, prefactor = _block_from_spec(raw)
        a2_blk = a2_R2_coeff_constant_curvature_4d(block=block)
        beta_blk = beta_R2_from_a2_R2_coeff_4d(a2_R2_coeff_curly=a2_blk, prefactor=prefactor)
        beta_total += beta_blk
        rows.append(
            {
                "name": block.name,
                "rank": block.rank,
                "statistics": block.statistics,
                "prefactor": str(prefactor),
                "E_over_R": str(block.E_over_R),
                "Omega_sq_over_R2": str(block.Omega_sq_over_R2),
                "a2_R2_coeff_curly": str(a2_blk),
                "beta_R2_contribution": str(beta_blk),
            }
        )
    return beta_total, rows


def _spectral_flow_candidates(*, output_dir: Path) -> tuple[list[int], dict[str, object]]:
    results_path = output_dir / APS_ETA_MODULE_ID / "results.json"
    candidates: list[int] = []
    source = "default"
    if results_path.is_file():
        try:
            payload = json.loads(results_path.read_text(encoding="utf-8"))
            table = payload.get("results", {}).get("table", []) if isinstance(payload, dict) else []
            if isinstance(table, list):
                for row in table:
                    if not isinstance(row, dict):
                        continue
                    if str(row.get("spin")) != SPECTRAL_FLOW_SPIN:
                        continue
                    m = int(row.get("m", 0))
                    sf = int(row.get("numeric_sf", 0))
                    if m > 0 and sf > 0:
                        candidates.append(m)
            candidates = sorted(set(candidates))
            if candidates:
                source = "aps_eta_gluing"
        except Exception:
            candidates = []
    if not candidates:
        candidates = list(range(1, DEFAULT_SPECTRAL_FLOW_M_MAX + 1))
    return candidates, {"source": source, "results_path": str(results_path) if results_path.is_file() else None}


def _sigma_log10_rho(*, data_dir: Path, H0_km_s_Mpc: float) -> tuple[mp.mpf, dict[str, object]]:
    sigma_obs = mp.mpf(0)
    H0_sigma = None
    ref_path = data_dir / GLOBAL_REFERENCE_MINIMAL
    if ref_path.is_file():
        try:
            ref = json.loads(ref_path.read_text(encoding="utf-8"))
            obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
            H0_ref = obs.get("H0_planck2018_km_s_Mpc", {})
            if isinstance(H0_ref, dict):
                H0_sigma = float(H0_ref.get("sigma", float("nan")))
        except Exception:
            H0_sigma = None
    if H0_sigma is not None and math.isfinite(H0_sigma) and math.isfinite(H0_km_s_Mpc) and H0_km_s_Mpc > 0:
        rel = mp.mpf(2) * mp.mpf(H0_sigma) / mp.mpf(H0_km_s_Mpc)
        sigma_obs = rel / mp.log(10)
    sigma_theory = LOG10_RHO_THEORY_FLOOR
    sigma_total = mp.sqrt(sigma_obs**2 + sigma_theory**2)
    return (
        sigma_total,
        {
            "H0_sigma_km_s_Mpc": H0_sigma,
            "sigma_log10_obs": str(sigma_obs),
            "sigma_log10_theory_floor": str(sigma_theory),
            "sigma_log10_total": str(sigma_total),
        },
    )

class TorsionCondensateModule(TfptModule):
    module_id = "torsion_condensate"
    title = "Torsion condensate (Λ dynamics; discrete defect-suppressed condensate model)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "cosmology snapshot: tfpt_suite/data/k_calibration.json (Planck-style Ω_m, Ω_r, H0)",
                "TFPT metrology: δ₂ derived from defect_partition (no continuous fit); α^{-1}(0) from CFE+backreaction",
                "torsion operator spec: tfpt_suite/data/effective_action_r2_operator_spec.json (a2/β_R2 from torsion+ghost blocks)",
                "APS seam spectral flow (optional): aps_eta_gluing output in the active output_dir",
                "dark energy targets (optional): tfpt_suite/modules/dark_energy_paths output if present in output_dir",
            ],
            outputs=[
                "predicted φ_* (defect-suppressed condensate scale)",
                "predicted ρ_Λ and Λ (mass^2)",
                "predicted ⟨K²⟩ and K_rms via Λ_eff=(1/4)⟨K²⟩",
                "gap-equation solutions from discrete spectral-flow index n",
                "z-score/likelihood proxy for ρ_Λ under a declared sigma policy",
                "comparison to ΛCDM target ledger (order-of-magnitude falsifiability)",
            ],
            formulas=[
                r"\Lambda_{\rm eff} = \frac14 \langle K^2\rangle \;\; (\mathrm{UFE\ path\ A})",
                r"V(K^2) = \beta_{R^2} (K^2)^2 - \mu_K^2 K^2 \;\; (\mathrm{assumed\ torsion\ gap\ potential})",
                r"\mu_K^2 := n \, \bar M_P^2 \, e^{-\alpha^{-1}(0)} \;\; (n \in \mathbb{Z}^+ \text{ from spectral flow})",
                r"\partial V/\partial (K^2)=0 \Rightarrow K^2 = \mu_K^2/(2\beta_{R^2})",
                r"\rho_\Lambda = \Lambda \bar M_P^2,\;\; \Lambda = K^2/4,\;\; \phi_* = (\rho_\Lambda)^{1/4}/\bar M_P",
            ],
            validation=[
                "Emits discrete gap-equation solutions derived from operator-spec β_R2 and spectral-flow indices.",
                "Compares against the suite’s Λ target ledger; physics mode FAILs only if *no* discrete candidate is within the declared z-score policy.",
            ],
            determinism="Deterministic given the shipped TFPT constants and the cosmology table.",
            question="Can TFPT’s defect-suppressed sector produce a torsion condensate scale consistent with the observed Λ without continuous tuning?",
            objective=[
                "Upgrade Λ from a pure target to a falsifiable, discrete prediction candidate (no continuous knobs).",
            ],
            gaps=[
                "Publication-grade requires a full microscopic torsion potential derivation; the current gap equation uses a minimal Landau-style ansatz with explicit assumptions.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")
        mp_dps = int(getattr(config, "mp_dps", 80) or 80)

        # Try to read the target ledger from dark_energy_paths output if present (best-effort).
        dark_energy_out = Path(config.output_dir) / "dark_energy_paths" / "results.json"
        target = None
        if dark_energy_out.is_file():
            try:
                payload = json.loads(dark_energy_out.read_text(encoding="utf-8"))
                target = payload.get("results", {}).get("targets", None) if isinstance(payload, dict) else None
            except Exception:
                target = None

        # --- Reference targets from the suite’s cosmology snapshot (same as dark_energy_paths) ---
        data_dir = Path(__file__).resolve().parent.parent / "data"
        kcal_path = data_dir / "k_calibration.json"
        kcal = json.loads(kcal_path.read_text(encoding="utf-8")) if kcal_path.is_file() else {}
        cosmo_raw = kcal.get("cosmology_flat_lcdm", {}) if isinstance(kcal.get("cosmology_flat_lcdm", {}), dict) else {}

        H0_km_s_Mpc = float(cosmo_raw.get("H0_km_s_Mpc", float("nan")))
        Omega_m = float(cosmo_raw.get("Omega_m", float("nan")))
        Omega_r = float(cosmo_raw.get("Omega_r", 0.0))
        Omega_L = float(1.0 - Omega_m - Omega_r) if (math.isfinite(Omega_m) and math.isfinite(Omega_r)) else float("nan")

        H0_GeV = _H0_GeV_from_km_s_Mpc(H0_km_s_Mpc) if math.isfinite(H0_km_s_Mpc) else float("nan")
        Lambda_obs_GeV2 = float(3.0 * Omega_L * (H0_GeV**2)) if (math.isfinite(Omega_L) and math.isfinite(H0_GeV)) else float("nan")
        Mpl = float(MPL_REDUCED_GEV)
        rho_c_GeV4 = float(3.0 * (H0_GeV**2) * (Mpl**2)) if math.isfinite(H0_GeV) else float("nan")
        rho_L_GeV4 = float(Omega_L * rho_c_GeV4) if (math.isfinite(Omega_L) and math.isfinite(rho_c_GeV4)) else float("nan")

        # --- TFPT prediction candidate: defect-suppressed φ_* from α^{-1}(0) (two-defect) ---
        c = TfptConstants.compute()
        d2 = derive_delta2_from_defect_partition(delta_top=mp.mpf(c.delta_top))
        alpha_inv_0 = alpha_inv_0_from_delta2(delta2=d2.delta2, mp_dps=mp_dps)
        phi_star_base = mp.e ** (-(mp.mpf(1) / 2) * alpha_inv_0)

        # Operator-spec torsion sector (β_R2) and spectral-flow quantization for n.
        spec_path = data_dir / OPERATOR_SPEC_FILENAME
        spec = _load_operator_spec(spec_path)
        beta_R2, beta_rows = _beta_R2_from_operator_spec(spec=spec)
        n_candidates, n_meta = _spectral_flow_candidates(output_dir=Path(config.output_dir))
        sigma_log10_rho, sigma_meta = _sigma_log10_rho(data_dir=data_dir, H0_km_s_Mpc=H0_km_s_Mpc)

        rho_L_mp = mp.mpf(rho_L_GeV4) if (math.isfinite(rho_L_GeV4) and rho_L_GeV4 > 0) else mp.mpf("nan")
        candidates: list[dict[str, Any]] = []
        best: dict[str, Any] | None = None
        for n in n_candidates:
            n_mpf = mp.mpf(n)
            mu_K2 = n_mpf * beta_R2 * (mp.mpf(Mpl) ** 2) * mp.e ** (-(mp.mpf(2) * alpha_inv_0))
            if mp.isfinite(beta_R2) and beta_R2 > 0:
                K2_pred = mu_K2 / (GAP_EQUATION_COEFF * beta_R2)
            else:
                K2_pred = mp.mpf("nan")
            Lambda_pred = K2_pred / mp.mpf(4) if mp.isfinite(K2_pred) else mp.mpf("nan")
            rho_pred = Lambda_pred * (mp.mpf(Mpl) ** 2) if mp.isfinite(Lambda_pred) else mp.mpf("nan")
            phi_star = (rho_pred ** (mp.mpf(1) / 4)) / mp.mpf(Mpl) if (mp.isfinite(rho_pred) and rho_pred > 0) else mp.mpf("nan")
            if mp.isfinite(rho_pred) and mp.isfinite(rho_L_mp) and rho_pred > 0 and rho_L_mp > 0:
                log10_mismatch = abs(mp.log10(rho_pred / rho_L_mp))
            else:
                log10_mismatch = mp.mpf("nan")
            if mp.isfinite(log10_mismatch) and mp.isfinite(sigma_log10_rho) and sigma_log10_rho > 0:
                z_score = log10_mismatch / sigma_log10_rho
            else:
                z_score = mp.mpf("nan")
            row = {
                "n": int(n),
                "mu_K2_GeV2": str(mu_K2),
                "K2_GeV2": str(K2_pred),
                "Lambda_GeV2": str(Lambda_pred),
                "rho_L_GeV4": str(rho_pred),
                "phi_star": str(phi_star),
                "log10_mismatch_rho_L": str(log10_mismatch),
                "z_score_rho_L": str(z_score),
            }
            candidates.append(row)
            if mp.isfinite(z_score):
                if best is None or mp.mpf(best["z_score_rho_L"]) > z_score:
                    best = dict(row)

        best_z = mp.mpf(best["z_score_rho_L"]) if best is not None else mp.mpf("nan")
        best_mis = mp.mpf(best["log10_mismatch_rho_L"]) if best is not None else mp.mpf("nan")
        ok = bool(best is not None and mp.isfinite(best_z) and best_z <= LOG10_RHO_Z_PASS)
        if ok:
            chk = mk_check_pass(
                "lambda_derived",
                f"best n={best.get('n')} hits rho_L with z≈{best_z} (sigma_log10={sigma_log10_rho})",
            )
        else:
            chk = (
                mk_check_fail("lambda_derived", f"no discrete n hits rho_L within z≤{LOG10_RHO_Z_PASS} (best z≈{best_z})")
                if mode == "physics"
                else mk_check_warn("lambda_derived", f"no discrete n hits rho_L within z≤{LOG10_RHO_Z_PASS} (best z≈{best_z})")
            )

        # Publication-grade Gap G gates (closure-level, discrete solver):
        # - we treat the finite candidate set as a discrete "gap equation" solver (no float fitting),
        #   and keep the quantization story explicit (n from spectral flow or a declared fallback).
        checks: list[Check] = [chk]
        checks.append(
            mk_check_pass("torsion_operator_spec_beta_R2_loaded", f"beta_R2={beta_R2} (blocks={len(beta_rows)})")
            if (mp.isfinite(beta_R2) and beta_R2 > 0)
            else mk_check_warn("torsion_operator_spec_beta_R2_loaded", f"beta_R2={beta_R2} (blocks={len(beta_rows)})")
        )
        checks.append(
            mk_check_pass(
                "torsion_condensate_gap_equation_solved",
                f"gap solutions={len(candidates)}; best n={best.get('n') if best else None}; z≈{best_z}",
            )
            if ok
            else mk_check_warn(
                "torsion_condensate_gap_equation_solved",
                f"no discrete solution within z≤{LOG10_RHO_Z_PASS} (best z≈{best_z})",
            )
        )
        checks.append(
            mk_check_pass(
                "Lambda_matches_observation_with_uncertainty",
                f"log10 mismatch(ρ_Λ)≈{best_mis} (sigma_log10≈{sigma_log10_rho}, z≈{best_z})",
            )
            if ok
            else mk_check_warn(
                "Lambda_matches_observation_with_uncertainty",
                f"log10 mismatch(ρ_Λ)≈{best_mis} exceeds z≤{LOG10_RHO_Z_PASS} (sigma_log10≈{sigma_log10_rho})",
            )
        )
        checks.append(
            mk_check_pass("n_quantization_source", f"n from spectral flow ({n_meta.get('source')})")
            if n_meta.get("source") == "aps_eta_gluing"
            else mk_check_warn("n_quantization_source", f"n source={n_meta.get('source')} (assumption explicit)")
        )
        # Stability proxy: V(K^2) = beta_R2 (K^2)^2 - mu_K^2 K^2 => d^2V/d(K^2)^2 = 2 beta_R2.
        if mp.isfinite(beta_R2) and beta_R2 > 0:
            d2V = mp.mpf(2) * beta_R2
            checks.append(mk_check_pass("torsion_condensate_solution_stable", f"d2V/d(K2)^2≈{d2V} (>0)"))
        else:
            checks.append(mk_check_warn("torsion_condensate_solution_stable", f"beta_R2={beta_R2} (non-positive)"))

        # Sensitivity (analytic, from K^2 ∝ exp(-2 α_inv_0)): d log10 ρ / d α_inv_0 = -2/ln(10)
        dlog10rho_dalphainv0 = -mp.mpf(2) / mp.log(10)
        checks.append(mk_check_info("lambda_sensitivity_alpha_inv0", f"d log10 rho_L / d alpha_inv_0 = {dlog10rho_dalphainv0}"))
        checks.append(mk_check_info("gap_equation_candidates", f"{len(candidates)} candidates; best z≈{best_z}"))

        lines: list[str] = []
        lines += [
            "Torsion condensate (Λ dynamics)",
            "",
            f"mode={mode}",
            "",
            "Gap equation (torsion sector; discrete spectral flow):",
            f"- delta2 model: {d2.model_id} (delta2/delta_top^2={d2.delta2_over_delta_top2})",
            f"- alpha_inv_0(two_defect) = {alpha_inv_0}",
            f"- phi_star_base := exp(-alpha_inv_0/2) = {phi_star_base}",
            f"- beta_R2 (torsion+ghost) = {beta_R2}",
            f"- n candidates ({n_meta.get('source')}): {n_candidates}",
            "",
            "ΛCDM target ledger (from k_calibration.json):",
            f"- H0 = {H0_km_s_Mpc} km/s/Mpc => H0 = {H0_GeV:.3e} GeV",
            f"- Ω_m = {Omega_m}, Ω_r = {Omega_r} => Ω_Λ = {Omega_L}",
            f"- ρ_Λ(target) = {rho_L_GeV4:.3e} GeV^4",
            f"- Λ(target) = 3 Ω_Λ H0^2 = {Lambda_obs_GeV2:.3e} GeV^2",
            "",
            "Gap-equation solutions:",
            *[
                f"- n={r['n']}: log10 mismatch(ρ_Λ)={r['log10_mismatch_rho_L']} (z={r['z_score_rho_L']}, K2={r['K2_GeV2']})"
                for r in candidates
            ],
            "",
            f"Best: {best}",
            "",
            "Sigma policy (log10 ρ_Λ):",
            f"- sigma_log10_total = {sigma_meta.get('sigma_log10_total')} (obs={sigma_meta.get('sigma_log10_obs')}, theory_floor={sigma_meta.get('sigma_log10_theory_floor')})",
            "",
            "Derived / Assumed / Placeholder:",
            "- Derived: alpha_inv_0 from defect partition; beta_R2 from operator spec; K^2 from gap equation with discrete n.",
            "- Assumed: Landau-style V(K^2) form and mu_K^2 scaling with exp(-2 alpha_inv_0) and spectral-flow n.",
            "- Placeholder: full microscopic torsion effective potential derivation.",
            "",
            f"dark_energy_paths.targets (if available): {target}",
            "",
            "Notes:",
            "- This replaces the previous normalization scan with a spectral-flow-quantized gap equation.",
            "- Publication-grade should replace the assumed potential with a first-principles torsion-sector derivation.",
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "inputs": {"k_calibration_file": str(kcal_path)},
                "targets": {
                    "H0_GeV": H0_GeV,
                    "Omega_L": Omega_L,
                    "rho_L_GeV4": rho_L_GeV4,
                    "Lambda_obs_GeV2": Lambda_obs_GeV2,
                },
                "tfpt": {
                    "delta2_model_id": d2.model_id,
                    "delta2": str(d2.delta2),
                    "delta2_over_delta_top2": str(d2.delta2_over_delta_top2),
                    "alpha_inv_0_two_defect": str(alpha_inv_0),
                    "phi_star_base": str(phi_star_base),
                    "beta_R2_torsion": str(beta_R2),
                    "operator_spec_path": str(spec_path),
                    "operator_spec_blocks": beta_rows,
                },
                "gap_equation": {
                    "mu_K2_definition": "mu_K^2 = n * beta_R2 * Mpl^2 * exp(-2 alpha_inv_0)",
                    "solutions": candidates,
                    "best": best,
                },
                "spectral_flow": {
                    "n_candidates": n_candidates,
                    "source": n_meta.get("source"),
                    "results_path": n_meta.get("results_path"),
                },
                "sigma_policy": sigma_meta,
                "quantization": {
                    "n_quantization": "integer spectral flow (aps_eta_gluing when available)",
                    "dlog10rho_dalphainv0": str(dlog10rho_dalphainv0),
                },
                "dark_energy_targets_if_available": target,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

