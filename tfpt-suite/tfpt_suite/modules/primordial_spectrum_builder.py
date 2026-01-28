from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.constants import TfptConstants
from tfpt_suite.cosmo_scale_map import CosmoScaleInputs, gev_to_mpc_inv, k_mpc_inv_from_k_hat, a0_over_a_transition
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass, mk_check_warn


def _interp_logx(*, x: float, x_grid: np.ndarray, y_grid: np.ndarray) -> float:
    """
    Interpolate y(x) on a log-spaced positive grid (clamped to endpoints).
    """
    if not np.isfinite(x) or x <= 0:
        return float("nan")
    xg = np.asarray(x_grid, dtype=float)
    yg = np.asarray(y_grid, dtype=float)
    if xg.size < 2:
        return float(yg[0]) if yg.size == 1 else float("nan")
    xmin = float(np.min(xg))
    xmax = float(np.max(xg))
    if x <= xmin:
        return float(yg[int(np.argmin(xg))])
    if x >= xmax:
        return float(yg[int(np.argmax(xg))])
    lx = float(np.log(x))
    lxg = np.log(xg)
    return float(np.interp(lx, lxg, yg))


class PrimordialSpectrumBuilderModule(TfptModule):
    module_id = "primordial_spectrum_builder"
    title = "Primordial spectrum builder (bounce T(k) injection → P_R(k), P_t(k) tables for CAMB)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "bounce transfer functions: bounce_perturbations output (k_grid, T_scalar, T_tensor, k_bounce_*_est_raw)",
                "k→Mpc^-1 mapping policy: tfpt_suite/data/k_calibration.json (a0/a_transition inputs)",
                "TFPT Starobinsky baseline: M/Mpl (constants) and N policy (from k_calibration.json assumptions)",
            ],
            outputs=[
                "k_grid_Mpc_inv (monotone log grid)",
                "P_R(k) and P_t(k) tables suitable for CAMBparams.set_initial_power_table",
                "explicit mapping metadata (a0/a_transition, k_bounce_s/t in k_hat=M units)",
            ],
            formulas=[
                r"P_\mathcal{R}(k) = P_{\mathcal{R},\mathrm{base}}(k)\,|T_s(k/k_{\mathrm{bounce},s})|^2",
                r"P_t(k) = P_{t,\mathrm{base}}(k)\,|T_t(k/k_{\mathrm{bounce},t})|^2",
                r"P_{\mathcal{R},\mathrm{base}}(k) = A_s (k/k_*)^{n_s-1},\;\; P_{t,\mathrm{base}}(k) = r\,A_s (k/k_*)^{n_t}",
            ],
            validation=[
                "bounce_feature_injection_wired is PASS when bounce outputs are present and P(k) tables are produced.",
                "tables are finite and monotone in k.",
            ],
            determinism="Deterministic given bounce outputs + declared cosmology mapping inputs.",
            question="Can we deterministically inject bounce transfer features into a primordial spectrum table consumable by a Boltzmann solver?",
            objective=[
                "Provide the missing bridge between bounce_perturbations and boltzmann_transfer (CAMB): a concrete P(k) table.",
            ],
            gaps=[
                "This uses a minimal multiplicative injection P_base*|T|^2. Publication-grade should document any windowing/smoothing policy and validate against numerical artifacts.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        # Load bounce outputs (required for this module’s core purpose).
        bounce_path = Path(config.output_dir) / "bounce_perturbations" / "results.json"
        if not bounce_path.is_file():
            return ModuleResult(
                results={"mode": mode, "bounce_results_file": str(bounce_path), "present": False},
                checks=[mk_check_warn("bounce_feature_injection_wired", f"missing bounce_perturburbations output at {bounce_path}")],
                report="\n".join(
                    [
                        "Primordial spectrum builder",
                        f"mode={mode}",
                        f"missing bounce output: {bounce_path}",
                    ]
                ),
                warnings=[],
            )

        bounce_payload = json.loads(bounce_path.read_text(encoding="utf-8"))
        bounce = bounce_payload.get("results", {}) if isinstance(bounce_payload, dict) else {}

        k_grid = np.asarray(bounce.get("k_grid", []), dtype=float)
        T_s = np.asarray(bounce.get("T_scalar", []), dtype=float)
        T_t = np.asarray(bounce.get("T_tensor", []), dtype=float)
        diag = bounce.get("diagnostics", {}) if isinstance(bounce.get("diagnostics", {}), dict) else {}
        k_bounce_s_raw = float(diag.get("k_bounce_s_est_raw", float("nan")))
        k_bounce_t_raw = float(diag.get("k_bounce_t_est_raw", float("nan")))

        # Load mapping policy (same source as boltzmann_transfer).
        data_dir = Path(__file__).resolve().parent.parent / "data"
        kcal_path = data_dir / "k_calibration.json"
        kcal = json.loads(kcal_path.read_text(encoding="utf-8"))
        asm = kcal.get("assumptions", {}) if isinstance(kcal.get("assumptions", {}), dict) else {}

        inp = CosmoScaleInputs(
            transition=str(asm.get("transition", "horizon_exit_of_pivot")),
            N_inflation_from_transition=float(asm.get("N_inflation_from_transition", 56.0)),
            N_reheat=float(asm.get("N_reheat", 0.0)),
            T_reheat_GeV=float(asm.get("T_reheat_GeV", 1.0e13)),
            g_star_s_reheat=float(asm.get("g_star_s_reheat", 120.0)),
            g_star_s_today=float(asm.get("g_star_s_today", 3.91)),
            T0_K=float(asm.get("T0_K", 2.7255)),
        )
        # Prefer the k_calibration expansion-budget estimate when available (brings k-scales into the declared
        # reheating/transition policy used by the suite).
        a0_over_atr = None
        note = "n/a"
        N_eff = float(inp.N_inflation_from_transition)
        kcal_out_path = Path(config.output_dir) / "k_calibration" / "results.json"
        if kcal_out_path.is_file():
            try:
                payload = json.loads(kcal_out_path.read_text(encoding="utf-8"))
                kres = payload.get("results", {}) if isinstance(payload, dict) else {}
                est = kres.get("expansion_budget_estimate", {}) if isinstance(kres.get("expansion_budget_estimate", {}), dict) else {}
                a0_est = float(est.get("a0_over_a_transition", float("nan")))
                N_est = float(est.get("N_inflation_from_transition", float("nan")))
                if math.isfinite(a0_est) and a0_est > 0:
                    a0_over_atr = a0_est
                    note = f"from k_calibration expansion_budget_estimate ({kcal_out_path})"
                if math.isfinite(N_est) and N_est > 0:
                    N_eff = N_est
            except Exception:
                pass
        if a0_over_atr is None:
            a0_over_atr, note = a0_over_a_transition(inp)
        if a0_over_atr is None or not math.isfinite(float(a0_over_atr)) or float(a0_over_atr) <= 0:
            return ModuleResult(
                results={"mode": mode, "bounce_results_file": str(bounce_path), "present": True, "a0_over_a_transition": a0_over_atr, "note": note},
                checks=[mk_check_warn("bounce_feature_injection_wired", f"a0/a_transition unavailable: {note}")],
                report="\n".join(
                    [
                        "Primordial spectrum builder",
                        f"mode={mode}",
                        f"bounce output: {bounce_path}",
                        f"a0/a_transition unavailable: {note}",
                    ]
                ),
                warnings=[],
            )

        # TFPT scale M and baseline primordial parameters.
        cst = TfptConstants.compute()
        M_GeV = float(float(cst.M_over_Mpl) * 2.435e18)
        N = int(round(float(N_eff)))
        n_s = float(1.0 - 2.0 / float(N))
        r = float(12.0 / (float(N) ** 2))
        n_t = float(-r / 8.0)
        A_s = float(((float(N) ** 2) / (24.0 * (math.pi**2))) * (float(cst.M_over_Mpl) ** 2))
        k_pivot = 0.05

        # Sanity: inputs present.
        ok_inputs = bool(
            k_grid.size == T_s.size == T_t.size
            and k_grid.size >= 3
            and math.isfinite(k_bounce_s_raw)
            and math.isfinite(k_bounce_t_raw)
            and k_bounce_s_raw > 0
            and k_bounce_t_raw > 0
        )

        # Choose a k grid (Mpc^-1) spanning both scalar and tensor bounce windows.
        # Map k_hat(M)=k/M using: k_hat(M)=k_hat(bounce)*k_bounce_raw.
        gev2mpc = float(gev_to_mpc_inv())
        kh_s_min = float(np.min(k_grid) * k_bounce_s_raw)
        kh_s_max = float(np.max(k_grid) * k_bounce_s_raw)
        kh_t_min = float(np.min(k_grid) * k_bounce_t_raw)
        kh_t_max = float(np.max(k_grid) * k_bounce_t_raw)
        k_s_min = float(k_mpc_inv_from_k_hat(k_hat=kh_s_min, M_GeV=M_GeV, a0_over_a_tr=float(a0_over_atr)))
        k_s_max = float(k_mpc_inv_from_k_hat(k_hat=kh_s_max, M_GeV=M_GeV, a0_over_a_tr=float(a0_over_atr)))
        k_t_min = float(k_mpc_inv_from_k_hat(k_hat=kh_t_min, M_GeV=M_GeV, a0_over_a_tr=float(a0_over_atr)))
        k_t_max = float(k_mpc_inv_from_k_hat(k_hat=kh_t_max, M_GeV=M_GeV, a0_over_a_tr=float(a0_over_atr)))
        k_min = float(min(k_s_min, k_t_min))
        k_max = float(max(k_s_max, k_t_max))

        # Ensure valid positive range for the table.
        k_min = max(k_min, 1e-6)
        k_max = max(k_max, k_min * 10.0)
        k_points = int(250)
        k_table = np.logspace(np.log10(k_min), np.log10(k_max), num=k_points)

        # Build transfer interpolation functions in k_hat(bounce) space.
        # For a given physical k, k_hat(M) = k * a0/a_tr / (M * GeV->Mpc^-1).
        def k_hat_M_from_k_mpc(k_mpc: float) -> float:
            return float((float(k_mpc) * float(a0_over_atr)) / (float(M_GeV) * float(gev2mpc)))

        pk_s = np.zeros_like(k_table)
        pk_t = np.zeros_like(k_table)
        for i, k_mpc in enumerate(k_table):
            khM = k_hat_M_from_k_mpc(float(k_mpc))
            khb_s = float(khM / k_bounce_s_raw)
            khb_t = float(khM / k_bounce_t_raw)
            Ts = _interp_logx(x=khb_s, x_grid=k_grid, y_grid=T_s)
            Tt = _interp_logx(x=khb_t, x_grid=k_grid, y_grid=T_t)
            # Baselines
            pr_base = float(A_s * (float(k_mpc) / k_pivot) ** (n_s - 1.0))
            pt_base = float((r * A_s) * (float(k_mpc) / k_pivot) ** (n_t))
            pk_s[i] = float(pr_base * (Ts * Ts))
            pk_t[i] = float(pt_base * (Tt * Tt))

        finite = bool(np.all(np.isfinite(k_table)) and np.all(np.isfinite(pk_s)) and np.all(np.isfinite(pk_t)))
        monotone = bool(np.all(np.diff(k_table) > 0))

        checks: list[Check] = []
        checks.append(mk_check_pass("bounce_feature_injection_wired", f"built P(k) tables from {bounce_path}") if ok_inputs and finite and monotone else mk_check_warn("bounce_feature_injection_wired", "failed to build finite monotone tables"))
        checks.append(mk_check_info("k_range", f"k range [{k_min:.3e}, {k_max:.3e}] Mpc^-1; points={k_points}"))

        report_lines: list[str] = [
            "Primordial spectrum builder (bounce injection)",
            f"mode={mode}",
            f"bounce input: {bounce_path}",
            f"k_calibration: {kcal_path}",
            "",
            "Mapping policy:",
            f"- a0/a_transition = {a0_over_atr} ({note})",
            f"- M ≈ {M_GeV:.3e} GeV",
            f"- k_bounce_s_raw (k_hat=M units) ≈ {k_bounce_s_raw}",
            f"- k_bounce_t_raw (k_hat=M units) ≈ {k_bounce_t_raw}",
            "",
            "Baseline primordial (TFPT/Starobinsky):",
            f"- N={N}, n_s={n_s:.6f}, A_s={A_s:.4e}, r={r:.6f}, n_t={n_t:.6g}, pivot={k_pivot} 1/Mpc",
            "",
            "Table:",
            f"- k_Mpc^-1 range [{k_min:.3e}, {k_max:.3e}] with {k_points} log points",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "bounce_results_file": str(bounce_path),
                "k_calibration_file": str(kcal_path),
                "mapping": {
                    "a0_over_a_transition": float(a0_over_atr),
                    "a0_over_a_transition_note": note,
                    "M_GeV": M_GeV,
                    "gev_to_mpc_inv": gev2mpc,
                },
                "bounce": {"k_bounce_s_est_raw": k_bounce_s_raw, "k_bounce_t_est_raw": k_bounce_t_raw, "k_grid_hat": [float(x) for x in k_grid.tolist()]},
                "baseline": {"N": N, "A_s": A_s, "n_s": n_s, "r": r, "n_t": n_t, "pivot_Mpc_inv": k_pivot},
                "table": {
                    "k_Mpc_inv": [float(x) for x in k_table.tolist()],
                    "P_R": [float(x) for x in pk_s.tolist()],
                    "P_t": [float(x) for x in pk_t.tolist()],
                },
            },
            checks=checks,
            report="\n".join(report_lines),
            warnings=[],
        )

