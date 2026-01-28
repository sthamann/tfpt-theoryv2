from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from tfpt_suite.cosmo_scale_map import (
    As_from_ln10_As,
    CosmoScaleInputs,
    MPL_REDUCED_GEV,
    N_reheat_from_rho_ratio,
    a0_over_a_transition as a0_over_a_transition_entropy,
    deltaN_from_rho_ratio,
    rho_end_GeV4_from_As_r,
    rho_reheat_GeV4,
    starobinsky_N_from_ns,
    starobinsky_r_from_N,
)
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass, mk_check_warn


DEFAULT_REHEAT_THRESHOLD_KEYS = ("MSigma", "MG8", "MNR1", "MNR2", "MNR3")
DEFAULT_W_REH = 0.0
DEFAULT_GSTAR_REHEAT = 120.0
DEFAULT_GSTAR_S_REHEAT = 120.0
DEFAULT_TREH_POLICY = "min_threshold"
DEFAULT_RADIATION_W = 1.0 / 3.0


def _load_thresholds(thr: dict[str, Any]) -> dict[str, float]:
    raw = thr.get("thresholds_GeV", {}) if isinstance(thr.get("thresholds_GeV", {}), dict) else {}
    thresholds: dict[str, float] = {}
    for key in DEFAULT_REHEAT_THRESHOLD_KEYS:
        if key in raw:
            try:
                val = float(raw[key])
            except Exception:
                continue
            if math.isfinite(val) and val > 0:
                thresholds[key] = val
    return thresholds


def _sorted_regimes(thresholds: dict[str, float], w_reh: float, T_reh: float) -> list[dict[str, float | str]]:
    ordered = sorted(thresholds.items(), key=lambda kv: kv[1], reverse=True)
    regimes: list[dict[str, float | str]] = []
    for idx, (name, val) in enumerate(ordered):
        T_high = val
        T_low = ordered[idx + 1][1] if idx + 1 < len(ordered) else T_reh
        w_eff = w_reh if T_low >= T_reh else DEFAULT_RADIATION_W
        regimes.append(
            {
                "label": f"{name}_to_{ordered[idx + 1][0]}" if idx + 1 < len(ordered) else f"{name}_to_Treh",
                "T_high_GeV": float(T_high),
                "T_low_GeV": float(T_low),
                "w_eff": float(w_eff),
            }
        )
    if T_reh > 0:
        regimes.append(
            {
                "label": "post_reheat_radiation",
                "T_high_GeV": float(T_reh),
                "T_low_GeV": 0.0,
                "w_eff": float(DEFAULT_RADIATION_W),
            }
        )
    return regimes


class CosmoThresholdHistoryModule(TfptModule):
    module_id = "cosmo_threshold_history"
    title = "Cosmology threshold history → derived reheating temperature and expansion"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "thresholds: tfpt_suite/data/rge_thresholds_v25.json",
                "reheating policy: tfpt_suite/data/cosmo_reheating_policy_v106.json",
                "external refs: tfpt_suite/data/global_reference_minimal.json (n_s, ln10_As)",
            ],
            outputs=[
                "threshold-derived T_reh, N_reh, and DeltaN",
                "entropy-based a0/a_transition estimate using derived reheating inputs",
                "explicit threshold regime table with effective w",
            ],
            formulas=[
                "T_reh := min(MSigma, MG8, MNR1..3) under the default threshold policy",
                "rho_reh = (pi^2/30) g* T_reh^4",
                "DeltaN = (1-3w)/(12(1+w)) ln(rho_reh/rho_end)",
                "N_reh = (1/(3(1+w))) ln(rho_end/rho_reh)",
            ],
            validation=[
                "Derived T_reh is sourced from TFPT thresholds when available.",
                "Derived N_reh and a0/a_transition are finite under the declared policy.",
            ],
            determinism="Deterministic given thresholds + declared reheating policy inputs.",
            question="Can the reheating temperature and expansion be derived from the TFPT threshold ladder instead of being a free input?",
            objective=[
                "Make reheating a deterministic function of TFPT thresholds (MSigma/MG8/MNR).",
                "Provide derived N_reh and a0/a_transition for k-calibration without manual tuning.",
            ],
            assumptions=[
                "Default threshold policy uses the minimum of (MSigma, MG8, MNR1..3) as T_reh.",
                "Use a single effective reheating w_reh from the v1.06 policy for the pre-radiation stage.",
                "Entropy conservation is used to estimate a0/a_transition once T_reh is fixed.",
            ],
            gaps=[
                "A full first-principles reheating history would derive w(t) and g*(T) from microphysics; this module encodes a deterministic policy using available thresholds.",
            ],
            references=[
                "tfpt_suite/data/rge_thresholds_v25.json",
                "tfpt_suite/data/cosmo_reheating_policy_v106.json",
                "tfpt_suite/cosmo_scale_map.py",
            ],
            maturity="policy-derived, threshold-anchored (deterministic but still simplified)",
        )

    def run(self, config) -> ModuleResult:
        data_dir = Path(__file__).resolve().parent.parent / "data"
        thr_path = data_dir / "rge_thresholds_v25.json"
        pol_path = data_dir / "cosmo_reheating_policy_v106.json"
        ref_path = data_dir / "global_reference_minimal.json"

        thr_raw = json.loads(thr_path.read_text(encoding="utf-8")) if thr_path.is_file() else {}
        pol_raw = json.loads(pol_path.read_text(encoding="utf-8")) if pol_path.is_file() else {}
        ref_raw = json.loads(ref_path.read_text(encoding="utf-8")) if ref_path.is_file() else {}

        thresholds = _load_thresholds(thr_raw if isinstance(thr_raw, dict) else {})
        threshold_values = list(thresholds.values())

        ass = pol_raw.get("assumptions", {}) if isinstance(pol_raw.get("assumptions", {}), dict) else {}
        w_reh = float(ass.get("w_reh", DEFAULT_W_REH))
        g_star = float(ass.get("g_star_reheat", DEFAULT_GSTAR_REHEAT))
        g_star_s = float(ass.get("g_star_s_reheat", DEFAULT_GSTAR_S_REHEAT))
        T_can = float(ass.get("T_reheat_GeV_canonical", 0.0))

        derived_from_thresholds = bool(threshold_values)
        if derived_from_thresholds:
            T_reh = min(threshold_values)
            T_source = min(thresholds, key=thresholds.get)
        else:
            T_reh = T_can
            T_source = "policy_fallback"

        obs = ref_raw.get("observables", {}) if isinstance(ref_raw.get("observables", {}), dict) else {}
        ns = float(obs.get("n_s_planck2018", {}).get("mean", float("nan"))) if isinstance(obs.get("n_s_planck2018", {}), dict) else float("nan")
        ln10_As = (
            float(obs.get("ln10_As_planck2018", {}).get("mean", float("nan")))
            if isinstance(obs.get("ln10_As_planck2018", {}), dict)
            else float("nan")
        )
        As = As_from_ln10_As(ln10_As)
        N_pivot = starobinsky_N_from_ns(ns)
        r = starobinsky_r_from_N(N_pivot)
        rho_end = rho_end_GeV4_from_As_r(As=As, r=r, c_end=float(ass.get("c_end", 0.35)), Mpl_reduced_GeV=MPL_REDUCED_GEV)
        rho_reh = rho_reheat_GeV4(T_reheat_GeV=T_reh, g_star=g_star)
        deltaN = deltaN_from_rho_ratio(w_reh=w_reh, rho_reh=rho_reh, rho_end=rho_end)
        N_reh = N_reheat_from_rho_ratio(w_reh=w_reh, rho_reh=rho_reh, rho_end=rho_end)

        a0_over_at, a0_note = None, "T_reh not available"
        if math.isfinite(T_reh) and T_reh > 0 and math.isfinite(N_pivot):
            cosmo_inp = CosmoScaleInputs(
                transition="horizon_exit_of_pivot",
                N_inflation_from_transition=float(N_pivot),
                N_reheat=float(N_reh),
                T_reheat_GeV=float(T_reh),
                g_star_s_reheat=float(g_star_s),
            )
            a0_over_at, a0_note = a0_over_a_transition_entropy(cosmo_inp)

        regimes = _sorted_regimes(thresholds, w_reh=w_reh, T_reh=T_reh)

        mode = str(getattr(config, "verification_mode", "engineering"))
        checks: list[Check] = []
        checks.append(mk_check_pass("thresholds_loaded", f"thresholds={list(thresholds.keys())}") if thresholds else mk_check_warn("thresholds_loaded", "no threshold values found"))
        checks.append(
            mk_check_pass("T_reh_from_thresholds", f"T_reh={T_reh:.3e} GeV (source={T_source})")
            if derived_from_thresholds
            else mk_check_warn("T_reh_from_thresholds", f"T_reh={T_reh:.3e} GeV (fallback={T_source})")
        )
        checks.append(
            mk_check_pass("N_reh_finite", f"N_reh={N_reh:.4f}, DeltaN={deltaN:.4f}")
            if math.isfinite(N_reh) and math.isfinite(deltaN)
            else mk_check_warn("N_reh_finite", f"N_reh={N_reh}, DeltaN={deltaN}")
        )
        if a0_over_at is not None and math.isfinite(float(a0_over_at)) and float(a0_over_at) > 0:
            checks.append(mk_check_pass("a0_over_a_transition_estimated", f"{a0_note}; a0/a_t={float(a0_over_at):.3e}"))
        else:
            checks.append(mk_check_warn("a0_over_a_transition_estimated", f"{a0_note}"))
        checks.append(mk_check_info("threshold_policy", f"policy={DEFAULT_TREH_POLICY}, threshold_keys={list(thresholds.keys())}"))

        report_lines = [
            "Threshold-driven reheating history (policy-derived)",
            f"mode = {mode}",
            f"threshold file: {thr_path}",
            f"policy file: {pol_path}",
            f"reference file: {ref_path}",
            "",
            "Thresholds used:",
            *[f"- {k} = {v:.3e} GeV" for k, v in thresholds.items()],
            "",
            "Derived reheating:",
            f"- T_reh = {T_reh:.3e} GeV (source={T_source})",
            f"- w_reh = {w_reh}",
            f"- g* = {g_star}, g*_s = {g_star_s}",
            f"- N_reh = {N_reh:.4f}",
            f"- DeltaN = {deltaN:.4f}",
            f"- a0/a_transition = {a0_over_at if a0_over_at is not None else 'n/a'}",
            "",
            "Regime summary (effective w):",
            *[f"- {row['label']}: T_high={row['T_high_GeV']:.3e} GeV, T_low={row['T_low_GeV']:.3e} GeV, w={row['w_eff']}" for row in regimes],
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "thresholds": thresholds,
                "reheating": {
                    "T_reheat_GeV": T_reh,
                    "T_reheat_source": T_source,
                    "N_reheat": N_reh,
                    "deltaN": deltaN,
                    "w_reh": w_reh,
                    "g_star_reheat": g_star,
                    "g_star_s_reheat": g_star_s,
                    "a0_over_a_transition": a0_over_at,
                    "a0_note": a0_note,
                },
                "pivot": {"n_s": ns, "ln10_As": ln10_As, "N_pivot": N_pivot, "r": r},
                "regimes": regimes,
                "policy": {"threshold_keys": list(thresholds.keys()), "T_reh_policy": DEFAULT_TREH_POLICY},
            },
            checks=checks,
            report="\n".join(report_lines),
            warnings=[],
        )
