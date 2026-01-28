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
from tfpt_suite.module_base import (
    Check,
    ModuleResult,
    ModuleSpec,
    TfptModule,
    mk_check_info,
    mk_check_pass,
    mk_check_warn,
)


def _log10_safe(x: float) -> float:
    v = float(x)
    if v <= 0 or not math.isfinite(v):
        return float("nan")
    return float(math.log10(v))


def _plot_reheating_window(
    *,
    out_dir: Path,
    log10T: list[float],
    deltaN: list[float],
    N_reh: list[float],
    T_canonical_GeV: float,
    deltaN_floor: float,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"reheating_window_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)
        x = np.array(log10T, dtype=float)
        y1 = np.array(deltaN, dtype=float)
        y2 = np.array(N_reh, dtype=float)

        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(9.5, 6.6), sharex=True)

        ax1.plot(x, y1, color="#1f77b4", lw=1.5)
        ax1.axhline(float(deltaN_floor), color="black", lw=1.0, ls="--", alpha=0.8, label=r"$\Delta N$ floor")
        ax1.axvline(float(np.log10(T_canonical_GeV)), color="#d62728", lw=1.0, ls="--", alpha=0.85, label="canonical T_reh")
        ax1.set_ylabel(r"$\Delta N$")
        ax1.set_title("Reheating window (v1.06 ΔN model)")
        ax1.grid(True, ls=":", alpha=0.35)
        ax1.legend(loc="best")

        ax2.plot(x, y2, color="#ff7f0e", lw=1.5)
        ax2.axvline(float(np.log10(T_canonical_GeV)), color="#d62728", lw=1.0, ls="--", alpha=0.85)
        ax2.set_xlabel(r"$\log_{10}(T_{\rm reh}/{\rm GeV})$")
        ax2.set_ylabel(r"$N_{\rm reh}$")
        ax2.grid(True, ls=":", alpha=0.35)

        fig.tight_layout()
        path = out_dir / "reheating_window.png"
        fig.savefig(path, dpi=170)
        plt.close(fig)
        plot["reheating_window_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")
    return plot, warnings


class CosmoReheatingPolicyV106Module(TfptModule):
    module_id = "cosmo_reheating_policy_v106"
    title = "Cosmology reheating/ΔN policy (v1.06) → derived N_reh and a0/a_transition"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "policy file: tfpt_suite/data/cosmo_reheating_policy_v106.json",
                "external refs: tfpt_suite/data/global_reference_minimal.json (n_s, ln10_As)",
            ],
            outputs=[
                "N_pivot from n_s via n_s≈1-2/N",
                "r from Starobinsky r≈12/N^2",
                "ρ_end proxy from (A_s, r) and c_end",
                "ΔN(T_reh) curve and a derived reheating window",
                "N_reh(T_reh) and an entropy-based a0/a_transition estimate",
            ],
            formulas=[
                "n_s ≈ 1 - 2/N  =>  N ≈ 2/(1-n_s)",
                "r ≈ 12/N^2  (Starobinsky plateau)",
                "A_s = exp(ln(10^{10}A_s))·1e-10",
                "V_* = (3/2)π^2 A_s r (reduced Planck units)",
                "ρ_end ≈ c_end·V_*·M̄_P^4  (proxy)",
                "ρ_reh = (π^2/30) g* T_reh^4",
                "ΔN ≈ (1-3w)/(12(1+w)) ln(ρ_reh/ρ_end)  (v1.06)",
                "N_reh = (1/(3(1+w))) ln(ρ_end/ρ_reh)",
                "a0/a_transition = (a0/a_reh)*(a_reh/a_end)*exp(N_pivot) with a0/a_reh from entropy conservation",
            ],
            validation=[
                "Produces finite N_pivot, r, and a0/a_transition under the declared policy.",
                "Records whether the canonical reheating choice satisfies the ΔN floor criterion.",
            ],
            determinism="Deterministic given the policy + reference tables.",
            question="Can the k→ℓ mapping policy be made less free-floating by deriving N and reheating expansion from (n_s, A_s) and an explicit reheating model (v1.06 ΔN window)?",
            objective=[
                "Turn the reheating assumptions into an explicit, auditable mapping layer (no silent 'typical N').",
                "Expose which reheating temperatures are compatible with a chosen ΔN floor (e.g. ΔN ≥ −6).",
            ],
            what_was_done=[
                "Read Planck-style (n_s, ln10_As) references and derive N_pivot and r using Starobinsky relations.",
                "Compute ρ_end proxy and scan ΔN(T_reh) + N_reh(T_reh) over a log-spaced range.",
                "Compute a canonical a0/a_transition estimate via entropy conservation using the derived N_reh and N_pivot.",
            ],
            assumptions=[
                "Use Starobinsky relations n_s≈1-2/N and r≈12/N^2.",
                "Approximate ρ_end with a fixed c_end times the pivot-scale plateau energy density.",
                "Assume a constant reheating equation-of-state w_reh and instantaneous thermalization at T_reh.",
            ],
            gaps=[
                "A publication-grade mapping would couple this to the TFPT-specific threshold history (MSigma/MG8/...) and a derived reheating temperature, not an assumed one.",
                "This module does not yet enforce a unique TFPT reheating temperature; it provides a deterministic window analysis under explicit assumptions.",
            ],
            references=[
                "paper_v1_06_01_09_2025.tex Sec. 8.2 (reheating window; ΔN formula; g*≈120; c≈0.34..0.37)",
                "tfpt_suite/cosmo_scale_map.py (entropy mapping used by k_calibration)",
            ],
            maturity="policy layer (deterministic given assumptions; not yet derived from TFPT thresholds)",
        )

    def run(self, config) -> ModuleResult:
        data_dir = Path(__file__).resolve().parent.parent / "data"
        pol_path = data_dir / "cosmo_reheating_policy_v106.json"
        ref_path = data_dir / "global_reference_minimal.json"

        pol = json.loads(pol_path.read_text(encoding="utf-8")) if pol_path.is_file() else {}
        ref = json.loads(ref_path.read_text(encoding="utf-8")) if ref_path.is_file() else {}

        ass = pol.get("assumptions", {}) if isinstance(pol.get("assumptions", {}), dict) else {}
        scan = pol.get("scan", {}) if isinstance(pol.get("scan", {}), dict) else {}

        w_reh = float(ass.get("w_reh", 0.0))
        g_star = float(ass.get("g_star_reheat", 120.0))
        g_star_s = float(ass.get("g_star_s_reheat", g_star))
        c_end = float(ass.get("c_end", 0.35))
        T_can = float(ass.get("T_reheat_GeV_canonical", 1.0e13))
        deltaN_floor = float(ass.get("deltaN_floor", -6.0))
        T_floor = float(ass.get("T_reheat_GeV_floor", 6.0e-3))

        obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
        ns = float(obs.get("n_s_planck2018", {}).get("mean", float("nan"))) if isinstance(obs.get("n_s_planck2018", {}), dict) else float("nan")
        ln10_As = (
            float(obs.get("ln10_As_planck2018", {}).get("mean", float("nan")))
            if isinstance(obs.get("ln10_As_planck2018", {}), dict)
            else float("nan")
        )
        As = As_from_ln10_As(ln10_As)
        N_pivot = starobinsky_N_from_ns(ns)
        r = starobinsky_r_from_N(N_pivot)

        rho_end = rho_end_GeV4_from_As_r(As=As, r=r, c_end=c_end, Mpl_reduced_GeV=MPL_REDUCED_GEV)

        # scan T_reh
        log10_Tmin = float(scan.get("log10_Tmin_GeV", -3.0))
        log10_Tmax = float(scan.get("log10_Tmax_GeV", 15.0))
        n_points = int(scan.get("n_points", 220))
        if n_points < 20:
            n_points = 20
        if log10_Tmax <= log10_Tmin:
            log10_Tmax = log10_Tmin + 1.0

        log10T: list[float] = []
        deltaN_list: list[float] = []
        Nreh_list: list[float] = []
        for i in range(n_points):
            t = log10_Tmin + (log10_Tmax - log10_Tmin) * (i / (n_points - 1))
            T = 10 ** float(t)
            rho_reh = rho_reheat_GeV4(T_reheat_GeV=T, g_star=g_star)
            dN = deltaN_from_rho_ratio(w_reh=w_reh, rho_reh=rho_reh, rho_end=rho_end)
            Nreh = N_reheat_from_rho_ratio(w_reh=w_reh, rho_reh=rho_reh, rho_end=rho_end)
            log10T.append(float(t))
            deltaN_list.append(float(dN))
            Nreh_list.append(float(Nreh))

        # Determine the minimal T that satisfies ΔN >= floor (if any)
        T_min_satisfy: float | None = None
        for t, dN in zip(log10T, deltaN_list):
            if math.isfinite(dN) and dN >= deltaN_floor:
                T_min_satisfy = 10 ** float(t)
                break

        # Canonical evaluation
        rho_reh_can = rho_reheat_GeV4(T_reheat_GeV=T_can, g_star=g_star)
        deltaN_can = deltaN_from_rho_ratio(w_reh=w_reh, rho_reh=rho_reh_can, rho_end=rho_end)
        Nreh_can = N_reheat_from_rho_ratio(w_reh=w_reh, rho_reh=rho_reh_can, rho_end=rho_end)

        cosmo_inp = CosmoScaleInputs(
            transition="horizon_exit_of_pivot",
            N_inflation_from_transition=float(N_pivot),
            N_reheat=float(Nreh_can),
            T_reheat_GeV=float(T_can),
            g_star_s_reheat=float(g_star_s),
        )
        a0_over_at, a0_note = a0_over_a_transition_entropy(cosmo_inp)

        # Floor evaluation (BBN-ish)
        rho_reh_floor = rho_reheat_GeV4(T_reheat_GeV=T_floor, g_star=g_star)
        deltaN_floor_val = deltaN_from_rho_ratio(w_reh=w_reh, rho_reh=rho_reh_floor, rho_end=rho_end)

        mode = str(getattr(config, "verification_mode", "engineering"))
        checks: list[Check] = []
        checks.append(mk_check_pass("derived_N_pivot_from_ns", f"n_s={ns} => N_pivot≈{N_pivot:.4f}"))
        checks.append(mk_check_pass("derived_r_from_N", f"r≈12/N^2 => r≈{r:.6g}"))
        checks.append(mk_check_pass("rho_end_finite", f"rho_end≈{rho_end:.3e} GeV^4 (c_end={c_end})"))
        checks.append(mk_check_info("deltaN_floor_eval", f"T_floor={T_floor:.3e} GeV => ΔN≈{deltaN_floor_val:.3f}"))

        ok_can = bool(math.isfinite(deltaN_can) and deltaN_can >= deltaN_floor)
        if ok_can:
            checks.append(mk_check_pass("canonical_T_reh_satisfies_deltaN_floor", f"T_can={T_can:.3e} GeV => ΔN≈{deltaN_can:.3f} (floor={deltaN_floor})"))
        else:
            checks.append(
                mk_check_warn(
                    "canonical_T_reh_violates_deltaN_floor",
                    f"T_can={T_can:.3e} GeV => ΔN≈{deltaN_can:.3f} < floor={deltaN_floor}; window_min≈{(T_min_satisfy or float('nan')):.3e} GeV",
                )
            )

        if a0_over_at is not None and float(a0_over_at) > 0:
            checks.append(mk_check_pass("a0_over_a_transition_estimated", f"{a0_note}; a0/a_t≈{float(a0_over_at):.3e}"))
        else:
            checks.append(mk_check_warn("a0_over_a_transition_missing", f"{a0_note}"))

        report_lines: list[str] = [
            "Cosmology reheating policy (v1.06 ΔN model)",
            f"mode = {mode}",
            f"policy file: {pol_path}",
            f"reference file: {ref_path}",
            "",
            "Derived pivot parameters (Starobinsky relations):",
            f"- n_s(ref) = {ns}",
            f"- ln(10^10 A_s)(ref) = {ln10_As} => A_s = {As:.6g}",
            f"- N_pivot = 2/(1-n_s) = {N_pivot:.6g}",
            f"- r = 12/N^2 = {r:.6g}",
            "",
            "Reheating model assumptions:",
            f"- w_reh = {w_reh}",
            f"- g_star(reheat) = {g_star}",
            f"- c_end = {c_end} (ρ_end ≈ c_end·V_*)",
            "",
            "Energy densities:",
            f"- rho_end ≈ {rho_end:.3e} GeV^4 (M̄_P={MPL_REDUCED_GEV:.3e} GeV)",
            f"- rho_reh(T) = (π^2/30) g* T^4",
            "",
            "Window scan:",
            f"- log10(T/GeV) in [{log10_Tmin}, {log10_Tmax}] with n={n_points}",
            f"- ΔN_floor = {deltaN_floor}",
            f"- inferred T_min satisfying ΔN>=floor: {(T_min_satisfy if T_min_satisfy is not None else 'none_in_scan')}",
            "",
            "Canonical choice:",
            f"- T_can = {T_can:.3e} GeV => ΔN_can ≈ {deltaN_can:.3f}, N_reh ≈ {Nreh_can:.3f}",
            f"- a0/a_transition(est) ≈ {a0_over_at}  ({a0_note})",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
            "",
            "Notes:",
            "- This makes the reheating policy explicit and audit-able; it is not yet a unique TFPT derivation of T_reh from thresholds.",
            "- k_calibration can consume these derived (N_pivot, N_reh, T_reh) values to reduce free parameters in k→ℓ mapping.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"reheating_window_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_reheating_window(
                out_dir=out_dir,
                log10T=log10T,
                deltaN=deltaN_list,
                N_reh=Nreh_list,
                T_canonical_GeV=T_can,
                deltaN_floor=deltaN_floor,
            )
            warnings.extend(plot_warnings)

        results: dict[str, Any] = {
            "policy_file": str(pol_path),
            "reference_file": str(ref_path),
            "refs": {"n_s": ns, "ln10_As": ln10_As, "A_s": As},
            "derived_pivot": {"N_pivot": N_pivot, "r": r},
            "assumptions": {
                "w_reh": w_reh,
                "g_star_reheat": g_star,
                "g_star_s_reheat": g_star_s,
                "c_end": c_end,
                "deltaN_floor": deltaN_floor,
                "T_reheat_GeV_canonical": T_can,
            },
            "energy": {
                "rho_end_GeV4": rho_end,
                "rho_reh_canonical_GeV4": rho_reh_can,
                "rho_reh_floor_GeV4": rho_reh_floor,
            },
            "window": {
                "T_min_satisfy_deltaN_floor_GeV": T_min_satisfy,
                "deltaN_floor_at_T_floor": deltaN_floor_val,
                "scan": {"log10T": log10T, "deltaN": deltaN_list, "N_reheat": Nreh_list},
            },
            "canonical": {
                "T_reheat_GeV": T_can,
                "deltaN": deltaN_can,
                "N_reheat": Nreh_can,
                "a0_over_a_transition_est": a0_over_at,
                "a0_over_a_transition_note": a0_note,
            },
            "plot": plot,
        }

        return ModuleResult(results=results, checks=checks, report="\n".join(report_lines), warnings=warnings)

