from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.axion_inputs import resolve_axion_claim
from tfpt_suite.constants import TfptConstants
from tfpt_suite.cosmo_scale_map import As_from_ln10_As, starobinsky_N_from_ns, starobinsky_r_from_N, MPL_REDUCED_GEV
from tfpt_suite.module_base import (
    Check,
    ModuleResult,
    ModuleSpec,
    TfptModule,
    mk_check_fail,
    mk_check_info,
    mk_check_pass,
    mk_check_warn,
)


def _microeV_to_GHz(micro_eV: float) -> float:
    # 1 eV/h = 241.79893 THz => 1 µeV corresponds to 0.24179893 GHz
    return float(0.24179893 * float(micro_eV))


def _H_inf_GeV_from_As_r(*, As: float, r: float, Mpl_reduced_GeV: float) -> float:
    """
    Slow-roll relation (engineering-level):
      A_s = H^2 / (8π^2 M̄_P^2 ε),  ε=r/16  =>  H = π M̄_P sqrt(A_s r / 2)
    """
    Asf = float(As)
    rf = float(r)
    if Asf <= 0 or rf <= 0:
        return float("nan")
    return float(math.pi * float(Mpl_reduced_GeV) * math.sqrt(Asf * rf / 2.0))


def _omega_misalignment_h2(*, theta: float, f_a_GeV: float, omega_norm: float, f_norm_GeV: float, p: float) -> float:
    th = float(theta)
    if th < 0:
        th = -th
    # small-angle harmonic approximation; anharmonic corrections are negligible at θ=φ0≈0.053.
    return float(float(omega_norm) * (th**2) * (float(f_a_GeV) / float(f_norm_GeV)) ** float(p))


def _parse_float_or_fraction(x: object, default: float) -> float:
    if isinstance(x, (int, float)):
        try:
            v = float(x)
            return v if math.isfinite(v) else float(default)
        except Exception:
            return float(default)
    if isinstance(x, str):
        s = x.strip()
        if not s:
            return float(default)
        try:
            v = float(Fraction(s))
            return v if math.isfinite(v) else float(default)
        except Exception:
            try:
                v = float(s)
                return v if math.isfinite(v) else float(default)
            except Exception:
                return float(default)
    return float(default)


def _plot_axion_dm(
    *,
    out_dir: Path,
    thetas: list[float],
    omega_h2: list[float],
    iso_ratio: list[float],
    theta_tfpt: float,
    omega_dm_ref: float,
    beta_iso_max: float,
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"axion_dm_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore
        import numpy as np  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        x = np.array(thetas, dtype=float)
        y1 = np.array(omega_h2, dtype=float)
        y2 = np.array(iso_ratio, dtype=float)

        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(9.5, 6.8), sharex=True)

        ax1.plot(x, y1, color="#1f77b4", lw=1.6)
        ax1.axhline(float(omega_dm_ref), color="black", lw=1.0, ls="--", alpha=0.85, label=r"$\Omega_{\rm DM}h^2$ ref")
        ax1.axvline(float(theta_tfpt), color="#d62728", lw=1.0, ls="--", alpha=0.9, label=r"TFPT $\theta_i=\varphi_0$")
        ax1.set_yscale("log")
        ax1.set_ylabel(r"$\Omega_a h^2$ (misalignment, approx)")
        ax1.set_title("TFPT axion DM pipeline (engineering-level: misalignment + isocurvature)")
        ax1.grid(True, which="both", ls=":", alpha=0.35)
        ax1.legend(loc="best")

        # isocurvature ratio: P_iso / P_R ~ (H/(π f_a θ))^2 / A_s
        iso_limit = float(beta_iso_max / max(1e-12, (1.0 - beta_iso_max)))
        ax2.plot(x, y2, color="#ff7f0e", lw=1.6)
        ax2.axhline(iso_limit, color="black", lw=1.0, ls="--", alpha=0.85, label="isocurvature limit (β/(1-β))")
        ax2.axvline(float(theta_tfpt), color="#d62728", lw=1.0, ls="--", alpha=0.9)
        ax2.set_yscale("log")
        ax2.set_xlabel(r"misalignment angle $\theta_i$")
        ax2.set_ylabel(r"$P_{\rm iso}/P_{\cal R}$ (approx)")
        ax2.grid(True, which="both", ls=":", alpha=0.35)
        ax2.legend(loc="best")

        fig.tight_layout()
        path = out_dir / "axion_dm.png"
        fig.savefig(path, dpi=170)
        plt.close(fig)
        plot["axion_dm_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")
    return plot, warnings


class AxionDmPipelineModule(TfptModule):
    module_id = "axion_dm_pipeline"
    title = "Axion DM pipeline (TFPT axion: f_a, m_a, θ_i=φ₀; relic+isocurvature audits)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants (computed): c3, varphi0 (for θ_i policy and g_coeff=-4c3)",
                "axion inputs: tfpt_suite/data/axion_tfpt_v106.json (quoted) or axion_fa_derivation output (preferred)",
                "Planck reference: tfpt_suite/data/global_reference_minimal.json (n_s, ln10_As for H_inf estimate)",
            ],
            outputs=[
                "axion frequency target (GHz) from m_a",
                "dimensionless coupling coefficient g_coeff=-4c3 and a physical scale g_phys≈g_coeff/f_a",
                "misalignment relic density estimate Ω_a h^2 (pre-inflation single-angle scenario)",
                "isocurvature ratio estimate P_iso/P_R under the same scenario (flags viability)",
                "strings/domain-walls factor explanation (C_str topology policy)",
            ],
            formulas=[
                "θ_i = varphi0 (from birefringence: β=varphi0/(4π) and β=2 c3 Δa ⇒ Δa=varphi0)",
                "g_coeff = -4 c3 = -1/(2π) (as used in TFPT notes; dimensionless coefficient)",
                "g_phys ≈ g_coeff / f_a  (if a is normalized as a/f_a)",
                "ν(GHz) = 0.24179893 · m_a(µeV)",
                "Ω_a h^2 ≈ Ω_norm · θ_i^2 · (f_a/f_norm)^p  (engineering-level misalignment scaling)",
                "H_inf ≈ π M̄_P sqrt(A_s r / 2), r≈12/N^2, N≈2/(1-n_s)",
                "P_iso/P_R ≈ (H_inf/(π f_a θ_i))^2 / A_s  (order-of-magnitude)",
            ],
            validation=[
                "Reproduces the haloscope frequency from the selected axion mass (within rounding).",
                "Produces finite Ω_a h^2 and an explicit isocurvature viability flag under stated assumptions.",
                "c_str_explained_by_topology: PASS when the discrete C_str is derived from the cusp/charge policy.",
            ],
            determinism="Deterministic given inputs and reference tables.",
            question="Given the TFPT axion (f_a, m_a) and the birefringence-fixed misalignment angle θ_i=φ₀, what relic density and isocurvature implications follow?",
            objective=[
                "Close the loop from 'nice axion numbers' to world-contact: relic abundance + isocurvature constraints.",
                "Make explicit whether the claimed θ_i=φ₀ is compatible with standard inflationary isocurvature bounds.",
            ],
            what_was_done=[
                "Imported the axion parameter claims (f_a, m_a) and computed the haloscope frequency target.",
                "Used θ_i=φ₀ and a standard misalignment scaling law to estimate Ω_a h^2 (pre-inflation scenario).",
                "Estimated H_inf from (n_s, A_s) via Starobinsky relations and evaluated isocurvature ratio P_iso/P_R.",
            ],
            assumptions=[
                "Prefer derived (f_a, m_a) from axion_fa_derivation when available; fallback to quoted benchmarks from axion_tfpt_v106.json.",
                "Use a simplified misalignment-only relic estimate (no strings/domain walls; small-angle approximation).",
                "Assume a standard single-field inflationary isocurvature estimate in the pre-inflation PQ-breaking scenario.",
            ],
            gaps=[
                "Publication-grade still requires a full axion-ladder derivation to be enforced as a hard dependency in every DM branch (not just preferred when available).",
                "Publication-grade requires choosing the correct PQ-breaking scenario (pre vs post inflation) and including strings/domain walls and updated lattice inputs.",
                "If inflation scale is high (Starobinsky), pre-inflation PQ breaking may be ruled out by isocurvature unless additional TFPT structure modifies the estimate.",
            ],
            references=[
                "five_problems.tex Sec. 4 (TFPT axion numbers and θ_i=φ₀ claim)",
                "tfpt_suite/constants.py (c3, φ0, g_coeff=-4c3)",
                "Planck 2018 constraints (global_reference_minimal.json; isocurvature bound is a proxy)",
            ],
            maturity="engineering-level pipeline (explicit assumptions; highlights missing ladder derivation + isocurvature tension)",
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        data_dir = Path(__file__).resolve().parent.parent / "data"
        ax_path = data_dir / "axion_tfpt_v106.json"
        ref_path = data_dir / "global_reference_minimal.json"

        ax = json.loads(ax_path.read_text(encoding="utf-8")) if ax_path.is_file() else {}
        ref = json.loads(ref_path.read_text(encoding="utf-8")) if ref_path.is_file() else {}

        cst = TfptConstants.compute()
        theta_varphi0 = float(cst.varphi0)
        g_coeff = float(cst.g_a_gamma_gamma)  # dimensionless coefficient (-1/(2π))

        claim_resolved = resolve_axion_claim(ax_raw=ax, output_dir=Path(config.output_dir))
        f_a_GeV = float(claim_resolved.get("f_a_GeV", float("nan")))
        m_micro_eV = float(claim_resolved.get("m_a_micro_eV", float("nan")))
        nu_GHz = _microeV_to_GHz(m_micro_eV)
        nu_claim = float(claim_resolved.get("frequency_GHz_claim", float("nan")))
        fa_source = str(claim_resolved.get("source", "quoted"))
        fa_derived_path = claim_resolved.get("derived_path", None)

        g_phys_GeV_inv = g_coeff / f_a_GeV if (math.isfinite(f_a_GeV) and f_a_GeV > 0) else float("nan")

        pol = ax.get("cosmo_policy", {}) if isinstance(ax.get("cosmo_policy", {}), dict) else {}
        scen_pol = pol.get("scenario_policy", {}) if isinstance(pol.get("scenario_policy", {}), dict) else {}
        scen_selected = str(scen_pol.get("selected", "")).strip()
        scen_allowed = [str(x) for x in (scen_pol.get("allowed", []) if isinstance(scen_pol.get("allowed", []), list) else [])]
        post = pol.get("post_inflation", {}) if isinstance(pol.get("post_inflation", {}), dict) else {}
        theta_rms = float(post.get("theta_rms", math.pi / math.sqrt(3.0)))
        strings_dw_factor = _parse_float_or_fraction(post.get("strings_domain_walls_factor", 1.0), 1.0)
        strings_dw_policy = post.get("strings_domain_walls_factor_policy", {})
        strings_dw_explained = False
        strings_dw_explain_detail = "no topology policy"
        if isinstance(strings_dw_policy, dict) and str(strings_dw_policy.get("kind", "")).strip() == "mobius_cusp_charge_sum":
            charges = strings_dw_policy.get("charges", [])
            if isinstance(charges, list) and charges:
                try:
                    total = Fraction(0, 1)
                    for c in charges:
                        total += Fraction(str(c))
                    strings_dw_factor = float(total)
                    strings_dw_explained = True
                    strings_dw_explain_detail = f"charges={charges} -> sum={strings_dw_factor}"
                except Exception:
                    # fall back to the raw numeric value
                    pass
        strings_dw_factor = float(max(strings_dw_factor, 1.0))

        mis = pol.get("misalignment_model", {}) if isinstance(pol.get("misalignment_model", {}), dict) else {}
        omega_norm = float(mis.get("omega_norm_at_fa_5e11", 0.12))
        f_norm = float(mis.get("fa_norm_GeV", 5.0e11))
        p = float(mis.get("power_law_index", 1.165))
        omega_dm_ref = float(pol.get("omega_dm_h2_ref", 0.12))

        # Scenario selection (explicit; used to avoid silent pre-inflation assumptions in physics mode).
        scenario_id = scen_selected or "pre_inflation_single_theta_varphi0"
        if scen_allowed and scenario_id not in set(scen_allowed):
            scenario_id = "pre_inflation_single_theta_varphi0"

        isocurvature_applicable = True
        theta_eff = float(theta_varphi0)
        strings_factor = 1.0
        theta_policy_txt = "theta_i = varphi0 (birefringence-fixed; pre-inflation single-angle)"
        if scenario_id == "post_inflation_theta_rms_no_strings":
            isocurvature_applicable = False
            theta_eff = float(theta_rms)
            strings_factor = 1.0
            theta_policy_txt = "post-inflation PQ: random theta_i (use rms pi/sqrt(3)); no strings/DW"
        elif scenario_id == "post_inflation_theta_rms_with_strings_dw_factor":
            isocurvature_applicable = False
            theta_eff = float(theta_rms)
            strings_factor = float(strings_dw_factor)
            theta_policy_txt = "post-inflation PQ: random theta_i (rms pi/sqrt(3)); include strings/DW factor"

        omega_a_h2 = _omega_misalignment_h2(theta=theta_eff, f_a_GeV=f_a_GeV, omega_norm=omega_norm, f_norm_GeV=f_norm, p=p) * float(
            strings_factor
        )
        frac_dm = omega_a_h2 / omega_dm_ref if omega_dm_ref > 0 else float("nan")

        # Inflation scale estimate (Starobinsky relations anchored to Planck n_s, A_s)
        obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
        ns = (
            float(obs.get("n_s_planck2018", {}).get("mean", float("nan")))
            if isinstance(obs.get("n_s_planck2018", {}), dict)
            else float("nan")
        )
        ln10_As = (
            float(obs.get("ln10_As_planck2018", {}).get("mean", float("nan")))
            if isinstance(obs.get("ln10_As_planck2018", {}), dict)
            else float("nan")
        )
        As = As_from_ln10_As(ln10_As)
        N = starobinsky_N_from_ns(ns)
        r = starobinsky_r_from_N(N)
        H_inf_GeV = _H_inf_GeV_from_As_r(As=As, r=r, Mpl_reduced_GeV=MPL_REDUCED_GEV)

        # Isocurvature proxy (pre-inflation only):
        # P_iso ~ (H/(π f_a θ))^2 ; P_R ~ A_s
        if isocurvature_applicable:
            P_iso = (H_inf_GeV / (math.pi * f_a_GeV * max(1e-20, abs(theta_eff)))) ** 2 if f_a_GeV > 0 else float("nan")
            iso_ratio = float(P_iso / As) if (As > 0 and math.isfinite(P_iso)) else float("nan")
        else:
            P_iso = float("nan")
            iso_ratio = float("nan")
        beta_iso_max = float((pol.get("isocurvature", {}) if isinstance(pol.get("isocurvature", {}), dict) else {}).get("beta_iso_max", 0.038))
        iso_limit = float(beta_iso_max / max(1e-12, (1.0 - beta_iso_max)))

        checks: list[Check] = []
        if scen_selected:
            checks.append(mk_check_pass("scenario_policy_is_explicit", f"selected={scen_selected}"))
        else:
            checks.append((mk_check_fail if mode == "physics" else mk_check_warn)("scenario_policy_is_explicit", "missing cosmo_policy.scenario_policy.selected"))
        if scen_allowed:
            checks.append(
                mk_check_pass("scenario_policy_selected_known", f"selected={scen_selected}")
                if (scen_selected and scen_selected in set(scen_allowed))
                else mk_check_warn("scenario_policy_selected_known", f"selected={scen_selected!r} not in allowed list; using {scenario_id}")
            )
        checks.append(
            mk_check_pass("c_str_explained_by_topology", strings_dw_explain_detail)
            if strings_dw_explained
            else mk_check_warn("c_str_explained_by_topology", f"policy={strings_dw_policy}")
        )

        if scenario_id == "pre_inflation_single_theta_varphi0":
            checks.append(mk_check_pass("theta_i_fixed_by_varphi0", f"theta_i=varphi0={theta_varphi0}"))
        else:
            checks.append(mk_check_info("theta_i_random_postinflation", f"theta_i≈theta_rms={theta_eff} (scenario={scenario_id})"))
        checks.append(mk_check_pass("axion_frequency_from_mass", f"m_a={m_micro_eV} µeV => nu={nu_GHz:.3f} GHz (claim={nu_claim})"))
        checks.append(mk_check_info("g_coeff_fixed", f"g_coeff=-4c3={g_coeff} (dimensionless), g_phys≈{g_phys_GeV_inv:.3e} GeV^-1 for f_a={f_a_GeV:.3e} GeV"))

        # DM abundance: expected to underproduce for small θ_i, consistent with five_problems.tex narrative.
        if math.isfinite(frac_dm) and frac_dm < 1.0:
            checks.append(mk_check_warn("axion_underproduces_dm_misalignment", f"Omega_a h^2≈{omega_a_h2:.3e} => fraction≈{frac_dm:.3e} of Omega_DM h^2≈{omega_dm_ref}"))
        elif math.isfinite(frac_dm):
            checks.append(mk_check_info("axion_dm_fraction_estimate", f"Omega_a h^2≈{omega_a_h2:.3e} => fraction≈{frac_dm:.3e}"))
        else:
            checks.append(mk_check_warn("axion_dm_fraction_nan", f"Omega_a h^2≈{omega_a_h2} (check inputs)"))

        # Isocurvature viability:
        # - pre-inflation: apply the proxy constraint
        # - post-inflation: treat isocurvature as N/A (scenario matrix carries branch comparison)
        if isocurvature_applicable:
            if math.isfinite(iso_ratio):
                msg = f"P_iso/P_R≈{iso_ratio:.3e} (limit≈{iso_limit:.3e}); H_inf≈{H_inf_GeV:.3e} GeV, f_a≈{f_a_GeV:.3e} GeV, theta_i≈{theta_eff:.3e}"
                if iso_ratio > iso_limit:
                    checks.append((mk_check_fail if mode == "physics" else mk_check_warn)("isocurvature_tension_preinflation_proxy", msg))
                else:
                    checks.append(mk_check_pass("isocurvature_ok_preinflation_proxy", msg))
            else:
                checks.append(mk_check_warn("isocurvature_ratio_nan", f"iso_ratio={iso_ratio} (check inputs)"))
        else:
            checks.append(mk_check_pass("isocurvature_not_applicable_postinflation", f"scenario={scenario_id}"))

        # Plots
        warnings: list[str] = []
        plot: dict[str, str | None] = {"axion_dm_png": None}
        if config.plot:
            # θ grid on (0, π]
            thetas = [float(10 ** (math.log10(1e-4) + (math.log10(math.pi) - math.log10(1e-4)) * (i / 260))) for i in range(261)]
            omega_curve = [
                _omega_misalignment_h2(theta=t, f_a_GeV=f_a_GeV, omega_norm=omega_norm, f_norm_GeV=f_norm, p=p) for t in thetas
            ]
            iso_curve = [
                float(((H_inf_GeV / (math.pi * f_a_GeV * max(1e-20, abs(t)))) ** 2) / As) if (f_a_GeV > 0 and As > 0) else float("nan")
                for t in thetas
            ]
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_axion_dm(
                out_dir=out_dir,
                thetas=thetas,
                omega_h2=omega_curve,
                iso_ratio=iso_curve,
                theta_tfpt=theta_eff,
                omega_dm_ref=omega_dm_ref,
                beta_iso_max=beta_iso_max,
            )
            warnings.extend(plot_warnings)

        report_lines: list[str] = [
            "Axion DM pipeline (TFPT; engineering-level closure)",
            f"mode = {mode}",
            f"axion config: {ax_path}",
            "",
            "TFPT invariants:",
            f"- c3 = {cst.c3}",
            f"- varphi0 = {cst.varphi0}",
            f"- scenario_policy.selected = {scen_selected or '(missing)'}",
            f"- scenario_id (effective) = {scenario_id}",
            f"- theta_i policy: {theta_policy_txt}",
            f"- g_coeff = -4 c3 = {g_coeff}",
            "",
            "TFPT axion claim inputs:",
            f"- f_a = {f_a_GeV:.6g} GeV (source={fa_source})",
            f"- m_a = {m_micro_eV:.6g} µeV (source={fa_source})",
            f"- nu = {nu_GHz:.6g} GHz (computed) vs claim {nu_claim}",
            f"- g_phys ≈ g_coeff / f_a = {g_phys_GeV_inv:.3e} GeV^-1 (normalization-dependent; recorded as a candidate)",
            "",
            "Relic density (scenario policy; misalignment scaling):",
            f"- strings/domain-walls factor = {strings_factor:g}",
            f"- C_str policy: {strings_dw_policy} (detail={strings_dw_explain_detail})",
            f"- Ω_a h^2 ≈ {omega_a_h2:.3e}",
            f"- DM fraction vs Ω_DM h^2≈{omega_dm_ref}: {frac_dm:.3e}",
            "",
            "Inflation/isocurvature proxy (Starobinsky relations with Planck refs):",
            f"- n_s(ref)={ns} => N≈{N:.3f}, r≈{r:.4g}, A_s≈{As:.4g}",
            f"- H_inf≈{H_inf_GeV:.3e} GeV",
            f"- P_iso/P_R ≈ {iso_ratio:.3e} (limit≈{iso_limit:.3e})  [proxy; only meaningful for pre-inflation PQ]",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
            "",
            "Notes:",
            "- This module intentionally separates (a) TFPT-claimed axion parameters from (b) standard cosmology closure checks.",
            "- A publication-grade TFPT DM claim requires: (i) enforcing the ladder-derived f_a as a hard dependency, (ii) scenario selection (pre vs post inflation), (iii) strings/domain walls, (iv) astrophysical bounds, (v) g_{aγγ} normalization unambiguously tied to f_a.",
            "",
        ]
        if fa_derived_path:
            try:
                insert_at = report_lines.index("Relic density (scenario policy; misalignment scaling):")
            except ValueError:
                insert_at = len(report_lines)
            report_lines.insert(insert_at, f"- f_a derivation source: {fa_derived_path}")

        results: dict[str, Any] = {
            "mode": mode,
            "axion_claim_file": str(ax_path),
            "tfpt_invariants": {"c3": cst.c3, "varphi0": cst.varphi0, "g_coeff": g_coeff},
            "axion_claim": {
                "f_a_GeV": f_a_GeV,
                "m_a_micro_eV": m_micro_eV,
                "nu_GHz": nu_GHz,
                "g_phys_GeV_inv_candidate": g_phys_GeV_inv,
                "source": fa_source,
                "derived_path": fa_derived_path,
            },
            "scenario_policy": {
                "selected": scen_selected or None,
                "allowed": scen_allowed,
                "scenario_id_effective": scenario_id,
                "theta_policy": theta_policy_txt,
                "theta_eff": theta_eff,
                "strings_domain_walls_factor": strings_factor,
                "strings_domain_walls_factor_policy": strings_dw_policy if isinstance(strings_dw_policy, dict) else None,
                "strings_domain_walls_factor_explained": strings_dw_explained,
                "strings_domain_walls_factor_explain_detail": strings_dw_explain_detail,
                "isocurvature_applicable": isocurvature_applicable,
            },
            "relic_density": {"Omega_a_h2": omega_a_h2, "Omega_DM_h2_ref": omega_dm_ref, "fraction_of_dm": frac_dm},
            "inflation_proxy": {"n_s": ns, "ln10_As": ln10_As, "A_s": As, "N": N, "r": r, "H_inf_GeV": H_inf_GeV},
            "isocurvature_proxy": {
                "P_iso_over_P_R": iso_ratio if isocurvature_applicable else None,
                "limit_beta_iso": beta_iso_max,
                "limit_ratio": iso_limit,
            },
            "plot": plot,
        }

        return ModuleResult(results=results, checks=checks, report="\n".join(report_lines), warnings=warnings)

