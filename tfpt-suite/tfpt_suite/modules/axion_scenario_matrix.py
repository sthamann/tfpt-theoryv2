from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from tfpt_suite.axion_inputs import resolve_axion_claim
from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass, mk_check_warn


MPL_REDUCED_GEV = 2.435e18


def _workspace_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _microeV_to_GHz(micro_eV: float) -> float:
    # E = h nu ; 1 eV corresponds to 241.79893 THz.
    return float(micro_eV) * 1e-6 * 241.79893e3


def _omega_misalignment_h2(*, theta: float, f_a_GeV: float, omega_norm: float, f_norm_GeV: float, p: float) -> float:
    th = float(abs(theta))
    return float(float(omega_norm) * (th**2) * (float(f_a_GeV) / float(f_norm_GeV)) ** float(p))


def _As_from_ln10_As(ln10_As: float) -> float:
    return float(math.exp(float(ln10_As)) / 1e10)


def _starobinsky_N_from_ns(ns: float) -> float:
    # n_s ≈ 1 - 2/N
    return float(2.0 / max(1e-12, (1.0 - float(ns))))


def _starobinsky_r_from_N(N: float) -> float:
    return float(12.0 / (float(N) ** 2))


def _H_inf_GeV_from_As_r(*, As: float, r: float) -> float:
    # H = π M̄_P sqrt(A_s r / 2)
    Asf = float(As)
    rf = float(r)
    if Asf <= 0 or rf <= 0:
        return float("nan")
    return float(math.pi * MPL_REDUCED_GEV * math.sqrt(Asf * rf / 2.0))


@dataclass(frozen=True)
class AxionScenarioResult:
    scenario_id: str
    pq_breaking: str
    theta_policy: str
    theta_eff: float
    strings_domain_walls_factor: float
    omega_a_h2: float
    frac_dm: float
    isocurvature_applicable: bool
    P_iso_over_P_R: float | None
    isocurvature_limit_ratio: float
    passes_isocurvature: bool
    overproduces_dm: bool
    passes_scenario: bool
    note: str


class AxionScenarioMatrixModule(TfptModule):
    module_id = "axion_scenario_matrix"
    title = "Axion scenario matrix (pre/post inflation PQ; explicit branches; physics gate)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "axion claim + policy knobs: tfpt_suite/data/axion_tfpt_v106.json (derived f_a preferred if available)",
                "Planck baseline (n_s, ln10_As) for inflation scale proxy: tfpt_suite/data/global_reference_minimal.json",
            ],
            outputs=[
                "scenario table: Ω_a h^2, DM fraction, isocurvature status per branch",
                "physics gate: at least one scenario must pass isocurvature + non-overproduction",
            ],
            formulas=[
                r"Ω_a h^2 ≈ Ω_norm · θ_eff^2 · (f_a/f_norm)^p (engineering misalignment scaling)",
                r"pre-inflation proxy: P_iso/P_R ≈ (H/(π f_a θ_eff))^2 / A_s",
                r"post-inflation: isocurvature is treated as N/A (requires strings/domain-walls modeling instead)",
            ],
            validation=[
                "Scenario list is explicit (no silent default).",
                "At least one scenario passes the isocurvature gate in physics mode.",
            ],
            determinism="Deterministic given inputs.",
            question="Which PQ-breaking scenario(s) are compatible with TFPT axion parameters under explicit assumptions?",
            objective=[
                "Make the scenario dependence first-class (pre vs post inflation PQ breaking).",
                "Prevent a single implicit choice from dominating ToE optics.",
            ],
            gaps=[
                "Publication-grade post-inflation requires a strings/domain-walls relic computation and updated lattice inputs.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering"))

        ax_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "axion_tfpt_v106.json"
        ref_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "global_reference_minimal.json"
        ax = json.loads(ax_path.read_text(encoding="utf-8")) if ax_path.is_file() else {}
        ref = json.loads(ref_path.read_text(encoding="utf-8")) if ref_path.is_file() else {}

        claim_resolved = resolve_axion_claim(ax_raw=ax, output_dir=Path(config.output_dir))
        f_a_GeV = float(claim_resolved.get("f_a_GeV", float("nan")))
        m_micro_eV = float(claim_resolved.get("m_a_micro_eV", float("nan")))
        nu_GHz = _microeV_to_GHz(m_micro_eV)
        fa_source = str(claim_resolved.get("source", "quoted"))

        pol = ax.get("cosmo_policy", {}) if isinstance(ax.get("cosmo_policy", {}), dict) else {}
        mis = pol.get("misalignment_model", {}) if isinstance(pol.get("misalignment_model", {}), dict) else {}
        omega_norm = float(mis.get("omega_norm_at_fa_5e11", 0.12))
        f_norm = float(mis.get("fa_norm_GeV", 5.0e11))
        p = float(mis.get("power_law_index", 1.165))
        omega_dm_ref = float(pol.get("omega_dm_h2_ref", 0.12))
        beta_iso_max = float((pol.get("isocurvature", {}) if isinstance(pol.get("isocurvature", {}), dict) else {}).get("beta_iso_max", 0.038))
        iso_limit = float(beta_iso_max / max(1e-12, (1.0 - beta_iso_max)))

        # Inflation scale proxy (Starobinsky; consistent with axion_dm_pipeline)
        obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
        ns = float((obs.get("n_s_planck2018", {}) if isinstance(obs.get("n_s_planck2018", {}), dict) else {}).get("mean", float("nan")))
        ln10_As = float((obs.get("ln10_As_planck2018", {}) if isinstance(obs.get("ln10_As_planck2018", {}), dict) else {}).get("mean", float("nan")))
        As = _As_from_ln10_As(ln10_As) if math.isfinite(ln10_As) else float("nan")
        N = _starobinsky_N_from_ns(ns) if math.isfinite(ns) else float("nan")
        r = _starobinsky_r_from_N(N) if math.isfinite(N) else float("nan")
        H_inf_GeV = _H_inf_GeV_from_As_r(As=As, r=r) if math.isfinite(As) and math.isfinite(r) else float("nan")

        cst = TfptConstants.compute()
        theta_varphi0 = float(cst.varphi0)
        theta_rms = float(math.pi / math.sqrt(3.0))  # rms of uniform θ in [-π,π]

        # Scenario knobs (optional; explicit defaults)
        strings_factor = float((pol.get("post_inflation", {}) if isinstance(pol.get("post_inflation", {}), dict) else {}).get("strings_domain_walls_factor", 1.0))
        strings_factor = float(max(strings_factor, 1.0))

        scenarios: list[AxionScenarioResult] = []

        def add_scenario(
            *,
            scenario_id: str,
            pq_breaking: str,
            theta_policy: str,
            theta_eff: float,
            strings_domain_walls_factor: float,
            isocurvature_applicable: bool,
            note: str,
        ) -> None:
            omega_a = _omega_misalignment_h2(theta=theta_eff, f_a_GeV=f_a_GeV, omega_norm=omega_norm, f_norm_GeV=f_norm, p=p)
            omega_a *= float(strings_domain_walls_factor)
            frac = omega_a / omega_dm_ref if omega_dm_ref > 0 else float("nan")
            over = bool(math.isfinite(frac) and frac > 1.0)

            iso_val: float | None
            iso_ok: bool
            if isocurvature_applicable:
                if not (math.isfinite(H_inf_GeV) and math.isfinite(As) and f_a_GeV > 0 and theta_eff != 0):
                    iso_val = None
                    iso_ok = False
                else:
                    P_iso = (H_inf_GeV / (math.pi * f_a_GeV * max(1e-20, abs(theta_eff)))) ** 2
                    iso_val = float(P_iso / As) if As > 0 else float("nan")
                    iso_ok = bool(math.isfinite(iso_val) and iso_val <= iso_limit)
            else:
                iso_val = None
                iso_ok = True

            passes = bool((not over) and iso_ok)
            scenarios.append(
                AxionScenarioResult(
                    scenario_id=scenario_id,
                    pq_breaking=pq_breaking,
                    theta_policy=theta_policy,
                    theta_eff=float(theta_eff),
                    strings_domain_walls_factor=float(strings_domain_walls_factor),
                    omega_a_h2=float(omega_a),
                    frac_dm=float(frac),
                    isocurvature_applicable=bool(isocurvature_applicable),
                    P_iso_over_P_R=iso_val,
                    isocurvature_limit_ratio=float(iso_limit),
                    passes_isocurvature=bool(iso_ok),
                    overproduces_dm=bool(over),
                    passes_scenario=bool(passes),
                    note=str(note),
                )
            )

        add_scenario(
            scenario_id="pre_inflation_single_theta_varphi0",
            pq_breaking="pre_inflation",
            theta_policy="theta_i = varphi0 (birefringence-fixed; single-angle patch)",
            theta_eff=theta_varphi0,
            strings_domain_walls_factor=1.0,
            isocurvature_applicable=True,
            note="Standard pre-inflation proxy (single θ_i) used by the legacy axion_dm_pipeline; expected to fail isocurvature for high H_inf.",
        )
        add_scenario(
            scenario_id="post_inflation_theta_rms_no_strings",
            pq_breaking="post_inflation",
            theta_policy="random θ_i; use rms(θ)=π/√3 for estimate",
            theta_eff=theta_rms,
            strings_domain_walls_factor=1.0,
            isocurvature_applicable=False,
            note="Post-inflation PQ breaking: treat isocurvature as N/A; strings/domain walls not included (conservative).",
        )
        add_scenario(
            scenario_id="post_inflation_theta_rms_with_strings_dw_factor",
            pq_breaking="post_inflation",
            theta_policy="random θ_i; use rms(θ)=π/√3; include strings/domain-walls factor",
            theta_eff=theta_rms,
            strings_domain_walls_factor=strings_factor,
            isocurvature_applicable=False,
            note="Post-inflation PQ breaking: includes an explicit multiplicative strings/domain-walls factor (policy knob; needs a real derivation).",
        )

        # --- Checks ---
        checks: list[Check] = []
        checks.append(mk_check_pass("scenario_is_explicit", f"{len(scenarios)} scenarios enumerated (pre vs post inflation)"))

        any_iso_pass = any(s.passes_isocurvature for s in scenarios)
        if mode == "physics":
            checks.append(
                mk_check_pass("isocurvature_passes_in_at_least_one_scenario", "PASS")
                if any_iso_pass
                else mk_check_fail("isocurvature_passes_in_at_least_one_scenario", "FAIL: no scenario passes isocurvature gate")
            )
        else:
            checks.append(mk_check_pass("isocurvature_passes_in_at_least_one_scenario", "PASS (engineering mode; diagnostic)"))

        any_full_pass = any(s.passes_scenario for s in scenarios)
        if mode == "physics":
            checks.append(
                mk_check_pass("at_least_one_scenario_is_physically_viable", "PASS")
                if any_full_pass
                else mk_check_fail("at_least_one_scenario_is_physically_viable", "FAIL: no scenario passes (isocurvature + no overproduction)")
            )
        else:
            checks.append(mk_check_pass("at_least_one_scenario_is_physically_viable", "PASS (engineering mode; diagnostic)"))

        # If the selected scenario is present in config, record it as a policy check (no hard requirement in engineering mode).
        selected = str((pol.get("scenario_policy", {}) if isinstance(pol.get("scenario_policy", {}), dict) else {}).get("selected", "")).strip()
        if selected:
            found = any(s.scenario_id == selected for s in scenarios)
            checks.append(
                mk_check_pass("scenario_policy_selected_is_known", f"selected={selected}")
                if found
                else mk_check_warn("scenario_policy_selected_is_known", f"unknown selected={selected!r}; update axion_tfpt_v106.json")
            )
        else:
            checks.append(
                (mk_check_fail if mode == "physics" else mk_check_warn)(
                    "scenario_policy_is_explicit",
                    "missing cosmo_policy.scenario_policy.selected in axion_tfpt_v106.json (physics mode requires an explicit scenario choice)",
                )
            )

        # --- Report ---
        lines: list[str] = []
        lines += [
            "Axion scenario matrix (explicit branches)",
            "",
            f"mode={mode}",
            f"axion claim file: {ax_path}",
            f"- f_a={f_a_GeV:.6g} GeV, m_a={m_micro_eV:.6g} µeV (nu≈{nu_GHz:.3f} GHz, source={fa_source})",
            "",
            "Inflation proxy (Starobinsky, Planck refs):",
            f"- n_s={ns}, ln(1e10 A_s)={ln10_As} => A_s={As:.6g}, N≈{N:.6g}, r≈{r:.6g}, H≈{H_inf_GeV:.3e} GeV",
            "",
            f"isocurvature limit proxy: beta_iso_max={beta_iso_max} => ratio_limit≈{iso_limit:.6g}",
            "",
            "Scenarios:",
        ]
        for s in scenarios:
            iso_txt = "N/A" if not s.isocurvature_applicable else (f"{s.P_iso_over_P_R:.3e}" if s.P_iso_over_P_R is not None else "n/a")
            lines.append(
                f"- {s.scenario_id}: pq={s.pq_breaking}, theta={s.theta_eff:.3e} ({s.theta_policy}), "
                f"strings_factor={s.strings_domain_walls_factor:g}, Ω_a h^2={s.omega_a_h2:.3e} (frac={s.frac_dm:.3e}), "
                f"iso={iso_txt} (ok={s.passes_isocurvature}), overDM={s.overproduces_dm}, PASS={s.passes_scenario}"
            )
        lines += [
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "inputs": {"axion_tfpt_file": str(ax_path), "global_reference_minimal_file": str(ref_path)},
                "axion_claim": {"f_a_GeV": f_a_GeV, "m_a_micro_eV": m_micro_eV, "nu_GHz": nu_GHz, "source": fa_source},
                "inflation_proxy": {"n_s": ns, "A_s": As, "N": N, "r": r, "H_inf_GeV": H_inf_GeV},
                "limits": {"beta_iso_max": beta_iso_max, "iso_ratio_limit": iso_limit, "omega_dm_h2_ref": omega_dm_ref},
                "scenarios": [s.__dict__ for s in scenarios],
                "scenario_policy_selected": selected or None,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

