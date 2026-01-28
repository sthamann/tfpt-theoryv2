from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from tfpt_suite.constants import TfptConstants
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


def _eta10_from_omega_b_h2(omega_b_h2: float) -> float:
    return float(273.9 * float(omega_b_h2))


def _kappa_eff_strong_washout(K: float) -> float:
    """
    Minimal deterministic efficiency proxy for thermal leptogenesis.

    For K>>1, a standard strong-washout scaling is roughly:
      κ ~ 0.3 / (K (ln K)^0.6)
    """
    Kf = float(K)
    if not math.isfinite(Kf) or Kf <= 0:
        return float("nan")
    if Kf < 3.0:
        return float(1.0 / max(2.0, Kf))  # very rough weak-washout proxy
    return float(0.3 / (Kf * (math.log(Kf) ** 0.6)))


class BaryogenesisMechanismModule(TfptModule):
    module_id = "baryogenesis_mechanism"
    title = "Baryogenesis mechanism (η_b; vanilla leptogenesis proxy from TFPT seesaw scales + CP suppression)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "seesaw scales: tfpt_suite/data/rge_thresholds_v25.json (MNR1..3; use M1=MNR1 for leptogenesis proxy)",
                "light neutrino masses input: tfpt_suite/data/flavor_texture_v24.json (m_i for Δm in Davidson–Ibarra bound)",
                "Planck reference: tfpt_suite/data/global_reference_minimal.json (Ω_b h^2 → η_b target)",
                "TFPT invariant: β_rad=varphi0/(4π) used as a discrete CP-suppression proxy",
            ],
            outputs=[
                "ε1 (CP asymmetry proxy) from Davidson–Ibarra bound × discrete CP suppression",
                "κ_eff (washout efficiency proxy) from K≈m3/m*",
                "η_b prediction and comparison to Planck η_b(Ω_b h^2)",
            ],
            formulas=[
                r"\epsilon_1^{\max} \approx \frac{3}{16\pi}\frac{M_1 (m_3-m_1)}{v^2} \;\; (\text{Davidson–Ibarra bound})",
                r"\epsilon_1 \approx \beta_{\rm rad}\,\epsilon_1^{\max} \;\; (\text{TFPT discrete CP-suppression proxy})",
                r"K \approx m_3/m_*,\; m_*\approx 1.08\times10^{-3}\,\mathrm{eV},\; \kappa\approx 0.3/[K(\ln K)^{0.6}]",
                r"\eta_b \approx 0.96\times 10^{-2}\,\epsilon_1\,\kappa \;\; (\text{order-of-magnitude leptogenesis proxy})",
                r"\eta_{10}\approx 273.9\,\Omega_b h^2,\;\; \eta_b=\eta_{10}\times10^{-10}",
            ],
            validation=[
                "Computes a deterministic η_b proxy (no continuous fit parameters) and checks consistency against the Planck Ω_b h^2 anchor.",
            ],
            determinism="Deterministic given the shipped input tables.",
            question="Is there an implemented, falsifiable baryogenesis mechanism that predicts η_b?",
            objective=[
                "Close the explicit baryogenesis placeholder with a concrete, auditable η_b proxy (vanilla leptogenesis).",
                "Keep publication-grade scope explicit: this is a deterministic proxy, not a full Boltzmann/leptogenesis solver.",
            ],
            gaps=[
                "Publication-grade requires a full flavor-aware leptogenesis computation (Boltzmann equations) and a rigorous CP source derived from the TFPT operator map.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")

        data_dir = Path(__file__).resolve().parent.parent / "data"
        thr_path = data_dir / "rge_thresholds_v25.json"
        tex_path = data_dir / "flavor_texture_v24.json"
        ref_path = data_dir / "global_reference_minimal.json"

        thr = json.loads(thr_path.read_text(encoding="utf-8")) if thr_path.is_file() else {}
        tex = json.loads(tex_path.read_text(encoding="utf-8")) if tex_path.is_file() else {}
        ref = json.loads(ref_path.read_text(encoding="utf-8")) if ref_path.is_file() else {}

        thresholds = thr.get("thresholds_GeV", {}) if isinstance(thr.get("thresholds_GeV", {}), dict) else {}
        M1_GeV = float(thresholds.get("MNR1", float("nan")))

        nm = (
            ((tex.get("neutrino_mechanism", {}) if isinstance(tex.get("neutrino_mechanism", {}), dict) else {}).get("neutrino_masses_input_eV", None))
        )
        m_eV = [float(x) for x in nm] if isinstance(nm, list) and len(nm) == 3 else [0.0, 0.0086, 0.05]
        m1_eV, m2_eV, m3_eV = float(m_eV[0]), float(m_eV[1]), float(m_eV[2])

        # Planck Ω_b h^2 -> η_b target
        obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
        omega_b = obs.get("omega_b_h2_planck2018", {}) if isinstance(obs.get("omega_b_h2_planck2018", {}), dict) else {}
        omega_b_mean = float(omega_b.get("mean", float("nan")))
        omega_b_sigma = float(omega_b.get("sigma", float("nan")))
        eta10_obs = _eta10_from_omega_b_h2(omega_b_mean) if math.isfinite(omega_b_mean) else float("nan")
        eta_b_obs = float(eta10_obs * 1e-10) if math.isfinite(eta10_obs) else float("nan")
        sigma_eta_obs = float(_eta10_from_omega_b_h2(omega_b_sigma) * 1e-10) if math.isfinite(omega_b_sigma) else float("nan")

        # Davidson–Ibarra bound (use v≈246 GeV, electroweak vev; engineering proxy).
        v_GeV = 246.0
        delta_m_eV = max(0.0, m3_eV - m1_eV)
        delta_m_GeV = float(delta_m_eV) * 1e-9
        eps1_max = (
            float((3.0 / (16.0 * math.pi)) * (float(M1_GeV) * float(delta_m_GeV) / (float(v_GeV) ** 2)))
            if (math.isfinite(M1_GeV) and M1_GeV > 0 and delta_m_GeV > 0)
            else float("nan")
        )

        # Discrete CP suppression proxy from TFPT invariants (β_rad = varphi0/(4π)).
        c = TfptConstants.compute()
        beta_rad = float(c.beta_rad)
        eps1 = float(beta_rad * eps1_max) if math.isfinite(eps1_max) else float("nan")

        # Washout proxy
        m_star_eV = 1.08e-3
        K = float(m3_eV / m_star_eV) if m_star_eV > 0 else float("nan")
        kappa = _kappa_eff_strong_washout(K)

        # η_b proxy
        eta_b_pred = float(0.96e-2 * eps1 * kappa) if (math.isfinite(eps1) and math.isfinite(kappa)) else float("nan")

        # Gate: require agreement within 0.5 dex (factor ~3) against the Planck η_b anchor.
        tol_dex = 0.5
        log10_ratio = float(math.log10(eta_b_pred / eta_b_obs)) if (eta_b_pred > 0 and eta_b_obs > 0) else float("nan")
        ok = bool(math.isfinite(log10_ratio) and abs(log10_ratio) <= tol_dex)

        # Also provide a conservative z-score using a large theory floor (proxy model).
        sigma_theory = float(0.5 * eta_b_obs) if (math.isfinite(eta_b_obs) and eta_b_obs > 0) else float("nan")
        sigma_eff2 = (sigma_eta_obs**2 if math.isfinite(sigma_eta_obs) else 0.0) + (sigma_theory**2 if math.isfinite(sigma_theory) else 0.0)
        z = float((eta_b_pred - eta_b_obs) / math.sqrt(sigma_eff2)) if (math.isfinite(eta_b_pred) and math.isfinite(eta_b_obs) and sigma_eff2 > 0) else float("nan")

        checks: list[Check] = []
        checks.append(mk_check_info("inputs_loaded", f"M1={M1_GeV:.3g} GeV, m_eV={m_eV}, beta_rad={beta_rad:.6g}"))
        checks.append(mk_check_info("leptogenesis_proxy", f"eps1_max≈{eps1_max:.3g}, eps1≈{eps1:.3g}, K≈{K:.3g}, kappa≈{kappa:.3g}"))
        checks.append(mk_check_info("eta_b_target", f"eta_b(obs)≈{eta_b_obs:.3g} from Ω_b h^2={omega_b_mean}"))
        if ok:
            checks.append(mk_check_pass("eta_b_match", f"eta_b≈{eta_b_pred:.3g} matches within {tol_dex} dex (log10 ratio≈{log10_ratio:.3g}, z≈{z:.3g})"))
        else:
            checks.append(
                (mk_check_fail if mode == "physics" else mk_check_warn)(
                    "eta_b_match",
                    f"eta_b≈{eta_b_pred:.3g} not within {tol_dex} dex of target {eta_b_obs:.3g} (log10 ratio≈{log10_ratio}, z≈{z})",
                )
            )

        lines = [
            "Baryogenesis mechanism (η_b; deterministic vanilla leptogenesis proxy)",
            "",
            f"mode={mode}",
            f"inputs: {thr_path} (MNR1), {tex_path} (m_i), {ref_path} (Ω_b h^2), TFPT β_rad",
            "",
            "Inputs:",
            f"- M1 = MNR1 ≈ {M1_GeV:.6g} GeV",
            f"- light neutrino masses (input) m_i = {m_eV} eV",
            f"- v ≈ {v_GeV} GeV (EW vev proxy)",
            "",
            "CP + washout proxies:",
            f"- ε1_max (Davidson–Ibarra) ≈ {eps1_max:.6g}",
            f"- β_rad (TFPT invariant) = {beta_rad:.6g} ⇒ ε1 ≈ β_rad·ε1_max ≈ {eps1:.6g}",
            f"- m* ≈ {m_star_eV:.3g} eV, K≈m3/m*≈{K:.3g} ⇒ κ_eff≈{kappa:.3g}",
            "",
            "η_b prediction:",
            f"- η_b ≈ 0.96e-2 · ε1 · κ_eff ≈ {eta_b_pred:.6g}",
            "",
            "η_b target (Planck Ω_b h^2 anchor):",
            f"- Ω_b h^2 = {omega_b_mean} ± {omega_b_sigma} ⇒ η10≈{eta10_obs:.6g} ⇒ η_b(obs)≈{eta_b_obs:.6g}",
            f"- log10(η_b/η_obs)≈{log10_ratio:.6g}, z(proxy)≈{z:.6g} (theory floor=50%)",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
            "",
            "Notes:",
            "- This is an engineering-level, deterministic leptogenesis proxy (no Boltzmann solver).",
            "- CP suppression uses β_rad as a TFPT-discrete proxy; a publication-grade derivation should tie ε1 to a topology→operator map.",
            "",
        ]

        results: dict[str, Any] = {
            "mode": mode,
            "inputs": {
                "rge_thresholds_file": str(thr_path),
                "flavor_texture_file": str(tex_path),
                "global_reference_minimal_file": str(ref_path),
                "M1_GeV": M1_GeV,
                "masses_input_eV": m_eV,
                "v_GeV": v_GeV,
            },
            "tfpt": {"beta_rad": beta_rad},
            "proxy": {
                "eps1_max": eps1_max,
                "eps1": eps1,
                "m_star_eV": m_star_eV,
                "K": K,
                "kappa_eff": kappa,
                "eta_b_pred": eta_b_pred,
            },
            "target": {
                "omega_b_h2_planck": {"mean": omega_b_mean, "sigma": omega_b_sigma},
                "eta_b_obs": eta_b_obs,
                "sigma_eta_obs": sigma_eta_obs,
            },
            "gate": {"tol_log10": tol_dex, "log10_ratio": log10_ratio, "z_proxy": z, "sigma_theory_floor": sigma_theory},
        }

        return ModuleResult(results=results, checks=checks, report="\n".join(lines), warnings=[])

