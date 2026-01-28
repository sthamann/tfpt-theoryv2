from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_info, mk_check_pass, mk_check_warn

CM_IN_GEV_INV = 5.0677307e13
FM_IN_GEV_INV = 5.0677307
GEV_TO_HZ = 2.417989242e23
K_MINIMAL_COUPLING = 0.75  # |b| ≈ (3/4)|S| (minimal coupling policy)
DEFAULT_PSD_SIGMA_REF_HZ = 1e-5
DEFAULT_PSD_F_REF_HZ = 1e-8
DEFAULT_PSD_GAMMA = 4.33
DEFAULT_PSD_F_SIGNAL_HZ = 1e-8
DEFAULT_PSD_BANDWIDTH_HZ = 1e-9

# Thermodynamic constants (for explicit magnetar polarization model).
KB_eV_per_K = 8.617333262e-5
MU_B_eV_per_T = 5.7883818060e-5


def _tfpt_M_GeV() -> float:
    cst = TfptConstants.compute()
    Mpl_red = 2.435e18
    return float(float(cst.M_over_Mpl) * Mpl_red)


def _n_cm3_to_GeV3(n_cm3: float) -> float:
    return float(float(n_cm3) * (1.0 / CM_IN_GEV_INV) ** 3)


def _n_fm3_to_GeV3(n_fm3: float) -> float:
    return float(float(n_fm3) * (1.0 / FM_IN_GEV_INV) ** 3)


def _sigma_from_frequency_psd(channel: dict[str, Any]) -> float:
    sigma_ref = float(channel.get("sigma_nu_ref_Hz", DEFAULT_PSD_SIGMA_REF_HZ))
    f_ref = float(channel.get("f_ref_Hz", DEFAULT_PSD_F_REF_HZ))
    gamma = float(channel.get("gamma", DEFAULT_PSD_GAMMA))
    f_signal = float(channel.get("f_signal_Hz", f_ref if f_ref > 0 else DEFAULT_PSD_F_SIGNAL_HZ))
    bandwidth = float(channel.get("bandwidth_Hz", DEFAULT_PSD_BANDWIDTH_HZ))
    if not (math.isfinite(sigma_ref) and math.isfinite(f_ref) and math.isfinite(gamma) and math.isfinite(f_signal) and math.isfinite(bandwidth)):
        return float("nan")
    if sigma_ref <= 0 or f_ref <= 0 or f_signal <= 0 or bandwidth <= 0:
        return float("nan")
    scale = (f_signal / f_ref) ** (-gamma / 2.0)
    return float(sigma_ref * scale * math.sqrt(bandwidth))


@dataclass(frozen=True)
class ChannelResult:
    channel_id: str
    source_id: str
    label: str
    delta_nu_Hz: float
    sigma_nu_Hz: float
    snr: float
    note: str


class TorsionFalsifiabilitySnrModule(TfptModule):
    module_id = "torsion_falsifiability_snr"
    title = "Torsion falsifiability (explicit source + noise model; SNR gate)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "noise + source model policy: tfpt_suite/data/torsion_falsifiability_noise_v1.json",
                "TFPT scale M from constants (R^2 scale) used as coupling-scale proxy",
                "minimal-coupling mapping constant k=3/4 (b_mu ≈ k S_mu)",
            ],
            outputs=[
                "SNR table for declared channels (lab + astro)",
                "go/no-go gate based on SNR ≥ 5 in at least one channel",
            ],
            formulas=[
                r"P_\mathrm{spin}=\tanh(\mu_B B/(k_B T)) \;\; (\mathrm{magnetar\ electron\ polarization\ proxy})",
                r"n_e = Y_e\, n_b,\;\;\rho_\mathrm{spin}\sim P\,s\,n_e,\;\; |S|\sim \rho_\mathrm{spin}/M_\mathrm{eff}^2,\;\; \Delta\nu\approx 2k|S|\cdot(\mathrm{GeV}\to\mathrm{Hz})",
                r"\sigma_\nu(f)=\sigma_{\nu,\mathrm{ref}}(f/f_\mathrm{ref})^{-\gamma/2}\sqrt{\Delta f}\;\;(\mathrm{frequency\ PSD\ proxy})",
                r"\mathrm{SNR} = |\Delta\nu|/\sigma_\nu",
            ],
            validation=[
                "Emits an explicit SNR gate for at least one channel (PASS in physics mode if SNR≥5).",
                "noise_psd_frequency_dependent: PASS if a frequency-dependent PSD noise model is declared.",
                "Keeps the source+noise assumptions as a machine-readable policy file.",
            ],
            determinism="Deterministic given the shipped policy JSON and TFPT constants.",
            question="Given an explicit source model and a noise model, is there at least one falsifiable torsion observable channel with SNR≥5?",
            objective=[
                "Upgrade torsion falsifiability from proxy thresholds to an explicit SNR calculation.",
                "Provide a concrete upgrade path to publication-grade likelihoods (replace σν proxies with real instrument PSD/likelihood).",
            ],
            gaps=[
                "This remains a proxy until a real dataset likelihood (timing residuals, polarimetry noise PSD, lab magnetometer systematics) is wired in.",
            ],
        )

    def run(self, config) -> ModuleResult:
        mode = str(getattr(config, "verification_mode", "engineering") or "engineering")
        data_dir = Path(__file__).resolve().parent.parent / "data"
        policy_path = data_dir / "torsion_falsifiability_noise_v1.json"
        policy = json.loads(policy_path.read_text(encoding="utf-8")) if policy_path.is_file() else {}

        # Load channel noise σν.
        channels = policy.get("channels", []) if isinstance(policy.get("channels", []), list) else []
        sigma_by_id: dict[str, float] = {}
        noise_by_id: dict[str, dict[str, Any]] = {}
        for ch in channels:
            if not isinstance(ch, dict) or "channel_id" not in ch:
                continue
            try:
                kind = str(ch.get("kind", "frequency_shift"))
                if kind == "frequency_psd":
                    sigma = _sigma_from_frequency_psd(ch)
                else:
                    sigma = float(ch.get("sigma_nu_Hz", float("nan")))
            except Exception:
                sigma = float("nan")
            channel_id = str(ch["channel_id"])
            sigma_by_id[channel_id] = sigma
            noise_by_id[channel_id] = {
                "kind": str(ch.get("kind", "frequency_shift")),
                "sigma_nu_Hz": sigma,
                "sigma_nu_ref_Hz": ch.get("sigma_nu_ref_Hz", None),
                "f_ref_Hz": ch.get("f_ref_Hz", None),
                "gamma": ch.get("gamma", None),
                "f_signal_Hz": ch.get("f_signal_Hz", None),
                "bandwidth_Hz": ch.get("bandwidth_Hz", None),
                "note": ch.get("note", None),
            }

        # --- Source model A: lab He-3 benchmark (kept for explicit "lab is too small" evidence) ---
        M_eff = _tfpt_M_GeV()
        spin_density_lab = float(0.7 * 0.5 * _n_cm3_to_GeV3(2.0e22))
        S_lab = float(abs(spin_density_lab) / (M_eff * M_eff)) if (math.isfinite(spin_density_lab) and M_eff > 0) else float("nan")
        dnu_lab = float(2.0 * (K_MINIMAL_COUPLING * S_lab) * GEV_TO_HZ) if math.isfinite(S_lab) else float("nan")

        # --- Source model B: magnetar electron polarization (derived from B,T via tanh) ---
        src_models = policy.get("source_models", []) if isinstance(policy.get("source_models", []), list) else []
        mag = next((s for s in src_models if isinstance(s, dict) and s.get("source_id") == "magnetar_electron_polarization_v1"), None)
        B_gauss = float(mag.get("B_gauss", 1e15) if isinstance(mag, dict) else 1e15)
        T_K = float(mag.get("T_kelvin", 1e6) if isinstance(mag, dict) else 1e6)
        nb_fm3 = float(mag.get("baryon_density_fm3", 0.16) if isinstance(mag, dict) else 0.16)
        Ye = float(mag.get("electron_fraction_Ye", 0.1) if isinstance(mag, dict) else 0.1)
        spin_per = float(mag.get("spin_per_particle", 0.5) if isinstance(mag, dict) else 0.5)

        B_T = float(B_gauss / 1.0e4)
        ratio = float((MU_B_eV_per_T * B_T) / (KB_eV_per_K * T_K)) if (T_K > 0 and math.isfinite(B_T)) else float("nan")
        pol_e = float(math.tanh(ratio)) if math.isfinite(ratio) else float("nan")
        ne_fm3 = float(max(0.0, Ye) * max(0.0, nb_fm3))
        ne_GeV3 = _n_fm3_to_GeV3(ne_fm3)
        spin_density_mag = float(pol_e * spin_per * ne_GeV3) if (math.isfinite(pol_e) and math.isfinite(ne_GeV3)) else float("nan")
        S_mag = float(abs(spin_density_mag) / (M_eff * M_eff)) if (math.isfinite(spin_density_mag) and M_eff > 0) else float("nan")
        dnu_mag = float(2.0 * (K_MINIMAL_COUPLING * S_mag) * GEV_TO_HZ) if math.isfinite(S_mag) else float("nan")

        rows: list[ChannelResult] = []

        def _snr(delta_nu: float, sigma: float) -> float:
            if not (math.isfinite(delta_nu) and math.isfinite(sigma) and sigma > 0):
                return float("nan")
            return float(abs(delta_nu) / sigma)

        # Channel: lab frequency proxy
        sig_lab = float(sigma_by_id.get("lab_frequency_shift_proxy", float("nan")))
        rows.append(
            ChannelResult(
                channel_id="lab_frequency_shift_proxy",
                source_id="lab_spin_fluid_He3",
                label="Lab (He-3 benchmark) — frequency shift proxy",
                delta_nu_Hz=dnu_lab,
                sigma_nu_Hz=sig_lab,
                snr=_snr(dnu_lab, sig_lab),
                note="Lab benchmark retained explicitly; expected to be far below sensitivity under TFPT-scale coupling.",
            )
        )

        # Channel: magnetar timing proxy
        sig_mag = float(sigma_by_id.get("magnetar_timing_proxy", float("nan")))
        rows.append(
            ChannelResult(
                channel_id="magnetar_timing_proxy",
                source_id="magnetar_electron_polarization_v1",
                label="Magnetar proxy (electron polarization) — timing/polarimetry frequency shift",
                delta_nu_Hz=dnu_mag,
                sigma_nu_Hz=sig_mag,
                snr=_snr(dnu_mag, sig_mag),
                note="Derived polarization P=tanh(μ_B B/(k_B T)) with n_e=Y_e n_b at nuclear density scale; explicit proxy σν.",
            )
        )

        snr_threshold = 5.0
        astro_snr = float(next((r.snr for r in rows if r.channel_id == "magnetar_timing_proxy"), float("nan")))
        lab_snr = float(next((r.snr for r in rows if r.channel_id == "lab_frequency_shift_proxy"), float("nan")))
        astro_ok = bool(math.isfinite(astro_snr) and astro_snr >= snr_threshold)
        any_ok = any(math.isfinite(r.snr) and r.snr >= snr_threshold for r in rows)

        checks: list[Check] = []
        checks.append(mk_check_pass("torsion_falsifiability_snr_table_present", f"rows={len(rows)} policy={policy_path}"))
        psd_channels = [cid for cid, meta in noise_by_id.items() if meta.get("kind") == "frequency_psd" and math.isfinite(float(meta.get("sigma_nu_Hz") or float("nan")))]
        checks.append(
            mk_check_pass("noise_psd_frequency_dependent", f"channels={psd_channels}")
            if psd_channels
            else mk_check_warn("noise_psd_frequency_dependent", "no frequency_psd channels declared")
        )
        checks.append(
            mk_check_pass("source_model_derived_not_benchmark", f"magnetar P=tanh(mu_B*B/(k_B*T)) => P≈{pol_e:.6g} (B={B_gauss:.3e} G, T={T_K:.3e} K, Ye={Ye})")
            if math.isfinite(pol_e)
            else mk_check_warn("source_model_derived_not_benchmark", f"magnetar polarization non-finite (B={B_gauss}, T={T_K})")
        )
        checks.append(
            mk_check_info("lab_channel_measurable_under_realistic_noise", f"SNR≈{lab_snr:.3e} (<{snr_threshold:g}; lab channel not measurable under current TFPT-scale coupling assumptions)")
        )
        checks.append(
            mk_check_pass("astro_channel_measurable_under_realistic_noise", f"SNR≈{astro_snr:.3g} (≥{snr_threshold:g})")
            if astro_ok
            else mk_check_warn("astro_channel_measurable_under_realistic_noise", f"SNR≈{astro_snr:.3g} (<{snr_threshold:g})")
        )
        if mode == "physics":
            checks.append(
                mk_check_pass("go_no_go_snr_ge_5", f"PASS: at least one channel has SNR≥{snr_threshold:g}")
                if any_ok
                else mk_check_warn("go_no_go_snr_ge_5", f"WARN: no channel reaches SNR≥{snr_threshold:g} (design incomplete)")
            )
        else:
            checks.append(mk_check_info("go_no_go_snr_ge_5", f"INFO: any_ok={any_ok} (threshold={snr_threshold:g})"))

        lines: list[str] = [
            "Torsion falsifiability (explicit source + noise model; SNR gate)",
            "",
            f"mode={mode}",
            f"policy={policy_path}",
            f"M_eff(TFPT_M)≈{M_eff:.3e} GeV, k=3/4, GeV→Hz={GEV_TO_HZ:.3e}",
            f"noise_policy: {noise_by_id}",
            "",
            "Derived source model (magnetar electron polarization):",
            f"- B={B_gauss:.3e} G ({B_T:.3e} T), T={T_K:.3e} K => μ_B B/(k_B T)≈{ratio:.3e} => P≈tanh(...)={pol_e:.6g}",
            f"- n_b≈{nb_fm3:.3g} fm^-3, Y_e≈{Ye:.3g} => n_e≈{ne_fm3:.3g} fm^-3 => {ne_GeV3:.3e} GeV^3",
            "",
            "Channels:",
        ]
        for r in rows:
            lines.append(f"- {r.channel_id}: Δν≈{r.delta_nu_Hz:.3e} Hz, σν≈{r.sigma_nu_Hz:.3e} Hz => SNR≈{r.snr:.3g} ({r.label})")
        lines += [
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "mode": mode,
                "policy_file": str(policy_path),
                "M_eff_GeV": M_eff,
                "snr_threshold": snr_threshold,
                "noise_policy": noise_by_id,
                "rows": [r.__dict__ for r in rows],
                "magnetar_source_model": {
                    "B_gauss": B_gauss,
                    "T_kelvin": T_K,
                    "ratio_muB_B_over_kBT": ratio,
                    "polarization": pol_e,
                    "baryon_density_fm3": nb_fm3,
                    "electron_fraction_Ye": Ye,
                    "n_e_fm3": ne_fm3,
                    "n_e_GeV3": ne_GeV3,
                    "spin_density_GeV3": spin_density_mag,
                    "S_abs_GeV": S_mag,
                    "delta_nu_Hz": dnu_mag,
                },
                "lab_source_model": {
                    "n_cm3": 2.0e22,
                    "polarization": 0.7,
                    "spin_per_particle": 0.5,
                    "spin_density_GeV3": spin_density_lab,
                    "S_abs_GeV": S_lab,
                    "delta_nu_Hz": dnu_lab,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

