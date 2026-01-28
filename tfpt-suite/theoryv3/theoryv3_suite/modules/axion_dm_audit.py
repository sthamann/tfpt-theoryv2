from __future__ import annotations

import json
from pathlib import Path

from mpmath import mp

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn
from theoryv3_suite.utils import coerce_float, ensure_ascii, load_tfpt_results


GHZ_PER_MICRO_EV = 0.24179893
OMEGA_DM_REF = 0.12


def _plot_axion_summary(
    *, out_dir: Path, omega_pred: float, omega_ref: float
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"axion_summary_png": None}
    warnings: list[str] = []
    try:
        import matplotlib.pyplot as plt  # type: ignore

        fig, ax = plt.subplots(figsize=(6.5, 3.2))
        ax.bar(["Omega_a h^2 pred", "Omega_DM h^2 ref"], [omega_pred, omega_ref], color=["#2f855a", "#c05621"])
        ax.set_ylabel("Omega h^2")
        ax.set_title("Axion DM relic fraction")
        ax.grid(True, axis="y", ls=":", alpha=0.4)

        fig.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "axion_summary.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["axion_summary_png"] = str(path)
    except Exception as exc:
        warnings.append(f"plot_generation_failed: {exc}")
    return plot, warnings


class AxionDmAuditModule(TfptModule):
    module_id = "axion_dm_audit"
    title = "Axion DM audit (frequency, relic fraction)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "axion_dm_pipeline outputs (preferred)",
                "axion_tfpt_v106.json (fallback)",
            ],
            outputs=[
                "nu_GHz from m_a",
                "Omega_a h^2",
                "theta_eff and strings/domain-walls factor (if available)",
            ],
            formulas=[
                "nu_GHz = 0.24179893 * m_a_micro_eV",
            ],
            validation=[
                "frequency matches the m_a conversion",
                "relic fraction is near Omega_DM reference",
            ],
            determinism="Deterministic given inputs.",
            question="Do the axion numbers map cleanly to a frequency and relic fraction?",
            objective=[
                "Expose the axion target frequency and relic fraction in a single audit.",
            ],
        )

    def run(self, config) -> ModuleResult:
        payload = load_tfpt_results("axion_dm_pipeline", prefer_physics=True)
        axion_claim: dict[str, float] = {}
        relic: dict[str, float] = {}
        scenario: dict[str, float | str] = {}

        if payload:
            axion_claim_raw = payload.get("results", {}).get("axion_claim", {})
            axion_claim = {
                "f_a_GeV": coerce_float(axion_claim_raw.get("f_a_GeV")),
                "m_a_micro_eV": coerce_float(axion_claim_raw.get("m_a_micro_eV")),
                "nu_GHz": coerce_float(axion_claim_raw.get("nu_GHz")),
            }
            relic_raw = payload.get("results", {}).get("relic_density", {})
            relic = {
                "Omega_a_h2": coerce_float(relic_raw.get("Omega_a_h2")),
                "fraction_of_dm": coerce_float(relic_raw.get("fraction_of_dm")),
            }
            scenario_raw = payload.get("results", {}).get("scenario_policy", {})
            scenario = {
                "theta_eff": coerce_float(scenario_raw.get("theta_eff")),
                "strings_domain_walls_factor": coerce_float(scenario_raw.get("strings_domain_walls_factor")),
                "scenario_id": str(scenario_raw.get("scenario_id_effective", "")),
            }
        else:
            data_dir = Path(__file__).resolve().parents[3] / "tfpt_suite" / "data"
            ax_path = data_dir / "axion_tfpt_v106.json"
            ax = json.loads(ax_path.read_text(encoding="utf-8"))
            axion_claim = {
                "f_a_GeV": coerce_float(ax.get("axion", {}).get("f_a_GeV")),
                "m_a_micro_eV": coerce_float(ax.get("axion", {}).get("m_a_micro_eV")),
            }
            axion_claim["nu_GHz"] = GHZ_PER_MICRO_EV * axion_claim["m_a_micro_eV"]
            relic = {"Omega_a_h2": float("nan"), "fraction_of_dm": float("nan")}
            scenario = {"theta_eff": float("nan"), "strings_domain_walls_factor": float("nan"), "scenario_id": "n/a"}

        nu_calc = GHZ_PER_MICRO_EV * axion_claim.get("m_a_micro_eV", float("nan"))
        nu_reported = axion_claim.get("nu_GHz", float("nan"))
        freq_ok = abs(nu_calc - nu_reported) < 1e-3 if (nu_calc == nu_calc and nu_reported == nu_reported) else False

        checks: list[Check] = []
        checks.append(
            mk_check_pass("frequency_matches", f"nu_calc={nu_calc}, nu_reported={nu_reported}")
            if freq_ok
            else mk_check_warn("frequency_matches", f"nu_calc={nu_calc}, nu_reported={nu_reported}")
        )
        if relic.get("Omega_a_h2", float("nan")) == relic.get("Omega_a_h2", float("nan")):
            checks.append(
                mk_check_pass("omega_a_near_ref", f"Omega_a_h2={relic['Omega_a_h2']}")
                if abs(relic["Omega_a_h2"] - OMEGA_DM_REF) <= 0.03
                else mk_check_warn("omega_a_near_ref", f"Omega_a_h2={relic['Omega_a_h2']}")
            )

        lines = [
            "Axion DM audit",
            "",
            f"f_a_GeV = {axion_claim.get('f_a_GeV')}",
            f"m_a_micro_eV = {axion_claim.get('m_a_micro_eV')}",
            f"nu_GHz (reported) = {nu_reported}",
            f"nu_GHz (calc) = {nu_calc}",
            "",
            f"Omega_a_h2 = {relic.get('Omega_a_h2')} (ref={OMEGA_DM_REF})",
            f"fraction_of_dm = {relic.get('fraction_of_dm')}",
            f"theta_eff = {scenario.get('theta_eff')}",
            f"strings_domain_walls_factor = {scenario.get('strings_domain_walls_factor')}",
            f"scenario_id = {scenario.get('scenario_id')}",
            "",
            "Checks:",
            *[
                f"- {c.check_id}: {str(c.severity or ('PASS' if c.passed else 'FAIL')).upper()} ({ensure_ascii(c.detail)})"
                for c in checks
            ],
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"axion_summary_png": None}
        if getattr(config, "plot", True) and relic.get("Omega_a_h2", float("nan")) == relic.get("Omega_a_h2", float("nan")):
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_axion_summary(out_dir=out_dir, omega_pred=relic["Omega_a_h2"], omega_ref=OMEGA_DM_REF)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "axion_claim": axion_claim,
                "relic": relic,
                "scenario": scenario,
                "nu_calc_GHz": nu_calc,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )
