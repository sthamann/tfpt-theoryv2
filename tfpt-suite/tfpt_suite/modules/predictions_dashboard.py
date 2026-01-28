from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.alpha_running import AlphaRunningInputs, alpha_bar5_MZ_from_alpha0, alpha_running_inputs_from_pdg
from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


@dataclass(frozen=True)
class Prediction:
    key: str
    label: str
    value: mp.mpf
    units: str
    status: str
    depends_on: list[str]


def _cabibbo_lambda_from_varphi0(varphi0: mp.mpf) -> mp.mpf:
    # Paper v2.5 (Cabibbo corollary): λ = sqrt(varphi0)*(1 - varphi0/2)
    return mp.sqrt(varphi0) * (mp.mpf(1) - mp.mpf(1) / 2 * varphi0)


def _starobinsky_triplet(*, M_over_Mpl: mp.mpf, N: int) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    # Paper v2.5 (Starobinsky): n_s=1-2/N, r=12/N^2, A_s≈N^2/(24π^2)*(M/Mpl)^2
    Nmp = mp.mpf(N)
    n_s = mp.mpf(1) - mp.mpf(2) / Nmp
    r = mp.mpf(12) / (Nmp**2)
    A_s = (Nmp**2) / (mp.mpf(24) * (mp.pi**2)) * (M_over_Mpl**2)
    return n_s, r, A_s


def _load_global_reference() -> dict[str, Any]:
    path = Path(__file__).resolve().parent.parent / "data" / "global_reference.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _z_score(*, pred: mp.mpf, mean: mp.mpf, sigma: mp.mpf, sigma_theory: mp.mpf = mp.mpf(0), sigma_floor: mp.mpf = mp.mpf(0)) -> mp.mpf:
    sig = mp.sqrt(sigma**2 + sigma_theory**2 + sigma_floor**2)
    if sig == 0:
        return mp.mpf("nan")
    return (pred - mean) / sig


def _ppm(*, pred: mp.mpf, ref: mp.mpf) -> mp.mpf:
    if ref == 0:
        return mp.mpf("nan")
    return (pred - ref) / ref * mp.mpf(1_000_000)


def _fmt_sci(x: mp.mpf | float, *, sig: int = 10) -> str:
    xf = float(x)
    if xf == 0.0:
        return "0"
    ax = abs(xf)
    if ax < 1e-3 or ax >= 1e4:
        return f"{xf:.{sig}e}"
    return f"{xf:.{sig}g}"


class PredictionsDashboardModule(TfptModule):
    module_id = "predictions_dashboard"
    title = "Predictions dashboard (paper-ready: 5–10 key numbers + uncertainties/dependencies)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants from paper v2.5 (c3, varphi0)",
                "Reference table (means/σ): tfpt_suite/data/global_reference.json",
                "Derived sectors: α (via alpha_precision_audit logic), Starobinsky R^2 (via M/Mpl), RG fingerprints (via two-loop gauge running table)",
            ],
            outputs=[
                "A compact predictions table (value + reference + z/ppm where applicable)",
                "Explicit dependencies (e.g. N for inflation; structural assumptions flags)",
            ],
            formulas=[
                "β_deg = (180/π) * varphi0/(4π)",
                "λ = sqrt(varphi0)*(1 - varphi0/2)",
                "Starobinsky: n_s=1-2/N, r=12/N^2, A_s≈N^2/(24π^2)*(M/Mpl)^2 with M/Mpl=sqrt(8π)c3^4",
            ],
            validation=[
                "All reported predictions are finite numbers",
                "Reference values are loaded and z-scores are computable when σ>0",
            ],
            determinism="Deterministic (algebraic; reference table is fixed input).",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()
        ref = _load_global_reference()
        obs = ref.get("observables", {})

        # --- Core predictions (parameter-free / structural) ---
        beta_deg = c.beta_deg
        cabibbo = _cabibbo_lambda_from_varphi0(c.varphi0)

        # --- α prediction (uses the same deterministic logic as alpha_precision_audit) ---
        # Import locally to avoid import-time cycles for registry loading.
        from tfpt_suite.modules.alpha_precision_audit import _fixed_point_alpha  # type: ignore

        sol_alpha = _fixed_point_alpha(c3=c.c3, b1=c.b1, varphi_tree=c.varphi0_tree, delta_top=c.delta_top, k=mp.mpf(2))
        sol_alpha_2def = _fixed_point_alpha(
            c3=c.c3, b1=c.b1, varphi_tree=c.varphi0_tree, delta_top=c.delta_top, k=mp.mpf(2), delta2=c.delta_top**2
        )

        # --- Inflation predictions (depends on N; we publish N=56 as benchmark + a small N-range) ---
        infl = []
        for N in (55, 56, 57):
            n_s, r, A_s = _starobinsky_triplet(M_over_Mpl=c.M_over_Mpl, N=N)
            infl.append({"N": N, "n_s": n_s, "r": r, "A_s": A_s})

        # Reference comparison helpers
        def get_ref(name: str) -> dict[str, Any]:
            return dict(obs.get(name, {}))

        alpha0_ref = get_ref("alpha_inv_0")
        alpha_bar5_ref = get_ref("alpha_bar5_inv_MZ")
        alpha_0_ref = get_ref("alpha_inv_0")
        beta_ref = get_ref("beta_deg")
        cab_ref = get_ref("cabibbo_lambda")
        ns_ref = get_ref("n_s")
        As_ref = get_ref("A_s")
        r_ref = get_ref("r_upper_95")

        # Secondary conversion model (external SM/QED running):
        pdg_path = Path(__file__).resolve().parent.parent / "data" / "alpha_running_pdg.json"
        pdg = json.loads(pdg_path.read_text(encoding="utf-8"))
        inputs = alpha_running_inputs_from_pdg(pdg=pdg)

        # compute secondary alpha_bar5(MZ) from TFPT alpha(0) using the two-defect refinement
        alpha0_pred = mp.mpf(1) / sol_alpha_2def.alpha_inv
        alpha_mz_pred, dal = alpha_bar5_MZ_from_alpha0(alpha0=alpha0_pred, inputs=inputs)
        alpha_bar5_inv_mz_pred = mp.mpf(1) / alpha_mz_pred

        # Compute reference comparisons (z, ppm) where applicable
        def compare_gaussian(pred: mp.mpf, ref_cfg: dict[str, Any]) -> dict[str, Any]:
            mean = mp.mpf(str(ref_cfg.get("mean", "nan")))
            sigma = mp.mpf(str(ref_cfg.get("sigma", "nan")))
            sigma_theory = mp.mpf(str(ref_cfg.get("sigma_theory", "0")))
            sigma_floor = mp.mpf(str(ref_cfg.get("sigma_floor", "0")))
            z = _z_score(pred=pred, mean=mean, sigma=sigma, sigma_theory=sigma_theory, sigma_floor=sigma_floor)
            return {
                "mean": mean,
                "sigma": sigma,
                "sigma_theory": sigma_theory,
                "sigma_floor": sigma_floor,
                "z": z,
                "ppm": _ppm(pred=pred, ref=mean),
            }

        comparisons: dict[str, Any] = {
            "alpha_bar5_inv_MZ_primary": compare_gaussian(alpha_bar5_inv_mz_pred, alpha_bar5_ref) if alpha_bar5_ref else {},
            "alpha_inv_0_selfconsistent": compare_gaussian(sol_alpha.alpha_inv, alpha0_ref) if alpha0_ref else {},
            "alpha_inv_0_two_defect": compare_gaussian(sol_alpha_2def.alpha_inv, alpha0_ref) if alpha0_ref else {},
            "beta_deg": compare_gaussian(beta_deg, beta_ref) if beta_ref else {},
            "cabibbo_lambda": compare_gaussian(cabibbo, cab_ref) if cab_ref else {},
            # inflation uses N=56 as the benchmark for ref comparisons
            "n_s_N56": compare_gaussian(infl[1]["n_s"], ns_ref) if ns_ref else {},
            "A_s_N56": compare_gaussian(infl[1]["A_s"], As_ref) if As_ref else {},
        }

        # r is typically compared to an upper bound
        r_cmp: dict[str, Any] = {}
        try:
            upper = mp.mpf(str(r_ref.get("upper")))
            r_pred = infl[1]["r"]
            r_cmp = {"upper_95": upper, "cl": r_ref.get("cl"), "passes": bool(r_pred < upper), "ratio": (r_pred / upper) if upper != 0 else mp.mpf("nan")}
        except Exception:
            r_cmp = {}

        preds: list[Prediction] = [
            Prediction(
                key="alpha_bar5_inv_MZ_primary",
                label="ᾱ^{(5)}(MZ)^{-1} (primary; MSbar-at-MZ via α(0)→MZ running)",
                value=alpha_bar5_inv_mz_pred,
                units="dimensionless",
                status="primary comparison observable under the MSbar-at-MZ policy (TFPT α(0) + external SM/QED running inputs)",
                depends_on=[
                    "TFPT α(0) (two-defect refinement)",
                    "Δα_had^(5)(MZ) (PDG input)",
                    "Δα_lept(MZ) (1-loop model)",
                ],
            ),
            Prediction(
                key="alpha_inv_0_selfconsistent",
                label="α^{-1}(0) (self-consistent, k=2; diagnostic)",
                value=sol_alpha.alpha_inv,
                units="dimensionless",
                status="diagnostic IR/on-shell quantity (CFE + double-cover backreaction; interpreted as α(0))",
                depends_on=["c3, varphi0_tree, delta_top, b1, backreaction exponent k=2"],
            ),
            Prediction(
                key="alpha_inv_0_two_defect",
                label="α^{-1}(0) (two-defect refinement, δ2=δ_top^2; diagnostic)",
                value=sol_alpha_2def.alpha_inv,
                units="dimensionless",
                status="diagnostic refinement (parameter-free next term template)",
                depends_on=["c3, varphi0_tree, delta_top, b1, k=2, δ2=δ_top^2"],
            ),
            Prediction(
                key="beta_deg",
                label="β (cosmic birefringence, degrees)",
                value=beta_deg,
                units="deg",
                status="structural (Δa_top = varphi0, n=1)",
                depends_on=["varphi0"],
            ),
            Prediction(
                key="cabibbo_lambda",
                label="λ (Cabibbo proxy)",
                value=cabibbo,
                units="dimensionless",
                status="structural (from varphi0)",
                depends_on=["varphi0"],
            ),
            Prediction(
                key="n_s_N56",
                label="n_s (Starobinsky R^2, N=56)",
                value=infl[1]["n_s"],
                units="dimensionless",
                status="derived (assumptions K1–K3; N choice)",
                depends_on=["c3 via M/Mpl", "N (e-folds)"],
            ),
            Prediction(
                key="r_N56",
                label="r (Starobinsky R^2, N=56)",
                value=infl[1]["r"],
                units="dimensionless",
                status="derived (assumptions K1–K3; N choice)",
                depends_on=["N (e-folds)"],
            ),
            Prediction(
                key="A_s_N56",
                label="A_s (Starobinsky R^2, N=56)",
                value=infl[1]["A_s"],
                units="dimensionless",
                status="derived (assumptions K1–K3; N choice)",
                depends_on=["c3 via M/Mpl", "N (e-folds)"],
            ),
        ]

        # Checks
        checks: list[Check] = []
        checks.append(
            Check(
                check_id="predictions_are_finite",
                passed=all(mp.isfinite(p.value) for p in preds),
                detail="all core prediction values are finite",
            )
        )
        checks.append(
            Check(
                check_id="reference_loaded",
                passed=bool(ref.get("vetted", False)) and bool(obs),
                detail="loaded tfpt_suite/data/global_reference.json",
            )
        )

        # Build report.txt
        lines: list[str] = []
        lines.append("TFPT predictions dashboard (paper v2.5)")
        lines.append("")
        lines.append("Key numbers (with status + dependencies):")
        for p in preds:
            lines.append(f"- {p.label}: {_fmt_sci(p.value)} [{p.units}]")
            lines.append(f"  status: {p.status}")
            lines.append(f"  depends on: {', '.join(p.depends_on)}")
            if p.key == "A_s_N56":
                lines.append(f"  (10^9 A_s = {float(p.value * mp.mpf('1e9')):.5f})")
        lines.append("")
        lines.append("Reference comparisons (from global_reference.json):")
        for k, v in comparisons.items():
            if not v:
                continue
            lines.append(
                f"- {k}: z={_fmt_sci(v.get('z', mp.mpf('nan')), sig=6)}, ppm={_fmt_sci(v.get('ppm', mp.mpf('nan')), sig=6)} "
                f"(mean={_fmt_sci(v.get('mean', mp.mpf('nan')))}, sigma={_fmt_sci(v.get('sigma', mp.mpf('nan')))})"
            )
        if r_cmp:
            lines.append("")
            lines.append(
                "r (N=56) vs 95% CL upper bound: "
                f"r={_fmt_sci(infl[1]['r'])}, upper={_fmt_sci(r_cmp.get('upper_95'))} -> PASS={r_cmp.get('passes')}, ratio={_fmt_sci(r_cmp.get('ratio'), sig=6)}"
            )
        lines.append("")
        lines.append("Inflation (N scan; report 10^9 A_s):")
        for row in infl:
            lines.append(
                f"- N={row['N']}: n_s={float(row['n_s']):.6f}, r={float(row['r']):.6f}, 10^9 A_s={float(row['A_s'] * mp.mpf('1e9')):.5f}"
            )
        report = "\n".join(lines) + "\n"

        return ModuleResult(
            results={
                "predictions": [
                    {"key": p.key, "label": p.label, "value": p.value, "units": p.units, "status": p.status, "depends_on": p.depends_on} for p in preds
                ],
                "comparisons": comparisons,
                "r_upper_bound_comparison": r_cmp,
                "inflation_scan": infl,
                "reference_source": "tfpt_suite/data/global_reference.json",
            },
            checks=checks,
            report=report,
            warnings=[],
        )

