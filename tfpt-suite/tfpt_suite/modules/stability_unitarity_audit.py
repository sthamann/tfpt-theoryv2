from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np

from tfpt_suite.conventions import gY_from_g1_gut
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.pyrate_pythonoutputs import get_pyrate_pythonoutput
from tfpt_suite.rge_pyrate_2loop import load_pyrate_beta_module, run_flavor_rge_2loop_thresholds
from tfpt_suite.pyrate_boundary_runner import sm_boundary_conditions_at_mt


def _workspace_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _relpath(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(_workspace_root()))
    except Exception:
        return str(path)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass(frozen=True)
class LambdaCrossing:
    crosses: bool
    mu_cross_GeV: Optional[float]
    note: str


class StabilityUnitarityAuditModule(TfptModule):
    module_id = "stability_unitarity_audit"
    title = "Stability & perturbativity audit (λ(μ) + coupling red flags; mt→μUV)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "SM inputs at MZ: tfpt_suite/data/sm_inputs_mz.json (boundary conditions at mt)",
                "RG thresholds: tfpt_suite/data/rge_thresholds_v25.json",
                "2-loop PyR@TE betas (SM and E8 σ+yN)",
            ],
            outputs=[
                "λ(mt), λ(μUV) and a coarse instability-scale estimate (λ=0 crossing, if any)",
                "perturbativity red-flag summary (max |g|, |y|, |λ| at mt and μUV)",
            ],
            formulas=[
                "Vacuum stability proxy: track λ(μ). If λ crosses 0 at some μ<μUV, the EW vacuum is metastable/unstable (needs lifetime analysis for publication).",
                "Perturbativity proxy: require all dimensionless couplings remain finite and below a conservative ceiling (e.g. 4π).",
            ],
            validation=[
                "RG integration succeeds and outputs are finite.",
                "If λ(mt)>0 and λ(μUV)<0, a bracketed λ=0 crossing scale is found via bisection (coarse).",
            ],
            determinism="Deterministic given inputs (no random sampling).",
        )

    def run(self, config) -> ModuleResult:
        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"
        thr_path = Path(__file__).resolve().parent.parent / "data" / "rge_thresholds_v25.json"
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        thr_raw = json.loads(thr_path.read_text(encoding="utf-8"))
        thresholds_GeV = dict(thr_raw.get("thresholds_GeV", {}))
        mu_uv_GeV = float(thr_raw.get("mu_uv_GeV", 1.0e16))

        bc = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_raw)
        mu_start_GeV = float(bc.mu_mt_GeV)
        g1_gut_mt = float(bc.route_2loop["g1_gut_mt"])
        gY_mt = float(bc.route_2loop["gY_mt"])
        g2_mt = float(bc.route_2loop["g2_mt"])
        g3_mt = float(bc.route_2loop["g3_mt"])
        lam_mt = float(bc.route_2loop["lambda_mt"])
        yt_mt = float(bc.route_2loop["yt_mt"])
        yb_mt = float(bc.route_2loop["yb_mt"])
        ytau_mt = float(bc.route_2loop["ytau_mt"])

        # Minimal diagonal Yukawa matrices (phase-1): keep only 3rd generation.
        Yu_mt = np.diag([0.0, 0.0, yt_mt]).astype(complex)
        Yd_mt = np.diag([0.0, 0.0, yb_mt]).astype(complex)
        Ye_mt = np.diag([0.0, 0.0, ytau_mt]).astype(complex)

        py_sm = get_pyrate_pythonoutput("sm_tfpt_2loop_v25")
        py_e8 = get_pyrate_pythonoutput("e8_sigma_yN_2loop")
        beta_sm = load_pyrate_beta_module(
            kind="sm_tfpt_2loop_v25",
            pythonoutput_dir=py_sm.pythonoutput_dir,
            model_name_expected=py_sm.model_name_expected,
            yaml_source=py_sm.yaml_source,
        )
        beta_e8 = load_pyrate_beta_module(
            kind="e8_sigma_yN_2loop",
            pythonoutput_dir=py_e8.pythonoutput_dir,
            model_name_expected=py_e8.model_name_expected,
            yaml_source=py_e8.yaml_source,
        )

        # Primary run (mt→μUV)
        rge = run_flavor_rge_2loop_thresholds(
            mu_start_GeV=mu_start_GeV,
            mu_end_GeV=mu_uv_GeV,
            thresholds_GeV=thresholds_GeV,
            g_start=(gY_mt, g2_mt, g3_mt),
            Yu_start=Yu_mt,
            Yd_start=Yd_mt,
            Ye_start=Ye_mt,
            lambda_start=float(lam_mt),
            yN_start=None,
            beta_sm=beta_sm,
            beta_e8=beta_e8,
            apply_sigma_threshold=True,
            apply_g8_delta_b3=True,
            delta_b3_g8=2.0,
            apply_matching=True,
            matching_loop_order=1,
            rtol=1e-8,
            atol=1e-10,
            method="DOP853",
        )
        lam_uv = float(rge["lambda_end"])
        Yu_uv = np.array(rge["Yu_end"], dtype=complex)
        Yd_uv = np.array(rge["Yd_end"], dtype=complex)
        Ye_uv = np.array(rge["Ye_end"], dtype=complex)
        gY_uv = float(rge["g_end_sm"]["gY"])
        g2_uv = float(rge["g_end_sm"]["g2"])
        g3_uv = float(rge["g_end_sm"]["g3"])

        def maxabs(mat: np.ndarray) -> float:
            return float(np.max(np.abs(np.array(mat, dtype=complex))))

        # Coarse λ=0 crossing scale (bisection in log μ), only if we have a sign change.
        crossing = LambdaCrossing(crosses=False, mu_cross_GeV=None, note="no sign change bracketed")
        if np.isfinite(lam_mt) and np.isfinite(lam_uv) and (lam_mt > 0) and (lam_uv < 0):
            lo = float(mu_start_GeV)
            hi = float(mu_uv_GeV)

            def lam_at(mu: float) -> float:
                out = run_flavor_rge_2loop_thresholds(
                    mu_start_GeV=mu_start_GeV,
                    mu_end_GeV=float(mu),
                    thresholds_GeV=thresholds_GeV,
                    g_start=(gY_mt, g2_mt, g3_mt),
                    Yu_start=Yu_mt,
                    Yd_start=Yd_mt,
                    Ye_start=Ye_mt,
                    lambda_start=float(lam_mt),
                    yN_start=None,
                    beta_sm=beta_sm,
                    beta_e8=beta_e8,
                    apply_sigma_threshold=True,
                    apply_g8_delta_b3=True,
                    delta_b3_g8=2.0,
                    apply_matching=True,
                    matching_loop_order=1,
                    rtol=2e-7,
                    atol=2e-9,
                    method="DOP853",
                )
                return float(out["lambda_end"])

            lam_lo = float(lam_mt)
            lam_hi = float(lam_uv)
            # 20 bisection iterations is plenty for a coarse log-scale.
            for _ in range(20):
                mid = math.sqrt(lo * hi)
                lam_mid = lam_at(mid)
                if not np.isfinite(lam_mid):
                    break
                if lam_mid > 0:
                    lo, lam_lo = mid, lam_mid
                else:
                    hi, lam_hi = mid, lam_mid
            crossing = LambdaCrossing(crosses=True, mu_cross_GeV=float(math.sqrt(lo * hi)), note="coarse bisection bracket")

        # Perturbativity proxy (very conservative ceiling)
        ceiling = float(4.0 * np.pi)
        max_dimless_mt = max(abs(gY_mt), abs(g2_mt), abs(g3_mt), abs(lam_mt), maxabs(Yu_mt), maxabs(Yd_mt), maxabs(Ye_mt))
        max_dimless_uv = max(abs(gY_uv), abs(g2_uv), abs(g3_uv), abs(lam_uv), maxabs(Yu_uv), maxabs(Yd_uv), maxabs(Ye_uv))

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="rg_run_succeeds",
                passed=True,
                detail=f"mt→μUV run completed; λ(mt)={lam_mt:.6g}, λ(μUV)={lam_uv:.6g}",
            )
        )
        pub = rge.get("publication_grade", {}) or {}
        thr_ok = bool(pub.get("threshold_matching_ok", False))
        blocked = pub.get("blocked_thresholds", []) or []
        checks.append(
            Check(
                check_id="threshold_matching_publication_grade",
                passed=bool(thr_ok),
                detail=f"threshold_matching_ok={thr_ok}, blocked_thresholds={blocked}",
            )
        )
        checks.append(
            Check(
                check_id="lambda_crossing_reported_if_needed",
                passed=bool((lam_mt <= 0) or (lam_uv >= 0) or (crossing.crosses and crossing.mu_cross_GeV is not None)),
                detail=f"crosses={crossing.crosses}, mu_cross≈{crossing.mu_cross_GeV} GeV ({crossing.note})",
            )
        )
        checks.append(
            Check(
                check_id="couplings_finite",
                passed=bool(np.isfinite(max_dimless_mt) and np.isfinite(max_dimless_uv)),
                detail=f"max|dimless|(mt)={max_dimless_mt:.6g}, max|dimless|(μUV)={max_dimless_uv:.6g}",
            )
        )
        checks.append(
            Check(
                check_id="perturbative_window",
                passed=bool(max_dimless_mt < ceiling and max_dimless_uv < ceiling),
                detail=f"ceiling=4π≈{ceiling:.6g}; max|dimless|(mt)={max_dimless_mt:.6g}, max|dimless|(μUV)={max_dimless_uv:.6g}",
            )
        )

        lines: list[str] = []
        lines += [
            "Stability & perturbativity audit (diagnostic; not a publication-grade vacuum lifetime analysis)",
            "",
            f"SM inputs: {_relpath(sm_path)} (sha256={_sha256_file(sm_path)})",
            f"thresholds: {_relpath(thr_path)} (sha256={_sha256_file(thr_path)})",
            "",
            "Boundary at mt (from sm_boundary_conditions_at_mt):",
            f"- mt={mu_start_GeV:.6g} GeV, g1_GUT(mt)={g1_gut_mt:.6g} ⇒ gY(mt)={gY_mt:.6g}, g2(mt)={g2_mt:.6g}, g3(mt)={g3_mt:.6g}",
            f"- yt(mt)={yt_mt:.6g}, yb(mt)={yb_mt:.6g}, yτ(mt)={ytau_mt:.6g}, λ(mt)={lam_mt:.6g}",
            "",
            "UV endpoint (mt→μUV, 2-loop PyR@TE engine with explicit thresholds):",
            f"- μUV={mu_uv_GeV:.3e} GeV, (gY,g2,g3)(μUV)=({gY_uv:.6g},{g2_uv:.6g},{g3_uv:.6g}), λ(μUV)={lam_uv:.6g}",
            "",
            "Vacuum stability proxy:",
            f"- λ(mt)={lam_mt:.6g}, λ(μUV)={lam_uv:.6g}",
            f"- λ=0 crossing: {crossing.mu_cross_GeV:.3e} GeV ({crossing.note})" if crossing.crosses and crossing.mu_cross_GeV is not None else "- λ=0 crossing: none bracketed",
            "",
            "Perturbativity proxy:",
            f"- ceiling=4π≈{ceiling:.6g}",
            f"- max|dimless|(mt)={max_dimless_mt:.6g}, max|dimless|(μUV)={max_dimless_uv:.6g}",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This is a red-flag detector (engineering). A publication-grade metastability statement requires a vacuum lifetime computation and careful matching/uncertainty propagation.",
        ]

        return ModuleResult(
            results={
                "inputs": {
                    "sm_inputs_file": _relpath(sm_path),
                    "thresholds_file": _relpath(thr_path),
                },
                "mt": {
                    "mu_GeV": mu_start_GeV,
                    "gY": gY_mt,
                    "g2": g2_mt,
                    "g3": g3_mt,
                    "lambda": lam_mt,
                    "yt": yt_mt,
                    "yb": yb_mt,
                    "ytau": ytau_mt,
                },
                "mu_uv": {
                    "mu_GeV": mu_uv_GeV,
                    "gY": gY_uv,
                    "g2": g2_uv,
                    "g3": g3_uv,
                    "lambda": lam_uv,
                    "Yu_maxabs": maxabs(Yu_uv),
                    "Yd_maxabs": maxabs(Yd_uv),
                    "Ye_maxabs": maxabs(Ye_uv),
                },
                "lambda_crossing": crossing.__dict__,
                "perturbativity": {"ceiling_4pi": ceiling, "max_dimless_mt": max_dimless_mt, "max_dimless_mu_uv": max_dimless_uv},
                "publication_grade": rge.get("publication_grade", None),
                "segments": rge.get("segments", []),
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

