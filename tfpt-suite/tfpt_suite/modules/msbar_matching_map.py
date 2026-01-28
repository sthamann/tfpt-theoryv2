from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.conventions import alpha_from_g, gY_from_g1_gut
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule

COV_EIGENVALUE_TOL = 1e-18
from tfpt_suite.pyrate_boundary_runner import sm_boundary_conditions_at_mt
from tfpt_suite.rge_sm import run_sm_gauge_only_2loop_thresholds
from tfpt_suite.sm_inputs import SmMzInputs, gauge_couplings_from_mz_inputs


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


def _mb_run_1loop_qcd(*, mb_mb_GeV: float, alpha_s_mu0: float, alpha_s_mu1: float, nf: int) -> float:
    """
    1-loop QCD running of a MSbar quark mass:
      m(μ1) = m(μ0) * (αs(μ1)/αs(μ0))^(12/(33-2 nf))

    This is used here as a deterministic *phase-1* matching primitive, explicitly documented
    and intended to be upgraded (higher loops + EW finite pieces) as the matching layer matures.
    """
    m0 = float(mb_mb_GeV)
    if m0 <= 0:
        raise ValueError("mb_mb_GeV must be positive")
    a0 = float(alpha_s_mu0)
    a1 = float(alpha_s_mu1)
    if a0 <= 0 or a1 <= 0:
        raise ValueError("alpha_s inputs must be positive for 1-loop mass running")
    nf_i = int(nf)
    if nf_i < 3:
        raise ValueError("nf must be >=3")
    expo = 12.0 / (33.0 - 2.0 * float(nf_i))
    return float(m0 * (a1 / a0) ** expo)


class MsbarMatchingMapModule(TfptModule):
    module_id = "msbar_matching_map"
    title = "MSbar matching map (PDG-style inputs → MSbar couplings/Yukawas; EFT cascade + uncertainty propagation)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "SM inputs at MZ: tfpt_suite/data/sm_inputs_mz.json",
                "lepton masses (for y_tau proxy): tfpt_suite/data/lepton_masses_pdg.json",
            ],
            outputs=[
                "gauge couplings at MZ and mt (2-loop gauge-only, with αs threshold matching at heavy-quark thresholds)",
                "yt(mt) from mt_pole→mt_MSbar(mt) (QCD 2-loop, as used by RG-dressed pipelines)",
                "yb(mt) derived from mb(mb) via LO QCD running (nf=5) (explicit assumption)",
                "yτ(mt) proxy from mτ/v (tree; explicit assumption)",
                "EFT cascade diagnostics: αs(μ) at selected IR scales (mc, mb, MZ, mt)",
                "uncertainty propagation: deterministic Monte Carlo over input priors (if σ fields are provided)",
            ],
            formulas=[
                "g2(MZ)=e/sinθW, gY(MZ)=e/cosθW, g1_GUT=√(5/3) gY",
                "αs = g3^2/(4π)",
                "yt(mt) = √2 · m_t^MSbar(mt) / v (m_t^MSbar from QCD 2-loop pole→MSbar)",
                "mb(mt) ≈ mb(mb) [αs(mt)/αs(mb)]^(12/(33-2nf)) (1-loop; nf=5)",
                "yb(mt)=√2 · mb(mt) / v",
            ],
            validation=[
                "all derived couplings are finite and positive where expected",
                "αs sensitivity changes g3(mt) and (yt,yb) (sanity)",
                "Monte Carlo (if enabled) produces finite mean/σ and is deterministic under SuiteConfig.seed",
            ],
            determinism="Deterministic given inputs (no fitter).",
        )

    def run(self, config) -> ModuleResult:
        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"
        lep_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        lep_raw = json.loads(lep_path.read_text(encoding="utf-8"))

        mu_mz = float(sm_raw["mu_GeV"])
        MZ = float(sm_raw.get("mu_GeV", 91.1876))
        mt_pole = float(sm_raw.get("mt_GeV", 172.76))
        mb_mb = float(sm_raw.get("mb_GeV", 4.18))  # interpreted as mb(mb) in GeV (PDG-ish)
        mc = float(sm_raw.get("mc_GeV", 1.27))
        mb_thr = float(sm_raw.get("mb_GeV", 4.18))
        mt_thr = float(sm_raw.get("mt_GeV", 172.76))
        v_ev = float(sm_raw.get("v_ev_GeV", 246.0))

        alpha_em_inv = float(sm_raw["alpha_em_inv"])
        sin2 = float(sm_raw["sin2_thetaW"])
        alpha_s_mz = float(sm_raw["alpha_s"])
        alpha_s_sigma = float(sm_raw.get("alpha_s_sigma", 0.0) or 0.0)
        alpha_em_inv_sigma = float(sm_raw.get("alpha_em_inv_sigma", 0.0) or 0.0)
        sin2_sigma = float(sm_raw.get("sin2_thetaW_sigma", 0.0) or 0.0)
        mt_sigma = float(sm_raw.get("mt_sigma_GeV", 0.0) or 0.0)
        mb_sigma = float(sm_raw.get("mb_sigma_GeV", 0.0) or 0.0)
        mc_sigma = float(sm_raw.get("mc_sigma_GeV", 0.0) or 0.0)
        mc_samples = int(sm_raw.get("matching_mc_samples", 0) or 0)

        inp = SmMzInputs(mu_GeV=mu_mz, alpha_em_inv=alpha_em_inv, sin2_thetaW=sin2, alpha_s=alpha_s_mz)
        g1_gut_mz, g2_mz, g3_mz = gauge_couplings_from_mz_inputs(inp)
        gY_mz = gY_from_g1_gut(g1_gut_mz)

        # Standard mt boundary used by the suite (includes αs matching and mt pole→MSbar)
        bc = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_raw)
        gY_mt = float(bc.route_2loop["gY_mt"])
        g2_mt = float(bc.route_2loop["g2_mt"])
        g3_mt = float(bc.route_2loop["g3_mt"])
        alpha_s_mt = float(bc.route_2loop["alpha_s_mt"])
        mt_msbar_mt = float(bc.route_2loop["mt_msbar_mt_GeV"])
        yt_mt = float(bc.route_2loop["yt_mt"])

        # For mb(mt) running, use nf=5 couplings. We take αs(mt) from the authoritative mt boundary
        # (MZ→mt via PyR@TE) and keep the below-MZ αs policy explicit.
        g_mb_nf5 = run_sm_gauge_only_2loop_thresholds(
            mu_start_GeV=MZ,
            mu_end_GeV=mb_mb,
            g_start=(g1_gut_mz, g2_mz, g3_mz),
            mc_GeV=mc,
            mb_GeV=mb_thr,
            mt_GeV=mt_thr,
            apply_alpha3_matching=False,
        )
        alpha_s_mt_nf5 = float(bc.route_2loop["alpha_s_mt"])
        alpha_s_mb_nf5 = float(alpha_from_g(float(g_mb_nf5[2])))

        mb_mt = _mb_run_1loop_qcd(mb_mb_GeV=mb_mb, alpha_s_mu0=alpha_s_mb_nf5, alpha_s_mu1=alpha_s_mt_nf5, nf=5)
        yb_mt = float(math.sqrt(2.0) * mb_mt / v_ev)

        # Tau Yukawa proxy (no QED running here yet)
        mtau = float(lep_raw["masses"]["tau"]["mean"])
        ytau_mt = float(math.sqrt(2.0) * mtau / v_ev)

        # EFT cascade diagnostics (QCD-focused): αs at selected IR scales (same running policy as rge_sm gauge-only thresholds)
        g_mc_nf4 = run_sm_gauge_only_2loop_thresholds(
            mu_start_GeV=MZ,
            mu_end_GeV=mc,
            g_start=(g1_gut_mz, g2_mz, g3_mz),
            mc_GeV=mc,
            mb_GeV=mb_thr,
            mt_GeV=mt_thr,
            apply_alpha3_matching=True,
        )
        alpha_s_mc = float(alpha_from_g(float(g_mc_nf4[2])))

        # αs(MZ) sensitivity (deterministic ±σ scan)
        sens = None
        if alpha_s_sigma and alpha_s_sigma > 0:
            rows: list[dict[str, float]] = []
            for tag, a_s in [
                ("minus_sigma", float(alpha_s_mz) - float(alpha_s_sigma)),
                ("central", float(alpha_s_mz)),
                ("plus_sigma", float(alpha_s_mz) + float(alpha_s_sigma)),
            ]:
                sm_sc = dict(sm_raw)
                sm_sc["alpha_s"] = float(a_s)
                bc_sc = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_sc)

                inp_sc = SmMzInputs(mu_GeV=mu_mz, alpha_em_inv=alpha_em_inv, sin2_thetaW=sin2, alpha_s=float(a_s))
                g1_gut_mz_sc, g2_mz_sc, g3_mz_sc = gauge_couplings_from_mz_inputs(inp_sc)
                g_mb_nf5_sc = run_sm_gauge_only_2loop_thresholds(
                    mu_start_GeV=MZ,
                    mu_end_GeV=mb_mb,
                    g_start=(g1_gut_mz_sc, g2_mz_sc, g3_mz_sc),
                    mc_GeV=mc,
                    mb_GeV=mb_thr,
                    mt_GeV=mt_thr,
                    apply_alpha3_matching=False,
                )
                a_mt_nf5_sc = float(bc_sc.route_2loop["alpha_s_mt"])
                a_mb_nf5_sc = float(alpha_from_g(float(g_mb_nf5_sc[2])))
                mb_mt_sc = _mb_run_1loop_qcd(mb_mb_GeV=mb_mb, alpha_s_mu0=a_mb_nf5_sc, alpha_s_mu1=a_mt_nf5_sc, nf=5)
                yb_mt_sc = float(math.sqrt(2.0) * mb_mt_sc / v_ev)
                rows.append(
                    {
                        "alpha_s_MZ": float(a_s),
                        "g3_MZ": float(g3_mz_sc),
                        "g3_mt": float(bc_sc.route_2loop["g3_mt"]),
                        "yt_mt": float(bc_sc.route_2loop["yt_mt"]),
                        "yb_mt_derived": float(yb_mt_sc),
                    }
                )
            # half-range as a simple symmetric sensitivity proxy
            yb_half = 0.5 * abs(rows[2]["yb_mt_derived"] - rows[0]["yb_mt_derived"])
            yt_half = 0.5 * abs(rows[2]["yt_mt"] - rows[0]["yt_mt"])
            sens = {"rows": rows, "half_range": {"yb_mt": float(yb_half), "yt_mt": float(yt_half)}}

        # Monte Carlo uncertainty propagation (deterministic, optional).
        #
        # This is intentionally conservative: if no σ fields are provided, we do not pretend to have an error budget.
        mc_prop = None
        if mc_samples and mc_samples > 0 and (
            alpha_s_sigma > 0 or alpha_em_inv_sigma > 0 or sin2_sigma > 0 or mt_sigma > 0 or mb_sigma > 0 or mc_sigma > 0
        ):
            rng = np.random.default_rng(int(config.seed))

            def _clip_pos(x: float, eps: float = 1e-12) -> float:
                return float(max(float(x), float(eps)))

            def _clip_sin2(x: float) -> float:
                return float(min(max(float(x), 1e-6), 1.0 - 1e-6))

            def _sample_gauss(mean: float, sigma: float) -> float:
                if sigma <= 0:
                    return float(mean)
                return float(rng.normal(loc=float(mean), scale=float(sigma)))

            rows_mc: list[dict[str, float]] = []
            for _ in range(int(mc_samples)):
                aem_inv_i = _clip_pos(_sample_gauss(alpha_em_inv, alpha_em_inv_sigma))
                sin2_i = _clip_sin2(_sample_gauss(sin2, sin2_sigma))
                as_i = _clip_pos(_sample_gauss(alpha_s_mz, alpha_s_sigma))
                mt_i = _clip_pos(_sample_gauss(mt_pole, mt_sigma))
                mb_i = _clip_pos(_sample_gauss(mb_mb, mb_sigma))
                mc_i = _clip_pos(_sample_gauss(mc, mc_sigma))

                sm_sc = dict(sm_raw)
                sm_sc["alpha_em_inv"] = float(aem_inv_i)
                sm_sc["sin2_thetaW"] = float(sin2_i)
                sm_sc["alpha_s"] = float(as_i)
                sm_sc["mt_GeV"] = float(mt_i)
                sm_sc["mb_GeV"] = float(mb_i)
                sm_sc["mc_GeV"] = float(mc_i)

                bc_sc = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_sc)
                rows_mc.append(
                    {
                        "alpha_s_MZ": float(as_i),
                        "alpha_em_inv_MZ": float(aem_inv_i),
                        "sin2_thetaW_MZ": float(sin2_i),
                        "mt_pole_GeV": float(mt_i),
                        "g3_mt": float(bc_sc.route_2loop["g3_mt"]),
                        "yt_mt": float(bc_sc.route_2loop["yt_mt"]),
                        "yb_mt": float(bc_sc.route_2loop["yb_mt"]),
                    }
                )

            def _mean_std(key: str) -> tuple[float, float]:
                arr = np.array([r[key] for r in rows_mc], dtype=float)
                if arr.size < 2:
                    return float(np.mean(arr)), 0.0
                return float(np.mean(arr)), float(np.std(arr, ddof=1))

            g3_mean, g3_std = _mean_std("g3_mt")
            yt_mean, yt_std = _mean_std("yt_mt")
            yb_mean, yb_std = _mean_std("yb_mt")
            mc_prop = {
                "samples": int(mc_samples),
                "seed": int(config.seed),
                "inputs_sigma": {
                    "alpha_s_sigma": float(alpha_s_sigma),
                    "alpha_em_inv_sigma": float(alpha_em_inv_sigma),
                    "sin2_thetaW_sigma": float(sin2_sigma),
                    "mt_sigma_GeV": float(mt_sigma),
                    "mb_sigma_GeV": float(mb_sigma),
                    "mc_sigma_GeV": float(mc_sigma),
                },
                "summary": {
                    "g3_mt": {"mean": g3_mean, "std": g3_std},
                    "yt_mt": {"mean": yt_mean, "std": yt_std},
                    "yb_mt": {"mean": yb_mean, "std": yb_std},
                },
                "rows_truncated": rows_mc[: min(len(rows_mc), 50)],
                "note": "rows are truncated to 50 samples to keep results.json and PDF size manageable",
            }

        # Linear covariance propagation (Jacobian; deterministic).
        lin_cov: dict[str, Any] | None = None
        try:
            in_labels = ["alpha_s_MZ", "alpha_em_inv_MZ", "sin2_thetaW_MZ", "mt_pole_GeV", "mb_GeV", "mc_GeV"]
            x0 = np.array([alpha_s_mz, alpha_em_inv, sin2, mt_pole, mb_mb, mc], dtype=float)
            sx = np.array([alpha_s_sigma, alpha_em_inv_sigma, sin2_sigma, mt_sigma, mb_sigma, mc_sigma], dtype=float)

            def _compute_outputs(sm_inputs: dict[str, Any]) -> dict[str, float]:
                bc_loc = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_inputs)
                return {
                    "gY_mt": float(bc_loc.route_2loop["gY_mt"]),
                    "g2_mt": float(bc_loc.route_2loop["g2_mt"]),
                    "g3_mt": float(bc_loc.route_2loop["g3_mt"]),
                    "yt_mt": float(bc_loc.route_2loop["yt_mt"]),
                    "yb_mt": float(bc_loc.route_2loop.get("yb_mt", yb_mt)),
                }

            out_labels = ["gY_mt", "g2_mt", "g3_mt", "yt_mt", "yb_mt"]
            y0 = _compute_outputs(dict(sm_raw))
            y0_vec = np.array([float(y0[k]) for k in out_labels], dtype=float)

            # Jacobian via symmetric finite differences. Use dp=sigma when available, otherwise a tiny relative dp.
            J = np.zeros((len(out_labels), len(in_labels)), dtype=float)
            for j, (lab, xj, sj) in enumerate(zip(in_labels, x0.tolist(), sx.tolist())):
                dp = float(sj) if (np.isfinite(sj) and sj > 0) else float(max(1e-8, 1e-6 * abs(float(xj))))
                if not (np.isfinite(dp) and dp > 0):
                    continue
                sm_p = dict(sm_raw)
                sm_m = dict(sm_raw)
                if lab == "alpha_s_MZ":
                    sm_p["alpha_s"] = float(xj + dp)
                    sm_m["alpha_s"] = float(xj - dp)
                elif lab == "alpha_em_inv_MZ":
                    sm_p["alpha_em_inv"] = float(xj + dp)
                    sm_m["alpha_em_inv"] = float(xj - dp)
                elif lab == "sin2_thetaW_MZ":
                    sm_p["sin2_thetaW"] = float(xj + dp)
                    sm_m["sin2_thetaW"] = float(xj - dp)
                elif lab == "mt_pole_GeV":
                    sm_p["mt_GeV"] = float(xj + dp)
                    sm_m["mt_GeV"] = float(xj - dp)
                elif lab == "mb_GeV":
                    sm_p["mb_GeV"] = float(xj + dp)
                    sm_m["mb_GeV"] = float(xj - dp)
                elif lab == "mc_GeV":
                    sm_p["mc_GeV"] = float(xj + dp)
                    sm_m["mc_GeV"] = float(xj - dp)
                else:
                    continue
                yp = _compute_outputs(sm_p)
                ym = _compute_outputs(sm_m)
                yp_vec = np.array([float(yp[k]) for k in out_labels], dtype=float)
                ym_vec = np.array([float(ym[k]) for k in out_labels], dtype=float)
                J[:, j] = (yp_vec - ym_vec) / (2.0 * dp)

            # Input covariance (diagonal for now; upgrade to full covariance when supplied).
            Cov_x = np.diag(sx * sx)
            Cov_y = J @ Cov_x @ J.T
            sig_y = np.sqrt(np.clip(np.diag(Cov_y), 0.0, np.inf))
            Corr_y = np.zeros_like(Cov_y)
            for i in range(Cov_y.shape[0]):
                for k in range(Cov_y.shape[1]):
                    denom = float(sig_y[i] * sig_y[k])
                    Corr_y[i, k] = float(Cov_y[i, k] / denom) if denom > 0 else float("nan")

            lin_cov = {
                "inputs": {
                    "labels": in_labels,
                    "means": [float(x) for x in x0.tolist()],
                    "sigmas": [float(s) for s in sx.tolist()],
                    "note": "assumed diagonal input covariance from sm_inputs_mz.json σ fields",
                },
                "outputs": {"labels": out_labels, "means": [float(x) for x in y0_vec.tolist()]},
                "jacobian": [[float(x) for x in row] for row in J.tolist()],
                "covariance": [[float(x) for x in row] for row in Cov_y.tolist()],
                "correlation": [[float(x) for x in row] for row in Corr_y.tolist()],
            }
        except Exception as e:
            lin_cov = {"error": str(e)}

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="inputs_loaded",
                passed=bool(sm_path.exists() and lep_path.exists()),
                detail=f"sm={_relpath(sm_path)} lep={_relpath(lep_path)}",
            )
        )
        checks.append(
            Check(
                check_id="gauge_couplings_finite",
                passed=bool(all(np.isfinite(x) and x > 0 for x in [g1_gut_mz, g2_mz, g3_mz, gY_mz, gY_mt, g2_mt, g3_mt])),
                detail=f"(gY,g2,g3)(MZ)=({gY_mz:.6g},{g2_mz:.6g},{g3_mz:.6g}), (gY,g2,g3)(mt)=({gY_mt:.6g},{g2_mt:.6g},{g3_mt:.6g})",
            )
        )
        checks.append(
            Check(
                check_id="top_msbar_positive",
                passed=bool(np.isfinite(mt_msbar_mt) and mt_msbar_mt > 0 and np.isfinite(yt_mt) and yt_mt > 0),
                detail=f"mt_MSbar(mt)={mt_msbar_mt:.6g} GeV, yt(mt)={yt_mt:.6g} (from sm_boundary_conditions_at_mt)",
            )
        )
        checks.append(
            Check(
                check_id="bottom_running_finite",
                passed=bool(np.isfinite(mb_mt) and mb_mt > 0 and np.isfinite(yb_mt) and yb_mt > 0),
                detail=f"mb(mb)={mb_mb:.6g} GeV → mb(mt)≈{mb_mt:.6g} GeV (1-loop nf=5), yb(mt)≈{yb_mt:.6g}",
            )
        )
        if sens is not None:
            checks.append(
                Check(
                    check_id="alpha_s_sensitivity_nonzero",
                    passed=bool(sens["half_range"]["yb_mt"] > 0 and sens["half_range"]["yt_mt"] > 0),
                    detail=f"half-range: Δyb_mt≈{sens['half_range']['yb_mt']:.3e}, Δyt_mt≈{sens['half_range']['yt_mt']:.3e} for αs(MZ)±σ",
                )
            )
        if mc_prop is not None:
            checks.append(
                Check(
                    check_id="matching_mc_present",
                    passed=bool(int(mc_prop.get("samples", 0)) >= 10),
                    detail=f"samples={mc_prop.get('samples')}, σ(yt_mt)≈{mc_prop['summary']['yt_mt']['std']:.3e}, σ(yb_mt)≈{mc_prop['summary']['yb_mt']['std']:.3e}",
                )
            )
        if isinstance(lin_cov, dict) and "covariance" in lin_cov:
            try:
                cov = np.array(lin_cov["covariance"], dtype=float)
                shape_ok = bool(np.all(np.isfinite(cov)) and cov.shape[0] == cov.shape[1] and cov.shape[0] >= 2)
                if shape_ok:
                    eigvals = np.linalg.eigvalsh(cov)
                    min_eig = float(np.min(eigvals))
                    ok_cov = bool(min_eig >= -COV_EIGENVALUE_TOL)
                else:
                    min_eig = float("nan")
                    ok_cov = False
            except Exception:
                min_eig = float("nan")
                ok_cov = False
            checks.append(
                Check(
                    check_id="covariance_propagated_end_to_end",
                    passed=ok_cov,
                    detail=(
                        f"linear Jacobian covariance computed (min_eig={min_eig:.3e})"
                        if ok_cov
                        else "linear Jacobian covariance computation failed"
                    ),
                )
            )

        lines: list[str] = []
        lines += [
            "MSbar matching map (explicit, deterministic; EFT-cascade diagnostics + uncertainty propagation hooks)",
            "",
            f"SM inputs: {_relpath(sm_path)} (sha256={_sha256_file(sm_path)})",
            f"lepton masses: {_relpath(lep_path)} (sha256={_sha256_file(lep_path)})",
            "",
            "Gauge couplings at MZ (from α_em(MZ), sin²θW(MZ), αs(MZ)):",
            f"- gY(MZ)={gY_mz:.8g}, g2(MZ)={g2_mz:.8g}, g3(MZ)={g3_mz:.8g}",
            f"- αs(MZ) input = {alpha_s_mz:.6g} ± {alpha_s_sigma:.3g}",
            "",
            "mt boundary (suite policy; includes αs threshold matching + mt pole→MSbar(mt) QCD 2-loop):",
            f"- gY(mt)={gY_mt:.8g}, g2(mt)={g2_mt:.8g}, g3(mt)={g3_mt:.8g}  (αs(mt)={alpha_s_mt:.6g})",
            f"- mt_pole={mt_pole:.6g} GeV → mt_MSbar(mt)≈{mt_msbar_mt:.6g} GeV ⇒ yt(mt)≈{yt_mt:.6g} (v={v_ev:.6g} GeV)",
            "",
            "Derived Yukawas (phase-1 placeholders for full matching):",
            f"- mb(mt)≈{mb_mt:.6g} GeV from mb(mb)={mb_mb:.6g} GeV via 1-loop QCD (nf=5) ⇒ yb(mt)≈{yb_mt:.6g}",
            f"- yτ(mt)≈{ytau_mt:.6g} from mτ={mtau:.6g} GeV (no QED running yet)",
            "",
            "EFT cascade (QCD-focused diagnostics; gauge-only running policy, no EW decoupling below MZ):",
            f"- αs(mc={mc:.6g} GeV) ≈ {alpha_s_mc:.6g} (incl. 2-loop finite αs matching at mb threshold)",
            f"- αs(mb={mb_mb:.6g} GeV) ≈ {alpha_s_mb_nf5:.6g}",
            f"- αs(MZ={MZ:.6g} GeV) = {alpha_s_mz:.6g} (input)",
            f"- αs(mt={mt_pole:.6g} GeV) ≈ {alpha_s_mt_nf5:.6g} (nf=5, no top-threshold matching)",
            "",
        ]
        if sens is not None:
            lines += [
                "αs(MZ) sensitivity (±σ):",
                *[
                    f"- {r['alpha_s_MZ']:.6g}: g3(MZ)={r['g3_MZ']:.6g}, g3(mt)={r['g3_mt']:.6g}, yt(mt)={r['yt_mt']:.6g}, yb(mt)={r['yb_mt_derived']:.6g}"
                    for r in sens["rows"]
                ],
                f"- half-range: Δyt≈{sens['half_range']['yt_mt']:.3e}, Δyb≈{sens['half_range']['yb_mt']:.3e}",
                "",
            ]
        if mc_prop is not None:
            lines += [
                "Uncertainty propagation (Monte Carlo; deterministic):",
                f"- samples = {mc_prop['samples']} (seed={mc_prop['seed']})",
                f"- yt(mt): mean={mc_prop['summary']['yt_mt']['mean']:.6g}, std={mc_prop['summary']['yt_mt']['std']:.3e}",
                f"- yb(mt): mean={mc_prop['summary']['yb_mt']['mean']:.6g}, std={mc_prop['summary']['yb_mt']['std']:.3e}",
                "",
            ]
        if isinstance(lin_cov, dict) and "covariance" in lin_cov and isinstance(lin_cov.get("outputs", {}), dict):
            try:
                diag = np.diag(np.array(lin_cov["covariance"], dtype=float))
                sigs = np.sqrt(np.clip(diag, 0.0, np.inf))
                out_labs = list(lin_cov["outputs"]["labels"])
                lines += [
                    "Uncertainty propagation (linear Jacobian):",
                    "- " + ", ".join([f"{lab}±{sig:.3e}" for lab, sig in zip(out_labs, sigs.tolist())]),
                    "",
                ]
            except Exception:
                pass

        lines += [
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module is a bookkeeping bridge: it makes the PDG-style input policy explicit and reproducible.",
            "- EFT cascade here is QCD-focused and does not implement a full electroweak decoupling policy below MZ (W/Z/H thresholds).",
            "- The bottom and tau Yukawas here use explicit approximations (LO QCD for mb, tree-level yτ). Upgrade to higher-loop running + finite EW/QCD pieces for publication-grade matching.",
        ]

        return ModuleResult(
            results={
                "inputs": {
                    "sm_inputs_file": _relpath(sm_path),
                    "sm_inputs_sha256": _sha256_file(sm_path),
                    "lepton_masses_file": _relpath(lep_path),
                    "lepton_masses_sha256": _sha256_file(lep_path),
                },
                "MZ": {
                    "mu_GeV": mu_mz,
                    "gY": gY_mz,
                    "g2": g2_mz,
                    "g3": g3_mz,
                    "alpha_s": alpha_s_mz,
                },
                "mt_boundary": {
                    "mu_GeV": float(bc.mu_mt_GeV),
                    "gY": gY_mt,
                    "g2": g2_mt,
                    "g3": g3_mt,
                    "alpha_s_mt": alpha_s_mt,
                    "mt_msbar_mt_GeV": mt_msbar_mt,
                    "yt_mt": yt_mt,
                },
                "derived_yukawas_mt": {
                    "v_ev_GeV": v_ev,
                    "mb_mt_GeV_1loop_qcd_nf5": mb_mt,
                    "yb_mt_1loop_qcd_nf5": yb_mt,
                    "ytau_mt_tree": ytau_mt,
                },
                "eft_cascade_qcd_diagnostics": {
                    "alpha_s_mc_GeV": mc,
                    "alpha_s_mc": alpha_s_mc,
                    "alpha_s_mb_GeV": mb_mb,
                    "alpha_s_mb_nf5": alpha_s_mb_nf5,
                    "alpha_s_mz_GeV": MZ,
                    "alpha_s_mz": alpha_s_mz,
                    "alpha_s_mt_GeV": mt_pole,
                    "alpha_s_mt_nf5": alpha_s_mt_nf5,
                    "finite_alpha3_matching_at_mb_threshold": True,
                },
                "alpha_s_sensitivity": sens,
                "matching_mc": mc_prop,
                "linear_covariance": lin_cov,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

