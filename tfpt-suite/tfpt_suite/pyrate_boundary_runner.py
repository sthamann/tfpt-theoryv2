from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.conventions import alpha_from_g, g1_gut_over_gY, gY_from_g1_gut
from tfpt_suite.pyrate_pythonoutputs import get_pyrate_pythonoutput
from tfpt_suite.rge_pyrate_2loop import load_pyrate_beta_module, run_flavor_rge_2loop_thresholds
from tfpt_suite.rge_sm import SmBoundaryAtMt, msbar_mass_run_qcd_1loop, pole_mass_to_msbar_mass_qcd_2loop, run_sm_gauge_only_2loop_thresholds
from tfpt_suite.sm_inputs import SmMzInputs, gauge_couplings_from_mz_inputs


def _workspace_root() -> Path:
    # .../wolfram_latex_attachments/tfpt-suite/tfpt_suite/<file>.py -> parents[2] is workspace root
    return Path(__file__).resolve().parents[2]


def _load_thresholds_v25() -> dict[str, float]:
    thr_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "rge_thresholds_v25.json"
    raw = json.loads(thr_path.read_text(encoding="utf-8"))
    thresholds = raw.get("thresholds_GeV", {}) if isinstance(raw, dict) else {}
    if not isinstance(thresholds, dict):
        return {}
    return {str(k): float(v) for k, v in thresholds.items()}


def _default_yukawa_seeds_mz() -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Deterministic MZ seeds for a short SM PyR@TE run.

    These are *not* a fit; they are a stable seed consistent with the existing
    `two_loop_rg_fingerprints` module, and the MZ→mt interval is small.
    """
    Yu = np.diag([0.0, 0.0, 0.94]).astype(complex)
    Yd = np.diag([0.0, 0.0, 0.017]).astype(complex)
    Ye = np.diag([0.0, 0.0, 0.010]).astype(complex)
    lam = 0.13
    return Yu, Yd, Ye, float(lam)


def sm_boundary_conditions_at_mt(*, sm_inputs_mz: dict[str, Any]) -> SmBoundaryAtMt:
    """
    Construct SM MSbar-ish boundary conditions at μ = m_t, using a *PyR@TE-driven* short RG run MZ→mt.

    This is the authoritative mt boundary for mt→UV modules:
    - gauge couplings (gY,g2,g3) are obtained by integrating PyR@TE beta functions from μ=MZ to μ=mt
      (SM model `sm_tfpt_2loop_v25`, loop-truncated as requested).
    - Yukawas and λ at mt are *matching inputs* (mass-derived, not RG solved here), kept explicit.

    We provide two routes:
    - route_2loop: PyR@TE betas truncated at 2 loops (default)
    - route_1loop: same betas truncated at 1 loop (diagnostic)
    """
    mu_mz = float(sm_inputs_mz["mu_GeV"])
    mt_pole = float(sm_inputs_mz.get("mt_GeV", 172.76))
    mc = float(sm_inputs_mz.get("mc_GeV", 1.27))
    mb = float(sm_inputs_mz.get("mb_GeV", 4.18))
    mH = float(sm_inputs_mz.get("mH_GeV", 125.25))

    alpha_em_inv = float(sm_inputs_mz["alpha_em_inv"])
    sin2 = float(sm_inputs_mz["sin2_thetaW"])
    alpha_s_mz = float(sm_inputs_mz["alpha_s"])

    v_ev = float(sm_inputs_mz.get("v_ev_GeV", 246.0))
    yb_mt_in = sm_inputs_mz.get("yb_mt", None)
    ytau_mt_in = sm_inputs_mz.get("ytau_mt", None)

    if "lambda_mt" in sm_inputs_mz:
        lam_mt = float(sm_inputs_mz["lambda_mt"])
        lam_source = "input"
    else:
        lam_mt = float((mH * mH) / (2.0 * v_ev * v_ev))
        lam_source = "tree_from_mH_and_v"

    # MZ gauge couplings (algebraic conversion, not running)
    inp = SmMzInputs(mu_GeV=mu_mz, alpha_em_inv=alpha_em_inv, sin2_thetaW=sin2, alpha_s=alpha_s_mz)
    g1_gut_mz, g2_mz, g3_mz = gauge_couplings_from_mz_inputs(inp)
    gY_mz = float(gY_from_g1_gut(g1_gut_mz))

    # PyR@TE beta modules (fail-fast via expected model names + SHA256 in the loader)
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
    thresholds_GeV = _load_thresholds_v25()

    Yu_seed, Yd_seed, Ye_seed, lam_seed = _default_yukawa_seeds_mz()

    def _run_pyrate_route(max_loop: int) -> dict[str, float]:
        out = run_flavor_rge_2loop_thresholds(
            mu_start_GeV=mu_mz,
            mu_end_GeV=mt_pole,
            thresholds_GeV=thresholds_GeV,
            g_start=(gY_mz, float(g2_mz), float(g3_mz)),
            Yu_start=Yu_seed,
            Yd_start=Yd_seed,
            Ye_start=Ye_seed,
            lambda_start=float(lam_seed),
            yN_start=None,
            beta_sm=beta_sm,
            beta_e8=beta_e8,
            apply_sigma_threshold=True,
            apply_g8_delta_b3=True,
            delta_b3_g8=2.0,
            apply_matching=False,
            apply_gravity_alpha3=False,
            max_loop=int(max_loop),
            rtol=1e-10,
            atol=1e-12,
            method="DOP853",
        )
        gY_mt = float(out["g_end_sm"]["gY"])
        g2_mt = float(out["g_end_sm"]["g2"])
        g3_mt = float(out["g_end_sm"]["g3"])
        return {"gY_mt": gY_mt, "g2_mt": g2_mt, "g3_mt": g3_mt}

    # Route A: PyR@TE-truncated 2-loop
    gA = _run_pyrate_route(2)
    gY_mt_2l = float(gA["gY_mt"])
    g2_mt_2l = float(gA["g2_mt"])
    g3_mt_2l = float(gA["g3_mt"])
    g1_gut_mt_2l = float(g1_gut_over_gY() * gY_mt_2l)
    alpha_s_mt_2l = float(alpha_from_g(g3_mt_2l))
    mt_msbar_mt_2l = float(pole_mass_to_msbar_mass_qcd_2loop(pole_mass_GeV=mt_pole, alpha_s_mu=alpha_s_mt_2l, n_light=5))
    yt_mt_2l = float(math.sqrt(2.0) * mt_msbar_mt_2l / v_ev)

    # Route B: PyR@TE-truncated 1-loop (diagnostic)
    gB = _run_pyrate_route(1)
    gY_mt_1l = float(gB["gY_mt"])
    g2_mt_1l = float(gB["g2_mt"])
    g3_mt_1l = float(gB["g3_mt"])
    g1_gut_mt_1l = float(g1_gut_over_gY() * gY_mt_1l)
    alpha_s_mt_1l = float(alpha_from_g(g3_mt_1l))
    mt_msbar_mt_1l = float(pole_mass_to_msbar_mass_qcd_2loop(pole_mass_GeV=mt_pole, alpha_s_mu=alpha_s_mt_1l, n_light=5))
    yt_mt_1l = float(math.sqrt(2.0) * mt_msbar_mt_1l / v_ev)

    # yb(mt), ytau(mt): explicit inputs or deterministic derivations.
    #
    # Important: keep below-MZ mass-running as an explicit EFT assumption; do NOT replace the mt→UV RG master.
    if yb_mt_in is not None:
        yb_mt = float(yb_mt_in)
        yb_source = "input"
    else:
        # αs(mb) from a below-MZ gauge-only policy (QCD-focused).
        g_mb = run_sm_gauge_only_2loop_thresholds(
            mu_start_GeV=mu_mz,
            mu_end_GeV=mb,
            g_start=(float(g1_gut_mz), float(g2_mz), float(g3_mz)),
            mc_GeV=mc,
            mb_GeV=mb,
            mt_GeV=mt_pole,
            apply_alpha3_matching=False,
        )
        alpha_s_mb = float(alpha_from_g(float(g_mb[2])))
        mb_mt = float(msbar_mass_run_qcd_1loop(m_mu0_GeV=mb, alpha_s_mu0=alpha_s_mb, alpha_s_mu1=alpha_s_mt_2l, nf=5))
        yb_mt = float(math.sqrt(2.0) * mb_mt / v_ev)
        yb_source = "derived_from_mb_mb_via_1loop_qcd_nf5_with_alpha_s_mb_below_MZ"

    if ytau_mt_in is not None:
        ytau_mt = float(ytau_mt_in)
        ytau_source = "input"
    else:
        # Deterministic proxy: yτ ≈ √2 mτ/v (no QED running here).
        try:
            lep_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "lepton_masses_pdg.json"
            lep = json.loads(lep_path.read_text(encoding="utf-8")) if lep_path.exists() else {}
            mtau = float(lep.get("masses", {}).get("tau", {}).get("mean", 1.77686))
        except Exception:
            mtau = 1.77686
        ytau_mt = float(math.sqrt(2.0) * mtau / v_ev)
        ytau_source = "derived_from_mtau_over_v_tree"

    route_2loop: dict[str, Any] = {
        "mu_mz_GeV": mu_mz,
        "mu_mt_GeV": mt_pole,
        "g1_gut_mt": g1_gut_mt_2l,
        "gY_mt": gY_mt_2l,
        "g2_mt": g2_mt_2l,
        "g3_mt": g3_mt_2l,
        "alpha_s_mt": alpha_s_mt_2l,
        "mt_msbar_mt_GeV": mt_msbar_mt_2l,
        "yt_mt": yt_mt_2l,
        "yb_mt": yb_mt,
        "yb_source": yb_source,
        "ytau_mt": ytau_mt,
        "ytau_source": ytau_source,
        "lambda_mt": lam_mt,
        "v_ev_GeV": v_ev,
        "lambda_source": lam_source,
        "rge_source": "pyrate_beta_integrator",
        "max_loop": 2,
    }
    route_1loop: dict[str, Any] = {
        "mu_mz_GeV": mu_mz,
        "mu_mt_GeV": mt_pole,
        "g1_gut_mt": g1_gut_mt_1l,
        "gY_mt": gY_mt_1l,
        "g2_mt": g2_mt_1l,
        "g3_mt": g3_mt_1l,
        "alpha_s_mt": alpha_s_mt_1l,
        "mt_msbar_mt_GeV": mt_msbar_mt_1l,
        "yt_mt": yt_mt_1l,
        "yb_mt": yb_mt,
        "yb_source": yb_source,
        "ytau_mt": ytau_mt,
        "ytau_source": ytau_source,
        "lambda_mt": lam_mt,
        "v_ev_GeV": v_ev,
        "lambda_source": lam_source,
        "rge_source": "pyrate_beta_integrator",
        "max_loop": 1,
    }

    diffs = {
        "gY_mt": float(route_2loop["gY_mt"] - route_1loop["gY_mt"]),
        "g2_mt": float(route_2loop["g2_mt"] - route_1loop["g2_mt"]),
        "g3_mt": float(route_2loop["g3_mt"] - route_1loop["g3_mt"]),
        "yt_mt": float(route_2loop["yt_mt"] - route_1loop["yt_mt"]),
    }

    return SmBoundaryAtMt(
        mu_mt_GeV=float(mt_pole),
        scheme="MSbar-ish inputs @ MZ; MZ→mt gauge running via PyR@TE beta functions (loop-truncated); mt pole→MSbar(mt) via QCD 2-loop",
        conventions={
            "hypercharge_definition": "Q = T3 + Y",
            "g1_gut_over_gY": f"{g1_gut_over_gY():.15g}",
            "g1_gut_definition": "g1_GUT = sqrt(5/3) * gY",
            "gY_definition": "gY ≡ g′ (SM hypercharge coupling)",
            "lambda_source": lam_source,
            "mz_yukawa_seeds": "Yu33=0.94, Yd33=0.017, Ye33=0.010, lambda=0.13 (deterministic seed for short MZ→mt run)",
        },
        route_2loop=route_2loop,
        route_1loop=route_1loop,
        diffs=diffs,
    )

