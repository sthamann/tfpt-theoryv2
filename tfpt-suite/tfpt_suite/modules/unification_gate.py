from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.conventions import g1_gut_over_gY, gY_from_g1_gut
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass, mk_check_warn
from tfpt_suite.pyrate_boundary_runner import sm_boundary_conditions_at_mt
from tfpt_suite.pyrate_pythonoutputs import get_pyrate_pythonoutput
from tfpt_suite.rge_pyrate_2loop import load_pyrate_beta_module, run_flavor_rge_2loop_thresholds
from tfpt_suite.sm_inputs import SmMzInputs, gauge_couplings_from_mz_inputs


@dataclass(frozen=True)
class UnificationPoint:
    log10_mu: float
    mu_GeV: float
    alphaY: float
    alpha1_gut: float
    alpha2: float
    alpha3: float
    mismatch_rel: float


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


def _load_two_loop_cfg() -> dict[str, object]:
    cfg_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "two_loop_rg_fingerprints.json"
    raw = json.loads(cfg_path.read_text(encoding="utf-8"))
    raw["_cfg_path"] = str(cfg_path)
    return raw


def _load_unification_policy() -> dict[str, object]:
    cfg_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "unification_gate_policy.json"
    if cfg_path.exists():
        raw = json.loads(cfg_path.read_text(encoding="utf-8"))
        raw["_cfg_loaded"] = True
    else:
        raw = {"_cfg_loaded": False}
    raw["_cfg_path"] = str(cfg_path)
    return raw


def _alpha_from_g(g: np.ndarray) -> np.ndarray:
    return (np.array(g, dtype=float) ** 2) / (4.0 * np.pi)


def _best_unification_point(*, log10_mu: np.ndarray, alpha1_gut: np.ndarray, alpha2: np.ndarray, alpha3: np.ndarray) -> UnificationPoint:
    a1 = np.array(alpha1_gut, dtype=float)
    a2 = np.array(alpha2, dtype=float)
    a3 = np.array(alpha3, dtype=float)
    t = np.array(log10_mu, dtype=float)

    # mismatch in absolute α units
    f = np.maximum.reduce([np.abs(a1 - a2), np.abs(a1 - a3), np.abs(a2 - a3)])
    abar = (a1 + a2 + a3) / 3.0
    f_rel = f / np.maximum(abar, 1e-30)

    i = int(np.argmin(f_rel))
    t0 = float(t[i])

    # Quadratic refinement around the grid minimum (if possible).
    t_ref = t0
    if 0 < i < (t.size - 1):
        xs = np.array([t[i - 1], t[i], t[i + 1]], dtype=float)
        ys = np.array([f_rel[i - 1], f_rel[i], f_rel[i + 1]], dtype=float)
        try:
            a, b, c = np.polyfit(xs, ys, 2)
            if a != 0:
                t_ref = float(-b / (2.0 * a))
                t_ref = float(min(max(t_ref, xs.min()), xs.max()))
        except Exception:
            t_ref = t0

    # Interpolate α values at refined t
    a1_ref = float(np.interp(t_ref, t, a1))
    a2_ref = float(np.interp(t_ref, t, a2))
    a3_ref = float(np.interp(t_ref, t, a3))
    f_ref = float(max(abs(a1_ref - a2_ref), abs(a1_ref - a3_ref), abs(a2_ref - a3_ref)))
    abar_ref = float((a1_ref + a2_ref + a3_ref) / 3.0)
    f_rel_ref = float(f_ref / max(abar_ref, 1e-30))

    return UnificationPoint(
        log10_mu=t_ref,
        mu_GeV=float(10.0**t_ref),
        alphaY=float((3.0 / 5.0) * a1_ref),
        alpha1_gut=a1_ref,
        alpha2=a2_ref,
        alpha3=a3_ref,
        mismatch_rel=f_rel_ref,
    )


class UnificationGateModule(TfptModule):
    module_id = "unification_gate"
    title = "Unification gate (explicit 2-loop PyR@TE RG check)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "SM inputs at MZ: tfpt_suite/data/sm_inputs_mz.json (for initial gauge couplings)",
                "unification policy: tfpt_suite/data/unification_gate_policy.json (mt boundary + G8 patch policy)",
                "threshold policy: tfpt_suite/data/rge_thresholds_v25.json (segment boundaries)",
                "beta sources: tfpt_suite/data/pyrate_pythonoutputs.json (SM + E8)",
                "optional runner patch: tfpt_suite/data/two_loop_rg_fingerprints.json (gravity α^3 patch config)",
            ],
            outputs=[
                "best-fit unification scale μ* (minimizes max pairwise |Δα|)",
                "relative mismatch at μ* (gate metric)",
                "PASS/FAIL gate for unification within a declared tolerance",
            ],
            formulas=[
                "α_i = g_i^2/(4π); α1_GUT=(5/3)αY",
                "mismatch(μ) = max(|α1-α2|, |α1-α3|, |α2-α3|) / ((α1+α2+α3)/3)",
            ],
            validation=[
                "PyR@TE RGE solve succeeds and produces finite α arrays",
                "gate is FAIL unless mismatch_rel is below the declared tolerance",
            ],
            determinism="Deterministic given config + PyR@TE model.",
            question="Do the gauge couplings unify at a single scale (2-loop, declared model, declared tolerance)?",
            objective=[
                "Replace vague 'near unification' statements by an explicit numeric gate.",
                "Force loop-order consistency: use a single PyR@TE 2-loop run (no parallel gauge-only runner).",
            ],
            what_was_done=[
                "Run the TFPT segment runner (SM→E8 at MSigma; Δb3 patch above MG8) with matching enabled (1-loop identity at thresholds).",
                "Sample the gauge-coupling trajectory on a log grid (plus exact segment boundaries).",
                "Compute the best-fit unification point μ* by minimizing the maximum pairwise α mismatch.",
                "Emit a PASS/FAIL gate based on a fixed mismatch tolerance.",
                "Apply mt-boundary and optional G8 2-loop patch policies if enabled.",
            ],
            assumptions=[
                "This gate is an EFT+threshold diagnostic: it follows the declarative TFPT segment policy (SM↔E8 switch, Δb3 patch) rather than a single monolithic model.",
                "Finite matching pieces beyond 1-loop identity are not assumed unless explicitly provided in matching_finite_delta_alpha.",
                "Yukawa/quartic initial values are deterministic seeds (same spirit as rg_fingerprints) or mt-boundary outputs (if enabled) to keep the gate reproducible.",
                "Optional 2-loop Δb3 patch assumes a Weyl adjoint G8 unless explicit Δb3(2-loop) is provided.",
            ],
            gaps=[
                "Publication-grade unification would require full scheme-consistent finite matching (EW/QED) and a derivation of the threshold policy from the microscopic TFPT ladder (discrete threshold candidates).",
            ],
        )

    def run(self, config) -> ModuleResult:
        cfg = _load_two_loop_cfg()
        cfg_path = Path(str(cfg.get("_cfg_path", "")))
        grav = cfg.get("gravity_alpha3_patch", {}) if isinstance(cfg.get("gravity_alpha3_patch", {}), dict) else {}
        grav_enabled = bool(grav.get("enabled", False))
        grav_kappa = grav.get("kappa_vector", [0.0, 0.0, 0.0])
        try:
            k1, k2, k3 = (float(grav_kappa[0]), float(grav_kappa[1]), float(grav_kappa[2]))
        except Exception:
            k1, k2, k3 = (0.0, 0.0, 0.0)
        alpha_def_u1 = str(grav.get("alpha_definition_for_U1", "alphaY"))

        policy = _load_unification_policy()
        policy_path = Path(str(policy.get("_cfg_path", "")))
        policy_loaded = bool(policy.get("_cfg_loaded", False))
        use_mt_boundary = bool(policy.get("use_mt_boundary", True))
        boundary_route = str(policy.get("boundary_route", "pyrate_2loop")).strip().lower()
        mu_end = float(policy.get("mu_end_GeV", 1.0e19))
        apply_g8_delta_b3_2loop = bool(policy.get("apply_g8_delta_b3_2loop", False))
        delta_b3_candidates_raw = policy.get("delta_b3_candidates", None)
        delta_b3_candidates: list[float] = []
        if isinstance(delta_b3_candidates_raw, list):
            for v in delta_b3_candidates_raw:
                try:
                    delta_b3_candidates.append(float(v))
                except Exception:
                    continue
        if not delta_b3_candidates:
            delta_b3_candidates = [0.0, 1.0, 2.0, 3.0, 4.0]

        sm_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "sm_inputs_mz.json"
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        sm_inp = SmMzInputs(
            mu_GeV=float(sm_raw["mu_GeV"]),
            alpha_em_inv=float(sm_raw["alpha_em_inv"]),
            sin2_thetaW=float(sm_raw["sin2_thetaW"]),
            alpha_s=float(sm_raw["alpha_s"]),
        )
        mu0 = float(sm_inp.mu_GeV)
        if mu0 <= 0:
            raise ValueError("Invalid SM input scale mu_GeV")

        boundary = {"mode": "MZ", "route": "algebraic", "mu_GeV": mu0}
        if use_mt_boundary:
            sm_bc = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_raw)
            route = sm_bc.route_2loop if boundary_route != "pyrate_1loop" else sm_bc.route_1loop
            mu0 = float(route["mu_mt_GeV"])
            gY_mz = float(route["gY_mt"])
            g2_mz = float(route["g2_mt"])
            g3_mz = float(route["g3_mt"])
            Yu_seed = np.diag([0.0, 0.0, float(route["yt_mt"])]).astype(complex)
            Yd_seed = np.diag([0.0, 0.0, float(route["yb_mt"])]).astype(complex)
            Ye_seed = np.diag([0.0, 0.0, float(route["ytau_mt"])]).astype(complex)
            lam_seed = float(route["lambda_mt"])
            boundary = {
                "mode": "mt",
                "route": "pyrate_2loop" if boundary_route != "pyrate_1loop" else "pyrate_1loop",
                "mu_GeV": mu0,
                "yt_mt": float(route["yt_mt"]),
                "yb_mt": float(route["yb_mt"]),
                "ytau_mt": float(route["ytau_mt"]),
                "lambda_mt": float(route["lambda_mt"]),
            }
        else:
            # Initial gauge couplings at μ=MZ (algebraic)
            g1_gut_mz, g2_mz, g3_mz = gauge_couplings_from_mz_inputs(sm_inp)
            gY_mz = float(gY_from_g1_gut(g1_gut_mz))
            Yu_seed = np.diag([0.0, 0.0, 0.94]).astype(complex)
            Yd_seed = np.diag([0.0, 0.0, 0.017]).astype(complex)
            Ye_seed = np.diag([0.0, 0.0, 0.010]).astype(complex)
            lam_seed = 0.13

        gut_norm = float(g1_gut_over_gY())
        g1_gut_mz = float(gut_norm * gY_mz)

        # Thresholds (segment runner)
        thresholds_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "rge_thresholds_v25.json"
        thresholds_raw = json.loads(thresholds_path.read_text(encoding="utf-8"))
        thresholds_GeV = dict(thresholds_raw.get("thresholds_GeV", {}))

        # Beta sources (fail-fast via PyR@TE registry)
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

        # Run a small discrete policy scan (no continuous tuning):
        # - Δb3 patch above MG8 on/off (TFPT paper v1.06 note)
        # - gravity α^3 patch on/off (runner-level; controlled via two_loop_rg_fingerprints.json)
        variants: list[dict[str, object]] = []
        kappa_nonzero = bool((abs(k1) > 0) or (abs(k2) > 0) or (abs(k3) > 0))
        gravity_toggle = [False, True] if kappa_nonzero else [bool(grav_enabled)]
        g8_two_loop_tag = "g8_2l=on" if apply_g8_delta_b3_2loop else "g8_2l=off"
        # Discrete Δb3 candidates (integer multiplicities). This is the minimal “TFPT-style discrete search”
        # knob for unification without introducing a continuous tuner (policy-driven if provided).
        for apply_grav in gravity_toggle:
            # MG8 patch OFF (no Δb3)
            variants.append(
                {
                    "label": f"policy_g8=off_{g8_two_loop_tag}_grav={'on' if apply_grav else 'off'}",
                    "apply_g8_delta_b3": False,
                    "delta_b3_g8": 0.0,
                    "apply_g8_delta_b3_2loop": bool(apply_g8_delta_b3_2loop),
                    "apply_gravity_alpha3": bool(apply_grav),
                }
            )
            # MG8 patch ON (scan Δb3 integers)
            for db3 in delta_b3_candidates:
                db3_tag = f"{float(db3):g}"
                variants.append(
                    {
                        "label": f"policy_g8=on_db3={db3_tag}_{g8_two_loop_tag}_grav={'on' if apply_grav else 'off'}",
                        "apply_g8_delta_b3": True,
                        "delta_b3_g8": float(db3),
                        "apply_g8_delta_b3_2loop": bool(apply_g8_delta_b3_2loop),
                        "apply_gravity_alpha3": bool(apply_grav),
                    }
                )

        rows: list[dict[str, object]] = []
        best = None
        best_row = None
        best_rge2 = None

        for v in variants:
            apply_g8_delta_b3 = bool(v["apply_g8_delta_b3"])
            apply_gravity_alpha3 = bool(v["apply_gravity_alpha3"])
            apply_g8_delta_b3_2loop = bool(v.get("apply_g8_delta_b3_2loop", False))
            delta_b3_g8 = float(v.get("delta_b3_g8", 2.0))
            rge2 = run_flavor_rge_2loop_thresholds(
                mu_start_GeV=float(mu0),
                mu_end_GeV=float(mu_end),
                thresholds_GeV=thresholds_GeV,
                g_start=(float(gY_mz), float(g2_mz), float(g3_mz)),
                Yu_start=Yu_seed,
                Yd_start=Yd_seed,
                Ye_start=Ye_seed,
                lambda_start=float(lam_seed),
                yN_start=None,
                beta_sm=beta_sm,
                beta_e8=beta_e8,
                apply_sigma_threshold=True,
                apply_g8_delta_b3=bool(apply_g8_delta_b3),
                delta_b3_g8=float(delta_b3_g8),
                apply_g8_delta_b3_2loop=bool(apply_g8_delta_b3_2loop),
                apply_matching=True,
                matching_loop_order=1,
                matching_finite_delta_alpha=None,
                apply_gravity_alpha3=bool(apply_gravity_alpha3),
                gravity_kappa=(float(k1), float(k2), float(k3)),
                gravity_c3=float(1.0 / (8.0 * np.pi)),
                alpha_definition_for_u1=str(alpha_def_u1),
                max_loop=2,
                rtol=1e-9,
                atol=1e-11,
                method="DOP853",
                return_trajectory=True,
                trajectory_n_points=260,
            )

            traj = rge2.get("trajectory", {}) if isinstance(rge2.get("trajectory", {}), dict) else {}
            t = np.array(traj.get("log10_mu", []), dtype=float)
            g_sm = traj.get("g_sm", {}) if isinstance(traj.get("g_sm", {}), dict) else {}
            g1_sm = np.array(g_sm.get("gY", []), dtype=float)
            g2 = np.array(g_sm.get("g2", []), dtype=float)
            g3 = np.array(g_sm.get("g3", []), dtype=float)
            if t.size == 0 or g1_sm.size == 0 or g2.size == 0 or g3.size == 0:
                raise RuntimeError("Unification gate: missing trajectory sampling output")

            alphaY = _alpha_from_g(g1_sm)
            alpha1_gut = (5.0 / 3.0) * alphaY
            alpha2 = _alpha_from_g(g2)
            alpha3 = _alpha_from_g(g3)
            b = _best_unification_point(log10_mu=t, alpha1_gut=alpha1_gut, alpha2=alpha2, alpha3=alpha3)

            rec = {
                "label": str(v.get("label", "")),
                "apply_g8_delta_b3": bool(apply_g8_delta_b3),
                "delta_b3_g8": float(delta_b3_g8),
                "apply_g8_delta_b3_2loop": bool(apply_g8_delta_b3_2loop),
                "apply_gravity_alpha3": bool(apply_gravity_alpha3),
                "best": b.__dict__,
            }
            rows.append(rec)
            if best is None or (np.isfinite(b.mismatch_rel) and float(b.mismatch_rel) < float(best.mismatch_rel)):
                best = b
                best_row = dict(rec)
                best_rge2 = rge2

        assert best is not None and best_row is not None and best_rge2 is not None

        # Gate tolerance: 1% mismatch by default; physics mode is stricter.
        tol = 1.0e-2 if str(config.verification_mode) != "physics" else 5.0e-3
        passed = bool(np.isfinite(best.mismatch_rel) and best.mismatch_rel < tol)

        detail = f"mismatch_rel={best.mismatch_rel:.6g} @ mu*={best.mu_GeV:.3e} GeV (tol={tol:g}, mode={config.verification_mode})"
        if passed:
            checks = [mk_check_pass("unification_gate", detail)]
        else:
            checks = [mk_check_fail("unification_gate", detail) if str(config.verification_mode) == "physics" else mk_check_warn("unification_gate", detail)]

        lines: list[str] = []
        cfg_sha = _sha256_file(cfg_path) if cfg_path.exists() else "missing"
        sm_sha = _sha256_file(sm_path) if sm_path.exists() else "missing"
        thr_sha = _sha256_file(thresholds_path) if thresholds_path.exists() else "missing"
        lines += [
            "Unification gate (explicit 2-loop PyR@TE RG check)",
            "",
            f"config: {_relpath(cfg_path)} (sha256={cfg_sha})",
            f"SM inputs: {_relpath(sm_path)} (sha256={sm_sha})",
            f"thresholds: {_relpath(thresholds_path)} (sha256={thr_sha})",
            f"policy: {_relpath(policy_path)} (loaded={policy_loaded})",
            "",
            "Runner policy:",
            "- segment switch: SM→E8 at MSigma; Δb3 patch above MG8",
            "- matching: enabled (1-loop identity at thresholds unless finite pieces are provided)",
            f"- gravity α^3 patch (runner-level): configured enabled={bool(grav_enabled)} kappa=({k1},{k2},{k3}) alpha_def_u1={alpha_def_u1}",
            f"- G8 Δb3 2-loop patch: {'enabled' if apply_g8_delta_b3_2loop else 'disabled'} (assumption=Weyl adjoint)",
            "",
            f"boundary: mode={boundary.get('mode')} route={boundary.get('route')} mu0={mu0:.6g} GeV",
            *(
                [
                    f"boundary mt seeds: yt={boundary.get('yt_mt'):.6g}, yb={boundary.get('yb_mt'):.6g}, ytau={boundary.get('ytau_mt'):.6g}, lambda={boundary.get('lambda_mt'):.6g}"
                ]
                if boundary.get("mode") == "mt"
                else []
            ),
            f"initial @ mu0: gY={gY_mz:.6g}, g2={g2_mz:.6g}, g3={g3_mz:.6g} (g1_GUT={g1_gut_mz:.6g}, g1_GUT/gY={g1_gut_mz/gY_mz:.15g} vs sqrt(5/3)={gut_norm:.15g})",
            "",
            "Discrete policy scan (best per variant):",
        ]
        for r in sorted(rows, key=lambda x: float(((x.get("best") or {}).get("mismatch_rel")) or float("inf"))):  # type: ignore[arg-type]
            b = r.get("best", {}) if isinstance(r.get("best", {}), dict) else {}
            lines.append(
                f"- {r.get('label')}: mismatch_rel={b.get('mismatch_rel')} @ mu*={b.get('mu_GeV')} GeV "
                f"(g8_delta_b3={r.get('apply_g8_delta_b3')}, delta_b3_g8={r.get('delta_b3_g8')}, "
                f"g8_delta_b3_2loop={r.get('apply_g8_delta_b3_2loop')}, gravity_alpha3={r.get('apply_gravity_alpha3')})"
            )
        lines += [
            "",
            "Best-fit unification point (min over grid + quadratic refinement):",
            f"- mu* = {best.mu_GeV:.6e} GeV (log10 mu* = {best.log10_mu:.6f})",
            f"- alpha1_GUT(mu*) = {best.alpha1_gut:.6g}",
            f"- alpha2(mu*)      = {best.alpha2:.6g}",
            f"- alpha3(mu*)      = {best.alpha3:.6g}",
            f"- mismatch_rel(mu*) = {best.mismatch_rel:.6g}",
            "",
            f"Gate: {'PASS' if passed else 'FAIL'} (tol={tol:g}, mode={config.verification_mode})",
        ]

        return ModuleResult(
            results={
                "config": {"path": _relpath(cfg_path), "sha256": cfg_sha},
                "sm_inputs": {"path": _relpath(sm_path), "sha256": sm_sha},
                "thresholds": {"path": _relpath(thresholds_path), "sha256": thr_sha},
                "policy": {
                    "path": _relpath(policy_path),
                    "loaded": bool(policy_loaded),
                    "use_mt_boundary": bool(use_mt_boundary),
                    "boundary_route": str(boundary_route),
                    "delta_b3_candidates": [float(x) for x in delta_b3_candidates],
                    "apply_g8_delta_b3_2loop": bool(apply_g8_delta_b3_2loop),
                    "mu_end_GeV": float(mu_end),
                },
                "boundary": dict(boundary),
                "runner": {
                    "segments": best_rge2.get("segments"),
                    "matching": best_rge2.get("matching"),
                    "gravity_alpha3_patch": best_rge2.get("gravity_alpha3_patch"),
                    "g8_patch": best_rge2.get("g8_patch"),
                },
                "policy_scan": {"variants": rows, "best_variant": best_row},
                "best_unification_point": best.__dict__,
                "tolerance": {"mismatch_rel_tol": float(tol), "mode": str(config.verification_mode)},
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

