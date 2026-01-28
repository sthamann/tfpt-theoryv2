from __future__ import annotations

import json
import math
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.constants import TfptConstants
from tfpt_suite.flavor_textures import left_unitary_from_yukawa
from tfpt_suite.mobius_z3_yukawa_generator import generate_quark_yukawas_mt
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.pyrate_pythonoutputs import get_pyrate_pythonoutput
from tfpt_suite.rge_pyrate_2loop import load_pyrate_beta_module, run_flavor_rge_2loop_thresholds
from tfpt_suite.pyrate_boundary_runner import sm_boundary_conditions_at_mt


def _plot_flavor_holdout(
    *,
    out_dir: Path,
    rows: list[dict[str, Any]],
    fit_keys: list[str],
    holdout_keys: list[str],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"flavor_holdout_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)
        if not rows:
            return plot, warnings

        best = rows[0]

        # Scatter: chi2_fit vs chi2_holdout, colored by complexity.
        xs = np.array([float(r.get("chi2_fit", float("nan"))) for r in rows], dtype=float)
        ys = np.array([float(r.get("chi2_holdout", float("nan"))) for r in rows], dtype=float)
        cs = np.array([int(r.get("candidate", {}).get("complexity", 0)) for r in rows], dtype=float)
        ds = [str(r.get("candidate", {}).get("delta_source", "?")) for r in rows]

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 4.8))

        markers = {"delta_star": "o", "tau_mu": "^"}
        # Normalize colors to complexity range.
        cmin = float(np.min(cs)) if cs.size else 0.0
        cmax = float(np.max(cs)) if cs.size else 1.0
        norm = (cs - cmin) / (cmax - cmin + 1e-12)
        cmap = plt.get_cmap("viridis")

        for i in range(len(rows)):
            if not (np.isfinite(xs[i]) and np.isfinite(ys[i])):
                continue
            ax1.scatter(
                [xs[i]],
                [ys[i]],
                s=70,
                marker=markers.get(ds[i], "o"),
                color=cmap(float(norm[i])),
                edgecolor="black",
                linewidth=0.6,
                alpha=0.9,
            )

        ax1.set_xlabel(r"$\chi^2_{\mathrm{fit}}$ (7 CKM entries)")
        ax1.set_ylabel(r"$\chi^2_{\mathrm{holdout}}$ (Vub,Vtd)")
        ax1.set_title("Discrete convention scan (holdout protocol)")
        ax1.grid(True, ls=":", alpha=0.35)

        # Annotate best point.
        bx = float(best["chi2_fit"])
        by = float(best["chi2_holdout"])
        label = best.get("candidate", {})
        ax1.scatter([bx], [by], s=160, marker="*", color="red", edgecolor="black", linewidth=0.8, zorder=5)
        ax1.text(
            bx,
            by,
            f" best\n {label}",
            fontsize=9,
            ha="left",
            va="bottom",
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
        )

        # Colorbar for complexity.
        sm = plt.cm.ScalarMappable(cmap=cmap)
        sm.set_array(cs)
        cbar = fig.colorbar(sm, ax=ax1, fraction=0.046, pad=0.04)
        cbar.set_label("complexity (proxy)")

        # Bar plot: holdout contributions for the best candidate.
        contrib_hold = best.get("contrib_holdout", [])
        labels_h = []
        vals_h = []
        for c in contrib_hold:
            try:
                labels_h.append(str(c["key"]))
                vals_h.append(float(c["chi2"]))
            except Exception:
                continue
        if labels_h:
            ax2.bar(labels_h, vals_h, color="#d62728", alpha=0.9)
            ax2.set_title("Holdout χ² contributions (best candidate)")
            ax2.set_ylabel(r"$\chi^2$ contribution")
            ax2.grid(True, axis="y", ls=":", alpha=0.35)
        else:
            ax2.text(0.5, 0.5, "no holdout contributions", ha="center", va="center")

        fig.tight_layout()
        path = out_dir / "flavor_holdout.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["flavor_holdout_png"] = str(path)

    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


@dataclass(frozen=True)
class Candidate:
    delta_source: str
    s13_mode: str
    delta_mode: str
    complexity: int


def _tfpt_suite_dir() -> Path:
    # .../tfpt-suite/unconventional/tfpt_unconventional/modules/<this_file>
    return Path(__file__).resolve().parents[3]


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _ckm_from_yukawas(Yu: np.ndarray, Yd: np.ndarray) -> np.ndarray:
    Uu = left_unitary_from_yukawa(Yu)
    Ud = left_unitary_from_yukawa(Yd)
    return (Uu.conj().T @ Ud).astype(complex)


def _Vabs(V: np.ndarray) -> dict[str, float]:
    return {
        k: float(np.abs(V[i, j]))
        for (k, (i, j)) in {
            "Vud": (0, 0),
            "Vus": (0, 1),
            "Vub": (0, 2),
            "Vcd": (1, 0),
            "Vcs": (1, 1),
            "Vcb": (1, 2),
            "Vtd": (2, 0),
            "Vts": (2, 1),
            "Vtb": (2, 2),
        }.items()
    }


def _jarlskog(V: np.ndarray) -> float:
    # One standard invariant: J = Im(Vud Vcs Vus* Vcd*)
    return float(np.imag(V[0, 0] * V[1, 1] * np.conj(V[0, 1]) * np.conj(V[1, 0])))


def _chi2(Vabs: dict[str, float], ref: dict[str, Any], *, keys: list[str]) -> tuple[float, list[dict[str, float | str]]]:
    contrib: list[dict[str, float | str]] = []
    chi2 = 0.0
    for k in keys:
        cfg = ref[k]
        mean = float(cfg["mean"])
        sigma = float(cfg["sigma"])
        pred = float(Vabs[k])
        c2 = ((pred - mean) / sigma) ** 2 if sigma > 0 else float("inf")
        chi2 += c2
        contrib.append({"key": k, "pred": pred, "mean": mean, "sigma": sigma, "chi2": c2})
    contrib.sort(key=lambda x: float(x["chi2"]), reverse=True)
    return float(chi2), contrib


def _complexity(delta_source: str, s13_mode: str, delta_mode: str) -> int:
    # Minimal, explicit complexity proxy (not physics):
    # - prefer simple π conventions; penalize "2π" and fixed phases
    base = 1
    base += 0 if delta_source == "delta_star" else 1
    base += 1 if s13_mode == "A_lam3_over_3" else 2
    if delta_mode in {"pi_times_delta", "pi_times_1_minus_delta"}:
        base += 1
    elif delta_mode == "2pi_times_delta":
        base += 2
    else:  # koide_pi_over_12
        base += 3
    return int(base)


def _perm_complexity(perm: list[int]) -> int:
    """
    Minimal inversion-count proxy for permutation complexity (0..3 for 3×3).
    """
    p = [int(x) for x in perm]
    inv = 0
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            inv += 1 if p[i] > p[j] else 0
    return int(inv)


def _pmns_sin2(theta_deg: float) -> float:
    th = math.radians(float(theta_deg))
    s = math.sin(th)
    return float(s * s)


def _pmns_chi2_split(
    *,
    angles_deg: dict[str, Any],
    ref: dict[str, Any],
    fit_keys: list[str],
    holdout_keys: list[str],
) -> tuple[float, float, float, list[dict[str, float | str]], list[dict[str, float | str]]]:
    """
    Compute PMNS χ² with a holdout split.

    Reference table uses sin²θij for angles and δCP in degrees.
    """
    pred = {
        "sin2_theta12": _pmns_sin2(float(angles_deg["theta12_deg"])),
        "sin2_theta13": _pmns_sin2(float(angles_deg["theta13_deg"])),
        "sin2_theta23": _pmns_sin2(float(angles_deg["theta23_deg"])),
        "delta_cp_deg": float(angles_deg["delta_cp_deg"]),
    }

    def contrib(key: str) -> dict[str, float | str]:
        cfg = ref[key]
        mean = float(cfg["mean"])
        sigma = float(cfg["sigma"])
        x = float(pred[key])
        c2 = ((x - mean) / sigma) ** 2 if sigma > 0 else float("inf")
        return {"key": key, "pred": x, "mean": mean, "sigma": sigma, "chi2": float(c2)}

    fit = [contrib(k) for k in fit_keys]
    hold = [contrib(k) for k in holdout_keys]
    chi2_fit = float(np.nansum([float(t["chi2"]) for t in fit]))
    chi2_hold = float(np.nansum([float(t["chi2"]) for t in hold]))
    return chi2_fit, chi2_hold, float(chi2_fit + chi2_hold), fit, hold


class FlavorHoldoutSearchModule(TfptModule):
    """
    Unconventional module (C): discrete phase-map / CKM construction search with a holdout principle.

    This module does NOT introduce continuous fitting.
    It evaluates a fixed grid of discrete conventions and reports:
      - “fit” χ² on a chosen subset of CKM entries,
      - “holdout” χ² on the remaining entries.

    This is meant to prevent silent convention drift and to keep the “no free parameters” story honest.
    """

    module_id = "ux_flavor_holdout_search"
    title = "Unconventional: flavor holdout search (discrete CKM phase-map conventions)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "CKM reference snapshot: tfpt_suite/data/ckm_reference.json (diagnostic, scale-labeled)",
                "SM mt boundary conditions: tfpt_suite/rge_sm.py (from sm_inputs_mz.json)",
                "PyR@TE 2-loop runner: tfpt_suite/rge_pyrate_2loop.py (mt→μ_ref short run)",
                "discrete convention space: {delta_source}×{s13_mode}×{delta_mode}",
            ],
            outputs=[
                "ranked candidates by χ²_fit (with χ²_holdout reported separately)",
                "top-k candidates with full per-entry contributions",
            ],
            formulas=[
                "holdout split: fit on 7 CKM entries; hold out the small elements Vub and Vtd",
                "no continuous tuning; only discrete convention choices",
            ],
            validation=[
                "all candidates are evaluated deterministically and yield finite χ² values",
                "report the best candidate under the holdout protocol",
            ],
            determinism="Deterministic (finite discrete scan; no randomness).",
        )

    def run(self, config) -> ModuleResult:
        suite_dir = _tfpt_suite_dir()
        ref_path = suite_dir / "tfpt_suite" / "data" / "ckm_reference.json"
        lep_path = suite_dir / "tfpt_suite" / "data" / "lepton_masses_pdg.json"
        pmns_ref_path = suite_dir / "tfpt_suite" / "data" / "pmns_reference.json"
        sm_path = suite_dir / "tfpt_suite" / "data" / "sm_inputs_mz.json"
        thr_path = suite_dir / "tfpt_suite" / "data" / "rge_thresholds_v25.json"

        ref_raw = _read_json(ref_path)
        ref = ref_raw["matrix_abs"]
        mu_ref = float(ref_raw["reference"]["reference_scale_GeV"])

        lep = _read_json(lep_path)
        mtau = float(lep["masses"]["tau"]["mean"])
        mmu = float(lep["masses"]["muon"]["mean"])
        R = float(np.sqrt(mtau / mmu))
        delta_M = float((R - 1.0) / (R + 1.0))

        c = TfptConstants.compute()
        delta_star = float(c.delta_star)

        # Boundary conditions at mt
        sm_inputs = _read_json(sm_path)
        bnd = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_inputs)
        mt = float(bnd.route_2loop["mu_mt_GeV"])
        g_start = (float(bnd.route_2loop["gY_mt"]), float(bnd.route_2loop["g2_mt"]), float(bnd.route_2loop["g3_mt"]))
        yt = float(bnd.route_2loop["yt_mt"])
        yb = float(bnd.route_2loop["yb_mt"])
        ytau = float(bnd.route_2loop["ytau_mt"])
        lam = float(bnd.route_2loop["lambda_mt"])

        Ye = np.diag([2e-6, 4e-4, ytau]).astype(complex)

        thresholds_GeV = dict(_read_json(thr_path).get("thresholds_GeV", {}))

        # Load beta sources once (even though mt→MZ stays in the SM segment).
        sm_model = get_pyrate_pythonoutput("sm_tfpt_2loop_v25")
        e8_model = get_pyrate_pythonoutput("e8_sigma_yN_2loop")
        beta_sm = load_pyrate_beta_module(
            kind="sm_tfpt_2loop_v25",
            pythonoutput_dir=sm_model.pythonoutput_dir,
            model_name_expected=sm_model.model_name_expected,
            yaml_source=sm_model.yaml_source,
        )
        beta_e8 = load_pyrate_beta_module(
            kind="e8_sigma_yN_2loop",
            pythonoutput_dir=e8_model.pythonoutput_dir,
            model_name_expected=e8_model.model_name_expected,
            yaml_source=e8_model.yaml_source,
        )

        # Candidate grid (discrete; no continuous fit).
        delta_sources = ["delta_star", "tau_mu"]
        s13_modes = ["A_lam3_over_3", "A_lam3_times_1_minus_delta"]
        delta_modes = ["pi_times_delta", "pi_times_1_minus_delta", "2pi_times_delta", "koide_pi_over_12"]
        candidates = [Candidate(ds, s13, dm, _complexity(ds, s13, dm)) for ds in delta_sources for s13 in s13_modes for dm in delta_modes]

        fit_keys = ["Vud", "Vus", "Vcd", "Vcs", "Vcb", "Vts", "Vtb"]
        holdout_keys = ["Vub", "Vtd"]

        rows: list[dict[str, Any]] = []
        for cand in candidates:
            delta_used = float(delta_star if cand.delta_source == "delta_star" else delta_M)
            Yu_mt, Yd_mt, meta = generate_quark_yukawas_mt(
                varphi0=float(c.varphi0),
                delta_used=float(delta_used),
                yt_mt=float(yt),
                yb_mt=float(yb),
                reference_scale_GeV=float(mt),
                s13_mode=cand.s13_mode,  # type: ignore[arg-type]
                delta_mode=cand.delta_mode,  # type: ignore[arg-type]
            )
            V_mt = _ckm_from_yukawas(Yu_mt, Yd_mt)

            # Run mt→mu_ref (down) with 2-loop runner (no thresholds in this range).
            rge = run_flavor_rge_2loop_thresholds(
                mu_start_GeV=float(mt),
                mu_end_GeV=float(mu_ref),
                thresholds_GeV=thresholds_GeV,
                g_start=g_start,
                Yu_start=Yu_mt,
                Yd_start=Yd_mt,
                Ye_start=Ye,
                lambda_start=float(lam),
                yN_start=None,
                beta_sm=beta_sm,
                beta_e8=beta_e8,
                apply_sigma_threshold=True,
                apply_g8_delta_b3=False,
                apply_matching=False,
                matching_loop_order=1,
                rtol=1e-9,
                atol=1e-11,
                method="DOP853",
            )
            Yu_ref = np.array(rge["Yu_end"], dtype=complex)
            Yd_ref = np.array(rge["Yd_end"], dtype=complex)
            V_ref = _ckm_from_yukawas(Yu_ref, Yd_ref)

            vabs_ref = _Vabs(V_ref)
            vabs_mt = _Vabs(V_mt)
            chi2_fit, contrib_fit = _chi2(vabs_ref, ref, keys=fit_keys)
            chi2_hold, contrib_hold = _chi2(vabs_ref, ref, keys=holdout_keys)
            chi2_total = float(chi2_fit + chi2_hold)

            unitarity_dev = float(np.max(np.abs(V_ref.conj().T @ V_ref - np.eye(3))))

            rows.append(
                {
                    "candidate": cand.__dict__,
                    "delta_M": delta_M,
                    "delta_star": delta_star,
                    "delta_used": delta_used,
                    "ckm_mt_abs": vabs_mt,
                    "ckm_ref_abs": vabs_ref,
                    "J_mt": _jarlskog(V_mt),
                    "J_ref": _jarlskog(V_ref),
                    "chi2_fit": chi2_fit,
                    "chi2_holdout": chi2_hold,
                    "chi2_total": chi2_total,
                    "contrib_fit_top": contrib_fit[:5],
                    "contrib_holdout": contrib_hold,
                    "unitarity_dev_ref": unitarity_dev,
                }
            )

        rows.sort(key=lambda r: (float(r["chi2_fit"]), int(r["candidate"]["complexity"]), float(r["chi2_holdout"])))
        top = rows[:8]

        checks = [
            Check(
                check_id="evaluated_all_candidates",
                passed=bool(len(rows) == len(candidates)),
                detail=f"candidates={len(candidates)}, rows={len(rows)}",
            ),
            Check(
                check_id="best_candidate_has_finite_chi2",
                passed=bool(np.isfinite(float(rows[0]["chi2_total"]))),
                detail=f"best chi2_fit={rows[0]['chi2_fit']:.3f}, chi2_holdout={rows[0]['chi2_holdout']:.3f}, total={rows[0]['chi2_total']:.3f}",
            ),
        ]

        # --- PMNS holdout (new): avoid silent convention shopping for δ_CP ---
        pmns_holdout: dict[str, Any] = {}
        try:
            from tfpt_suite.config import SuiteConfig
            from tfpt_suite.modules.pmns_full_pipeline import PmnsFullPipelineModule

            pmns_ref_raw = _read_json(pmns_ref_path)
            pmns_refs = {
                "normal": dict(pmns_ref_raw.get("normal_ordering", {})),
                "inverted": dict(pmns_ref_raw.get("inverted_ordering", {})),
            }
            pmns_fit_keys = ["sin2_theta12", "sin2_theta13", "sin2_theta23"]
            pmns_hold_keys = ["delta_cp_deg"]

            # Run once in a temp dir to obtain the discrete convention candidates (6 permutations).
            with tempfile.TemporaryDirectory() as tmp:
                cfg_pmns = SuiteConfig(
                    output_dir=Path(tmp),
                    mp_dps=int(getattr(config, "mp_dps", 60) or 60),
                    seed=int(getattr(config, "seed", 0) or 0),
                    overwrite=True,
                    plot=False,
                    verification_mode="engineering",
                )
                pmns_mod = PmnsFullPipelineModule()
                pmns_res = pmns_mod.run_and_write(config=cfg_pmns)
                pmns_payload = json.loads((Path(tmp) / "pmns_full_pipeline" / "results.json").read_text(encoding="utf-8"))

            pmns_mt = pmns_payload.get("results", {}).get("pmns_mt", {}) if isinstance(pmns_payload, dict) else {}
            cand_list = pmns_mt.get("candidates", []) if isinstance(pmns_mt.get("candidates", []), list) else []

            pmns_rows: list[dict[str, Any]] = []
            for cand in cand_list:
                if not isinstance(cand, dict):
                    continue
                perm = [int(x) for x in (cand.get("perm", []) if isinstance(cand.get("perm", []), list) else [])]
                ang = cand.get("angles_deg", {}) if isinstance(cand.get("angles_deg", {}), dict) else {}
                if not ang:
                    continue

                # Compute holdout split for both orderings; choose best ordering by total χ².
                best = None
                for ordering, ref_o in pmns_refs.items():
                    if not isinstance(ref_o, dict):
                        continue
                    chi2_fit, chi2_hold, chi2_total, contrib_fit, contrib_hold = _pmns_chi2_split(
                        angles_deg=ang, ref=ref_o, fit_keys=pmns_fit_keys, holdout_keys=pmns_hold_keys
                    )
                    row = {
                        "perm": perm,
                        "ordering": ordering,
                        "complexity": _perm_complexity(perm),
                        "angles_deg": ang,
                        "chi2_fit": float(chi2_fit),
                        "chi2_holdout": float(chi2_hold),
                        "chi2_total": float(chi2_total),
                        "contrib_fit": contrib_fit,
                        "contrib_holdout": contrib_hold,
                    }
                    if best is None or float(row["chi2_total"]) < float(best["chi2_total"]):
                        best = row
                if best is not None:
                    pmns_rows.append(best)

            pmns_rows.sort(key=lambda r: (float(r["chi2_fit"]), int(r["complexity"]), float(r["chi2_holdout"])))
            pmns_top = pmns_rows[:6]

            pmns_holdout = {
                "reference": {"file": str(pmns_ref_path), "fit_keys": pmns_fit_keys, "holdout_keys": pmns_hold_keys},
                "rows": pmns_rows,
                "top": pmns_top,
            }

            checks.append(
                Check(
                    check_id="pmns_holdout_evaluated",
                    passed=bool(len(pmns_rows) >= 1),
                    detail=f"pmns_candidates={len(pmns_rows)} (from 6 permutation scan; best chi2_fit={pmns_rows[0]['chi2_fit']:.3f} chi2_holdout={pmns_rows[0]['chi2_holdout']:.3f})",
                )
            )
        except Exception as e:
            pmns_holdout = {"error": str(e), "note": "PMNS holdout is best-effort; rerun with pmns_full_pipeline dependencies available."}
            checks.append(Check(check_id="pmns_holdout_evaluated", passed=False, detail=f"failed: {e}"))

        if isinstance(pmns_holdout.get("top", None), list) and pmns_holdout.get("top"):
            pmns_best_line = f"- best PMNS candidate: {pmns_holdout.get('top', [])[0]}"
        else:
            pmns_best_line = f"- pmns_holdout: {pmns_holdout}"

        lines: list[str] = []
        lines += [
            "Unconventional: flavor holdout search (discrete CKM conventions)",
            "",
            "Holdout protocol:",
            f"- fit keys: {fit_keys}",
            f"- holdout keys: {holdout_keys}",
            "",
            f"Reference: {ref_path} (mu_ref={mu_ref} GeV)",
            f"Boundary: mt={mt} GeV; delta_M={delta_M:.6g}; delta_star={delta_star:.6g}",
            "",
            "Top candidates (sorted by chi2_fit, then complexity, then chi2_holdout):",
        ]
        for r in top:
            cand = r["candidate"]
            lines.append(
                f"- {cand}  chi2_fit={r['chi2_fit']:.3f}  chi2_holdout={r['chi2_holdout']:.3f}  chi2_total={r['chi2_total']:.3f}  unitarity_dev={r['unitarity_dev_ref']:.3e}"
            )
        lines += [
            "",
            "Best candidate details (per-entry contributions):",
            f"- candidate: {rows[0]['candidate']}",
            f"- contrib_fit_top: {rows[0]['contrib_fit_top']}",
            f"- contrib_holdout: {rows[0]['contrib_holdout']}",
            "",
            "PMNS holdout (new; avoids silent δCP convention shopping):",
            f"- reference: {pmns_ref_path}",
            f"- fit keys: sin2_theta12,sin2_theta13,sin2_theta23; holdout key: delta_cp_deg",
            f"- status: {'OK' if ('rows' in pmns_holdout and pmns_holdout['rows']) else 'FAILED/EMPTY'}",
            pmns_best_line,
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module is diagnostic. It does not claim an operator-level topology→phase derivation.",
            "- It enforces a holdout split to reduce the risk of 'convention shopping' that overfits the full CKM table.",
            "- PMNS holdout is implemented as a best-effort scan over the discrete PMNS permutation conventions emitted by pmns_full_pipeline.",
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"flavor_holdout_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_flavor_holdout(out_dir=out_dir, rows=rows, fit_keys=fit_keys, holdout_keys=holdout_keys)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "reference": {"file": str(ref_path), "mu_ref_GeV": mu_ref, "fit_keys": fit_keys, "holdout_keys": holdout_keys},
                "deltas": {"delta_M": delta_M, "delta_star": delta_star},
                "scan": {"candidates": [c.__dict__ for c in candidates], "rows": rows, "top": top},
                "pmns_holdout": pmns_holdout,
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

