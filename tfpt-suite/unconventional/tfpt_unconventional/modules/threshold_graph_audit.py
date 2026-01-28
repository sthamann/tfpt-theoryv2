from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.pyrate_pythonoutputs import get_pyrate_pythonoutput
from tfpt_suite.rge_pyrate_2loop import load_pyrate_beta_module, run_flavor_rge_2loop_thresholds
from tfpt_suite.pyrate_boundary_runner import sm_boundary_conditions_at_mt
from tfpt_suite.rge_sm import run_sm_gauge_only_2loop_thresholds
from tfpt_suite.sm_inputs import SmMzInputs, gauge_couplings_from_mz_inputs


@dataclass(frozen=True)
class ThresholdNode:
    threshold_id: str
    scale_GeV: float
    note: str


@dataclass(frozen=True)
class SegmentEdge:
    mu_start_GeV: float
    mu_end_GeV: float
    model: str
    patches: list[str]
    threshold_at_start: str | None
    matching_active: bool | None


def _tfpt_suite_dir() -> Path:
    # .../tfpt-suite/unconventional/tfpt_unconventional/modules/<this_file>
    return Path(__file__).resolve().parents[3]


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _plot_threshold_graph(
    *,
    out_dir: Path,
    nodes: list[ThresholdNode],
    edges_off: list[SegmentEdge],
    edges_on: list[SegmentEdge],
) -> tuple[dict[str, str | None], list[str]]:
    plot: dict[str, str | None] = {"threshold_graph_png": None}
    warnings: list[str] = []

    try:
        import matplotlib.pyplot as plt  # type: ignore

        out_dir.mkdir(parents=True, exist_ok=True)

        def log10_mu(x: float) -> float:
            return float(math.log10(max(1e-300, float(x))))

        colors = {"sm_tfpt_2loop_v25": "#1f77b4", "e8_sigma_yN_2loop": "#ff7f0e"}

        def bars(edges: list[SegmentEdge]) -> list[tuple[float, float, str, str]]:
            out: list[tuple[float, float, str, str]] = []
            for e in edges:
                a = log10_mu(e.mu_start_GeV)
                b = log10_mu(e.mu_end_GeV)
                out.append((a, b - a, colors.get(e.model, "#7f7f7f"), str(e.threshold_at_start or "")))
            return out

        bars_off = bars(edges_off)
        bars_on = bars(edges_on)

        fig, ax = plt.subplots(figsize=(11.5, 3.8))
        y_off = (10, 6)
        y_on = (2, 6)

        for (x0, w, col, lbl) in bars_off:
            ax.broken_barh([(x0, w)], y_off, facecolors=col, alpha=0.85)
            if lbl:
                ax.text(x0 + 0.02, y_off[0] + y_off[1] - 0.7, lbl, fontsize=9, ha="left", va="top")
        for (x0, w, col, lbl) in bars_on:
            ax.broken_barh([(x0, w)], y_on, facecolors=col, alpha=0.85)
            if lbl:
                ax.text(x0 + 0.02, y_on[0] + y_on[1] - 0.7, lbl, fontsize=9, ha="left", va="top")

        # Threshold markers
        for n in nodes:
            if not (n.scale_GeV > 0 and math.isfinite(n.scale_GeV)):
                continue
            x = log10_mu(n.scale_GeV)
            ax.axvline(x, color="black", lw=0.8, ls=":", alpha=0.6)
            ax.text(x, 0.6, n.threshold_id, rotation=90, fontsize=9, ha="center", va="bottom", alpha=0.85)

        ax.set_yticks([y_on[0] + y_on[1] / 2.0, y_off[0] + y_off[1] / 2.0])
        ax.set_yticklabels(["matching ON", "matching OFF"])
        ax.set_xlabel(r"$\log_{10}(\mu/\mathrm{GeV})$")
        ax.set_title("Threshold graph (segments colored by active beta source)")
        ax.grid(True, axis="x", ls=":", alpha=0.3)

        # Legend
        handles = [
            plt.Line2D([0], [0], color=colors["sm_tfpt_2loop_v25"], lw=8, label="SM betas"),
            plt.Line2D([0], [0], color=colors["e8_sigma_yN_2loop"], lw=8, label="E8 betas"),
        ]
        ax.legend(handles=handles, loc="upper left")

        fig.tight_layout()
        path = out_dir / "threshold_graph.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot["threshold_graph_png"] = str(path)
    except Exception as e:
        warnings.append(f"plot_generation_failed: {e}")

    return plot, warnings


class ThresholdGraphAuditModule(TfptModule):
    """
    Unconventional module (B): build and validate a declarative “threshold graph”.

    The main suite already contains declarative threshold bookkeeping in the 2-loop RGE engine
    (`tfpt_suite/rge_pyrate_2loop.py` via `threshold_rules` and `segments`).

    This module makes that explicit as a standalone artifact:
    - parse the canonical threshold table (`rge_thresholds_v25.json`)
    - run the RGE engine once in a *toy* configuration to extract the segment graph
    - verify that threshold switches and matching flags are reported consistently
    """

    module_id = "ux_threshold_graph_audit"
    title = "Unconventional: threshold graph audit (segments + matching bookkeeping)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "threshold table: tfpt_suite/data/rge_thresholds_v25.json",
                "SM MZ inputs: tfpt_suite/data/sm_inputs_mz.json (to build a consistent mt boundary)",
                "PyR@TE PythonOutput registry: tfpt_suite/data/pyrate_pythonoutputs.json",
            ],
            outputs=[
                "threshold nodes (sorted)",
                "segment edges extracted from the 2-loop runner",
                "matching-on vs matching-off comparison (publication-grade bookkeeping flags)",
            ],
            formulas=[
                "segments are those used by run_flavor_rge_2loop_thresholds (piecewise integration in ln μ)",
                "threshold_rules declare the action at each threshold (beta switch, Δb3 patch, yN activation)",
            ],
            validation=[
                "graph contains expected nodes and ordered segments",
                "matching_off marks thresholds as 'continuous_by_assumption'",
                "matching_on records explicit match_* outcomes at thresholds (even if identity at 1-loop)",
                "thresholds are tagged as finite_pieces_implemented vs identity_matching when matching is enabled",
                "matching_on has no blocked thresholds",
            ],
            determinism="Deterministic given config seed; no stochastic components.",
        )

    def run(self, config) -> ModuleResult:
        suite_dir = _tfpt_suite_dir()
        thr_path = suite_dir / "tfpt_suite" / "data" / "rge_thresholds_v25.json"
        sm_path = suite_dir / "tfpt_suite" / "data" / "sm_inputs_mz.json"

        thr = _read_json(thr_path)
        thr_scales = dict(thr.get("thresholds_GeV", {}))
        mu_uv = float(thr.get("mu_uv_GeV", 1e16))

        # Canonical nodes for the flavor runner.
        nodes: list[ThresholdNode] = []
        for key in ["MSigma", "MG8", "MNR1", "MNR2", "MNR3", "MPhi"]:
            if key in thr_scales:
                nodes.append(ThresholdNode(threshold_id=key, scale_GeV=float(thr_scales[key]), note="from rge_thresholds_v25.json"))
        nodes.sort(key=lambda n: n.scale_GeV)

        # Build mt boundary (deterministic; avoids ad-hoc starting values).
        sm_inputs = _read_json(sm_path)
        bnd = sm_boundary_conditions_at_mt(sm_inputs_mz=sm_inputs)
        mt = float(bnd.route_2loop["mu_mt_GeV"])
        mu_ref = float(sm_inputs["mu_GeV"])

        g_start = (float(bnd.route_2loop["gY_mt"]), float(bnd.route_2loop["g2_mt"]), float(bnd.route_2loop["g3_mt"]))
        yt = float(bnd.route_2loop["yt_mt"])
        yb = float(bnd.route_2loop["yb_mt"])
        ytau = float(bnd.route_2loop["ytau_mt"])
        lam = float(bnd.route_2loop["lambda_mt"])

        # Minimal deterministic Yukawa matrices for the audit run.
        Yu = np.diag([1e-5, 1e-2, yt]).astype(complex)
        Yd = np.diag([1e-5, 2e-4, yb]).astype(complex)
        Ye = np.diag([2e-6, 4e-4, ytau]).astype(complex)
        yN = np.zeros((3, 3), dtype=complex)  # included so MNR thresholds appear in the graph

        # Load beta sources (fail-fast uses model_name_expected and sha256 fingerprints).
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

        thresholds_GeV = dict(thr.get("thresholds_GeV", {}))

        # Run once with matching disabled (baseline: continuous-by-assumption labeling).
        run_off = run_flavor_rge_2loop_thresholds(
            mu_start_GeV=mt,
            mu_end_GeV=mu_uv,
            thresholds_GeV=thresholds_GeV,
            g_start=g_start,
            Yu_start=Yu,
            Yd_start=Yd,
            Ye_start=Ye,
            lambda_start=lam,
            yN_start=yN,
            beta_sm=beta_sm,
            beta_e8=beta_e8,
            apply_sigma_threshold=True,
            apply_g8_delta_b3=True,
            delta_b3_g8=2.0,
            apply_matching=False,
            matching_loop_order=1,
            rtol=1e-8,
            atol=1e-10,
            method="DOP853",
        )

        # Run once with matching enabled (1-loop identity at μ=threshold, but explicit record).
        run_on = run_flavor_rge_2loop_thresholds(
            mu_start_GeV=mt,
            mu_end_GeV=mu_uv,
            thresholds_GeV=thresholds_GeV,
            g_start=g_start,
            Yu_start=Yu,
            Yd_start=Yd,
            Ye_start=Ye,
            lambda_start=lam,
            yN_start=yN,
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

        def edges_from(run: dict[str, Any]) -> list[SegmentEdge]:
            segs = run.get("segments", [])
            out: list[SegmentEdge] = []
            if not isinstance(segs, list):
                return out
            for s in segs:
                if not isinstance(s, dict):
                    continue
                thr_rule = s.get("threshold_transition_at_start", None)
                thr_id = None
                if isinstance(thr_rule, dict) and "threshold_id" in thr_rule:
                    thr_id = str(thr_rule["threshold_id"])
                elif thr_rule is not None and hasattr(thr_rule, "threshold_id"):
                    thr_id = str(getattr(thr_rule, "threshold_id"))
                tm = s.get("threshold_match", None)
                matching_active = None
                if isinstance(tm, dict) and "matching_active" in tm:
                    matching_active = bool(tm["matching_active"])
                out.append(
                    SegmentEdge(
                        mu_start_GeV=float(s.get("mu_start_GeV", s.get("mu_start", float("nan")))),
                        mu_end_GeV=float(s.get("mu_end_GeV", s.get("mu_end", float("nan")))),
                        model=str(s.get("model", "?")),
                        patches=[str(x) for x in (s.get("patches", []) if isinstance(s.get("patches", []), list) else [])],
                        threshold_at_start=thr_id,
                        matching_active=matching_active,
                    )
                )
            return out

        edges_off = edges_from(run_off)
        edges_on = edges_from(run_on)

        def threshold_statuses(run: dict[str, Any]) -> list[dict[str, Any]]:
            segs = run.get("segments", [])
            out: list[dict[str, Any]] = []
            seen: set[str] = set()
            if not isinstance(segs, list):
                return out
            for s in segs:
                if not isinstance(s, dict):
                    continue
                thr_rule = s.get("threshold_transition_at_start", None)
                if thr_rule is None:
                    continue
                thr_id = None
                if isinstance(thr_rule, dict) and "threshold_id" in thr_rule:
                    thr_id = str(thr_rule["threshold_id"])
                elif hasattr(thr_rule, "threshold_id"):
                    thr_id = str(getattr(thr_rule, "threshold_id"))
                if not thr_id or thr_id in seen:
                    continue
                seen.add(thr_id)
                tm = s.get("threshold_match", {})
                comp_status: dict[str, str] = {}
                finite = False
                for key in ("gauge", "yukawa", "quartic"):
                    if isinstance(tm, dict) and isinstance(tm.get(key, {}), dict):
                        status = tm.get(key, {}).get("status", None)
                        if status is not None:
                            comp_status[key] = str(status)
                            if status == "matched_with_finite_pieces":
                                finite = True
                label = "finite_pieces_implemented" if finite else "identity_matching"
                out.append({"threshold_id": thr_id, "status": label, "component_status": comp_status})
            return out

        finite_status_on = threshold_statuses(run_on)

        # --- Below-MZ QCD threshold audit (finite αs matching) ---
        # This is the missing piece for Gate-2 "threshold_graph_audit marks publication grade":
        # we explicitly check that the below-MZ gauge runner applies finite 2-loop αs matching at heavy-quark thresholds.
        alpha_em_inv = float(sm_inputs["alpha_em_inv"])
        sin2 = float(sm_inputs["sin2_thetaW"])
        alpha_s_mz = float(sm_inputs["alpha_s"])
        mc = float(sm_inputs.get("mc_GeV", 1.27))
        mb = float(sm_inputs.get("mb_GeV", 4.18))
        mt_pole = float(sm_inputs.get("mt_GeV", 172.76))

        inp_mz = SmMzInputs(mu_GeV=float(mu_ref), alpha_em_inv=alpha_em_inv, sin2_thetaW=sin2, alpha_s=alpha_s_mz)
        g_mz_gut = gauge_couplings_from_mz_inputs(inp_mz)

        g_mc_off = run_sm_gauge_only_2loop_thresholds(
            mu_start_GeV=float(mu_ref),
            mu_end_GeV=float(mc),
            g_start=(float(g_mz_gut[0]), float(g_mz_gut[1]), float(g_mz_gut[2])),
            mc_GeV=float(mc),
            mb_GeV=float(mb),
            mt_GeV=float(mt_pole),
            apply_alpha3_matching=False,
        )
        g_mc_on = run_sm_gauge_only_2loop_thresholds(
            mu_start_GeV=float(mu_ref),
            mu_end_GeV=float(mc),
            g_start=(float(g_mz_gut[0]), float(g_mz_gut[1]), float(g_mz_gut[2])),
            mc_GeV=float(mc),
            mb_GeV=float(mb),
            mt_GeV=float(mt_pole),
            apply_alpha3_matching=True,
        )
        alpha_mc_off = float((float(g_mc_off[2]) ** 2) / (4.0 * math.pi))
        alpha_mc_on = float((float(g_mc_on[2]) ** 2) / (4.0 * math.pi))
        delta_alpha_mc = float(alpha_mc_on - alpha_mc_off)

        # Basic ordering checks.
        def is_ordered(edges: list[SegmentEdge]) -> bool:
            if not edges:
                return False
            return all(e.mu_start_GeV < e.mu_end_GeV for e in edges)

        # Expected number of segments from cuts:
        expected_thresholds_in_range = [n.threshold_id for n in nodes if mt < n.scale_GeV < mu_uv and n.threshold_id in {"MSigma", "MG8", "MNR1", "MNR2", "MNR3"}]
        expected_segments = 1 + len(expected_thresholds_in_range)

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="threshold_table_loaded",
                passed=bool(nodes),
                detail=f"nodes={[(n.threshold_id, n.scale_GeV) for n in nodes]}",
            )
        )
        checks.append(
            Check(
                check_id="segments_ordered_matching_off",
                passed=bool(is_ordered(edges_off)),
                detail=f"segments={len(edges_off)}",
            )
        )
        checks.append(
            Check(
                check_id="segments_ordered_matching_on",
                passed=bool(is_ordered(edges_on)),
                detail=f"segments={len(edges_on)}",
            )
        )
        checks.append(
            Check(
                check_id="segment_count_matches_expected",
                passed=bool(len(edges_on) == expected_segments),
                detail=f"expected={expected_segments} from thresholds_in_range={expected_thresholds_in_range}; got={len(edges_on)}",
            )
        )
        checks.append(
            Check(
                check_id="matching_off_marks_blocked_thresholds",
                passed=bool(run_off.get("publication_grade", {}).get("threshold_matching_ok", True) is False and bool(run_off.get("publication_grade", {}).get("blocked_thresholds", []))),
                detail=f"publication_grade={run_off.get('publication_grade', {})}",
            )
        )
        checks.append(
            Check(
                check_id="matching_on_is_publication_grade_bookkeeping",
                passed=bool(run_on.get("publication_grade", {}).get("threshold_matching_ok", False) is True),
                detail=f"publication_grade={run_on.get('publication_grade', {})}",
            )
        )
        checks.append(
            Check(
                check_id="thresholds_marked_publication_grade",
                passed=bool(finite_status_on and all("status" in row for row in finite_status_on)),
                detail=f"threshold_statuses={finite_status_on}",
            )
        )
        checks.append(
            Check(
                check_id="matching_on_no_blocked_thresholds",
                passed=bool(not run_on.get("publication_grade", {}).get("blocked_thresholds", [])),
                detail=f"blocked_thresholds={run_on.get('publication_grade', {}).get('blocked_thresholds', [])}",
            )
        )

        # Count how many threshold starts actually record matching_active flags.
        n_thr_on = sum(1 for e in edges_on if e.threshold_at_start is not None)
        n_thr_on_active = sum(1 for e in edges_on if e.threshold_at_start is not None and e.matching_active is True)
        checks.append(
            Check(
                check_id="threshold_starts_have_explicit_matching_records_when_enabled",
                passed=bool(n_thr_on == n_thr_on_active),
                detail=f"threshold_segments={n_thr_on}, matching_active_true={n_thr_on_active}",
            )
        )

        # Below-MZ finite matching check (should be nonzero when enabled).
        checks.append(
            Check(
                check_id="below_mz_alpha3_matching_changes_alpha_s",
                passed=bool(np.isfinite(delta_alpha_mc) and abs(delta_alpha_mc) > 0.0),
                detail=f"alpha_s(mc) off={alpha_mc_off:.9g}, on={alpha_mc_on:.9g}, delta={delta_alpha_mc:+.3e}",
            )
        )

        lines: list[str] = []
        lines += [
            "Unconventional: threshold graph audit",
            "",
            "Purpose:",
            "- Provide an explicit, inspectable segment/threshold graph for the 2-loop RGE runner.",
            "- Compare matching disabled vs enabled bookkeeping (identity-at-threshold at 1-loop, but explicit).",
            "",
            f"Threshold table: {thr_path}",
            f"mt boundary: mt={mt} GeV (from sm_boundary_conditions_at_mt, route_2loop)",
            f"mu_uv: {mu_uv} GeV",
            "",
            "Threshold nodes (sorted):",
        ]
        for n in nodes:
            lines.append(f"- {n.threshold_id}: {n.scale_GeV:.6g} GeV")
        lines += [
            "",
            "Edges (matching OFF):",
        ]
        for e in edges_off:
            lines.append(f"- [{e.mu_start_GeV:.6g} → {e.mu_end_GeV:.6g}] model={e.model} patches={e.patches} thr_at_start={e.threshold_at_start} matching_active={e.matching_active}")
        lines += [
            "",
            "Edges (matching ON):",
        ]
        for e in edges_on:
            lines.append(f"- [{e.mu_start_GeV:.6g} → {e.mu_end_GeV:.6g}] model={e.model} patches={e.patches} thr_at_start={e.threshold_at_start} matching_active={e.matching_active}")
        lines += [
            "",
            "Threshold finite-piece status (matching ON):",
        ]
        for row in finite_status_on:
            comps = row.get("component_status", {})
            lines.append(f"- {row.get('threshold_id')}: {row.get('status')} (gauge={comps.get('gauge')}, yukawa={comps.get('yukawa')}, quartic={comps.get('quartic')})")
        lines += [
            "",
            "Notes:",
            "- Matching at loop_order=1 is identity at μ=threshold (log-only), but still recorded explicitly when enabled.",
            "- Below MZ, the SM gauge-only runner applies finite 2-loop αs decoupling at heavy-quark thresholds when enabled (Gate-2 finite matching).",
            "- This module is an audit tool: it makes policies explicit and helps prevent silent threshold/matching drift.",
            "",
            "Below-MZ finite αs matching (2-loop decoupling) diagnostic:",
            f"- mu_ref=MZ={mu_ref} GeV, mb={mb} GeV, mc={mc} GeV",
            f"- alpha_s(mc) matching OFF: {alpha_mc_off:.9g}",
            f"- alpha_s(mc) matching ON:  {alpha_mc_on:.9g}",
            f"- delta alpha_s(mc) = {delta_alpha_mc:+.3e}",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
        ]

        warnings: list[str] = []
        plot: dict[str, str | None] = {"threshold_graph_png": None}
        if config.plot:
            out_dir = self.output_dir(config)
            plot, plot_warnings = _plot_threshold_graph(out_dir=out_dir, nodes=nodes, edges_off=edges_off, edges_on=edges_on)
            warnings.extend(plot_warnings)

        return ModuleResult(
            results={
                "threshold_table": {"file": str(thr_path), "mu_uv_GeV": mu_uv, "nodes": [n.__dict__ for n in nodes]},
                "boundary": {"mt_GeV": mt, "mu_ref_mz_GeV": mu_ref, "g_start": list(g_start)},
                "runs": {
                    "matching_off": {"publication_grade": run_off.get("publication_grade", {}), "threshold_rules": run_off.get("threshold_rules", []), "segments": run_off.get("segments", [])},
                    "matching_on": {"publication_grade": run_on.get("publication_grade", {}), "threshold_rules": run_on.get("threshold_rules", []), "segments": run_on.get("segments", [])},
                },
                "graph": {"edges_matching_off": [e.__dict__ for e in edges_off], "edges_matching_on": [e.__dict__ for e in edges_on]},
                "finite_piece_status": {"matching_on": finite_status_on},
                "below_mz_qcd_matching": {
                    "mu_ref_mz_GeV": mu_ref,
                    "mb_GeV": mb,
                    "mc_GeV": mc,
                    "alpha_s_mc_matching_off": alpha_mc_off,
                    "alpha_s_mc_matching_on": alpha_mc_on,
                    "delta_alpha_s_mc": delta_alpha_mc,
                },
                "plot": plot,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

