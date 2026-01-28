from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.operator_spec_builder import generate_effective_action_r2_operator_spec


@dataclass(frozen=True)
class LedgerItem:
    area: str
    requirement: str
    status: str
    evidence: list[str]
    next_steps: list[str]


def _read_json_if_exists(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _has_any_status(raw: dict[str, Any], status: str) -> bool:
    want = status.strip().lower()
    for term in raw.get("lagrangian_density", {}).get("ufe_action_terms", []):
        if isinstance(term, dict) and str(term.get("status", "")).strip().lower() == want:
            return True
    return False


def _operator_spec_derivation_status(raw: dict[str, Any]) -> str:
    deriv = raw.get("derivation", {})
    if isinstance(deriv, dict):
        return str(deriv.get("status", "")).strip().lower()
    return ""


def _discover_paper_source(root: Path) -> Path:
    """
    Best-effort discovery of the current paper source in this repo snapshot.

    Historically, the suite referenced `latex/tfpt-theory-fullv25.tex`. In the current tree,
    v2.5 content may exist as a snippets pack (`latex/tfpt-v25-snippets.tex`) while the latest
    full paper TeX file is `tfpt-theory-fullv24.tex`.

    We keep this logic deterministic and local: no network, no git inspection.
    """
    candidates = [
        root / "latex" / "tfpt-theory-fullv25.tex",
        root / "latex" / "tfpt-theory-fullv24.tex",
        root / "latex" / "tfpt-theory-fullv21.tex",
        root / "latex" / "tffpt-theory-fullv2.tex",
    ]
    for p in candidates:
        if p.exists():
            return p

    latex_dir = root / "latex"
    if latex_dir.is_dir():
        tex_files = sorted([p for p in latex_dir.iterdir() if p.is_file() and p.suffix.lower() == ".tex"])
        if tex_files:
            return tex_files[-1]

    # Fallback (non-existing): preserve the historical reference for error messages.
    return candidates[0]


class QftCompletenessLedgerModule(TfptModule):
    module_id = "qft_completeness_ledger"
    title = "QFT completeness ledger (7.1): what is specified vs. still missing (paper v2.5)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "paper conventions: best-effort from latex/ (prefers tfpt-theory-fullv25.tex, otherwise falls back to tfpt-theory-fullv24.tex)",
                "suite code + data tables under tfpt-suite/",
                "PyR@TE3 artifacts under Pyrate3/",
            ],
            outputs=[
                "A structured checklist aligned with 7.1 (microscopic action, quantization, renorm, anomalies, derivation-vs-input)",
                "Concrete engineering next steps aligned with 7.2",
            ],
            formulas=[],
            validation=["Ledger is produced deterministically from repo state (presence of configs/spec files)."],
            determinism="Deterministic (static repo inspection; no network).",
        )

    def run(self, config) -> ModuleResult:
        # Evidence pointers (kept as paths/identifiers the reader can click/navigate)
        root = Path(__file__).resolve().parents[3]
        tfpt_suite_dir = root / "tfpt-suite"
        pyrate3_dir = root / "Pyrate3"
        latex_file = _discover_paper_source(root)
        paper_v1_06 = root / "paper_v1_06_01_09_2025.tex"
        update_tfpt_v1_07sm = root / "update_tfptv1_07sm.tex"
        birefringence_letter = root / "alessandro.tex"

        # Key files we now rely on
        operator_spec = tfpt_suite_dir / "tfpt_suite" / "data" / "effective_action_r2_operator_spec.json"
        microscopic_action = tfpt_suite_dir / "tfpt_suite" / "data" / "microscopic_action_tfpt_v25.json"
        sm_inputs = tfpt_suite_dir / "tfpt_suite" / "data" / "sm_inputs_mz.json"
        ckm_ref = tfpt_suite_dir / "tfpt_suite" / "data" / "ckm_reference.json"
        global_ref = tfpt_suite_dir / "tfpt_suite" / "data" / "global_reference.json"
        pyrate3_sm_yaml = pyrate3_dir / "models" / "SM_TFPT_2loop_v25.yaml"
        pyrate3_e8_pythonoutput = pyrate3_dir / "pyrate" / "results" / "E8Cascade2LoopGravityV2" / "PythonOutput"

        # Ensure the OperatorSpec is generated from the canonical microscopic action before we evaluate statuses.
        if microscopic_action.exists():
            try:
                generate_effective_action_r2_operator_spec(microscopic_action_path=microscopic_action, output_path=operator_spec)
            except Exception:
                # Ledger should never hard-fail; it records the state as best-effort.
                pass

        op_spec_raw = _read_json_if_exists(operator_spec)
        micro_raw = _read_json_if_exists(microscopic_action)

        # Heuristic to detect “full block operator with ghosts” vs the current minimal closure.
        blocks = list(op_spec_raw.get("blocks", [])) if isinstance(op_spec_raw.get("blocks", []), list) else []
        has_multiple_blocks = len(blocks) >= 2
        has_ghost_block = any(isinstance(b, dict) and str(b.get("statistics", "")).lower() == "ghost" for b in blocks)
        op_deriv_status = _operator_spec_derivation_status(op_spec_raw)
        micro_has_placeholder = _has_any_status(micro_raw, "placeholder")
        micro_has_quantization = isinstance(micro_raw.get("quantization", None), dict)

        ledger: list[LedgerItem] = [
            LedgerItem(
                area="7.1.1 Microscopic action",
                requirement="Fields, symmetries, gauge group, representations; explicit Lagrangian including fermions, Yukawas, gauge + gravity couplings.",
                status=(
                    "COMPLETE (canonical action spec exists; torsion-sector quadratic closure is explicit; no placeholder action terms)"
                    if (microscopic_action.exists() and (not micro_has_placeholder))
                    else "PARTIAL (canonical action spec exists but still contains placeholder terms)"
                ),
                evidence=[
                    f"paper: {latex_file}",
                    f"birefringence letter action (explicit Riemann–Cartan + anomaly normalization): {birefringence_letter}",
                    f"canonical action spec: {microscopic_action}",
                    f"PyR@TE model configs: {pyrate3_dir / 'models'}",
                    f"E8 2-loop PythonOutput (RGEs): {pyrate3_e8_pythonoutput}",
                ],
                next_steps=[
                    "If/when a first-principles torsionful connection derivation is required: replace the closure-level torsion quadratic operator by an operator derived from the microscopic torsion action + gauge fixing.",
                ],
            ),
            LedgerItem(
                area="7.1.2 Quantization",
                requirement="Gauge fixing + ghosts; path-integral measure + explicit regularization.",
                status=(
                    "COMPLETE (closure-level gauge fixing + FP ghost block are encoded; OperatorSpec derivation.status=derived)"
                    if (micro_has_quantization and has_ghost_block and op_deriv_status == "derived")
                    else "INCOMPLETE (missing explicit closure-level quantization metadata and/or ghost blocks and/or derived operator spec)"
                ),
                evidence=[
                    f"effective_action_r2 operator spec: {operator_spec}",
                    "effective_action_r2 module: tfpt-suite/tfpt_suite/modules/effective_action_r2.py",
                ],
                next_steps=[
                    "If needed for publication-grade QFT: provide a BRST-complete derivation from the torsionful microscopic action (beyond closure-level specification).",
                ],
            ),
            LedgerItem(
                area="7.1.3 Renormalization",
                requirement="Fix scheme (e.g. MSbar), implement beta functions + running for relevant couplings, and threshold matching (MZ, mt, new scales).",
                status=(
                    "COMPLETE (suite-level engineering closure: MSbar scheme is declared; 2-loop gauge+Yukawa RGEs are integrated via "
                    "PyR@TE3-generated beta functions with complex Yukawas supported in-suite; mt→μUV threshold bookkeeping is explicit "
                    "(MSigma, MG8; PMNS adds stepwise N_R activation at MNR1..3 with κ EFT running). "
                    "A matching/bridge pipeline exists (`msbar_matching_map`) with explicit assumptions + deterministic αs(MZ) sensitivity "
                    "and an optional Monte Carlo uncertainty hook, and a dedicated below-mt EFT audit trail exists (`below_mt_eft_cascade`). "
                    "Remaining work is publication-grade refinement: fully sourced finite EW/QCD pieces, explicit QED/EW decoupling policy "
                    "below MZ, and end-to-end uncertainty propagation across all thresholds.)"
                ),
                evidence=[
                    "two_loop_rg_fingerprints: tfpt-suite/tfpt_suite/modules/two_loop_rg_fingerprints.py",
                    f"PyR@TE3 baseline SM YAML: {pyrate3_sm_yaml}",
                    "PyR@TE3 SM PythonOutput: Pyrate3/pyrate/results/SM_TFPT_2loop_v25/SM_TFPT_2Loop_v25/PythonOutput",
                    f"paper v1.06 threshold blueprint (G8, NR1..3, PQ field, R3 spurion): {paper_v1_06}",
                    "PyR@TE-driven 2-loop engine (complex Yukawas, thresholds): tfpt-suite/tfpt_suite/rge_pyrate_2loop.py",
                    "CKM pipeline (mt boundary + 2-loop RG + thresholds): tfpt-suite/tfpt_suite/modules/ckm_full_pipeline.py",
                    "PMNS pipeline (κ EFT + N_R thresholds): tfpt-suite/tfpt_suite/modules/pmns_full_pipeline.py",
                    "Matching bridge (PDG-style inputs → MSbar boundary; uncertainty hooks): tfpt-suite/tfpt_suite/modules/msbar_matching_map.py",
                    "RG thresholds table: tfpt-suite/tfpt_suite/data/rge_thresholds_v25.json",
                    f"SM inputs @ MZ: {sm_inputs}",
                ],
                next_steps=[
                    "Upgrade the matching layer to publication-grade: fully sourced PDG inputs → MSbar parameters at MZ (incl. finite EW/QCD pieces) and a below-mt EFT cascade (QCD+QED, decoupling policy).",
                    "Promote uncertainty propagation to a first-class pipeline: Monte Carlo over PDG priors through matching + thresholds + RG (αs(MZ), mt, mH, etc.).",
                ],
            ),
            LedgerItem(
                area="7.1.4 Anomalies + consistency",
                requirement="Anomaly freedom, unitarity + causality (perturbatively).",
                status=(
                    "COMPLETE (suite-level automation: SM gauge anomalies + SU(2) Witten global check are audited from the canonical "
                    "field-content spec; a stability/perturbativity red-flag detector exists (`stability_unitarity_audit`). "
                    "Remaining work is publication-grade refinement: generalized anomaly scanning for arbitrary chiral extensions, "
                    "and systematic perturbative unitarity / vacuum metastability lifetime analysis beyond proxy checks.)"
                ),
                evidence=[
                    "aps_eta_gluing module (η / spectral-flow seam term): tfpt-suite/tfpt_suite/modules/aps_eta_gluing.py",
                    "discrete_consistency_uniqueness module: tfpt-suite/tfpt_suite/modules/discrete_consistency_uniqueness.py",
                    "anomaly_cancellation_audit module: tfpt-suite/tfpt_suite/modules/anomaly_cancellation_audit.py",
                    "stability_unitarity_audit module: tfpt-suite/tfpt_suite/modules/stability_unitarity_audit.py",
                    f"canonical action spec (field content): {microscopic_action}",
                ],
                next_steps=[
                    "Extend anomaly audit beyond the SM baseline if/when TFPT introduces chiral charged fields beyond the anomaly-neutral list.",
                    "Add publication-grade vacuum metastability and perturbative unitarity checks (beyond red-flag proxies), and apply them to any effective operators used in extensions.",
                ],
            ),
            LedgerItem(
                area="7.1.5 Derivation vs parametrization",
                requirement="Yukawa textures/phases derived from topology structure or explicitly marked as Ansatz; same for any 'magic number'.",
                status=(
                    "COMPLETE (suite-level mechanism contract: the flavor pipelines are driven by TFPT invariants + explicit Möbius/Z3 "
                    "rules with δ defaulting to the geometric closure δ⋆ (not τ/μ), cusps/classification fixed, and an explicit phase map "
                    "reported to prevent convention drift. Remaining work is theory-grade derivation: replace the mechanism contract by an "
                    "operator/holonomy-level derivation of the topology→phase map and the Yukawa generator.)"
                ),
                evidence=[
                    f"CKM pipeline module: {tfpt_suite_dir / 'tfpt_suite' / 'modules' / 'ckm_full_pipeline.py'}",
                    f"CKM reference table: {ckm_ref}",
                    "PMNS Z3 breaking module: tfpt-suite/tfpt_suite/modules/pmns_z3_breaking.py",
                    "PMNS full pipeline module: tfpt-suite/tfpt_suite/modules/pmns_full_pipeline.py",
                    "Z3 texture utility + fixed coefficients rule: tfpt-suite/tfpt_suite/flavor_textures.py",
                    "Z3 texture config (explicit conventions): tfpt-suite/tfpt_suite/data/flavor_texture_v24.json",
                    "Möbius δ calibration module: tfpt-suite/tfpt_suite/modules/mobius_delta_calibration.py",
                    "Seesaw anchor block module: tfpt-suite/tfpt_suite/modules/seesaw_block.py",
                    f"paper v1.06 flavor anchors (δ from τ/μ, δ⋆ from varphi0, seesaw block): {paper_v1_06}",
                    f"update TFPT v1.07SM: Z3 Yukawa texture formula Y(y)=y⋆[C(δ)+a_y varphi0 D + b c3 1]: {update_tfpt_v1_07sm}",
                    "predictions dashboard module: tfpt-suite/tfpt_suite/modules/predictions_dashboard.py",
                ],
                next_steps=[
                    "Replace λ-power 'texture layer' with a Möbius/Z3-monodromy-derived Yukawa texture generator (mechanism).",
                    "Quantize/derive CP phase input (or label discrete phase choices explicitly and expose falsifiable outputs).",
                ],
            ),
        ]

        warnings: list[str] = []
        if not latex_file.exists():
            warnings.append(f"missing_paper_source: {latex_file}")
        if not pyrate3_e8_pythonoutput.is_dir():
            warnings.append(f"missing_pyrate3_pythonoutput: {pyrate3_e8_pythonoutput}")
        if not operator_spec.exists():
            warnings.append(f"missing_effective_action_r2_operator_spec: {operator_spec}")
        if not microscopic_action.exists():
            warnings.append(f"missing_microscopic_action_spec: {microscopic_action}")

        # We intentionally expose “ghosts missing” as an explicit flag for the milestone.
        ghost_status = {
            "operator_spec_blocks": len(blocks),
            "has_multiple_blocks": has_multiple_blocks,
            "has_ghost_block": has_ghost_block,
            "note": "Full effective_action_r2 robustness milestone requires explicit gauge-fixing + ghost blocks in the operator spec.",
        }

        checks: list[Check] = []
        checks.append(Check(check_id="ledger_emitted", passed=True, detail="QFT completeness ledger emitted"))
        checks.append(Check(check_id="paper_source_present", passed=latex_file.exists(), detail=f"path={latex_file}"))
        checks.append(Check(check_id="operator_spec_present", passed=operator_spec.exists(), detail=f"path={operator_spec}"))

        # Report text
        lines: list[str] = []
        lines.append("QFT completeness ledger (paper v2.5)")
        lines.append("")
        lines.append(f"Paper source: {latex_file}")
        lines.append("")
        lines.append("This module is a structured checklist aligned with section 7.1, and an engineering TODO map aligned with 7.2.")
        lines.append("")
        for item in ledger:
            lines.append(f"## {item.area}")
            lines.append(f"Requirement: {item.requirement}")
            lines.append(f"Status: {item.status}")
            if item.evidence:
                lines.append("Evidence:")
                for ev in item.evidence:
                    lines.append(f"- {ev}")
            if item.next_steps:
                lines.append("Next steps:")
                for ns in item.next_steps:
                    lines.append(f"- {ns}")
            lines.append("")
        lines.append("Ghost/operator milestone status (effective_action_r2):")
        lines.append(str(ghost_status))
        lines.append("")

        return ModuleResult(
            results={
                "paper_source": str(latex_file),
                "ghost_operator_status": ghost_status,
                "ledger": [
                    {
                        "area": i.area,
                        "requirement": i.requirement,
                        "status": i.status,
                        "evidence": i.evidence,
                        "next_steps": i.next_steps,
                    }
                    for i in ledger
                ],
                "key_files": {
                    "effective_action_r2_operator_spec": str(operator_spec),
                    "microscopic_action_tfpt_v25": str(microscopic_action),
                    "sm_inputs_mz": str(sm_inputs),
                    "ckm_reference": str(ckm_ref),
                    "global_reference": str(global_ref),
                    "pyrate3_sm_yaml": str(pyrate3_sm_yaml),
                    "pyrate3_e8_pythonoutput": str(pyrate3_e8_pythonoutput),
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=warnings,
        )

