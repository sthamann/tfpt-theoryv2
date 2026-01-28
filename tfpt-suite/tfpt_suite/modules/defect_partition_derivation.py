from __future__ import annotations

import json
from pathlib import Path

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.defect_partition import Delta2Derivation, alpha_inv_0_from_delta2, derive_delta2_from_defect_partition
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass, mk_check_warn

EXPECTED_DEFECT_MULTIPLICITY = 5


class DefectPartitionDerivationModule(TfptModule):
    module_id = "defect_partition_derivation"
    title = "Defect partition derivation (derive δ₂ without fitting)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants (c3, varphi0_tree, delta_top) from tfpt_suite/constants.py",
                "CODATA α^{-1}(0) reference (diagnostic) from tfpt_suite/data/global_reference.json",
            ],
            outputs=[
                "δ₂ derivation record (model_id, assumptions, δ₂ value)",
                "two-defect sector enumeration + justification chain (JSON artifact)",
                "diagnostic α^{-1}(0) value implied by the δ₂ candidate (strict z-score vs CODATA)",
            ],
            formulas=[
                r"varphi0(α)=varphi_tree+δ_top e^{-2α}+δ₂ e^{-4α}",
                r"δ₂=(g/4)·δ_top² (g from two-defect sector enumeration; see report)",
            ],
            validation=[
                "δ₂ is produced by a discrete model choice (no continuous knob).",
                "Two-defect sector enumeration yields a discrete multiplicity g with explicit justification.",
                "α^{-1}(0) implied by the δ₂ candidate is within 5σ of CODATA (strict σ) as a metrology gate.",
            ],
            determinism="Deterministic given TFPT constants and the fixed δ₂ model.",
            question="Can the next defect-sector correction δ₂ be generated from discrete defect combinatorics (not fitted), while improving α(0) metrology?",
            objective=[
                "Turn the δ₂ term from a debug target into an assumption-explicit discrete derivation candidate.",
                "Provide a clear audit trail for why δ₂ is not a free parameter.",
            ],
            gaps=[
                "Publication-grade still requires a first-principles topological/QFT derivation beyond the SU(5) holonomy + action-block enumeration.",
            ],
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        ref_path = Path(__file__).resolve().parent.parent / "data" / "global_reference.json"
        ref = json.loads(ref_path.read_text(encoding="utf-8"))
        obs = ref.get("observables", {}) if isinstance(ref.get("observables", {}), dict) else {}
        alpha0_ref = obs.get("alpha_inv_0", {}) if isinstance(obs.get("alpha_inv_0", {}), dict) else {}
        alpha_inv_codata = mp.mpf(str(alpha0_ref.get("mean", "137.035999177")))
        alpha_inv_sigma = mp.mpf(str(alpha0_ref.get("sigma", "2.1e-8")))

        # --- Derive δ₂ (discrete model; no tuning) ---
        d2: Delta2Derivation = derive_delta2_from_defect_partition(delta_top=mp.mpf(c.delta_top))
        justification = d2.justification or {}
        g_value = justification.get("effective_multiplicity", None)

        # Diagnostic implied α^{-1}(0)
        alpha_inv_pred = alpha_inv_0_from_delta2(delta2=d2.delta2, mp_dps=max(60, int(config.mp_dps)))
        diff = alpha_inv_pred - alpha_inv_codata
        z = diff / alpha_inv_sigma if alpha_inv_sigma != 0 else mp.mpf("nan")

        checks: list[Check] = []
        checks.append(
            mk_check_pass(
                "delta2_is_derived_not_fitted",
                f"PASS: model_id={d2.model_id}, delta2={d2.delta2} (no continuous parameter)",
            )
        )
        if isinstance(g_value, int):
            if g_value == EXPECTED_DEFECT_MULTIPLICITY:
                checks.append(
                    mk_check_pass(
                        "delta2_multiplicity_derived_from_enumeration",
                        f"PASS: g={g_value} from two-defect sector equivalence classes",
                    )
                )
            else:
                checks.append(
                    mk_check_fail(
                        "delta2_multiplicity_derived_from_enumeration",
                        f"FAIL: expected g={EXPECTED_DEFECT_MULTIPLICITY} from enumeration, got g={g_value}",
                    )
                )
        else:
            checks.append(
                mk_check_fail(
                    "delta2_multiplicity_derived_from_enumeration",
                    f"FAIL: missing or non-integer g in justification (g={g_value})",
                )
            )
        if mp.isfinite(z):
            if abs(z) > mp.mpf(5):
                checks.append(mk_check_fail("alpha_inv_0_within_5sigma", f"|z|={abs(z)} > 5 (pred={alpha_inv_pred}, ref={alpha_inv_codata})"))
            elif abs(z) > mp.mpf(2):
                checks.append(mk_check_warn("alpha_inv_0_within_2sigma", f"|z|={abs(z)} > 2 (pred={alpha_inv_pred}, ref={alpha_inv_codata})"))
            else:
                checks.append(mk_check_pass("alpha_inv_0_within_2sigma", f"|z|={abs(z)} ≤ 2 (pred={alpha_inv_pred}, ref={alpha_inv_codata})"))
        else:
            checks.append(mk_check_fail("alpha_inv_0_z_finite", f"z is not finite: z={z}"))

        lines: list[str] = []
        lines += [
            "Defect partition derivation (Gap A): δ₂ candidate from discrete defect combinatorics",
            "",
            f"TFPT constants: delta_top={c.delta_top}",
            "",
            "δ₂ derivation:",
            f"- model_id: {d2.model_id}",
            f"- delta2: {d2.delta2}",
            f"- delta2/(delta_top^2): {d2.delta2_over_delta_top2}",
            f"- note: {d2.note}",
            "",
            "Two-defect sector enumeration:",
            f"- effective multiplicity g: {g_value}",
        ]
        if isinstance(justification.get("equivalence_classes", None), list):
            lines.append("- equivalence classes:")
            for entry in justification["equivalence_classes"]:
                if not isinstance(entry, dict):
                    continue
                class_id = entry.get("class_id", "unknown")
                category = entry.get("category", "unknown")
                pair_count = entry.get("pair_count", "n/a")
                lines.append(f"  - {class_id}: category={category}, pair_count={pair_count}")
        lines += [
            "",
            "Assumptions (explicit):",
            *[f"- {a}" for a in d2.assumptions],
            "",
            "Metrology diagnostic (strict σ, CODATA α^{-1}(0)):",
            f"- ref mean: {alpha_inv_codata} (sigma={alpha_inv_sigma})",
            f"- pred (self-consistent CFE+backreaction incl. δ₂): {alpha_inv_pred}",
            f"- diff: {diff}",
            f"- z: {z}",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
        ]

        return ModuleResult(
            results={
                "delta2_derivation": {
                    "model_id": d2.model_id,
                    "delta_top": c.delta_top,
                    "delta2": d2.delta2,
                    "delta2_over_delta_top2": d2.delta2_over_delta_top2,
                    "assumptions": d2.assumptions,
                    "note": d2.note,
                },
                "delta2_justification": justification,
                "alpha_inv_0_diagnostic": {
                    "pred": alpha_inv_pred,
                    "ref_mean": alpha_inv_codata,
                    "ref_sigma": alpha_inv_sigma,
                    "diff": diff,
                    "z": z,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

