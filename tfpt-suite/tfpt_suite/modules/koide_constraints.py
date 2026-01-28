from __future__ import annotations

import json
import math
from pathlib import Path

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_pass, mk_check_warn


class KoideConstraintsModule(TfptModule):
    module_id = "koide_constraints"
    title = "Koide constraints (charged leptons; diagnostic docking check)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=["lepton masses: tfpt_suite/data/lepton_masses_pdg.json"],
            outputs=[
                "Koide Q for charged leptons",
                "deviation from 2/3 (diagnostic; not a TFPT derivation yet)",
            ],
            formulas=[
                r"Q := (m_e+m_\mu+m_\tau)/(\sqrt{m_e}+\sqrt{m_\mu}+\sqrt{m_\tau})^2",
                r"Koide (empirical): Q \approx 2/3",
            ],
            validation=["Reports Q and deviation; check is a diagnostic PASS/WARN only (not a ToE gate)."],
            determinism="Deterministic given input table.",
            question="Do charged lepton pole masses satisfy the Koide relation (diagnostic)?",
            objective=[
                "Provide an explicit, citeable Koide diagnostic as a docking point for future topology/mass-ratio derivations.",
            ],
            gaps=[
                "This module does not derive Koide from TFPT; it only records the diagnostic value and deviation.",
            ],
        )

    def run(self, config) -> ModuleResult:
        lep_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        lep = json.loads(lep_path.read_text(encoding="utf-8"))
        me = float(lep["masses"]["electron"]["mean"])
        mm = float(lep["masses"]["muon"]["mean"])
        mt = float(lep["masses"]["tau"]["mean"])

        denom = (math.sqrt(me) + math.sqrt(mm) + math.sqrt(mt)) ** 2
        Q = float((me + mm + mt) / denom) if denom > 0 else float("nan")
        target = 2.0 / 3.0
        dev = float(Q - target) if math.isfinite(Q) else float("nan")

        checks: list[Check] = []
        # Diagnostic-only: Koide is empirical; keep this as PASS with deviation recorded.
        tol = 1e-4
        if math.isfinite(dev) and abs(dev) <= tol:
            checks.append(mk_check_pass("koide_leptons_close_to_2_over_3", f"Q={Q:.12g}, Q-2/3={dev:.3e} (|dev|≤{tol})"))
        else:
            checks.append(mk_check_warn("koide_leptons_close_to_2_over_3", f"Q={Q}, Q-2/3={dev} (tol={tol})"))

        lines: list[str] = []
        lines += [
            "Koide constraints (charged leptons; diagnostic)",
            "",
            f"input: {lep_path}",
            f"m_e={me:.12g} GeV, m_mu={mm:.12g} GeV, m_tau={mt:.12g} GeV",
            "",
            f"Q = {Q:.12g}",
            f"Q - 2/3 = {dev:.3e}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "input_file": str(lep_path),
                "masses_GeV": {"electron": me, "muon": mm, "tau": mt},
                "koide": {"Q": Q, "target": target, "deviation": dev},
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

