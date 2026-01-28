from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


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


def _ledger_entry(*, value: float, status: str, note: str, source: str | None = None) -> dict[str, object]:
    return {
        "value": float(value),
        "status": str(status),
        "note": str(note),
        "source": str(source) if source else None,
    }


class MassSpectrumMinimalModule(TfptModule):
    module_id = "mass_spectrum_minimal"
    title = "Mass spectrum minimal (ledger: derived vs input)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "SM inputs at MZ: tfpt_suite/data/sm_inputs_mz.json",
                "Lepton masses: tfpt_suite/data/lepton_masses_pdg.json",
            ],
            outputs=[
                "mass ledger with explicit status tags (input/derived/placeholder)",
                "minimal v and Higgs λ bookkeeping tags",
            ],
            formulas=[
                "This module is intentionally a ledger, not a fitter: it records what is treated as input vs what is derived elsewhere.",
            ],
            validation=[
                "Every numeric entry carries an explicit status tag.",
                "No silent 'magic' values: placeholders must be marked as such.",
            ],
            determinism="Deterministic given inputs.",
            question="Which masses/scales are inputs vs derived quantities in the current suite?",
            objective=[
                "Provide a minimal ToE scope map without pretending missing derivations are solved.",
                "Make it reviewer-proof which quantities are 'input' vs 'derived'.",
            ],
        )

    def run(self, config) -> ModuleResult:
        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"
        lep_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        lep_raw = json.loads(lep_path.read_text(encoding="utf-8"))

        v_ev = float(sm_raw.get("v_ev_GeV", 246.0))
        mH = float(sm_raw.get("mH_GeV", 125.25))
        lam_tree = float((mH * mH) / (2.0 * v_ev * v_ev))

        masses: dict[str, dict[str, object]] = {}
        masses["SM_inputs_mz"] = {
            "MZ_GeV": _ledger_entry(value=float(sm_raw.get("mu_GeV", 91.1876)), status="input", note="declared boundary scale μ=MZ", source=_relpath(sm_path)),
            "mW_GeV": _ledger_entry(value=float(sm_raw.get("mW_GeV", 80.379)), status="input", note="PDG-ish pole mass (threshold bookkeeping)", source=_relpath(sm_path)),
            "mH_GeV": _ledger_entry(value=mH, status="input", note="PDG-ish Higgs pole mass (used for λ_tree)", source=_relpath(sm_path)),
            "mc_GeV": _ledger_entry(value=float(sm_raw.get("mc_GeV", 1.27)), status="input", note="threshold proxy (scheme-dependent)", source=_relpath(sm_path)),
            "mb_GeV": _ledger_entry(value=float(sm_raw.get("mb_GeV", 4.18)), status="input", note="interpreted as mb(mb) in GeV (scheme-dependent)", source=_relpath(sm_path)),
            "mt_GeV": _ledger_entry(value=float(sm_raw.get("mt_GeV", 172.76)), status="input", note="top threshold proxy (pole vs MSbar is scheme-dependent)", source=_relpath(sm_path)),
        }
        masses["leptons_pdg"] = {
            "electron_GeV": _ledger_entry(
                value=float(lep_raw["masses"]["electron"]["mean"]),
                status="input",
                note="PDG lepton pole mass (used for QED thresholds / Yukawa proxies)",
                source=_relpath(lep_path),
            ),
            "muon_GeV": _ledger_entry(
                value=float(lep_raw["masses"]["muon"]["mean"]),
                status="input",
                note="PDG lepton pole mass (used for Möbius δ calibration)",
                source=_relpath(lep_path),
            ),
            "tau_GeV": _ledger_entry(
                value=float(lep_raw["masses"]["tau"]["mean"]),
                status="input",
                note="PDG lepton pole mass (used for Möbius δ calibration / Yukawa proxy)",
                source=_relpath(lep_path),
            ),
        }
        masses["ew_derived"] = {
            "v_ev_GeV": _ledger_entry(value=v_ev, status="input", note="suite convention; can be derived from G_F in a publication-grade pass", source=_relpath(sm_path)),
            "lambda_tree_from_mH_v": _ledger_entry(value=lam_tree, status="derived", note="λ_tree = mH^2/(2 v^2) (explicit approximation)", source=_relpath(sm_path)),
        }

        # Placeholder slots for a publication-grade ToE spectrum derivation.
        masses["placeholders"] = {
            "neutrino_mass_scale_eV": _ledger_entry(
                value=0.05,
                status="placeholder",
                note="placeholder (order-of-magnitude): requires κ/yN/MR reconstruction and a unique selection rule",
                source=None,
            ),
            "proton_mass_GeV": _ledger_entry(
                value=0.938272,
                status="placeholder",
                note="placeholder (PDG): deriving hadron masses from TFPT/QCD completeness is tracked as a separate gap",
                source=None,
            ),
        }

        # Validation: every ledger entry must carry a valid status tag.
        allowed_status = {"input", "derived", "placeholder"}
        bad: list[str] = []
        for block_name, block in masses.items():
            for key, entry in dict(block).items():
                st = str(entry.get("status", ""))
                if st not in allowed_status:
                    bad.append(f"{block_name}.{key} (status={st!r})")

        checks = [
            Check(
                check_id="ledger_status_tags_present",
                passed=bool(len(bad) == 0),
                detail="PASS" if len(bad) == 0 else f"Invalid/missing status tags: {bad[:10]}",
            )
        ]

        lines: list[str] = []
        lines += [
            "Mass spectrum minimal (ledger: derived vs input)",
            "",
            f"SM inputs: {_relpath(sm_path)} (sha256={_sha256_file(sm_path)})",
            f"lepton masses: {_relpath(lep_path)} (sha256={_sha256_file(lep_path)})",
            "",
            "Ledger (value + status):",
        ]
        for block_name in ["SM_inputs_mz", "leptons_pdg", "ew_derived", "placeholders"]:
            lines.append(f"- {block_name}:")
            for key, entry in masses[block_name].items():
                lines.append(f"  - {key}: {entry.get('value')}  [{entry.get('status')}]  ({entry.get('note')})")
        lines += [
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "inputs": {
                    "sm_inputs_file": _relpath(sm_path),
                    "sm_inputs_sha256": _sha256_file(sm_path),
                    "lepton_masses_file": _relpath(lep_path),
                    "lepton_masses_sha256": _sha256_file(lep_path),
                },
                "ledger": masses,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

