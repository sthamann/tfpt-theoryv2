from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


@dataclass(frozen=True)
class Fermion:
    name: str
    generations: int
    su3_rep: str
    su3_dim: int
    su2_rep: str
    su2_dim: int
    Y: float


def _load_field_content(path: Path) -> tuple[list[Fermion], dict[str, Any]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    # Canonical source: microscopic action spec (fields.sm.*).
    sm = raw.get("fields", {}).get("sm", {})
    fermions_raw = []
    fermions_raw += list(sm.get("fermions_left_handed_weyl", []))
    fermions_raw += list(sm.get("optional_anomaly_neutral_extensions", []))
    out: list[Fermion] = []
    for f in fermions_raw:
        if not isinstance(f, dict):
            continue
        reps = f.get("reps", {})
        dims = f.get("dims", {})
        out.append(
            Fermion(
                name=str(f["name"]),
                generations=int(f.get("generations", 1)),
                su3_rep=str(reps.get("SU3c", "singlet")),
                su3_dim=int(dims.get("SU3c", 1)),
                su2_rep=str(reps.get("SU2L", "singlet")),
                su2_dim=int(dims.get("SU2L", 1)),
                Y=float(reps.get("U1Y", 0.0)),
            )
        )
    return out, raw


def _rep_constants(field_content_raw: dict[str, Any]) -> dict[str, dict[str, float]]:
    rep = field_content_raw.get("conventions", {}).get("representation_data", {})
    su3_T2 = dict(rep.get("SU3_T2", {})) or {"fund": 0.5, "anti": 0.5, "singlet": 0.0, "adj": 3.0}
    su3_A3 = dict(rep.get("SU3_A3", {})) or {"fund": 1.0, "anti": -1.0, "singlet": 0.0, "adj": 0.0}
    su2_T2 = dict(rep.get("SU2_T2", {})) or {"doublet": 0.5, "triplet": 2.0, "singlet": 0.0}
    return {
        "SU3_T2": {k: float(v) for k, v in su3_T2.items()},
        "SU3_A3": {k: float(v) for k, v in su3_A3.items()},
        "SU2_T2": {k: float(v) for k, v in su2_T2.items()},
    }


def _sum_anomalies(fermions: list[Fermion], rep: dict[str, dict[str, float]]) -> dict[str, float]:
    su3_T2 = rep["SU3_T2"]
    su3_A3 = rep["SU3_A3"]
    su2_T2 = rep["SU2_T2"]

    A_U1_3 = 0.0
    A_SU2_2_U1 = 0.0
    A_SU3_2_U1 = 0.0
    A_grav_2_U1 = 0.0
    A_SU3_3 = 0.0
    n_su2_doublets = 0

    for f in fermions:
        mult = int(f.generations)
        d3 = int(f.su3_dim)
        d2 = int(f.su2_dim)
        Y = float(f.Y)

        A_U1_3 += mult * d3 * d2 * (Y**3)
        A_grav_2_U1 += mult * d3 * d2 * Y

        # SU(2)^2-U(1): sum d(SU3)*T2(SU2)*Y
        t2_su2 = float(su2_T2.get(f.su2_rep, 0.0))
        A_SU2_2_U1 += mult * d3 * t2_su2 * Y

        # SU(3)^2-U(1): sum d(SU2)*T2(SU3)*Y
        t2_su3 = float(su3_T2.get(f.su3_rep, 0.0))
        A_SU3_2_U1 += mult * d2 * t2_su3 * Y

        # SU(3)^3: sum d(SU2)*A3(SU3)
        a3_su3 = float(su3_A3.get(f.su3_rep, 0.0))
        A_SU3_3 += mult * d2 * a3_su3

        # Witten SU(2) global anomaly counts SU(2) doublets (including color multiplicity)
        if f.su2_rep == "doublet":
            n_su2_doublets += mult * d3

    return {
        "U1Y_cubed": A_U1_3,
        "SU2_sq_U1Y": A_SU2_2_U1,
        "SU3_sq_U1Y": A_SU3_2_U1,
        "grav_sq_U1Y": A_grav_2_U1,
        "SU3_cubed": A_SU3_3,
        "su2_doublets_count": float(n_su2_doublets),
    }


class AnomalyCancellationAuditModule(TfptModule):
    module_id = "anomaly_cancellation_audit"
    title = "Anomaly cancellation audit (SM + TFPT anomaly-neutral extensions)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "Microscopic action file: tfpt_suite/data/microscopic_action_tfpt_v25.json (canonical SM left-handed Weyl fermions, SM hypercharge normalization)",
            ],
            outputs=[
                "Gauge anomaly coefficients (U(1)^3, SU(2)^2-U(1), SU(3)^2-U(1), grav^2-U(1), SU(3)^3)",
                "Global SU(2) Witten check (even number of doublets)",
            ],
            formulas=[
                "U(1)^3: Σ d3 d2 Y^3",
                "SU(2)^2-U(1): Σ d3 T2(SU2) Y",
                "SU(3)^2-U(1): Σ d2 T2(SU3) Y",
                "grav^2-U(1): Σ d3 d2 Y",
                "SU(3)^3: Σ d2 A3(SU3)",
            ],
            validation=[
                "SM anomalies cancel exactly (up to floating error from decimals used for Y).",
                "SU(2) Witten anomaly: even number of doublets (counting color + generations).",
            ],
            determinism="Deterministic given the fixed field-content JSON.",
        )

    def run(self, config) -> ModuleResult:
        path = Path(__file__).resolve().parent.parent / "data" / "microscopic_action_tfpt_v25.json"
        fermions, raw = _load_field_content(path)
        rep = _rep_constants(raw)
        sums = _sum_anomalies(fermions, rep)

        # Use tight tolerance; Ys are decimals, so expect tiny float error.
        tol = 1e-12
        checks: list[Check] = []
        checks.append(Check("u1_cubed_cancels", abs(sums["U1Y_cubed"]) < tol, f"U(1)^3 = {sums['U1Y_cubed']}"))
        checks.append(Check("su2_sq_u1_cancels", abs(sums["SU2_sq_U1Y"]) < tol, f"SU(2)^2-U(1) = {sums['SU2_sq_U1Y']}"))
        checks.append(Check("su3_sq_u1_cancels", abs(sums["SU3_sq_U1Y"]) < tol, f"SU(3)^2-U(1) = {sums['SU3_sq_U1Y']}"))
        checks.append(Check("grav_sq_u1_cancels", abs(sums["grav_sq_U1Y"]) < tol, f"grav^2-U(1) = {sums['grav_sq_U1Y']}"))
        checks.append(Check("su3_cubed_cancels", abs(sums["SU3_cubed"]) < tol, f"SU(3)^3 = {sums['SU3_cubed']}"))

        # Witten: even number of SU(2) doublets (counting color/generations)
        n_doublets = int(round(sums["su2_doublets_count"]))
        checks.append(
            Check(
                "su2_witten_even_doublets",
                (n_doublets % 2) == 0,
                f"SU(2) doublets count (color×gen) = {n_doublets}",
            )
        )

        lines: list[str] = []
        lines.append("Anomaly cancellation audit (SM + TFPT anomaly-neutral extensions)")
        lines.append("")
        lines.append(f"Field content file: {path}")
        lines.append("")
        lines.append("Anomaly sums (SM hypercharge normalization):")
        lines.append(f"- U(1)_Y^3            = {sums['U1Y_cubed']}")
        lines.append(f"- SU(2)_L^2-U(1)_Y    = {sums['SU2_sq_U1Y']}")
        lines.append(f"- SU(3)_c^2-U(1)_Y    = {sums['SU3_sq_U1Y']}")
        lines.append(f"- grav^2-U(1)_Y       = {sums['grav_sq_U1Y']}")
        lines.append(f"- SU(3)_c^3           = {sums['SU3_cubed']}")
        lines.append("")
        lines.append(f"SU(2) Witten global check: #doublets = {n_doublets} (even required) -> {'PASS' if (n_doublets%2)==0 else 'FAIL'}")
        lines.append("")
        lines.append("Checks:")
        for chk in checks:
            lines.append(f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})")
        lines.append("")

        return ModuleResult(
            results={
                "field_content_file": str(path),
                "anomaly_sums": sums,
                "tolerance": tol,
                "n_fermions_listed": len(fermions),
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

