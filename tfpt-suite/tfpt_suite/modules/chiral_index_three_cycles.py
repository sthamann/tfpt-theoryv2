from __future__ import annotations

import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.module_base import (
    Check,
    ModuleResult,
    ModuleSpec,
    TfptModule,
    mk_check_info,
    mk_check_pass,
    mk_check_warn,
)


def _parse_fraction(value: object) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, (int,)):
        return Fraction(int(value), 1)
    if isinstance(value, float):
        # interpret floats as exact rationals only when safe-ish
        return Fraction(value).limit_denominator(10_000)
    s = str(value).strip()
    if "/" in s:
        n, d = s.split("/", 1)
        return Fraction(int(n.strip()), int(d.strip()))
    return Fraction(int(s), 1)


def _angle_mod_2pi(angle: Any) -> Any:
    two_pi = 2 * mp.pi
    # normalize into [0, 2π)
    a = angle % two_pi
    if a < 0:
        a += two_pi
    return a


@dataclass(frozen=True)
class _ChargePhase:
    label: str
    Y: Fraction
    cycle: str
    nu: int
    theta_rad: Any
    theta_over_pi: Any


class ChiralIndexThreeCyclesModule(TfptModule):
    module_id = "chiral_index_three_cycles"
    title = "Chiral index on orientable double cover (3 boundary cycles → 3 families; Wilson-line phase atoms)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "configuration: tfpt_suite/data/chiral_index_three_cycles.json (fluxes ν1,ν2,νT and example hypercharges)",
                "theory reference: paper_v1_06_01_09_2025.tex Appendix J (index theorem sketch + Wilson lines)",
            ],
            outputs=[
                "index IndD = ν1 + ν2 + νT (family number under the minimal flux choice)",
                "Gauss–Bonnet boundary-cycle bookkeeping: 2π + 2π + 2π = 6π (two physical boundaries + seam Γ)",
                "discrete U(1) Wilson-line phase atoms exp(i·2π·Y·νi) for representative SM hypercharges",
            ],
            formulas=[
                "Ind D_(M~) = (1/2π) ∫_{M~} F = ν1 + ν2 + νT  (Appendix J sketch)",
                "Gauss–Bonnet on orientable double cover with seam: K_total = 2π+2π+2π = 6π  (Appendix D/J sketch)",
                "U(1) Wilson line: ∮_Ci A = 2π νi ⇒ phase for charge Y is exp(i·Y·∮A)=exp(i·2π·Y·νi)",
            ],
            validation=[
                "index is integer and reproduces 3 families for the minimal flux choice (ν1,ν2,νT)=(1,1,1)",
                "phase-atom set is nontrivial for fractional hypercharges",
            ],
            determinism="Deterministic given the JSON config and mpmath precision.",
            question="Do three boundary cycles (two physical + seam Γ) naturally yield three chiral families via an index theorem, and do Wilson lines along these cycles provide a discrete phase source?",
            objective=[
                "Make the Appendix J mechanism audit-able: IndD=3 from the minimal integer flux choice.",
                "Provide a concrete, discrete 'phase atom' set from U(1) Wilson lines as a docking point for a future topology→CP-phase map.",
            ],
            what_was_done=[
                "Parsed a small explicit config (ν1,ν2,νT and example SM hypercharges).",
                "Computed IndD and the resulting family count for the default flux choice.",
                "Computed Wilson-line phase angles for representative SM hypercharges and enumerated unique phase atoms.",
            ],
            assumptions=[
                "Use the schematic index relation IndD = ν1+ν2+νT as stated in paper_v1_06_01_09_2025.tex Appendix J.",
                "Treat Wilson-line holonomies for U(1)_Y as exp(i·2π·Y·νi) with ∮A=2πνi.",
                "This module exposes Wilson-line phase atoms as a docking point. The discrete δ/δ_CP candidate map is implemented in `topology_phase_map`; a publication-grade operator-level derivation remains future work.",
            ],
            gaps=[
                "Publication-grade operator/holonomy-level derivation of topology→(δ, δ_CP) remains open (beyond the shipped finite candidate map).",
                "Enforce a unique, mechanism-derived selection rule for δ/δ_CP (beyond finite candidate enumeration + downstream χ² gates).",
            ],
            references=[
                "paper_v1_06_01_09_2025.tex Appendix D/J (6π seam + index sketch)",
                "update_tfptv1_07sm.tex (Z3 flavor architecture and phase conventions)",
            ],
            maturity="mechanism scaffold (topology→family number; phase docking point; not yet a full flavor derivation)",
        )

    def run(self, config) -> ModuleResult:
        data_dir = Path(__file__).resolve().parent.parent / "data"
        cfg_path = data_dir / "chiral_index_three_cycles.json"
        raw = json.loads(cfg_path.read_text(encoding="utf-8"))

        fluxes = raw.get("default_fluxes", {}) if isinstance(raw.get("default_fluxes", {}), dict) else {}
        nu1 = int(fluxes.get("nu1", 1))
        nu2 = int(fluxes.get("nu2", 1))
        nuT = int(fluxes.get("nuT", 1))
        indD = int(nu1 + nu2 + nuT)

        charges_in = raw.get("example_charges_one_generation", [])
        charges: list[tuple[str, Fraction]] = []
        for e in charges_in:
            if not isinstance(e, dict):
                continue
            label = str(e.get("label", "")).strip() or "unnamed"
            Y = _parse_fraction(e.get("Y", "0"))
            charges.append((label, Y))

        phase_rows: list[_ChargePhase] = []
        two_pi = 2 * mp.pi

        for label, Y in charges:
            for cycle, nu in [("C1", nu1), ("C2", nu2), ("CT", nuT)]:
                theta = two_pi * mp.mpf(Y.numerator) * mp.mpf(nu) / mp.mpf(Y.denominator)
                theta_mod = _angle_mod_2pi(theta)
                phase_rows.append(
                    _ChargePhase(
                        label=label,
                        Y=Y,
                        cycle=cycle,
                        nu=nu,
                        theta_rad=theta_mod,
                        theta_over_pi=theta_mod / mp.pi,
                    )
                )

        # Unique “phase atoms” (angles modulo 2π), represented as θ/π with limited denom for readability.
        atoms: list[Fraction] = []
        for row in phase_rows:
            # Convert θ/π to a rational when possible (we generated it from rationals anyway).
            # θ/π = 2*Y*nu mod 2
            atoms.append(Fraction(2 * row.Y * row.nu).limit_denominator(12))

        # Reduce modulo 2 (since angles are modulo 2π).
        atoms_mod2 = [Fraction(a.numerator % (2 * a.denominator), a.denominator).limit_denominator(12) for a in atoms]
        atoms_unique = sorted(set(atoms_mod2), key=lambda x: (float(x), x.denominator, x.numerator))

        has_nontrivial = any(a != 0 for a in atoms_unique)
        minimal_flux_choice = (nu1, nu2, nuT) == (1, 1, 1)

        checks: list[Check] = []
        checks.append(mk_check_pass("index_is_integer", f"IndD={indD} from nu=(nu1,nu2,nuT)=({nu1},{nu2},{nuT})"))
        if minimal_flux_choice:
            checks.append(mk_check_pass("minimal_flux_gives_three_families", "nu=(1,1,1) => IndD=3 families"))
        else:
            checks.append(mk_check_info("nonminimal_flux_choice_used", f"nu=({nu1},{nu2},{nuT}) => IndD={indD}"))
        if has_nontrivial:
            checks.append(
                mk_check_pass(
                    "wilson_line_phase_atoms_nontrivial",
                    f"distinct theta/pi atoms (mod 2): {[str(a) for a in atoms_unique]}",
                )
            )
        else:
            checks.append(mk_check_warn("wilson_line_phase_atoms_trivial", "all phase atoms are 0 (unexpected for fractional Y)"))

        # Explicit gap marker (prevents “everything green” optics in physics mode).
        checks.append(
            mk_check_info(
                "topology_to_cp_phase_map",
                "Wilson-line phase atoms are provided here; the shipped finite δ/δ_CP candidate map is implemented in `topology_phase_map` (operator-level derivation remains future work).",
            )
        )

        report_lines: list[str] = [
            "Chiral index on orientable double cover (three cycles)",
            f"config: {cfg_path}",
            f"mp.dps = {mp.dps}",
            "",
            "Flux choice:",
            f"- (nu1,nu2,nuT)=({nu1},{nu2},{nuT})",
            f"- IndD = nu1 + nu2 + nuT = {indD}",
            "",
            "Gauss–Bonnet / seam bookkeeping (sketch):",
            "- two physical boundary cycles + one seam Γ contribute 2π each ⇒ total 6π (used in φ_tree normalization in the main theory)",
            "",
            "Wilson-line phase atoms (U(1)_Y holonomy model):",
            "- phase for charge Y on cycle Ci: exp(i·2π·Y·νi)",
            f"- unique atoms (θ/π mod 2, limited denom): {[str(a) for a in atoms_unique]}",
            "",
            "Per-charge / per-cycle angles (θ/π):",
        ]
        for row in phase_rows:
            report_lines.append(f"- {row.label:24s} Y={row.Y!s:>6s}  {row.cycle} (nu={row.nu}): theta/pi={row.theta_over_pi}")
        report_lines += [
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
            "",
            "Notes:",
            "- This module is a docking point: it makes the Appendix J mechanism explicit and machine-readable.",
            "- The next step is to connect these discrete holonomies to the Möbius/Z3 deformation parameter δ and the CKM/PMNS CP phases (operator/holonomy-level).",
            "",
        ]

        results: dict[str, Any] = {
            "config_file": str(cfg_path),
            "fluxes": {"nu1": nu1, "nu2": nu2, "nuT": nuT},
            "index": {"IndD": indD, "interpretation": "family number under the Appendix J sketch"},
            "gauss_bonnet": {"cycles": ["C1", "C2", "CT"], "total_boundary_curvature": "6π (2π+2π+2π)"},
            "wilson_line_phase_atoms": {
                "atoms_theta_over_pi_mod2": [str(a) for a in atoms_unique],
                "u1_convention": raw.get("u1_convention", {}),
                "example_charges_one_generation": raw.get("example_charges_one_generation", []),
                "per_charge_cycles": [
                    {
                        "label": row.label,
                        "Y": str(row.Y),
                        "cycle": row.cycle,
                        "nu": row.nu,
                        "theta_over_pi": str(row.theta_over_pi),
                    }
                    for row in phase_rows
                ],
            },
        }

        return ModuleResult(results=results, checks=checks, report="\n".join(report_lines), warnings=[])

