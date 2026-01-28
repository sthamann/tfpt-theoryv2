from __future__ import annotations

import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass, mk_check_warn


def _workspace_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _parse_fraction(x: object) -> Fraction:
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x, 1)
    if isinstance(x, str):
        s = x.strip()
        if "/" in s:
            a, b = s.split("/", 1)
            return Fraction(int(a.strip()), int(b.strip()))
        return Fraction(int(s), 1)
    raise TypeError(f"Unsupported fraction input: {type(x)}")


@dataclass(frozen=True)
class PhaseAtom:
    theta_over_pi_mod2: str
    theta_rad: float
    cycle_id: str


@dataclass(frozen=True)
class HolonomyClass:
    class_id: str
    cycle_id: str
    nu: int
    atoms_theta_over_pi_mod2: list[str]


class TopologyPhaseMapModule(TfptModule):
    module_id = "topology_phase_map"
    title = "Topology phase map (Wilson-line atoms → discrete δ / δ_CP candidates)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "Wilson-line / flux config: tfpt_suite/data/chiral_index_three_cycles.json",
                "lepton masses: tfpt_suite/data/lepton_masses_pdg.json (for δ_M anchor)",
                "TFPT invariants: tfpt_suite/constants.py (δ⋆ anchor from varphi0)",
            ],
            outputs=[
                "finite candidate set for δ (texture deformation) and δ_CP (CP phase) values",
                "explicit phase atoms (θ/π mod 2) derived from the U(1) holonomy model",
            ],
            formulas=[
                r"Wilson-line phase: exp(i·2π·Y·ν_i) ⇒ θ/π = 2·Y·ν_i (mod 2)",
                r"δ_M = (sqrt(m_tau/m_mu)-1)/(sqrt(m_tau/m_mu)+1) (Möbius anchor)",
                r"δ⋆ = 3/5 + varphi0/6 (suite/paper anchor)",
            ],
            validation=[
                "phase_map_is_discrete: no continuous tuning parameters appear in the candidate set",
                "delta_star_anchor_present: δ⋆ must be included as a baseline δ candidate",
                "gauge_relabeling_invariant: phase atoms are stored mod 2 (U(1) relabeling)",
                "cycle_permutation_covariant: candidate set invariant under cycle relabeling",
                "complex_conjugation_consistent: δ_CP branches include conjugates",
            ],
            determinism="Deterministic given inputs (finite enumeration).",
            question="Given discrete Wilson-line phase atoms, what discrete δ / δ_CP candidates arise under a minimal mapping policy?",
            objective=[
                "Provide the explicit discrete candidate set needed by downstream joint flavor scans (CKM+PMNS).",
                "Keep this module assumption-explicit: it is a docking map, not yet a full operator derivation.",
            ],
            gaps=[
                "A publication-grade topology→phase derivation must connect these atoms to the actual Yukawa operator structure (holonomy/APS η input → δ, δ_CP).",
            ],
        )

    def run(self, config) -> ModuleResult:
        cfg_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "chiral_index_three_cycles.json"
        lep_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "lepton_masses_pdg.json"
        raw = json.loads(cfg_path.read_text(encoding="utf-8"))
        lep = json.loads(lep_path.read_text(encoding="utf-8"))

        fluxes = raw.get("default_fluxes", {}) if isinstance(raw.get("default_fluxes", {}), dict) else {}
        nu1 = int(fluxes.get("nu1", 1))
        nu2 = int(fluxes.get("nu2", 1))
        nuT = int(fluxes.get("nuT", 1))

        charges_in = raw.get("example_charges_one_generation", [])
        charges: list[Fraction] = []
        for e in charges_in:
            if not isinstance(e, dict):
                continue
            charges.append(_parse_fraction(e.get("Y", "0")))

        # Compute phase atoms per cycle (θ/π mod 2) from Wilson-line holonomies.
        def atoms_for_nu(nu: int) -> list[Fraction]:
            return [
                Fraction((2 * Y * nu)).limit_denominator(12)
                for Y in charges
            ]

        cycle_nus = {"C1": nu1, "C2": nu2, "CT": nuT}
        atoms_by_cycle: dict[str, list[Fraction]] = {}
        for cycle_id, nu in cycle_nus.items():
            atoms_raw = atoms_for_nu(nu)
            atoms_by_cycle[cycle_id] = [
                Fraction(a.numerator % (2 * a.denominator), a.denominator).limit_denominator(12) for a in atoms_raw
            ]

        atoms_union = [a for atoms in atoms_by_cycle.values() for a in atoms]
        atoms_unique = sorted(set(atoms_union), key=lambda x: (float(x), x.denominator, x.numerator))

        phase_atoms: list[PhaseAtom] = []
        for cycle_id, atoms in atoms_by_cycle.items():
            for a in sorted(set(atoms), key=lambda x: (float(x), x.denominator, x.numerator)):
                theta = float(math.pi * float(a))  # θ = π (θ/π)
                phase_atoms.append(PhaseAtom(theta_over_pi_mod2=str(a), theta_rad=theta, cycle_id=cycle_id))

        holonomy_classes: list[HolonomyClass] = []
        for cycle_id, atoms in atoms_by_cycle.items():
            atoms_sorted = sorted(set(atoms), key=lambda x: (float(x), x.denominator, x.numerator))
            class_id = f"{cycle_id}:{','.join(str(a) for a in atoms_sorted)}"
            holonomy_classes.append(
                HolonomyClass(
                    class_id=class_id,
                    cycle_id=cycle_id,
                    nu=cycle_nus[cycle_id],
                    atoms_theta_over_pi_mod2=[str(a) for a in atoms_sorted],
                )
            )

        # δ anchors (discrete candidate set): δ⋆ (TFPT) and δ_M (Möbius-from-leptons).
        c = TfptConstants.compute()
        delta_star = mp.mpf(c.delta_star)
        mtau = mp.mpf(str(lep["masses"]["tau"]["mean"]))
        mmu = mp.mpf(str(lep["masses"]["muon"]["mean"]))
        R = mp.sqrt(mtau / mmu)
        delta_M = (R - 1) / (R + 1)

        delta_candidates: list[dict[str, object]] = [
            {"label": "delta_star", "delta": float(delta_star), "source": "TfptConstants.delta_star"},
            {"label": "delta_M_from_tau_mu", "delta": float(delta_M), "source": "mobius_delta_calibration formula"},
        ]

        # δ_CP candidates from phase atoms (include conjugate branch 2π-θ).
        delta_cp_candidates: list[dict[str, object]] = []
        for a in atoms_unique:
            theta = float(math.pi * float(a))
            delta_cp_candidates.append({"theta_over_pi_mod2": str(a), "delta_cp_rad": float(theta), "branch": "theta"})
            if theta != 0.0:
                delta_cp_candidates.append(
                    {"theta_over_pi_mod2": str(a), "delta_cp_rad": float((2.0 * math.pi) - theta), "branch": "conjugate"}
                )

        # Explicit holonomy-class → (δ, δ_CP) map.
        holonomy_map: list[dict[str, object]] = []
        for cls in holonomy_classes:
            cls_atoms = [Fraction(x) for x in cls.atoms_theta_over_pi_mod2]
            cls_delta_cp: list[dict[str, object]] = []
            for a in cls_atoms:
                theta = float(math.pi * float(a))
                cls_delta_cp.append({"theta_over_pi_mod2": str(a), "delta_cp_rad": float(theta), "branch": "theta"})
                if theta != 0.0:
                    cls_delta_cp.append(
                        {"theta_over_pi_mod2": str(a), "delta_cp_rad": float((2.0 * math.pi) - theta), "branch": "conjugate"}
                    )
            holonomy_map.append(
                {
                    "class_id": cls.class_id,
                    "cycle_id": cls.cycle_id,
                    "nu": cls.nu,
                    "delta_candidates": delta_candidates,
                    "delta_cp_candidates": cls_delta_cp,
                }
            )

        # Pair candidates (finite, explicit; downstream scans can apply additional penalties/filters).
        pairs: list[dict[str, object]] = []
        for d in delta_candidates:
            for ph in delta_cp_candidates:
                pairs.append(
                    {
                        "delta_label": str(d["label"]),
                        "delta": float(d["delta"]),
                        "delta_cp_rad": float(ph["delta_cp_rad"]),
                        "theta_over_pi_mod2": str(ph["theta_over_pi_mod2"]),
                        "branch": str(ph["branch"]),
                    }
                )

        checks: list[Check] = []
        checks.append(mk_check_pass("phase_map_is_discrete", f"|atoms|={len(atoms_unique)}, |pairs|={len(pairs)} (finite enumeration)"))
        checks.append(
            mk_check_pass("delta_star_anchor_present", f"delta_star={float(delta_star)}")
            if any(d.get("label") == "delta_star" for d in delta_candidates)
            else mk_check_fail("delta_star_anchor_present", "delta_star missing from delta candidates")
        )
        # Gauge relabeling invariance: all atoms are stored mod 2 (U(1) reparam).
        relabel_ok = all(
            a == Fraction(a.numerator % (2 * a.denominator), a.denominator).limit_denominator(12) for a in atoms_unique
        )
        checks.append(
            mk_check_pass("gauge_relabeling_invariant", "atoms stored mod 2 (U(1) relabeling)")
            if relabel_ok
            else mk_check_fail("gauge_relabeling_invariant", "found atoms outside mod-2 canonicalization")
        )
        # Cycle permutation covariance: union of atoms invariant under relabeling cycles.
        perm_ok = True
        base = set(atoms_unique)
        for perm in ((nu1, nu2, nuT), (nu1, nuT, nu2), (nu2, nu1, nuT), (nu2, nuT, nu1), (nuT, nu1, nu2), (nuT, nu2, nu1)):
            atoms_perm: list[Fraction] = []
            for nu in perm:
                atoms_perm.extend(atoms_for_nu(nu))
            atoms_perm_mod2 = {
                Fraction(a.numerator % (2 * a.denominator), a.denominator).limit_denominator(12) for a in atoms_perm
            }
            if atoms_perm_mod2 != base:
                perm_ok = False
                break
        checks.append(
            mk_check_pass("cycle_permutation_covariant", "atom set invariant under cycle relabeling")
            if perm_ok
            else mk_check_fail("cycle_permutation_covariant", "atom set changes under cycle relabeling")
        )
        # Complex conjugation consistency: each nonzero atom yields both branches in δ_CP candidates.
        branches_by_atom: dict[str, set[str]] = {}
        for entry in delta_cp_candidates:
            key = str(entry.get("theta_over_pi_mod2"))
            branch = str(entry.get("branch"))
            branches_by_atom.setdefault(key, set()).add(branch)
        conj_ok = True
        for a in atoms_unique:
            if a == 0:
                continue
            branches = branches_by_atom.get(str(a), set())
            if not {"theta", "conjugate"}.issubset(branches):
                conj_ok = False
                break
        checks.append(
            mk_check_pass("complex_conjugation_consistent", "conjugate branches present for all atoms")
            if conj_ok
            else mk_check_fail("complex_conjugation_consistent", "missing conjugate atom for some phases")
        )
        if not atoms_unique:
            checks.append(mk_check_warn("phase_atoms_empty", "no phase atoms produced (check config)"))

        lines: list[str] = []
        lines += [
            "Topology phase map (docking module): Wilson-line atoms -> discrete δ / δ_CP candidates",
            "",
            f"config: {cfg_path}",
            f"default fluxes: (nu1,nu2,nuT)=({nu1},{nu2},{nuT})",
            "",
            "Phase atoms by cycle (θ/π mod 2):",
            *[f"- {cid}: {[str(a) for a in atoms_by_cycle[cid]]}" for cid in ("C1", "C2", "CT")],
            "",
            "Holonomy classes:",
            *[f"- {cls.class_id}" for cls in holonomy_classes],
            "",
            "δ anchors:",
            f"- delta_star = {float(delta_star)}",
            f"- delta_M(from tau/mu) = {float(delta_M)}",
            "",
            "Candidate set summary:",
            f"- delta_candidates={len(delta_candidates)}, delta_cp_candidates={len(delta_cp_candidates)}, pairs={len(pairs)}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "config_file": str(cfg_path),
                "fluxes": {"nu1": nu1, "nu2": nu2, "nuT": nuT},
                "phase_atoms": [a.__dict__ for a in phase_atoms],
                "holonomy_classes": [cls.__dict__ for cls in holonomy_classes],
                "holonomy_map": holonomy_map,
                "delta_candidates": delta_candidates,
                "delta_cp_candidates": delta_cp_candidates,
                "pairs": pairs,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

