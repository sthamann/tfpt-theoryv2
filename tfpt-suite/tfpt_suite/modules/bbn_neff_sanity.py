from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule, mk_check_fail, mk_check_pass


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


@dataclass(frozen=True)
class GStarPoint:
    T_GeV: float
    g_star_energy: float
    note: str


def _g_star_energy_bbn(*, T_GeV: float, m_e_GeV: float) -> float:
    """
    Minimal g_* (energy density) for the MeV-era plasma relevant to BBN:
    photons + e± + 3 neutrinos.

    Convention:
    - T is the photon temperature.
    - neutrinos are assumed to decouple before e± annihilation; for T < m_e we apply T_ν/T = (4/11)^(1/3).
    """
    T = float(T_GeV)
    if T <= 0:
        raise ValueError("T_GeV must be positive")

    g_gamma = 2.0
    g_e = 4.0  # e- and e+
    g_nu = 6.0  # 3 species × (ν + anti-ν), treated as effectively massless
    factor_f = 7.0 / 8.0

    include_e = bool(T >= float(m_e_GeV))
    if include_e:
        Tnu_over_T = 1.0
    else:
        Tnu_over_T = float((4.0 / 11.0) ** (1.0 / 3.0))

    g = g_gamma
    if include_e:
        g += factor_f * g_e
    g += factor_f * g_nu * (Tnu_over_T**4)
    return float(g)


class BbnNeffSanityModule(TfptModule):
    module_id = "bbn_neff_sanity"
    title = "BBN / N_eff sanity (MeV-era g_* ledger; no full BBN simulation)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "Lepton masses: tfpt_suite/data/lepton_masses_pdg.json (electron mass threshold)",
            ],
            outputs=[
                "g_*(T) energy-density ledger around the e± annihilation threshold (MeV era)",
                "N_eff baseline marker (3.046; placeholder for full BBN/CMB likelihood)",
            ],
            formulas=[
                r"g_* = 2 + (7/8) g_e + (7/8) g_\nu (T_\nu/T)^4, with T_\nu/T = 1 above m_e and (4/11)^{1/3} below",
            ],
            validation=[
                "g_*(T≈1 MeV) ≈ 10.75 (photons + e± + 3 neutrinos)",
                "g_*(T≪m_e) ≈ 3.36 (photons + reheated neutrinos)",
            ],
            determinism="Deterministic given inputs.",
            question="Is the MeV-era thermodynamic bookkeeping consistent (g_* before/after e± annihilation; N_eff marker)?",
            objective=[
                "Provide a minimal BBN/N_eff sanity check without pretending to run a full BBN code.",
                "Make the temperature conventions explicit (T is photon temperature; neutrino reheating applied below m_e).",
            ],
            gaps=[
                "A publication-grade BBN module would integrate weak rates, nuclear network, and use a likelihood for (D/H, Y_p, N_eff).",
            ],
        )

    def run(self, config) -> ModuleResult:
        lep_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        lep_raw = json.loads(lep_path.read_text(encoding="utf-8"))
        m_e = float(lep_raw["masses"]["electron"]["mean"])

        # MeV-era temperature grid (photon temperature).
        Ts = np.array([1.0e-4, 3.0e-4, 5.0e-4, 1.0e-3, 3.0e-3, 1.0e-2], dtype=float)  # GeV
        pts: list[GStarPoint] = []
        for T in Ts:
            g = _g_star_energy_bbn(T_GeV=float(T), m_e_GeV=m_e)
            note = "above m_e (e± present; Tν=T)" if float(T) >= m_e else "below m_e (e± annihilated; Tν/T=(4/11)^(1/3))"
            pts.append(GStarPoint(T_GeV=float(T), g_star_energy=float(g), note=note))

        # Reference sanity points
        g_1MeV = float(_g_star_energy_bbn(T_GeV=1.0e-3, m_e_GeV=m_e))
        g_0p1MeV = float(_g_star_energy_bbn(T_GeV=1.0e-4, m_e_GeV=m_e))

        checks: list[Check] = []
        if abs(g_1MeV - 10.75) < 0.05:
            checks.append(mk_check_pass("g_star_at_1MeV", f"g*(1 MeV)={g_1MeV:.6g} (expected ≈10.75)"))
        else:
            checks.append(mk_check_fail("g_star_at_1MeV", f"g*(1 MeV)={g_1MeV:.6g} (expected ≈10.75)"))

        if abs(g_0p1MeV - 3.36) < 0.05:
            checks.append(mk_check_pass("g_star_below_me", f"g*(0.1 MeV)={g_0p1MeV:.6g} (expected ≈3.36)"))
        else:
            checks.append(mk_check_fail("g_star_below_me", f"g*(0.1 MeV)={g_0p1MeV:.6g} (expected ≈3.36)"))

        # N_eff marker (not derived here)
        N_eff = 3.046
        checks.append(mk_check_pass("neff_marker_present", f"N_eff marker={N_eff} (standard; placeholder for full likelihood)"))

        lines: list[str] = []
        lines += [
            "BBN / N_eff sanity (MeV-era g_* ledger; no full BBN simulation)",
            "",
            f"lepton masses: {_relpath(lep_path)} (sha256={_sha256_file(lep_path)})",
            f"electron mass threshold: m_e={m_e:.6g} GeV",
            "",
            "g_*(T) (energy density; T is photon temperature):",
        ]
        for p in pts:
            lines.append(f"- T={p.T_GeV:.3e} GeV: g*={p.g_star_energy:.6g}  ({p.note})")
        lines += [
            "",
            f"N_eff marker (standard): {N_eff}",
            "",
            "Checks:",
            *[f"- {c.check_id}: {'PASS' if c.passed else 'FAIL'} ({c.detail})" for c in checks],
        ]

        return ModuleResult(
            results={
                "inputs": {"lepton_masses_file": _relpath(lep_path), "lepton_masses_sha256": _sha256_file(lep_path)},
                "electron_mass_GeV": m_e,
                "g_star_energy_points": [p.__dict__ for p in pts],
                "neff_marker": N_eff,
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

