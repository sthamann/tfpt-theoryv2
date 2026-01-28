from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.conventions import alpha_from_g, gY_from_g1_gut
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule
from tfpt_suite.rge_sm import msbar_mass_run_qcd_1loop, run_sm_gauge_only_2loop_thresholds
from tfpt_suite.sm_inputs import SmMzInputs, gauge_couplings_from_mz_inputs


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
class AlphaSPoint:
    mu_GeV: float
    alpha_s_nf5_above: float | None
    alpha_s_matched: float
    nf_effective: int


class BelowMtEftCascadeModule(TfptModule):
    module_id = "below_mt_eft_cascade"
    title = "Below-mt EFT cascade (QCD thresholds + running masses; explicit bookkeeping)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "SM inputs at MZ: tfpt_suite/data/sm_inputs_mz.json (MZ, α_em_inv, sin²θW, αs, mc, mb, mt)",
                "Lepton pole masses: tfpt_suite/data/lepton_masses_pdg.json (for QED running thresholds e,μ,τ)",
                "Below-MZ QED policy: tfpt_suite/data/below_mz_policy.json (explicit charged-fermion thresholds)",
            ],
            outputs=[
                "αs(μ) across mc/mb thresholds (above+below; 2-loop matching for αs)",
                "mb(μ) at selected IR scales (LO QCD running, piecewise nf)",
                "α_em(μ) below MZ via a declared QED EFT policy (EW integrated out at MZ; 1-loop QED running with charged-fermion thresholds)",
                "diagnostic table suitable as an EFT-cascade audit trail",
            ],
            formulas=[
                "αs = g3^2/(4π)",
                "αs matching at quark thresholds: 2-loop MSbar decoupling (implemented in rge_sm.run_sm_gauge_only_2loop_thresholds)",
                "mb running (LO): m(μ1)=m(μ0)[αs(μ1)/αs(μ0)]^(12/(33-2nf)) (piecewise nf)",
                "QED (1-loop, MSbar): dα/d ln μ = (2/(3π)) α^2 Σ_f N_c Q_f^2  ⇒  1/α(μ1)=1/α(μ0) - (2/(3π))Σ ln(μ1/μ0)",
            ],
            validation=[
                "αs values are finite/positive and show the expected monotonic increase toward the IR",
                "mb(μ) is finite/positive for all reported scales",
                "α_em decreases toward the IR under the QED EFT policy (since QED β>0 in the UV)",
            ],
            determinism="Deterministic given the SM input table (no fitter).",
        )

    def run(self, config) -> ModuleResult:
        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"
        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        lep_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        lep_raw = json.loads(lep_path.read_text(encoding="utf-8"))
        policy_path = Path(__file__).resolve().parent.parent / "data" / "below_mz_policy.json"
        policy_raw = json.loads(policy_path.read_text(encoding="utf-8"))

        MZ = float(sm_raw.get("mu_GeV", 91.1876))
        alpha_em_inv = float(sm_raw["alpha_em_inv"])
        sin2 = float(sm_raw["sin2_thetaW"])
        alpha_s_mz = float(sm_raw["alpha_s"])

        mc = float(sm_raw.get("mc_GeV", 1.27))
        mb = float(sm_raw.get("mb_GeV", 4.18))
        mt = float(sm_raw.get("mt_GeV", 172.76))

        m_e = float(lep_raw["masses"]["electron"]["mean"])
        m_mu = float(lep_raw["masses"]["muon"]["mean"])
        m_tau = float(lep_raw["masses"]["tau"]["mean"])

        inp = SmMzInputs(mu_GeV=MZ, alpha_em_inv=alpha_em_inv, sin2_thetaW=sin2, alpha_s=alpha_s_mz)
        g1_gut_mz, g2_mz, g3_mz = gauge_couplings_from_mz_inputs(inp)
        gY_mz = float(gY_from_g1_gut(float(g1_gut_mz)))
        alpha_em_mz = 1.0 / float(alpha_em_inv)

        # Helper: g3(μ) from MZ→μ with optional αs matching at heavy-quark thresholds.
        def g3_at(mu_end: float, *, apply_matching: bool) -> float:
            g = run_sm_gauge_only_2loop_thresholds(
                mu_start_GeV=MZ,
                mu_end_GeV=float(mu_end),
                g_start=(float(g1_gut_mz), float(g2_mz), float(g3_mz)),
                mc_GeV=mc,
                mb_GeV=mb,
                mt_GeV=mt,
                apply_alpha3_matching=bool(apply_matching),
            )
            return float(g[2])

        # αs values:
        # - "above" at mb/mc: run without matching (keeps nf=5 all the way down; gives αs just above threshold)
        # - "matched" at mb/mc: run with matching (gives αs just below threshold in the lower-nf EFT)
        alpha_s_mt_nf5 = float(alpha_from_g(g3_at(mt, apply_matching=False)))
        alpha_s_mb_nf5_above = float(alpha_from_g(g3_at(mb, apply_matching=False)))
        alpha_s_mb_nf4_below = float(alpha_from_g(g3_at(mb, apply_matching=True)))
        alpha_s_mc_nf4_above = float(alpha_from_g(g3_at(mc, apply_matching=False)))  # nf≈4 above mc when mb matching applied? (policy: no matching => nf=5 drift)
        alpha_s_mc_nf3_below = float(alpha_from_g(g3_at(mc, apply_matching=True)))

        alpha_s_2gev = float(alpha_from_g(g3_at(2.0, apply_matching=True)))

        # mb running (LO, piecewise):
        # Treat mb(mb) as the MSbar running mass at μ=mb in the nf=5 theory (PDG-ish input).
        mb_mb = float(sm_raw.get("mb_GeV", 4.18))

        # Upward (nf=5): mb(mb) -> mb(MZ), mb(mt)
        mb_MZ = msbar_mass_run_qcd_1loop(m_mu0_GeV=mb_mb, alpha_s_mu0=alpha_s_mb_nf5_above, alpha_s_mu1=alpha_s_mz, nf=5)
        mb_mt = msbar_mass_run_qcd_1loop(m_mu0_GeV=mb_mb, alpha_s_mu0=alpha_s_mb_nf5_above, alpha_s_mu1=alpha_s_mt_nf5, nf=5)

        # Downward below mb (nf=4): approximate continuity of mb at μ=mb (LO) and run to 2 GeV.
        mb_2gev = msbar_mass_run_qcd_1loop(m_mu0_GeV=mb_mb, alpha_s_mu0=alpha_s_mb_nf4_below, alpha_s_mu1=alpha_s_2gev, nf=4)

        # --- QED policy below MZ (EW integrated out) ---
        # We run α_em in a QED EFT with charged fermions as active degrees of freedom.
        # One-loop MSbar decoupling is log-only, so matching at μ=m_f is continuous at this order.
        qed_policy = policy_raw.get("qed_running", {}) if isinstance(policy_raw, dict) else {}
        fields = qed_policy.get("active_fields", []) if isinstance(qed_policy, dict) else []
        if not isinstance(fields, list) or not fields:
            raise ValueError("below_mz_policy.json missing qed_running.active_fields")
        fermions = []
        for entry in fields:
            if not isinstance(entry, dict):
                continue
            fid = str(entry.get("id", "")).strip()
            if not fid:
                raise ValueError("below_mz_policy.json contains a field with empty id")
            mass = float(entry.get("mass_GeV", 0.0))
            nc = int(entry.get("Nc", 1))
            charge = float(entry.get("charge", 0.0))
            fermions.append((fid, mass, nc, charge))

        def sum_q2(mu_GeV: float) -> float:
            mu = float(mu_GeV)
            s = 0.0
            for _, m, nc, q in fermions:
                if float(m) <= 0.0 or mu >= float(m):
                    s += float(nc) * float(q) ** 2
            return float(s)

        def alpha_em_run_1loop(*, mu_start_GeV: float, alpha_start: float, mu_end_GeV: float) -> float:
            mu0 = float(mu_start_GeV)
            mu1 = float(mu_end_GeV)
            if mu0 <= 0 or mu1 <= 0:
                raise ValueError("mu must be positive")
            if alpha_start <= 0:
                raise ValueError("alpha_start must be positive")
            if mu0 == mu1:
                return float(alpha_start)
            direction = "up" if mu1 > mu0 else "down"
            lo, hi = (mu0, mu1) if mu0 < mu1 else (mu1, mu0)
            # breakpoints = endpoints + all fermion masses in between (excluding 0-mass EFT fields)
            cuts = [mu0, mu1] + [float(m) for _, m, _, _ in fermions if float(m) > 0.0 and lo < float(m) < hi]
            cuts = sorted(set(cuts))
            if direction == "down":
                cuts = list(reversed(cuts))
            a = float(alpha_start)
            for a_mu, b_mu in zip(cuts[:-1], cuts[1:]):
                mu_mid = float(math.sqrt(float(a_mu) * float(b_mu)))
                s_q2 = float(sum_q2(mu_mid))
                b = (2.0 / (3.0 * math.pi)) * s_q2
                inv_a_mu = 1.0 / a
                inv_a_next = inv_a_mu - b * math.log(float(b_mu) / float(a_mu))
                if inv_a_next <= 0:
                    raise RuntimeError("QED 1-loop running hit a non-positive 1/alpha (Landau pole / invalid range)")
                a = 1.0 / inv_a_next
            return float(a)

        alpha_em_mb = alpha_em_run_1loop(mu_start_GeV=MZ, alpha_start=alpha_em_mz, mu_end_GeV=mb)
        alpha_em_2gev = alpha_em_run_1loop(mu_start_GeV=MZ, alpha_start=alpha_em_mz, mu_end_GeV=2.0)
        alpha_em_mc = alpha_em_run_1loop(mu_start_GeV=MZ, alpha_start=alpha_em_mz, mu_end_GeV=mc)

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="below_MZ_policy_explicit_and_applied",
                passed=True,
                detail=f"policy={_relpath(policy_path)}; fields={[f[0] for f in fermions]}",
            )
        )
        checks.append(
            Check(
                check_id="alpha_s_finite_positive",
                passed=bool(all(np.isfinite(x) and x > 0 for x in [alpha_s_mz, alpha_s_mt_nf5, alpha_s_mb_nf5_above, alpha_s_mb_nf4_below, alpha_s_2gev])),
                detail=f"αs(MZ)={alpha_s_mz:.6g}, αs(mt,nf5)={alpha_s_mt_nf5:.6g}, αs(mb^+)={alpha_s_mb_nf5_above:.6g}, αs(mb^-)={alpha_s_mb_nf4_below:.6g}, αs(2GeV)={alpha_s_2gev:.6g}",
            )
        )
        checks.append(
            Check(
                check_id="alpha_s_ir_increases",
                passed=bool(alpha_s_2gev > alpha_s_mb_nf4_below > alpha_s_mz),
                detail=f"αs(2GeV)={alpha_s_2gev:.6g} > αs(mb^- )={alpha_s_mb_nf4_below:.6g} > αs(MZ)={alpha_s_mz:.6g}",
            )
        )
        checks.append(
            Check(
                check_id="mb_running_finite_positive",
                passed=bool(all(np.isfinite(x) and x > 0 for x in [mb_mb, mb_MZ, mb_mt, mb_2gev])),
                detail=f"mb(mb)={mb_mb:.6g} → mb(MZ)={mb_MZ:.6g}, mb(mt)={mb_mt:.6g}, mb(2GeV)={mb_2gev:.6g}",
            )
        )
        checks.append(
            Check(
                check_id="alpha_em_finite_positive",
                passed=bool(all(np.isfinite(x) and x > 0 for x in [alpha_em_mz, alpha_em_mb, alpha_em_2gev, alpha_em_mc])),
                detail=f"αem(MZ)={alpha_em_mz:.9g}, αem(mb)={alpha_em_mb:.9g}, αem(2GeV)={alpha_em_2gev:.9g}, αem(mc)={alpha_em_mc:.9g}",
            )
        )
        checks.append(
            Check(
                check_id="alpha_em_ir_decreases",
                passed=bool(alpha_em_mz > alpha_em_mb > alpha_em_2gev > 0 and alpha_em_2gev > alpha_em_mc),
                detail=f"αem(MZ)={alpha_em_mz:.9g} > αem(mb)={alpha_em_mb:.9g} > αem(2GeV)={alpha_em_2gev:.9g} > αem(mc)={alpha_em_mc:.9g}",
            )
        )

        lines: list[str] = []
        lines += [
            "Below-mt EFT cascade (QCD thresholds + running masses; explicit bookkeeping)",
            "",
            f"SM inputs: {_relpath(sm_path)} (sha256={_sha256_file(sm_path)})",
            f"lepton masses: {_relpath(lep_path)} (sha256={_sha256_file(lep_path)})",
            f"below-MZ policy: {_relpath(policy_path)} (sha256={_sha256_file(policy_path)})",
            "",
            "Gauge initialization at MZ:",
            f"- MZ={MZ:.6g} GeV, αs(MZ)={alpha_s_mz:.6g} (input), αem(MZ)={alpha_em_mz:.9g} (input), gY(MZ)={gY_mz:.6g}",
            "",
            "αs across thresholds (policy: 2-loop gauge-only running with optional 2-loop αs matching at thresholds):",
            f"- αs(mt={mt:.6g} GeV) [nf=5, no top matching] ≈ {alpha_s_mt_nf5:.6g}",
            f"- αs(mb={mb:.6g} GeV)^+ [nf=5, above threshold] ≈ {alpha_s_mb_nf5_above:.6g}",
            f"- αs(mb={mb:.6g} GeV)^- [matched down, below threshold] ≈ {alpha_s_mb_nf4_below:.6g}",
            f"- αs(mc={mc:.6g} GeV)^- [matched down] ≈ {alpha_s_mc_nf3_below:.6g}",
            f"- αs(2 GeV) [matched down] ≈ {alpha_s_2gev:.6g}",
            "",
            "αem below MZ (policy: EW integrated out at MZ; QED 1-loop running with charged-fermion thresholds):",
            f"- Σ_f N_c Q_f^2 at MZ: {sum_q2(MZ):.6g} (active: e,μ,τ,u,d,s,c,b)",
            f"- αem(mb={mb:.6g} GeV) ≈ {alpha_em_mb:.9g}",
            f"- αem(2 GeV) ≈ {alpha_em_2gev:.9g}",
            f"- αem(mc={mc:.6g} GeV) ≈ {alpha_em_mc:.9g} (includes τ threshold at mτ={m_tau:.6g} GeV)",
            "",
            "mb running (LO QCD; explicit approximation):",
            f"- mb(mb) input = {mb_mb:.6g} GeV",
            f"- mb(MZ) (nf=5 LO) ≈ {mb_MZ:.6g} GeV",
            f"- mb(mt) (nf=5 LO) ≈ {mb_mt:.6g} GeV",
            f"- mb(2 GeV) (nf=4 LO; continuity at mb assumed) ≈ {mb_2gev:.6g} GeV",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- This module is an audit trail for the below-mt EFT cascade. It is intentionally explicit about what is matched and what is still approximated.",
            "- Publication-grade upgrades: higher-loop QCD running for masses, explicit mass matching at thresholds, and (if needed) higher-loop QED matching/decoupling beyond the 1-loop identity-at-μ=m_f policy used here.",
        ]

        return ModuleResult(
            results={
                "inputs": {
                    "sm_inputs_file": _relpath(sm_path),
                    "sm_inputs_sha256": _sha256_file(sm_path),
                    "lepton_masses_file": _relpath(lep_path),
                    "lepton_masses_sha256": _sha256_file(lep_path),
                    "below_mz_policy_file": _relpath(policy_path),
                    "below_mz_policy_sha256": _sha256_file(policy_path),
                },
                "alpha_s": {
                    "MZ": {"mu_GeV": MZ, "alpha_s": alpha_s_mz},
                    "mt_nf5": {"mu_GeV": mt, "alpha_s": alpha_s_mt_nf5},
                    "mb_above_nf5": {"mu_GeV": mb, "alpha_s": alpha_s_mb_nf5_above},
                    "mb_below_nf4": {"mu_GeV": mb, "alpha_s": alpha_s_mb_nf4_below},
                    "mc_below_nf3": {"mu_GeV": mc, "alpha_s": alpha_s_mc_nf3_below},
                    "mu_2_GeV": {"mu_GeV": 2.0, "alpha_s": alpha_s_2gev},
                },
                "alpha_em": {
                    "MZ": {"mu_GeV": MZ, "alpha_em": alpha_em_mz},
                    "mb": {"mu_GeV": mb, "alpha_em": alpha_em_mb},
                    "mu_2_GeV": {"mu_GeV": 2.0, "alpha_em": alpha_em_2gev},
                    "mc": {"mu_GeV": mc, "alpha_em": alpha_em_mc},
                },
                "qed_policy": {
                    "ew_integrated_out_below_MZ": True,
                    "beta_1loop": "dα/d ln μ = (2/(3π)) α^2 Σ_f N_c Q_f^2",
                    "charged_fermions": [
                        {"name": n, "mass_GeV": float(m), "Nc": int(nc), "Q": float(q)} for (n, m, nc, q) in fermions
                    ],
                },
                "mb_running_lo": {
                    "mb_mb_GeV": mb_mb,
                    "mb_MZ_GeV": mb_MZ,
                    "mb_mt_GeV": mb_mt,
                    "mb_2GeV_GeV": mb_2gev,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

