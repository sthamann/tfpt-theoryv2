from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from tfpt_suite.constants import TfptConstants
from tfpt_suite.mobius_z3_yukawa_generator import _ckm_from_yukawas, generate_quark_yukawas_mt
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


class MobiusZ3YukawaGeneratorModule(TfptModule):
    module_id = "mobius_z3_yukawa_generator"
    title = "Möbius/Z3 Yukawa generator (v1.07SM-style CKM + Möbius mass hierarchies)"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "TFPT invariants (c3,varphi0,delta_star) via constants.py",
                "flavor config: tfpt_suite/data/flavor_texture_v24.json (delta_source policy)",
                "lepton masses: tfpt_suite/data/lepton_masses_pdg.json (δ_M from τ/μ; diagnostic)",
                "CKM reference: tfpt_suite/data/ckm_reference.json (for χ² diagnostic)",
                "SM inputs: tfpt_suite/data/sm_inputs_mz.json (yt(mt), yb(mt) targets)",
            ],
            outputs=[
                "Yu(mt), Yd(mt) matrices (complex)",
                "CKM |V_ij|(mt) from diagonalization and χ² vs reference (diagnostic)",
                "explicit generator metadata (δ used, CKM construction modes, Möbius ratio predictions)",
            ],
            formulas=[
                "δ_M = (sqrt(mτ/mμ)-1)/(sqrt(mτ/mμ)+1) (diagnostic, if delta_source=delta_star then δ_* is used instead)",
                "δ_* = 3/5 + varphi0/6 (TFPT geometric closure)",
                "λ = sqrt(varphi0) * (1 - varphi0/2)",
                "A = 5/6 (from Z3 cusp slopes 2 and 1)",
                "s23 = A λ^2; s13 = A λ^3 (1-δ); δ_CP = π(1-δ) (explicit convention)",
                "hierarchies from Möbius relations: ms/md=(M1(δ))^2, mb/ms=(M1(δ)(1+δ))^2, mc/mu=(M2/3(δ))^2, mt/mc=(2/3/(2/3-δ))^2",
            ],
            validation=[
                "Generated CKM is unitary (numerical) and χ² is finite",
            ],
            determinism="Deterministic given input tables; no optimizer/scan.",
        )

    def run(self, config) -> ModuleResult:
        c = TfptConstants.compute()

        flavor_cfg_path = Path(__file__).resolve().parent.parent / "data" / "flavor_texture_v24.json"
        lepton_masses_path = Path(__file__).resolve().parent.parent / "data" / "lepton_masses_pdg.json"
        ckm_ref_path = Path(__file__).resolve().parent.parent / "data" / "ckm_reference.json"
        sm_path = Path(__file__).resolve().parent.parent / "data" / "sm_inputs_mz.json"

        flavor_cfg = json.loads(flavor_cfg_path.read_text(encoding="utf-8"))
        tex = flavor_cfg.get("yukawa_texture", {}) if isinstance(flavor_cfg.get("yukawa_texture", {}), dict) else {}
        delta_source = str(tex.get("delta_source", "tau_mu")).strip()

        lep = json.loads(lepton_masses_path.read_text(encoding="utf-8"))
        mtau = float(lep["masses"]["tau"]["mean"])
        mmu = float(lep["masses"]["muon"]["mean"])
        R = float(np.sqrt(mtau / mmu))
        delta_M = float((R - 1.0) / (R + 1.0))
        delta_star = float(c.delta_star)
        if delta_source == "tau_mu":
            delta_used = delta_M
        elif delta_source in ("delta_star", "delta_*", "delta_star_varphi0"):
            delta_used = delta_star
        else:
            raise ValueError(f"Unsupported delta_source in {flavor_cfg_path}: {delta_source}")

        sm_raw = json.loads(sm_path.read_text(encoding="utf-8"))
        yt_mt = float(sm_raw.get("yt_mt_target", 0.949758))  # fall back to typical mt estimate
        yb_mt = float(sm_raw.get("yb_mt", 0.017))

        Yu, Yd, meta = generate_quark_yukawas_mt(varphi0=float(c.varphi0), delta_used=delta_used, yt_mt=yt_mt, yb_mt=yb_mt)
        V = _ckm_from_yukawas(Yu, Yd)
        Vabs = {k: float(np.abs(V[i, j])) for k, (i, j) in {
            "Vud": (0, 0), "Vus": (0, 1), "Vub": (0, 2),
            "Vcd": (1, 0), "Vcs": (1, 1), "Vcb": (1, 2),
            "Vtd": (2, 0), "Vts": (2, 1), "Vtb": (2, 2),
        }.items()}

        ref = json.loads(ckm_ref_path.read_text(encoding="utf-8"))["matrix_abs"]
        chi2 = 0.0
        contrib = []
        for k, cfg in ref.items():
            mean = float(cfg["mean"])
            sigma = float(cfg["sigma"])
            pred = float(Vabs[k])
            c2 = ((pred - mean) / sigma) ** 2
            chi2 += c2
            contrib.append({"key": k, "pred": pred, "mean": mean, "sigma": sigma, "chi2": c2})
        contrib.sort(key=lambda t: t["chi2"], reverse=True)

        unitarity_dev = float(np.max(np.abs(V.conj().T @ V - np.eye(3))))

        checks = [
            Check(check_id="ckm_unitarity", passed=bool(unitarity_dev < 1e-10), detail=f"max|V†V-I|={unitarity_dev:.3e}"),
            Check(check_id="ckm_chi2_finite", passed=bool(np.isfinite(chi2) and chi2 >= 0), detail=f"chi2={chi2:.6g}"),
        ]

        lines: list[str] = []
        lines += [
            "Möbius/Z3 Yukawa generator (v1.07SM-style CKM + Möbius mass hierarchies)",
            "",
            f"flavor config file: {_relpath(flavor_cfg_path)} (sha256={_sha256_file(flavor_cfg_path)})",
            f"lepton masses file: {_relpath(lepton_masses_path)} (sha256={_sha256_file(lepton_masses_path)})",
            f"CKM reference file: {_relpath(ckm_ref_path)} (sha256={_sha256_file(ckm_ref_path)})",
            f"SM inputs file: {_relpath(sm_path)} (sha256={_sha256_file(sm_path)})",
            "",
            f"delta_source = {delta_source}",
            f"delta_M (from tau/mu) = {delta_M}",
            f"delta_star (from geometry) = {delta_star}",
            f"delta_used (generator input) = {delta_used}",
            f"generator meta: {meta}",
            "",
            "CKM |V_ij| (mt):",
        ]
        for k in ["Vud","Vus","Vub","Vcd","Vcs","Vcb","Vtd","Vts","Vtb"]:
            lines.append(f"- {k} = {Vabs[k]}")
        lines += [
            "",
            f"chi2 = {chi2}",
            "largest chi2 contributions (top 5):",
        ]
        for t in contrib[:5]:
            lines.append(f"- {t}")
        lines += [
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
        ]

        return ModuleResult(
            results={
                "delta_source": delta_source,
                "delta_M": delta_M,
                "delta_star": delta_star,
                "delta_used": delta_used,
                "yukawas_mt": {"Yu": Yu.tolist(), "Yd": Yd.tolist()},
                "ckm_mt_abs": Vabs,
                "ckm_mt_unitarity_dev": unitarity_dev,
                "chi2": chi2,
                "chi2_contributions": contrib,
                "meta": {
                    "scheme": meta.scheme,
                    "reference_scale_GeV": meta.reference_scale_GeV,
                    "phase_mode": meta.phase_mode,
                    "ckm": meta.ckm.__dict__,
                    "ratios": meta.ratios,
                },
            },
            checks=checks,
            report="\n".join(lines),
            warnings=[],
        )

