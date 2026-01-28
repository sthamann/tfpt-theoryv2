from __future__ import annotations

import argparse
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _suite_root() -> Path:
    # tfpt-suite/tfpt_suite/export_prediction_ledger.py -> tfpt-suite/
    return Path(__file__).resolve().parents[1]


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _relpath(path: Path, *, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root))
    except Exception:
        return str(path)


def _tex_escape(s: str) -> str:
    return (
        s.replace("\\", "\\textbackslash{}")
        .replace("_", "\\_")
        .replace("%", "\\%")
        .replace("&", "\\&")
        .replace("#", "\\#")
        .replace("{", "\\{")
        .replace("}", "\\}")
        .replace("^", "\\textasciicircum{}")
        .replace("~", "\\textasciitilde{}")
    )


def _num(x: Any) -> str:
    if x is None:
        return "---"
    try:
        return f"\\num{{{x}}}"
    except Exception:
        return "---"


_NUM_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def _tex_cell(x: Any) -> str:
    """
    Best-effort TeX cell renderer:
    - numeric -> \\num{...}
    - non-numeric -> TeX-escaped text
    """
    if x is None:
        return "---"
    if isinstance(x, (int, float)):
        return _num(x)
    if isinstance(x, str):
        s = x.strip()
        if _NUM_RE.match(s):
            return _num(s)
        return _tex_escape(s)
    return _tex_escape(str(x))


def _ref_with_sigma(*, mean: Any, sigma: Any) -> str:
    if mean is None:
        return "---"
    if sigma is None:
        return _tex_cell(mean)
    return rf"{_tex_cell(mean)} $\pm$ {_tex_cell(sigma)}"


def export_prediction_ledger_outputs(*, write_json: Path, write_tex: Path) -> tuple[Path, Path]:
    """
    Export a reviewer-friendly prediction ledger from suite outputs.

    This is intentionally conservative: it only includes a curated set of headline
    observables whose sources are already tracked by the suite.
    """
    suite_root = _suite_root()

    # --- Load sources (skip gracefully if missing) ---
    sources: dict[str, Any] = {}
    def _maybe_load(key: str, rel: str) -> None:
        p = (suite_root / Path(rel)).resolve()
        if p.is_file():
            sources[key] = {"path": p, "payload": _load_json(p)}

    _maybe_load("defect_g5", "theoryv3/out/defect_partition_g5_audit/results.json")
    _maybe_load("global", "out/global_consistency_test/results.json")
    _maybe_load("omega_b", "out/omega_b_conjecture_scan/results.json")
    _maybe_load("dark_energy_exp", "theoryv3/out/dark_energy_exponential_audit/results.json")
    _maybe_load("torsion_cond", "out_physics/torsion_condensate/results.json")
    _maybe_load("axion_dm", "theoryv3/out/axion_dm_audit/results.json")
    _maybe_load("pmns_z3", "out/pmns_z3_breaking/results.json")
    _maybe_load("pmns_ref", "tfpt_suite/data/pmns_reference.json")

    rows: list[dict[str, Any]] = []

    def _add_row(
        *,
        group: str,
        observable_tex: str,
        formula_tex: str,
        pred: Any,
        ref_mean: Any | None = None,
        ref_sigma: Any | None = None,
        z: Any | None = None,
        status: str,
        module_id: str,
        output_rel: str,
    ) -> None:
        rows.append(
            {
                "group": group,
                "observable_tex": observable_tex,
                "formula_tex": formula_tex,
                "pred": pred,
                "ref_mean": ref_mean,
                "ref_sigma": ref_sigma,
                "z": z,
                "status": status,
                "module_id": module_id,
                "output_rel": output_rel,
            }
        )

    # --- Core, parameter-free / kernel-level ---
    if "defect_g5" in sources:
        payload = sources["defect_g5"]["payload"]
        out_rel = _relpath(sources["defect_g5"]["path"], root=suite_root)
        alpha = (payload.get("results", {}) or {}).get("alpha_inv_0", {})
        if isinstance(alpha, dict):
            _add_row(
                group="Core (parameter-free)",
                observable_tex=r"\(\alpha^{-1}(0)\)",
                formula_tex=r"\(\mathrm{CFE} + \varphi_0(\alpha)\) with \(\delta_2=(g/4)\delta_{\mathrm{top}}^2\) and \(g=5\) [Kernel+Grammar]",
                pred=alpha.get("pred", None),
                ref_mean=alpha.get("ref", None),
                ref_sigma=alpha.get("sigma", None),
                z=alpha.get("z", None),
                status="P",
                module_id=str(payload.get("module_id", "defect_partition_g5_audit")),
                output_rel=out_rel,
            )

    # g_{aγγ} is purely algebraic from kernel; include as a ledger entry.
    gagg = -1.0 / (2.0 * math.pi)
    _add_row(
        group="Core (parameter-free)",
        observable_tex=r"\(g_{a\gamma\gamma}\)",
        formula_tex=r"\(g_{a\gamma\gamma}=-4c_3=-\frac{1}{2\pi}\) [Kernel]",
        pred=str(gagg),
        ref_mean=None,
        ref_sigma=None,
        z=None,
        status="P",
        module_id="(algebraic)",
        output_rel="paper kernel (no suite output file)",
    )

    # --- Global consistency terms (comparison layer + core observables) ---
    if "global" in sources:
        payload = sources["global"]["payload"]
        out_rel = _relpath(sources["global"]["path"], root=suite_root)
        terms = (payload.get("results", {}) or {}).get("terms", [])
        if isinstance(terms, list):
            term_by_name = {t.get("name"): t for t in terms if isinstance(t, dict)}

            def _term_row(name: str, *, group: str, observable_tex: str, formula_tex: str, status: str) -> None:
                t = term_by_name.get(name, None)
                if not isinstance(t, dict):
                    return
                _add_row(
                    group=group,
                    observable_tex=observable_tex,
                    formula_tex=formula_tex,
                    pred=t.get("pred", None),
                    ref_mean=t.get("mean", None),
                    ref_sigma=t.get("sigma_exp", None),
                    z=t.get("z_total", t.get("z_exp", None)),
                    status=status,
                    module_id=str(payload.get("module_id", "global_consistency_test")),
                    output_rel=out_rel,
                )

            _term_row(
                "beta_deg",
                group="Core (parameter-free)",
                observable_tex=r"\(\beta_{\mathrm{deg}}\)",
                formula_tex=r"\(\beta_{\mathrm{deg}}=\frac{180}{\pi}\frac{\varphi_0}{4\pi}\) [Kernel]",
                status="P",
            )
            _term_row(
                "cabibbo_lambda",
                group="Core (parameter-free)",
                observable_tex=r"\(\lambda\equiv |V_{us}|\)",
                formula_tex=r"\(\lambda=\sqrt{\varphi_0}\left(1-\frac{1}{2}\varphi_0\right)\) [Kernel]",
                status="P",
            )
            _term_row(
                "n_s",
                group="Conditional dynamics",
                observable_tex=r"\(n_s\)",
                formula_tex=r"Starobinsky \(R^2\) completion: \(n_s=1-\frac{2}{N}\) (benchmark \(N=56\)) [Conditional]",
                status="conditional",
            )
            _term_row(
                "A_s",
                group="Conditional dynamics",
                observable_tex=r"\(A_s\)",
                formula_tex=r"Starobinsky \(R^2\) completion: \(A_s\simeq \frac{N^2}{24\pi^2}(M/\bar M_{\rm Pl})^2\) [Conditional]",
                status="conditional",
            )
            _term_row(
                "alpha_bar5_inv_MZ",
                group="Comparison layer",
                observable_tex=r"\(\overline{\alpha}^{(5)}(M_Z)^{-1}\)",
                formula_tex=r"QED running + matching (declared policy) from \(\alpha(0)\) [Policy]",
                status="comparison",
            )

        # r is a predicted scalar here; suite treats only an upper bound proxy in this module.
        preds = (payload.get("results", {}) or {}).get("predictions", {})
        r_proxy = (payload.get("results", {}) or {}).get("r_bound_proxy", {})
        if isinstance(preds, dict) and isinstance(r_proxy, dict):
            r_pred = preds.get("r", None)
            r_upper = r_proxy.get("r_upper_95", None)
            if r_pred is not None:
                _add_row(
                    group="Conditional dynamics",
                    observable_tex=r"\(r\)",
                    formula_tex=r"Starobinsky \(R^2\) completion: \(r=\frac{12}{N^2}\) (benchmark \(N=56\)); compared to an upper bound [Conditional]",
                    pred=r_pred,
                    ref_mean=f"< {r_upper}" if r_upper is not None else None,
                    ref_sigma=None,
                    z=None,
                    status="conditional",
                    module_id=str(payload.get("module_id", "global_consistency_test")),
                    output_rel=out_rel,
                )

    # --- PMNS: closed sub-block from identity + TM1; references from NuFIT ---
    if "pmns_z3" in sources and "pmns_ref" in sources:
        pmns = sources["pmns_z3"]["payload"]
        pmns_ref = sources["pmns_ref"]["payload"]
        out_rel = _relpath(sources["pmns_z3"]["path"], root=suite_root)

        baseline = (pmns.get("results", {}) or {}).get("baseline", {})
        normal = (pmns_ref.get("normal_ordering", {}) or {})
        if isinstance(baseline, dict) and isinstance(normal, dict):
            theta13 = float(baseline.get("theta13_deg", 0.0))
            theta12 = float(baseline.get("theta12_deg", 0.0))
            sin2_13 = math.sin(math.radians(theta13)) ** 2
            sin2_12 = math.sin(math.radians(theta12)) ** 2

            ref_13 = normal.get("sin2_theta13", {})
            ref_12 = normal.get("sin2_theta12", {})
            if isinstance(ref_13, dict) and isinstance(ref_12, dict):
                mean13 = float(ref_13.get("mean", 0.0))
                sig13 = float(ref_13.get("sigma", 1.0))
                mean12 = float(ref_12.get("mean", 0.0))
                sig12 = float(ref_12.get("sigma", 1.0))
                z13 = (sin2_13 - mean13) / sig13 if sig13 else None
                z12 = (sin2_12 - mean12) / sig12 if sig12 else None

                _add_row(
                    group="Core (parameter-free)",
                    observable_tex=r"\(\sin^2\theta_{13}\)",
                    formula_tex=r"Identity + TM1 [Kernel+Grammar]",
                    pred=f"{sin2_13:.8f}",
                    ref_mean=f"{mean13}",
                    ref_sigma=f"{sig13}",
                    z=f"{z13:.6f}" if z13 is not None else None,
                    status="P",
                    module_id=str(pmns.get("module_id", "pmns_z3_breaking")),
                    output_rel=out_rel,
                )
                _add_row(
                    group="Core (parameter-free)",
                    observable_tex=r"\(\sin^2\theta_{12}\)",
                    formula_tex=r"TM1 sum rule from \(\sin^2\theta_{13}\) [Kernel+Grammar]",
                    pred=f"{sin2_12:.8f}",
                    ref_mean=f"{mean12}",
                    ref_sigma=f"{sig12}",
                    z=f"{z12:.6f}" if z12 is not None else None,
                    status="P",
                    module_id=str(pmns.get("module_id", "pmns_z3_breaking")),
                    output_rel=out_rel,
                )

    # --- Candidates / open interfaces ---
    if "omega_b" in sources:
        payload = sources["omega_b"]["payload"]
        out_rel = _relpath(sources["omega_b"]["path"], root=suite_root)
        ident = (payload.get("results", {}) or {}).get("identity", {})
        ref = (payload.get("results", {}) or {}).get("reference", {})
        if isinstance(ident, dict) and isinstance(ref, dict):
            pred = ident.get("omega_b_pred", None)
            mean = ref.get("omega_b_ref", None)
            sigma = ref.get("sigma_omega_b_ref", None)
            z = None
            try:
                if pred is not None and mean is not None and sigma is not None:
                    z = (float(pred) - float(mean)) / float(sigma)
            except Exception:
                z = None
            _add_row(
                group="Candidates / open interfaces",
                observable_tex=r"\(\Omega_b\)",
                formula_tex=r"\(\Omega_b=(4\pi-1)\beta_{\rm rad}\) under explicit sector-counting assumptions [Candidate]",
                pred=pred,
                ref_mean=mean,
                ref_sigma=sigma,
                z=f"{z:.6f}" if isinstance(z, float) else None,
                status="C (conditional)",
                module_id=str(payload.get("module_id", "omega_b_conjecture_scan")),
                output_rel=out_rel,
            )

    if "dark_energy_exp" in sources:
        payload = sources["dark_energy_exp"]["payload"]
        out_rel = _relpath(sources["dark_energy_exp"]["path"], root=suite_root)
        best = (payload.get("results", {}) or {}).get("best", {})
        targets = (payload.get("results", {}) or {}).get("targets", {})
        if isinstance(best, dict) and isinstance(targets, dict):
            _add_row(
                group="Candidates / open interfaces",
                observable_tex=r"\(\rho_\Lambda\) (normalization audit)",
                formula_tex=r"\(\phi_\*\!=n\,e^{-\alpha^{-1}(0)/2}\), \(\rho_\Lambda=(\bar M_{\rm Pl}\phi_\*)^4\); best \(n=1/2\) with log10 mismatch [Candidate]",
                pred=best.get("rho_L_GeV4", None),
                ref_mean=targets.get("rho_L_target", None),
                ref_sigma=None,
                z=None,
                status="C (PASS audit)",
                module_id=str(payload.get("module_id", "dark_energy_exponential_audit")),
                output_rel=out_rel,
            )

    if "torsion_cond" in sources:
        payload = sources["torsion_cond"]["payload"]
        out_rel = _relpath(sources["torsion_cond"]["path"], root=suite_root)
        best = ((payload.get("results", {}) or {}).get("gap_equation", {}) or {}).get("best", {})
        targets = (payload.get("results", {}) or {}).get("targets", {})
        if isinstance(best, dict) and isinstance(targets, dict):
            _add_row(
                group="Candidates / open interfaces",
                observable_tex=r"\(\rho_\Lambda\) (torsion condensate closure)",
                formula_tex=r"Gap-equation completion with discrete \(n\) (spectral flow); compare to cosmology target [Candidate]",
                pred=best.get("rho_L_GeV4", None),
                ref_mean=targets.get("rho_L_GeV4", None),
                ref_sigma=None,
                z=best.get("z_score_rho_L", None),
                status="C (FAIL gate)",
                module_id=str(payload.get("module_id", "torsion_condensate")),
                output_rel=out_rel,
            )

    if "axion_dm" in sources:
        payload = sources["axion_dm"]["payload"]
        out_rel = _relpath(sources["axion_dm"]["path"], root=suite_root)
        claim = (payload.get("results", {}) or {}).get("axion_claim", {})
        relic = (payload.get("results", {}) or {}).get("relic", {})
        if isinstance(claim, dict):
            _add_row(
                group="Conditional dynamics",
                observable_tex=r"\(f_a\) [GeV]",
                formula_tex=r"PQ cascade block (axion claim) [Conditional]",
                pred=claim.get("f_a_GeV", None),
                ref_mean=None,
                ref_sigma=None,
                z=None,
                status="conditional",
                module_id=str(payload.get("module_id", "axion_dm_audit")),
                output_rel=out_rel,
            )
            _add_row(
                group="Conditional dynamics",
                observable_tex=r"\(m_a\) [\(\mu\)eV]",
                formula_tex=r"PQ cascade block (axion claim) [Conditional]",
                pred=claim.get("m_a_micro_eV", None),
                ref_mean=None,
                ref_sigma=None,
                z=None,
                status="conditional",
                module_id=str(payload.get("module_id", "axion_dm_audit")),
                output_rel=out_rel,
            )
            _add_row(
                group="Conditional dynamics",
                observable_tex=r"\(\nu\) [GHz]",
                formula_tex=r"\(\nu_{\rm GHz}\approx 0.24179893\, m_a[\mu\mathrm{eV}]\) [Conditional]",
                pred=claim.get("nu_GHz", None),
                ref_mean=None,
                ref_sigma=None,
                z=None,
                status="conditional",
                module_id=str(payload.get("module_id", "axion_dm_audit")),
                output_rel=out_rel,
            )
        if isinstance(relic, dict):
            _add_row(
                group="Conditional dynamics",
                observable_tex=r"\(\Omega_a h^2\)",
                formula_tex=r"Axion relic accounting (scenario-dependent) [Conditional]",
                pred=relic.get("Omega_a_h2", None),
                ref_mean=None,
                ref_sigma=None,
                z=None,
                status="conditional",
                module_id=str(payload.get("module_id", "axion_dm_audit")),
                output_rel=out_rel,
            )

    # Sort rows by group then observable
    rows = sorted(rows, key=lambda r: (str(r.get("group", "")), str(r.get("observable_tex", ""))))

    now = datetime.now(tz=timezone.utc).isoformat()
    out_json = {"generated_at_utc": now, "rows": rows}
    write_json.parent.mkdir(parents=True, exist_ok=True)
    write_json.write_text(json.dumps(out_json, indent=2, sort_keys=False) + "\n", encoding="utf-8")

    # --- TeX output ---
    lines: list[str] = []
    lines.append("% AUTO-GENERATED by tfpt_suite/export_prediction_ledger.py — do not edit by hand.")
    lines.append(f"% generated_at_utc: {now}")
    lines.append("% sources: curated suite outputs (global_consistency_test, defect_partition_g5_audit, etc.)")
    lines.append(r"\providecommand{\codepath}[1]{\texttt{\detokenize{#1}}}")
    lines.append("")
    lines.append(r"\begin{longtable}{@{}p{0.15\textwidth}p{0.27\textwidth}rrrp{0.10\textwidth}p{0.18\textwidth}@{}}")
    lines.append(r"\caption{Prediction ledger (auto-generated from TFPT suite outputs).}\label{tab:prediction-ledger}\\")
    lines.append(r"\toprule")
    lines.append(r"Observable & Formula / dependencies & Pred. & Ref. & $z$ & Status & Suite evidence\\")
    lines.append(r"\midrule")
    lines.append(r"\endfirsthead")
    lines.append(r"\toprule")
    lines.append(r"Observable & Formula / dependencies & Pred. & Ref. & $z$ & Status & Suite evidence\\")
    lines.append(r"\midrule")
    lines.append(r"\endhead")
    lines.append(r"\midrule")
    lines.append(r"\multicolumn{7}{r}{\small Continued on next page.}\\")
    lines.append(r"\midrule")
    lines.append(r"\endfoot")
    lines.append(r"\bottomrule")
    lines.append(r"\endlastfoot")

    last_group: str | None = None
    for r in rows:
        group = str(r.get("group", ""))
        if group != last_group:
            gtex = _tex_escape(group)
            lines.append(rf"\multicolumn{{7}}{{@{{}}l}}{{\textbf{{{gtex}}}}}\\")
            last_group = group

        obs = str(r.get("observable_tex", ""))
        formula = str(r.get("formula_tex", ""))
        pred = r.get("pred", None)
        mean = r.get("ref_mean", None)
        sigma = r.get("ref_sigma", None)
        z = r.get("z", None)
        status = _tex_escape(str(r.get("status", "")))
        module_id = _tex_escape(str(r.get("module_id", "")))
        out_rel = str(r.get("output_rel", ""))
        out_cell = rf"\texttt{{{module_id}}}\\ \codepath{{{out_rel}}}"

        pred_cell = _tex_cell(pred)
        ref_cell = _ref_with_sigma(mean=mean, sigma=sigma) if mean is not None else "---"
        z_cell = _tex_cell(z)

        lines.append(rf"{obs} & {formula} & {pred_cell} & {ref_cell} & {z_cell} & {status} & {out_cell}\\")

    lines.append(r"\end{longtable}")
    lines.append("")

    write_tex.parent.mkdir(parents=True, exist_ok=True)
    write_tex.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return write_json, write_tex


def main() -> None:
    parser = argparse.ArgumentParser(description="Export a prediction ledger (JSON + TeX) from TFPT suite outputs.")
    parser.add_argument(
        "--write-json",
        default="out/predictions_ledger.json",
        help="Path (relative to tfpt-suite/) to write the JSON output (default: out/predictions_ledger.json).",
    )
    parser.add_argument(
        "--write-tex",
        default="out/predictions_ledger.tex",
        help="Path (relative to tfpt-suite/) to write the TeX output (default: out/predictions_ledger.tex).",
    )
    args = parser.parse_args()

    suite_root = _suite_root()
    out_json, out_tex = export_prediction_ledger_outputs(
        write_json=(suite_root / Path(str(args.write_json))).resolve(),
        write_tex=(suite_root / Path(str(args.write_tex))).resolve(),
    )
    print(str(out_json))
    print(str(out_tex))


if __name__ == "__main__":
    main()

