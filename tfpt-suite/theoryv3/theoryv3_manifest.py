from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _suite_root() -> Path:
    # tfpt-suite/theoryv3/theoryv3_manifest.py -> tfpt-suite/
    return Path(__file__).resolve().parents[1]


def _repo_root() -> Path:
    return _suite_root().parent


def _relpath_from_suite(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(_suite_root()))
    except Exception:
        return str(path)


def _resolve_artifact_path(path_str: str) -> Path:
    p = Path(path_str)
    if p.is_absolute():
        return p
    s = str(path_str)
    if s.startswith("tfpt-suite/"):
        return (_repo_root() / p).resolve()
    return (_suite_root() / p).resolve()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _infer_check_severity(check: dict[str, Any]) -> str:
    sev = check.get("severity")
    if isinstance(sev, str) and sev.strip():
        return sev.strip().upper()
    return "PASS" if bool(check.get("passed")) else "FAIL"


def _counts_from_checks(checks: list[Any]) -> tuple[int, int, int, int]:
    total = 0
    passc = 0
    warnc = 0
    failc = 0
    for c in checks:
        if not isinstance(c, dict):
            continue
        total += 1
        sev = _infer_check_severity(c)
        if sev == "FAIL":
            failc += 1
        elif sev == "WARN":
            warnc += 1
        else:
            passc += 1
    return total, passc, warnc, failc


def _can_load_image(path: Path) -> bool:
    try:
        from PIL import Image  # type: ignore

        with Image.open(path) as img:
            img.verify()
        return True
    except Exception:
        pass

    try:
        import matplotlib.image as mpimg  # type: ignore

        _ = mpimg.imread(str(path))
        return True
    except Exception:
        return False


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


THEORYV3_SIGNAL: dict[str, str] = {
    "seed_invariants_audit": "Seed invariants (c3, φ0, δ_top, β)",
    "g5_origin_audit": "Holonomy degeneracy g=5",
    "g5_crosslink_audit": "Crosslink consistency (g=5 anchors)",
    "defect_partition_g5_audit": "Two-defect partition (δ2=(5/4)δ_top^2)",
    "alpha_backreaction_sensitivity_audit": "α backreaction sensitivity (k-grid)",
    "constant_factory_audit": "Constant factory audit (derived outputs)",
    "dark_energy_exponential_audit": "Dark energy exponential candidate",
    "dark_energy_norm_half_origin_audit": "Dark energy normalization n=1/2 scan",
    "flavor_pattern_audit": "Flavor pattern audit (texture invariants)",
    "pmns_tm1_audit": "PMNS TM1 audit",
    "yukawa_exponent_index_audit": "Yukawa exponent-index mapping",
    "yukawa_index_mapping_audit": "Yukawa index mapping",
    "axion_dm_audit": "Axion DM audit",
    "baryon_consistency_audit": "Baryon consistency / H0 audit",
}


@dataclass(frozen=True)
class Entry:
    module_id: str
    mode: str
    checks_total: int
    pass_count: int
    warn_count: int
    fail_count: int
    plot_expected: int
    plot_present: int
    plot_missing: int
    plot_invalid: int
    out_dir: str


def scan_theoryv3_out(out_root: Path) -> list[Entry]:
    out: list[Entry] = []
    if not out_root.exists():
        return out
    for child in sorted(out_root.iterdir(), key=lambda p: p.name):
        if not child.is_dir():
            continue
        results_path = child / "results.json"
        if not results_path.exists():
            continue
        payload = _load_json(results_path)
        meta_path = child / "meta.json"
        mode = "unknown"
        plot_enabled = True
        if meta_path.exists():
            meta = _load_json(meta_path)
            cfg = meta.get("config", {}) if isinstance(meta.get("config", {}), dict) else {}
            mode = str(cfg.get("verification_mode", mode))
            plot_enabled = bool(cfg.get("plot", True))
        checks = payload.get("checks", [])
        checks_total, passc, warnc, failc = _counts_from_checks(checks if isinstance(checks, list) else [])
        plot_block = payload.get("plot", {})
        plot_files: dict[str, str | None] = {}
        if isinstance(plot_block, dict):
            for k, v in plot_block.items():
                if isinstance(k, str):
                    plot_files[k] = v if isinstance(v, str) else None

        plot_expected = len(plot_files)
        plot_present = 0
        plot_missing = 0
        plot_invalid = 0
        if plot_enabled and plot_files:
            for _, path_str in plot_files.items():
                if not path_str:
                    plot_missing += 1
                    continue
                p = _resolve_artifact_path(path_str)
                if not p.exists():
                    plot_missing += 1
                    continue
                if p.suffix.lower() in (".png", ".jpg", ".jpeg"):
                    if _can_load_image(p):
                        plot_present += 1
                    else:
                        plot_invalid += 1
                else:
                    plot_present += 1

        out.append(
            Entry(
                module_id=str(payload.get("module_id", child.name)),
                mode=mode,
                checks_total=checks_total,
                pass_count=passc,
                warn_count=warnc,
                fail_count=failc,
                plot_expected=plot_expected,
                plot_present=plot_present,
                plot_missing=plot_missing,
                plot_invalid=plot_invalid,
                out_dir=_relpath_from_suite(child),
            )
        )
    return sorted(out, key=lambda e: e.module_id)


def write_theoryv3_manifest(*, entries: list[Entry], write_dir: Path) -> tuple[Path, Path, Path]:
    write_dir.mkdir(parents=True, exist_ok=True)
    json_path = write_dir / "theoryv3_manifest.json"
    tex_summary_path = write_dir / "theoryv3_manifest_summary.tex"
    tex_table_path = write_dir / "theoryv3_manifest_table.tex"

    totals = {
        "modules_total": len(entries),
        "checks_total": sum(e.checks_total for e in entries),
        "checks_pass_total": sum(e.pass_count for e in entries),
        "checks_warn_total": sum(e.warn_count for e in entries),
        "checks_fail_total": sum(e.fail_count for e in entries),
        "modules_with_warn": sum(1 for e in entries if e.warn_count > 0),
        "modules_with_fail": sum(1 for e in entries if e.fail_count > 0),
        "plots_expected_total": sum(e.plot_expected for e in entries),
        "plots_present_total": sum(e.plot_present for e in entries),
        "plots_missing_total": sum(e.plot_missing for e in entries),
        "plots_invalid_total": sum(e.plot_invalid for e in entries),
    }

    json_payload = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "out_root": _relpath_from_suite(write_dir),
        "totals": totals,
        "modules": [e.__dict__ for e in entries],
    }
    json_path.write_text(json.dumps(json_payload, indent=2, sort_keys=True), encoding="utf-8")

    # Summary TeX macros
    s: list[str] = []
    s.append("% AUTO-GENERATED by tfpt-suite/theoryv3/theoryv3_manifest.py — do not edit by hand.")
    s.append(f"% generated_at_utc: {json_payload['generated_at_utc']}")
    s.append("")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreeModulesTotal}}{{{int(totals['modules_total'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreeChecksTotal}}{{{int(totals['checks_total'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreeChecksPassTotal}}{{{int(totals['checks_pass_total'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreeChecksWarnTotal}}{{{int(totals['checks_warn_total'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreeChecksFailTotal}}{{{int(totals['checks_fail_total'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreeModulesWithWarn}}{{{int(totals['modules_with_warn'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreeModulesWithFail}}{{{int(totals['modules_with_fail'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreePlotsExpectedTotal}}{{{int(totals['plots_expected_total'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreePlotsPresentTotal}}{{{int(totals['plots_present_total'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreePlotsMissingTotal}}{{{int(totals['plots_missing_total'])}}}")
    s.append(f"\\newcommand{{\\TFPTTheoryVThreePlotsInvalidTotal}}{{{int(totals['plots_invalid_total'])}}}")
    s.append("")
    tex_summary_path.write_text("\n".join(s) + "\n", encoding="utf-8")

    # Evidence table (compact)
    t: list[str] = []
    t.append("% AUTO-GENERATED by tfpt-suite/theoryv3/theoryv3_manifest.py — do not edit by hand.")
    t.append("\\begin{table}[H]")
    t.append("\\centering")
    t.append("\\small")
    t.append("\\begin{tabular}{@{}lllcc@{}}")
    t.append("\\toprule")
    t.append("Module & Signal & Result & Checks & Plots OK\\\\")
    t.append("\\midrule")
    for e in entries:
        mid = _tex_escape(e.module_id)
        signal = _tex_escape(THEORYV3_SIGNAL.get(e.module_id, e.module_id))
        if e.fail_count > 0:
            result = "FAIL"
        elif e.warn_count > 0:
            result = "WARN"
        else:
            result = "PASS"
        checks = f"{e.pass_count}/{e.checks_total}"
        plots_ok = "yes" if (e.plot_missing == 0 and e.plot_invalid == 0) else "no"
        t.append(f"\\texttt{{{mid}}} & {signal} & {result} & {checks} & {plots_ok}\\\\")
    t.append("\\bottomrule")
    t.append("\\end{tabular}")
    t.append("\\caption{theoryv3 compressed evidence layer (auto-generated from \\,\\texttt{tfpt-suite/theoryv3/out/}).}")
    t.append("\\end{table}")
    tex_table_path.write_text("\n".join(t) + "\n", encoding="utf-8")

    return json_path, tex_summary_path, tex_table_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate theoryv3 manifest JSON + TeX evidence table.")
    parser.add_argument(
        "--out-root",
        default="theoryv3/out",
        help="theoryv3 output root (relative to tfpt-suite/) to scan (default: theoryv3/out).",
    )
    args = parser.parse_args()

    out_root = (_suite_root() / Path(str(args.out_root))).resolve()
    entries = scan_theoryv3_out(out_root)
    write_theoryv3_manifest(entries=entries, write_dir=out_root)


if __name__ == "__main__":
    main()

