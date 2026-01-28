from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def _workspace_root() -> Path:
    # tfpt-suite/tfpt_suite/suite_manifest.py -> tfpt-suite/
    #
    # Note: This script is typically executed from inside the `tfpt-suite/` folder
    # (where `tfpt_suite` is importable). Therefore we treat `tfpt-suite/` as the
    # workspace root so default paths like `out/` and `out_physics/` resolve
    # correctly and match the CLI runner defaults.
    return Path(__file__).resolve().parents[1]


def _repo_root() -> Path:
    # repo root is one level above tfpt-suite/
    return _workspace_root().parent


def _resolve_artifact_path(path_str: str) -> Path:
    """
    Resolve an artifact path emitted by modules.

    We accept three common styles:
    - absolute paths
    - paths relative to tfpt-suite/ (e.g. "out/module/file.png")
    - paths relative to repo root that include the "tfpt-suite/" prefix
      (e.g. "tfpt-suite/out_physics/module/file.png")
    """
    p = Path(path_str)
    if p.is_absolute():
        return p
    s = str(path_str)
    if s.startswith("tfpt-suite/"):
        return (_repo_root() / p).resolve()
    return (_workspace_root() / p).resolve()


def _relpath(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(_workspace_root()))
    except Exception:
        return str(path)


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _infer_check_severity(check: dict[str, Any]) -> str:
    sev = check.get("severity")
    if isinstance(sev, str) and sev.strip():
        return sev.strip().upper()
    return "PASS" if bool(check.get("passed")) else "FAIL"


def _counts_from_checks(checks: Iterable[dict[str, Any]]) -> tuple[int, int, int, int]:
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


_DATA_JSON_RE = re.compile(r"tfpt_suite/data/([A-Za-z0-9_.-]+\.json)")


def _extract_reference_keys(spec: dict[str, Any]) -> list[str]:
    keys: set[str] = set()
    for field in ("inputs", "references"):
        entries = spec.get(field, [])
        if not isinstance(entries, list):
            continue
        for s in entries:
            if not isinstance(s, str):
                continue
            for m in _DATA_JSON_RE.findall(s):
                keys.add(m)
            # Some modules mention reference files without the prefix.
            for known in ("global_reference.json", "global_reference_minimal.json", "references.json"):
                if known in s:
                    keys.add(known)
    return sorted(keys)


def _flatten_plot_files(plot_block: Any) -> dict[str, str | None]:
    if not isinstance(plot_block, dict):
        return {}
    out: dict[str, str | None] = {}
    for k, v in plot_block.items():
        if isinstance(k, str):
            out[k] = v if isinstance(v, str) else None
    return out


def _can_load_image(path: Path) -> bool:
    # Prefer Pillow if present; fall back to matplotlib.image if available.
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
    # Minimal escaping for safe inline TeX (we keep it conservative: module IDs, filenames, check IDs).
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


@dataclass(frozen=True)
class ModuleManifest:
    module_id: str
    mode: str
    out_dir: str
    checks_total: int
    pass_count: int
    warn_count: int
    fail_count: int
    plot_files: dict[str, str | None]
    plot_expected: int
    plot_present: int
    plot_missing: int
    plot_invalid: int
    reference_keys: list[str]
    warnings: list[str]


def scan_output_root(out_root: Path) -> list[ModuleManifest]:
    if not out_root.exists():
        return []

    manifests: list[ModuleManifest] = []
    for child in sorted(out_root.iterdir(), key=lambda p: p.name):
        if not child.is_dir():
            continue
        # Skip unconventional outputs inside conventional roots.
        if child.name == "unconventional":
            continue
        results_path = child / "results.json"
        if not results_path.exists():
            continue

        results = _load_json(results_path)
        module_id = str(results.get("module_id", child.name))
        checks = results.get("checks", [])
        checks_total, passc, warnc, failc = _counts_from_checks(checks if isinstance(checks, list) else [])
        spec = results.get("spec", {}) if isinstance(results.get("spec", {}), dict) else {}
        ref_keys = _extract_reference_keys(spec)
        warnings_list = results.get("warnings", [])
        warnings_strs = [str(w) for w in warnings_list] if isinstance(warnings_list, list) else []

        meta_path = child / "meta.json"
        mode = "unknown"
        plot_enabled = True
        if meta_path.exists():
            meta = _load_json(meta_path)
            cfg = meta.get("config", {}) if isinstance(meta.get("config", {}), dict) else {}
            mode = str(cfg.get("verification_mode", mode))
            plot_enabled = bool(cfg.get("plot", True))
        else:
            # fallback heuristic
            mode = "physics" if "out_physics" in str(out_root) else "engineering"

        plot_block = _flatten_plot_files(results.get("plot"))
        # Plot verification for declared plots in results.json
        plot_expected = len(plot_block)
        plot_present = 0
        plot_missing = 0
        plot_invalid = 0
        if plot_enabled and plot_block:
            for _, path_str in plot_block.items():
                if not path_str:
                    plot_missing += 1
                    continue
                p = _resolve_artifact_path(path_str)
                if not p.exists():
                    plot_missing += 1
                    continue
                # Only validate image load for raster images; PDF is treated as existence-only.
                if p.suffix.lower() in (".png", ".jpg", ".jpeg"):
                    if _can_load_image(p):
                        plot_present += 1
                    else:
                        plot_invalid += 1
                else:
                    plot_present += 1

        manifests.append(
            ModuleManifest(
                module_id=module_id,
                mode=mode,
                out_dir=_relpath(child),
                checks_total=checks_total,
                pass_count=passc,
                warn_count=warnc,
                fail_count=failc,
                plot_files=plot_block,
                plot_expected=plot_expected,
                plot_present=plot_present,
                plot_missing=plot_missing,
                plot_invalid=plot_invalid,
                reference_keys=ref_keys,
                warnings=warnings_strs,
            )
        )

    return manifests


def build_suite_manifest(*, out_roots: list[Path]) -> dict[str, Any]:
    modules: list[ModuleManifest] = []
    for r in out_roots:
        modules.extend(scan_output_root(r))

    # Sort for stable output.
    modules = sorted(modules, key=lambda m: (m.mode, m.module_id))

    ids_all = sorted({m.module_id for m in modules})
    ids_eng = sorted({m.module_id for m in modules if m.mode == "engineering"})
    ids_phy = sorted({m.module_id for m in modules if m.mode == "physics"})

    def _totals_for(subset: list[ModuleManifest]) -> dict[str, int]:
        return {
            "modules_total": len(subset),
            "checks_total": sum(m.checks_total for m in subset),
            "modules_with_warn": sum(1 for m in subset if m.warn_count > 0),
            "modules_with_fail": sum(1 for m in subset if m.fail_count > 0),
            "plots_expected_total": sum(m.plot_expected for m in subset),
            "plots_present_total": sum(m.plot_present for m in subset),
            "plots_missing_total": sum(m.plot_missing for m in subset),
            "plots_invalid_total": sum(m.plot_invalid for m in subset),
        }

    total_checks = sum(m.checks_total for m in modules)
    total_warn_modules = sum(1 for m in modules if m.warn_count > 0)
    total_fail_modules = sum(1 for m in modules if m.fail_count > 0)
    total_plot_expected = sum(m.plot_expected for m in modules)
    total_plot_present = sum(m.plot_present for m in modules)
    total_plot_missing = sum(m.plot_missing for m in modules)
    total_plot_invalid = sum(m.plot_invalid for m in modules)

    return {
        "schema_version": 1,
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "out_roots": [_relpath(p) for p in out_roots],
        "totals": {
            "modules_total": len(modules),
            "module_ids_unique_total": len(ids_all),
            "module_ids_unique_engineering": len(ids_eng),
            "module_ids_unique_physics": len(ids_phy),
            "checks_total": total_checks,
            "modules_with_warn": total_warn_modules,
            "modules_with_fail": total_fail_modules,
            "plots_expected_total": total_plot_expected,
            "plots_present_total": total_plot_present,
            "plots_missing_total": total_plot_missing,
            "plots_invalid_total": total_plot_invalid,
            "engineering": _totals_for([m for m in modules if m.mode == "engineering"]),
            "physics": _totals_for([m for m in modules if m.mode == "physics"]),
        },
        "modules": [
            {
                "module_id": m.module_id,
                "mode": m.mode,
                "out_dir": m.out_dir,
                "checks_total": m.checks_total,
                "pass_count": m.pass_count,
                "warn_count": m.warn_count,
                "fail_count": m.fail_count,
                "plot_files": {k: (_relpath(_resolve_artifact_path(v)) if v else None) for k, v in m.plot_files.items()},
                "plot_expected": m.plot_expected,
                "plot_present": m.plot_present,
                "plot_missing": m.plot_missing,
                "plot_invalid": m.plot_invalid,
                "reference_keys": m.reference_keys,
            }
            for m in modules
        ],
    }


def write_suite_manifest_outputs(*, manifest: dict[str, Any], out_dir: Path) -> tuple[Path, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "suite_manifest.json"
    tex_summary_path = out_dir / "suite_manifest_summary.tex"
    tex_table_path = out_dir / "suite_manifest_table.tex"

    json_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    totals = manifest.get("totals", {}) if isinstance(manifest.get("totals", {}), dict) else {}
    modules = manifest.get("modules", []) if isinstance(manifest.get("modules", []), list) else []

    # Build problem lists (WARN/FAIL and plot issues).
    warn_mods: list[dict[str, Any]] = [m for m in modules if isinstance(m, dict) and int(m.get("warn_count", 0)) > 0]
    fail_mods: list[dict[str, Any]] = [m for m in modules if isinstance(m, dict) and int(m.get("fail_count", 0)) > 0]
    plot_bad: list[dict[str, Any]] = [
        m
        for m in modules
        if isinstance(m, dict) and (int(m.get("plot_missing", 0)) > 0 or int(m.get("plot_invalid", 0)) > 0)
    ]

    def _row_module(m: dict[str, Any]) -> str:
        mid = _tex_escape(str(m.get("module_id", "")))
        mode = _tex_escape(str(m.get("mode", "")))
        outp = _tex_escape(str(m.get("out_dir", "")))
        ctot = int(m.get("checks_total", 0))
        w = int(m.get("warn_count", 0))
        f = int(m.get("fail_count", 0))
        return f"\\item \\texttt{{{mid}}} ({mode}): checks={ctot}, WARN={w}, FAIL={f} (out: \\texttt{{{outp}}})"

    lines: list[str] = []
    lines.append("% AUTO-GENERATED by tfpt_suite.suite_manifest.py — do not edit by hand.")
    lines.append(f"% generated_at_utc: {_tex_escape(str(manifest.get('generated_at_utc','')))}")
    lines.append("")
    lines.append(f"\\newcommand{{\\TFPTSuiteModulesTotal}}{{{int(totals.get('modules_total', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuiteUniqueModulesTotal}}{{{int(totals.get('module_ids_unique_total', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuiteUniqueModulesEngineering}}{{{int(totals.get('module_ids_unique_engineering', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuiteUniqueModulesPhysics}}{{{int(totals.get('module_ids_unique_physics', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuiteChecksTotal}}{{{int(totals.get('checks_total', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuiteModulesWithWarn}}{{{int(totals.get('modules_with_warn', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuiteModulesWithFail}}{{{int(totals.get('modules_with_fail', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuitePlotsExpectedTotal}}{{{int(totals.get('plots_expected_total', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuitePlotsPresentTotal}}{{{int(totals.get('plots_present_total', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuitePlotsMissingTotal}}{{{int(totals.get('plots_missing_total', 0))}}}")
    lines.append(f"\\newcommand{{\\TFPTSuitePlotsInvalidTotal}}{{{int(totals.get('plots_invalid_total', 0))}}}")
    lines.append("")
    lines.append("% Problem modules (WARN/FAIL) and plot inventory issues:")
    lines.append("\\newcommand{\\TFPTSuiteProblemModules}{%")
    lines.append("\\begin{itemize}")
    if fail_mods:
        lines.append("\\item \\textbf{FAIL modules:}")
        lines.append("\\begin{itemize}")
        for m in fail_mods:
            lines.append(_row_module(m))
        lines.append("\\end{itemize}")
    if warn_mods:
        lines.append("\\item \\textbf{WARN modules:}")
        lines.append("\\begin{itemize}")
        for m in warn_mods:
            lines.append(_row_module(m))
        lines.append("\\end{itemize}")
    if plot_bad:
        lines.append("\\item \\textbf{Plot issues (missing/invalid):}")
        lines.append("\\begin{itemize}")
        for m in plot_bad:
            mid = _tex_escape(str(m.get("module_id", "")))
            mode = _tex_escape(str(m.get("mode", "")))
            miss = int(m.get("plot_missing", 0))
            inv = int(m.get("plot_invalid", 0))
            outp = _tex_escape(str(m.get("out_dir", "")))
            lines.append(f"\\item \\texttt{{{mid}}} ({mode}): missing={miss}, invalid={inv} (out: \\texttt{{{outp}}})")
        lines.append("\\end{itemize}")
    if not (fail_mods or warn_mods or plot_bad):
        lines.append("\\item \\textbf{No WARN/FAIL modules and no plot issues detected.}")
    lines.append("\\end{itemize}")
    lines.append("}%")
    lines.append("")
    tex_summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Optional compact table (only WARN/FAIL or plot-issue modules; keeps the paper readable).
    table_lines: list[str] = []
    table_lines.append("% AUTO-GENERATED by tfpt_suite.suite_manifest.py — do not edit by hand.")
    table_lines.append("\\begin{table}[H]")
    table_lines.append("\\centering")
    table_lines.append("\\small")
    table_lines.append("\\begin{tabular}{@{}llrrrrrr@{}}")
    table_lines.append("\\toprule")
    table_lines.append("Module & Mode & Checks & PASS & WARN & FAIL & Plots(exp) & Plots(ok)\\\\")
    table_lines.append("\\midrule")
    focus = sorted({(m.get("module_id", ""), m.get("mode", "")) for m in (warn_mods + fail_mods + plot_bad)}, key=lambda x: (x[1], x[0]))
    idx = {(m.get("module_id", ""), m.get("mode", "")): m for m in modules if isinstance(m, dict)}
    for key in focus:
        m = idx.get(key)
        if not m:
            continue
        mid = _tex_escape(str(m.get("module_id", "")))
        mode = _tex_escape(str(m.get("mode", "")))
        ctot = int(m.get("checks_total", 0))
        pc = int(m.get("pass_count", 0))
        wc = int(m.get("warn_count", 0))
        fc = int(m.get("fail_count", 0))
        pexp = int(m.get("plot_expected", 0))
        pok = int(m.get("plot_present", 0))
        table_lines.append(f"\\texttt{{{mid}}} & {mode} & {ctot} & {pc} & {wc} & {fc} & {pexp} & {pok}\\\\")
    if not focus:
        table_lines.append("\\multicolumn{8}{l}{No problem modules to report.}\\\\")
    table_lines.append("\\bottomrule")
    table_lines.append("\\end{tabular}")
    table_lines.append("\\caption{Suite manifest: compact table of WARN/FAIL modules and plot issues (auto-generated).}")
    table_lines.append("\\end{table}")
    tex_table_path.write_text("\n".join(table_lines) + "\n", encoding="utf-8")

    return json_path, tex_summary_path, tex_table_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a conventional-suite manifest from out/ and out_physics/ artifacts.")
    parser.add_argument(
        "--out-roots",
        default="out,out_physics",
        help="Comma-separated output roots to scan (default: out,out_physics).",
    )
    parser.add_argument(
        "--write-dir",
        default="out",
        help="Directory to write suite_manifest.{json,tex} outputs (default: out).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit nonzero if any plot issues (missing/invalid) are detected.",
    )
    args = parser.parse_args()

    roots = [(_workspace_root() / Path(p.strip())).resolve() for p in str(args.out_roots).split(",") if p.strip()]
    manifest = build_suite_manifest(out_roots=roots)
    _, _, _ = write_suite_manifest_outputs(manifest=manifest, out_dir=(_workspace_root() / Path(str(args.write_dir))).resolve())

    totals = manifest.get("totals", {}) if isinstance(manifest.get("totals", {}), dict) else {}
    has_plot_issues = int(totals.get("plots_missing_total", 0)) > 0 or int(totals.get("plots_invalid_total", 0)) > 0
    if args.strict and has_plot_issues:
        raise SystemExit(2)


if __name__ == "__main__":
    main()

