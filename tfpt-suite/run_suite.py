from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

from tfpt_suite.config import SuiteConfig
from tfpt_suite.modules.registry import get_module_registry
from tfpt_suite.report_builder import build_tfpt_test_results_pdf
from tfpt_suite.export_reference_ledger import export_reference_ledger_tex
from tfpt_suite.export_prediction_ledger import export_prediction_ledger_outputs
from tfpt_suite.suite_manifest import build_suite_manifest, write_suite_manifest_outputs


def _parse_modules_arg(value: str) -> list[str]:
    value = value.strip()
    if not value:
        return []
    return [part.strip() for part in value.split(",") if part.strip()]


def _iter_selected_modules(all_ids: Iterable[str], selected: list[str] | None) -> list[str]:
    all_ids_list = list(all_ids)
    if not selected:
        return all_ids_list

    unknown = [m for m in selected if m not in all_ids_list]
    if unknown:
        raise SystemExit(f"Unknown module(s): {', '.join(unknown)}. Use `list-modules` to see available modules.")

    return selected


def main() -> None:
    parser = argparse.ArgumentParser(prog="tfpt-suite", description="TFPT modular verification suite (paper v2.4).")

    parser.add_argument(
        "--output-dir",
        default=str((Path(__file__).resolve().parent / "out").resolve()),
        help="Output directory for artifacts (default: tfpt-suite/out).",
    )
    parser.add_argument(
        "--mode",
        choices=["engineering", "physics"],
        default="engineering",
        help=(
            "Verification semantics (default: engineering). "
            "'engineering' emphasizes deterministic execution and explicit assumptions; "
            "'physics' upgrades large deviations / missing-derivation flags to WARN/FAIL."
        ),
    )
    parser.add_argument("--mp-dps", type=int, default=80, help="mpmath decimal digits of precision (default: 80).")
    parser.add_argument("--seed", type=int, default=0, help="Deterministic RNG seed (default: 0).")
    parser.add_argument("--no-plot", action="store_true", help="Disable plot generation.")
    parser.add_argument("--no-overwrite", action="store_true", help="Do not overwrite existing outputs.")

    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list-modules", help="List available modules.")

    run_all = sub.add_parser("run-all", help="Run all modules.")
    run_all.add_argument("--modules", default="", help="Comma-separated subset of modules to run.")

    run = sub.add_parser("run", help="Run a subset of modules.")
    run.add_argument("--modules", required=True, help="Comma-separated module IDs to run.")

    report = sub.add_parser("build-report", help="Build a central PDF report from module outputs.")
    report.add_argument(
        "--pdf-path",
        default=str((Path(__file__).resolve().parent / "tfpt-test-results.pdf").resolve()),
        help="Output PDF path (default: tfpt-suite/tfpt-test-results.pdf).",
    )

    manifest_cmd = sub.add_parser("build-manifest", help="Build suite manifest JSON + TeX summaries from out/ and out_physics/.")
    manifest_cmd.add_argument(
        "--out-roots",
        default="out,out_physics",
        help="Comma-separated output roots to scan (default: out,out_physics).",
    )
    manifest_cmd.add_argument(
        "--write-dir",
        default="out",
        help="Directory to write suite_manifest.{json,tex} outputs (default: out).",
    )
    manifest_cmd.add_argument(
        "--strict",
        action="store_true",
        help="Exit nonzero if any plot issues (missing/invalid) are detected.",
    )

    ref_cmd = sub.add_parser("export-reference-ledger", help="Export paper_reference_ledger.tex from tfpt_suite/data/references.json.")
    ref_cmd.add_argument(
        "--write-tex",
        default="out/paper_reference_ledger.tex",
        help="Path (relative to tfpt-suite/) to write the TeX reference ledger (default: out/paper_reference_ledger.tex).",
    )

    pred_cmd = sub.add_parser("export-prediction-ledger", help="Export predictions_ledger.{json,tex} from suite outputs (curated).")
    pred_cmd.add_argument(
        "--write-json",
        default="out/predictions_ledger.json",
        help="Path (relative to tfpt-suite/) to write the JSON prediction ledger (default: out/predictions_ledger.json).",
    )
    pred_cmd.add_argument(
        "--write-tex",
        default="out/predictions_ledger.tex",
        help="Path (relative to tfpt-suite/) to write the TeX prediction ledger (default: out/predictions_ledger.tex).",
    )

    args = parser.parse_args()

    registry = get_module_registry()
    all_ids = list(registry.keys())

    if args.command == "list-modules":
        for module_id in all_ids:
            module = registry[module_id]
            print(f"{module_id}: {module.title}")
        return

    if args.command == "build-report":
        verification_mode = getattr(args, "mode", "engineering")
        pdf_path = build_tfpt_test_results_pdf(
            out_dir=Path(args.output_dir), 
            pdf_path=Path(args.pdf_path),
            verification_mode=verification_mode
        )
        print(str(pdf_path))
        return

    if args.command == "build-manifest":
        suite_root = Path(__file__).resolve().parent
        roots = [(suite_root / Path(p.strip())).resolve() for p in str(args.out_roots).split(",") if p.strip()]
        manifest = build_suite_manifest(out_roots=roots)
        write_suite_manifest_outputs(manifest=manifest, out_dir=(suite_root / Path(str(args.write_dir))).resolve())
        totals = manifest.get("totals", {}) if isinstance(manifest.get("totals", {}), dict) else {}
        has_plot_issues = int(totals.get("plots_missing_total", 0)) > 0 or int(totals.get("plots_invalid_total", 0)) > 0
        if bool(getattr(args, "strict", False)) and has_plot_issues:
            raise SystemExit(2)
        return

    if args.command == "export-reference-ledger":
        suite_root = Path(__file__).resolve().parent
        out = export_reference_ledger_tex(write_path=(suite_root / Path(str(args.write_tex))).resolve())
        print(str(out))
        return

    if args.command == "export-prediction-ledger":
        suite_root = Path(__file__).resolve().parent
        out_json, out_tex = export_prediction_ledger_outputs(
            write_json=(suite_root / Path(str(args.write_json))).resolve(),
            write_tex=(suite_root / Path(str(args.write_tex))).resolve(),
        )
        print(str(out_json))
        print(str(out_tex))
        return

    selected = _parse_modules_arg(getattr(args, "modules", ""))
    module_ids = _iter_selected_modules(all_ids, selected)

    config = SuiteConfig(
        output_dir=Path(args.output_dir),
        mp_dps=args.mp_dps,
        seed=args.seed,
        overwrite=not args.no_overwrite,
        plot=not args.no_plot,
        verification_mode=str(args.mode),
    )

    for module_id in module_ids:
        module = registry[module_id]
        module.run_and_write(config=config)


if __name__ == "__main__":
    main()

