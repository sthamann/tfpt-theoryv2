from __future__ import annotations

import argparse
import sys
import unittest
from pathlib import Path
from typing import Iterable

TFPT_SUITE_ROOT = Path(__file__).resolve().parents[1]
if str(TFPT_SUITE_ROOT) not in sys.path:
    sys.path.insert(0, str(TFPT_SUITE_ROOT))

from tfpt_suite.config import SuiteConfig
from theoryv3_suite.modules.registry import get_module_registry
from theoryv3_suite.runtime import apply_numeric_context
from report_builder import build_theoryv3_report
from validator import validate_artifacts


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
        raise SystemExit(f"Unknown module(s): {', '.join(unknown)}. Use list-modules to see available modules.")
    return selected


def _run_tests() -> int:
    suite = unittest.defaultTestLoader.discover(
        str(Path(__file__).resolve().parent / "tests"),
        pattern="test_*.py",
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


def _require_out_physics() -> Path:
    out_physics = TFPT_SUITE_ROOT / "out_physics"
    if not out_physics.is_dir():
        raise SystemExit(f"Missing required out_physics directory: {out_physics}")
    return out_physics


def main() -> None:
    parser = argparse.ArgumentParser(prog="theoryv3", description="TFPT theoryv3 analysis suite.")

    parser.add_argument(
        "--output-dir",
        default=str((Path(__file__).resolve().parent / "out").resolve()),
        help="Output directory for artifacts (default: tfpt-suite/theoryv3/out).",
    )
    parser.add_argument(
        "--mode",
        choices=["engineering", "physics"],
        default="engineering",
        help="Verification semantics (default: engineering).",
    )
    parser.add_argument("--mp-dps", type=int, default=80, help="mpmath decimal digits of precision (default: 80).")
    parser.add_argument("--seed", type=int, default=0, help="Deterministic RNG seed (default: 0).")
    parser.add_argument("--no-plot", action="store_true", help="Disable plot generation.")
    parser.add_argument("--no-overwrite", action="store_true", help="Do not overwrite existing outputs.")

    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list-modules", help="List available theoryv3 modules.")

    run_all = sub.add_parser("run-all", help="Run all theoryv3 modules.")
    run_all.add_argument("--modules", default="", help="Comma-separated subset of modules to run.")

    run = sub.add_parser("run", help="Run a subset of theoryv3 modules.")
    run.add_argument("--modules", required=True, help="Comma-separated module IDs to run.")

    report = sub.add_parser("build-report", help="Build theoryv3 PDF report from outputs.")
    report.add_argument(
        "--pdf-path",
        default=str((Path(__file__).resolve().parent / "theoryv3-analysis.pdf").resolve()),
        help="Output PDF path (default: tfpt-suite/theoryv3/theoryv3-analysis.pdf).",
    )

    sub.add_parser("run-tests", help="Run theoryv3 unit tests.")
    sub.add_parser("validate-artifacts", help="Validate theoryv3 results.json against schema.")

    args = parser.parse_args()

    registry = get_module_registry()
    all_ids = list(registry.keys())

    if args.command == "list-modules":
        for module_id in all_ids:
            module = registry[module_id]
            print(f"{module_id}: {module.title}")
        return

    if args.command == "build-report":
        _require_out_physics()
        pdf_path = build_theoryv3_report(out_dir=Path(args.output_dir), pdf_path=Path(args.pdf_path))
        print(str(pdf_path))
        return

    if args.command == "run-tests":
        sys.exit(_run_tests())

    if args.command == "validate-artifacts":
        _require_out_physics()
        errors = validate_artifacts(Path(args.output_dir))
        if errors:
            raise SystemExit("\n".join(errors))
        print("OK")
        return

    _require_out_physics()
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

    apply_numeric_context(config)
    for module_id in module_ids:
        module = registry[module_id]
        module.run_and_write(config=config)

    errors = validate_artifacts(Path(args.output_dir))
    if errors:
        raise SystemExit("\n".join(errors))


if __name__ == "__main__":
    main()
