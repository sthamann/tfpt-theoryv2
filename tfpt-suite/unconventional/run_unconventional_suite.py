from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable


# Ensure `tfpt_suite` is importable when running from repo root.
TFPT_SUITE_DIR = Path(__file__).resolve().parents[1]
if str(TFPT_SUITE_DIR) not in sys.path:
    sys.path.insert(0, str(TFPT_SUITE_DIR))


from tfpt_suite.config import SuiteConfig  # noqa: E402
from tfpt_suite.report_builder import build_tfpt_test_results_pdf  # noqa: E402
from tfpt_unconventional.modules.registry import get_unconventional_module_registry  # noqa: E402


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
        raise SystemExit(
            f"Unknown module(s): {', '.join(unknown)}. Use `list-modules` to see available unconventional modules."
        )

    return selected


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="tfpt-unconventional-suite",
        description="TFPT unconventional auxiliary suite (search/audit tooling; see tfpt-suite/unconventional/README.md).",
    )

    parser.add_argument(
        "--output-dir",
        default=str((TFPT_SUITE_DIR / "out" / "unconventional").resolve()),
        help="Output directory for unconventional artifacts (default: tfpt-suite/out/unconventional).",
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
    sub.add_parser("list-modules", help="List available unconventional modules.")

    run_all = sub.add_parser("run-all", help="Run all unconventional modules.")
    run_all.add_argument("--modules", default="", help="Comma-separated subset of modules to run.")

    run = sub.add_parser("run", help="Run a subset of unconventional modules.")
    run.add_argument("--modules", required=True, help="Comma-separated module IDs to run.")

    report = sub.add_parser("build-report", help="Build a central PDF report from unconventional module outputs.")
    report.add_argument(
        "--pdf-path",
        default=str((Path(__file__).resolve().parent / "unconventional-test-results.pdf").resolve()),
        help="Output PDF path (default: tfpt-suite/unconventional/unconventional-test-results.pdf).",
    )

    args = parser.parse_args()

    registry = get_unconventional_module_registry()
    all_ids = list(registry.keys())

    if args.command == "list-modules":
        for module_id in all_ids:
            module = registry[module_id]
            print(f"{module_id}: {module.title}")
        return

    if args.command == "build-report":
        pdf_path = build_tfpt_test_results_pdf(out_dir=Path(args.output_dir), pdf_path=Path(args.pdf_path))
        print(str(pdf_path))
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

