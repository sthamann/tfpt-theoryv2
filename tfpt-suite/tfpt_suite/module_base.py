from __future__ import annotations

import random
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from mpmath import mp

from tfpt_suite.config import SuiteConfig
from tfpt_suite.io import build_run_metadata, write_json, write_text
from tfpt_suite.rg_authority import rg_module_context


@dataclass(frozen=True)
class ModuleSpec:
    name: str
    module_id: str
    inputs: list[str]
    outputs: list[str]
    formulas: list[str]
    validation: list[str]
    determinism: str
    # Optional narrative fields (rendered in PDF when present).
    question: str | None = None
    objective: list[str] = field(default_factory=list)
    what_was_done: list[str] = field(default_factory=list)
    assumptions: list[str] = field(default_factory=list)
    gaps: list[str] = field(default_factory=list)
    references: list[str] = field(default_factory=list)
    maturity: str | None = None


@dataclass(frozen=True)
class Check:
    check_id: str
    passed: bool
    detail: str
    # Optional severity label (PASS/WARN/FAIL/INFO). If omitted, PASS/FAIL is inferred from `passed`.
    severity: str | None = None


SEVERITY_PASS = "PASS"
SEVERITY_WARN = "WARN"
SEVERITY_FAIL = "FAIL"
SEVERITY_INFO = "INFO"


def check_severity(c: Check) -> str:
    if c.severity:
        return str(c.severity)
    return SEVERITY_PASS if bool(c.passed) else SEVERITY_FAIL


def mk_check_pass(check_id: str, detail: str) -> Check:
    return Check(check_id=check_id, passed=True, detail=detail, severity=SEVERITY_PASS)


def mk_check_warn(check_id: str, detail: str) -> Check:
    return Check(check_id=check_id, passed=True, detail=detail, severity=SEVERITY_WARN)


def mk_check_fail(check_id: str, detail: str) -> Check:
    return Check(check_id=check_id, passed=False, detail=detail, severity=SEVERITY_FAIL)


def mk_check_info(check_id: str, detail: str) -> Check:
    return Check(check_id=check_id, passed=True, detail=detail, severity=SEVERITY_INFO)


@dataclass(frozen=True)
class ModuleResult:
    results: dict[str, Any]
    checks: list[Check]
    report: str
    warnings: list[str]


class TfptModule(ABC):
    module_id: str
    title: str

    def output_dir(self, config: SuiteConfig) -> Path:
        return config.output_dir / self.module_id

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[],
            outputs=[],
            formulas=[],
            validation=[],
            determinism="Deterministic given config seed + precision.",
        )

    @abstractmethod
    def run(self, config: SuiteConfig) -> ModuleResult:  # pragma: no cover
        raise NotImplementedError

    def run_and_write(self, *, config: SuiteConfig) -> ModuleResult:
        mp.dps = int(config.mp_dps)

        random.seed(int(config.seed))
        try:
            import numpy as np  # type: ignore

            np.random.seed(int(config.seed))
        except Exception:
            pass

        out_dir = self.output_dir(config)
        out_dir.mkdir(parents=True, exist_ok=True)

        with rg_module_context(self.module_id):
            result = self.run(config)

        meta = build_run_metadata(module_id=self.module_id, module_title=self.title, config=config)
        write_json(out_dir / "meta.json", meta, overwrite=config.overwrite)
        schema_version = 1
        plot_block = None
        if isinstance(result.results, dict) and "plot" in result.results:
            plot_block = result.results.get("plot")
        spec = self.spec()
        refs_from_results: list[str] = []
        if isinstance(result.results, dict):
            refs_block = result.results.get("references")
            if isinstance(refs_block, dict):
                refs_from_results.extend([str(k) for k in refs_block.keys()])
            elif isinstance(refs_block, list):
                refs_from_results.extend([str(k) for k in refs_block])
            ref_single = result.results.get("reference")
            if isinstance(ref_single, dict) and ref_single.get("dataset_id"):
                refs_from_results.append(str(ref_single.get("dataset_id")))
        if refs_from_results and not spec.references:
            spec = replace(spec, references=sorted(set(refs_from_results)))

        write_json(
            out_dir / "results.json",
            {
                "schema_version": schema_version,
                "module_id": self.module_id,
                "spec": spec,
                "results": result.results,
                "plot": plot_block,
                "checks": result.checks,
                "warnings": result.warnings,
            },
            overwrite=config.overwrite,
        )
        write_text(out_dir / "report.txt", result.report, overwrite=config.overwrite)

        return result

