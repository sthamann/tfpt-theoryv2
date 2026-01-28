from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _schema_path() -> Path:
    return Path(__file__).resolve().parent / "schema" / "results.schema.json"


def _constant_schema_path() -> Path:
    return Path(__file__).resolve().parent / "schema" / "constant_ledger.schema.json"


def _load_schema() -> dict[str, Any]:
    path = _schema_path()
    return json.loads(path.read_text(encoding="utf-8"))


def _load_constant_schema() -> dict[str, Any] | None:
    path = _constant_schema_path()
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _basic_validate(payload: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    required = ["schema_version", "module_id", "spec", "results", "checks", "warnings"]
    for key in required:
        if key not in payload:
            errors.append(f"missing key: {key}")
    if payload.get("schema_version") != 1:
        errors.append(f"schema_version must be 1 (got {payload.get('schema_version')})")
    if not isinstance(payload.get("spec"), dict):
        errors.append("spec must be object")
    if not isinstance(payload.get("results"), dict):
        errors.append("results must be object")
    if not isinstance(payload.get("checks"), list):
        errors.append("checks must be list")
    if not isinstance(payload.get("warnings"), list):
        errors.append("warnings must be list")

    spec = payload.get("spec", {})
    if isinstance(spec, dict):
        for key in ("module_id", "inputs", "formulas", "determinism", "validation"):
            if key not in spec:
                errors.append(f"spec missing key: {key}")
    checks = payload.get("checks", [])
    if isinstance(checks, list):
        for idx, chk in enumerate(checks):
            if not isinstance(chk, dict):
                errors.append(f"checks[{idx}] must be object")
                continue
            for key in ("check_id", "passed", "detail"):
                if key not in chk:
                    errors.append(f"checks[{idx}] missing key: {key}")
    return errors


def _validate_with_jsonschema(payload: dict[str, Any], schema: dict[str, Any]) -> list[str] | None:
    try:
        import jsonschema  # type: ignore
    except Exception:
        return None
    try:
        jsonschema.validate(instance=payload, schema=schema)
        return []
    except Exception as exc:
        return [str(exc)]


def _basic_validate_constant(entry: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    required = ["name", "definition", "primitives", "fixed_point", "sensitivity", "value_pred", "status"]
    for key in required:
        if key not in entry:
            errors.append(f"missing key: {key}")
    if not isinstance(entry.get("primitives"), list):
        errors.append("primitives must be list")
    if not isinstance(entry.get("sensitivity"), dict):
        errors.append("sensitivity must be object")
    return errors


def _validate_constant_entry(entry: dict[str, Any], schema: dict[str, Any] | None) -> list[str]:
    if schema is None:
        return _basic_validate_constant(entry)
    schema_errors = _validate_with_jsonschema(entry, schema)
    if schema_errors is None:
        return _basic_validate_constant(entry)
    return schema_errors


def validate_artifacts(out_dir: Path) -> list[str]:
    schema = _load_schema()
    constant_schema = _load_constant_schema()
    errors: list[str] = []
    for child in sorted(out_dir.iterdir()):
        if not child.is_dir():
            continue
        results_path = child / "results.json"
        if not results_path.is_file():
            continue
        payload = json.loads(results_path.read_text(encoding="utf-8"))
        schema_errors = _validate_with_jsonschema(payload, schema)
        if schema_errors is None:
            schema_errors = _basic_validate(payload)
        if schema_errors:
            for err in schema_errors:
                errors.append(f"{results_path}: {err}")
        spec = payload.get("spec", {})
        spec_refs = spec.get("references", []) if isinstance(spec, dict) else []
        results_refs = payload.get("results", {}).get("references", {}) if isinstance(payload.get("results", {}), dict) else {}
        if isinstance(results_refs, dict) and results_refs and not spec_refs:
            errors.append(f"{results_path}: spec.references empty while results.references populated")
        if payload.get("module_id") == "constant_factory_audit":
            ledger = payload.get("results", {}).get("ledger", [])
            if not isinstance(ledger, list):
                errors.append(f"{results_path}: ledger must be list")
            else:
                for idx, entry in enumerate(ledger):
                    if not isinstance(entry, dict):
                        errors.append(f"{results_path}: ledger[{idx}] must be object")
                        continue
                    entry_errors = _validate_constant_entry(entry, constant_schema)
                    for err in entry_errors:
                        errors.append(f"{results_path}: ledger[{idx}] {err}")
            groups = payload.get("results", {}).get("groups", [])
            if isinstance(groups, list):
                for g_idx, group in enumerate(groups):
                    items = group.get("items", []) if isinstance(group, dict) else []
                    if not isinstance(items, list):
                        continue
                    for i_idx, item in enumerate(items):
                        if not isinstance(item, dict):
                            continue
                        if item.get("comparison") and not item.get("reference") and not item.get("reference_info"):
                            errors.append(
                                f"{results_path}: groups[{g_idx}].items[{i_idx}] comparison without reference"
                            )
    return errors
