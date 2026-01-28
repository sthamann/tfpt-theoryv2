from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _ledger_path() -> Path:
    return Path(__file__).resolve().parent / "data" / "references.json"


def load_reference_ledger() -> dict[str, Any]:
    path = _ledger_path()
    if not path.is_file():
        raise FileNotFoundError(f"Reference ledger not found: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def get_dataset(dataset_id: str) -> dict[str, Any]:
    ledger = load_reference_ledger()
    datasets = ledger.get("datasets", {})
    if not isinstance(datasets, dict):
        raise KeyError("Invalid reference ledger: datasets is not a dict")
    if dataset_id not in datasets:
        raise KeyError(f"Dataset not found in reference ledger: {dataset_id}")
    return datasets[dataset_id]


def get_value_sigma(dataset_id: str) -> tuple[float | None, float | None]:
    ds = get_dataset(dataset_id)
    value = ds.get("value", None)
    sigma = ds.get("sigma", None)
    return value, sigma
