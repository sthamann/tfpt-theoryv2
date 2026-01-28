from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable

from mpmath import mp


Z_TOL_DEFAULT = mp.mpf("2")
REL_ERROR_DEFAULT = mp.mpf("0.02")
LOG10_DEX_TOL_DEFAULT = mp.mpf("0.5")

def tfpt_suite_root() -> Path:
    # .../tfpt-suite/theoryv3/theoryv3_suite/utils.py -> parents[2] == tfpt-suite
    return Path(__file__).resolve().parents[2]


def theoryv3_root() -> Path:
    return Path(__file__).resolve().parents[1]


def workspace_root() -> Path:
    return Path(__file__).resolve().parents[3]


def read_json_if_exists(path: Path) -> dict[str, Any] | None:
    try:
        if path.is_file():
            return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return None


def load_tfpt_results(module_id: str, *, prefer_physics: bool = True) -> dict[str, Any] | None:
    root = tfpt_suite_root()
    out_physics = root / "out_physics" / module_id / "results.json"
    out_engineering = root / "out" / module_id / "results.json"

    if prefer_physics:
        for path in (out_physics, out_engineering):
            payload = read_json_if_exists(path)
            if payload is not None:
                return payload
    else:
        for path in (out_engineering, out_physics):
            payload = read_json_if_exists(path)
            if payload is not None:
                return payload
    return None


def extract_nested(payload: dict[str, Any], keys: Iterable[str], default: Any = None) -> Any:
    cur: Any = payload
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def coerce_float(value: Any, *, default: float = float("nan")) -> float:
    if isinstance(value, (int, float)):
        try:
            v = float(value)
            return v if math.isfinite(v) else float(default)
        except Exception:
            return float(default)
    if isinstance(value, str):
        s = value.strip()
        if not s:
            return float(default)
        try:
            v = float(s)
            return v if math.isfinite(v) else float(default)
        except Exception:
            return float(default)
    return float(default)


def coerce_mpf(value: Any, *, default: mp.mpf | None = None) -> mp.mpf:
    if default is None:
        default = mp.mpf("nan")
    if isinstance(value, mp.mpf):
        return value
    if isinstance(value, (int, float)):
        return mp.mpf(value)
    if isinstance(value, str):
        s = value.strip()
        if not s:
            return mp.mpf(default)
        try:
            return mp.mpf(s)
        except Exception:
            return mp.mpf(default)
    return mp.mpf(default)


def safe_log10(value: float) -> float:
    if value <= 0 or not math.isfinite(value):
        return float("nan")
    return float(math.log10(value))


def ensure_ascii(text: str) -> str:
    # Replace non-ASCII characters with '?'
    out = []
    for ch in str(text):
        if ord(ch) <= 127:
            out.append(ch)
        else:
            out.append("?")
    return "".join(out)
