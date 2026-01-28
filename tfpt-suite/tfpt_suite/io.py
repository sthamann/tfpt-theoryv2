from __future__ import annotations

import json
import platform
import sys
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _jsonify(obj: Any) -> Any:
    # dataclasses
    if is_dataclass(obj):
        return _jsonify(asdict(obj))

    # pathlib
    if isinstance(obj, Path):
        return str(obj)

    # mpmath (mpf/mpc) -> string to preserve precision
    if obj.__class__.__module__.startswith("mpmath"):
        try:
            return str(obj)
        except Exception:
            return repr(obj)

    # numpy scalars/arrays (avoid hard dependency at import time)
    if obj.__class__.__module__.startswith("numpy"):
        try:
            import numpy as np  # type: ignore

            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, (np.integer, np.floating)):
                return obj.item()
        except Exception:
            return repr(obj)

    # containers
    if isinstance(obj, dict):
        return {str(k): _jsonify(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonify(v) for v in obj]

    # basic types
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj

    return repr(obj)


def write_json(path: Path, payload: Any, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        return
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as f:
        json.dump(_jsonify(payload), f, indent=2, sort_keys=True)
        f.write("\n")


def write_text(path: Path, text: str, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        return
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


def build_run_metadata(*, module_id: str, module_title: str, config: Any) -> dict[str, Any]:
    return {
        "module": {"id": module_id, "title": module_title},
        "timestamp_utc": now_utc_iso(),
        "python": {"version": sys.version, "executable": sys.executable},
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "config": config,
    }

