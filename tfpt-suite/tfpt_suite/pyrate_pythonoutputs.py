from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


def _workspace_root() -> Path:
    # .../wolfram_latex_attachments/tfpt-suite/tfpt_suite/<file>.py -> parents[2] is workspace root
    return Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class PyratePythonOutputModel:
    kind: str
    model_name_expected: str
    pythonoutput_dir: Path
    yaml_source: Path | None


def load_pyrate_pythonoutputs() -> dict[str, PyratePythonOutputModel]:
    """
    Load the central PyR@TE3 PythonOutput registry used by TFPT runners.

    File: `tfpt_suite/data/pyrate_pythonoutputs.json`
    """
    cfg_path = _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / "pyrate_pythonoutputs.json"
    raw = json.loads(cfg_path.read_text(encoding="utf-8"))
    models = raw.get("models", {}) if isinstance(raw, dict) else {}
    if not isinstance(models, dict) or not models:
        raise ValueError(f"Invalid pyrate_pythonoutputs.json (missing/empty 'models'): {cfg_path}")

    out: dict[str, PyratePythonOutputModel] = {}
    for kind, cfg in models.items():
        if not isinstance(cfg, dict):
            continue
        exp = str(cfg.get("model_name_expected", "")).strip()
        py_dir = str(cfg.get("pythonoutput_dir", "")).strip()
        if not exp or not py_dir:
            raise ValueError(f"Invalid entry for kind={kind!r} in {cfg_path}: need model_name_expected + pythonoutput_dir")
        yaml_src = str(cfg.get("yaml_source", "")).strip()
        out[str(kind)] = PyratePythonOutputModel(
            kind=str(kind),
            model_name_expected=exp,
            pythonoutput_dir=_workspace_root() / py_dir,
            yaml_source=(_workspace_root() / yaml_src) if yaml_src else None,
        )

    if not out:
        raise ValueError(f"No usable entries found in {cfg_path}")
    return out


def get_pyrate_pythonoutput(kind: str) -> PyratePythonOutputModel:
    models = load_pyrate_pythonoutputs()
    if kind not in models:
        raise KeyError(f"PyR@TE PythonOutput kind not found in config: {kind!r} (available: {sorted(models.keys())})")
    return models[kind]

