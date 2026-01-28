from __future__ import annotations

import json
import os
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from pathlib import Path


def _workspace_root() -> Path:
    # .../wolfram_latex_attachments/tfpt-suite/tfpt_suite/<file>.py -> parents[2] is workspace root
    return Path(__file__).resolve().parents[2]


def _data_path(name: str) -> Path:
    return _workspace_root() / "tfpt-suite" / "tfpt_suite" / "data" / name


@dataclass(frozen=True)
class RgAuthorityPolicy:
    above_MZ: str
    models: dict[str, str]
    below_MZ: str
    allowed_1loop_modules: list[str]
    note: str | None = None


_POLICY_CACHE: RgAuthorityPolicy | None = None


def _policy_cached() -> RgAuthorityPolicy:
    global _POLICY_CACHE
    if _POLICY_CACHE is None:
        _POLICY_CACHE = load_rg_authority_policy()
    return _POLICY_CACHE


_CURRENT_MODULE_ID: ContextVar[str | None] = ContextVar("tfpt_rg_current_module_id", default=None)


def current_module_id() -> str | None:
    return _CURRENT_MODULE_ID.get()


@contextmanager
def rg_module_context(module_id: str):
    """
    Context marker for RG-authority checks.

    The suite sets this while a module executes so legacy helpers can be permitted for explicitly
    whitelisted "below-MZ EFT" modules (see `rg_authority.json`).
    """
    token = _CURRENT_MODULE_ID.set(str(module_id))
    try:
        yield
    finally:
        _CURRENT_MODULE_ID.reset(token)


def load_rg_authority_policy() -> RgAuthorityPolicy:
    path = _data_path("rg_authority.json")
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError(f"Invalid rg_authority.json (expected object): {path}")
    above = str(raw.get("above_MZ", "")).strip()
    below = str(raw.get("below_MZ", "")).strip()
    models = raw.get("models", {}) if isinstance(raw.get("models", {}), dict) else {}
    allowed = raw.get("allowed_1loop_modules", []) if isinstance(raw.get("allowed_1loop_modules", []), list) else []
    note = raw.get("note", None)
    if not above or not below or not models:
        raise ValueError(f"Invalid rg_authority.json (missing above_MZ/below_MZ/models): {path}")
    return RgAuthorityPolicy(
        above_MZ=above,
        models={str(k): str(v) for k, v in dict(models).items()},
        below_MZ=below,
        allowed_1loop_modules=[str(x) for x in allowed],
        note=str(note) if isinstance(note, str) else None,
    )


def mz_scale_GeV(*, fallback: float = 91.1876) -> float:
    """
    Canonical MZ scale used by the RG-authority policy.

    We read the shared suite input `sm_inputs_mz.json` so the policy tracks the suite’s chosen boundary scale.
    """
    path = _data_path("sm_inputs_mz.json")
    if not path.exists():
        return float(fallback)
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
        return float(raw.get("mu_GeV", fallback))
    except Exception:
        return float(fallback)


def legacy_sm_rge_above_mz_is_allowed() -> bool:
    """
    Escape hatch for debugging only.

    The TFPT task policy requires PyR@TE betas as the single source of truth above MZ; legacy SM RG helpers
    (`rge_sm.py`) are blocked above MZ by default. Set this env var to *temporarily* allow legacy helpers:
      TFPT_ALLOW_LEGACY_SM_RGE_ABOVE_MZ=1
    """
    return str(os.environ.get("TFPT_ALLOW_LEGACY_SM_RGE_ABOVE_MZ", "")).strip() in ("1", "true", "TRUE", "yes", "YES")


def enforce_no_legacy_sm_rge_above_mz(*, mu_start_GeV: float, mu_end_GeV: float, caller: str) -> None:
    """
    Guardrail: legacy SM RG helpers must not be used above MZ.
    """
    if legacy_sm_rge_above_mz_is_allowed():
        return
    pol = _policy_cached()
    mz = mz_scale_GeV()
    if float(mu_start_GeV) > mz or float(mu_end_GeV) > mz:
        mod = current_module_id()
        if mod is not None and str(mod) in set(pol.allowed_1loop_modules):
            return
        raise RuntimeError(
            "RG authority violation: legacy SM RG helper used above MZ. "
            f"caller={caller}, module_id={mod!r}, mu_start_GeV={mu_start_GeV}, mu_end_GeV={mu_end_GeV}, MZ={mz}. "
            "Use the PyR@TE-driven runner (tfpt_suite/rge_pyrate_2loop.py) or the MZ→mt boundary helper "
            "(tfpt_suite/pyrate_boundary_runner.py). For debugging only, set TFPT_ALLOW_LEGACY_SM_RGE_ABOVE_MZ=1."
        )

