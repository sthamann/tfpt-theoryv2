from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def resolve_axion_claim(*, ax_raw: dict[str, Any], output_dir: Path | None = None) -> dict[str, Any]:
    claim = ax_raw.get("tfpt_claim", {}) if isinstance(ax_raw.get("tfpt_claim", {}), dict) else {}
    f_a = claim.get("f_a_GeV", float("nan"))
    m_a = claim.get("m_a_micro_eV", float("nan"))
    nu = claim.get("frequency_GHz_claim", float("nan"))
    source = "quoted"
    derived_path = None

    if output_dir is not None:
        candidate = Path(output_dir) / "axion_fa_derivation" / "results.json"
        if candidate.is_file():
            try:
                payload = json.loads(candidate.read_text(encoding="utf-8"))
                derived = payload.get("results", {}).get("derived", {}) if isinstance(payload, dict) else {}
                if isinstance(derived, dict) and derived.get("f_a_GeV", None) is not None:
                    f_a = float(derived.get("f_a_GeV", f_a))
                    if derived.get("m_a_micro_eV", None) is not None:
                        m_a = float(derived.get("m_a_micro_eV", m_a))
                    source = "derived"
                    derived_path = str(candidate)
            except Exception:
                pass

    return {
        "f_a_GeV": f_a,
        "m_a_micro_eV": m_a,
        "frequency_GHz_claim": nu,
        "source": source,
        "derived_path": derived_path,
    }
