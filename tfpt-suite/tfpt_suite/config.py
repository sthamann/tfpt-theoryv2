from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SuiteConfig:
    output_dir: Path
    mp_dps: int = 80
    seed: int = 0
    overwrite: bool = True
    plot: bool = True
    # Verification semantics:
    # - engineering: "plumbing" + deterministic calculations; avoid false negatives from known-conditional closures
    # - physics:     treat large deviations / missing-derivation flags as FAIL/WARN to prevent "everything green" optics
    verification_mode: str = "engineering"

