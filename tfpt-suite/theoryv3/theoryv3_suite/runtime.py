from __future__ import annotations

import random
from typing import Any

from mpmath import mp


def apply_numeric_context(config: Any) -> None:
    mp.dps = int(getattr(config, "mp_dps", 80))
    seed = int(getattr(config, "seed", 0))
    random.seed(seed)
    try:
        import numpy as np  # type: ignore

        np.random.seed(seed)
    except Exception:
        pass
