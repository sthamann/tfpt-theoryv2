#!/usr/bin/env python3
"""Quick grid-scan for E8-Cascade gauge-unification.

Passt sich automatisch an deine CSV-Header an ("mu_GeV", "log10_mu",
"g1_GUT", …). Misst den Gauge-Spread bei µ⋆ (Default 1×10¹⁶ GeV).
"""

from __future__ import annotations

import argparse
import math
import sys
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

CURRENT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CURRENT_DIR / "src"))

from e8cascade.solver_patch import SolverPatched as E8CascadeSolver
from e8cascade.analysis import load_results  # type: ignore

###############################################################################
# Helpers
###############################################################################

MU_ALIASES = [
    ("mu_GeV", lambda s: s),
    ("mu", lambda s: s),
    ("scale", lambda s: s),
    ("log10_mu", lambda s: 10 ** s),
    ("log10mu", lambda s: 10 ** s),
]

G1_ALIASES = ["g1_GUT", "g1_SM", "g1", "g_1", "alpha1", "a1", "α1"]
G2_ALIASES = ["g2", "g_2", "alpha2", "a2", "α2"]
G3_ALIASES = ["g3", "g_3", "alpha3", "a3", "α3"]


def logspace(min_val: float, max_val: float, steps: int) -> list[float]:
    return list(np.logspace(math.log10(min_val), math.log10(max_val), int(steps), endpoint=True))


def pick_column(df: pd.DataFrame, candidates: list[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"Keine passende Spalte gefunden für {candidates}")


def get_mu_series(df: pd.DataFrame) -> pd.Series:
    for col, fn in MU_ALIASES:
        if col in df.columns:
            return fn(df[col].astype(float))
    raise KeyError("µ-Spalte nicht gefunden (mu_GeV/mu/scale/log10_mu)")


def spread_at(df: pd.DataFrame, mu_star: float) -> float:
    mu_series = get_mu_series(df)
    idx = (mu_series - mu_star).abs().idxmin()

    g1_col = pick_column(df, G1_ALIASES)
    g2_col = pick_column(df, G2_ALIASES)
    g3_col = pick_column(df, G3_ALIASES)

    alphas = df.loc[idx, [g1_col, g2_col, g3_col]]
    return max(alphas) - min(alphas)

###############################################################################
# Main
###############################################################################

def main() -> None:  # noqa: D401
    p = argparse.ArgumentParser(description="Grid-Scan für Gauge-Unifikation")
    p.add_argument("--gut-normalisation", action="store_true")
    p.add_argument("--gravity", action="store_true")
    p.add_argument("--m-sigma", nargs=3, type=float, default=[1e7, 1e10, 4])
    p.add_argument("--m-nr",    nargs=3, type=float, default=[1e11, 1e14, 4])
    #  ▸▸ cR  – beliebig viele Einzelwerte  (wie bisher)
    p.add_argument("--cR",  nargs="*", type=float, default=[1.0],
                   help="Liste einzelner cR-Faktoren (g1,g2)")
    #  ▸▸ cR3 – entweder drei Zahlen MIN MAX N  ➜ log-grid
    #           oder beliebig viele Einzelwerte wie bei --cR
    p.add_argument("--cR3", nargs="+", type=float, metavar=("…"),
                   default=[1.0],
                   help="Einzelwerte   »0.3 1 3«  oder MIN MAX N  »0.3 3 5«")
    p.add_argument("--mu-star", type=float, default=1e16)
    p.add_argument("--out",     type=str, default="scan_results")
    args = p.parse_args()

    m_sigma_grid = logspace(*args.m_sigma)
    m_nr_grid    = logspace(*args.m_nr)
    # -------- cR-Raster -----------------------------------------------
    cR_grid = list(map(float, args.cR))

    # Wenn genau drei Zahlen → als Log-Gitter interpretieren
    if len(args.cR3) == 3:
        cR3_min, cR3_max, n_cR3 = args.cR3
        cR3_grid = logspace(cR3_min, cR3_max, n_cR3)
    else:
        cR3_grid = list(map(float, args.cR3))

    outdir = Path(args.out).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []

    for MSigma, MNR, cR, cR3 in product(m_sigma_grid, m_nr_grid,cR_grid,cR3_grid):
        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            solver = E8CascadeSolver(
                model_package="PythonOutput",
                model_module="E8Cascade2LoopGravity",
                results_dir=tmp_path,
                enable_gravity_portal=args.gravity,
            )
            if hasattr(solver, "cR_factor"):
                solver.cR_factor = cR
            if hasattr(solver, "cR3_factor"):
                solver.cR3_factor = cR3
            if args.gut_normalisation and hasattr(solver, "rescale_hypercharge"):
                solver.rescale_hypercharge()
            thresh = {"MSigma": float(MSigma), "MNR": float(MNR)}
            if hasattr(solver, "set_thresholds"):
                solver.set_thresholds(thresh)
            else:
                solver.thresholds = thresh  # type: ignore[attr-defined]

            print(f"• Run MSigma={MSigma:.2e}, MNR={MNR:.2e}, "f"cR={cR}, cR3={cR3}")
            solver.solve()

            df = load_results(tmp_path)
            try:
                spr = spread_at(df, args.mu_star)
            except KeyError as e:
                print(f"[warn] {e} – überspringe Kombination.")
                continue
            rows.append(dict(MSigma=MSigma, MNR=MNR,
                             cR=cR, cR3=cR3, spread=spr))

    if not rows:
        print("Keine gültigen Datensätze – bitte Header prüfen!")
        sys.exit(1)

    res = pd.DataFrame(rows).sort_values("spread", ignore_index=True)
    csv_path = outdir / "best_points.csv"
    res.to_csv(csv_path, index=False)
    print("\n✓ Scan fertig – Top 10:")
    print(res.head(10).to_string(index=False))
    print(f"\nVollständige Tabelle: {csv_path}\n")


if __name__ == "__main__":
    main()
