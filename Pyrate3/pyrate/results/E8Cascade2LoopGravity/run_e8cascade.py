#!/usr/bin/env python3
"""
E8 Cascade Model Runner (fixed)
==============================
Single entry point for running the E8 Cascade RGE analysis with
GUT-normalized hypercharge, explicit mass thresholds and optional
gravity portal.

Usage examples::

    # vanilla run with defaults
    python run_e8cascade.py

    # enable gravity portal with cR factor 2
    python run_e8cascade.py --gravity --cR-factor 2

    # turn on GUT normalisation and custom thresholds
    python run_e8cascade.py --gut-normalisation \
        --m-sigma 1e8 --m-nr 1e11
"""

import sys
import math
from pathlib import Path

# Add src to Python path so that "import e8cascade" works when launched from
# the project root.
sys.path.insert(0, str(Path(__file__).parent / "src"))

from e8cascade.solver_patch import SolverPatched as E8CascadeSolver
from e8cascade.analysis import (
    load_results,
    check_perturbativity,
    plot_stability_analysis,
)

###############################################################################
# Helpers
###############################################################################

class TeeOutput:
    """Duplicate stdout/stderr to console *and* a logfile."""

    def __init__(self, file_path: Path):
        self.terminal = sys.stdout
        self.log = open(file_path, "w", encoding="utf-8")

    def write(self, message: str):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):  # pragma: no cover – required for `python -u`
        self.terminal.flush()
        self.log.flush()

    def close(self):
        self.log.close()


###############################################################################
# Main entry point
###############################################################################

def main() -> None:  # noqa: D401 (imperative mood)
    """Command-line wrapper around :class:`~e8cascade.solver.E8CascadeSolver`."""

    # ------------------------------------------------------------------ CLI --
    import argparse

    parser = argparse.ArgumentParser(
        prog="run_e8cascade.py",
        description="E8 Cascade two-loop RGE solver with optional gravity portal",
    )

    parser.add_argument(
        "--gravity",
        action="store_true",
        help="Enable gravity portal corrections (disabled by default)",
    )
    parser.add_argument(
        "--cR-factor",
        type=float,
        default=1.0,
        metavar="FACTOR",
        help="Scale gravity coefficients by this factor (default 1.0)",
    )
    parser.add_argument(
        "--cR3-factor",
        type=float,
        default=1.0,
        metavar="FACTOR",
        help="Extra scaling for the g3 gravity term (default 1.0)",
    )
    parser.add_argument(
        "--gut-normalisation",
        action="store_true",
        help="Rescale g1 and β(g1) to GUT normalisation (√3⁄5, 3⁄5)",
    )
    parser.add_argument(
        "--m-sigma",
        type=float,
        default=1e9,
        metavar="GEV",
        help="Threshold MSigma in GeV (default 1e9)",
    )
    parser.add_argument(
        "--m-nr",
        type=float,
        default=1e12,
        metavar="GEV",
        help="Threshold MNR in GeV (default 1e12)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results",
        metavar="DIR",
        help="Output directory (default ./results)",
    )

    args = parser.parse_args()

    # ----------------------------------------------------------- filesystem --
    outdir = Path(args.output).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------- logging --
    log_file = outdir / "results.txt"
    tee = TeeOutput(log_file)
    old_stdout, old_stderr = sys.stdout, sys.stderr
    sys.stdout = tee  # type: ignore
    sys.stderr = tee  # type: ignore

    try:
        # ----------------------------------------------------------- banner --
        print("=" * 72)
        print("E8 CASCADE MODEL RGE SOLVER – FIXED EDITION")
        print("=" * 72)
        print(f"Gravity portal         : {'ENABLED' if args.gravity else 'disabled'}")
        if args.gravity:
            print(f"  ↳ cR  (g1,g2) factor : {args.cR_factor}")
            print(f"  ↳ cR3 (g3)   factor : {args.cR3_factor}")
        print(f"GUT normalisation      : {'ON' if args.gut_normalisation else 'off'}")
        print(f"Thresholds (GeV)       : MSigma = {args.m_sigma:.3g}, MNR = {args.m_nr:.3g}")
        print(f"Output directory       : {outdir}")
        print("=" * 72, flush=True)

        # --------------------------------------------------------- solver --
        solver = E8CascadeSolver(
            model_package="PythonOutput",          # directory with PyR@TE output
            model_module="E8Cascade2LoopGravity",  # module inside that dir
            results_dir=outdir,
            enable_gravity_portal=args.gravity,
        )

        # gravity scaling --------------------------------------------------
        if args.cR_factor != 1.0:
            if hasattr(solver, "cR_factor"):
                solver.cR_factor = args.cR_factor
            else:
                print("[warn] solver lacks attribute 'cR_factor' – ignoring cR-factor.")
        
        if hasattr(solver, "cR_factor"):
            solver.cR_factor  = args.cR_factor
        if hasattr(solver, "cR3_factor"):
            solver.cR3_factor = args.cR3_factor

        # GUT normalisation -----------------------------------------------
        if args.gut_normalisation:
            _apply_gut_normalisation(solver)

        # thresholds -------------------------------------------------------
        _apply_thresholds(solver, args.m_sigma, args.m_nr)

        # ------------------------------------------------------------ steps --
        print("\nStep 1 / 3  –  Solving RGEs …", flush=True)
        solver.solve()

        print("\nStep 2 / 3  –  Generating plots …", flush=True)
        solver.make_plots()

        print("\nStep 3 / 3  –  Post-analysis …", flush=True)
        df = load_results(outdir)
        check_perturbativity(df)
        plot_stability_analysis(df, outdir / "stability_analysis.png")

        # ----------------------------------------------------------- finish --
        print("\n" + "=" * 72)
        print("✓ Run complete! All results saved to:", outdir)
        print("=" * 72)

        print("\nGenerated files:")
        for f in sorted(outdir.glob("*")):
            if f.is_file():
                print("  •", f.name)

    finally:
        # ensure stdout/stderr are restored even if an exception bubbles up
        sys.stdout = old_stdout  # type: ignore
        sys.stderr = old_stderr  # type: ignore
        tee.close()


###############################################################################
# Implementation details
###############################################################################

def _apply_gut_normalisation(solver) -> None:
    """Rescale **g1** and **β(g1)** to the conventional GUT normalisation.

    If the solver provides a public convenience method ``rescale_hypercharge``
    it is used. Otherwise we try a best-effort direct modification of common
    attribute names (``gauge_couplings``, ``beta_funcs`` …). Fallback is a loud
    warning so the user can patch the solver class.
    """

    try:
        # preferred – built-in helper
        solver.rescale_hypercharge()  # type: ignore[attr-defined]
        print("[info] Hypercharge rescaled via solver.rescale_hypercharge().")
        return
    except AttributeError:
        pass  # we will attempt a manual patch

    # heuristic: search for g1 in typical containers
    done = False

    if hasattr(solver, "gauge_couplings"):
        gc = solver.gauge_couplings  # type: ignore[attr-defined]
        if isinstance(gc, dict) and "g1" in gc:
            gc["g1"] *= math.sqrt(3 / 5)
            done = True

    if hasattr(solver, "beta_funcs"):
        bf = solver.beta_funcs  # type: ignore[attr-defined]
        if isinstance(bf, dict) and "g1" in bf:
            if isinstance(bf["g1"], (float, int)):
                bf["g1"] *= 3 / 5
            elif isinstance(bf["g1"], (list, tuple)):
                bf["g1"] = [(3 / 5) * coeff for coeff in bf["g1"]]
            done = True

    if done:
        print("[info] Hypercharge rescaled manually (√3⁄5, 3⁄5).")
    else:
        print(
            "[warn] Could not locate g1 coupling/β-function – GUT normalisation "
            "*not* applied.\n       Patch your solver or rename attributes so that "
            "they contain 'g1'.",
            file=sys.stderr,
        )


def _apply_thresholds(solver, m_sigma: float, m_nr: float) -> None:
    """Attach threshold info to the solver.

    Solver API is not fixed – we try ``set_thresholds`` first, then fall back to
    a simple attribute. The heavy lifting (disabling β-functions below the
    scale) must be done inside the solver implementation. This function only
    injects the numbers so that *something* is there to read.
    """

    thresholds = {"MSigma": float(m_sigma), "MNR": float(m_nr)}

    if hasattr(solver, "set_thresholds"):
        solver.set_thresholds(thresholds)  # type: ignore[attr-defined]
        print("[info] Thresholds passed via solver.set_thresholds().")
    else:
        solver.thresholds = thresholds  # type: ignore[attr-defined]
        print("[info] Thresholds attached as solver.thresholds dict.")


###############################################################################
if __name__ == "__main__":
    main()
