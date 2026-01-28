#!/usr/bin/env python3
"""
E8 Cascade Model Runner
=======================
Single entry point for running the E8 Cascade RGE analysis.

Usage:
    python run_e8cascade.py [--gravity]
"""
import sys
from pathlib import Path

# Add src to Python path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from e8cascade.solver import E8CascadeSolver
from e8cascade.analysis import load_results, check_perturbativity, plot_stability_analysis


def main():
    """Main entry point"""
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(description='E8 Cascade RGE Solver')
    parser.add_argument('--gravity', action='store_true',
                        help='Enable gravity portal corrections')
    parser.add_argument('--cR-factor', type=float, default=1.0,
                        help='Multiply gravity coefficients by this factor (default: 1.0)')
    parser.add_argument('--output', type=str, default='results',
                        help='Output directory (default: results)')
    
    args = parser.parse_args()
    
    # Setup paths
    outdir = Path(args.output)
    outdir.mkdir(exist_ok=True)
    
    # Setup logging to file
    log_file = outdir / 'results.txt'
    
    class TeeOutput:
        """Output to both console and file"""
        def __init__(self, file_path):
            self.terminal = sys.stdout
            self.log = open(file_path, 'w')
            
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
            self.log.flush()
            
        def flush(self):
            self.terminal.flush()
            self.log.flush()
            
        def close(self):
            self.log.close()
    
    # Redirect stdout to both console and file
    tee = TeeOutput(log_file)
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = tee
    sys.stderr = tee  # Also capture errors
    
    # Initialize solver
    print("=" * 60)
    print("E8 CASCADE MODEL RGE SOLVER")
    print("=" * 60)
    print(f"Gravity portal: {'ENABLED' if args.gravity else 'DISABLED'}")
    if args.gravity and args.cR_factor != 1.0:
        print(f"Gravity factor: {args.cR_factor}×")
    print(f"Output directory: {outdir.resolve()}")
    print("=" * 60)
    
    solver = E8CascadeSolver(
        model_package="PythonOutput",          # Directory with PyR@TE output
        model_module="E8Cascade2LoopGravityV2",  # Module name
        results_dir=outdir,
        enable_gravity_portal=args.gravity
    )
    
    # Apply cR factor if specified
    if args.cR_factor != 1.0:
        solver.cR_factor = args.cR_factor
    
    # Solve RGEs
    print("\nStep 1: Solving RGEs...")
    solver.solve()
    
    # Generate plots
    print("\nStep 2: Generating plots...")
    solver.make_plots()
    
    # Additional analysis
    print("\nStep 3: Running analysis...")
    df = load_results(outdir)
    
    # Check perturbativity
    check_perturbativity(df)
    
    # Stability analysis
    plot_stability_analysis(df, outdir / "stability_analysis.png")
    
    print("\n" + "=" * 60)
    print(f"✓ Run complete! All results saved to: {outdir.resolve()}")
    print("=" * 60)
    
    # Print file list
    print("\nGenerated files:")
    for f in sorted(outdir.glob("*")):
        if f.is_file():
            print(f"  - {f.name}")
    
    # Restore original stdout/stderr and close log file
    sys.stdout = old_stdout
    sys.stderr = old_stderr
    tee.close()


if __name__ == "__main__":
    main() 