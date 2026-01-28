#!/usr/bin/env python3
"""
E8 Cascade Model Runner - 2-Loop Version
==========================================

Runs the E8 Cascade gauge unification model with 2-loop RGEs.
Generates plots and saves results as CSV files.

Usage:
    python run_e8cascade.py

Output:
    - results/gauge_couplings.csv
    - results/gauge_running.png
    - results/unification.png
    - results/stability_analysis.png
"""

import sys
import os
from pathlib import Path

def main():
    print("=" * 60)
    print("E8 Cascade Model - 2-Loop Analysis")
    print("=" * 60)
    
    # Setup paths
    current_dir = Path(__file__).parent
    results_dir = current_dir / "results"
    pyrate_dir = current_dir / "E8Cascade2LoopGravity" / "PythonOutput"
    
    # Add PyR@TE output to path
    sys.path.insert(0, str(pyrate_dir))
    
    try:
        from e8cascade.solver import E8CascadeSolver
        from e8cascade.analysis import plot_stability_analysis, check_perturbativity
        
        print(f"✓ Using PyR@TE output from: {pyrate_dir}")
        print(f"✓ Results will be saved to: {results_dir}")
        
        # Initialize solver
        solver = E8CascadeSolver(
            model_package=str(pyrate_dir),
            model_module="E8Cascade2LoopGravity", 
            results_dir=str(results_dir),
            enable_gravity_portal=False  # 2-loop version without gravity
        )
        
        print("\n🔄 Solving 2-Loop RGE system...")
        solver.solve()
        
        print("\n📊 Generating plots...")
        solver.make_plots()
        
        # Load results for analysis
        import pandas as pd
        df = pd.read_csv(results_dir / "gauge_couplings.csv")
        
        # Additional analyses
        print("\n🔍 Running stability analysis...")
        plot_stability_analysis(df, output_path=results_dir / "stability_analysis.png")
        
        print("\n⚖️  Checking perturbativity...")
        check_perturbativity(df)
        
        print("\n" + "=" * 60)
        print("✅ E8 Cascade 2-Loop Analysis Complete!")
        print("=" * 60)
        print(f"📁 Results saved in: {results_dir}")
        print("\nGenerated files:")
        print("  • gauge_couplings.csv - Numerical RGE solutions")
        print("  • gauge_running.png - Coupling evolution plot")
        print("  • unification.png - Gauge unification plot") 
        print("  • stability_analysis.png - Detailed analysis plots")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
