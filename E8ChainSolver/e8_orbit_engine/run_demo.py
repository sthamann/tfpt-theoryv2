#!/usr/bin/env python3
"""
Standalone demo script for E8 Orbit Engine.
Runs without full installation.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

# Basic test with mock data if pandas not available
try:
    import pandas as pd
    import numpy as np
    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False
    print("WARNING: pandas/numpy not installed. Using mock data.")

def run_demo():
    """Run a basic demonstration of the engine."""
    
    if not HAS_DEPS:
        print("\n" + "="*60)
        print("E8 ORBIT ENGINE - Mock Demo (Install dependencies for full version)")
        print("="*60)
        
        # Mock data for demonstration
        print("\n1. Mock orbit chain (first 10 orbits):")
        orbits = [
            ("0 (adjoint)", 248, 0),
            ("A4+A1", 206, 42),
            ("D5", 200, 48),
            ("A5", 196, 52),
            ("A4+2A1", 192, 56),
            ("2A3", 188, 60),
            ("D4+A1", 184, 64),
            ("A4", 180, 68),
            ("D4(a1)+A1", 176, 72),
            ("A3+2A1", 172, 76)
        ]
        
        for i, (label, D, orbit_dim) in enumerate(orbits):
            print(f"  n={i:2d}: {label:15s} D={D:3d}  orbit_dim={orbit_dim:3d}")
        
        print("\n2. Expected quadratic fit:")
        print("  γ(n) ≈ 0.834 + 0.108n + 0.0105627n²")
        
        print("\n3. Key theoretical relations:")
        print(f"  c₃ = 1/(8π) = {1/(8*3.14159):.8f}")
        print(f"  γ₂ = γ₀/(8π²) = 0.834/(8π²) = {0.834/(8*3.14159**2):.8f}")
        print(f"  Expected: γ₂ ≈ 0.0105627")
        
        print("\n4. Calibration-free test:")
        print("  Δ³ln(φ) = -2γ₂ = -0.0211254 (constant)")
        
        print("\n" + "="*60)
        print("Install dependencies with: pip3 install pandas numpy scipy matplotlib")
        print("Then run: python3 -m src.e8_orbit_engine.cli fit data/nilpotent_orbits.csv")
        return
    
    # Full demo with actual engine
    from e8_orbit_engine import load_orbits, build_chain, fit_gamma
    from e8_orbit_engine.verify import comprehensive_verification
    from e8_orbit_engine.fit import analyze_fit_quality
    
    print("\n" + "="*60)
    print("E8 ORBIT ENGINE - Full Demo")
    print("="*60)
    
    # Load data
    data_path = Path(__file__).parent / 'data' / 'nilpotent_orbits.csv'
    print(f"\n1. Loading orbit data from {data_path.name}...")
    df = load_orbits(data_path)
    print(f"   ✓ Loaded {len(df)} orbits")
    print(f"   D range: {df['D'].min():.0f} - {df['D'].max():.0f}")
    
    # Build chain
    print("\n2. Building monotonic chain...")
    chain = build_chain(df, max_len=20)
    print(f"   ✓ Built chain with {len(chain)} orbits")
    
    # Fit gamma
    print("\n3. Fitting quadratic γ(n)...")
    results = fit_gamma(chain, gamma0_target=0.834)
    
    print(f"\n   Fitted parameters:")
    print(f"   γ₀ = {results['gamma0']:.6f} (target: 0.834)")
    print(f"   γ₁ = {results['gamma1']:.6f} (target: 0.108)")
    print(f"   γ₂ = {results['gamma2']:.8f} (target: 0.0105627)")
    print(f"\n   R² = {results['r_squared']:.6f}")
    
    # Verify
    print("\n4. Running verification tests...")
    verification = comprehensive_verification(chain, results)
    
    if verification['all_tests_passed']:
        print("   ✓ All verification tests passed!")
    else:
        print("   ⚠ Some tests need attention")
    
    # Quality analysis
    quality = analyze_fit_quality(results)
    print(f"\n5. Quality assessment: {quality['quality'].upper()}")
    
    # Show deviations
    print("\n   Parameter deviations from theoretical values:")
    for param in ['gamma0', 'gamma1', 'gamma2']:
        dev = quality['deviations'][param]
        print(f"   {param}: {dev['deviation_pct']:.2f}%")
    
    print("\n" + "="*60)
    print("✓ Demo complete!")
    print("\nFor full analysis with plots and HTML report, run:")
    print("  python3 -m src.e8_orbit_engine.cli fit data/nilpotent_orbits.csv")
    print("="*60)

if __name__ == "__main__":
    run_demo()
