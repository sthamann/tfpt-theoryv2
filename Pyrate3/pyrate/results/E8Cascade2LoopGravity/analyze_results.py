#!/usr/bin/env python3
"""Quick analysis of E8 Cascade RGE results"""

import pandas as pd
import numpy as np
import sys

def analyze_results(csv_file='results/gauge_couplings.csv'):
    """Analyze RGE running results"""
    
    # Load data
    df = pd.read_csv(csv_file)
    
    print("=" * 60)
    print("E8 CASCADE RGE ANALYSIS")
    print("=" * 60)
    
    # 1. Initial conditions at M_Z
    print("\n1. Initial conditions at M_Z:")
    mz_row = df.iloc[0]
    print(f"   g₁(M_Z) = {mz_row['g1_SM']:.4f} (SM), {mz_row['g1_GUT']:.4f} (GUT)")
    print(f"   g₂(M_Z) = {mz_row['g2']:.4f}")
    print(f"   g₃(M_Z) = {mz_row['g3']:.4f}")
    print(f"   α₁⁻¹ = {mz_row['alpha1_inv_GUT']:.2f}")
    print(f"   α₂⁻¹ = {mz_row['alpha2_inv']:.2f}")
    print(f"   α₃⁻¹ = {mz_row['alpha3_inv']:.2f}")
    
    # 2. Unification analysis
    print("\n2. Gauge unification:")
    diff = (df['alpha1_inv_GUT'] - df['alpha2_inv']).abs()
    uni_idx = diff.idxmin()
    uni_row = df.iloc[uni_idx]
    
    print(f"   Closest approach at μ = {uni_row['mu_GeV']:.2e} GeV")
    print(f"   α₁⁻¹(GUT) = {uni_row['alpha1_inv_GUT']:.2f}")
    print(f"   α₂⁻¹ = {uni_row['alpha2_inv']:.2f}")
    print(f"   |Δ| = {diff.iloc[uni_idx]:.3f}")
    
    # Check if truly unified (within 1%)
    relative_diff = diff.iloc[uni_idx] / uni_row['alpha1_inv_GUT'] * 100
    if relative_diff < 1:
        print(f"   ✓ UNIFIED within {relative_diff:.1f}%!")
    else:
        print(f"   ✗ Not unified ({relative_diff:.1f}% difference)")
    
    # 3. Perturbativity check
    print("\n3. Perturbativity:")
    g_max = df[['g1_GUT', 'g2', 'g3']].max().max()
    alpha_max = g_max**2 / (4 * np.pi)
    
    if g_max < np.sqrt(4 * np.pi):
        print(f"   ✓ All couplings perturbative (g_max = {g_max:.3f} < {np.sqrt(4*np.pi):.3f})")
    else:
        print(f"   ✗ Non-perturbative! (g_max = {g_max:.3f})")
    print(f"   Max α = {alpha_max:.3f}")
    
    # 4. Running analysis
    print("\n4. Running characteristics:")
    
    # Where does g2 beta function change sign?
    if 'g2' in df.columns and len(df) > 1:
        g2_diff = np.diff(df['g2'].values)
        sign_changes = np.where(np.diff(np.sign(g2_diff)))[0]
        if len(sign_changes) > 0:
            for idx in sign_changes:
                print(f"   g₂ beta changes sign at μ ≈ {df.iloc[idx]['mu_GeV']:.2e} GeV")
    
    # 5. Threshold scales (if visible in data)
    print("\n5. Visible thresholds:")
    
    # Look for kinks in beta functions (2nd derivative)
    for coupling in ['g1_GUT', 'g2', 'g3']:
        if len(df) > 10:
            second_diff = np.diff(np.diff(df[coupling].values))
            threshold_candidates = np.where(np.abs(second_diff) > 0.001)[0]
            if len(threshold_candidates) > 0:
                print(f"   Possible thresholds in {coupling}:")
                for idx in threshold_candidates[:3]:  # Show max 3
                    print(f"     - μ ≈ {df.iloc[idx]['mu_GeV']:.2e} GeV")
    
    # 6. Gravity portal effect (if multiple files to compare)
    print("\n6. Gravity portal analysis:")
    print("   Run with different --cR-factor values to see the effect")
    print("   Example: python run_e8cascade.py --gravity --cR-factor 100")
    
    print("\n" + "=" * 60)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        analyze_results(sys.argv[1])
    else:
        analyze_results() 