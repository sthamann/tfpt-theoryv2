"""Analysis tools for E8 Cascade model"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def load_results(results_dir):
    """Load saved results from CSV"""
    results_dir = Path(results_dir)
    csv_path = results_dir / "gauge_couplings.csv"
    
    if not csv_path.exists():
        raise FileNotFoundError(f"No results found at {csv_path}")
    
    return pd.read_csv(csv_path)


def find_unification_scale(df, threshold=1.0):
    """Find the scale where g1 and g2 are closest"""
    # Use GUT normalization for g1
    diff = np.abs(df['alpha1_inv_GUT'] - df['alpha2_inv'])
    idx_min = diff.idxmin()
    
    return {
        'mu_GeV': df.loc[idx_min, 'mu_GeV'],
        'log10_mu': df.loc[idx_min, 'log10_mu'],
        'alpha1_inv_GUT': df.loc[idx_min, 'alpha1_inv_GUT'],
        'alpha2_inv': df.loc[idx_min, 'alpha2_inv'],
        'difference': diff.iloc[idx_min],
        'unified': diff.iloc[idx_min] < threshold
    }


def plot_stability_analysis(df, output_path=None):
    """Analyze the stability of running"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    mu = df['mu_GeV']
    
    # Plot 1: Gauge coupling evolution
    ax = axes[0, 0]
    ax.plot(mu, df['g1_GUT'], 'b-', label='g₁ (GUT)')
    ax.plot(mu, df['g2'], 'r-', label='g₂')
    ax.plot(mu, df['g3'], 'g-', label='g₃')
    ax.set_xscale('log')
    ax.set_xlabel('μ [GeV]')
    ax.set_ylabel('Gauge couplings')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Gauge Coupling Evolution')
    
    # Plot 2: Beta functions (numerical derivative)
    ax = axes[0, 1]
    dt = np.diff(df['log10_mu'].values)
    # Convert from d/d(log10 mu) to d/d(ln mu) by multiplying by ln(10)
    beta_g1 = np.diff(df['g1_GUT'].values) / dt * np.log(10)
    beta_g2 = np.diff(df['g2'].values) / dt * np.log(10)
    beta_g3 = np.diff(df['g3'].values) / dt * np.log(10)
    
    ax.plot(mu[1:], beta_g1, 'b-', label='β(g₁)')
    ax.plot(mu[1:], beta_g2, 'r-', label='β(g₂)')
    ax.plot(mu[1:], beta_g3, 'g-', label='β(g₃)')
    ax.set_xscale('log')
    ax.set_xlabel('μ [GeV]')
    ax.set_ylabel('Beta functions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Numerical Beta Functions')
    
    # Plot 3: Unification measure
    ax = axes[1, 0]
    diff_12 = np.abs(df['alpha1_inv_GUT'] - df['alpha2_inv'])
    diff_13 = np.abs(df['alpha1_inv_GUT'] - df['alpha3_inv'])
    diff_23 = np.abs(df['alpha2_inv'] - df['alpha3_inv'])
    
    ax.plot(mu, diff_12, 'b-', label='|α₁⁻¹ - α₂⁻¹|')
    ax.plot(mu, diff_13, 'r-', label='|α₁⁻¹ - α₃⁻¹|')
    ax.plot(mu, diff_23, 'g-', label='|α₂⁻¹ - α₃⁻¹|')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('μ [GeV]')
    ax.set_ylabel('Coupling differences')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Gauge Unification Measure')
    
    # Plot 4: Running of alpha_EM
    ax = axes[1, 1]
    # α_EM = (g_Y² cos²θ_W + g2² sin²θ_W)/(4π)
    # At tree level: sin²θ_W ≈ g_Y²/(g_Y² + g2²)
    # Use SM normalization for EM coupling!
    sin2_theta = df['g1_SM']**2 / (df['g1_SM']**2 + df['g2']**2)
    cos2_theta = 1 - sin2_theta
    alpha_em = (df['g1_SM']**2 * cos2_theta + df['g2']**2 * sin2_theta) / (4*np.pi)
    
    ax.plot(mu, 1/alpha_em, 'k-', linewidth=2)
    ax.axhline(y=137.035999, color='r', linestyle='--', label='α⁻¹(M_Z) exp')
    ax.set_xscale('log')
    ax.set_xlabel('μ [GeV]')
    ax.set_ylabel('α_EM⁻¹')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Electromagnetic Coupling')
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"✓ Stability analysis saved to {output_path}")
    else:
        plt.show()
    
    plt.close()


def check_perturbativity(df, threshold=4*np.pi):
    """Check if couplings remain perturbative"""
    violations = []
    
    # Check both SM and GUT normalized g1
    for col, name in [('g1_SM', 'g1(SM)'), ('g1_GUT', 'g1(GUT)'), ('g2', 'g2'), ('g3', 'g3')]:
        if col in df.columns:
            mask = df[col] > threshold
            if mask.any():
                first_violation = df.loc[mask].iloc[0]
                violations.append({
                    'coupling': name,
                    'mu_GeV': first_violation['mu_GeV'],
                    'value': first_violation[col]
                })
    
    if violations:
        print("⚠️  Perturbativity violations detected:")
        for v in violations:
            print(f"   {v['coupling']} > {threshold:.2f} at μ = {v['mu_GeV']:.2e} GeV")
    else:
        print("✓ All couplings remain perturbative")
    
    return violations


def analyze_gravity_effect(results_dir_vanilla, results_dir_gravity):
    """Compare runs with and without gravity portal"""
    df_vanilla = load_results(results_dir_vanilla)
    df_gravity = load_results(results_dir_gravity)
    
    # Find unification scales
    unif_vanilla = find_unification_scale(df_vanilla)
    unif_gravity = find_unification_scale(df_gravity)
    
    print("Gravity Portal Effect Analysis:")
    print("=" * 50)
    print(f"Vanilla: Closest approach at μ = {unif_vanilla['mu_GeV']:.2e} GeV")
    print(f"         |Δ| = {unif_vanilla['difference']:.2f}")
    print(f"Gravity: Closest approach at μ = {unif_gravity['mu_GeV']:.2e} GeV")
    print(f"         |Δ| = {unif_gravity['difference']:.2f}")
    print(f"\nImprovement: {(1 - unif_gravity['difference']/unif_vanilla['difference'])*100:.1f}%")
    
    # Plot comparison
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.plot(df_vanilla['mu_GeV'], df_vanilla['alpha1_inv_GUT'], 'b--', 
            alpha=0.5, label='α₁⁻¹ (GUT, vanilla)')
    ax.plot(df_gravity['mu_GeV'], df_gravity['alpha1_inv_GUT'], 'b-', 
            label='α₁⁻¹ (GUT, gravity)')
    
    ax.plot(df_vanilla['mu_GeV'], df_vanilla['alpha2_inv'], 'r--', 
            alpha=0.5, label='α₂⁻¹ (vanilla)')
    ax.plot(df_gravity['mu_GeV'], df_gravity['alpha2_inv'], 'r-', 
            label='α₂⁻¹ (gravity)')
    
    ax.set_xscale('log')
    ax.set_xlabel('μ [GeV]', fontsize=14)
    ax.set_ylabel('α⁻¹', fontsize=14)
    ax.set_title('Effect of Gravity Portal on Gauge Unification', fontsize=16)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show() 