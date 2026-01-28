"""
Critical tests for the normalization fix.
These tests verify that the first step anchoring is correct.
"""

import pytest
import numpy as np
import pandas as pd
from e8_orbit_engine.fit import fit_gamma
from e8_orbit_engine.verify import phi_third_diff_from_fit


def test_first_step_anchoring():
    """Test that s₀ is correctly anchored at ln(248/206)."""
    # Create a chain starting at 206 (A4+A1)
    chain = pd.DataFrame({
        'label': ['A4+A1', 'D5', 'A5'],
        'D': [206.0, 190.0, 175.0],
        'orbit_dim': [42, 58, 73]
    })
    
    results = fit_gamma(chain, include_adjoint=True)
    
    # Check s₀ against theoretical value
    s0_expected = np.log(248.0 / 206.0)  # 0.1855525773754012
    s0_actual = results['s0_true']
    
    assert abs(s0_actual - s0_expected) < 1e-10, \
        f"s₀ = {s0_actual:.10f}, expected {s0_expected:.10f}"
    
    # Check that λ is correctly calculated
    lam_expected = 0.834 / s0_expected
    lam_actual = results['lam']
    assert abs(lam_actual - lam_expected) < 1e-6


def test_topological_constraint():
    """Test that γ₂ = γ₀/(8π²) to better than 1%."""
    # Create theoretical chain
    gamma0 = 0.834
    gamma1 = 0.108
    gamma2 = 0.0105627
    
    n = np.arange(15)
    gamma = gamma0 + gamma1 * n + gamma2 * n**2
    
    # Build D values from gamma
    D0 = 248
    ln_D = np.zeros(15)
    ln_D[0] = np.log(206)  # Start at A4+A1
    
    s0_true = np.log(248/206)
    lam = gamma0 / s0_true
    s = gamma / lam
    
    for i in range(14):
        ln_D[i+1] = ln_D[i] - s[i+1]  # Use s[i+1] because s[0] is 248->206
    
    D = np.exp(ln_D)
    
    chain = pd.DataFrame({
        'label': [f'Orbit_{i}' for i in range(15)],
        'D': D,
        'orbit_dim': 248 - D
    })
    
    results = fit_gamma(chain, include_adjoint=True)
    
    # Check topological constraint
    gamma2_fitted = results['gamma2']
    gamma2_theoretical = results['gamma0'] / (8 * np.pi**2)
    
    deviation = abs(gamma2_fitted - gamma2_theoretical) / gamma2_theoretical
    
    assert deviation < 0.01, \
        f"γ₂ deviation {deviation*100:.2f}% exceeds 1% tolerance"
    
    # Check exact value
    assert abs(gamma2_fitted - 0.0105627) < 0.0001, \
        f"γ₂ = {gamma2_fitted:.8f}, expected 0.0105627"


def test_calibration_free_third_difference():
    """Test that Δ³ln(φ) = -0.0211254668 to 1e-8."""
    # Use exact theoretical parameters
    gamma0 = 0.834
    gamma1 = 0.108
    gamma2 = 0.0105627
    
    result = phi_third_diff_from_fit(gamma0, gamma1, gamma2, N=30)
    
    # Check against paper value
    target = -0.0211254668
    d3_phi = result['d3_phi']
    
    # All values should be very close to target
    for i, val in enumerate(d3_phi):
        assert abs(val - target) < 1e-8, \
            f"Δ³ln(φ)[{i}] = {val:.10f}, expected {target:.10f}"
    
    # Mean should be exact
    mean_d3 = np.mean(d3_phi) if len(d3_phi) > 0 else 0
    assert abs(mean_d3 - target) < 1e-8, \
        f"Mean Δ³ln(φ) = {mean_d3:.10f}, expected {target:.10f}"


def test_complete_workflow_with_fix():
    """Test the complete workflow with the normalization fix."""
    # Create a realistic E8 chain
    # Using known orbits with correct centralizer dimensions
    orbits = [
        ('A4+A1', 206),
        ('D5', 200),
        ('A5', 196),
        ('A4+2A1', 192),
        ('2A3', 188),
        ('D4+A1', 184),
        ('A4', 180),
        ('D4(a1)+A1', 176),
        ('A3+2A1', 172),
        ('D4', 168)
    ]
    
    chain = pd.DataFrame(orbits, columns=['label', 'D'])
    chain['orbit_dim'] = 248 - chain['D']
    
    # Fit with fix
    results = fit_gamma(chain, include_adjoint=True)
    
    # Verify all three critical assertions
    
    # 1. First step check
    s0_expected = np.log(248.0 / 206.0)
    assert abs(results['s0_true'] - s0_expected) < 1e-10
    
    # 2. Topological constraint
    gamma2_theoretical = results['gamma0'] / (8 * np.pi**2)
    deviation = abs(results['gamma2'] - gamma2_theoretical) / gamma2_theoretical
    assert deviation < 0.02  # 2% tolerance for real data
    
    # 3. Third difference (using fitted parameters)
    d3_result = phi_third_diff_from_fit(
        results['gamma0'], 
        results['gamma1'], 
        results['gamma2']
    )
    
    # Should be close to -2*gamma2
    expected_d3 = -2 * results['gamma2']
    mean_d3 = np.mean(d3_result['d3_phi']) if len(d3_result['d3_phi']) > 0 else 0
    assert abs(mean_d3 - expected_d3) < 1e-6


def test_backward_compatibility():
    """Test that include_adjoint=False preserves old behavior."""
    chain = pd.DataFrame({
        'label': ['O1', 'O2', 'O3'],
        'D': [200.0, 180.0, 160.0],
        'orbit_dim': [48, 68, 88]
    })
    
    # Old behavior
    results_old = fit_gamma(chain, include_adjoint=False)
    
    # Should use first step as 200->180
    s0_old = np.log(200.0 / 180.0)
    assert abs(results_old['steps'][0] - s0_old) < 1e-10
