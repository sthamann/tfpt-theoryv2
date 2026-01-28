"""
Tests for fitting module.
"""

import pytest
import numpy as np
import pandas as pd
from e8_orbit_engine.fit import (
    fit_gamma, 
    quad, 
    gamma_scale_from_first_step,
    analyze_fit_quality
)


def test_quad_function():
    """Test quadratic function."""
    n = np.array([0, 1, 2, 3])
    result = quad(n, 1.0, 0.5, 0.1)
    expected = np.array([1.0, 1.6, 2.4, 3.4])
    np.testing.assert_array_almost_equal(result, expected)


def test_gamma_scale():
    """Test scaling calculation."""
    s0 = 0.1852  # Approximate ln(248/206)
    lam = gamma_scale_from_first_step(s0, gamma0_target=0.834)
    assert abs(lam - 4.5) < 0.1  # Should be around 4.5


def test_fit_synthetic_data():
    """Test fitting on synthetic quadratic data."""
    # Create perfect quadratic data
    n = np.arange(20)
    true_a, true_b, true_c = 0.834, 0.108, 0.0105627
    
    # Create synthetic D values that give the right pattern
    gamma_true = quad(n, true_a, true_b, true_c)
    
    # Work backwards: if gamma = lambda * s, and s = ln(D_n) - ln(D_{n+1})
    # Then ln(D_n) is a cubic with appropriate coefficients
    ln_D = np.zeros(21)
    ln_D[0] = np.log(248)  # Start at adjoint
    
    # Use approximate lambda
    lam = 4.5
    s = gamma_true / lam
    
    for i in range(20):
        ln_D[i+1] = ln_D[i] - s[i]
    
    D = np.exp(ln_D)
    
    # Create chain DataFrame
    chain = pd.DataFrame({
        'label': [f'Orbit_{i}' for i in range(21)],
        'D': D,
        'orbit_dim': 248 - D
    })
    
    # Fit
    results = fit_gamma(chain, gamma0_target=0.834)
    
    # Check results
    assert abs(results['gamma0'] - true_a) < 0.01
    assert abs(results['gamma1'] - true_b) < 0.01
    assert abs(results['gamma2'] - true_c) < 0.001
    assert results['r_squared'] > 0.99


def test_fit_real_orbits():
    """Test fitting on actual E8 orbit data."""
    # Create realistic chain based on known orbits
    # Using key orbits from the paper (starts at A4+A1, not adjoint)
    orbits_data = [
        ('A4+A1', 206),  # Key step (adjoint 248 will be added by include_adjoint)
        ('D5', 200),
        ('A5', 196),
        ('A4+2A1', 192),
        ('2A3', 188),
        ('D4+A1', 184),
        ('A4', 180),
        ('D4(a1)+A1', 176),
        ('A3+2A1', 172),
        ('D4', 168),
        ('D4(a1)', 166),
        ('A3+A1', 164),
        ('2A2+A1', 162),
        ('2A2', 156),
        ('A2+3A1', 154),
        ('A3', 148),
        ('A2+2A1', 146),
        ('A2+A1', 136),
        ('4A1', 128)
    ]
    
    chain = pd.DataFrame(orbits_data, columns=['label', 'D'])
    chain['orbit_dim'] = 248 - chain['D']
    
    # Fit with the corrected anchoring
    results = fit_gamma(chain, gamma0_target=0.834, include_adjoint=True)
    
    # With the fix, values should be much closer to theory
    assert abs(results['gamma0'] - 0.834) < 0.01  # Within 1%
    assert abs(results['gamma1'] - 0.108) < 0.01  # Within ~10%
    assert abs(results['gamma2'] - 0.0105627) < 0.001  # Within ~10%
    
    # Check topological constraint
    theoretical_gamma2 = results['gamma0'] / (8 * np.pi**2)
    assert abs(results['gamma2'] - theoretical_gamma2) / theoretical_gamma2 < 0.1


def test_analyze_fit_quality():
    """Test fit quality analysis."""
    # Create mock fit results
    fit_results = {
        'gamma0': 0.835,
        'gamma1': 0.107,
        'gamma2': 0.0106,
        'r_squared': 0.995,
        'residuals': np.random.normal(0, 0.001, 20)
    }
    
    quality = analyze_fit_quality(fit_results)
    
    assert 'deviations' in quality
    assert 'topological_check' in quality
    assert 'quality' in quality
    
    # Should be excellent or good quality
    assert quality['quality'] in ['excellent', 'good']


def test_weighted_fit():
    """Test weighted fitting option."""
    # Create chain with some noise
    np.random.seed(42)
    n = np.arange(15)
    D = 248 * np.exp(-0.2 * n) + np.random.normal(0, 2, 15)
    D = np.maximum(D, 10)  # Ensure positive
    
    chain = pd.DataFrame({
        'label': [f'O_{i}' for i in range(15)],
        'D': D,
        'orbit_dim': 248 - D
    })
    
    # Fit with and without weighting
    results_unweighted = fit_gamma(chain, weighted=False)
    results_weighted = fit_gamma(chain, weighted=True)
    
    # Both should give reasonable fits
    assert results_unweighted['r_squared'] > 0.8
    assert results_weighted['r_squared'] > 0.8
    
    # Weighted should emphasize early terms more
    # (This is a soft test - exact behavior depends on data)
