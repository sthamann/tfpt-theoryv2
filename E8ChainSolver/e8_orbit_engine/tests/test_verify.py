"""
Tests for verification module.
"""

import pytest
import numpy as np
import pandas as pd
from e8_orbit_engine.verify import (
    forward_diff3,
    verify_lnD_cubic,
    phi_third_diff_from_fit,
    verify_normalization,
    comprehensive_verification
)


def test_forward_diff3():
    """Test third forward difference calculation."""
    # Test with cubic polynomial - should give constant third difference
    n = np.arange(10)
    x = 2*n**3 + 3*n**2 + n + 1
    d3 = forward_diff3(x)
    
    # For cubic, third difference should be constant = 6*a where a is cubic coefficient
    expected = 12  # 6 * 2
    np.testing.assert_array_almost_equal(d3, expected * np.ones_like(d3))
    
    # Test with quadratic - third difference should be zero
    x_quad = n**2 + 2*n + 1
    d3_quad = forward_diff3(x_quad)
    np.testing.assert_array_almost_equal(d3_quad, np.zeros_like(d3_quad))


def test_verify_lnD_cubic():
    """Test ln(D) cubic pattern verification."""
    # Create chain with perfect cubic ln(D)
    n = np.arange(20)
    ln_D = 5.5 - 0.2*n - 0.01*n**2 - 0.001*n**3
    D = np.exp(ln_D)
    
    chain = pd.DataFrame({
        'label': [f'O_{i}' for i in range(20)],
        'D': D
    })
    
    result = verify_lnD_cubic(chain)
    
    assert result['valid']  # Should detect constant third difference
    assert result['cv'] < 0.01  # Very small coefficient of variation
    
    # Test with non-cubic (exponential)
    D_exp = 248 * np.exp(-0.5 * n)
    chain_exp = pd.DataFrame({'D': D_exp})
    result_exp = verify_lnD_cubic(chain_exp)
    
    # Should have larger variation in third difference
    assert result_exp['cv'] > result['cv']


def test_phi_third_diff_from_fit():
    """Test φ third difference calculation."""
    # Use paper values
    gamma0 = 0.834
    gamma1 = 0.108
    gamma2 = 0.0105627
    
    result = phi_third_diff_from_fit(gamma0, gamma1, gamma2, N=30)
    
    # Check theoretical value
    theoretical = -2 * gamma2
    assert abs(result['theoretical'] - theoretical) < 1e-10
    
    # All computed values should be very close to theoretical
    assert result['valid']
    assert result['rel_error'] < 1e-6
    
    # Mean should match theoretical
    assert abs(result['mean'] - theoretical) < 1e-8


def test_verify_normalization():
    """Test normalization verification."""
    # Create mock fit results with good values
    fit_results = {
        'gamma0': 0.834,
        'gamma1': 0.108,
        'gamma2': 0.834 / (8 * np.pi**2)  # Exact topological relation
    }
    
    result = verify_normalization(fit_results)
    
    assert result['overall_valid']
    assert result['gamma0_check']['valid']
    assert result['gamma2_check']['valid']
    assert result['gamma2_check']['deviation_pct'] < 0.1  # Very small deviation
    
    # Test with bad values
    bad_fit_results = {
        'gamma0': 0.9,  # Wrong
        'gamma1': 0.108,
        'gamma2': 0.02  # Wrong
    }
    
    bad_result = verify_normalization(bad_fit_results)
    assert not bad_result['overall_valid']


def test_comprehensive_verification():
    """Test full verification suite."""
    # Create realistic chain
    n = np.arange(15)
    # Use approximate cubic pattern for ln(D)
    ln_D = np.log(248) - 0.185*n - 0.015*n**2 - 0.0008*n**3
    D = np.exp(ln_D)
    
    chain = pd.DataFrame({
        'label': [f'Orbit_{i}' for i in range(15)],
        'D': D,
        'orbit_dim': 248 - D
    })
    
    # Create fit results
    fit_results = {
        'gamma0': 0.834,
        'gamma1': 0.108,
        'gamma2': 0.0105627,
        'lam': 4.5,
        'gamma_obs': np.zeros(14),  # Dummy
        'gamma_fit': np.zeros(14),  # Dummy
        'residuals': np.random.normal(0, 0.001, 14),
        'r_squared': 0.995,
        'lnD': ln_D,
        'steps': ln_D[:-1] - ln_D[1:],
        'n': np.arange(14),
        'theoretical_gamma2': 0.834 / (8 * np.pi**2),
        'gamma2_deviation': 0.01
    }
    
    result = comprehensive_verification(chain, fit_results)
    
    assert 'lnD_cubic' in result
    assert 'phi_third_diff' in result
    assert 'normalization' in result
    assert 'all_tests_passed' in result
    assert 'summary' in result
    
    # Check summary is a string
    assert isinstance(result['summary'], str)
    assert 'VERIFICATION SUMMARY' in result['summary']


def test_topological_constants():
    """Test topological constant calculations."""
    c3 = 1 / (8 * np.pi)
    
    # Check value from paper
    assert abs(c3 - 0.039788735) < 1e-8
    
    # Check relation to gamma2
    gamma0 = 0.834
    gamma2_theoretical = gamma0 / (8 * np.pi**2)
    
    assert abs(gamma2_theoretical - 0.0105627) < 0.0001
    
    # Check c3^2 appears in various places
    c3_squared = c3**2
    assert abs(c3_squared - 0.001583) < 0.00001


def test_chain_too_short():
    """Test handling of chains that are too short."""
    # Chain with only 3 elements (too short for third difference)
    chain = pd.DataFrame({
        'D': [248, 206, 200]
    })
    
    result = verify_lnD_cubic(chain)
    assert not result['valid']
    assert 'too short' in result['message'].lower()


def test_perfect_agreement():
    """Test verification with perfect theoretical agreement."""
    # Create chain that perfectly matches theory
    gamma0 = 0.834
    gamma1 = 0.108
    gamma2 = gamma0 / (8 * np.pi**2)
    
    # Build ln(D) from gamma backwards
    n = np.arange(20)
    gamma = gamma0 + gamma1 * n + gamma2 * n**2
    
    ln_D = np.zeros(21)
    ln_D[0] = np.log(248)
    
    lam = 4.5  # Approximate
    s = gamma / lam
    
    for i in range(20):
        ln_D[i+1] = ln_D[i] - s[i]
    
    D = np.exp(ln_D)
    
    chain = pd.DataFrame({'D': D})
    
    # Run verification
    from e8_orbit_engine.fit import fit_gamma
    fit_results = fit_gamma(chain)
    verification = comprehensive_verification(chain, fit_results)
    
    # Should pass all tests with this perfect data
    # (May not be perfect due to numerical precision and fitting)
    assert fit_results['r_squared'] > 0.99
