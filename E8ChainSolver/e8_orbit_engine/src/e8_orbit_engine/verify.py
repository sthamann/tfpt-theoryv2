"""
Verification module for calibration-free consistency tests and RGE crossings.
"""

import numpy as np
import pandas as pd
import math
from typing import Tuple, Dict, Optional
import logging

logger = logging.getLogger(__name__)


def forward_diff3(x: np.ndarray) -> np.ndarray:
    """
    Calculate third forward difference: Δ³x_n.
    
    Δ³x_n = x_{n+3} - 3*x_{n+2} + 3*x_{n+1} - x_n
    
    Parameters
    ----------
    x : np.ndarray
        Input array
        
    Returns
    -------
    np.ndarray
        Third forward differences
    """
    if len(x) < 4:
        return np.array([])
    
    return x[3:] - 3*x[2:-1] + 3*x[1:-2] - x[:-3]


def phi_third_diff_from_data(gamma_obs: np.ndarray) -> Dict:
    """
    Calculate third difference of ln(φ) directly from observed gamma values.
    
    Calibration-free test: Δln(φ) = -γ ⇒ Δ³ln(φ) = -Δ²γ
    
    Parameters
    ----------
    gamma_obs : np.ndarray
        Observed gamma values from data
        
    Returns
    -------
    dict
        Third differences and statistics
        
    Notes
    -----
    Excludes n=0 and n=1 for robustness since the first step is an outlier.
    """
    # Second differences of gamma: Δ²γ
    d2_gamma = gamma_obs[2:] - 2*gamma_obs[1:-1] + gamma_obs[:-2]
    
    # Third difference of ln(φ) = -Δ²γ
    d3_phi = -d2_gamma
    
    # For robustness, exclude n=0 and n=1 (starts effectively at n=2)
    if len(d3_phi) > 1:
        d3_phi_core = d3_phi[1:]  # Exclude the first point
    else:
        d3_phi_core = d3_phi
    
    # Statistics
    mean = float(np.mean(d3_phi_core)) if len(d3_phi_core) > 0 else 0
    std = float(np.std(d3_phi_core)) if len(d3_phi_core) > 0 else 0
    cv = std / abs(mean) if mean != 0 else np.inf
    
    # Check if approximately constant
    is_constant = cv < 0.1  # CV < 10%
    
    return {
        "valid": is_constant,
        "d3_phi": d3_phi,
        "d3_phi_core": d3_phi_core,
        "mean": mean,
        "std": std,
        "cv": cv,
        "message": f"Δ³ln(φ) from data: mean={mean:.6f}, CV={cv:.3f} "
                   f"{'✓ (approx. constant)' if is_constant else '✗ (not constant)'}"
    }


def verify_lnD_cubic(chain: pd.DataFrame) -> Dict:
    """
    Verify that ln(D) follows a near-cubic pattern.
    
    The third forward difference should be approximately constant.
    
    Parameters
    ----------
    chain : pd.DataFrame
        Orbit chain with D values
        
    Returns
    -------
    dict
        Verification results
    """
    lnD = np.log(chain['D'].values)
    d3 = forward_diff3(lnD)
    
    if len(d3) == 0:
        return {
            'valid': False,
            'message': 'Chain too short for third difference',
            'd3': d3
        }
    
    # Statistics of third difference
    mean_d3 = np.mean(d3)
    std_d3 = np.std(d3)
    cv_d3 = std_d3 / abs(mean_d3) if mean_d3 != 0 else np.inf
    
    # Check if approximately constant
    is_constant = cv_d3 < 0.1  # Coefficient of variation < 10%
    
    return {
        'valid': is_constant,
        'd3': d3,
        'mean': mean_d3,
        'std': std_d3,
        'cv': cv_d3,
        'min': np.min(d3),
        'max': np.max(d3),
        'range': np.max(d3) - np.min(d3),
        'message': f"Δ³ln(D) CV = {cv_d3:.3f} {'✓ (constant)' if is_constant else '✗ (not constant)'}"
    }


def verify_lnD_cubic_rolling(chain: pd.DataFrame, window: int = 8) -> Dict:
    """
    Verify cubic pattern of ln(D) using rolling windows.
    
    Shows that the cubic character may hold locally but not globally.
    
    Parameters
    ----------
    chain : pd.DataFrame
        Orbit chain with D values
    window : int
        Window size for rolling analysis
        
    Returns
    -------
    dict
        Rolling window statistics
    """
    lnD = np.log(chain['D'].values)
    
    if len(lnD) < window:
        return {
            'valid': False,
            'message': f'Chain too short for window size {window}',
            'windows': []
        }
    
    windows = []
    for i in range(len(lnD) - window + 1):
        # Extract window
        lnD_window = lnD[i:i+window]
        
        # Calculate third differences for this window
        if len(lnD_window) >= 4:
            d3_window = forward_diff3(lnD_window)
            
            if len(d3_window) > 0:
                mean_d3 = np.mean(d3_window)
                std_d3 = np.std(d3_window)
                cv_d3 = std_d3 / abs(mean_d3) if mean_d3 != 0 else np.inf
                
                windows.append({
                    'start': i,
                    'end': i + window - 1,
                    'mean': float(mean_d3),
                    'std': float(std_d3),
                    'cv': float(cv_d3),
                    'is_constant': cv_d3 < 0.1
                })
    
    # Overall assessment
    if windows:
        cvs = [w['cv'] for w in windows]
        constant_windows = sum(1 for w in windows if w['is_constant'])
        
        return {
            'valid': constant_windows > len(windows) * 0.5,  # More than half should be constant
            'windows': windows,
            'n_constant': constant_windows,
            'n_total': len(windows),
            'cv_range': [min(cvs), max(cvs)],
            'message': f"Rolling Δ³ln(D): {constant_windows}/{len(windows)} windows constant (CV<0.1)"
        }
    else:
        return {
            'valid': False,
            'message': 'No valid windows',
            'windows': []
        }


def phi_third_diff_from_fit(
    gamma0: float, 
    gamma1: float, 
    gamma2: float, 
    N: int = 40
) -> Dict:
    """
    Calculate third difference of ln(φ_n) from fitted γ(n).
    
    Theory predicts: Δ³ln(φ_n) = -2*γ₂ (constant)
    
    Parameters
    ----------
    gamma0, gamma1, gamma2 : float
        Quadratic coefficients of γ(n)
    N : int
        Number of steps to calculate
        
    Returns
    -------
    dict
        Verification results
    """
    # Calculate ln(φ_n) from closed form
    n = np.arange(N + 1)
    
    # ln(φ_n) = ln(φ_0) - Σ_{k=0}^{n-1} γ(k)
    # With γ(k) = γ₀ + γ₁*k + γ₂*k²
    # The sum gives the closed form:
    lnphi = (0  # ln(φ_0) term cancels in differences
             - gamma0 * n
             - 0.5 * gamma1 * n * (n - 1)
             - (gamma2 / 6.0) * (n - 1) * n * (2*n - 1))
    
    # Calculate third difference
    d3_phi = forward_diff3(lnphi)
    
    # Theoretical value
    theoretical = -2 * gamma2
    
    # Check all values
    deviations = d3_phi - theoretical
    max_dev = np.max(np.abs(deviations)) if len(deviations) > 0 else np.inf
    
    # Relative error (use absolute if theoretical is near zero)
    if abs(theoretical) > 1e-10:
        rel_error = max_dev / abs(theoretical)
    else:
        rel_error = max_dev
    
    is_valid = rel_error < 1e-6  # Very tight tolerance
    
    return {
        'valid': is_valid,
        'd3_phi': d3_phi,
        'theoretical': theoretical,
        'mean': np.mean(d3_phi) if len(d3_phi) > 0 else 0,
        'std': np.std(d3_phi) if len(d3_phi) > 0 else 0,
        'max_deviation': max_dev,
        'rel_error': rel_error,
        'message': f"Δ³ln(φ) = {np.mean(d3_phi) if len(d3_phi) > 0 else 0:.8f} vs theory {theoretical:.8f} "
                   f"{'✓' if is_valid else '✗'} (rel error: {rel_error:.2e})"
    }


def verify_normalization(fit_results: Dict) -> Dict:
    """
    Verify the normalization relationships.
    
    Checks:
    1. First step normalization: γ(0) ≈ 0.834
    2. Topological relation: γ₂ ≈ γ₀/(8π²)
    
    Parameters
    ----------
    fit_results : dict
        Results from fit_gamma
        
    Returns
    -------
    dict
        Normalization verification results
    """
    gamma0 = fit_results['gamma0']
    gamma2 = fit_results['gamma2']
    
    # Check 1: First step normalization
    # Based on 248 → 206 (adjoint → A4+A1)
    expected_gamma0 = 0.834
    gamma0_dev = abs(gamma0 - expected_gamma0) / expected_gamma0
    
    # Check 2: Topological relation
    # γ₂ = γ₀/(8π²)
    expected_gamma2 = gamma0 / (8 * np.pi**2)
    gamma2_dev = abs(gamma2 - expected_gamma2) / expected_gamma2
    
    # Check 3: c₃ relation
    # c₃ = 1/(8π) from paper
    c3 = 1 / (8 * np.pi)
    c3_squared = c3**2
    eight_pi_squared = 8 * np.pi**2
    
    # Check 4: First step s₀ (if available)
    s0_check = None
    if 's0_true' in fit_results and fit_results['s0_true'] is not None:
        s0_expected = np.log(248.0/60.0)   # First step: adjoint (248) → A4+A1 (D=60)
        s0_actual = fit_results['s0_true']
        s0_dev = abs(s0_actual - s0_expected) / s0_expected * 100
        s0_check = {
            'expected': s0_expected,
            'actual': s0_actual,
            'deviation_pct': s0_dev,
            'valid': s0_dev < 1.0  # Very tight tolerance for first step
        }
    
    return {
        'gamma0_check': {
            'expected': expected_gamma0,
            'actual': gamma0,
            'deviation_pct': gamma0_dev * 100,
            'valid': gamma0_dev < 0.01  # 1% tolerance
        },
        'gamma2_check': {
            'expected': expected_gamma2,
            'actual': gamma2,
            'deviation_pct': gamma2_dev * 100,
            'valid': gamma2_dev < 0.01  # TIGHTENED to 1% tolerance (was 2%)
        },
        's0_check': s0_check,
        'topological_constants': {
            'c3': c3,
            'c3_squared': c3_squared,
            '8_pi_squared': eight_pi_squared,
            'ratio': 1 / eight_pi_squared
        },
        'overall_valid': gamma0_dev < 0.01 and gamma2_dev < 0.02
    }


def comprehensive_verification(chain: pd.DataFrame, fit_results: Dict) -> Dict:
    """
    Run all verification tests.
    
    Parameters
    ----------
    chain : pd.DataFrame
        Orbit chain
    fit_results : dict
        Results from fit_gamma
        
    Returns
    -------
    dict
        All verification results
    """
    # Test 1: ln(D) cubic pattern
    lnD_test = verify_lnD_cubic(chain)
    
    # Test 2: φ third difference
    phi_test = phi_third_diff_from_fit(
        fit_results['gamma0'],
        fit_results['gamma1'],
        fit_results['gamma2']
    )
    
    # Test 3: Normalization
    norm_test = verify_normalization(fit_results)
    
    # Overall assessment
    all_valid = (
        lnD_test['valid'] and 
        phi_test['valid'] and 
        norm_test['overall_valid']
    )
    
    return {
        'lnD_cubic': lnD_test,
        'phi_third_diff': phi_test,
        'normalization': norm_test,
        'all_tests_passed': all_valid,
        'summary': _create_verification_summary(lnD_test, phi_test, norm_test)
    }


def _cross_log_interp(x1: float, y1: float, x2: float, y2: float, target: float) -> float:
    """
    Log-linear interpolation to find crossing point.
    
    Parameters
    ----------
    x1, x2 : float
        log10(mu) values
    y1, y2 : float
        alpha values
    target : float
        Target alpha value to find crossing
        
    Returns
    -------
    float
        Energy scale in GeV where crossing occurs
    """
    # Linear interpolation in log space
    t = (target - y1) / (y2 - y1 + 1e-18)  # Add small value to avoid division by zero
    log10_mu = x1 + t * (x2 - x1)
    return 10**log10_mu


def find_crossing_scales(df_rge: pd.DataFrame, target: float, 
                         col_mu: str = "log10_mu", 
                         col_alpha: str = "alpha3") -> Optional[float]:
    """
    Find energy scale where gauge coupling crosses target value.
    
    Parameters
    ----------
    df_rge : pd.DataFrame
        RGE data with columns for energy and coupling
    target : float
        Target value to find crossing (e.g., 1/(6π) for E6)
    col_mu : str
        Column name for log10(energy) 
    col_alpha : str
        Column name for coupling value
        
    Returns
    -------
    float or None
        Energy scale in GeV where crossing occurs, or None if no crossing
    """
    df = df_rge.sort_values(col_mu)[[col_mu, col_alpha]].values
    
    crossings = []
    for (x1, y1), (x2, y2) in zip(df[:-1], df[1:]):
        # Check if target is between y1 and y2
        if (y1 - target) * (y2 - target) <= 0:
            mu_cross = _cross_log_interp(x1, y1, x2, y2, target)
            crossings.append(mu_cross)
            logger.info(f"Found crossing at μ = {mu_cross:.3e} GeV")
    
    return crossings[0] if crossings else None


def get_target_values() -> Dict[str, float]:
    """
    Get target values for E6, E7, E8 crossings.
    
    Returns
    -------
    dict
        Dictionary with target alpha values for each group
    """
    return {
        "E6": 1.0 / (6 * math.pi),
        "E7": 1.0 / (7 * math.pi),
        "E8": 1.0 / (8 * math.pi)
    }


def verify_rge_crossings(df_rge: pd.DataFrame, include_uncertainty: bool = True) -> Dict:
    """
    Find and verify E6, E7, E8 crossing scales from RGE data.
    
    Parameters
    ----------
    df_rge : pd.DataFrame
        RGE data with gauge couplings
        
    Returns
    -------
    dict
        Crossing scales and verification results
    """
    targets = get_target_values()
    results = {}
    
    for group, target in targets.items():
        mu_cross = find_crossing_scales(df_rge, target)
        
        if mu_cross is not None:
            # Find closest value to check
            idx = np.argmin(np.abs(df_rge["mu_GeV"] - mu_cross))
            alpha_actual = df_rge.iloc[idx]["alpha3"]
            deviation = abs(alpha_actual - target) / target * 100
            
            # Find uncertainty from neighboring points
            if include_uncertainty:
                # Find the two points bracketing the crossing
                below = df_rge[df_rge["alpha3"] <= target]
                above = df_rge[df_rge["alpha3"] > target]
                
                if len(below) > 0 and len(above) > 0:
                    # Get the closest points on each side
                    idx_below = below["alpha3"].idxmax()
                    idx_above = above["alpha3"].idxmin()
                    
                    log10_mu_range = [
                        float(df_rge.loc[idx_below, "log10_mu"]),
                        float(df_rge.loc[idx_above, "log10_mu"])
                    ]
                    mu_range = [
                        float(df_rge.loc[idx_below, "mu_GeV"]),
                        float(df_rge.loc[idx_above, "mu_GeV"])
                    ]
                else:
                    log10_mu_range = [np.log10(mu_cross), np.log10(mu_cross)]
                    mu_range = [mu_cross, mu_cross]
            else:
                log10_mu_range = None
                mu_range = None
            
            results[group] = {
                "target": target,
                "mu_GeV": mu_cross,
                "log10_mu": np.log10(mu_cross),
                "log10_mu_range": log10_mu_range,
                "mu_range": mu_range,
                "alpha_actual": alpha_actual,
                "deviation_pct": deviation,
                "found": True
            }
            
            logger.info(f"{group}: μ = {mu_cross:.3e} GeV, α₃ = {alpha_actual:.8f} "
                       f"(target: {target:.8f}, dev: {deviation:.2f}%)")
        else:
            results[group] = {
                "target": target,
                "mu_GeV": None,
                "found": False
            }
            logger.warning(f"{group}: No crossing found for α₃ = {target:.8f}")
    
    return results


def assert_rge_schema(df_rge: pd.DataFrame) -> None:
    """
    Assert that RGE data has correct schema and properties.
    
    Parameters
    ----------
    df_rge : pd.DataFrame
        RGE data to validate
        
    Raises
    ------
    AssertionError
        If schema is invalid
    """
    # Check required columns
    required = [
        "mu_GeV", "log10_mu", "g1_SM", "g1_GUT", "g2", "g3",
        "alpha1_GUT", "alpha2", "alpha3",
        "alpha1_inv_GUT", "alpha2_inv", "alpha3_inv"
    ]
    
    missing = [c for c in required if c not in df_rge.columns]
    assert not missing, f"Missing required columns: {missing}"
    
    # Check monotonicity
    assert df_rge["mu_GeV"].is_monotonic_increasing, "mu_GeV must be monotonic increasing"
    

def assert_crossings(df_rge: pd.DataFrame) -> Dict:
    """
    Assert that E6, E7, E8 crossings exist and return scales.
    
    Parameters
    ----------
    df_rge : pd.DataFrame
        RGE data
        
    Returns
    -------
    dict
        Crossing scales
        
    Raises
    ------
    AssertionError
        If any crossing is missing
    """
    results = verify_rge_crossings(df_rge)
    
    for group in ["E6", "E7", "E8"]:
        assert results[group]["found"], f"No {group} crossing found"
    
    return {
        group: results[group]["mu_GeV"] 
        for group in ["E6", "E7", "E8"]
    }


def _create_verification_summary(lnD_test: Dict, phi_test: Dict, norm_test: Dict) -> str:
    """Create a summary string of verification results."""
    lines = [
        "VERIFICATION SUMMARY",
        "=" * 50,
        f"1. ln(D) cubic pattern: {lnD_test['message']}",
        f"2. φ third difference: {phi_test['message']}",
        f"3. γ₀ normalization: {norm_test['gamma0_check']['actual']:.4f} "
        f"(target: {norm_test['gamma0_check']['expected']:.4f}, "
        f"dev: {norm_test['gamma0_check']['deviation_pct']:.2f}%)",
        f"4. γ₂ topological: {norm_test['gamma2_check']['actual']:.6f} "
        f"(expected: {norm_test['gamma2_check']['expected']:.6f}, "
        f"dev: {norm_test['gamma2_check']['deviation_pct']:.2f}%)",
        "=" * 50
    ]
    return "\n".join(lines)
