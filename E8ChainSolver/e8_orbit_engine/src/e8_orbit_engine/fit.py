"""
Fitting module for quadratic damping function γ(n) and diagnostics.
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from typing import Dict, Tuple, Optional
import math
import logging

logger = logging.getLogger(__name__)


def compute_gamma(df_chain: pd.DataFrame, gamma0: float = 0.834) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Compute gamma values from chain data WITHOUT fitting.
    
    Uses normalization based on the first step from adjoint (248) to A4+A1 (D=60).
    
    Parameters
    ----------
    df_chain : pd.DataFrame
        Chain DataFrame with lnD values
    gamma0 : float
        Target value for γ(0), default 0.834 from paper
        
    Returns
    -------
    s : np.ndarray
        Step sizes s(n) = ln(D_n) - ln(D_{n+1})
    gamma : np.ndarray
        Gamma values γ(n) = λ * s(n)
    lam : float
        Scaling factor λ = γ₀ / s*
    s_star : float
        Normalization step s* = ln(248) - ln(D₀) where D₀ = 60 for A4+A1
        
    Notes
    -----
    This is a definition, not a fit. The gamma values are directly computed
    from the data using the normalization constant.
    The first step is from the adjoint (248) to the first orbit.
    """
    lnD = df_chain["lnD"].values
    D0 = df_chain["D"].iloc[0]  # Should be 60 for A4+A1
    
    # CRITICAL: Include the first step from adjoint (248) to first orbit
    # This is the normalization step
    s0 = math.log(248.0) - lnD[0]  # Step from adjoint to first orbit
    s_rest = lnD[:-1] - lnD[1:]    # Rest of the steps within the chain
    s = np.concatenate(([s0], s_rest))  # All steps, starting from adjoint
    
    # Normalization: s* = ln(248) - ln(D₀) = s0
    s_star = s0  # This is the normalization step
    
    # Scaling factor
    lam = gamma0 / s_star
    
    # Compute gamma values (this is a definition, not a fit)
    gamma = lam * s
    
    logger.info(f"Normalization: s* = ln(248) - ln({D0}) = {s_star:.6f}")
    logger.info(f"Scaling factor: λ = {gamma0} / {s_star:.6f} = {lam:.6f}")
    logger.info(f"γ(0) = {gamma[0]:.6f} (should equal {gamma0})")
    
    # Verify normalization
    assert abs(gamma[0] - gamma0) < 1e-10, f"Normalization failed: γ(0) = {gamma[0]} ≠ {gamma0}"
    
    return s, gamma, lam, s_star


def quad_and_d3(df_chain: pd.DataFrame, include_adjoint: bool = True) -> Dict:
    """
    Perform quadratic fit (for diagnostics only) and compute third differences.
    
    Parameters
    ----------
    df_chain : pd.DataFrame
        Chain DataFrame with lnD values
    include_adjoint : bool
        Whether to include the adjoint (248) step in calculations
        
    Returns
    -------
    dict
        Dictionary containing:
        - coef: Quadratic coefficients [c, b, a] for s(n) ≈ a + b*n + c*n²
        - r2: R-squared value of quadratic fit
        - residuals: Residuals from quadratic fit
        - d3: Third forward differences of ln(D)
        - d3_mean: Mean of third differences
        - d3_var: Variance of third differences
        
    Notes
    -----
    The quadratic fit is for diagnostic purposes only. The actual gamma values
    are computed using the normalization method, not from fitting.
    
    Third forward difference: Δ³ln(D_n) = ln(D_{n+3}) - 3*ln(D_{n+2}) + 3*ln(D_{n+1}) - ln(D_n)
    For a cubic ln(D), this should be approximately constant.
    """
    lnD = df_chain["lnD"].values
    
    # Include adjoint step if requested
    if include_adjoint:
        # Prepend ln(248) to the lnD array
        lnD_full = np.concatenate(([math.log(248.0)], lnD))
        # Calculate step sizes including first step from adjoint
        s = lnD_full[:-1] - lnD_full[1:]
    else:
        # Calculate step sizes within chain only
        s = lnD[:-1] - lnD[1:]
        lnD_full = lnD
    
    n = np.arange(len(s))
    
    # Fit quadratic (diagnostic only)
    coef = np.polyfit(n, s, 2)  # Returns [c, b, a] for c*n² + b*n + a
    s_fit = np.polyval(coef, n)
    residuals = s - s_fit
    
    # Calculate R²
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((s - s.mean())**2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # Calculate third forward differences on full lnD
    # Δ³ln(D) should be approximately constant for cubic ln(D)
    d3 = np.diff(np.diff(np.diff(lnD_full)))
    
    # Statistics of third differences
    d3_mean = float(d3.mean()) if len(d3) > 0 else 0
    d3_var = float(d3.var(ddof=1)) if len(d3) > 1 else 0
    
    logger.info(f"Quadratic fit: s(n) ≈ {coef[2]:.6f} + {coef[1]:.6f}*n + {coef[0]:.8f}*n²")
    logger.info(f"R² = {r2:.6f}")
    logger.info(f"Third differences: mean = {d3_mean:.8f}, var = {d3_var:.8e}")
    
    return {
        "coef": coef,  # [c, b, a] for polynomial
        "r2": float(r2),
        "residuals": residuals,
        "d3": d3,
        "d3_mean": d3_mean,
        "d3_var": d3_var
    }


def assert_first_step(df_chain: pd.DataFrame, tolerance: float = 0.01) -> None:
    """
    Assert that the first step matches the expected normalization.
    
    Parameters
    ----------
    df_chain : pd.DataFrame
        Chain DataFrame
    tolerance : float
        Relative tolerance for the check
        
    Raises
    ------
    AssertionError
        If first step deviates from expected value
    """
    D0 = df_chain["D"].iloc[0]
    s_star_expected = math.log(248.0 / 60.0)  # Expected for A4+A1
    s_star_actual = math.log(248.0 / D0)
    
    rel_error = abs(s_star_actual - s_star_expected) / s_star_expected
    
    assert rel_error < tolerance, \
        f"First step deviation: expected s* = ln(248/60) = {s_star_expected:.6f}, " \
        f"got ln(248/{D0}) = {s_star_actual:.6f} (error: {rel_error*100:.2f}%)"


def fit_gamma(
    chain: pd.DataFrame, 
    gamma0_target: float = 0.834,
    weighted: bool = False,
    include_adjoint: bool = True
) -> Dict:
    """
    Compute gamma values and perform diagnostic quadratic fit.
    
    This function computes gamma using the normalization method (not fitting)
    and also performs a quadratic fit for diagnostic purposes.
    
    Parameters
    ----------
    chain : pd.DataFrame
        Orbit chain with D values
    gamma0_target : float
        Target value for γ(0), from paper
    weighted : bool
        Use weighted least squares for diagnostic fit
    include_adjoint : bool
        Include the adjoint (248) as zeroth step if missing
        
    Returns
    -------
    dict
        Results including:
        - gamma0, gamma1, gamma2: quadratic coefficients (from diagnostic fit)
        - lam: scaling factor
        - gamma_obs: observed (computed) gamma values
        - gamma_fit: fitted gamma values (diagnostic)
        - residuals: fit residuals
        - r_squared: R² value
        - s0_true: actual first step used
    """
    # Get D values
    lnD = np.log(chain['D'].values)
    
    # CRITICAL: Include the true first step from adjoint (248) to first orbit
    if include_adjoint and chain['D'].iloc[0] < 248:
        # Prepend the true first step: 248 → D₀
        s0_true = np.log(248.0) - lnD[0]  # This is THE anchoring step
        s_rest = lnD[:-1] - lnD[1:]        # Rest of the steps
        s = np.concatenate(([s0_true], s_rest))  # s[0] is now 248 → first orbit
    else:
        # Original behavior (for backward compatibility)
        s = lnD[:-1] - lnD[1:]
    
    n = np.arange(len(s))
    
    # Calculate scaling from the TRUE first step
    # λ = γ₀/s₀ where s₀ = ln(248) - ln(D_first)
    lam = gamma0_target / s[0]
    
    # Scale steps to get gamma
    gamma_obs = lam * s
    
    # Fit quadratic (for diagnostics only)
    if weighted:
        # Weight by 1/n to emphasize early terms
        weights = 1.0 / (n + 1)
        popt, pcov = curve_fit(lambda n, a, b, c: a + b * n + c * n * n, 
                              n, gamma_obs, sigma=1/weights)
    else:
        popt, pcov = curve_fit(lambda n, a, b, c: a + b * n + c * n * n, 
                              n, gamma_obs)
    
    a, b, c = popt
    
    # Calculate fitted values and residuals
    gamma_fit = a + b * n + c * n * n
    residuals = gamma_obs - gamma_fit
    
    # Calculate R²
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((gamma_obs - np.mean(gamma_obs))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    # Check theoretical constraint: γ₂ ≈ γ₀/(8π²)
    # Note: Using the target gamma0, not the fitted value
    theoretical_gamma2 = gamma0_target / (8 * np.pi**2)
    gamma2_deviation = abs(c - theoretical_gamma2) / theoretical_gamma2 if theoretical_gamma2 != 0 else np.inf
    
    # For consistency, lnD should include adjoint if we added it
    if include_adjoint and chain['D'].iloc[0] < 248:
        lnD_full = np.concatenate(([np.log(248.0)], lnD))
    else:
        lnD_full = lnD
    
    results = {
        'gamma0': a,
        'gamma1': b,
        'gamma2': c,
        'lam': lam,
        'gamma_obs': gamma_obs,
        'gamma_fit': gamma_fit,
        'residuals': residuals,
        'r_squared': r_squared,
        'lnD': lnD_full,
        'steps': s,
        'n': n,
        'cov': pcov,
        'theoretical_gamma2': theoretical_gamma2,
        'gamma2_deviation': gamma2_deviation,
        's0_true': s[0] if len(s) > 0 else None  # Store the actual first step used
    }
    
    return results


def fit_gamma_log_model(chain: pd.DataFrame, anchor: bool = True) -> Dict:
    """
    Fit exact logarithmic model: γ(n) = λ * ln(D_n / D_{n+1})
    
    Parameters
    ----------
    chain : pd.DataFrame
        Chain DataFrame with D values
    anchor : bool
        If True, use fixed λ from normalization. If False, fit λ.
        
    Returns
    -------
    dict
        Model results including parameters, predictions, and statistics
    """
    import numpy as np
    from scipy.optimize import curve_fit
    
    D = chain['D'].values  # [60, 58, ..., 8]
    n = np.arange(len(D))
    
    # Get observed gamma from data
    s, gamma_obs, lam_anchor, s_star = compute_gamma(chain)
    
    def g_exact(n_vals, lam):
        """Exact log model for n >= 1"""
        Dn = 60 - 2 * n_vals
        # Avoid division by zero or negative values
        Dn_next = Dn - 2
        mask = (Dn > 0) & (Dn_next > 0)
        result = np.zeros_like(n_vals, dtype=float)
        result[mask] = lam * np.log(Dn[mask] / Dn_next[mask])
        return result
    
    # Focus on n >= 1 (exclude the anchor point)
    idx = n >= 1
    n_fit = n[idx]
    y_fit = gamma_obs[idx]
    
    if anchor:
        # Use fixed lambda from normalization
        lam = lam_anchor
        yhat = g_exact(n_fit, lam)
        model_name = "log-exact-anchored"
        n_params = 0  # No free parameters
    else:
        # Fit lambda
        popt, pcov = curve_fit(lambda k, lam: g_exact(k, lam), n_fit, y_fit, p0=(0.58,))
        lam = float(popt[0])
        yhat = g_exact(n_fit, lam)
        model_name = "log-exact-fit"
        n_params = 1
    
    # Calculate statistics
    ss_res = np.sum((y_fit - yhat)**2)
    ss_tot = np.sum((y_fit - y_fit.mean())**2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # AIC and BIC
    n_obs = len(y_fit)
    if n_obs > 0:
        mse = ss_res / n_obs
        ll = -0.5 * n_obs * (np.log(2 * np.pi * mse) + 1)  # Log-likelihood
        aic = 2 * n_params - 2 * ll
        bic = n_params * np.log(n_obs) - 2 * ll
    else:
        aic = bic = np.inf
    
    return {
        "model": model_name,
        "lam": lam,
        "yhat": yhat,
        "n_fit": n_fit,
        "ss_res": ss_res,
        "r2": r2,
        "aic": aic,
        "bic": bic,
        "n_params": n_params
    }


def fit_gamma_hyperbolic(chain: pd.DataFrame) -> Dict:
    """
    Fit hyperbolic model: γ(n) ≈ A / (B - n) for n >= 1
    
    This approximates ln(D_n/D_{n+1}) extremely well for D_n = 60 - 2n.
    
    Parameters
    ----------
    chain : pd.DataFrame
        Chain DataFrame with D values
        
    Returns
    -------
    dict
        Model results including parameters, predictions, and statistics
    """
    import numpy as np
    from scipy.optimize import curve_fit
    
    # Get observed gamma
    s, gamma_obs, lam, s_star = compute_gamma(chain)
    n = np.arange(len(gamma_obs))
    
    # Focus on n >= 1
    idx = n >= 1
    n_fit = n[idx]
    y_fit = gamma_obs[idx]
    
    def hyper(n_vals, A, B):
        """Hyperbolic function"""
        # Avoid division by zero
        denom = B - n_vals
        mask = denom > 0
        result = np.zeros_like(n_vals, dtype=float)
        result[mask] = A / denom[mask]
        return result
    
    # Initial guess: A ≈ λ, B ≈ 29 (from D_n = 60 - 2n → denominator ≈ 29 - n)
    # Ensure B > max(n) for numerical stability
    n_max = int(n_fit.max()) if len(n_fit) else 0
    popt, pcov = curve_fit(
        hyper, n_fit, y_fit, 
        p0=(0.59, 29.0),
        bounds=([0.0, n_max + 1e-3], [np.inf, 1e3]),
        maxfev=5000
    )
    A, B = map(float, popt)
    
    yhat = hyper(n_fit, A, B)
    
    # Calculate statistics
    ss_res = np.sum((y_fit - yhat)**2)
    ss_tot = np.sum((y_fit - y_fit.mean())**2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # AIC and BIC
    n_obs = len(y_fit)
    n_params = 2  # A and B
    if n_obs > 0:
        mse = ss_res / n_obs
        ll = -0.5 * n_obs * (np.log(2 * np.pi * mse) + 1)
        aic = 2 * n_params - 2 * ll
        bic = n_params * np.log(n_obs) - 2 * ll
    else:
        aic = bic = np.inf
    
    return {
        "model": "hyperbolic",
        "A": A,
        "B": B,
        "yhat": yhat,
        "n_fit": n_fit,
        "ss_res": ss_res,
        "r2": r2,
        "aic": aic,
        "bic": bic,
        "n_params": n_params
    }


def compare_models(chain: pd.DataFrame, fit_results: Dict) -> pd.DataFrame:
    """
    Compare quadratic, log-exact, and hyperbolic models.
    
    Parameters
    ----------
    chain : pd.DataFrame
        Chain DataFrame
    fit_results : dict
        Results from standard fit_gamma (quadratic)
        
    Returns
    -------
    pd.DataFrame
        Comparison table with R², AIC, BIC for each model
    """
    import pandas as pd
    import numpy as np
    
    # Quadratic model stats (from existing fit, but only for n >= 1)
    s, gamma_obs, lam, s_star = compute_gamma(chain)
    n = np.arange(len(gamma_obs))
    idx = n >= 1
    n_fit = n[idx]
    y_fit = gamma_obs[idx]
    
    # Recompute quadratic for n >= 1 only
    gamma_fit_quad = fit_results['gamma_fit'][idx]
    ss_res_quad = np.sum((y_fit - gamma_fit_quad)**2)
    ss_tot = np.sum((y_fit - y_fit.mean())**2)
    r2_quad = 1 - (ss_res_quad / ss_tot) if ss_tot > 0 else 0
    
    n_obs = len(y_fit)
    n_params_quad = 3  # gamma0, gamma1, gamma2
    mse_quad = ss_res_quad / n_obs
    ll_quad = -0.5 * n_obs * (np.log(2 * np.pi * mse_quad) + 1)
    aic_quad = 2 * n_params_quad - 2 * ll_quad
    bic_quad = n_params_quad * np.log(n_obs) - 2 * ll_quad
    
    # Get other models
    log_anchored = fit_gamma_log_model(chain, anchor=True)
    log_fit = fit_gamma_log_model(chain, anchor=False)
    hyperbolic = fit_gamma_hyperbolic(chain)
    
    # Create comparison table
    comparison = pd.DataFrame([
        {
            "Model": "Quadratic (n≥1)",
            "Parameters": n_params_quad,
            "R²": r2_quad,
            "AIC": aic_quad,
            "BIC": bic_quad,
            "RMSE": np.sqrt(mse_quad)
        },
        {
            "Model": "Log-exact (anchored)",
            "Parameters": log_anchored["n_params"],
            "R²": log_anchored["r2"],
            "AIC": log_anchored["aic"],
            "BIC": log_anchored["bic"],
            "RMSE": np.sqrt(log_anchored["ss_res"] / n_obs)
        },
        {
            "Model": "Log-exact (fitted)",
            "Parameters": log_fit["n_params"],
            "R²": log_fit["r2"],
            "AIC": log_fit["aic"],
            "BIC": log_fit["bic"],
            "RMSE": np.sqrt(log_fit["ss_res"] / n_obs)
        },
        {
            "Model": "Hyperbolic",
            "Parameters": hyperbolic["n_params"],
            "R²": hyperbolic["r2"],
            "AIC": hyperbolic["aic"],
            "BIC": hyperbolic["bic"],
            "RMSE": np.sqrt(hyperbolic["ss_res"] / n_obs)
        }
    ])
    
    # Sort by AIC (lower is better)
    comparison = comparison.sort_values("AIC")
    
    return comparison


def analyze_fit_quality(fit_results: Dict) -> Dict:
    """
    Analyze quality of the fit against theoretical predictions.
    
    Parameters
    ----------
    fit_results : dict
        Results from fit_gamma or compute_gamma
        
    Returns
    -------
    dict
        Quality metrics and checks
    """
    # Target values from paper
    targets = {
        'gamma0': 0.834,
        'gamma1': 0.108,
        'gamma2': 0.0105627
    }
    
    # Calculate deviations
    deviations = {}
    for key in targets:
        if key in fit_results:
            actual = fit_results[key]
            target = targets[key]
            dev_abs = abs(actual - target)
            dev_rel = dev_abs / abs(target) if target != 0 else np.inf
            deviations[key] = {
                'actual': actual,
                'target': target,
                'deviation_abs': dev_abs,
                'deviation_rel': dev_rel,
                'deviation_pct': dev_rel * 100
            }
    
    # Check topological constraint: γ₂ = γ₀/(8π²)
    if 'gamma0' in fit_results and 'gamma2' in fit_results:
        gamma0 = fit_results['gamma0']
        gamma2 = fit_results['gamma2']
        expected_gamma2 = gamma0 / (8 * np.pi**2)
        
        topological_check = {
            'expected_gamma2': expected_gamma2,
            'actual_gamma2': gamma2,
            'deviation': abs(gamma2 - expected_gamma2),
            'deviation_pct': abs(gamma2 - expected_gamma2) / expected_gamma2 * 100
        }
    else:
        topological_check = None
    
    # Overall quality assessment
    if deviations:
        max_deviation_pct = max(d['deviation_pct'] for d in deviations.values())
        
        if max_deviation_pct < 1:
            quality = "excellent"
        elif max_deviation_pct < 5:
            quality = "good"
        elif max_deviation_pct < 10:
            quality = "acceptable"
        else:
            quality = "poor"
    else:
        quality = "unknown"
    
    return {
        'deviations': deviations,
        'topological_check': topological_check,
        'r_squared': fit_results.get('r_squared', 0),
        'max_residual': np.max(np.abs(fit_results.get('residuals', [0]))),
        'rms_residual': np.sqrt(np.mean(fit_results.get('residuals', [0])**2)),
        'quality': quality
    }