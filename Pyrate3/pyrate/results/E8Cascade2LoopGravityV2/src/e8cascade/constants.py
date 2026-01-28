"""Physical constants and model parameters for E8 Cascade"""
import numpy as np

# Physical scales
M_Pl = 1.22e19   # GeV - Planck mass
M_GUT = 2e17     # GeV - E8 breaking scale
M_Z = 91.1876    # GeV - Z boson mass

# Threshold scales
M_SIGMA = 1e3    # GeV - SigmaF triplet mass
M_NR = 1e15      # GeV - Right-handed neutrino mass  
M_PHI = 1e16     # GeV - PQ scalar mass

# Beta coefficients at different energy regimes
# SM values (below all thresholds)
BETA_SM = {
    'g1': 41/6,    # U(1)_Y with GUT normalization
    'g2': -11/6,   # SU(2)_L
    'g3': -7       # SU(3)_C
}

# E8 cascade values (above all thresholds)
# These include contributions from SigmaF triplet, extra quarks, etc.
BETA_E8 = {
    'g1': 50/6,    # Enhanced by extra fields
    'g2': 2,       # Positive due to SigmaF triplet!
    'g3': -3       # Less negative due to extra colored fields
}

# Intermediate values (for threshold matching)
BETA_ABOVE_SIGMA = {
    'g1': 44/6,    # Partial enhancement
    'g2': 1/2,     # Already turning positive
    'g3': -6       # Slightly less negative
}

BETA_ABOVE_NR = {
    'g1': 47/6,    
    'g2': 3/2,     
    'g3': -4       
}

def gravity_coefficient(b_i, M_GUT=M_GUT, M_Pl=M_Pl):
    """
    Compute gravity portal coefficient
    c_Ri = b_i/(16π²) * (M_GUT/M_Pl)²
    
    Note: Keeps the sign of b_i!
    """
    return b_i / (16 * np.pi**2) * (M_GUT / M_Pl)**2

def get_beta_coefficients(mu):
    """
    Get appropriate beta coefficients for given energy scale
    Implements threshold matching
    """
    if mu < M_SIGMA:
        return BETA_SM
    elif mu < M_NR:
        return BETA_ABOVE_SIGMA
    elif mu < M_PHI:
        return BETA_ABOVE_NR
    else:
        return BETA_E8 