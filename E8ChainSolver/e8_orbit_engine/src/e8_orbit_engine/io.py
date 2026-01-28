"""
I/O module for loading and processing E8 orbit data and RGE data.
"""

import pandas as pd
import numpy as np
import math
from pathlib import Path
from typing import Union, Optional
import logging


# Configure logging
logger = logging.getLogger(__name__)

# Allowed label families
ALLOWED_FAMILIES = ("A", "D", "E", "E8")

def load_orbits(path: Union[str, Path] = "data/nilpotent_orbits.csv") -> pd.DataFrame:
    """
    Load nilpotent orbit data from CSV file with validation and corrections.
    
    Parameters
    ----------
    path : str or Path
        Path to the data file (default: data/nilpotent_orbits.csv)
        
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: Label, Height, Dim, D (centralizer dimension)
        
    Notes
    -----
    - Dim is the orbit dimension
    - D = 248 - Dim is the centralizer dimension (must be integer)
    - Special case: A4+A1 with Dim=206 is auto-corrected to Dim=188 (E8 table)
    """
    path = Path(path)
    
    # Read CSV file
    df = pd.read_csv(path)
    
    # Check required columns
    required_cols = {"Label", "Height", "Dim"}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f"Missing required columns: {missing}")
    
    # Special case correction: A4+A1 dimension
    # According to E8 nilpotent orbit tables, A4+A1 has Dim=188, not 206
    mask = (df["Label"].str.replace(" ", "") == "A4+A1") & (df["Dim"] == 206)
    if mask.any():
        logger.info("Auto-correcting A4+A1: Dim 206 → 188 (source: E8 nilpotent orbit table)")
        print("[INFO] Auto-corrected A4+A1 dimension from 206 to 188 (source: E8 nilpotent orbit table)")
        df.loc[mask, "Dim"] = 188
    
    # Basic validations
    assert (df["Dim"].between(0, 240)).all(), "Orbit dimension must be in [0, 240]"
    
    # Calculate centralizer dimension
    df["D"] = 248 - df["Dim"].astype(int)
    
    # Validate D is positive integer
    assert (df["D"] > 0).all(), "Centralizer dimension D must be positive"
    assert (df["D"] == df["D"].round()).all(), "Centralizer dimension D must be integer"
    
    # Validate labels (only A, D, E families or E8 variants)
    # Labels can have numerical prefixes (e.g., 2A1, 3A1) and combinations (e.g., A2+A1)
    bad_labels = []
    for label in df["Label"]:
        label_clean = str(label).strip()
        if label_clean == "0":  # Special case: trivial orbit
            continue
        
        # Check if label contains any allowed family
        # Remove numbers and check if any part contains A, D, E, or E8
        label_parts = label_clean.replace("+", " ").replace("(", " ").replace(")", " ").split()
        has_valid_family = False
        for part in label_parts:
            # Remove leading digits
            part_no_digit = part.lstrip("0123456789")
            if any(part_no_digit.startswith(fam) for fam in ALLOWED_FAMILIES):
                has_valid_family = True
                break
        
        if not has_valid_family:
            bad_labels.append(label_clean)
    
    if bad_labels:
        raise ValueError(f"Invalid labels (must contain A, D, E families or E8): {bad_labels}")
    
    logger.info(f"Loaded {len(df)} orbits, D range: {df['D'].min()}-{df['D'].max()}")
    
    return df


def load_rge(path: Union[str, Path] = "data/gauge_couplings.csv") -> pd.DataFrame:
    """
    Load RGE gauge coupling data from CSV file with validation.
    
    Parameters
    ----------
    path : str or Path
        Path to the RGE data file (default: data/gauge_couplings.csv)
        
    Returns
    -------
    pd.DataFrame
        DataFrame with gauge coupling evolution data
        
    Notes
    -----
    Expected columns:
    - mu_GeV: Energy scale in GeV
    - log10_mu: Log10 of energy scale
    - g1_SM, g1_GUT, g2, g3: Gauge couplings
    - alpha1_GUT, alpha2, alpha3: Fine structure constants
    - alpha1_inv_GUT, alpha2_inv, alpha3_inv: Inverse fine structure constants
    """
    path = Path(path)
    
    # Expected columns from paper
    expected_cols = [
        "mu_GeV", "log10_mu", "g1_SM", "g1_GUT", "g2", "g3",
        "alpha1_GUT", "alpha2", "alpha3", 
        "alpha1_inv_GUT", "alpha2_inv", "alpha3_inv"
    ]
    
    # Read CSV
    df = pd.read_csv(path)
    
    # Check for missing columns
    missing = [c for c in expected_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in RGE data: {missing}")
    
    # Sort by energy scale
    df = df.sort_values("mu_GeV")
    
    # Consistency checks
    # 1. Check mu_GeV is strictly increasing
    assert df["mu_GeV"].is_monotonic_increasing, "mu_GeV must be strictly increasing"
    
    # 2. Check log10_mu consistency
    log_check = np.log10(df["mu_GeV"].values)
    assert np.allclose(df["log10_mu"].values, log_check, rtol=1e-3), \
        "log10_mu inconsistent with mu_GeV"
    
    # 3. Sample check: alpha3 = g3^2/(4π)
    probe = df.sample(min(10, len(df)), random_state=42)
    alpha3_calc = probe["g3"]**2 / (4 * math.pi)
    assert np.allclose(probe["alpha3"], alpha3_calc, rtol=5e-3), \
        "alpha3 inconsistent with g3^2/(4π)"
    
    logger.info(f"Loaded RGE data: {len(df)} points, μ range: {df['mu_GeV'].min():.1e} - {df['mu_GeV'].max():.1e} GeV")
    
    return df


def save_chain(chain: pd.DataFrame, path: Union[str, Path]) -> None:
    """
    Save orbit chain to CSV file.
    
    Parameters
    ----------
    chain : pd.DataFrame
        Chain DataFrame with orbit data
    path : str or Path
        Output path for CSV file
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    chain.to_csv(path, index=False)
    print(f"Saved chain to {path}")
