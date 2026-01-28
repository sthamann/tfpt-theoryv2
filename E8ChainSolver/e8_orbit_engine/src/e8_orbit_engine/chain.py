"""
Chain building module for constructing monotonic E8 orbit sequences.
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Tuple
import logging

logger = logging.getLogger(__name__)


# Baseline E8 orbit order from the nilpotent orbit table
# This is a known valid chain, used as baseline for comparison
# Source: Nilpotent Orbits in Type E8 table
BASELINE_E8_ORDER = [
    "A4+A1",      # Dim=188, D=60
    "D5(a1)",     # Dim=190, D=58
    "A4+2A1",     # Dim=192, D=56
    "A4+A2",      # Dim=194, D=54
    "A5",         # Dim=196, D=52
    "D4+A2",      # Dim=198, D=50
    "D5",         # Dim=200, D=48
    "A5+A1",      # Dim=202, D=46
    "D5(a1)+A2",  # Dim=204, D=44
    "E6(a3)+A1",  # Dim=206, D=42
    "E8(a7)",     # Dim=208, D=40
    "A6",         # Dim=210, D=38
    "A6+A1",      # Dim=212, D=36
    "D5+A2",      # Dim=214, D=34
    "E6",         # Dim=216, D=32
    "A7",         # Dim=218, D=30
    "E7(a3)",     # Dim=220, D=28
    "D7(a1)",     # Dim=222, D=26
    "E7(a2)",     # Dim=224, D=24
    "D7",         # Dim=226, D=22
    "E7(a1)",     # Dim=228, D=20
    "E8(b4)",     # Dim=230, D=18
    "E8(a4)",     # Dim=232, D=16
    "E8(a3)",     # Dim=234, D=14
    "E8(a2)",     # Dim=236, D=12
    "E8(a1)",     # Dim=238, D=10
    "E8"          # Dim=240, D=8
]

def build_chain(df_orbits: pd.DataFrame, use_baseline: bool = False, custom_order: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Build the E8 orbit chain.
    
    Parameters
    ----------
    df_orbits : pd.DataFrame
        DataFrame with orbit data (Label, Dim, D columns)
    use_baseline : bool
        If True, use the baseline predefined order
    custom_order : List[str], optional
        Custom list of labels to use for the chain
        
    Returns
    -------
    pd.DataFrame
        Chain DataFrame with columns: n, label, dim_orbit, D, lnD
        
    Notes
    -----
    If neither use_baseline nor custom_order is provided, will search for
    optimal chain using beam search algorithm.
    """
    # Determine which order to use
    if custom_order is not None:
        order = custom_order
        logger.info("Using custom chain order")
    elif use_baseline:
        order = BASELINE_E8_ORDER
        logger.info("Using baseline E8 chain order")
    else:
        # Use beam search to find optimal chain
        logger.info("Searching for optimal chain...")
        from .chain_search import beam_search_chains
        
        # Search for chains
        chains = beam_search_chains(df_orbits, beam_width=32, top_k=1)
        
        if not chains:
            logger.warning("No valid chains found, falling back to baseline")
            order = BASELINE_E8_ORDER
        else:
            best_chain = chains[0]
            order = best_chain['labels']
            logger.info(f"Found optimal chain with cost {best_chain['cost']:.4f}")
            logger.info(f"CV(Δ³ln(D)): {best_chain['metrics']['cv_d3']:.4f}")
    
    # Set index to Label for easy lookup
    df = df_orbits.set_index("Label")
    
    rows = []
    for i, lbl in enumerate(order):
        if lbl not in df.index:
            logger.warning(f"Orbit {lbl} not found in data, skipping")
            continue
            
        dim = int(df.loc[lbl, "Dim"])
        D = int(248 - dim)
        rows.append({
            "n": i,
            "label": lbl,
            "dim_orbit": dim,
            "D": D
        })
    
    # Create chain DataFrame
    chain = pd.DataFrame(rows)
    
    # Validate monotonicity
    assert chain["dim_orbit"].is_monotonic_increasing, "Orbit dimension must be strictly increasing"
    assert chain["D"].is_monotonic_decreasing, "Centralizer dimension D must be strictly decreasing"
    
    # Add ln(D) for calculations
    chain["lnD"] = np.log(chain["D"].astype(float))
    
    logger.info(f"Built E8 chain with {len(chain)} orbits")
    logger.info(f"D range: {chain['D'].max()} → {chain['D'].min()}")
    
    return chain


def assert_integer_D(df_chain: pd.DataFrame) -> None:
    """
    Assert that all D values are integers.
    
    Parameters
    ----------
    df_chain : pd.DataFrame
        Chain DataFrame with D column
    """
    D_values = df_chain["D"].values
    assert np.all(D_values == D_values.astype(int)), "All D values must be integers"
    

def assert_monotone(df_chain: pd.DataFrame) -> None:
    """
    Assert that the chain is monotonic.
    
    Parameters
    ----------
    df_chain : pd.DataFrame
        Chain DataFrame with dim_orbit and D columns
    """
    assert df_chain["dim_orbit"].is_monotonic_increasing, "Orbit dimension must be increasing"
    assert df_chain["D"].is_monotonic_decreasing, "D must be decreasing"
    

def assert_allowed_labels(df_chain: pd.DataFrame) -> None:
    """
    Assert that all labels are from allowed families.
    
    Parameters
    ----------
    df_chain : pd.DataFrame
        Chain DataFrame with label column
    """
    ALLOWED = ("A", "D", "E", "E8")
    bad_labels = []
    for label in df_chain["label"]:
        # Check if label contains any allowed family
        # Labels can have numerical prefixes and combinations
        label_parts = label.replace("+", " ").replace("(", " ").replace(")", " ").split()
        has_valid_family = False
        for part in label_parts:
            # Remove leading digits
            part_no_digit = part.lstrip("0123456789")
            if any(part_no_digit.startswith(fam) for fam in ALLOWED):
                has_valid_family = True
                break
        
        if not has_valid_family:
            bad_labels.append(label)
    
    if bad_labels:
        raise ValueError(f"Invalid labels found: {bad_labels}")


def validate_chain(chain: pd.DataFrame) -> Tuple[bool, List[str]]:
    """
    Validate that chain satisfies monotonicity and other requirements.
    
    Returns
    -------
    valid : bool
        Whether chain is valid
    messages : list of str
        Validation messages
    """
    messages = []
    valid = True
    
    try:
        # Check integer D values
        assert_integer_D(chain)
        messages.append("✓ D values are all integers")
    except AssertionError as e:
        messages.append(f"✗ {str(e)}")
        valid = False
    
    try:
        # Check monotonicity
        assert_monotone(chain)
        messages.append("✓ Chain is monotonic (dim_orbit ↑, D ↓)")
    except AssertionError as e:
        messages.append(f"✗ {str(e)}")
        valid = False
    
    try:
        # Check allowed labels
        assert_allowed_labels(chain)
        messages.append("✓ All labels are from allowed families (A, D, E, E8)")
    except ValueError as e:
        messages.append(f"✗ {str(e)}")
        valid = False
    
    # Check reasonable range
    if len(chain) < 5:
        messages.append(f"✗ Chain too short: {len(chain)} < 5")
        valid = False
    else:
        messages.append(f"✓ Chain length: {len(chain)} orbits")
    
    # Check for valid D values
    D_values = chain['D'].values
    if np.any(D_values <= 0) or np.any(D_values > 248):
        messages.append("✗ Invalid D values (must be in (0, 248])")
        valid = False
    else:
        messages.append(f"✓ D range: {D_values.max()} → {D_values.min()}")
    
    # Summary
    if valid:
        messages.append("✅ Chain validation PASSED")
    else:
        messages.append("❌ Chain validation FAILED")
    
    return valid, messages
