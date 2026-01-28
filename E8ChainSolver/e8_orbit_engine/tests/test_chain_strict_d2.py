"""
Test strict ΔD=2 chain properties for regression safety.
"""
import numpy as np
from src.e8_orbit_engine.io import load_orbits
from src.e8_orbit_engine.chain_search import beam_search_chains
from src.e8_orbit_engine.chain import build_chain
from src.e8_orbit_engine.fit import compute_gamma, assert_first_step


def test_strict_d2_chain_has_27_steps_and_no_jumps():
    """Test that strict ΔD=2 mode produces a 27-step chain with no large jumps."""
    df = load_orbits("data/nilpotent_orbits.csv")
    chains = beam_search_chains(
        df, beam_width=64, top_k=1, allowed_steps=(2,), verbose=False
    )
    assert chains, "No chain found"
    best = chains[0]
    
    # Check chain properties
    assert len(best["labels"]) == 27, f"Expected 27 steps, got {len(best['labels'])}"
    assert best["metrics"]["n_jumps"] == 0, f"Expected 0 jumps, got {best['metrics']['n_jumps']}"
    assert best["D_values"][0] == 60, f"Expected start D=60, got {best['D_values'][0]}"
    assert best["D_values"][-1] == 8, f"Expected end D=8, got {best['D_values'][-1]}"
    
    # Verify ΔD=2 throughout
    D_values = best["D_values"]
    for i in range(len(D_values) - 1):
        delta_D = D_values[i] - D_values[i+1]
        assert delta_D == 2, f"Step {i}: ΔD={delta_D}, expected 2"


def test_normalization_anchor_gamma0():
    """Test that gamma normalization produces γ₀=0.834 exactly."""
    df = load_orbits("data/nilpotent_orbits.csv")
    
    # Use baseline chain for this test (or strict ΔD=2)
    chain = build_chain(df, use_baseline=True)
    s, gamma, lam, s_star = compute_gamma(chain)
    
    # Check γ₀ anchor
    assert np.isclose(gamma[0], 0.834, rtol=0, atol=1e-12), \
        f"γ₀ = {gamma[0]}, expected 0.834"
    
    # Check first step normalization
    assert_first_step(chain, tolerance=1e-6)
    
    # Verify s* calculation
    expected_s_star = np.log(248.0) - np.log(60.0)
    assert np.isclose(s_star, expected_s_star, rtol=1e-10), \
        f"s* = {s_star}, expected {expected_s_star}"


def test_a4_a1_has_corrected_dimension():
    """Test that A4+A1 has corrected dimension 188 (D=60)."""
    df = load_orbits("data/nilpotent_orbits.csv")
    
    # Find A4+A1 orbit
    a4_a1 = df[df['Label'].str.replace(' ', '') == 'A4+A1']
    assert len(a4_a1) > 0, "A4+A1 orbit not found"
    
    # Check corrected dimension
    dim = a4_a1.iloc[0]['Dim']
    assert dim == 188, f"A4+A1 dimension = {dim}, expected 188 (corrected from 206)"
    
    # Check corresponding D value
    D = 248 - dim
    assert D == 60, f"A4+A1 has D={D}, expected 60"


def test_strict_d2_chain_is_monotonic():
    """Test that strict ΔD=2 chain has monotonic properties."""
    df = load_orbits("data/nilpotent_orbits.csv")
    chains = beam_search_chains(
        df, beam_width=32, top_k=1, allowed_steps=(2,), verbose=False
    )
    assert chains, "No chain found"
    best = chains[0]
    
    # Check D monotonically decreasing
    D_values = best["D_values"]
    for i in range(len(D_values) - 1):
        assert D_values[i] > D_values[i+1], \
            f"D not decreasing at step {i}: {D_values[i]} -> {D_values[i+1]}"
    
    # Check dimensions monotonically increasing
    dims = best["dims"]
    for i in range(len(dims) - 1):
        assert dims[i] < dims[i+1], \
            f"Dimensions not increasing at step {i}: {dims[i]} -> {dims[i+1]}"


def test_baseline_chain_properties():
    """Test that baseline chain has expected properties."""
    df = load_orbits("data/nilpotent_orbits.csv")
    chain = build_chain(df, use_baseline=True)
    
    # Check length
    assert len(chain) == 27, f"Baseline chain length = {len(chain)}, expected 27"
    
    # Check D range
    assert chain['D'].iloc[0] == 60, f"Start D = {chain['D'].iloc[0]}, expected 60"
    assert chain['D'].iloc[-1] == 8, f"End D = {chain['D'].iloc[-1]}, expected 8"
    
    # Check all ΔD = 2
    D_diffs = -chain['D'].diff().dropna()
    assert all(D_diffs == 2), f"Not all ΔD=2 in baseline chain: {list(D_diffs.unique())}"
    
    # Check first orbit is A4+A1
    assert chain['label'].iloc[0] == 'A4+A1', \
        f"First orbit = {chain['label'].iloc[0]}, expected A4+A1"
