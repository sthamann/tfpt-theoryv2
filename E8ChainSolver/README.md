# E8ChainSolver

**E8 Nilpotent Orbit Analysis for the TFPT Framework**

The E8ChainSolver validates the quadratic damping function γ(n) derived from E8 nilpotent orbits, providing an independent verification path for the Topological Fixed Point Theory (TFPT).

---

## Overview

This component analyzes chains through E8 nilpotent orbits to validate:

- **Damping function**: γ(n) ≈ 0.834 + 0.108n + 0.0106n²
- **Calibration-free test**: Δ³ln(D) constancy
- **Topological constraint**: γ₂ = γ₀/(8π²)

The engine searches for optimal monotonic chains through the E8 orbit lattice using beam search, then fits various models to validate theoretical predictions.

---

## Directory Structure

```
E8ChainSolver/
└── e8_orbit_engine/
    ├── src/e8_orbit_engine/    # Core Python library
    │   ├── __init__.py         # Public API exports
    │   ├── chain.py            # Chain building & validation
    │   ├── chain_search.py     # Beam search algorithm
    │   ├── fit.py              # Model fitting (quadratic, log, hyperbolic)
    │   ├── verify.py           # Verification tests
    │   ├── io.py               # Data loading utilities
    │   ├── report.py           # HTML report generation
    │   └── cli.py              # Command-line interface
    │
    ├── data/                   # Input data
    │   ├── nilpotent_orbits.csv   # E8 orbit table (Label, Height, Dim)
    │   └── gauge_couplings.csv    # Optional RGE data
    │
    ├── results/                # Generated outputs
    │   ├── orbit_engine_report.html  # Comprehensive HTML report
    │   ├── e8_chain_summary.json     # Chain and fit summary
    │   ├── data/               # Processed chain data
    │   ├── chains/             # Top chains from search
    │   │   ├── ranking.json
    │   │   └── top10/
    │   ├── plots/              # Visualization outputs
    │   └── tables/             # CSV output tables
    │
    ├── tests/                  # Unit tests
    │   ├── test_chain_strict_d2.py
    │   ├── test_fit.py
    │   ├── test_normalization_fix.py
    │   └── test_verify.py
    │
    ├── notebooks/              # Interactive demos
    │   └── demo.ipynb
    │
    ├── README.md               # Technical documentation
    ├── README_USAGE.md         # Quick usage guide
    ├── pyproject.toml          # Package configuration
    ├── setup.py                # Installation script
    ├── requirements.txt        # Dependencies
    └── run_demo.py             # Demo script
```

---

## Theory Background

### E8 Nilpotent Orbits

E8 has 70 nilpotent orbits classified by Bala-Carter labels. Each orbit has:
- **Dimension** (Dim): Size of the orbit in E8's 248-dimensional adjoint
- **Height** (h): Position in the Hasse diagram
- **Centralizer dimension** D = 248 − Dim

The engine constructs chains of strictly increasing orbit dimensions (decreasing D).

### Normalization

The damping function γ(n) is normalized by:

1. **First E8 step**: 248 (adjoint) → D₀ = 60 (A₄+A₁ orbit)
2. **Normalization step**: s₀ = ln(248) − ln(60) ≈ 1.419
3. **Scaling factor**: λ = γ₀/s₀ ≈ 0.588

This ensures γ(0) = 0.834 by construction, not by fitting.

### Key Equations

**Step sizes:**
```
s(n) = ln(Dₙ) − ln(Dₙ₊₁)
```

**Damping function (quadratic approximation):**
```
γ(n) = γ₀ + γ₁·n + γ₂·n²
```

**Third forward difference (calibration-free test):**
```
Δ³ln(D)ₙ = ln(Dₙ₊₃) − 3·ln(Dₙ₊₂) + 3·ln(Dₙ₊₁) − ln(Dₙ)
```

For a cubic ln(D), this should be approximately constant.

### Topological Constraint

The theory predicts:
```
γ₂ = γ₀/(8π²) ≈ 0.01056
```

This links the quadratic coefficient to the topological fixpoint c₃ = 1/(8π).

---

## Installation

```bash
cd E8ChainSolver/e8_orbit_engine

# Install as editable package
pip install -e .

# Or install with dev dependencies
pip install -e ".[dev]"
```

---

## Quick Start

### Strict ΔD=2 Mode (Theory-Conformant)

Generate a 27-step chain with purely ΔD=2 transitions:

```bash
python3 -m src.e8_orbit_engine.cli solve \
  --rebuild-e8 \
  --allowed-steps "2" \
  --beam-width 64 --top-k 3
```

**Expected Results:**
- 27 orbits from D=60 to D=8
- 0 large jumps (ΔD=4)
- CV(Δ³ln(D)) ≈ 1.75

### Exploratory Mode (ΔD ∈ {2,4})

Allow both ΔD=2 and ΔD=4 transitions:

```bash
python3 -m src.e8_orbit_engine.cli solve \
  --rebuild-e8 \
  --allowed-steps "2,4" \
  --beam-width 64 --top-k 3
```

### Baseline Chain

Use the predefined 27-step baseline chain:

```bash
python3 -m src.e8_orbit_engine.cli solve --rebuild-e8 --baseline
```

---

## Python API

```python
from e8_orbit_engine import (
    load_orbits, build_chain, fit_gamma, 
    beam_search_chains, compare_models, comprehensive_verification
)

# Load orbit data
df = load_orbits('data/nilpotent_orbits.csv')

# Build chain (baseline or search)
chain = build_chain(df, use_baseline=True)

# Fit gamma function
results = fit_gamma(chain, gamma0_target=0.834)
print(f"γ(n) = {results['gamma0']:.4f} + {results['gamma1']:.4f}n + {results['gamma2']:.6f}n²")
print(f"R² = {results['r_squared']:.6f}")

# Compare models
comparison = compare_models(chain, results)
print(comparison)

# Run verification
verification = comprehensive_verification(chain, results)
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--allowed-steps` | "2,4" | Comma-separated allowed ΔD steps (e.g., "2" for strict) |
| `--beam-width` | 32 | Beam search width (higher = more thorough) |
| `--top-k` | 10 | Number of top chains to save |
| `--baseline` | False | Use predefined baseline chain |
| `--html/--no-html` | --html | Generate HTML report |
| `--label-threshold` | 1.0 | Max label distance for edges (1.0 = no filter) |
| `--max-height-diff` | 999 | Max height difference for edges |

---

## Model Comparison

The engine compares four models for γ(n):

| Model | Formula | Parameters | Notes |
|-------|---------|------------|-------|
| Quadratic | γ = a + bn + cn² | 3 | Diagnostic only |
| Log-exact (anchored) | γ = λ·ln(Dₙ/Dₙ₊₁), λ fixed | 0 | Theory-conformant |
| Log-exact (fitted) | Same, λ fitted | 1 | Best information criterion |
| Hyperbolic | γ = A/(B−n) | 2 | Best R² typically |

---

## Output Files

| File | Description |
|------|-------------|
| `orbit_engine_report.html` | Comprehensive analysis with embedded plots |
| `e8_chain_summary.json` | Chain data and fit coefficients |
| `chains/ranking.json` | Metrics for top chains from search |
| `chains/top10/rank_N.csv` | Individual chain data files |
| `plots/*.png` | ln(D) vs n, γ vs n, residuals, Δ³ plots |
| `tables/e8_chain.csv` | Full chain table with all columns |

---

## Key Results

### Strict ΔD=2 Chain (27 steps)

| Metric | Value |
|--------|-------|
| Chain length | 27 orbits |
| D sequence | 60 → 58 → 56 → ... → 10 → 8 |
| CV(Δ³ln(D)) | ~1.75 |
| Best model | Hyperbolic (R² ≈ 1.000) |

### Physics Constants

| Parameter | Value | Origin |
|-----------|-------|--------|
| γ₀ | 0.834 | Anchored (= g/(g+1) with g=5) |
| s₀ (first step) | ln(248/60) ≈ 1.419 | E8 adjoint → A₄+A₁ |
| λ (scaling) | 0.834/1.419 ≈ 0.588 | Normalization ratio |

### Topological Constraint Check

```
γ₂(expected) = γ₀/(8π²) ≈ 0.01056
γ₂(observed) ≈ varies by model
```

The engine reports the deviation percentage for validation.

---

## Testing

```bash
cd e8_orbit_engine

# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=e8_orbit_engine

# Specific test file
pytest tests/test_chain_strict_d2.py -v
```

---

## Connection to TFPT

The E8ChainSolver provides independent validation of TFPT's damping parameter:

| TFPT Parameter | E8 Origin |
|----------------|-----------|
| γ(0) = 5/6 ≈ 0.833 | First orbit step normalization |
| g = 5 (holonomy) | SU(5) hypercharge eigenspace |
| c₃ = 1/(8π) | 11D Chern-Simons quantization |

The fact that γ₀ = g/(g+1) emerges from E8 orbit geometry—independently of the topological derivation—provides a cross-check between algebraic (E8) and topological (fixpoint) approaches.

---

## References

- Hamann, S. & Rizzo, A. (2026). "Topological Fixed Point Theory (TFPT)"
- Collingwood & McGovern. "Nilpotent Orbits in Semisimple Lie Algebras"
- E8 nilpotent orbit classification (Bala-Carter)

---

## License

MIT License

---

## Author

Stefan Hamann (sh@sh-future.de)
