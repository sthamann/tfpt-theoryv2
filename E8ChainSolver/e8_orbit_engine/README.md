# E8 Orbit Engine

A data-driven engine to derive the quadratic damping function γ(n) from nilpotent E8 orbits and validate calibration-free consistency tests.

## Theory Background

Based on the paper "Von Topologie zu Dynamik: Die Ordnung hinter α und den Naturkonstanten" by Stefan Hamann (2025).

The engine validates:
- γ(n) ≈ 0.834 + 0.108n + 0.0105627n²
- Calibration-free test: Δ³ln(φ) = -2γ₂
- Topological constraint: γ₂ = γ₀/(8π²)

## Installation

```bash
# Clone the repository
cd e8_orbit_engine

# Install dependencies
pip install -e .

# Or with development dependencies
pip install -e ".[dev]"
```

## Quick Start

```bash
# Run analysis on orbit data
e8-orbit fit data/nilpotent_orbits.csv

# Get help
e8-orbit --help

# Run with custom parameters
e8-orbit fit data/nilpotent_orbits.csv --max-len 20 --gamma0 0.834 --output results/
```

## Usage

### Command Line

```bash
# Basic fit
e8-orbit fit data/nilpotent_orbits.csv

# Verify with custom parameters
e8-orbit verify chain.csv --gamma0 0.834 --gamma1 0.108 --gamma2 0.0105627

# Show info
e8-orbit info
```

### Python API

```python
from e8_orbit_engine import load_orbits, build_chain, fit_gamma

# Load data
df = load_orbits('data/nilpotent_orbits.csv')

# Build chain
chain = build_chain(df, max_len=25)

# Fit gamma
results = fit_gamma(chain, gamma0_target=0.834)

print(f"γ(n) = {results['gamma0']:.4f} + {results['gamma1']:.4f}n + {results['gamma2']:.6f}n²")
print(f"R² = {results['r_squared']:.6f}")
```

### Jupyter Notebook

See `notebooks/demo.ipynb` for an interactive demonstration.

## Project Structure

```
e8_orbit_engine/
├── src/e8_orbit_engine/
│   ├── io.py          # Data loading
│   ├── chain.py       # Chain building
│   ├── fit.py         # Quadratic fitting
│   ├── verify.py      # Verification tests
│   ├── report.py      # Report generation
│   └── cli.py         # Command line interface
├── data/              # Orbit data files
├── tests/             # Unit tests
├── notebooks/         # Demo notebooks
├── reports/           # Generated reports
└── tables/            # Output tables
```

## Key Features

1. **Data-Driven Analysis**: Derives γ(n) directly from E8 orbit dimensions
2. **Calibration-Free Tests**: Validates Δ³ln(φ) = -2γ₂ without adjustable parameters
3. **Topological Verification**: Confirms γ₂ = γ₀/(8π²) relation
4. **Comprehensive Reporting**: Generates HTML reports with visualizations
5. **Robust Chain Building**: Handles monotonic orbit sequences

## Mathematical Foundation

### Normalization
- First E8 orbit step: D₀ = 248 (adjoint) → D₁ = 206 (A₄+A₁)
- s₀ = ln(248/206) anchors the normalization
- λ = 0.834/s₀ scales to target γ₀

### Quadratic Form
- γ(n) = γ₀ + γ₁n + γ₂n²
- Target: γ₀ ≈ 0.834, γ₁ ≈ 0.108, γ₂ ≈ 0.0105627

### Topological Constraint
- c₃ = 1/(8π) (topological fixpoint)
- γ₂ = γ₀/(8π²) (derived relation)

### Calibration-Free Test
- From γ(n) quadratic: Δ³ln(φₙ) = -2γ₂ (constant)
- Independent of overall scale

## Testing

```bash
# Run tests
pytest tests/

# With coverage
pytest tests/ --cov=e8_orbit_engine

# Specific test
pytest tests/test_fit.py::test_fit_real_orbits
```

## Output

The engine generates:
- **HTML Report**: Comprehensive analysis with plots
- **CSV Tables**: Orbit chain with fitted values
- **Plots**: ln(D), γ(n), residuals, third differences
- **Verification**: All consistency checks

## References

- Hamann, S. (2025). "Von Topologie zu Dynamik: Die Ordnung hinter α und den Naturkonstanten"
- E8 nilpotent orbits from Lie theory
- Topological fixpoint c₃ = 1/(8π)

## License

MIT License

## Author

Stefan Hamann (stefan.hamann@example.com)
