# E8 Cascade Model RGE Solver

This package provides a clean implementation of RGE solving for the E8 cascade model with gravity portal corrections.

## Recent Updates (v0.2.0)

### Major Fixes Implemented:

1. **Sign-preserving gravity coefficients**: Removed `abs()` calls that destroyed physical information
2. **Robust parameter indexing**: No longer assumes g1,g2,g3 are at indices 0,1,2
3. **Threshold-aware beta functions**: Beta coefficients now change at M_Σ, M_NR, M_Φ thresholds
4. **Corrected beta derivatives**: Now properly scaled with ln(10) factor for literature comparison
5. **Enhanced CLI**: Added `--cR-factor` for parameter studies

### Physical Constants Module

Created `src/e8cascade/constants.py` with:
- Threshold scales: M_Σ = 10³ GeV, M_NR = 10¹⁵ GeV, M_Φ = 10¹⁶ GeV
- Beta coefficients for each energy regime (SM → E8)
- Sign-preserving gravity coefficient calculation

## Usage

```bash
# Standard run without gravity portal
python run_e8cascade.py

# Enable gravity portal with default coefficients
python run_e8cascade.py --gravity

# Parameter study with enhanced gravity (100× factor)
python run_e8cascade.py --gravity --cR-factor 100

# Custom output directory
python run_e8cascade.py --output my_results/
```

## Results

### Without Gravity Portal
- Gauge unification at μ = 6.31×10¹⁵ GeV
- |α₁⁻¹ - α₂⁻¹| = 0.08 (excellent!)
- All couplings remain perturbative

### With Gravity Portal (100× enhanced)
- Unification shifts to μ = 5.01×10¹⁵ GeV
- |α₁⁻¹ - α₂⁻¹| = 0.10 (still very good)
- Demonstrates gravity effects on running

## Output Files

- `gauge_couplings.csv`: Full RGE solution data
- `gauge_running.png`: Evolution of gauge couplings
- `unification.png`: Detailed view near unification scale
- `stability_analysis.png`: Perturbativity and beta function analysis

## Known Issues

1. PyR@TE-generated RGEs still use SM beta coefficients
2. LaTeX rendering issue with "lambda_" in plots (cosmetic only)
3. Need proper model file with SigmaF declared as Dirac fermion

## Next Steps

1. **P1**: Regenerate RGEs with proper E8 field content and thresholds
2. **P3**: Implement adaptive ODE solver parameters
3. **P4**: Create proper Python package structure
4. **P5**: Fix plot LaTeX issues

## Technical Details

The solver implements threshold matching by dynamically adjusting beta coefficients:
- Below M_Σ: Standard Model running
- Above M_Σ: SigmaF triplet contributions
- Above M_NR: Right-handed neutrino effects  
- Above M_Φ: Full E8 cascade spectrum

Gravity portal corrections follow: c_Ri = b_i/(16π²) × (M_GUT/M_Pl)² 